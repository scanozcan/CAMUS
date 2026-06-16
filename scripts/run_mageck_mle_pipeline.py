#!/usr/bin/env python3
"""
CAMUS - Pipeline with MAGeCK MLE for multi-condition comparison.

This pipeline runs MAGeCK MLE (Maximum Likelihood Estimation) to compare
multiple conditions simultaneously against a control group, estimating
individual effects for each condition.

Complete pipeline:
1. Demultiplex pooled FASTQ files by cell type
2. Combine counts into matrices per cell type
3. Run MAGeCK MLE for each cell type (all conditions vs control)
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from camus.config.experiment import ExperimentConfig
from camus.demux.demux_hash_table import FastDemultiplexer
from camus.wrappers.mageck_wrapper import run_mle_multi_condition


def find_fastq(fastq_dir: Path, sample_name: str, read: str):
    """Return the FASTQ path for a sample/read, accepting .fastq or .fastq.gz."""
    for suffix in (".fastq", ".fastq.gz", ".fq", ".fq.gz"):
        candidate = fastq_dir / f"{sample_name}_{read}{suffix}"
        if candidate.exists():
            return candidate
    return None


def resolve_control_sgrna(config, demux_dir: Path):
    """Resolve a control sgRNA list (config path, else auto-detected demux list)."""
    if getattr(config, 'control_sgrna', None):
        p = Path(config.control_sgrna)
        if p.exists():
            return str(p)
        print(f"  Warning: configured control_sgrna not found: {p}")
    matches = sorted(demux_dir.glob('*_control_sgrnas.txt'))
    if matches:
        return str(matches[0])
    return None


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='CAMUS - MAGeCK MLE Multi-Condition Pipeline'
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to experiment YAML configuration'
    )

    parser.add_argument(
        '--data-dir',
        type=str,
        required=True,
        help='Directory containing test data (library files + FASTQ)'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='pipeline_outputs_mle',
        help='Output directory (default: pipeline_outputs_mle)'
    )

    parser.add_argument(
        '--conda-env',
        type=str,
        default=None,
        help='Conda environment with MAGeCK installed (fallback if not on PATH)'
    )

    parser.add_argument(
        '--mageck-binary',
        type=str,
        default=None,
        help='Explicit path to the mageck executable (highest priority)'
    )

    parser.add_argument(
        '--normalization',
        type=str,
        default='median',
        choices=['median', 'total', 'control'],
        help='MAGeCK normalization method (default: median)'
    )

    parser.add_argument(
        '--fdr',
        type=float,
        default=0.05,
        help='FDR threshold for hit flagging (default: 0.05)'
    )

    parser.add_argument(
        '--no-control-sgrna',
        action='store_true',
        help='Do not pass non-targeting controls to MAGeCK even if detected'
    )

    parser.add_argument(
        '--skip-demux',
        action='store_true',
        help='Skip demultiplexing step (use existing count files)'
    )

    parser.add_argument(
        '--use-existing-counts',
        type=str,
        help='Path to existing counts directory (skips demux and count steps)'
    )

    return parser.parse_args()


def combine_sample_counts(
    demux_dir: Path,
    sample_names: list,
    cell_types: list,
    output_file: Path
) -> pd.DataFrame:
    """
    Combine demux count files into a single count matrix.
    """
    # Load first sample to get sgRNA/Gene template
    first_sample = sample_names[0]
    first_file = demux_dir / f"{first_sample}_count.txt"

    if not first_file.exists():
        raise FileNotFoundError(f"Missing count file: {first_file}")

    template_df = pd.read_csv(first_file, sep='\t')

    # Start with sgRNA and Gene columns
    matrix_df = template_df[['sgRNA', 'Gene']].copy()
    matrix_df['in_library'] = 1

    # Add counts for each sample
    for sample_name in sample_names:
        count_file = demux_dir / f"{sample_name}_count.txt"

        if not count_file.exists():
            print(f"  Warning: Missing {count_file.name}, using zeros")
            for cell_type in cell_types:
                matrix_df[f"{sample_name}_{cell_type}"] = 0
            continue

        count_df = pd.read_csv(count_file, sep='\t')

        # Add columns for each cell type in this sample
        for cell_type in cell_types:
            if cell_type in count_df.columns:
                matrix_df[f"{sample_name}_{cell_type}"] = count_df[cell_type].values
            else:
                matrix_df[f"{sample_name}_{cell_type}"] = 0

    # Save combined matrix
    matrix_df.to_csv(output_file, sep='\t', index=False)
    return matrix_df


def extract_celltype_matrix(
    combined_matrix: pd.DataFrame,
    cell_type: str,
    sample_names: list,
    output_file: Path
) -> pd.DataFrame:
    """
    Extract counts for a specific cell type into MAGeCK format.
    """
    # Start with sgRNA, Gene, in_library
    celltype_df = combined_matrix[['sgRNA', 'Gene', 'in_library']].copy()

    # Add sample columns for this cell type
    for sample_name in sample_names:
        col_name = f"{sample_name}_{cell_type}"
        if col_name in combined_matrix.columns:
            celltype_df[sample_name] = combined_matrix[col_name].values
        else:
            celltype_df[sample_name] = 0

    # Save
    celltype_df.to_csv(output_file, sep='\t', index=False)
    return celltype_df


def main():
    """Main execution."""
    args = parse_args()

    print("=" * 80)
    print("CAMUS - MAGeCK MLE MULTI-CONDITION PIPELINE")
    print("=" * 80)

    # Load config
    print(f"\nLoading configuration: {args.config}")
    config = ExperimentConfig.from_yaml(args.config)

    # Setup directories
    data_dir = Path(args.data_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True, parents=True)

    demux_dir = output_dir / 'demux'
    counts_dir = output_dir / 'counts'
    results_dir = output_dir / 'results_mle'

    demux_dir.mkdir(exist_ok=True)
    counts_dir.mkdir(exist_ok=True)
    results_dir.mkdir(exist_ok=True)

    # Get all sample names
    all_samples = []
    sample_groups = [config.control] + config.conditions
    for sample_group in sample_groups:
        for rep_idx in range(1, sample_group.replicates + 1):
            all_samples.append(f"{sample_group.name}_Rep{rep_idx}")

    # ========================================================================
    # OPTION: Use existing counts
    # ========================================================================
    if args.use_existing_counts:
        print("\n⊳ Using existing counts from:", args.use_existing_counts)
        counts_dir = Path(args.use_existing_counts)

    # ========================================================================
    # STEP 1: DEMULTIPLEXING (if needed)
    # ========================================================================
    elif not args.skip_demux:
        print("\n" + "=" * 80)
        print("STEP 1: DEMULTIPLEXING BY CELL TYPE")
        print("=" * 80)

        # Initialize demuxer from the experiment read geometry
        g = config.geometry
        barcode_csv = data_dir / "celltype_barcodes.csv"
        demuxer = FastDemultiplexer(
            barcode_csv=str(barcode_csv),
            barcode_start=g.barcode_start,
            barcode_length=g.barcode_length,
            grna_start=g.grna_start,
            grna_length=g.grna_length,
            max_barcode_mismatches=g.max_barcode_mismatches,
            allow_grna_mismatch=g.allow_grna_mismatch,
            min_qual=g.min_qual,
            grna_anchor=g.grna_anchor,
            anchor_max_offset=g.anchor_max_offset,
        )

        # Process all samples
        fastq_dir = data_dir / 'fastq_files'
        for sample_name in all_samples:
            r1_file = find_fastq(fastq_dir, sample_name, 'R1')
            r2_file = find_fastq(fastq_dir, sample_name, 'R2')

            if r1_file is None or r2_file is None:
                print(f"\nWarning: Missing FASTQ files for {sample_name} "
                      f"(looked for .fastq/.fastq.gz), skipping...")
                continue

            print(f"\nDemultiplexing {sample_name}...")

            library_file = data_dir / "combined_library.txt"
            output_prefix = demux_dir / sample_name

            demuxer.process(
                read1_file=str(r1_file),
                read2_file=str(r2_file),
                library_file=str(library_file),
                output_prefix=str(output_prefix)
            )

        print("\n✓ Demultiplexing complete!")

        # ====================================================================
        # STEP 2: CREATE COUNT MATRICES
        # ====================================================================
        print("\n" + "=" * 80)
        print("STEP 2: CREATING CELL-TYPE-SPECIFIC COUNT MATRICES")
        print("=" * 80)

        cell_type_names = [ct.name for ct in config.cell_types]

        print("\nCombining all sample counts...")
        combined_matrix_file = counts_dir / "all_samples_combined.txt"

        combined_matrix = combine_sample_counts(
            demux_dir=demux_dir,
            sample_names=all_samples,
            cell_types=cell_type_names,
            output_file=combined_matrix_file
        )

        print("\nExtracting cell-type-specific matrices...")
        for cell_type in config.cell_types:
            matrix_file = counts_dir / f"{cell_type.name}_count_matrix.txt"

            print(f"\n  [{cell_type.name}]")
            extract_celltype_matrix(
                combined_matrix=combined_matrix,
                cell_type=cell_type.name,
                sample_names=all_samples,
                output_file=matrix_file
            )

        print("\n✓ Count matrix creation complete!")

    # ========================================================================
    # STEP 3: RUN MAGeCK MLE
    # ========================================================================
    print("\n" + "=" * 80)
    print("STEP 3: RUNNING MAGeCK MLE ANALYSIS")
    print("=" * 80)
    print("\nMLE compares ALL conditions simultaneously vs control.")
    print("Each condition's effect is estimated independently.\n")

    # Resolve a non-targeting control sgRNA list for control normalization.
    control_sgrna = None
    if not args.no_control_sgrna and demux_dir.exists():
        control_sgrna = resolve_control_sgrna(config, demux_dir)
        if control_sgrna:
            print(f"Using control sgRNA list: {control_sgrna}\n")

    analysis_count = 0

    # Prepare condition groups
    control_sample_names = [
        f"Control_Rep{i}" for i in range(1, config.control.replicates + 1)
    ]

    condition_groups = {}
    for condition in config.conditions:
        condition_sample_names = [
            f"{condition.name}_Rep{i}" for i in range(1, condition.replicates + 1)
        ]
        condition_groups[condition.name] = condition_sample_names

    # Run MLE for each cell type
    for cell_type in config.cell_types:
        print(f"\n[{cell_type.name}] Running MLE...")

        matrix_file = counts_dir / f"{cell_type.name}_count_matrix.txt"

        if not matrix_file.exists():
            print(f"  Error: Count matrix not found: {matrix_file}")
            continue

        output_prefix = results_dir / f"{cell_type.name}_MLE"

        print(f"  Control samples: {', '.join(control_sample_names)}")
        for cond_name, cond_samples in condition_groups.items():
            print(f"  {cond_name}: {', '.join(cond_samples)}")

        try:
            results = run_mle_multi_condition(
                count_matrix_file=str(matrix_file),
                control_samples=control_sample_names,
                condition_groups=condition_groups,
                output_prefix=str(output_prefix),
                conda_env=args.conda_env,
                mageck_binary=args.mageck_binary,
                normalization=args.normalization,
                control_sgrna=control_sgrna,
                fdr_threshold=args.fdr,
            )

            # Save parsed gene results
            gene_output = results_dir / f"{cell_type.name}_MLE_gene_results.txt"
            results['gene_results'].to_csv(gene_output, sep='\t', index=False)

            print(f"  ✓ Gene results saved: {gene_output.name}")
            print(f"  ✓ Design matrix saved: {Path(results['design_matrix']).name}")

            analysis_count += 1

        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            continue

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "=" * 80)
    print("✓ MLE PIPELINE COMPLETE!")
    print("=" * 80)

    print(f"\nCompleted {analysis_count} cell type analyses")
    print(f"Output directory: {output_dir}")

    print("\nGenerated files:")
    if not args.use_existing_counts and not args.skip_demux:
        print(f"  Demux: {len(list(demux_dir.glob('*_count.txt')))} count files")
        print(f"  Counts: {len(list(counts_dir.glob('*_count_matrix.txt')))} count matrices")
    print(f"  Results: {len(list(results_dir.glob('*.gene_summary.txt')))} MAGeCK MLE gene summaries")
    print(f"  Results: {len(list(results_dir.glob('*_gene_results.txt')))} parsed gene results")
    print(f"  Design matrices: {len(list(results_dir.glob('*_design_matrix.txt')))} files")

    print("\nMLE Results:")
    print("  Each gene has estimates for:")
    for condition in config.conditions:
        print(f"    - {condition.name}: beta (effect size), z-score, p-value, FDR")
    print("\n  Conditions are compared simultaneously, sharing information")
    print("  for better variance estimation.\n")

    print("Next steps:")
    print("  - Review gene_summary.txt files for condition-specific effects")
    print("  - Compare MLE vs RRA results to validate findings")
    print("  - Use MLE for final publication (more statistically rigorous)")


if __name__ == '__main__':
    main()
