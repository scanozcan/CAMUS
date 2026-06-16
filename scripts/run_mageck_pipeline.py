#!/usr/bin/env python3
"""
CAMUS - Pipeline with MAGeCK RRA wrapper.

Complete pipeline:
1. Demultiplex pooled FASTQ files by cell type (generates counts directly)
2. Combine counts into matrices per cell type
3. Run MAGeCK test (RRA) for each comparison
"""

import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from camus.config.experiment import ExperimentConfig
from camus.demux.demux_hash_table import FastDemultiplexer
from camus.wrappers.mageck_wrapper import run_celltype_comparison


def find_fastq(fastq_dir: Path, sample_name: str, read: str):
    """Return the FASTQ path for a sample/read, accepting .fastq or .fastq.gz.

    ``read`` is 'R1' or 'R2'. Returns a Path or None if neither exists.
    """
    for suffix in (".fastq", ".fastq.gz", ".fq", ".fq.gz"):
        candidate = fastq_dir / f"{sample_name}_{read}{suffix}"
        if candidate.exists():
            return candidate
    return None


def resolve_control_sgrna(config, demux_dir: Path):
    """Resolve a control sgRNA list for MAGeCK control normalization.

    Priority: explicit path in config -> auto-detected list written by the
    demuxer (``*_control_sgrnas.txt``). Returns a path string or None.
    """
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
        description='CAMUS - MAGeCK RRA Pipeline'
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
        default='pipeline_outputs',
        help='Output directory (default: pipeline_outputs)'
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
        help='FDR threshold for MAGeCK and hit flagging (default: 0.05)'
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

    return parser.parse_args()


def combine_sample_counts(
    demux_dir: Path,
    sample_names: list,
    cell_types: list,
    output_file: Path
) -> pd.DataFrame:
    """
    Combine demux count files into a single count matrix.

    Args:
        demux_dir: Directory with demux count files
        sample_names: List of sample names
        cell_types: List of cell type names
        output_file: Path to save combined matrix

    Returns:
        Combined count matrix DataFrame
    """
    # Load first sample to get sgRNA/Gene template
    first_sample = sample_names[0]
    first_file = demux_dir / f"{first_sample}_count.txt"

    if not first_file.exists():
        raise FileNotFoundError(f"Missing count file: {first_file}")

    template_df = pd.read_csv(first_file, sep='\t')

    # Start with sgRNA and Gene columns
    matrix_df = template_df[['sgRNA', 'Gene']].copy()
    matrix_df['in_library'] = 1  # All sgRNAs are in library

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

    print(f"  Combined {len(sample_names)} samples × {len(cell_types)} cell types")
    print(f"  Output: {output_file.name}")

    return matrix_df


def extract_celltype_matrix(
    combined_matrix: pd.DataFrame,
    cell_type: str,
    sample_names: list,
    output_file: Path,
    library_file: Path = None
) -> pd.DataFrame:
    """
    Extract counts for a specific cell type into MAGeCK format.

    IMPORTANT: Filters to only include sgRNAs from this cell type's library,
    removing sgRNAs from other cell types that have all-zero counts.

    Args:
        combined_matrix: Combined count matrix
        cell_type: Cell type to extract
        sample_names: List of sample names
        output_file: Path to save cell-type-specific matrix
        library_file: Optional path to cell-type-specific library file for filtering

    Returns:
        Cell-type-specific count matrix
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

    # Set in_library based on whether sgRNA is in this cell type's library
    if library_file is not None and library_file.exists():
        # Load cell-type-specific library (no header, first column is sgRNA)
        lib_df = pd.read_csv(library_file, sep='\t', header=None, names=['sgRNA', 'Gene'])
        lib_sgrnas = set(lib_df['sgRNA'].values)

        # Set in_library: 1 if in this library, 0 otherwise
        celltype_df['in_library'] = celltype_df['sgRNA'].isin(lib_sgrnas).astype(int)

        # Filter to only keep sgRNAs from this library
        celltype_df = celltype_df[celltype_df['in_library'] == 1].copy()

        n_in_library = (celltype_df['in_library'] == 1).sum()
        n_filtered = len(combined_matrix) - len(celltype_df)
        print(f"    Filtered to {n_in_library:,} sgRNAs from library (removed {n_filtered:,} from other cell types)")
    else:
        # Fallback: filter out rows where all counts are zero
        count_cols = [col for col in celltype_df.columns if col not in ['sgRNA', 'Gene', 'in_library']]
        row_sums = celltype_df[count_cols].sum(axis=1)
        celltype_df = celltype_df[row_sums > 0].copy()

        n_filtered = len(combined_matrix) - len(celltype_df)
        print(f"    Filtered to {len(celltype_df):,} sgRNAs with non-zero counts (removed {n_filtered:,} all-zero rows)")

    # Save
    celltype_df.to_csv(output_file, sep='\t', index=False)

    return celltype_df


def write_group_barcode_csv(cell_types, output_path: Path) -> Path:
    """Write a Celltype,Barcode CSV for the demuxer from CellTypeConfig objects.

    Restricting the demuxer to only a group's cell-type barcodes means pure
    groups (e.g. a single-cell-type culture) don't get spurious near-empty lanes
    for cell types they don't contain.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write("Celltype,Barcode\n")
        for ct in cell_types:
            f.write(f"{ct.name},{ct.barcode}\n")
    return output_path


def run_group(group, config, data_dir: Path, output_root: Path, args) -> int:
    """Demux + count + RRA for one group, compared against its OWN control baseline.

    Each group writes a self-contained ``<output_root>/<group>/{demux,counts,results}``
    subtree (the top level itself in flat mode, where ``group.name`` is empty), and
    runs ``{control}_vs_{condition}`` per cell type using the group's own control
    samples. Returns the number of completed comparisons.
    """
    group_out = output_root / group.name if group.name else output_root
    demux_dir = group_out / 'demux'
    counts_dir = group_out / 'counts'
    results_dir = group_out / 'results'
    for d in (demux_dir, counts_dir, results_dir):
        d.mkdir(parents=True, exist_ok=True)

    cell_types = config.get_group_cell_types(group)
    cell_type_names = [ct.name for ct in cell_types]
    all_samples = group.get_sample_names()
    label = group.name or config.name

    print("\n" + "#" * 80)
    print(f"GROUP: {label}   cell types: {', '.join(cell_type_names)}")
    print(f"  control: {group.control.name} ({group.control.replicates} reps)"
          f"   conditions: {', '.join(c.name for c in group.conditions)}")
    print("#" * 80)

    # ---- STEP 1: demultiplex (restricted to this group's barcodes) ----
    if not args.skip_demux:
        print("\n[1] Demultiplexing by cell type ...")
        g = config.geometry
        barcode_csv = write_group_barcode_csv(
            cell_types, demux_dir / "celltype_barcodes.csv")
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
        fastq_dir = data_dir / 'fastq_files'
        library_file = data_dir / "combined_library.txt"
        for sample_name in all_samples:
            r1_file = find_fastq(fastq_dir, sample_name, 'R1')
            r2_file = find_fastq(fastq_dir, sample_name, 'R2')
            if r1_file is None or r2_file is None:
                print(f"  Warning: missing FASTQ for {sample_name} "
                      f"(looked for .fastq/.fastq.gz), skipping...")
                continue
            print(f"  Demultiplexing {sample_name} ...")
            demuxer.process(
                read1_file=str(r1_file),
                read2_file=str(r2_file),
                library_file=str(library_file),
                output_prefix=str(demux_dir / sample_name),
            )
    else:
        print("\n[1] Skipping demultiplexing (using existing files)")

    # ---- STEP 2: per-cell-type count matrices ----
    print("\n[2] Building cell-type count matrices ...")
    combined_matrix = combine_sample_counts(
        demux_dir=demux_dir,
        sample_names=all_samples,
        cell_types=cell_type_names,
        output_file=counts_dir / "all_samples_combined.txt",
    )
    for cell_type in cell_types:
        print(f"  [{cell_type.name}]")
        extract_celltype_matrix(
            combined_matrix=combined_matrix,
            cell_type=cell_type.name,
            sample_names=all_samples,
            output_file=counts_dir / f"{cell_type.name}_count_matrix.txt",
            library_file=data_dir / f"{cell_type.name.lower()}_library.txt",
        )

    # ---- STEP 3: MAGeCK RRA, each condition vs THIS group's own control ----
    print("\n[3] Running MAGeCK RRA ...")
    control_sgrna = None
    if not args.no_control_sgrna:
        control_sgrna = resolve_control_sgrna(config, demux_dir)
        if control_sgrna:
            print(f"  Control sgRNA list: {control_sgrna}")
        else:
            print("  No control sgRNA list found; proceeding without --control-sgrna.")

    control_samples = group.control_sample_names()
    comparison_count = 0
    for cell_type in cell_types:
        matrix_file = counts_dir / f"{cell_type.name}_count_matrix.txt"
        if not matrix_file.exists():
            print(f"  Error: Count matrix not found: {matrix_file}")
            continue
        for condition in group.conditions:
            comparison_name = f"{group.control.name}_vs_{condition.name}"
            treatment_samples = group.condition_sample_names(condition)
            output_prefix = results_dir / f"{comparison_name}_{cell_type.name}"

            print(f"\n  [{cell_type.name}] {comparison_name}")
            print(f"    control:   {', '.join(control_samples)}")
            print(f"    treatment: {', '.join(treatment_samples)}")
            try:
                results = run_celltype_comparison(
                    count_matrix_file=str(matrix_file),
                    control_samples=control_samples,
                    treatment_samples=treatment_samples,
                    output_prefix=str(output_prefix),
                    conda_env=args.conda_env,
                    mageck_binary=args.mageck_binary,
                    normalization=args.normalization,
                    control_sgrna=control_sgrna,
                    fdr_threshold=args.fdr,
                )
                gene_output = results_dir / f"{comparison_name}_{cell_type.name}_gene_results.txt"
                results['gene_results'].to_csv(gene_output, sep='\t', index=False)
                print(f"    ✓ {gene_output.name}")
                comparison_count += 1
            except Exception as e:
                print(f"    ERROR: {e}")
                import traceback
                traceback.print_exc()
                continue

    return comparison_count


def main():
    """Main execution."""
    args = parse_args()

    print("=" * 80)
    print("CAMUS - MAGeCK RRA PIPELINE")
    print("=" * 80)

    # Load config
    print(f"\nLoading configuration: {args.config}")
    config = ExperimentConfig.from_yaml(args.config)

    data_dir = Path(args.data_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True, parents=True)

    groups = config.get_groups()
    grouped = bool(config.groups)
    print(f"\n{len(groups)} group(s): "
          f"{', '.join(g.name or '(flat)' for g in groups)}")

    total = 0
    for group in groups:
        total += run_group(group, config, data_dir, output_dir, args)

    # ========================================================================
    # SUMMARY
    # ========================================================================
    print("\n" + "=" * 80)
    print("✓ PIPELINE COMPLETE!")
    print("=" * 80)
    print(f"\nCompleted {total} comparisons across {len(groups)} group(s)")
    print(f"Output directory: {output_dir}")
    if grouped:
        print("Per-group results under <output>/<group>/results/")

    print("\nNext step (visualization):")
    print("  python scripts/post_analysis_visualization.py \\")
    if grouped:
        print(f"    --results-dir {output_dir} \\   # root containing per-group subdirs")
    else:
        print(f"    --results-dir {output_dir / 'results'} \\")
    print(f"    --config {args.config} \\")
    print(f"    --output visualizations")


if __name__ == '__main__':
    main()
