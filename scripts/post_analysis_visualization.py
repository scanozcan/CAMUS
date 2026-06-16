#!/usr/bin/env python3
"""
Comprehensive post-analysis visualization for CAMUS (MAGeCK wrapper).

Generates publication-quality plots and reports:
- Volcano plots (per cell type, per condition)
- QC plots (PCA, correlations, coverage)
- Cross-condition comparisons (Venn diagrams, concordance)
- Gene heatmaps
- Text summary reports

Adapted for MAGeCK output format.
"""

import sys
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from camus.config.experiment import ExperimentConfig
from camus.visualization import (
    VolcanoPlotter,
    QCPlotter,
    ComparisonPlotter,
    HeatmapPlotter,
    ReportGenerator
)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='CAMUS - Post-Analysis Visualization (MAGeCK)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help='Directory containing MAGeCK analysis results'
    )

    parser.add_argument(
        '--config',
        type=str,
        required=True,
        help='Path to experiment YAML configuration'
    )

    parser.add_argument(
        '--output',
        type=str,
        default='visualizations',
        help='Output directory for visualizations (default: visualizations)'
    )

    parser.add_argument(
        '--fdr',
        type=float,
        default=0.05,
        help='FDR threshold (default: 0.05)'
    )

    return parser.parse_args()


def load_mageck_results(results_dir: Path, config: ExperimentConfig) -> Dict[str, List[Dict]]:
    """
    Load MAGeCK analysis results (parsed gene_results.txt files).

    Args:
        results_dir: Directory with result files
        config: Experiment configuration

    Returns:
        Dictionary mapping comparison names to result lists
    """
    all_results = {}

    # Find all gene result files
    for condition in config.conditions:
        comparison_name = f"{config.control.name}_vs_{condition.name}"
        results_list = []

        for cell_type in config.cell_types:
            # Look for parsed gene-level results from MAGeCK wrapper
            gene_file = results_dir / f"{comparison_name}_{cell_type.name}_gene_results.txt"

            if not gene_file.exists():
                print(f"  Warning: Missing {gene_file.name}")
                continue

            try:
                gene_results = pd.read_csv(gene_file, sep='\t')

                # Verify required columns exist
                required_cols = ['Gene', 'FDR', 'logFC', 'Direction', 'Num_sgRNAs']
                missing_cols = [col for col in required_cols if col not in gene_results.columns]

                if missing_cols:
                    print(f"  Warning: {gene_file.name} missing columns: {missing_cols}")
                    continue

                # Convert Direction to match v3 format (lowercase)
                if 'Direction' in gene_results.columns:
                    gene_results['Direction'] = gene_results['Direction'].str.lower()

                results_list.append({
                    'cell_type': cell_type.name,
                    'gene_results': gene_results
                })

            except Exception as e:
                print(f"  Error loading {gene_file.name}: {e}")

        if results_list:
            all_results[comparison_name] = results_list

    return all_results


def load_count_matrices(counts_dir: Path, cell_types: List) -> Dict[str, pd.DataFrame]:
    """
    Load count matrices for QC plots.

    Args:
        counts_dir: Directory with count files
        cell_types: List of cell type objects

    Returns:
        Dictionary mapping cell type names to count DataFrames
    """
    count_matrices = {}

    for cell_type in cell_types:
        count_file = counts_dir / f"{cell_type.name}_count_matrix.txt"

        if count_file.exists():
            try:
                count_df = pd.read_csv(count_file, sep='\t')
                # Convert count columns (all except first 3: sgRNA, Gene, in_library) to numeric
                count_cols = count_df.columns[3:]
                for col in count_cols:
                    count_df[col] = pd.to_numeric(count_df[col], errors='coerce')
                count_matrices[cell_type.name] = count_df
            except Exception as e:
                print(f"  Error loading {count_file.name}: {e}")

    return count_matrices


def load_celltype_proportions(demux_dir: Path) -> Dict[str, pd.DataFrame]:
    """
    Load cell type proportion data from demux outputs.

    Args:
        demux_dir: Directory with demux output files

    Returns:
        Dictionary mapping sample names to celltype proportion DataFrames
    """
    celltype_data = {}

    if not demux_dir.exists():
        return celltype_data

    # Find all celltype_proportions.txt files
    for prop_file in demux_dir.glob('*_celltype_proportions.txt'):
        sample_name = prop_file.stem.replace('_celltype_proportions', '')
        try:
            df = pd.read_csv(prop_file, sep='\t')
            celltype_data[sample_name] = df
        except Exception as e:
            print(f"  Warning: Could not load {prop_file.name}: {e}")

    return celltype_data


def visualize_unit(config, results_dir, output_dir, fdr):
    """Generate all visualizations for one analysis unit (a flat run or one group).

    ``results_dir`` is the MAGeCK results directory; sibling ``counts/`` and
    ``demux/`` dirs are read for QC and cell-type proportions. ``config`` is a
    view exposing ``.control``, ``.conditions`` and ``.cell_types`` for this unit.
    Returns the number of plots generated (0 if no results found).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)

    results_dir = Path(results_dir)
    counts_dir = results_dir.parent / 'counts'

    # Load results
    print("\nLoading MAGeCK analysis results...")
    all_results = load_mageck_results(results_dir, config)

    if not all_results:
        print(f"  WARNING: No results found in {results_dir}")
        return 0

    print(f"  Loaded {len(all_results)} comparisons")

    # Load count matrices for QC
    print("\nLoading count matrices...")
    count_matrices = load_count_matrices(counts_dir, config.cell_types)
    print(f"  Loaded {len(count_matrices)} count matrices")

    # Load cell type proportions from demux
    demux_dir = results_dir.parent / 'demux'
    print("\nLoading cell type proportions...")
    celltype_data = load_celltype_proportions(demux_dir)
    print(f"  Loaded {len(celltype_data)} sample cell type distributions")

    # Initialize plotters (using v3 visualization modules)
    volcano_plotter = VolcanoPlotter(output_dir, fdr_threshold=fdr)
    qc_plotter = QCPlotter(output_dir)
    comparison_plotter = ComparisonPlotter(output_dir, fdr_threshold=fdr)
    heatmap_plotter = HeatmapPlotter(output_dir, fdr_threshold=fdr)
    report_generator = ReportGenerator(output_dir, fdr_threshold=fdr)

    print("\n" + "=" * 80)
    print("GENERATING VISUALIZATIONS")
    print("=" * 80)

    plot_count = 0

    # 1. Volcano plots
    print("\n[1] Generating volcano plots...")
    for comparison_name, results_list in all_results.items():
        for result in results_list:
            volcano_plotter.plot_volcano(
                result['gene_results'],
                result['cell_type'],
                comparison_name
            )
            plot_count += 1

        # Multi-panel volcano
        volcano_plotter.plot_multi_panel_volcano(
            all_results,
            comparison_name,
            [ct.name for ct in config.cell_types]
        )
        plot_count += 1

    print(f"  Generated {plot_count} volcano plots")

    # 2. QC plots
    if count_matrices:
        print("\n[2] Generating QC plots...")
        qc_count = 0

        for cell_type, count_df in count_matrices.items():
            # Coverage distribution
            sample_cols = [c for c in count_df.columns if c not in ['sgRNA', 'Gene', 'in_library']]
            if len(sample_cols) > 0:
                # Group samples by condition (baseline = this unit's control name)
                sample_groups = {'Control': [], 'Conditions': []}
                for sample in sample_cols:
                    if config.control.name in sample:
                        sample_groups['Control'].append(sample)
                    else:
                        sample_groups['Conditions'].append(sample)

                # Coverage distribution
                qc_plotter.plot_coverage_distribution(
                    count_df[sample_cols],
                    sample_groups
                )
                qc_count += 1

                # Count distributions
                qc_plotter.plot_count_distributions(count_df[sample_cols], cell_type)
                qc_count += 1

                # Replicate correlations
                if len(sample_cols) > 1:
                    qc_plotter.plot_replicate_correlations(
                        count_df[sample_cols],
                        sample_groups,
                        title=f"{cell_type} Replicate Correlations"
                    )
                    qc_count += 1

        print(f"  Generated {qc_count} QC plots")
        plot_count += qc_count

    # 2b. Cell type distribution plots
    if celltype_data:
        print("\n[2b] Generating cell type distribution plots...")
        celltype_count = 0

        # Get list of cell types
        cell_type_names = [ct.name for ct in config.cell_types]

        # Cell type proportions bar plot
        qc_plotter.plot_celltype_proportions(celltype_data, cell_type_names)
        celltype_count += 1

        # Cell type pie charts
        qc_plotter.plot_celltype_pie_charts(celltype_data, cell_type_names)
        celltype_count += 1

        # Replicate-aware cell type plot
        qc_plotter.plot_celltype_replicate_aware(celltype_data, cell_type_names)
        celltype_count += 1

        print(f"  Generated {celltype_count} cell type distribution plots")
        plot_count += celltype_count

    # 3. Cross-condition comparisons
    print("\n[3] Generating cross-condition comparisons...")
    comp_count = 0

    for cell_type in config.cell_types:
        # Collect hits per condition
        hits_per_condition = {}

        for comparison_name, results_list in all_results.items():
            for result in results_list:
                if result['cell_type'] == cell_type.name:
                    gene_results = result['gene_results']
                    sig_genes = set(
                        gene_results[gene_results['FDR'] <= fdr]['Gene'].values
                    )
                    condition_name = comparison_name.replace(
                        f"{config.control.name}_vs_", '')
                    hits_per_condition[condition_name] = sig_genes

        if len(hits_per_condition) >= 2:
            # Venn diagram
            comparison_plotter.plot_venn_diagram(
                hits_per_condition,
                cell_type.name
            )
            comp_count += 1

            # LogFC comparisons (pairwise)
            if len(hits_per_condition) >= 2:
                conditions = list(hits_per_condition.keys())
                for i in range(len(conditions)):
                    for j in range(i+1, len(conditions)):
                        cond1 = conditions[i]
                        cond2 = conditions[j]

                        # Find results for these conditions
                        results1 = None
                        results2 = None

                        for comp_name, results_list in all_results.items():
                            if cond1 in comp_name:
                                for r in results_list:
                                    if r['cell_type'] == cell_type.name:
                                        results1 = r['gene_results']
                            if cond2 in comp_name:
                                for r in results_list:
                                    if r['cell_type'] == cell_type.name:
                                        results2 = r['gene_results']

                        if results1 is not None and results2 is not None:
                            comparison_plotter.plot_logfc_comparison(
                                results1, results2, cond1, cond2, cell_type.name
                            )
                            comp_count += 1

    print(f"  Generated {comp_count} comparison plots")
    plot_count += comp_count

    # 4. Heatmaps
    print("\n[4] Generating heatmaps...")
    heatmap_count = 0

    # Reorganize data for cell type comparison heatmaps
    # Format: {cell_type: {condition: results_df}}
    results_by_celltype = {}
    for comparison_name, results_list in all_results.items():
        for result in results_list:
            cell_type = result['cell_type']
            if cell_type not in results_by_celltype:
                results_by_celltype[cell_type] = {}
            results_by_celltype[cell_type][comparison_name] = result['gene_results']

    # Cell-type-specific heatmaps for each condition
    for comparison_name in all_results.keys():
        heatmap_plotter.plot_celltype_heatmap(results_by_celltype, comparison_name, top_n=40)
        heatmap_count += 1

    # Overall clustered heatmap
    all_comparison_results = {}
    for comparison_name, results_list in all_results.items():
        for result in results_list:
            key = f"{comparison_name}_{result['cell_type']}"
            all_comparison_results[key] = result['gene_results']

    if all_comparison_results:
        heatmap_plotter.plot_clustered_heatmap(all_comparison_results, min_conditions=2)
        heatmap_count += 1

    # Dotplot showing logFC and FDR
    heatmap_plotter.plot_dotplot(results_by_celltype, top_n=50)
    heatmap_count += 1

    print(f"  Generated {heatmap_count} heatmaps")
    plot_count += heatmap_count

    # 5. Generate reports
    print("\n[5] Generating text reports...")

    # Summary report
    report_generator.generate_summary_report(all_results, config)

    # QC report
    if count_matrices:
        report_generator.generate_qc_report(
            count_matrices,
            [ct.name for ct in config.cell_types]
        )

    print("  Generated summary and QC reports")

    print(f"\n  Unit complete: {plot_count} plots, "
          f"{len(list(output_dir.glob('*.png')))} PNGs in {output_dir}")
    return plot_count


def cross_group_comparisons(config, results_root, output_dir, fdr):
    """Compare the SAME cell type ACROSS groups (e.g. prolif alone vs co-culture).

    The per-group visualizations only compare conditions *within* a group. This
    step answers the cross-group question — e.g. "prolif in the pure G1 culture
    vs prolif demultiplexed from the G3 co-culture". For each cell type present
    in >=2 groups it writes, under ``<output>/cross_group/``:
      - a logFC concordance scatter per group pair (cell-autonomous hits sit on
        the diagonal; context-dependent / rescue genes fall off it),
      - a hit-overlap Venn across groups,
      - ``CROSS_GROUP_REPORT.md`` flagging genes significant in one group only.
    """
    import itertools

    results_root = Path(results_root)
    out = Path(output_dir) / "cross_group"
    out.mkdir(parents=True, exist_ok=True)

    # Gather {cell_type: {group_label: gene_results_df}}
    per_celltype: Dict[str, Dict[str, pd.DataFrame]] = {}
    for group in config.get_groups():
        for ct in config.get_group_cell_types(group):
            for condition in group.conditions:
                comp = f"{group.control.name}_vs_{condition.name}"
                gene_file = (results_root / group.name / 'results'
                             / f"{comp}_{ct.name}_gene_results.txt")
                if not gene_file.exists():
                    continue
                label = (group.name if len(group.conditions) == 1
                         else f"{group.name}_{condition.name}")
                per_celltype.setdefault(ct.name, {})[label] = pd.read_csv(gene_file, sep='\t')

    plotter = ComparisonPlotter(out, fdr_threshold=fdr)
    n_plots = 0
    report = ["# Cross-group comparisons (same cell type across groups)\n",
              f"FDR threshold: {fdr}\n"]

    for ct, by_label in per_celltype.items():
        labels = list(by_label.keys())
        if len(labels) < 2:
            continue  # nothing to compare for a cell type in a single group
        report.append(f"\n## {ct}\n")

        # Hit overlap across all groups containing this cell type
        hits = {lab: set(df[df['FDR'] <= fdr]['Gene']) for lab, df in by_label.items()}
        plotter.plot_venn_diagram(hits, ct)
        n_plots += 1

        # Pairwise logFC concordance + context-dependent gene report
        for a, b in itertools.combinations(labels, 2):
            m = by_label[a][['Gene', 'logFC', 'FDR']].merge(
                by_label[b][['Gene', 'logFC', 'FDR']], on='Gene', suffixes=(f'_{a}', f'_{b}'))
            if m.empty:
                plotter.plot_logfc_comparison(by_label[a], by_label[b], a, b, ct)
                n_plots += 1
                continue
            r = (np.corrcoef(m[f'logFC_{a}'], m[f'logFC_{b}'])[0, 1]
                 if len(m) > 2 else float('nan'))
            sig_a = m[f'FDR_{a}'] <= fdr
            sig_b = m[f'FDR_{b}'] <= fdr
            context = m[sig_a ^ sig_b].copy()  # significant in exactly one group
            context['abs_delta'] = (m[f'logFC_{a}'] - m[f'logFC_{b}']).abs()
            context = context.sort_values('abs_delta', ascending=False)

            # Label genes with a real logFC swing (off the diagonal), not FDR noise.
            label_genes = context[context['abs_delta'] >= 1.0]['Gene'].tolist()
            plotter.plot_logfc_comparison(
                by_label[a], by_label[b], a, b, ct, label_genes=label_genes)
            n_plots += 1
            report.append(
                f"- **{a} vs {b}**: logFC Pearson r = {r:.3f}; "
                f"{int((sig_a & sig_b).sum())} shared hits, "
                f"{len(context)} context-dependent (significant in one group only)")
            if len(context):
                report.append(f"\n  | gene | logFC {a} (FDR) | logFC {b} (FDR) |")
                report.append("  |---|---|---|")
                for _, row in context.head(15).iterrows():
                    report.append(
                        f"  | {row['Gene']} | {row[f'logFC_{a}']:+.2f} ({row[f'FDR_{a}']:.2g}) "
                        f"| {row[f'logFC_{b}']:+.2f} ({row[f'FDR_{b}']:.2g}) |")

    (out / "CROSS_GROUP_REPORT.md").write_text("\n".join(report) + "\n")
    print(f"\n[cross-group] {n_plots} plots + CROSS_GROUP_REPORT.md in {out}")
    return n_plots


def main():
    """Main execution — group-aware over the experiment config."""
    args = parse_args()

    print("=" * 80)
    print("CAMUS - POST-ANALYSIS VISUALIZATION (MAGeCK)")
    print("=" * 80)

    print(f"\nLoading configuration: {args.config}")
    config = ExperimentConfig.from_yaml(args.config)

    base_out = Path(args.output)
    total = 0

    if config.groups:
        # Grouped mode: --results-dir is the pipeline OUTPUT ROOT containing
        # one <group>/ subtree per group. Visualize each group into its own dir.
        root = Path(args.results_dir)
        for group in config.get_groups():
            unit = ExperimentConfig(
                name=f"{config.name}_{group.name}",
                control=group.control,
                conditions=group.conditions,
                cell_types=config.get_group_cell_types(group),
                geometry=config.geometry,
            )
            print("\n" + "#" * 80)
            print(f"GROUP: {group.name}")
            print("#" * 80)
            total += visualize_unit(
                unit, root / group.name / 'results', base_out / group.name, args.fdr)

        # Cross-group: compare the same cell type across groups (alone vs co-culture)
        print("\n" + "#" * 80)
        print("CROSS-GROUP COMPARISONS (same cell type across groups)")
        print("#" * 80)
        total += cross_group_comparisons(config, root, base_out, args.fdr)
    else:
        # Flat mode: --results-dir is the results directory itself (legacy layout).
        total += visualize_unit(config, Path(args.results_dir), base_out, args.fdr)

    print("\n" + "=" * 80)
    print("✓ VISUALIZATION COMPLETE!")
    print("=" * 80)
    print(f"Total plots generated: {total}")
    print(f"Output directory: {base_out}")


if __name__ == '__main__':
    main()
