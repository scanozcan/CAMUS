#!/usr/bin/env python3
"""
MLE-specific post-analysis visualization for CAMUS.

Generates publication-quality plots for MAGeCK MLE results:
- Volcano plots (beta vs -log10(p-value))
- Heatmaps (beta scores across conditions)
- Dotplots (condition-specific effects)
- Cross-condition comparisons
- QC plots (reused from RRA)

Key difference from RRA: Uses BETA SCORES instead of logFC.
"""

import sys
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Dict, List, Tuple
from matplotlib.gridspec import GridSpec

# Add parent to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from camus.config.experiment import ExperimentConfig

try:
    from camus.visualization import QCPlotter
except ImportError:
    print("Warning: Cannot import QC plotting module")
    QCPlotter = None


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='CAMUS - MLE Post-Analysis Visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        '--results-dir',
        type=str,
        required=True,
        help='Directory containing MLE analysis results'
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
        default='visualizations_mle',
        help='Output directory for visualizations (default: visualizations_mle)'
    )

    parser.add_argument(
        '--fdr',
        type=float,
        default=0.05,
        help='FDR threshold (default: 0.05)'
    )

    return parser.parse_args()


def load_mle_results(results_dir: Path, config: ExperimentConfig) -> Dict[str, pd.DataFrame]:
    """
    Load MAGeCK MLE analysis results.

    Args:
        results_dir: Directory with MLE result files
        config: Experiment configuration

    Returns:
        Dictionary mapping cell type names to MLE result DataFrames
    """
    mle_results = {}

    for cell_type in config.cell_types:
        # MLE results have one file per cell type with all conditions
        gene_file = results_dir / f"{cell_type.name}_MLE_gene_results.txt"

        if not gene_file.exists():
            print(f"  Warning: Missing {gene_file.name}")
            continue

        try:
            gene_results = pd.read_csv(gene_file, sep='\t')

            # Verify MLE format (should have |beta, |p-value, |fdr columns)
            beta_cols = [col for col in gene_results.columns if col.endswith('|beta')]

            if not beta_cols:
                print(f"  Warning: {gene_file.name} doesn't have beta columns (not MLE format?)")
                continue

            mle_results[cell_type.name] = gene_results
            print(f"  Loaded {cell_type.name}: {len(gene_results)} genes, {len(beta_cols)} conditions")

        except Exception as e:
            print(f"  Error loading {gene_file.name}: {e}")

    return mle_results


def extract_condition_names(mle_df: pd.DataFrame) -> List[str]:
    """
    Extract condition names from MLE result DataFrame.

    Args:
        mle_df: MLE result DataFrame

    Returns:
        List of condition names
    """
    beta_cols = [col for col in mle_df.columns if col.endswith('|beta')]
    condition_names = [col.replace('|beta', '') for col in beta_cols]
    return condition_names


def plot_mle_volcano(
    mle_df: pd.DataFrame,
    condition: str,
    cell_type: str,
    output_path: Path,
    fdr_threshold: float = 0.05
):
    """
    Create volcano plot for MLE results (beta vs -log10(p-value)).

    Args:
        mle_df: MLE result DataFrame
        condition: Condition name
        cell_type: Cell type name
        output_path: Path to save plot
        fdr_threshold: FDR threshold for significance
    """
    # Get columns for this condition
    beta_col = f"{condition}|beta"
    pval_col = f"{condition}|p-value"
    fdr_col = f"{condition}|fdr"
    sig_col = f"{condition}|significant"

    if beta_col not in mle_df.columns or pval_col not in mle_df.columns:
        print(f"    Warning: Missing columns for {condition}")
        return

    # Prepare data
    plot_df = mle_df[['Gene', beta_col, pval_col, fdr_col]].copy()
    plot_df.columns = ['Gene', 'beta', 'p_value', 'FDR']

    # Add -log10(p-value)
    plot_df['neg_log10_pval'] = -np.log10(plot_df['p_value'].replace(0, 1e-300))

    # Classify significance using DUAL thresholds:
    # - FDR < 0.05: High confidence (dark colors)
    # - p-value < 0.05 (but FDR > 0.05): Suggestive (light colors)
    # - Neither: Not significant (gray)

    plot_df['FDR_sig'] = plot_df['FDR'] <= fdr_threshold
    plot_df['Pval_sig'] = plot_df['p_value'] <= fdr_threshold  # Use same threshold for p-value

    # Classify direction with confidence levels
    plot_df['Direction'] = 'Not Significant'

    # High confidence (FDR < 0.05)
    plot_df.loc[(plot_df['FDR_sig']) & (plot_df['beta'] < 0), 'Direction'] = 'Depleted (FDR<0.05)'
    plot_df.loc[(plot_df['FDR_sig']) & (plot_df['beta'] > 0), 'Direction'] = 'Enriched (FDR<0.05)'

    # Suggestive (p-value < 0.05 but FDR > 0.05)
    plot_df.loc[(~plot_df['FDR_sig']) & (plot_df['Pval_sig']) & (plot_df['beta'] < 0), 'Direction'] = 'Depleted (p<0.05)'
    plot_df.loc[(~plot_df['FDR_sig']) & (plot_df['Pval_sig']) & (plot_df['beta'] > 0), 'Direction'] = 'Enriched (p<0.05)'

    # Count significant genes
    n_depleted_fdr = (plot_df['Direction'] == 'Depleted (FDR<0.05)').sum()
    n_enriched_fdr = (plot_df['Direction'] == 'Enriched (FDR<0.05)').sum()
    n_depleted_p = (plot_df['Direction'] == 'Depleted (p<0.05)').sum()
    n_enriched_p = (plot_df['Direction'] == 'Enriched (p<0.05)').sum()
    n_total_fdr = n_depleted_fdr + n_enriched_fdr
    n_total_p = n_depleted_p + n_enriched_p

    # Create plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Color map with dual threshold system
    colors = {
        'Depleted (FDR<0.05)': '#2166ac',      # Dark Blue
        'Enriched (FDR<0.05)': '#b2182b',      # Dark Red
        'Depleted (p<0.05)': '#92c5de',        # Light Blue
        'Enriched (p<0.05)': '#f4a582',        # Light Red
        'Not Significant': '#cccccc'           # Gray
    }

    # Plot points (plot in order: not sig first, then suggestive, then FDR-sig on top)
    plot_order = ['Not Significant', 'Depleted (p<0.05)', 'Enriched (p<0.05)',
                  'Depleted (FDR<0.05)', 'Enriched (FDR<0.05)']

    for direction in plot_order:
        subset = plot_df[plot_df['Direction'] == direction]
        if len(subset) > 0:
            ax.scatter(
                subset['beta'],
                subset['neg_log10_pval'],
                c=colors[direction],
                label=f"{direction} ({len(subset)})",
                alpha=0.6 if 'Not' in direction else 0.8,
                s=20,
                edgecolors='none'
            )

    # Add threshold lines
    ax.axhline(y=-np.log10(fdr_threshold), color='black', linestyle='--', linewidth=0.8, alpha=0.5,
               label=f'p-value = {fdr_threshold}')
    ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8, alpha=0.3)

    # Labels
    ax.set_xlabel('Beta Score (Effect Size)', fontsize=12, fontweight='bold')
    ax.set_ylabel('-log10(p-value)', fontsize=12, fontweight='bold')

    # Build title with counts
    title_parts = [
        f'{cell_type} - {condition}',
        'MLE Volcano Plot',
        f'FDR<0.05: {n_total_fdr} ({n_depleted_fdr} depleted, {n_enriched_fdr} enriched)',
        f'p<0.05: {n_total_p} additional ({n_depleted_p} depleted, {n_enriched_p} enriched)'
    ]
    ax.set_title('\n'.join(title_parts), fontsize=12, fontweight='bold')

    # Legend
    ax.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_mle_multi_condition_volcano(
    mle_df: pd.DataFrame,
    conditions: List[str],
    cell_type: str,
    output_path: Path,
    fdr_threshold: float = 0.05
):
    """
    Create multi-panel volcano plot for all conditions.

    Args:
        mle_df: MLE result DataFrame
        conditions: List of condition names
        cell_type: Cell type name
        output_path: Path to save plot
        fdr_threshold: FDR threshold
    """
    n_conditions = len(conditions)
    n_cols = min(3, n_conditions)
    n_rows = (n_conditions + n_cols - 1) // n_cols

    fig = plt.figure(figsize=(7 * n_cols, 6 * n_rows))
    gs = GridSpec(n_rows, n_cols, figure=fig)

    for idx, condition in enumerate(conditions):
        row = idx // n_cols
        col = idx % n_cols
        ax = fig.add_subplot(gs[row, col])

        # Get columns for this condition
        beta_col = f"{condition}|beta"
        pval_col = f"{condition}|p-value"
        fdr_col = f"{condition}|fdr"
        sig_col = f"{condition}|significant"

        if beta_col not in mle_df.columns or pval_col not in mle_df.columns:
            continue

        # Prepare data
        plot_df = mle_df[['Gene', beta_col, pval_col, fdr_col]].copy()
        plot_df.columns = ['Gene', 'beta', 'p_value', 'FDR']
        plot_df['neg_log10_pval'] = -np.log10(plot_df['p_value'].replace(0, 1e-300))

        # Classify
        if sig_col in mle_df.columns:
            plot_df['Significant'] = mle_df[sig_col]
        else:
            plot_df['Significant'] = plot_df['FDR'].apply(
                lambda x: 'significant' if x <= fdr_threshold else 'not_significant'
            )

        plot_df['Direction'] = 'Not Significant'
        plot_df.loc[(plot_df['Significant'] == 'significant') & (plot_df['beta'] < 0), 'Direction'] = 'Depleted'
        plot_df.loc[(plot_df['Significant'] == 'significant') & (plot_df['beta'] > 0), 'Direction'] = 'Enriched'

        # Count
        n_depleted = (plot_df['Direction'] == 'Depleted').sum()
        n_enriched = (plot_df['Direction'] == 'Enriched').sum()

        # Plot
        colors = {
            'Depleted': '#2166ac',
            'Enriched': '#b2182b',
            'Not Significant': '#cccccc'
        }

        for direction in ['Not Significant', 'Depleted', 'Enriched']:
            subset = plot_df[plot_df['Direction'] == direction]
            ax.scatter(
                subset['beta'],
                subset['neg_log10_pval'],
                c=colors[direction],
                alpha=0.6,
                s=15,
                edgecolors='none'
            )

        # Threshold lines
        ax.axhline(y=-np.log10(fdr_threshold), color='black', linestyle='--', linewidth=0.8, alpha=0.5)
        ax.axvline(x=0, color='black', linestyle='-', linewidth=0.8, alpha=0.3)

        # Labels
        ax.set_xlabel('Beta Score', fontsize=10)
        ax.set_ylabel('-log10(p-value)', fontsize=10)
        ax.set_title(
            f'{condition}\n{n_depleted} depleted, {n_enriched} enriched',
            fontsize=11,
            fontweight='bold'
        )
        ax.grid(True, alpha=0.3, linestyle='--')

    fig.suptitle(
        f'{cell_type} - MLE Multi-Condition Volcano Plots',
        fontsize=14,
        fontweight='bold',
        y=0.995
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_beta_heatmap(
    mle_df: pd.DataFrame,
    conditions: List[str],
    cell_type: str,
    output_path: Path,
    fdr_threshold: float = 0.05,
    top_n: int = 50
):
    """
    Create heatmap of beta scores for top significant genes.

    Args:
        mle_df: MLE result DataFrame
        conditions: List of condition names
        cell_type: Cell type name
        output_path: Path to save plot
        fdr_threshold: FDR threshold
        top_n: Number of top genes to show
    """
    # Get significant genes in any condition
    sig_genes = set()
    for condition in conditions:
        sig_col = f"{condition}|significant"
        if sig_col in mle_df.columns:
            sig_genes.update(
                mle_df[mle_df[sig_col] == 'significant']['Gene'].values
            )

    if not sig_genes:
        print(f"    No significant genes for heatmap in {cell_type}")
        return

    # Filter to significant genes
    sig_df = mle_df[mle_df['Gene'].isin(sig_genes)].copy()

    # Extract beta scores
    beta_matrix = []
    for condition in conditions:
        beta_col = f"{condition}|beta"
        if beta_col in sig_df.columns:
            beta_matrix.append(sig_df[beta_col].values)

    if not beta_matrix:
        return

    beta_matrix = np.array(beta_matrix).T
    beta_df = pd.DataFrame(
        beta_matrix,
        index=sig_df['Gene'].values,
        columns=conditions
    )

    # Sort by mean absolute beta
    beta_df['mean_abs_beta'] = beta_df.abs().mean(axis=1)
    beta_df = beta_df.sort_values('mean_abs_beta', ascending=False)
    beta_df = beta_df.drop('mean_abs_beta', axis=1)

    # Take top N
    if len(beta_df) > top_n:
        beta_df = beta_df.iloc[:top_n]

    # Create heatmap
    fig, ax = plt.subplots(figsize=(8, max(10, len(beta_df) * 0.3)))

    sns.heatmap(
        beta_df,
        cmap='RdBu_r',
        center=0,
        vmin=-1.5,
        vmax=1.5,
        cbar_kws={'label': 'Beta Score'},
        linewidths=0.5,
        linecolor='white',
        ax=ax
    )

    ax.set_title(
        f'{cell_type} - Top {len(beta_df)} Significant Genes\nBeta Scores Across Conditions',
        fontsize=13,
        fontweight='bold',
        pad=20
    )
    ax.set_xlabel('Condition', fontsize=11, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=11, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_condition_comparison_dotplot(
    mle_df: pd.DataFrame,
    conditions: List[str],
    cell_type: str,
    output_path: Path,
    fdr_threshold: float = 0.05,
    top_n: int = 30
):
    """
    Create dotplot showing beta scores and significance across conditions.

    Args:
        mle_df: MLE result DataFrame
        conditions: List of condition names
        cell_type: Cell type name
        output_path: Path to save plot
        fdr_threshold: FDR threshold
        top_n: Number of top genes to show
    """
    # Get significant genes
    sig_genes = set()
    for condition in conditions:
        sig_col = f"{condition}|significant"
        if sig_col in mle_df.columns:
            sig_genes.update(
                mle_df[mle_df[sig_col] == 'significant']['Gene'].values
            )

    if not sig_genes:
        print(f"    No significant genes for dotplot in {cell_type}")
        return

    sig_df = mle_df[mle_df['Gene'].isin(sig_genes)].copy()

    # Sort by mean absolute beta
    beta_cols = [f"{c}|beta" for c in conditions]
    sig_df['mean_abs_beta'] = sig_df[beta_cols].abs().mean(axis=1)
    sig_df = sig_df.sort_values('mean_abs_beta', ascending=False)

    # Take top N
    if len(sig_df) > top_n:
        sig_df = sig_df.iloc[:top_n]

    # Prepare plot data
    plot_data = []
    for _, row in sig_df.iterrows():
        for condition in conditions:
            beta_col = f"{condition}|beta"
            fdr_col = f"{condition}|fdr"

            if beta_col in row and fdr_col in row:
                plot_data.append({
                    'Gene': row['Gene'],
                    'Condition': condition,
                    'Beta': row[beta_col],
                    'FDR': row[fdr_col],
                    'neg_log10_FDR': -np.log10(max(row[fdr_col], 1e-300))
                })

    if not plot_data:
        return

    plot_df = pd.DataFrame(plot_data)

    # Create dotplot
    fig, ax = plt.subplots(figsize=(8, max(8, len(sig_df) * 0.4)))

    # Create scatter plot
    scatter = ax.scatter(
        plot_df['Condition'],
        plot_df['Gene'],
        c=plot_df['Beta'],
        s=plot_df['neg_log10_FDR'] * 20,  # Size by significance
        cmap='RdBu_r',
        vmin=-1.5,
        vmax=1.5,
        alpha=0.8,
        edgecolors='black',
        linewidths=0.5
    )

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
    cbar.set_label('Beta Score', fontsize=10, fontweight='bold')

    # Labels
    ax.set_xlabel('Condition', fontsize=11, fontweight='bold')
    ax.set_ylabel('Gene', fontsize=11, fontweight='bold')
    ax.set_title(
        f'{cell_type} - Top {len(sig_df)} Significant Genes\n'
        f'Dot size = significance (-log10 FDR)',
        fontsize=13,
        fontweight='bold',
        pad=20
    )

    # Grid
    ax.grid(True, alpha=0.3, linestyle='--', axis='y')

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_cross_celltype_comparison(
    all_mle_results: Dict[str, pd.DataFrame],
    condition: str,
    output_path: Path,
    fdr_threshold: float = 0.05
):
    """
    Compare beta scores for a condition across cell types.

    Args:
        all_mle_results: Dict mapping cell types to MLE result DataFrames
        condition: Condition name to compare
        output_path: Path to save plot
        fdr_threshold: FDR threshold
    """
    cell_types = list(all_mle_results.keys())

    if len(cell_types) < 2:
        return

    # Collect beta scores for this condition
    beta_data = {}
    for cell_type, mle_df in all_mle_results.items():
        beta_col = f"{condition}|beta"
        if beta_col in mle_df.columns:
            beta_data[cell_type] = mle_df.set_index('Gene')[beta_col]

    if len(beta_data) < 2:
        return

    # Merge on gene
    merged_df = pd.DataFrame(beta_data)
    merged_df = merged_df.dropna()

    if len(merged_df) == 0:
        return

    # Create pairwise scatter plots
    n_celltypes = len(cell_types)
    fig, axes = plt.subplots(n_celltypes, n_celltypes, figsize=(12, 12))

    for i, ct1 in enumerate(cell_types):
        for j, ct2 in enumerate(cell_types):
            ax = axes[i, j] if n_celltypes > 1 else axes

            if i == j:
                # Diagonal: histogram
                if ct1 in merged_df.columns:
                    ax.hist(merged_df[ct1], bins=30, alpha=0.7, color='steelblue')
                    ax.set_ylabel('Frequency', fontsize=8)
                    ax.set_title(ct1, fontsize=9, fontweight='bold')
                    ax.grid(True, alpha=0.3)
            else:
                # Off-diagonal: scatter
                if ct1 in merged_df.columns and ct2 in merged_df.columns:
                    ax.scatter(
                        merged_df[ct2],
                        merged_df[ct1],
                        alpha=0.5,
                        s=10,
                        color='steelblue'
                    )

                    # Add diagonal line
                    lims = [
                        min(merged_df[ct2].min(), merged_df[ct1].min()),
                        max(merged_df[ct2].max(), merged_df[ct1].max())
                    ]
                    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)

                    # Correlation
                    corr = merged_df[ct2].corr(merged_df[ct1])
                    ax.text(
                        0.05, 0.95,
                        f'r = {corr:.3f}',
                        transform=ax.transAxes,
                        fontsize=8,
                        verticalalignment='top',
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.7)
                    )

                    ax.grid(True, alpha=0.3)

            # Labels
            if j == 0:
                ax.set_ylabel(ct1 if i != j else 'Frequency', fontsize=8)
            if i == n_celltypes - 1:
                ax.set_xlabel(ct2 if i != j else 'Beta Score', fontsize=8)

    fig.suptitle(
        f'{condition} - Beta Score Comparison Across Cell Types',
        fontsize=14,
        fontweight='bold',
        y=0.995
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def generate_summary_report(
    all_mle_results: Dict[str, pd.DataFrame],
    conditions: List[str],
    output_path: Path,
    fdr_threshold: float = 0.05
):
    """
    Generate text summary report of MLE results.

    Args:
        all_mle_results: Dict mapping cell types to MLE DataFrames
        conditions: List of condition names
        output_path: Path to save report
        fdr_threshold: FDR threshold
    """
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CAMUS - MLE RESULTS SUMMARY\n")
        f.write("=" * 80 + "\n\n")

        f.write(f"FDR Threshold: {fdr_threshold}\n")
        f.write(f"Cell Types: {len(all_mle_results)}\n")
        f.write(f"Conditions: {len(conditions)}\n\n")

        f.write("=" * 80 + "\n")
        f.write("SIGNIFICANT GENES PER CONDITION\n")
        f.write("=" * 80 + "\n\n")

        # Table header
        header = f"{'Cell Type':<15} | "
        for condition in conditions:
            header += f"{condition:<12} | "
        f.write(header + "\n")
        f.write("-" * len(header) + "\n")

        # Counts for each cell type
        for cell_type, mle_df in all_mle_results.items():
            row = f"{cell_type:<15} | "
            for condition in conditions:
                sig_col = f"{condition}|significant"
                if sig_col in mle_df.columns:
                    n_sig = (mle_df[sig_col] == 'significant').sum()
                    row += f"{n_sig:<12} | "
                else:
                    row += f"{'N/A':<12} | "
            f.write(row + "\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("BETA SCORE STATISTICS\n")
        f.write("=" * 80 + "\n\n")

        for cell_type, mle_df in all_mle_results.items():
            f.write(f"\n[{cell_type}]\n")
            f.write("-" * 40 + "\n")

            for condition in conditions:
                beta_col = f"{condition}|beta"
                sig_col = f"{condition}|significant"

                if beta_col in mle_df.columns and sig_col in mle_df.columns:
                    sig_df = mle_df[mle_df[sig_col] == 'significant']

                    if len(sig_df) > 0:
                        depleted = sig_df[sig_df[beta_col] < 0]
                        enriched = sig_df[sig_df[beta_col] > 0]

                        f.write(f"\n{condition}:\n")
                        f.write(f"  Total significant: {len(sig_df)}\n")
                        f.write(f"  Depleted (beta < 0): {len(depleted)}\n")
                        if len(depleted) > 0:
                            f.write(f"    Mean beta: {depleted[beta_col].mean():.4f}\n")
                            f.write(f"    Range: [{depleted[beta_col].min():.4f}, {depleted[beta_col].max():.4f}]\n")

                        f.write(f"  Enriched (beta > 0): {len(enriched)}\n")
                        if len(enriched) > 0:
                            f.write(f"    Mean beta: {enriched[beta_col].mean():.4f}\n")
                            f.write(f"    Range: [{enriched[beta_col].min():.4f}, {enriched[beta_col].max():.4f}]\n")

        f.write("\n" + "=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")


def main():
    """Main execution."""
    args = parse_args()

    print("=" * 80)
    print("CAMUS - MLE POST-ANALYSIS VISUALIZATION")
    print("=" * 80)

    # Load config
    print(f"\nLoading configuration: {args.config}")
    config = ExperimentConfig.from_yaml(args.config)

    # Setup output directory
    output_dir = Path(args.output)
    output_dir.mkdir(exist_ok=True, parents=True)

    results_dir = Path(args.results_dir)

    # Load MLE results
    print("\nLoading MAGeCK MLE analysis results...")
    all_mle_results = load_mle_results(results_dir, config)

    if not all_mle_results:
        print("ERROR: No MLE results found!")
        return

    print(f"  Loaded {len(all_mle_results)} cell types")

    # Get condition names from first result
    first_result = list(all_mle_results.values())[0]
    conditions = extract_condition_names(first_result)
    print(f"  Detected {len(conditions)} conditions: {', '.join(conditions)}")

    print("\n" + "=" * 80)
    print("GENERATING VISUALIZATIONS")
    print("=" * 80)

    plot_count = 0

    # 1. Volcano plots (per condition per cell type)
    print("\n[1] Generating MLE volcano plots...")
    for cell_type, mle_df in all_mle_results.items():
        for condition in conditions:
            output_path = output_dir / f"{cell_type}_{condition}_volcano.png"
            plot_mle_volcano(
                mle_df,
                condition,
                cell_type,
                output_path,
                args.fdr
            )
            plot_count += 1

        # Multi-condition volcano for each cell type
        output_path = output_dir / f"{cell_type}_all_conditions_volcano.png"
        plot_mle_multi_condition_volcano(
            mle_df,
            conditions,
            cell_type,
            output_path,
            args.fdr
        )
        plot_count += 1

    print(f"  Generated {plot_count} volcano plots")

    # 2. Beta score heatmaps
    print("\n[2] Generating beta score heatmaps...")
    heatmap_count = 0
    for cell_type, mle_df in all_mle_results.items():
        output_path = output_dir / f"{cell_type}_beta_heatmap.png"
        plot_beta_heatmap(
            mle_df,
            conditions,
            cell_type,
            output_path,
            args.fdr,
            top_n=50
        )
        heatmap_count += 1

    print(f"  Generated {heatmap_count} heatmaps")

    # 3. Dotplots
    print("\n[3] Generating condition comparison dotplots...")
    dotplot_count = 0
    for cell_type, mle_df in all_mle_results.items():
        output_path = output_dir / f"{cell_type}_beta_dotplot.png"
        plot_condition_comparison_dotplot(
            mle_df,
            conditions,
            cell_type,
            output_path,
            args.fdr,
            top_n=30
        )
        dotplot_count += 1

    print(f"  Generated {dotplot_count} dotplots")

    # 4. Cross-cell-type comparisons
    print("\n[4] Generating cross-cell-type comparisons...")
    cross_count = 0
    for condition in conditions:
        output_path = output_dir / f"{condition}_cross_celltype_comparison.png"
        plot_cross_celltype_comparison(
            all_mle_results,
            condition,
            output_path,
            args.fdr
        )
        cross_count += 1

    print(f"  Generated {cross_count} cross-cell-type plots")

    # 5. Summary report
    print("\n[5] Generating summary report...")
    report_path = output_dir / "mle_summary_report.txt"
    generate_summary_report(
        all_mle_results,
        conditions,
        report_path,
        args.fdr
    )
    print(f"  Report saved: {report_path.name}")

    # Summary
    total_plots = plot_count + heatmap_count + dotplot_count + cross_count
    print("\n" + "=" * 80)
    print("✓ VISUALIZATION COMPLETE!")
    print("=" * 80)
    print(f"\nGenerated {total_plots} plots + 1 report")
    print(f"Output directory: {output_dir}")

    print("\nGenerated files:")
    print(f"  Volcano plots: {plot_count}")
    print(f"  Beta heatmaps: {heatmap_count}")
    print(f"  Dotplots: {dotplot_count}")
    print(f"  Cross-cell-type: {cross_count}")
    print(f"  Summary report: 1")


if __name__ == '__main__':
    main()
