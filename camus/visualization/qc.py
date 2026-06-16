"""QC plots: coverage, count distributions, replicate correlations, and
cell-type composition from the demultiplexer proportions."""

import numpy as np
import pandas as pd

from ._base import BasePlotter, plt


class QCPlotter(BasePlotter):
    def __init__(self, output_dir, fdr_threshold: float = 0.05):
        super().__init__(output_dir, fdr_threshold)

    # ---- count-level QC ----------------------------------------------------

    def plot_coverage_distribution(self, counts_df, sample_groups):
        """Histogram of per-sgRNA total counts (log10) across samples."""
        fig, ax = plt.subplots(figsize=(8, 6))
        totals = counts_df.sum(axis=1).replace(0, np.nan).dropna()
        ax.hist(np.log10(totals + 1), bins=50, color="steelblue", alpha=0.8)
        ax.set_xlabel("log10(total counts per sgRNA + 1)", fontweight="bold")
        ax.set_ylabel("Number of sgRNAs", fontweight="bold")
        ax.set_title("Coverage Distribution", fontweight="bold")
        ax.grid(True, alpha=0.3, ls="--")
        return self._save(fig, "coverage_distribution.png")

    def plot_count_distributions(self, counts_df, cell_type):
        """Boxplot of log10 counts per sample for a cell type."""
        fig, ax = plt.subplots(figsize=(max(6, 1.2 * counts_df.shape[1]), 6))
        data = [np.log10(counts_df[c].values + 1) for c in counts_df.columns]
        ax.boxplot(data, labels=list(counts_df.columns), showfliers=False)
        ax.set_ylabel("log10(count + 1)", fontweight="bold")
        ax.set_title(f"{cell_type} - Count Distributions", fontweight="bold")
        ax.tick_params(axis="x", rotation=90)
        ax.grid(True, alpha=0.3, ls="--", axis="y")
        return self._save(fig, f"{cell_type}_count_distributions.png")

    def plot_replicate_correlations(self, counts_df, sample_groups, title="Replicate Correlations"):
        """Heatmap of Spearman correlations between samples."""
        log_df = np.log10(counts_df.astype(float) + 1)
        corr = log_df.corr(method="spearman")
        fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(corr)), max(5, 0.7 * len(corr))))
        im = ax.imshow(corr.values, vmin=0, vmax=1, cmap="viridis")
        ax.set_xticks(range(len(corr)))
        ax.set_yticks(range(len(corr)))
        ax.set_xticklabels(corr.columns, rotation=90, fontsize=8)
        ax.set_yticklabels(corr.index, fontsize=8)
        for i in range(len(corr)):
            for j in range(len(corr)):
                ax.text(j, i, f"{corr.values[i, j]:.2f}", ha="center", va="center",
                        color="white" if corr.values[i, j] < 0.7 else "black", fontsize=6)
        fig.colorbar(im, ax=ax, label="Spearman r")
        ax.set_title(title, fontweight="bold")
        return self._save(fig, "replicate_correlations.png")

    # ---- cell-type composition QC -----------------------------------------

    @staticmethod
    def _composition_matrix(celltype_data, cell_type_names):
        """Build a samples x cell_types matrix of read counts from demux
        proportion tables ({sample: df[CellType, Read_Count, ...]})."""
        rows = {}
        for sample, df in celltype_data.items():
            d = dict(zip(df["CellType"], df["Read_Count"]))
            rows[sample] = [float(d.get(ct, 0)) for ct in cell_type_names]
        mat = pd.DataFrame(rows, index=cell_type_names).T
        return mat.sort_index()

    def plot_celltype_proportions(self, celltype_data, cell_type_names):
        """Stacked bar of cell-type proportions per sample."""
        mat = self._composition_matrix(celltype_data, cell_type_names)
        props = mat.div(mat.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)

        fig, ax = plt.subplots(figsize=(max(7, 1.0 * len(mat)), 6))
        bottom = np.zeros(len(props))
        for ct in cell_type_names:
            ax.bar(props.index, props[ct].values, bottom=bottom, label=ct)
            bottom += props[ct].values
        ax.set_ylabel("Proportion of assigned reads", fontweight="bold")
        ax.set_title("Cell-Type Proportions per Sample", fontweight="bold")
        ax.tick_params(axis="x", rotation=90)
        ax.legend(loc="upper right", fontsize=8)
        return self._save(fig, "celltype_proportions.png")

    def plot_celltype_pie_charts(self, celltype_data, cell_type_names):
        """One pie per sample of cell-type composition."""
        mat = self._composition_matrix(celltype_data, cell_type_names)
        n = len(mat)
        ncol = min(4, n) or 1
        nrow = (n + ncol - 1) // ncol
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 4 * nrow), squeeze=False)
        for idx, sample in enumerate(mat.index):
            ax = axes[idx // ncol][idx % ncol]
            vals = mat.loc[sample].values
            if vals.sum() > 0:
                ax.pie(vals, labels=cell_type_names, autopct="%1.0f%%", textprops={"fontsize": 7})
            ax.set_title(sample, fontsize=9, fontweight="bold")
        for idx in range(n, nrow * ncol):
            axes[idx // ncol][idx % ncol].axis("off")
        fig.suptitle("Cell-Type Composition per Sample", fontweight="bold")
        return self._save(fig, "celltype_pie_charts.png")

    def plot_celltype_replicate_aware(self, celltype_data, cell_type_names):
        """Grouped bar of per-cell-type read counts, grouped by condition prefix."""
        mat = self._composition_matrix(celltype_data, cell_type_names)
        fig, ax = plt.subplots(figsize=(max(8, 1.0 * len(mat)), 6))
        x = np.arange(len(mat))
        width = 0.8 / max(1, len(cell_type_names))
        for i, ct in enumerate(cell_type_names):
            ax.bar(x + i * width, mat[ct].values, width=width, label=ct)
        ax.set_xticks(x + 0.4)
        ax.set_xticklabels(mat.index, rotation=90)
        ax.set_ylabel("Assigned reads", fontweight="bold")
        ax.set_title("Cell-Type Read Counts per Sample (replicate-aware)", fontweight="bold")
        ax.legend(fontsize=8)
        return self._save(fig, "celltype_replicate_aware.png")
