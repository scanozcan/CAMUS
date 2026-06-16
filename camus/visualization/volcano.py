"""Volcano plots for MAGeCK RRA gene-level results."""

import numpy as np

from ._base import BasePlotter, plt, DEPLETED, ENRICHED, NEUTRAL


class VolcanoPlotter(BasePlotter):
    """Volcano plots of logFC vs -log10(FDR) for RRA gene results."""

    def _classify(self, df):
        df = df.copy()
        if "logFC" not in df.columns:
            df["logFC"] = 0.0
        df["neg_log10_fdr"] = -np.log10(df["FDR"].clip(lower=1e-300))
        df["call"] = "Not Significant"
        sig = df["FDR"] <= self.fdr_threshold
        df.loc[sig & (df["logFC"] < 0), "call"] = "Depleted"
        df.loc[sig & (df["logFC"] > 0), "call"] = "Enriched"
        return df

    def plot_volcano(self, gene_results, cell_type, comparison_name):
        """Single volcano for one cell type / comparison."""
        df = self._classify(gene_results)
        colors = {"Depleted": DEPLETED, "Enriched": ENRICHED, "Not Significant": NEUTRAL}

        fig, ax = plt.subplots(figsize=(8, 7))
        for call in ["Not Significant", "Depleted", "Enriched"]:
            sub = df[df["call"] == call]
            ax.scatter(sub["logFC"], sub["neg_log10_fdr"], s=18,
                       c=colors[call], alpha=0.6 if call == "Not Significant" else 0.85,
                       edgecolors="none", label=f"{call} ({len(sub)})")

        ax.axhline(-np.log10(self.fdr_threshold), color="black", ls="--", lw=0.8, alpha=0.5)
        ax.axvline(0, color="black", ls="-", lw=0.8, alpha=0.3)
        ax.set_xlabel("log2 Fold Change", fontweight="bold")
        ax.set_ylabel("-log10(FDR)", fontweight="bold")
        ax.set_title(f"{cell_type} - {comparison_name}\nRRA Volcano", fontweight="bold")
        ax.legend(loc="upper right", frameon=True)
        ax.grid(True, alpha=0.3, ls="--")
        return self._save(fig, f"{comparison_name}_{cell_type}_volcano.png")

    def plot_multi_panel_volcano(self, all_results, comparison_name, cell_type_names):
        """One panel per cell type for a given comparison."""
        results_list = all_results.get(comparison_name, [])
        by_ct = {r["cell_type"]: r["gene_results"] for r in results_list}
        cts = [c for c in cell_type_names if c in by_ct]
        if not cts:
            return None

        n = len(cts)
        fig, axes = plt.subplots(1, n, figsize=(6 * n, 6), squeeze=False)
        colors = {"Depleted": DEPLETED, "Enriched": ENRICHED, "Not Significant": NEUTRAL}

        for ax, ct in zip(axes[0], cts):
            df = self._classify(by_ct[ct])
            for call in ["Not Significant", "Depleted", "Enriched"]:
                sub = df[df["call"] == call]
                ax.scatter(sub["logFC"], sub["neg_log10_fdr"], s=14,
                           c=colors[call], alpha=0.6 if call == "Not Significant" else 0.85,
                           edgecolors="none")
            n_dep = int((df["call"] == "Depleted").sum())
            n_enr = int((df["call"] == "Enriched").sum())
            ax.axhline(-np.log10(self.fdr_threshold), color="black", ls="--", lw=0.8, alpha=0.5)
            ax.axvline(0, color="black", ls="-", lw=0.8, alpha=0.3)
            ax.set_xlabel("log2 FC")
            ax.set_ylabel("-log10(FDR)")
            ax.set_title(f"{ct}\n{n_dep} depleted, {n_enr} enriched", fontweight="bold")
            ax.grid(True, alpha=0.3, ls="--")

        fig.suptitle(f"{comparison_name} - Volcano by Cell Type", fontweight="bold")
        return self._save(fig, f"{comparison_name}_all_celltypes_volcano.png")
