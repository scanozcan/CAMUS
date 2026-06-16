"""Heatmaps and dotplots summarizing logFC / significance across cell types
and conditions."""

import numpy as np
import pandas as pd

from ._base import BasePlotter, plt

try:
    import seaborn as sns
    _HAVE_SNS = True
except Exception:
    _HAVE_SNS = False


class HeatmapPlotter(BasePlotter):

    def _top_genes(self, frames, top_n):
        """Pick top_n genes by best (min) FDR across the provided result frames."""
        best = {}
        for df in frames:
            if df is None or "FDR" not in df.columns:
                continue
            for gene, fdr in zip(df["Gene"], df["FDR"]):
                best[gene] = min(best.get(gene, np.inf), float(fdr))
        if not best:
            return []
        ordered = sorted(best.items(), key=lambda kv: kv[1])
        return [g for g, _ in ordered[:top_n]]

    def plot_celltype_heatmap(self, results_by_celltype, comparison_name, top_n=40):
        """logFC heatmap (top genes x cell types) for one comparison."""
        frames = {ct: comp.get(comparison_name)
                  for ct, comp in results_by_celltype.items()
                  if comp.get(comparison_name) is not None}
        if not frames:
            return None

        genes = self._top_genes(list(frames.values()), top_n)
        if not genes:
            return None

        mat = pd.DataFrame(index=genes)
        for ct, df in frames.items():
            lfc = df.set_index("Gene")["logFC"] if "logFC" in df.columns else pd.Series(dtype=float)
            mat[ct] = [lfc.get(g, np.nan) for g in genes]

        fig, ax = plt.subplots(figsize=(max(5, 1.2 * mat.shape[1]), max(6, 0.25 * len(genes))))
        if _HAVE_SNS:
            sns.heatmap(mat, cmap="RdBu_r", center=0, ax=ax,
                        cbar_kws={"label": "logFC"}, linewidths=0.3, linecolor="white")
        else:
            im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r")
            ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90)
            ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
            fig.colorbar(im, ax=ax, label="logFC")
        ax.set_title(f"{comparison_name}\nTop {len(genes)} genes - logFC by cell type", fontweight="bold")
        return self._save(fig, f"{comparison_name}_celltype_heatmap.png")

    def plot_clustered_heatmap(self, all_comparison_results, min_conditions=2):
        """Clustered logFC heatmap across all comparison x cell-type columns,
        restricted to genes significant in >= min_conditions columns."""
        columns = list(all_comparison_results.keys())
        # significance count per gene
        sig_count = {}
        lfc_by_col = {}
        for col, df in all_comparison_results.items():
            if df is None or "Gene" not in df.columns:
                continue
            lfc_by_col[col] = df.set_index("Gene")["logFC"] if "logFC" in df.columns else pd.Series(dtype=float)
            if "FDR" in df.columns:
                for gene, fdr in zip(df["Gene"], df["FDR"]):
                    if float(fdr) <= self.fdr_threshold:
                        sig_count[gene] = sig_count.get(gene, 0) + 1

        genes = [g for g, c in sig_count.items() if c >= min_conditions]
        if len(genes) < 2:
            return None

        mat = pd.DataFrame(index=genes)
        for col in columns:
            s = lfc_by_col.get(col, pd.Series(dtype=float))
            mat[col] = [s.get(g, np.nan) for g in genes]
        mat = mat.fillna(0.0)

        if _HAVE_SNS:
            try:
                g = sns.clustermap(mat, cmap="RdBu_r", center=0, figsize=(
                    max(6, 1.0 * mat.shape[1]), max(6, 0.25 * len(genes))),
                    cbar_kws={"label": "logFC"})
                path = self.output_dir / "clustered_heatmap.png"
                g.savefig(path, dpi=200, bbox_inches="tight")
                plt.close(g.fig)
                return path
            except Exception:
                pass
        # Fallback unclustered
        fig, ax = plt.subplots(figsize=(max(6, 1.0 * mat.shape[1]), max(6, 0.25 * len(genes))))
        im = ax.imshow(mat.values, aspect="auto", cmap="RdBu_r")
        ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=90, fontsize=7)
        ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
        fig.colorbar(im, ax=ax, label="logFC")
        ax.set_title(f"Genes significant in >= {min_conditions} comparisons", fontweight="bold")
        return self._save(fig, "clustered_heatmap.png")

    def plot_dotplot(self, results_by_celltype, top_n=50):
        """Dotplot: color = logFC, size = -log10(FDR), across cell types for
        the union of top genes."""
        # Flatten to {(gene, cell_type): (logFC, FDR)}
        frames = {}
        for ct, comp in results_by_celltype.items():
            # merge all comparisons for this cell type by taking best FDR per gene
            merged = None
            for df in comp.values():
                if df is None:
                    continue
                keep = df[["Gene", "logFC", "FDR"]].copy() if {"Gene", "logFC", "FDR"}.issubset(df.columns) else None
                if keep is None:
                    continue
                merged = keep if merged is None else pd.concat([merged, keep])
            if merged is not None and not merged.empty:
                merged = merged.sort_values("FDR").drop_duplicates("Gene")
                frames[ct] = merged.set_index("Gene")

        if not frames:
            return None

        genes = self._top_genes([f.reset_index() for f in frames.values()], top_n)
        if not genes:
            return None

        cts = list(frames.keys())
        xs, ys, colors, sizes = [], [], [], []
        for xi, ct in enumerate(cts):
            f = frames[ct]
            for yi, g in enumerate(genes):
                if g in f.index:
                    xs.append(xi); ys.append(yi)
                    colors.append(float(f.loc[g, "logFC"]))
                    sizes.append(-np.log10(max(float(f.loc[g, "FDR"]), 1e-300)) * 20 + 5)

        fig, ax = plt.subplots(figsize=(max(5, 1.4 * len(cts)), max(6, 0.25 * len(genes))))
        sc = ax.scatter(xs, ys, c=colors, s=sizes, cmap="RdBu_r", vmin=-2, vmax=2,
                        edgecolors="black", linewidths=0.3)
        ax.set_xticks(range(len(cts))); ax.set_xticklabels(cts, rotation=30)
        ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
        fig.colorbar(sc, ax=ax, label="logFC")
        ax.set_title("Top genes - logFC (color) & significance (size)", fontweight="bold")
        ax.grid(True, alpha=0.3, ls="--")
        return self._save(fig, "dotplot_logFC_FDR.png")
