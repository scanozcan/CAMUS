"""Cross-condition comparison plots: hit-overlap (Venn) and logFC concordance."""

import numpy as np

from ._base import BasePlotter, plt

try:
    from matplotlib_venn import venn2, venn3
    _HAVE_VENN = True
except Exception:  # matplotlib-venn is optional
    _HAVE_VENN = False


class ComparisonPlotter(BasePlotter):

    def plot_venn_diagram(self, hits_per_condition, cell_type):
        """Venn (2-3 sets) or fallback overlap bar of significant-gene sets.

        hits_per_condition: {condition_name: set(genes)}
        """
        conditions = list(hits_per_condition.keys())
        sets = [set(hits_per_condition[c]) for c in conditions]
        fig, ax = plt.subplots(figsize=(7, 6))

        if _HAVE_VENN and len(sets) == 2:
            venn2(sets, set_labels=conditions, ax=ax)
        elif _HAVE_VENN and len(sets) == 3:
            venn3(sets, set_labels=conditions, ax=ax)
        else:
            # Fallback: bar chart of set sizes + shared-in-all count.
            sizes = [len(s) for s in sets]
            shared = set.intersection(*sets) if sets else set()
            labels = list(conditions) + ["Shared(all)"]
            values = sizes + [len(shared)]
            ax.bar(labels, values, color="steelblue")
            ax.set_ylabel("Significant genes")
            ax.tick_params(axis="x", rotation=30)

        ax.set_title(f"{cell_type} - Hit Overlap Across Conditions", fontweight="bold")
        return self._save(fig, f"{cell_type}_venn_diagram.png")

    def plot_logfc_comparison(self, results1, results2, cond1, cond2, cell_type,
                              label_genes=None):
        """Scatter of per-gene logFC for two conditions in the same cell type.

        ``label_genes`` (optional): an iterable of gene names to annotate on the
        plot. Used to mark off-diagonal (context-dependent) genes; points for
        labelled genes are also drawn larger/red so they stand out.
        """
        m = results1[["Gene", "logFC"]].merge(
            results2[["Gene", "logFC"]], on="Gene", suffixes=(f"_{cond1}", f"_{cond2}")
        )
        if m.empty:
            return None

        x = m[f"logFC_{cond2}"].values
        y = m[f"logFC_{cond1}"].values

        fig, ax = plt.subplots(figsize=(7, 7))
        ax.scatter(x, y, s=12, alpha=0.5, color="steelblue", edgecolors="none")
        lims = [np.nanmin([x.min(), y.min()]), np.nanmax([x.max(), y.max()])]
        ax.plot(lims, lims, "k--", lw=1, alpha=0.6)
        if len(m) > 2:
            r = np.corrcoef(x, y)[0, 1]
            ax.text(0.05, 0.95, f"r = {r:.3f}", transform=ax.transAxes,
                    va="top", bbox=dict(boxstyle="round", fc="white", alpha=0.7))

        # Highlight + annotate requested genes (e.g. context-dependent / rescue).
        if label_genes:
            label_set = set(label_genes)
            hl = m[m["Gene"].isin(label_set)]
            if not hl.empty:
                hx = hl[f"logFC_{cond2}"].values
                hy = hl[f"logFC_{cond1}"].values
                ax.scatter(hx, hy, s=45, color="crimson", edgecolors="black",
                           linewidths=0.5, zorder=5)
                # Stagger label offsets (with connector lines) so nearby labels
                # don't overlap when several genes cluster together.
                offsets = [(8, 8), (8, -12), (-40, 8), (-40, -12)]
                for i, (gene, gx, gy) in enumerate(zip(hl["Gene"].values, hx, hy)):
                    dx, dy = offsets[i % len(offsets)]
                    ax.annotate(gene, (gx, gy), xytext=(dx, dy),
                                textcoords="offset points", fontsize=8,
                                fontweight="bold", color="crimson", zorder=6,
                                arrowprops=dict(arrowstyle="-", color="crimson",
                                                lw=0.5, alpha=0.6))

        ax.set_xlabel(f"logFC ({cond2})", fontweight="bold")
        ax.set_ylabel(f"logFC ({cond1})", fontweight="bold")
        ax.set_title(f"{cell_type} - {cond1} vs {cond2} logFC", fontweight="bold")
        ax.grid(True, alpha=0.3, ls="--")
        return self._save(fig, f"{cell_type}_{cond1}_vs_{cond2}_logfc.png")
