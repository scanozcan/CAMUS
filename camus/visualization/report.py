"""Text reports summarizing RRA hits and count-matrix QC."""

import numpy as np

from ._base import BasePlotter


class ReportGenerator(BasePlotter):

    def generate_summary_report(self, all_results, config):
        """Write ANALYSIS_SUMMARY.txt with hit counts per comparison/cell type."""
        path = self.output_dir / "ANALYSIS_SUMMARY.txt"
        with open(path, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("CAMUS - RRA ANALYSIS SUMMARY\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Experiment: {getattr(config, 'name', 'NA')}\n")
            f.write(f"FDR threshold: {self.fdr_threshold}\n")
            f.write(f"Comparisons: {len(all_results)}\n\n")

            for comparison_name, results_list in all_results.items():
                f.write("-" * 80 + "\n")
                f.write(f"{comparison_name}\n")
                f.write("-" * 80 + "\n")
                for r in results_list:
                    df = r["gene_results"]
                    ct = r["cell_type"]
                    sig = df[df["FDR"] <= self.fdr_threshold] if "FDR" in df.columns else df.iloc[0:0]
                    dep = (sig["Direction"].str.lower() == "depleted").sum() if "Direction" in sig.columns else 0
                    enr = (sig["Direction"].str.lower() == "enriched").sum() if "Direction" in sig.columns else 0
                    f.write(f"  {ct:<16} total genes: {len(df):>6} | "
                            f"significant: {len(sig):>5} (depleted {dep}, enriched {enr})\n")
                    if "Gene" in sig.columns and len(sig) > 0:
                        top = sig.sort_values("FDR").head(10)["Gene"].tolist()
                        f.write(f"      top hits: {', '.join(map(str, top))}\n")
                f.write("\n")

            f.write("=" * 80 + "\n")
            f.write("NOTE: each cell type x comparison is an independent test; consider\n")
            f.write("multiple-testing burden across the full grid when interpreting hits.\n")
            f.write("=" * 80 + "\n")
        return path

    def generate_qc_report(self, count_matrices, cell_type_names):
        """Write QC_REPORT.txt with per-cell-type coverage statistics."""
        path = self.output_dir / "QC_REPORT.txt"
        with open(path, "w") as f:
            f.write("=" * 80 + "\n")
            f.write("CAMUS - COUNT MATRIX QC REPORT\n")
            f.write("=" * 80 + "\n\n")

            for ct in cell_type_names:
                df = count_matrices.get(ct)
                if df is None:
                    continue
                sample_cols = [c for c in df.columns if c not in ("sgRNA", "Gene", "in_library")]
                f.write(f"[{ct}]\n")
                f.write(f"  sgRNAs: {len(df):,}\n")
                f.write(f"  samples: {len(sample_cols)}\n")
                for c in sample_cols:
                    vals = df[c].astype(float).values
                    total = vals.sum()
                    zeros = int((vals == 0).sum())
                    gini = _gini(vals)
                    f.write(f"    {c:<22} reads={total:>12,.0f}  "
                            f"zero-count sgRNAs={zeros:>5}  Gini={gini:.3f}\n")
                f.write("\n")

            f.write("Gini coefficient of 0 = perfectly even library representation;\n")
            f.write("higher values indicate skew (a few sgRNAs dominate).\n")
        return path


def _gini(x):
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x)
    if n == 0 or x.sum() == 0:
        return 0.0
    cum = np.cumsum(x)
    return float((n + 1 - 2 * (cum / cum[-1]).sum()) / n)
