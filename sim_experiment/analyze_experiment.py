#!/usr/bin/env python3
"""
Analyze the simulated CAMUS experiment.

For each (group, cell type) build an initial-vs-D21 count matrix from the CAMUS
demux outputs, run MAGeCK RRA via the CAMUS wrapper (non-targeting controls
passed for normalization), then:

  * compare group3 (demultiplexed) vs group1 / group2 (pure) per cell type
    - logFC correlation (Pearson + Spearman)
    - hit concordance at FDR <= threshold (direction-aware)
  * score recovery against the simulation ground truth (precision/recall/F1)
  * write scatter figures and a markdown report.

Comparisons:
    prolif   : G1_prolif  vs  G3_mixed (prolif lane)
    high_cin : G2_high_cin vs  G3_mixed (high_cin lane)
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from camus.wrappers.mageck_wrapper import MAGeCKWrapper

TIMEPOINTS = ["initial", "D21"]
N_REPS = 3
# (group, cell_type, label) units to test
UNITS = [
    ("G1_prolif", "prolif", "prolif_pure"),
    ("G2_high_cin", "high_cin", "highcin_pure"),
    ("G3_mixed", "prolif", "prolif_mixed"),
    ("G3_mixed", "high_cin", "highcin_mixed"),
]
CAT_COLORS = {
    "shared_essential": "#6a3d9a", "prolif_essential": "#1f78b4",
    "highcin_essential": "#e31a1c", "prolif_enriched": "#33a02c",
    "highcin_enriched": "#ff7f00", "prolif_rescue": "#000000",
    "highcin_rescue": "#000000", "neutral": "#cccccc",
}
# which ground-truth column applies to each unit's context
CONTEXT_OF = {"prolif_pure": "alone", "highcin_pure": "alone",
              "prolif_mixed": "co", "highcin_mixed": "co"}


def build_matrix(demux_dir, group, cell_type, out_path):
    """Combine the cell_type column across the group's 6 samples into a MAGeCK
    count matrix (sgRNA, Gene, <samples>)."""
    cols = {}
    base = None
    for tp in TIMEPOINTS:
        for rep in range(1, N_REPS + 1):
            sample = f"{group}_{tp}_Rep{rep}"
            cf = os.path.join(demux_dir, f"{sample}_count.txt")
            df = pd.read_csv(cf, sep="\t")
            if base is None:
                base = df[["sgRNA", "Gene"]].copy()
            cols[f"{tp}_Rep{rep}"] = df[cell_type].values
    mat = base.copy()
    for name, vals in cols.items():
        mat[name] = vals
    mat.to_csv(out_path, sep="\t", index=False)
    return mat, list(cols.keys())


def run_mageck(wrapper, matrix_path, sample_cols, out_prefix, control_sgrna, fdr):
    control = [c for c in sample_cols if c.startswith("initial")]
    treat = [c for c in sample_cols if c.startswith("D21")]
    out = wrapper.run_test(
        count_matrix=matrix_path, control_samples=control, treatment_samples=treat,
        output_prefix=out_prefix, normalization="control" if control_sgrna else "median",
        control_sgrna=control_sgrna, fdr_threshold=fdr,
    )
    genes = wrapper.parse_gene_results(out["gene_summary"], fdr_threshold=fdr)
    return genes[["Gene", "logFC", "FDR", "Direction", "Significance", "Num_sgRNAs"]]


def compare(pure, mixed, truth, ct, fdr, out_png):
    """Compare pure vs mixed gene results for one cell type."""
    p = pure.set_index("Gene"); m = mixed.set_index("Gene")
    common = [g for g in p.index if g in m.index and g != "Non-Targeting"]
    p = p.loc[common]; m = m.loc[common]

    pear = float(np.corrcoef(p["logFC"], m["logFC"])[0, 1])
    spear = float(pd.Series(p["logFC"].values).corr(pd.Series(m["logFC"].values), method="spearman"))

    pure_hits = set(p.index[p["Significance"] == "significant"])
    mixed_hits = set(m.index[m["Significance"] == "significant"])
    inter = pure_hits & mixed_hits
    union = pure_hits | mixed_hits
    jaccard = len(inter) / len(union) if union else float("nan")
    # direction-concordant among shared hits
    dir_ok = sum(1 for g in inter if p.loc[g, "Direction"] == m.loc[g, "Direction"])

    # scatter colored by ground-truth category (rescue genes drawn large/black)
    cat = truth.set_index("Gene")["category"]
    fig, ax = plt.subplots(figsize=(7, 7))
    for c, color in CAT_COLORS.items():
        gs = [g for g in common if cat.get(g, "neutral") == c]
        if gs:
            is_rescue = "rescue" in c
            ax.scatter(p.loc[gs, "logFC"], m.loc[gs, "logFC"],
                       s=90 if is_rescue else 22, c=color,
                       marker="*" if is_rescue else "o",
                       label=f"{c} ({len(gs)})", alpha=0.9,
                       edgecolors="black" if is_rescue else "none",
                       zorder=5 if is_rescue else 2)
    lims = [min(p["logFC"].min(), m["logFC"].min()) - 0.3,
            max(p["logFC"].max(), m["logFC"].max()) + 0.3]
    ax.plot(lims, lims, "k--", lw=1, alpha=0.6)
    ax.set_xlim(lims); ax.set_ylim(lims)
    ax.set_xlabel(f"logFC pure ({ct})", fontweight="bold")
    ax.set_ylabel(f"logFC mixed/demuxed ({ct})", fontweight="bold")
    ax.set_title(f"{ct}: pure vs demultiplexed\nPearson r={pear:.3f}, Spearman={spear:.3f}, "
                 f"hit Jaccard={jaccard:.2f}", fontweight="bold")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3, ls="--")
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return {"pearson": pear, "spearman": spear, "pure_hits": len(pure_hits),
            "mixed_hits": len(mixed_hits), "shared_hits": len(inter),
            "jaccard": jaccard, "dir_concordant": dir_ok}


def recovery(genes, truth, ct, context, fdr):
    """Precision/recall/F1 of significant calls vs ground-truth movers (|lfc|>=1)
    for the appropriate context ('alone' for pure groups, 'co' for co-culture)."""
    g = genes[genes["Gene"] != "Non-Targeting"].set_index("Gene")
    t = truth.set_index("Gene")[f"{ct}_{context}_lfc"]
    true_movers = set(t.index[t.abs() >= 1.0])
    called = set(g.index[g["Significance"] == "significant"])
    tp = len(called & true_movers)
    prec = tp / len(called) if called else float("nan")
    rec = tp / len(true_movers) if true_movers else float("nan")
    f1 = (2 * prec * rec / (prec + rec)) if (prec and rec and prec + rec > 0) else float("nan")
    corr = float(np.corrcoef(g["logFC"], t.loc[g.index])[0, 1])
    return {"true_movers": len(true_movers), "called": len(called), "tp": tp,
            "precision": prec, "recall": rec, "f1": f1, "lfc_corr_truth": corr}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True)
    ap.add_argument("--demux", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mageck-binary", default=None)
    ap.add_argument("--conda-env", default=None)
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    truth = pd.read_csv(os.path.join(args.data, "ground_truth.csv"))

    # Non-targeting control list for control normalization
    control_file = os.path.join(args.out, "control_sgrnas.txt")
    lib = pd.read_csv(os.path.join(args.data, "combined_library.txt"), sep="\t",
                      header=None, names=["seq", "gene"])
    lib[lib["gene"] == "Non-Targeting"]["seq"].to_csv(control_file, index=False, header=False)

    wrapper = MAGeCKWrapper(conda_env=args.conda_env, mageck_binary=args.mageck_binary)

    results = {}
    for group, ct, label in UNITS:
        print(f"\n=== {label}: {group} / {ct} ===")
        mat_path = os.path.join(args.out, f"{label}_matrix.txt")
        mat, cols = build_matrix(args.demux, group, ct, mat_path)
        genes = run_mageck(wrapper, mat_path, cols,
                           os.path.join(args.out, f"{label}_RRA"), control_file, args.fdr)
        genes.to_csv(os.path.join(args.out, f"{label}_gene_results.txt"), sep="\t", index=False)
        results[label] = genes

    # Comparisons + recovery
    report = []
    report.append("# CAMUS simulated experiment — results\n")
    report.append("Comparison of demultiplexed group-3 results against the pure "
                  "single-cell-type screens (group 1 / group 2), and recovery vs ground truth.\n")

    rec_rows = []
    for label in results:
        ct = "prolif" if "prolif" in label else "high_cin"
        rec_rows.append((label, recovery(results[label], truth, ct, CONTEXT_OF[label], args.fdr)))

    report.append("## Recovery vs ground truth (|true logFC| >= 1)\n")
    report.append("| unit | true movers | called sig | TP | precision | recall | F1 | logFC corr(truth) |")
    report.append("|---|---|---|---|---|---|---|---|")
    for label, r in rec_rows:
        report.append(f"| {label} | {r['true_movers']} | {r['called']} | {r['tp']} | "
                      f"{r['precision']:.2f} | {r['recall']:.2f} | {r['f1']:.2f} | {r['lfc_corr_truth']:.3f} |")
    report.append("")

    report.append("## Demultiplexed (group 3) vs pure (group 1/2)\n")
    report.append("| cell type | Pearson r | Spearman | pure hits | mixed hits | shared | Jaccard | dir-concordant |")
    report.append("|---|---|---|---|---|---|---|---|")
    comps = [("prolif", "prolif_pure", "prolif_mixed"),
             ("high_cin", "highcin_pure", "highcin_mixed")]
    for ct, pure_label, mixed_label in comps:
        c = compare(results[pure_label], results[mixed_label], truth, ct, args.fdr,
                    os.path.join(args.out, f"{ct}_pure_vs_mixed.png"))
        report.append(f"| {ct} | {c['pearson']:.3f} | {c['spearman']:.3f} | {c['pure_hits']} | "
                      f"{c['mixed_hits']} | {c['shared_hits']} | {c['jaccard']:.2f} | {c['dir_concordant']} |")
    report.append("")
    report.append("Figures: `prolif_pure_vs_mixed.png`, `high_cin_pure_vs_mixed.png` "
                  "(points colored by ground-truth category; rescue genes drawn as black "
                  "stars; dashed line = identity).\n")

    # Context-dependent (rescue) genes: depleted ALONE, rescued in co-culture.
    report.append("## Context-dependent (rescue) genes\n")
    report.append("These genes are essential in pure culture but rescued when the two cell "
                  "types grow together (group 3) — a non-cell-autonomous interaction only "
                  "visible because CAMUS separates the co-culture by barcode.\n")
    report.append("| gene | cell type | pure logFC | pure sig | mixed logFC | mixed sig |")
    report.append("|---|---|---|---|---|---|")
    resc = truth[truth["category"].str.contains("rescue")]
    label_for = {"prolif": ("prolif_pure", "prolif_mixed"),
                 "high_cin": ("highcin_pure", "highcin_mixed")}
    for _, row in resc.iterrows():
        gene = row["Gene"]
        ct = "prolif" if row["category"] == "prolif_rescue" else "high_cin"
        pure_l, mixed_l = label_for[ct]
        pp = results[pure_l].set_index("Gene")
        mm = results[mixed_l].set_index("Gene")
        if gene in pp.index and gene in mm.index:
            report.append(
                f"| {gene} | {ct} | {pp.loc[gene,'logFC']:.2f} | {pp.loc[gene,'Significance']} "
                f"| {mm.loc[gene,'logFC']:.2f} | {mm.loc[gene,'Significance']} |")
    report.append("\nExpected pattern: **significant/depleted in pure, not-significant/~0 in mixed.**\n")

    with open(os.path.join(args.out, "REPORT.md"), "w") as f:
        f.write("\n".join(report))
    print("\nWrote", os.path.join(args.out, "REPORT.md"))
    print("\n".join(report))


if __name__ == "__main__":
    main()
