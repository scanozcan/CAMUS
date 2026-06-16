#!/usr/bin/env python3
"""
Simulate a two-cell-type, multiplexed CRISPR screen for CAMUS.

Cell types : prolif (proliferative), high_cin (high chromosomal instability)
Library    : 100 secretome / surface-receptor genes x 6 sgRNAs + 500 non-targeting
             controls = 1100 sgRNAs. Both cell types share the SAME library but
             carry different CAMUS cell-type barcodes.
Groups     : G1 = prolif only, G2 = high_cin only, G3 = both pooled together.
Timepoints : initial, D21      Replicates: 3
=> 3 groups x 2 timepoints x 3 reps = 18 samples.

Ground truth: each cell type has its own depleted (essential) and enriched genes,
plus a shared-essential set; non-targeting controls and neutral genes do not move.
At D21, abundances are scaled by 2**logFC (per cell type); the "initial" timepoint
reflects the input library. In G3 the read budget is split 50/50 between cell types.

Reads are written as paired FASTQ.gz with CAMUS's default geometry:
    R1 = 12 bp random prefix + 20 bp sgRNA            (>= grna_start 12 + 20)
    R2 = 22 bp random prefix + 8 bp cell-type barcode (>= barcode_start 22 + 8)
A small per-base substitution rate is injected into the sgRNA and barcode windows
to exercise CAMUS's 1-mismatch recovery and barcode error-correction.

This writes generation across samples in parallel (one worker per sample).
"""

import argparse
import csv
import gzip
import os
import random
from multiprocessing import Pool

import numpy as np

# ----------------------------------------------------------------------------
# Fixed design constants
# ----------------------------------------------------------------------------
CELL_TYPES = ["prolif", "high_cin"]
BARCODES = {"prolif": "ATGCAGGG", "high_cin": "GTTGCAGC"}  # Hamming distance 6

GROUPS = {
    "G1_prolif": ["prolif"],
    "G2_high_cin": ["high_cin"],
    "G3_mixed": ["prolif", "high_cin"],
}
TIMEPOINTS = ["initial", "D21"]
N_REPS = 3

SGRNAS_PER_GENE = 6
N_NTC = 500
GRNA_LEN = 20
R1_PREFIX = 12          # gRNA starts at 12
R2_PREFIX = 22          # barcode starts at 22
ERR_RATE_GRNA = 0.04    # fraction of reads with a 1bp sgRNA-window error
ERR_RATE_BC = 0.04      # fraction of reads with a 1bp barcode-window error
GLOBAL_SEED = 20260613

BASES = np.frombuffer(b"ACGT", dtype=np.uint8)

# 100 secretome + cell-surface-receptor gene symbols
GENES = [
    # secreted factors / cytokines / growth factors / ECM (curated)
    "VEGFA", "VEGFC", "EGF", "FGF1", "FGF2", "PDGFA", "PDGFB", "TGFB1", "TGFB2",
    "TGFB3", "IGF1", "IGF2", "HGF", "BMP2", "BMP4", "BMP7", "WNT3A", "WNT5A",
    "IL1B", "IL6", "CXCL8", "IL10", "IL11", "IL15", "IL18", "TNF", "LIF",
    "CSF1", "CSF2", "CSF3", "CXCL1", "CXCL2", "CXCL10", "CXCL12", "CCL2",
    "CCL5", "CCL20", "ANGPT1", "ANGPT2", "ANGPTL4", "SPP1", "MMP1", "MMP2",
    "MMP9", "MMP14", "TIMP1", "TIMP2", "LOX", "LOXL2", "SERPINE1", "THBS1",
    "FN1", "POSTN", "TNC", "CCN2", "CCN1", "DKK1", "SFRP1", "GREM1", "NOG",
    "INHBA", "FST", "WISP1", "AREG", "EREG", "HBEGF", "TGFA", "NRG1", "SHH",
    "DLL4", "JAG1", "GDF15",
    # cell-surface receptors
    "EGFR", "ERBB2", "ERBB3", "MET", "FGFR1", "FGFR2", "IGF1R", "PDGFRA",
    "PDGFRB", "KDR", "FLT1", "NOTCH1", "NOTCH2", "FZD7", "LRP6", "TGFBR1",
    "TGFBR2", "BMPR1A", "ACVR1", "ITGAV", "ITGB1", "ITGA5", "CD44", "CXCR4",
    "IL6R", "TNFRSF1A", "EPHA2", "EPHB2",
]
assert len(GENES) == 100, f"expected 100 genes, got {len(GENES)}"


# ----------------------------------------------------------------------------
# Library + ground truth
# ----------------------------------------------------------------------------
def rand_seq(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def build_library_and_truth(seed=GLOBAL_SEED):
    """Return (sgrnas, truth) where sgrnas is a list of (seq, gene, sgid) and
    truth maps gene -> {category, prolif_lfc, high_cin_lfc}."""
    rng = random.Random(seed)
    nprng = np.random.default_rng(seed)

    # Assign categories to the 100 genes
    genes = GENES[:]
    rng.shuffle(genes)
    cats = {}
    i = 0
    def take(n):
        nonlocal i
        chunk = genes[i:i + n]
        i += n
        return chunk
    for g in take(10):
        cats[g] = "shared_essential"
    for g in take(12):
        cats[g] = "prolif_essential"
    for g in take(12):
        cats[g] = "highcin_essential"
    for g in take(5):
        cats[g] = "prolif_enriched"
    for g in take(5):
        cats[g] = "highcin_enriched"
    for g in take(2):
        cats[g] = "prolif_rescue"     # depleted in prolif ALONE, rescued in co-culture
    for g in take(2):
        cats[g] = "highcin_rescue"    # depleted in high_cin ALONE, rescued in co-culture
    for g in genes[i:]:
        cats[g] = "neutral"

    def dep():
        return float(nprng.normal(-3.0, 0.4))
    def enr():
        return float(nprng.normal(2.0, 0.3))

    # Each gene gets an effect for the "alone" context (pure culture, groups 1/2)
    # and the "co" context (co-culture, group 3). They differ only for rescue genes:
    # a rescue knockout is lethal alone but neutral in co-culture, because the other
    # cell type supplies the missing secreted factor in trans.
    truth = {}
    for g in GENES:
        c = cats[g]
        pa = pc = ha = hc = 0.0  # prolif_alone, prolif_co, highcin_alone, highcin_co
        if c == "shared_essential":
            d = dep(); pa = pc = d; ha = hc = d
        elif c == "prolif_essential":
            d = dep(); pa = pc = d
        elif c == "highcin_essential":
            d = dep(); ha = hc = d
        elif c == "prolif_enriched":
            e = enr(); pa = pc = e
        elif c == "highcin_enriched":
            e = enr(); ha = hc = e
        elif c == "prolif_rescue":
            pa = dep(); pc = 0.0          # depleted alone, rescued together
        elif c == "highcin_rescue":
            ha = dep(); hc = 0.0
        truth[g] = {"category": c,
                    "prolif_alone": pa, "prolif_co": pc,
                    "high_cin_alone": ha, "high_cin_co": hc}

    # Build sgRNAs (unique random 20-mers)
    seen = set()
    sgrnas = []  # (seq, gene, sgid)
    for g in GENES:
        for k in range(1, SGRNAS_PER_GENE + 1):
            s = rand_seq(rng, GRNA_LEN)
            while s in seen:
                s = rand_seq(rng, GRNA_LEN)
            seen.add(s)
            sgrnas.append((s, g, f"{g}_sg{k}"))
    for k in range(1, N_NTC + 1):
        s = rand_seq(rng, GRNA_LEN)
        while s in seen:
            s = rand_seq(rng, GRNA_LEN)
        seen.add(s)
        sgrnas.append((s, "Non-Targeting", f"NTC_{k:04d}"))

    return sgrnas, truth


def write_assets(data_dir, sgrnas, truth):
    os.makedirs(data_dir, exist_ok=True)
    # CAMUS library: TSV  sequence <TAB> gene   (no header)
    with open(os.path.join(data_dir, "combined_library.txt"), "w") as f:
        for seq, gene, _ in sgrnas:
            f.write(f"{seq}\t{gene}\n")
    # Barcodes CSV
    with open(os.path.join(data_dir, "celltype_barcodes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Celltype", "Barcode"])
        for ct in CELL_TYPES:
            w.writerow([ct, BARCODES[ct]])
    # sgRNA id map (sequence -> sgid, gene)  for readable reporting
    with open(os.path.join(data_dir, "sgrna_id_map.tsv"), "w") as f:
        f.write("sequence\tsgID\tGene\n")
        for seq, gene, sgid in sgrnas:
            f.write(f"{seq}\t{sgid}\t{gene}\n")
    # Ground truth (separate effects for pure "alone" vs "co"-culture context)
    with open(os.path.join(data_dir, "ground_truth.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Gene", "category",
                    "prolif_alone_lfc", "prolif_co_lfc",
                    "high_cin_alone_lfc", "high_cin_co_lfc"])
        for g in GENES:
            t = truth[g]
            w.writerow([g, t["category"],
                        f"{t['prolif_alone']:.4f}", f"{t['prolif_co']:.4f}",
                        f"{t['high_cin_alone']:.4f}", f"{t['high_cin_co']:.4f}"])


# ----------------------------------------------------------------------------
# Read generation
# ----------------------------------------------------------------------------
def base_abundance(n_sgrnas, seed):
    rng = np.random.default_rng(seed)
    return rng.lognormal(mean=0.0, sigma=0.35, size=n_sgrnas)


def sgrna_weights(base, lfc_per_sg):
    w = base * np.power(2.0, lfc_per_sg)
    return w / w.sum()


def _seq_matrix(sgrna_seqs):
    """Encode sgRNA sequences as a (n, 20) uint8 matrix of base codes 0..3."""
    lut = {65: 0, 67: 1, 71: 2, 84: 3}  # A C G T
    mat = np.zeros((len(sgrna_seqs), GRNA_LEN), dtype=np.uint8)
    for i, s in enumerate(sgrna_seqs):
        for j, ch in enumerate(s.encode()):
            mat[i, j] = lut[ch]
    return mat


def _build_records(read_idx, sg_mat, barcode_codes, rng, chunk=500_000):
    """Yield (r1_bytes, r2_bytes) blocks for the given per-read sgRNA indices."""
    r1_len = 3 + (R1_PREFIX + GRNA_LEN) + 3 + (R1_PREFIX + GRNA_LEN) + 1
    seqlen1 = R1_PREFIX + GRNA_LEN  # 32
    r2_seqlen = R2_PREFIX + len(barcode_codes)  # 30

    n = len(read_idx)
    for start in range(0, n, chunk):
        idx = read_idx[start:start + chunk]
        m = len(idx)

        # ---- R1 sequence (uint8 codes) ----
        s1 = np.empty((m, seqlen1), dtype=np.uint8)
        s1[:, :R1_PREFIX] = rng.integers(0, 4, size=(m, R1_PREFIX), dtype=np.uint8)
        s1[:, R1_PREFIX:] = sg_mat[idx]
        # inject 1bp sgRNA-window errors
        err = rng.random(m) < ERR_RATE_GRNA
        if err.any():
            rows = np.nonzero(err)[0]
            pos = rng.integers(R1_PREFIX, seqlen1, size=rows.size)
            s1[rows, pos] = rng.integers(0, 4, size=rows.size, dtype=np.uint8)

        # ---- R2 sequence ----
        s2 = np.empty((m, r2_seqlen), dtype=np.uint8)
        s2[:, :R2_PREFIX] = rng.integers(0, 4, size=(m, R2_PREFIX), dtype=np.uint8)
        s2[:, R2_PREFIX:] = barcode_codes
        err2 = rng.random(m) < ERR_RATE_BC
        if err2.any():
            rows = np.nonzero(err2)[0]
            pos = rng.integers(R2_PREFIX, r2_seqlen, size=rows.size)
            s2[rows, pos] = rng.integers(0, 4, size=rows.size, dtype=np.uint8)

        yield (_encode_fastq(s1), _encode_fastq(s2))


def _encode_fastq(seq_codes):
    """Convert an (m, L) base-code matrix to a FASTQ byte block.
    Header '@r', '+' separator, quality all 'I' (Phred 40)."""
    m, L = seq_codes.shape
    reclen = 3 + L + 3 + L + 1  # '@r\n' seq '\n+\n' qual '\n'
    rec = np.empty((m, reclen), dtype=np.uint8)
    rec[:, 0] = 64   # @
    rec[:, 1] = 114  # r
    rec[:, 2] = 10   # \n
    rec[:, 3:3 + L] = BASES[seq_codes]
    off = 3 + L
    rec[:, off] = 10      # \n
    rec[:, off + 1] = 43  # +
    rec[:, off + 2] = 10  # \n
    rec[:, off + 3:off + 3 + L] = 73  # 'I'
    rec[:, off + 3 + L] = 10  # \n
    return rec.tobytes()


def simulate_sample(task):
    """Generate one sample's R1/R2 fastq.gz. task is a dict."""
    fastq_dir = task["fastq_dir"]
    sample = task["sample"]
    cell_types = task["cell_types"]
    reads_total = task["reads_total"]
    sg_mat = task["sg_mat"]
    base = task["base"]
    truth_lfc = task["truth_lfc"]          # {ct: {'alone': arr, 'co': arr}} per sgRNA
    timepoint = task["timepoint"]
    coculture = task["coculture"]          # True for group 3 (cells grown together)
    seed = task["seed"]
    ctx = "co" if coculture else "alone"

    rng = np.random.default_rng(seed)
    r1_path = os.path.join(fastq_dir, f"{sample}_R1.fastq.gz")
    r2_path = os.path.join(fastq_dir, f"{sample}_R2.fastq.gz")

    reads_per_ct = reads_total // len(cell_types)
    n_sg = sg_mat.shape[0]

    with gzip.open(r1_path, "wb", compresslevel=1) as f1, \
            gzip.open(r2_path, "wb", compresslevel=1) as f2:
        for ct in cell_types:
            barcode_codes = np.frombuffer(
                BARCODES[ct].encode().translate(bytes.maketrans(b"ACGT", bytes([0, 1, 2, 3]))),
                dtype=np.uint8,
            )
            if timepoint == "initial":
                lfc_sg = rng.normal(0.0, 0.05, size=n_sg)
            else:
                gene_lfc = truth_lfc[ct][ctx]
                lfc_sg = gene_lfc + rng.normal(0.0, 0.35, size=n_sg)
            p = sgrna_weights(base, lfc_sg)
            counts = rng.multinomial(reads_per_ct, p)
            read_idx = np.repeat(np.arange(n_sg), counts)
            rng.shuffle(read_idx)
            for r1b, r2b in _build_records(read_idx, sg_mat, barcode_codes, rng):
                f1.write(r1b)
                f2.write(r2b)

    return f"{sample}: {reads_total:,} reads ({'+'.join(cell_types)})"


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="Simulate CAMUS two-cell-type screen")
    ap.add_argument("--out", required=True, help="output experiment directory")
    ap.add_argument("--reads-per-sample", type=int, default=5_000_000)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--samples", default=None,
                    help="comma-separated subset of sample names to generate (default: all)")
    ap.add_argument("--assets-only", action="store_true",
                    help="only write library/barcodes/ground truth, no fastq")
    args = ap.parse_args()
    only = set(args.samples.split(",")) if args.samples else None

    data_dir = args.out
    fastq_dir = os.path.join(data_dir, "fastq_files")
    os.makedirs(fastq_dir, exist_ok=True)

    print("Building library and ground truth...")
    sgrnas, truth = build_library_and_truth()
    write_assets(data_dir, sgrnas, truth)
    print(f"  {len(sgrnas)} sgRNAs ({len(GENES)} genes x {SGRNAS_PER_GENE} + {N_NTC} NTC)")

    if args.assets_only:
        print("Assets written; skipping fastq generation.")
        return

    sg_seqs = [s for s, _, _ in sgrnas]
    genes_of_sg = [g for _, g, _ in sgrnas]
    sg_mat = _seq_matrix(sg_seqs)
    base = base_abundance(len(sgrnas), GLOBAL_SEED + 1)

    # Per-cell-type, per-context, per-sgRNA ground-truth logFC vectors.
    def lfc_vec(key):
        return np.array([truth.get(g, {}).get(key, 0.0) if g != "Non-Targeting" else 0.0
                         for g in genes_of_sg])
    truth_lfc = {
        "prolif": {"alone": lfc_vec("prolif_alone"), "co": lfc_vec("prolif_co")},
        "high_cin": {"alone": lfc_vec("high_cin_alone"), "co": lfc_vec("high_cin_co")},
    }

    tasks = []
    sidx = 0
    for group, cts in GROUPS.items():
        for tp in TIMEPOINTS:
            for rep in range(1, N_REPS + 1):
                sample = f"{group}_{tp}_Rep{rep}"
                if only is not None and sample not in only:
                    sidx += 1
                    continue
                tasks.append({
                    "fastq_dir": fastq_dir,
                    "sample": sample,
                    "cell_types": cts,
                    "reads_total": args.reads_per_sample,
                    "sg_mat": sg_mat,
                    "base": base,
                    "truth_lfc": truth_lfc,
                    "timepoint": tp,
                    "coculture": len(cts) > 1,
                    "seed": GLOBAL_SEED + 1000 + sidx,
                })
                sidx += 1

    print(f"Generating {len(tasks)} samples x {args.reads_per_sample:,} reads "
          f"using {args.workers} workers...")
    with open(os.path.join(data_dir, "_progress.log"), "w") as plog:
        plog.write(f"start {len(tasks)} samples\n")
        plog.flush()
        with Pool(args.workers) as pool:
            for msg in pool.imap_unordered(simulate_sample, tasks):
                print("  done", msg, flush=True)
                plog.write("done " + msg + "\n")
                plog.flush()
    print("All samples generated.")


if __name__ == "__main__":
    main()
