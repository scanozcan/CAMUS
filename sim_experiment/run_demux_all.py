#!/usr/bin/env python3
"""
Demultiplex all simulated samples with CAMUS.

Runs camus.demux.FastDemultiplexer on every sample's paired FASTQ. Groups 1/2
contain a single cell type (the other cell type simply receives ~0 reads); group 3
is genuinely separated into prolif and high_cin by the cell-type barcode.

Produces, per sample, <prefix>_count.txt with columns: sgRNA, Gene, high_cin, prolif.
Runs samples in parallel (one worker per sample).
"""

import argparse
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from camus.demux.demux_hash_table import FastDemultiplexer

GROUPS = {
    "G1_prolif": ["prolif"],
    "G2_high_cin": ["high_cin"],
    "G3_mixed": ["prolif", "high_cin"],
}
TIMEPOINTS = ["initial", "D21"]
N_REPS = 3


def _demux_one(args):
    data_dir, out_dir, sample = args
    fastq_dir = os.path.join(data_dir, "fastq_files")
    r1 = os.path.join(fastq_dir, f"{sample}_R1.fastq.gz")
    r2 = os.path.join(fastq_dir, f"{sample}_R2.fastq.gz")
    if not (os.path.exists(r1) and os.path.exists(r2)):
        return f"{sample}: MISSING fastq"
    demux = FastDemultiplexer(
        barcode_csv=os.path.join(data_dir, "celltype_barcodes.csv"),
        barcode_start=22, barcode_length=8,
        grna_start=12, grna_length=20,
        max_barcode_mismatches=1, allow_grna_mismatch=True,
        min_qual=0,
    )
    stats = demux.process(r1, r2, os.path.join(data_dir, "combined_library.txt"),
                          os.path.join(out_dir, sample))
    assigned = {k: v for k, v in stats["assigned"].items()}
    return (f"{sample}: total={stats['total']:,} assigned={sum(assigned.values()):,} "
            f"{assigned} invalid={stats['invalid_grna']:,} amb_bc={stats['ambiguous_barcode']:,}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="experiment data dir (with fastq_files/)")
    ap.add_argument("--out", required=True, help="demux output dir")
    ap.add_argument("--workers", type=int, default=4)
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    samples = [f"{g}_{tp}_Rep{r}"
               for g in GROUPS for tp in TIMEPOINTS for r in range(1, N_REPS + 1)]
    tasks = [(args.data, args.out, s) for s in samples]

    with open(os.path.join(args.out, "_demux_progress.log"), "w") as plog:
        with Pool(args.workers) as pool:
            for msg in pool.imap_unordered(_demux_one, tasks):
                print("done", msg, flush=True)
                plog.write("done " + msg + "\n"); plog.flush()
    print("All samples demultiplexed.")


if __name__ == "__main__":
    main()
