# CAMUS simulated experiment

A two-cell-type, multiplexed CRISPR screen that demonstrates what CAMUS is for:
recovering per-cell-type screen results from a pooled co-culture, including a
**non-cell-autonomous (rescue) interaction** that is only visible once the
co-culture is demultiplexed.

The paired `fastq.gz` are already generated under `data/fastq_files/` (18 samples,
5,000,000 read pairs each, ~2.3 GB). You can re-run any step yourself.

## Design

Cell types (same library, different CAMUS barcodes):

| cell type | meaning | barcode |
|---|---|---|
| `prolif` | proliferative | `ATGCAGGG` |
| `high_cin` | high chromosomal instability | `GTTGCAGC` |

Library: 100 secretome / surface-receptor genes × 6 sgRNAs + 500 non-targeting
controls = 1,100 sgRNAs.

Groups × timepoints × replicates = 18 samples:

| group | contents | demux? |
|---|---|---|
| `G1_prolif` | prolif only | trivial (one cell type) |
| `G2_high_cin` | high_cin only | trivial |
| `G3_mixed` | both, pooled together | **separated by barcode** |

Timepoints: `initial`, `D21`. Replicates: 3. In `G3` the 5M reads are split ~50/50
between the two cell types.

### Ground truth (`data/ground_truth.csv`)

Each gene has a logFC for the **alone** context (pure culture, groups 1/2) and the
**co**-culture context (group 3):

| category | genes | effect |
|---|---|---|
| shared_essential | 10 | depleted in both cell types (alone & co) |
| prolif_essential | 12 | depleted in prolif (alone & co) |
| highcin_essential | 12 | depleted in high_cin (alone & co) |
| prolif_enriched | 5 | enriched in prolif |
| highcin_enriched | 5 | enriched in high_cin |
| **prolif_rescue** | 2 (IGF1, TGFA) | **depleted in prolif ALONE, rescued in co-culture** |
| **highcin_rescue** | 2 (BMP4, DKK1) | **depleted in high_cin ALONE, rescued in co-culture** |
| neutral | 52 | no effect |

The rescue genes model a secreted factor whose knockout is lethal in pure culture
but rescued in trans by the neighbouring cell type in co-culture.

## Run it

MAGeCK must be on your `PATH` (or pass its location). Then, from the repo root:

```bash
# full run (re-simulates, demultiplexes, runs MAGeCK, compares)
bash sim_experiment/run_experiment.sh

# the reads already exist, so to skip simulation just run the last two steps:
python3 sim_experiment/run_demux_all.py  --data sim_experiment/data --out sim_experiment/demux
python3 sim_experiment/analyze_experiment.py \
    --data sim_experiment/data --demux sim_experiment/demux --out sim_experiment/results
```

If your MAGeCK lives in a conda env: add `--conda-env <env>` (or set `CONDA_ENV=<env>`
for `run_experiment.sh`); for an explicit binary use `--mageck-binary /path/to/mageck`.

Environment knobs for `run_experiment.sh`: `READS`, `WORKERS`, `FDR`, `MAGECK_BIN`,
`CONDA_ENV`.

## What the analysis does

For each populated (group, cell type) it builds an initial-vs-D21 count matrix and
runs MAGeCK RRA (non-targeting controls used for normalization), then compares:

- `prolif`  : `G1_prolif` (pure) vs `G3_mixed` prolif lane (demultiplexed)
- `high_cin`: `G2_high_cin` (pure) vs `G3_mixed` high_cin lane (demultiplexed)

Outputs in `sim_experiment/results/`:

- `REPORT.md` — recovery vs ground truth (precision/recall/F1, logFC correlation),
  pure-vs-demuxed concordance, and a dedicated **context-dependent (rescue) genes**
  table.
- `prolif_pure_vs_mixed.png`, `high_cin_pure_vs_mixed.png` — logFC scatter coloured
  by ground-truth category; rescue genes drawn as black stars. Cell-autonomous hits
  sit on the diagonal; **rescue genes fall off it** (depleted in pure, ~0 in mixed).
- `*_gene_results.txt`, `*_matrix.txt` — per-unit MAGeCK results and count matrices.

### Expected result

Cell-autonomous essential/enriched genes are recovered concordantly whether the
cell type was screened pure or demultiplexed from the co-culture (high logFC
correlation, hits on the diagonal). The four rescue genes are the exception by
design: significant and depleted in pure culture, but not significant (~0) in the
demultiplexed co-culture — the interaction CAMUS exists to reveal.

## Files

```
sim_experiment/
  simulate_experiment.py   library + ground truth + parallel fastq.gz generator
  run_demux_all.py         CAMUS demultiplexing of all 18 samples
  analyze_experiment.py    MAGeCK RRA + group comparison + report/figures
  run_experiment.sh        one-command end-to-end runner
  data/
    fastq_files/           36 fastq.gz (already generated)
    combined_library.txt   sgRNA<TAB>Gene
    celltype_barcodes.csv  Celltype,Barcode
    ground_truth.csv       per-gene alone/co logFC + category
    sgrna_id_map.tsv       sequence -> sgID, Gene
```
