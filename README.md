# CAMUS — CRISPR Analysis of Multiplexed Screens

Cell-type-aware demultiplexing and MAGeCK analysis for pooled CRISPR screens.

CAMUS takes paired-end FASTQs from a pooled CRISPR screen in which each cell
type carries a short cell-type barcode, splits the reads by cell type, builds a
per-cell-type sgRNA count matrix, and runs the official
[MAGeCK](https://sourceforge.net/p/mageck/wiki/Home/) (Li et al., *Genome
Biology* 2014) in either RRA (pairwise) or MLE (multi-condition) mode. It then
produces volcano plots, QC plots, heatmaps, and text reports.

## Install

```bash
# Recommended: conda env that also installs MAGeCK from bioconda
conda env create -f environment.yml
conda activate camus

# Or, if MAGeCK is already available, install just the Python package
pip install -e .
```

CAMUS locates MAGeCK automatically: an explicit `--mageck-binary`, then
`mageck` on `PATH`, then `conda run -n <env> mageck` if you pass `--conda-env`.
It fails early with the detected version if MAGeCK cannot be found.

## Input layout

```
data_dir/
  celltype_barcodes.csv        # columns: Celltype,Barcode
  combined_library.txt         # TSV: sgRNA_sequence <TAB> Gene  (no header)
  <celltype>_library.txt       # per-cell-type sub-libraries (lowercased name)
  fastq_files/
    <Sample>_R1.fastq[.gz]     # R1 carries the gRNA
    <Sample>_R2.fastq[.gz]     # R2 carries the cell-type barcode
```

Sample names come from the config: `Control_Rep1`, `Condition1_Rep1`, etc.
Both plain `.fastq` and gzipped `.fastq.gz` are accepted.

## Configure

Copy `config/experiment_template.yaml` and edit it. The `geometry` block defines
the read layout and matching parameters, so you can adapt to your vector without
touching code:

```yaml
experiment:
  name: My_Screen
  control: { name: Control, replicates: 3 }
  conditions:
    - { name: Condition1, replicates: 3 }
  cell_types:
    - { name: prolif,   barcode: ATGCAGGG }
    - { name: high_cin, barcode: GTTGCAGC }
  geometry:
    barcode_start: 22        # where the barcode sits in R2 (0-based)
    barcode_length: 8
    grna_start: 12           # fixed-offset gRNA start in R1
    grna_length: 20
    max_barcode_mismatches: 1
    min_qual: 0              # set >0 to filter low-quality barcode/gRNA windows
    grna_anchor: null        # e.g. a scaffold seq to locate staggered gRNAs
    anchor_max_offset: 40
  control_sgrna: null        # path to a control-sgRNA list (auto-detected if null)
```

If your reads are staggered or scaffolded so a fixed `grna_start` doesn't work,
set `grna_anchor` to the constant sequence immediately upstream of the gRNA and
CAMUS will search for it instead of using a fixed offset.

## Grouped experiments (per-group baselines)

By default the pipeline compares one shared `control` against each `condition`. If
your design has several **groups that each carry their own baseline** — e.g. cells
screened *alone* vs. the same cells *co-cultured*, each with its own `initial`
timepoint — use a `groups:` block instead. `cell_types` then becomes a global
catalogue (name → barcode) that groups reference by name; each group declares its
own `control` and `conditions` and is analysed against **its own** baseline:

```yaml
experiment:
  name: My_CoCulture_Screen
  cell_types:                      # global catalogue
    - { name: prolif,   barcode: ATGCAGGG }
    - { name: high_cin, barcode: GTTGCAGC }
  groups:
    - name: G1_prolif               # pure culture
      cell_types: [prolif]
      control:    { name: initial, replicates: 3 }
      conditions: [ { name: D21, replicates: 3 } ]
    - name: G3_mixed                # co-culture, demuxed into both lanes
      cell_types: [prolif, high_cin]
      control:    { name: initial, replicates: 3 }
      conditions: [ { name: D21, replicates: 3 } ]
  geometry: { ... }
```

FASTQ sample names are `{group}_{name}_Rep{i}` (e.g. `G1_prolif_initial_Rep1_R1.fastq.gz`),
so each group's reads are picked up automatically. The pipeline writes a self-contained
`<output>/<group>/{demux,counts,results}/` subtree per group and runs
`{control}_vs_{condition}` per cell type against that group's own baseline. Configs
*without* a `groups:` block behave exactly as before.

The visualization step is group-aware: it produces per-group plots under
`<output>/<group>/` **and** a `cross_group/` folder that compares the *same cell type
across groups* (e.g. `prolif` alone vs. co-cultured) — a logFC concordance scatter with
off-diagonal (context-dependent) genes labelled, a hit-overlap Venn, and a
`CROSS_GROUP_REPORT.md`. A worked end-to-end example lives in `sim_experiment/`.

## Run

```bash
# RRA: each condition vs control, per cell type
python scripts/run_mageck_pipeline.py \
    --config config/experiment_template.yaml \
    --data-dir data_dir \
    --output pipeline_outputs

# MLE: all conditions jointly vs control, per cell type
python scripts/run_mageck_mle_pipeline.py \
    --config config/experiment_template.yaml \
    --data-dir data_dir \
    --output pipeline_outputs_mle
```

Non-targeting controls are detected from the library by gene name (e.g.
`control`, `non-targeting`) and passed to MAGeCK as `--control-sgrna` for
control-based normalization. Disable with `--no-control-sgrna`.

## Visualize

```bash
python scripts/post_analysis_visualization.py \
    --results-dir pipeline_outputs/results \
    --config config/experiment_template.yaml \
    --output visualizations

python scripts/post_analysis_visualization_mle.py \
    --results-dir pipeline_outputs_mle/results_mle \
    --config config/experiment_template.yaml \
    --output visualizations_mle
```

## Demultiplex only

```bash
python -m camus.demux.demux_hash_table \
    --read1 R1.fastq.gz --read2 R2.fastq.gz \
    --library combined_library.txt \
    --celltype-barcodes celltype_barcodes.csv \
    --output-prefix out/Sample \
    --grna-anchor GTTTAAGAGC --min-qual 20
```

## Tests

```bash
pip install pytest
pytest
```

The suite runs on tiny synthetic FASTQs with known answers and covers the
ambiguous-1-mismatch-variant and barcode-tie fixes, quality filtering, anchor
search, gzip input, config parsing, and the MAGeCK wrapper plumbing.

## Validation status & next step

CAMUS has so far been exercised on simulated data and the included unit tests.
Before trusting it on a new platform, validate against a published screen with
known ground truth — e.g. confirm recovery of core-essential genes (Hart et al.)
on a public essentiality dataset, and, for the demux premise, compare barcode
assignments and unassigned/doublet rates against a real hashed/multiplexed study.

## Layout

```
CAMUS/
  camus/
    config/        experiment + library configuration (with read geometry)
    demux/         hash-table demultiplexer
    counting/      optional standalone sgRNA counting
    wrappers/      MAGeCK RRA/MLE wrapper + MAGeCK discovery
    visualization/ volcano, QC, comparison, heatmap, report
  scripts/         RRA + MLE pipelines, RRA + MLE visualization
  config/          experiment_template.yaml
  tests/           pytest suite
  environment.yml  conda env incl. pinned MAGeCK
  pyproject.toml
```
