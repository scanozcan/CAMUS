"""Tests for the MAGeCK wrapper that don't require MAGeCK to be installed."""

import pandas as pd
import pytest

from camus.wrappers.mageck_wrapper import (
    MAGeCKWrapper,
    MAGeCKNotFoundError,
    resolve_mageck_invocation,
)


def test_resolve_explicit_binary(tmp_path):
    fake = tmp_path / "mageck"
    fake.write_text("#!/bin/sh\necho 0.5.9\n")
    fake.chmod(0o755)
    inv = resolve_mageck_invocation(mageck_binary=str(fake))
    assert inv == [str(fake)]


def test_resolve_missing_raises(monkeypatch):
    # No mageck on PATH, no conda env -> clear error.
    monkeypatch.setattr("shutil.which", lambda name: None)
    with pytest.raises(MAGeCKNotFoundError):
        resolve_mageck_invocation()


def test_design_matrix(tmp_path):
    w = MAGeCKWrapper(verify=False)
    out = tmp_path / "design.txt"
    w.generate_design_matrix(
        control_samples=["Control_Rep1", "Control_Rep2"],
        condition_groups={
            "Cond1": ["Cond1_Rep1", "Cond1_Rep2"],
            "Cond2": ["Cond2_Rep1", "Cond2_Rep2"],
        },
        output_file=str(out),
    )
    df = pd.read_csv(out, sep="\t")
    assert list(df.columns) == ["Samples", "baseline", "Cond1", "Cond2"]
    assert (df["baseline"] == 1).all()
    row = df[df["Samples"] == "Cond1_Rep1"].iloc[0]
    assert row["Cond1"] == 1 and row["Cond2"] == 0


def test_parse_gene_results_threshold(tmp_path):
    # Minimal MAGeCK-style gene_summary with one clear hit and one non-hit.
    summary = tmp_path / "x.gene_summary.txt"
    summary.write_text(
        "id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\tneg|lfc\t"
        "pos|score\tpos|p-value\tpos|fdr\tpos|rank\tpos|lfc\n"
        "GeneA\t6\t0.001\t0.0001\t0.01\t1\t-2.0\t0.9\t0.9\t0.95\t900\t0.1\n"
        "GeneB\t6\t0.4\t0.4\t0.5\t500\t-0.1\t0.4\t0.4\t0.5\t500\t0.1\n"
    )
    w = MAGeCKWrapper(verify=False)
    df = w.parse_gene_results(str(summary), fdr_threshold=0.05)
    a = df[df["Gene"] == "GeneA"].iloc[0]
    b = df[df["Gene"] == "GeneB"].iloc[0]
    assert a["Significance"] == "significant"
    assert a["Direction"] == "depleted"
    assert b["Significance"] == "not_significant"
