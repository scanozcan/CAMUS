"""Tests for the CAMUS demultiplexer correctness fixes.

Geometry for these tests: barcode at R2[0:8], gRNA at R1[0:4].

Library:
    AAAA  GeneA
    AGGA  GeneB
    TTTT  control_1   (gene name marks it a non-targeting control)

Key constructed cases:
    AAAC  -> unique 1-mismatch of AAAA only            (should match GeneA)
    AAGA  -> reachable from AAAA (pos2) AND AGGA (pos1) (ambiguous -> dropped)

Barcodes (differ by 2 -> an in-between observed barcode is a tie):
    Type1  AAAAAAAA
    Type2  AAAAAATT
    AAAAAAAT -> distance 1 to BOTH -> ambiguous tie (dropped)
"""

import os
import pytest

from camus.demux.demux_hash_table import FastDemultiplexer

LIBRARY = "AAAA\tGeneA\nAGGA\tGeneB\nTTTT\tcontrol_1\n"
BARCODES = "Celltype,Barcode\nType1,AAAAAAAA\nType2,AAAAAATT\n"


def _write(path, text):
    with open(path, "w") as f:
        f.write(text)


def _fastq_record(idx, seq, qual=None):
    qual = qual or ("I" * len(seq))
    return f"@read{idx}\n{seq}\n+\n{qual}\n"


def _make_inputs(tmp_path, r1_reads, r2_reads):
    lib = tmp_path / "library.txt"
    bc = tmp_path / "barcodes.csv"
    r1 = tmp_path / "R1.fastq"
    r2 = tmp_path / "R2.fastq"
    _write(lib, LIBRARY)
    _write(bc, BARCODES)
    _write(r1, "".join(_fastq_record(i, s) if isinstance(s, str)
                        else _fastq_record(i, s[0], s[1]) for i, s in enumerate(r1_reads)))
    _write(r2, "".join(_fastq_record(i, s) if isinstance(s, str)
                        else _fastq_record(i, s[0], s[1]) for i, s in enumerate(r2_reads)))
    return str(lib), str(bc), str(r1), str(r2)


def _demux(bc, **kw):
    return FastDemultiplexer(
        barcode_csv=bc, barcode_start=0, barcode_length=8,
        grna_start=0, grna_length=4, max_barcode_mismatches=1,
        **kw,
    )


def test_exact_and_unique_mismatch_assign(tmp_path):
    # AAAA exact + Type1; AAAC unique 1-mismatch + Type1
    lib, bc, r1, r2 = _make_inputs(
        tmp_path,
        r1_reads=["AAAA", "AAAC"],
        r2_reads=["AAAAAAAA", "AAAAAAAA"],
    )
    d = _demux(bc)
    stats = d.process(r1, r2, lib, str(tmp_path / "out"))
    assert stats["exact_match"] == 1
    assert stats["one_mismatch"] == 1
    assert stats["assigned"]["Type1"] == 2
    assert stats["invalid_grna"] == 0


def test_ambiguous_variant_dropped(tmp_path):
    # AAGA is a 1-mismatch variant of BOTH AAAA and AGGA -> must be dropped.
    lib, bc, r1, r2 = _make_inputs(
        tmp_path, r1_reads=["AAGA"], r2_reads=["AAAAAAAA"]
    )
    d = _demux(bc)
    stats = d.process(r1, r2, lib, str(tmp_path / "out"))
    assert d.n_ambiguous_variants > 0
    assert stats["invalid_grna"] == 1
    assert sum(stats["assigned"].values()) == 0


def test_barcode_tie_rejected(tmp_path):
    # Barcode AAAAAAAT is distance 1 from both Type1 and Type2 -> ambiguous.
    lib, bc, r1, r2 = _make_inputs(
        tmp_path, r1_reads=["AAAA"], r2_reads=["AAAAAAAT"]
    )
    d = _demux(bc)
    stats = d.process(r1, r2, lib, str(tmp_path / "out"))
    assert stats["ambiguous_barcode"] == 1
    assert sum(stats["assigned"].values()) == 0


def test_quality_filter(tmp_path):
    # Same good read, but gRNA bases have low quality ('#': Phred 2).
    lib, bc, r1, r2 = _make_inputs(
        tmp_path,
        r1_reads=[("AAAA", "####")],
        r2_reads=[("AAAAAAAA", "IIIIIIII")],
    )
    d = _demux(bc, min_qual=20)
    stats = d.process(r1, r2, lib, str(tmp_path / "out"))
    assert stats["low_quality"] == 1
    assert sum(stats["assigned"].values()) == 0


def test_control_sgrna_detection_and_file(tmp_path):
    lib, bc, r1, r2 = _make_inputs(
        tmp_path, r1_reads=["TTTT"], r2_reads=["AAAAAAAA"]
    )
    out = str(tmp_path / "out")
    d = _demux(bc)
    d.process(r1, r2, lib, out)
    control_file = out + "_control_sgrnas.txt"
    assert os.path.exists(control_file)
    with open(control_file) as f:
        controls = f.read().split()
    assert "TTTT" in controls


def test_anchor_based_extraction(tmp_path):
    # Anchor 'GG' precedes the 4bp gRNA; R1 has a staggered prefix.
    lib, bc, r1, r2 = _make_inputs(
        tmp_path,
        r1_reads=["NNNGGAAAA"],   # anchor GG at idx 3, gRNA AAAA follows
        r2_reads=["AAAAAAAA"],
    )
    d = FastDemultiplexer(
        barcode_csv=bc, barcode_start=0, barcode_length=8,
        grna_length=4, max_barcode_mismatches=1,
        grna_anchor="GG", anchor_max_offset=10,
    )
    stats = d.process(r1, r2, lib, str(tmp_path / "out"))
    assert stats["assigned"]["Type1"] == 1
    assert stats["exact_match"] == 1


def test_gzip_input(tmp_path):
    import gzip
    lib = tmp_path / "library.txt"
    bc = tmp_path / "barcodes.csv"
    _write(lib, LIBRARY)
    _write(bc, BARCODES)
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    with gzip.open(r1, "wt") as f:
        f.write(_fastq_record(0, "AAAA"))
    with gzip.open(r2, "wt") as f:
        f.write(_fastq_record(0, "AAAAAAAA"))
    d = _demux(str(bc))
    stats = d.process(str(r1), str(r2), str(lib), str(tmp_path / "out"))
    assert stats["assigned"]["Type1"] == 1
