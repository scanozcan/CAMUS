"""Tests for experiment configuration parsing and validation."""

import textwrap
import pytest

from camus.config.experiment import ExperimentConfig, ReadGeometryConfig


BASE_YAML = """
experiment:
  name: Test_Screen
  control:
    name: Control
    replicates: 2
  conditions:
    - name: Cond1
      replicates: 2
  cell_types:
    - name: TypeA
      barcode: AAAAAAAA
    - name: TypeB
      barcode: GGGGGGGG
  geometry:
    barcode_start: 10
    barcode_length: 8
    grna_start: 0
    grna_length: 20
    max_barcode_mismatches: 1
    min_qual: 15
  library_config: default
"""


def _write_yaml(tmp_path, text):
    p = tmp_path / "exp.yaml"
    p.write_text(textwrap.dedent(text))
    return str(p)


def test_geometry_parsed(tmp_path):
    cfg = ExperimentConfig.from_yaml(_write_yaml(tmp_path, BASE_YAML))
    assert isinstance(cfg.geometry, ReadGeometryConfig)
    assert cfg.geometry.barcode_start == 10
    assert cfg.geometry.min_qual == 15
    assert cfg.geometry.grna_length == 20
    assert cfg.validate() == []


def test_default_geometry_when_absent(tmp_path):
    yaml_no_geom = BASE_YAML.split("  geometry:")[0] + "  library_config: default\n"
    cfg = ExperimentConfig.from_yaml(_write_yaml(tmp_path, yaml_no_geom))
    assert cfg.geometry.barcode_start == 22  # default
    assert cfg.geometry.grna_start == 12


def test_barcode_length_validation(tmp_path):
    bad = BASE_YAML.replace("barcode: AAAAAAAA", "barcode: AAA")  # 3bp vs expected 8
    # from_yaml validates and raises on invalid config
    with pytest.raises(ValueError):
        ExperimentConfig.from_yaml(_write_yaml(tmp_path, bad))


def test_comparisons_and_samples(tmp_path):
    cfg = ExperimentConfig.from_yaml(_write_yaml(tmp_path, BASE_YAML))
    samples = cfg.get_all_sample_names()
    assert "Control_Rep1" in samples and "Cond1_Rep2" in samples
    comps = cfg.get_comparisons()
    assert len(comps) == 1
    assert comps[0].name == "Control_vs_Cond1"
