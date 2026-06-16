#!/usr/bin/env python3
"""
Experiment configuration for CAMUS.

Provides YAML-based configuration for multi-condition experiments, including
the read geometry used by the demultiplexer.
"""

import yaml
from dataclasses import dataclass, field
from typing import List, Dict, Optional
from pathlib import Path


@dataclass
class ReadGeometryConfig:
    """Read geometry + matching parameters for the demultiplexer.

    All positions are 0-based. ``barcode_length`` also defines the expected
    cell-type barcode length used for validation.
    """
    barcode_start: int = 22
    barcode_length: int = 8
    grna_start: int = 12
    grna_length: int = 20
    max_barcode_mismatches: int = 1
    allow_grna_mismatch: bool = True
    min_qual: int = 0
    grna_anchor: Optional[str] = None
    anchor_max_offset: int = 40

    def validate(self) -> List[str]:
        errors = []
        for attr in ('barcode_start', 'grna_start', 'min_qual', 'anchor_max_offset'):
            if getattr(self, attr) < 0:
                errors.append(f"geometry.{attr} must be >= 0")
        if self.barcode_length < 1:
            errors.append("geometry.barcode_length must be >= 1")
        if self.grna_length < 1:
            errors.append("geometry.grna_length must be >= 1")
        if self.max_barcode_mismatches < 0:
            errors.append("geometry.max_barcode_mismatches must be >= 0")
        return errors


@dataclass
class CellTypeConfig:
    """Configuration for a single cell type."""
    name: str
    barcode: str
    barcode_length: int = 8  # expected length, taken from read geometry

    def validate(self) -> List[str]:
        """Validate cell type configuration."""
        errors = []
        if not self.name:
            errors.append("Cell type name cannot be empty")
        if not self.barcode or len(self.barcode) != self.barcode_length:
            errors.append(
                f"Barcode must be {self.barcode_length}bp, got: {self.barcode}")
        return errors


@dataclass
class ConditionConfig:
    """Configuration for a single experimental condition."""
    name: str
    replicates: int = 3

    def get_sample_names(self) -> List[str]:
        """Get sample names for this condition."""
        return [f"{self.name}_Rep{i+1}" for i in range(self.replicates)]

    def validate(self) -> List[str]:
        """Validate condition configuration."""
        errors = []
        if not self.name:
            errors.append("Condition name cannot be empty")
        if self.replicates < 1:
            errors.append(f"Replicates must be >= 1, got: {self.replicates}")
        return errors


@dataclass
class ControlConfig:
    """Configuration for control samples."""
    name: str = "Control"
    replicates: int = 3

    def get_sample_names(self) -> List[str]:
        """Get sample names for control."""
        return [f"{self.name}_Rep{i+1}" for i in range(self.replicates)]

    def validate(self) -> List[str]:
        """Validate control configuration."""
        errors = []
        if not self.name:
            errors.append("Control name cannot be empty")
        if self.replicates < 1:
            errors.append(f"Replicates must be >= 1, got: {self.replicates}")
        return errors


@dataclass
class GroupConfig:
    """An experimental group with its own matched baseline.

    A group has its own control (e.g. an ``initial`` timepoint) and one or more
    conditions (e.g. ``D21``), analysed over a subset of the experiment's cell
    types. Each group is compared against its OWN control rather than a single
    global baseline. FASTQ sample names are ``{group}_{name}_Rep{i}`` so a group
    named ``G1_prolif`` with control ``initial`` maps to ``G1_prolif_initial_Rep1``.
    """
    name: str
    control: ControlConfig
    conditions: List['ConditionConfig']
    cell_types: List[str] = field(default_factory=list)  # names; empty => all

    def _prefix(self) -> str:
        return f"{self.name}_" if self.name else ""

    def control_sample_names(self) -> List[str]:
        return [f"{self._prefix()}{self.control.name}_Rep{i+1}"
                for i in range(self.control.replicates)]

    def condition_sample_names(self, condition: 'ConditionConfig') -> List[str]:
        return [f"{self._prefix()}{condition.name}_Rep{i+1}"
                for i in range(condition.replicates)]

    def get_sample_names(self) -> List[str]:
        names = self.control_sample_names()
        for condition in self.conditions:
            names.extend(self.condition_sample_names(condition))
        return names

    def validate(self) -> List[str]:
        errors = []
        if not self.name:
            errors.append("Group name cannot be empty")
        errors.extend([f"Group '{self.name}' control: {e}"
                       for e in self.control.validate()])
        if not self.conditions:
            errors.append(f"Group '{self.name}': at least one condition required")
        for condition in self.conditions:
            errors.extend([f"Group '{self.name}' condition: {e}"
                           for e in condition.validate()])
        return errors


@dataclass
class Comparison:
    """A single control vs condition comparison (optionally within a group)."""
    control_name: str
    condition_name: str
    control_samples: List[str]
    condition_samples: List[str]
    group: Optional[str] = None

    @property
    def name(self) -> str:
        """Get comparison name (group-prefixed when part of a group)."""
        base = f"{self.control_name}_vs_{self.condition_name}"
        return f"{self.group}__{base}" if self.group else base


@dataclass
class ExperimentConfig:
    """Master experiment configuration."""
    name: str
    control: ControlConfig
    conditions: List[ConditionConfig]
    cell_types: List[CellTypeConfig]
    library_config: Optional[str] = 'default'  # Library preset name
    geometry: ReadGeometryConfig = field(default_factory=ReadGeometryConfig)
    control_sgrna: Optional[str] = None  # path to non-targeting control sgRNA list
    groups: List[GroupConfig] = field(default_factory=list)  # optional grouped design

    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'ExperimentConfig':
        """
        Load experiment configuration from YAML file.

        Args:
            yaml_path: Path to YAML configuration file

        Returns:
            ExperimentConfig instance

        Raises:
            FileNotFoundError: If YAML file doesn't exist
            yaml.YAMLError: If YAML is invalid
            ValueError: If configuration is invalid

        Examples:
            >>> config = ExperimentConfig.from_yaml('config/experiment.yaml')
            >>> print(config.name)
            Multi_Condition_Screen
        """
        yaml_path = Path(yaml_path)
        if not yaml_path.exists():
            raise FileNotFoundError(f"Config file not found: {yaml_path}")

        with open(yaml_path, 'r') as f:
            data = yaml.safe_load(f)

        # Parse experiment section
        exp_data = data.get('experiment', {})

        # Parse control
        control_data = exp_data.get('control', {})
        control = ControlConfig(
            name=control_data.get('name', 'Control'),
            replicates=control_data.get('replicates', 3)
        )

        # Parse conditions
        conditions = []
        for cond_data in exp_data.get('conditions', []):
            condition = ConditionConfig(
                name=cond_data.get('name'),
                replicates=cond_data.get('replicates', 3)
            )
            conditions.append(condition)

        # Parse read geometry (optional; sensible defaults)
        geom_data = exp_data.get('geometry', {}) or {}
        geometry = ReadGeometryConfig(
            barcode_start=geom_data.get('barcode_start', 22),
            barcode_length=geom_data.get('barcode_length', 8),
            grna_start=geom_data.get('grna_start', 12),
            grna_length=geom_data.get('grna_length', 20),
            max_barcode_mismatches=geom_data.get('max_barcode_mismatches', 1),
            allow_grna_mismatch=geom_data.get('allow_grna_mismatch', True),
            min_qual=geom_data.get('min_qual', 0),
            grna_anchor=geom_data.get('grna_anchor', None),
            anchor_max_offset=geom_data.get('anchor_max_offset', 40),
        )

        # Parse cell types (barcode length taken from geometry)
        cell_types = []
        for ct_data in exp_data.get('cell_types', []):
            cell_type = CellTypeConfig(
                name=ct_data.get('name'),
                barcode=ct_data.get('barcode'),
                barcode_length=geometry.barcode_length,
            )
            cell_types.append(cell_type)

        # Parse optional groups (each with its own control baseline + conditions)
        groups = []
        for g_data in exp_data.get('groups', []) or []:
            g_ctrl = g_data.get('control', {}) or {}
            group = GroupConfig(
                name=g_data.get('name'),
                control=ControlConfig(
                    name=g_ctrl.get('name', 'initial'),
                    replicates=g_ctrl.get('replicates', 3),
                ),
                conditions=[
                    ConditionConfig(
                        name=c.get('name'),
                        replicates=c.get('replicates', 3),
                    )
                    for c in (g_data.get('conditions', []) or [])
                ],
                cell_types=list(g_data.get('cell_types', []) or []),
            )
            groups.append(group)

        # Create config
        config = cls(
            name=exp_data.get('name', 'Experiment'),
            control=control,
            conditions=conditions,
            cell_types=cell_types,
            library_config=exp_data.get('library_config', 'default'),
            geometry=geometry,
            control_sgrna=exp_data.get('control_sgrna', None),
            groups=groups,
        )

        # Validate
        errors = config.validate()
        if errors:
            raise ValueError(f"Invalid configuration:\n" + "\n".join(f"  - {e}" for e in errors))

        return config

    def get_all_sample_names(self) -> List[str]:
        """Get all sample names (control + all conditions)."""
        samples = self.control.get_sample_names()
        for condition in self.conditions:
            samples.extend(condition.get_sample_names())
        return samples

    def get_comparisons(self) -> List[Comparison]:
        """
        Generate all Control vs Condition comparisons.

        Returns:
            List of Comparison objects
        """
        comparisons = []
        for condition in self.conditions:
            comparison = Comparison(
                control_name=self.control.name,
                condition_name=condition.name,
                control_samples=self.control.get_sample_names(),
                condition_samples=condition.get_sample_names()
            )
            comparisons.append(comparison)
        return comparisons

    def get_groups(self) -> List[GroupConfig]:
        """Return the experiment's groups.

        In grouped mode, returns the configured groups. In flat mode, synthesizes
        a single unnamed group from the top-level control/conditions so callers
        can treat both layouts uniformly. The unnamed group adds no name prefix,
        so flat-mode sample names and comparison names are unchanged.
        """
        if self.groups:
            return self.groups
        return [GroupConfig(
            name="",
            control=self.control,
            conditions=self.conditions,
            cell_types=[ct.name for ct in self.cell_types],
        )]

    def get_group_cell_types(self, group: GroupConfig) -> List['CellTypeConfig']:
        """Resolve a group's cell-type names against the global catalog.

        An empty ``group.cell_types`` means all cell types.
        """
        catalog = {ct.name: ct for ct in self.cell_types}
        if not group.cell_types:
            return list(self.cell_types)
        return [catalog[name] for name in group.cell_types if name in catalog]

    def get_cell_type_barcodes(self) -> Dict[str, str]:
        """Get mapping of cell type to barcode."""
        return {ct.name: ct.barcode for ct in self.cell_types}

    def validate(self) -> List[str]:
        """
        Validate experiment configuration.

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        # Validate name
        if not self.name:
            errors.append("Experiment name cannot be empty")

        # Validate control (always checked; in grouped mode it is a harmless default)
        errors.extend([f"Control: {e}" for e in self.control.validate()])

        # Require either top-level conditions (flat mode) or groups (grouped mode)
        if not self.conditions and not self.groups:
            errors.append("At least one condition or group must be specified")

        for i, condition in enumerate(self.conditions):
            cond_errors = condition.validate()
            errors.extend([f"Condition {i+1}: {e}" for e in cond_errors])

        # Check for duplicate condition names
        condition_names = [c.name for c in self.conditions]
        if len(condition_names) != len(set(condition_names)):
            errors.append("Duplicate condition names detected")

        # Validate groups (grouped mode)
        group_names = [g.name for g in self.groups]
        if len(group_names) != len(set(group_names)):
            errors.append("Duplicate group names detected")
        catalog_names = {ct.name for ct in self.cell_types}
        for group in self.groups:
            errors.extend(group.validate())
            for ct_name in group.cell_types:
                if ct_name not in catalog_names:
                    errors.append(
                        f"Group '{group.name}': cell type '{ct_name}' not in cell_types catalog")

        # Validate cell types
        if not self.cell_types:
            errors.append("At least one cell type must be specified")

        for i, cell_type in enumerate(self.cell_types):
            ct_errors = cell_type.validate()
            errors.extend([f"Cell type {i+1}: {e}" for e in ct_errors])

        # Check for duplicate cell type names or barcodes
        ct_names = [ct.name for ct in self.cell_types]
        if len(ct_names) != len(set(ct_names)):
            errors.append("Duplicate cell type names detected")

        ct_barcodes = [ct.barcode for ct in self.cell_types]
        if len(ct_barcodes) != len(set(ct_barcodes)):
            errors.append("Duplicate cell type barcodes detected")

        # Validate read geometry
        errors.extend([f"Geometry: {e}" for e in self.geometry.validate()])

        return errors

    def get_summary(self) -> str:
        """
        Get human-readable summary of experiment configuration.

        Returns:
            Formatted string with experiment details
        """
        summary = []
        summary.append("=" * 60)
        summary.append(f"EXPERIMENT: {self.name}")
        summary.append("=" * 60)
        summary.append("")
        summary.append(f"Control: {self.control.name} ({self.control.replicates} replicates)")
        summary.append("")
        summary.append(f"Conditions ({len(self.conditions)}):")
        for condition in self.conditions:
            summary.append(f"  - {condition.name} ({condition.replicates} replicates)")
        summary.append("")
        summary.append(f"Cell Types ({len(self.cell_types)}):")
        for cell_type in self.cell_types:
            summary.append(f"  - {cell_type.name}: {cell_type.barcode}")
        summary.append("")
        summary.append(f"Total Samples: {len(self.get_all_sample_names())}")
        summary.append(f"Total Comparisons: {len(self.get_comparisons())}")
        summary.append("")
        summary.append("Comparisons:")
        for comparison in self.get_comparisons():
            summary.append(f"  - {comparison.name}")
        summary.append("=" * 60)

        return "\n".join(summary)


def create_template_yaml(output_path: str):
    """
    Create a template experiment YAML file.

    Args:
        output_path: Path where template should be created
    """
    template = {
        'experiment': {
            'name': 'Multi_Condition_Screen',
            'control': {
                'name': 'Control',
                'replicates': 3
            },
            'conditions': [
                {'name': 'Condition1', 'replicates': 3},
                {'name': 'Condition2', 'replicates': 3},
                {'name': 'Condition3', 'replicates': 3}
            ],
            'cell_types': [
                {'name': 'Keratinocyte', 'barcode': 'ATGCAGGG'},
                {'name': 'Fibroblast', 'barcode': 'GTTGCAGC'},
                {'name': 'Endothelial', 'barcode': 'ATAGCACG'}
            ],
            'geometry': {
                'barcode_start': 22,
                'barcode_length': 8,
                'grna_start': 12,
                'grna_length': 20,
                'max_barcode_mismatches': 1,
                'allow_grna_mismatch': True,
                'min_qual': 0,
                'grna_anchor': None,
                'anchor_max_offset': 40,
            },
            'library_config': 'default'
        }
    }

    with open(output_path, 'w') as f:
        yaml.dump(template, f, default_flow_style=False, sort_keys=False)

    print(f"✓ Created template configuration: {output_path}")


if __name__ == '__main__':
    """Test experiment configuration."""
    import tempfile

    print("Testing Experiment Configuration Module")
    print("=" * 60)

    # Create temporary template
    with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as f:
        temp_path = f.name
        create_template_yaml(temp_path)

    # Load and validate
    print("\nLoading configuration...")
    config = ExperimentConfig.from_yaml(temp_path)

    # Show summary
    print("\n" + config.get_summary())

    # Validate
    print("\nValidation:")
    errors = config.validate()
    if errors:
        print("  Errors found:")
        for error in errors:
            print(f"    ✗ {error}")
    else:
        print("  ✓ Configuration is valid!")

    # Clean up
    import os
    os.unlink(temp_path)

    print("\n" + "=" * 60)
    print("✓ Experiment configuration module working correctly!")
