#!/usr/bin/env python3
"""
Library configuration for CAMUS.

Defines library sizes, gene categories, and sgRNA specifications.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class LibraryConfig:
    """Configuration for CRISPR library design."""

    # Library size parameters
    n_cell_specific_genes: int = 400  # Genes unique to each cell type
    n_common_essential: int = 25  # Genes shared across all cell types
    n_controls: int = 500  # Non-targeting controls
    sgrnas_per_gene: int = 6  # sgRNAs targeting each gene

    # Gene category breakdown (must sum to n_cell_specific_genes)
    n_condition1_specific: int = 30  # Genes depleted only in Condition1
    n_condition2_specific: int = 30  # Genes depleted only in Condition2
    n_celltype_essential: int = 50  # Genes depleted in both conditions
    n_resistance_per_condition: int = 10  # Genes enriched per condition (resistance)
    # Remaining genes are non-essential: 400 - 30 - 30 - 50 - (10*N) = variable

    # Cell types
    cell_types: List[str] = None

    def __post_init__(self):
        """Validate configuration after initialization."""
        if self.cell_types is None:
            self.cell_types = ['Keratinocyte', 'Fibroblast', 'Endothelial']

        # Validate gene counts
        total_categorized = (
            self.n_condition1_specific +
            self.n_condition2_specific +
            self.n_celltype_essential
        )

        if total_categorized > self.n_cell_specific_genes:
            raise ValueError(
                f"Categorized genes ({total_categorized}) exceed "
                f"total cell-specific genes ({self.n_cell_specific_genes})"
            )

    @property
    def n_non_essential(self) -> int:
        """Calculate number of non-essential genes."""
        return (
            self.n_cell_specific_genes -
            self.n_condition1_specific -
            self.n_condition2_specific -
            self.n_celltype_essential
        )

    @property
    def total_sgrnas_per_celltype(self) -> int:
        """Total sgRNAs in each cell type's library."""
        return (
            self.n_cell_specific_genes * self.sgrnas_per_gene +  # Cell-specific
            self.n_common_essential * self.sgrnas_per_gene +     # Common essential
            self.n_controls                                       # Controls
        )

    @property
    def total_sgrnas_combined(self) -> int:
        """Total sgRNAs across all cell types (with overlap for common)."""
        cell_specific_total = (
            self.n_cell_specific_genes * self.sgrnas_per_gene * len(self.cell_types)
        )
        common_total = self.n_common_essential * self.sgrnas_per_gene
        controls_total = self.n_controls * len(self.cell_types)

        return cell_specific_total + common_total + controls_total

    def get_library_summary(self) -> str:
        """
        Get human-readable summary of library configuration.

        Returns:
            Formatted string with library details
        """
        summary = []
        summary.append("=" * 60)
        summary.append("LIBRARY CONFIGURATION")
        summary.append("=" * 60)
        summary.append("")
        summary.append(f"Cell Types: {', '.join(self.cell_types)}")
        summary.append(f"sgRNAs per Gene: {self.sgrnas_per_gene}")
        summary.append("")
        summary.append("Per Cell Type Library:")
        summary.append(f"  Cell-specific genes:       {self.n_cell_specific_genes} genes")
        summary.append(f"    ├─ Condition1-specific:   {self.n_condition1_specific} genes")
        summary.append(f"    ├─ Condition2-specific:   {self.n_condition2_specific} genes")
        summary.append(f"    ├─ Cell-type essential:   {self.n_celltype_essential} genes")
        summary.append(f"    └─ Non-essential:         {self.n_non_essential} genes")
        summary.append(f"  Common essential genes:    {self.n_common_essential} genes (shared)")
        summary.append(f"  Non-targeting controls:    {self.n_controls} sgRNAs")
        summary.append("")
        summary.append("sgRNA Counts:")
        summary.append(f"  Per cell type:             {self.total_sgrnas_per_celltype:,} sgRNAs")
        summary.append(f"  Combined total:            {self.total_sgrnas_combined:,} sgRNAs")
        summary.append("")
        summary.append("=" * 60)

        return "\n".join(summary)

    def validate(self) -> List[str]:
        """
        Validate library configuration.

        Returns:
            List of validation errors (empty if valid)
        """
        errors = []

        # Check gene counts
        if self.n_cell_specific_genes < 1:
            errors.append("n_cell_specific_genes must be >= 1")

        if self.n_common_essential < 1:
            errors.append("n_common_essential must be >= 1")

        if self.n_controls < 1:
            errors.append("n_controls must be >= 1")

        if self.sgrnas_per_gene < 1:
            errors.append("sgrnas_per_gene must be >= 1")

        # Check category breakdown
        total_categorized = (
            self.n_condition1_specific +
            self.n_condition2_specific +
            self.n_celltype_essential
        )

        if total_categorized > self.n_cell_specific_genes:
            errors.append(
                f"Sum of categorized genes ({total_categorized}) exceeds "
                f"total ({self.n_cell_specific_genes})"
            )

        if self.n_non_essential < 0:
            errors.append(f"Non-essential count is negative: {self.n_non_essential}")

        # Check cell types
        if not self.cell_types or len(self.cell_types) < 1:
            errors.append("At least one cell type must be specified")

        return errors


# Predefined library configurations

DEFAULT_CONFIG = LibraryConfig(
    n_cell_specific_genes=400,
    n_common_essential=25,
    n_controls=500,
    sgrnas_per_gene=6,
    n_condition1_specific=30,
    n_condition2_specific=30,
    n_celltype_essential=50
)

SMALL_TEST_CONFIG = LibraryConfig(
    n_cell_specific_genes=120,
    n_common_essential=10,
    n_controls=100,
    sgrnas_per_gene=10,
    n_condition1_specific=50,
    n_condition2_specific=50,
    n_celltype_essential=10
)

LARGE_CONFIG = LibraryConfig(
    n_cell_specific_genes=1000,
    n_common_essential=50,
    n_controls=1000,
    sgrnas_per_gene=6,
    n_condition1_specific=100,
    n_condition2_specific=100,
    n_celltype_essential=100
)


def get_config(name: str = 'default') -> LibraryConfig:
    """
    Get predefined library configuration by name.

    Args:
        name: Configuration name ('default', 'small', 'large')

    Returns:
        LibraryConfig instance

    Raises:
        ValueError: If configuration name not found

    Examples:
        >>> config = get_config('default')
        >>> print(config.total_sgrnas_per_celltype)
        3050
    """
    configs = {
        'default': DEFAULT_CONFIG,
        'small': SMALL_TEST_CONFIG,
        'large': LARGE_CONFIG
    }

    if name.lower() not in configs:
        raise ValueError(
            f"Unknown configuration: {name}. "
            f"Available: {', '.join(configs.keys())}"
        )

    return configs[name.lower()]


if __name__ == '__main__':
    """Test library configuration."""
    print("Testing Library Configuration Module")
    print("")

    # Test default config
    config = get_config('default')
    print(config.get_library_summary())

    # Validate
    errors = config.validate()
    if errors:
        print("\nValidation Errors:")
        for error in errors:
            print(f"  ✗ {error}")
    else:
        print("\n✓ Configuration is valid!")

    # Show all presets
    print("\n" + "=" * 60)
    print("Available Presets:")
    print("=" * 60)
    for name in ['default', 'small', 'large']:
        cfg = get_config(name)
        print(f"\n{name.upper()}:")
        print(f"  Total sgRNAs per cell type: {cfg.total_sgrnas_per_celltype:,}")
        print(f"  Combined total: {cfg.total_sgrnas_combined:,}")
