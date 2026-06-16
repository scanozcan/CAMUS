"""Shared helpers for CAMUS visualization (non-interactive backend, colors)."""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless / pipeline-safe
import matplotlib.pyplot as plt  # noqa: E402

# Consistent color scheme
DEPLETED = "#2166ac"   # blue
ENRICHED = "#b2182b"   # red
NEUTRAL = "#cccccc"    # gray


class BasePlotter:
    def __init__(self, output_dir, fdr_threshold: float = 0.05):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.fdr_threshold = fdr_threshold

    def _save(self, fig, filename: str):
        path = self.output_dir / filename
        fig.savefig(path, dpi=200, bbox_inches="tight")
        plt.close(fig)
        return path
