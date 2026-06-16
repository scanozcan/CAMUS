"""CAMUS visualization package.

Self-contained, publication-oriented plotting for MAGeCK RRA results plus
QC and reporting. No external CAMUS version is required (this replaces the
old dependency on the v3 ``visualization`` package).
"""

from .volcano import VolcanoPlotter
from .qc import QCPlotter
from .comparison import ComparisonPlotter
from .heatmaps import HeatmapPlotter
from .report import ReportGenerator

__all__ = [
    "VolcanoPlotter",
    "QCPlotter",
    "ComparisonPlotter",
    "HeatmapPlotter",
    "ReportGenerator",
]
