"""
Analysis submodule for binding_analysis.

Provides tools for:
- Survival curve analysis
- Histogram analysis
- Exponential kinetic fitting
- High-level unified interface (ResidenceAnalysis)
"""

from .survival import SurvivalAnalysis
from .histogram import HistogramAnalysis
from .models import SingleExponentialModel
from .analysis import ResidenceAnalysis

__all__ = [
    "SurvivalAnalysis",
    "HistogramAnalysis",
    "SingleExponentialModel",
    "ResidenceAnalysis",
]
