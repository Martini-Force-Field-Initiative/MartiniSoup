"""
binding_analysis
================
Tools for analyzing metabolite–protein residence times extracted from
molecular dynamics simulations.

This package provides:
    • Efficient storage of binding events per atom and molecule
    • Aggregation tools for molecule and molecule-type statistics
    • Histogram, survival curve, hazard function, and bootstrap analyses
    • LMfit-based kinetic model fitting (exponential, Weibull, etc.)
"""

from .data_structures import BoundState, MetaboliteResidences
from .analysis import ResidenceAnalysis
from .models import KineticModels
from . import utils

__all__ = [
    "BoundState",
    "MetaboliteResidences",
    "ResidenceAnalysis",
    "KineticModels",
    "utils",
]