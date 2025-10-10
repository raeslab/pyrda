"""
pyrda: Python implementation of Distance-Based Redundancy Analysis

A Python package for performing distance-based redundancy analysis (db-RDA)
as practiced in the RaesLab, providing a streamlined workflow for identifying
and visualizing the impact of metadata variables on community composition.
"""

from .core import get_cumulative, get_r2, get_r2_single, prepare_plot_data, rldbrda
from .plotting import plot_dbrda
from .validation import validate_input

__version__ = "0.1.0"
__author__ = "Claude Code"
__email__ = "noreply@anthropic.com"

__all__ = [
    "rldbrda",
    "get_r2",
    "get_r2_single",
    "get_cumulative",
    "prepare_plot_data",
    "plot_dbrda",
    "validate_input",
]
