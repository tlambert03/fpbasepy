"""fpicker: Fluorophore Panel Optimizer.

Jointly selects fluorophores AND designs the imaging configuration
(laser lines, emission bands, multi-round assignments) that maximizes
computational unmixability.
"""

from ._library import (
    DEFAULT_FP_NAMES,
    fluorophore_from_arrays,
    load_fluorophore,
    load_library,
)
from ._matrix import (
    build_measurement_matrix,
    compute_condition_number,
    compute_diagnostics,
    compute_sigma_min,
)
from ._models import (
    WAVELENGTH_GRID,
    FluorophoreRecord,
    HardwareConfig,
    PanelResult,
    RoundConfig,
    make_autofluorescence,
)
from ._optimizer import optimize_panel
from ._output import format_results_table, plot_measurement_matrix, plot_panel

__all__ = [
    "DEFAULT_FP_NAMES",
    "FluorophoreRecord",
    "HardwareConfig",
    "PanelResult",
    "RoundConfig",
    "WAVELENGTH_GRID",
    "build_measurement_matrix",
    "compute_condition_number",
    "compute_diagnostics",
    "compute_sigma_min",
    "fluorophore_from_arrays",
    "format_results_table",
    "load_fluorophore",
    "load_library",
    "make_autofluorescence",
    "optimize_panel",
    "plot_measurement_matrix",
    "plot_panel",
]
