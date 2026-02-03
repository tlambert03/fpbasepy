"""Data models for the fluorophore panel optimizer."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np

# Common wavelength grid: 350-850nm at 1nm resolution
WAVELENGTH_GRID = np.arange(350, 851, dtype=float)


@dataclass
class FluorophoreRecord:
    """A fluorophore with interpolated spectra on a common wavelength grid.

    Spectra are interpolated onto WAVELENGTH_GRID (350-850nm, 1nm resolution).
    """

    name: str
    ex_spectrum: np.ndarray  # shape (n_wavelengths,), normalized 0-1
    em_spectrum: np.ndarray  # shape (n_wavelengths,), normalized 0-1
    ext_coeff: float  # M^-1 cm^-1
    qy: float  # 0-1

    @property
    def brightness(self) -> float:
        """Product of extinction coefficient and quantum yield."""
        return self.ext_coeff * self.qy

    @property
    def ex_peak(self) -> float:
        """Wavelength of peak excitation (nm)."""
        return float(WAVELENGTH_GRID[np.argmax(self.ex_spectrum)])

    @property
    def em_peak(self) -> float:
        """Wavelength of peak emission (nm)."""
        return float(WAVELENGTH_GRID[np.argmax(self.em_spectrum)])

    def __repr__(self) -> str:
        return (
            f"FluorophoreRecord({self.name!r}, "
            f"ex_peak={self.ex_peak:.0f}nm, em_peak={self.em_peak:.0f}nm, "
            f"EC={self.ext_coeff:.0f}, QY={self.qy:.2f}, "
            f"brightness={self.brightness:.0f})"
        )


@dataclass
class HardwareConfig:
    """Hardware constraints for the imaging system."""

    detector_mode: Literal["discrete", "spectral"] = "discrete"
    max_lasers_per_round: int = 8
    laser_range: tuple[float, float] = (400.0, 700.0)
    n_detectors: int = 5  # D bands per round (discrete mode)
    n_rounds: int = 1  # R rounds
    min_band_width: float = 10.0  # nm
    band_range: tuple[float, float] = (400.0, 800.0)
    # Spectral detector mode settings
    spectral_start: float = 410.0  # nm
    spectral_bin_width: float = 10.0  # nm
    spectral_n_bins: int = 32

    @property
    def total_channels(self) -> int:
        """Total number of detection channels across all rounds."""
        if self.detector_mode == "discrete":
            return self.n_rounds * self.n_detectors
        return self.n_rounds * self.spectral_n_bins

    @property
    def detectors_per_round(self) -> int:
        if self.detector_mode == "discrete":
            return self.n_detectors
        return self.spectral_n_bins


@dataclass
class RoundConfig:
    """Configuration for one imaging round."""

    laser_wavelengths: np.ndarray  # active laser lines (nm)
    laser_powers: np.ndarray  # relative power (0-1) per line
    band_edges: list[tuple[float, float]] | None = None  # discrete mode only
    target_fluorophores: list[int] = field(default_factory=list)


@dataclass
class PanelResult:
    """One solution from the Pareto front."""

    fluorophore_indices: list[int]  # indices into the original library
    fluorophores: list[FluorophoreRecord]
    round_configs: list[RoundConfig]
    measurement_matrix: np.ndarray  # shape (R*D, N)
    sigma_min: float
    condition_number: float
    brightness_values: np.ndarray
    min_brightness: float
    total_brightness: float

    def summary(self) -> str:
        """Human-readable summary of the panel."""
        lines = []
        names = [fp.name for fp in self.fluorophores]
        lines.append(f"Panel (N={len(self.fluorophores)}): {', '.join(names)}")
        lines.append("")

        for r, rc in enumerate(self.round_configs):
            target_names = [self.fluorophores[i].name for i in rc.target_fluorophores]
            lines.append(f"Round {r + 1} -- targets: {', '.join(target_names)}")
            lasers = ", ".join(f"{lam:.0f}nm" for lam in rc.laser_wavelengths)
            lines.append(f"  Lasers: {lasers}")
            if rc.band_edges:
                bands = ", ".join(
                    f"[{lo:.0f}-{hi:.0f}]" for lo, hi in rc.band_edges
                )
                lines.append(f"  Bands: {bands}")
            lines.append("")

        lines.append(
            f"sigma_min = {self.sigma_min:.4f} | "
            f"cond = {self.condition_number:.2f} | "
            f"brightness range: {self.min_brightness:.1f}-"
            f"{max(self.brightness_values):.1f}"
        )
        return "\n".join(lines)


def make_autofluorescence(
    peak: float = 500.0,
    width: float = 80.0,
    ex_peak: float = 380.0,
    ex_width: float = 50.0,
    ext_coeff: float = 5000.0,
    qy: float = 0.1,
) -> FluorophoreRecord:
    """Create a generic autofluorescence spectrum.

    Parameters
    ----------
    peak : float
        Emission peak wavelength (nm).
    width : float
        Emission Gaussian sigma (nm).
    ex_peak : float
        Excitation peak wavelength (nm).
    ex_width : float
        Excitation Gaussian sigma (nm).
    ext_coeff : float
        Effective extinction coefficient.
    qy : float
        Effective quantum yield.
    """
    em = np.exp(-0.5 * ((WAVELENGTH_GRID - peak) / width) ** 2)
    ex = np.exp(-0.5 * ((WAVELENGTH_GRID - ex_peak) / ex_width) ** 2)
    return FluorophoreRecord(
        name="Autofluorescence",
        ex_spectrum=ex,
        em_spectrum=em,
        ext_coeff=ext_coeff,
        qy=qy,
    )
