"""Load fluorophore library from FPbase via fpbasepy."""

from __future__ import annotations

from typing import Sequence

import numpy as np

from ._models import WAVELENGTH_GRID, FluorophoreRecord


def _interpolate_spectrum(
    data: list[tuple[float, float]], grid: np.ndarray
) -> np.ndarray:
    """Interpolate spectrum data onto a common wavelength grid.

    Parameters
    ----------
    data : list of (wavelength, intensity) tuples
        Raw spectrum data from fpbasepy.
    grid : np.ndarray
        Target wavelength grid.

    Returns
    -------
    np.ndarray
        Interpolated spectrum values, clipped to [0, 1].
    """
    if not data:
        return np.zeros_like(grid)
    wavelengths, intensities = zip(*sorted(data))
    interp = np.interp(grid, wavelengths, intensities, left=0.0, right=0.0)
    return np.clip(interp, 0.0, 1.0)


def load_fluorophore(name: str) -> FluorophoreRecord | None:
    """Load a single fluorophore from FPbase.

    Returns None if the fluorophore lacks required data (spectra, EC, QY).
    """
    import fpbase

    try:
        fp = fpbase.get_fluorophore(name)
    except (ValueError, Exception):
        return None

    state = fp.default_state
    if state is None:
        return None

    if state.ext_coeff is None or state.qy is None:
        return None

    ex_spec = state.excitation_spectrum
    em_spec = state.emission_spectrum
    if ex_spec is None or em_spec is None:
        return None

    ex_data = _interpolate_spectrum(ex_spec.data, WAVELENGTH_GRID)
    em_data = _interpolate_spectrum(em_spec.data, WAVELENGTH_GRID)

    # Skip if spectra are essentially empty
    if ex_data.max() < 0.01 or em_data.max() < 0.01:
        return None

    return FluorophoreRecord(
        name=fp.name,
        ex_spectrum=ex_data,
        em_spectrum=em_data,
        ext_coeff=state.ext_coeff,
        qy=state.qy,
    )


def load_library(
    names: Sequence[str] | None = None,
    *,
    min_brightness: float = 0.0,
    require_spectra: bool = True,
) -> list[FluorophoreRecord]:
    """Load a library of fluorophores from FPbase.

    Parameters
    ----------
    names : sequence of str, optional
        Specific fluorophore names to load. If None, loads common FPs.
    min_brightness : float
        Minimum brightness (EC * QY) threshold.
    require_spectra : bool
        If True, skip fluorophores without full spectra.

    Returns
    -------
    list of FluorophoreRecord
        Loaded fluorophores, sorted by emission peak.
    """
    if names is None:
        names = DEFAULT_FP_NAMES

    library = []
    for name in names:
        rec = load_fluorophore(name)
        if rec is None:
            continue
        if rec.brightness < min_brightness:
            continue
        library.append(rec)

    # Sort by emission peak
    library.sort(key=lambda fp: fp.em_peak)
    return library


def fluorophore_from_arrays(
    name: str,
    ex_wavelengths: np.ndarray,
    ex_intensities: np.ndarray,
    em_wavelengths: np.ndarray,
    em_intensities: np.ndarray,
    ext_coeff: float,
    qy: float,
) -> FluorophoreRecord:
    """Create a FluorophoreRecord from raw arrays.

    Useful for custom dyes or measured spectra not in FPbase.
    """
    ex_data = list(zip(ex_wavelengths, ex_intensities))
    em_data = list(zip(em_wavelengths, em_intensities))
    return FluorophoreRecord(
        name=name,
        ex_spectrum=_interpolate_spectrum(ex_data, WAVELENGTH_GRID),
        em_spectrum=_interpolate_spectrum(em_data, WAVELENGTH_GRID),
        ext_coeff=ext_coeff,
        qy=qy,
    )


# A curated set of commonly used fluorescent proteins spanning the visible spectrum
DEFAULT_FP_NAMES = [
    # Blue
    "mTagBFP2",
    "EBFP2",
    # Cyan
    "mTurquoise2",
    "mCerulean3",
    # Green
    "EGFP",
    "mNeonGreen",
    "mClover3",
    # Yellow
    "mVenus",
    "EYFP",
    "mCitrine",
    # Orange
    "mKO2",
    "mOrange2",
    "CyOFP1",
    # Red
    "mScarlet",
    "mCherry",
    "tdTomato",
    "mRuby3",
    # Far-red
    "mMaroon1",
    "mCardinal",
    "miRFP670nano",
    "miRFP680",
    "miRFP713",
]
