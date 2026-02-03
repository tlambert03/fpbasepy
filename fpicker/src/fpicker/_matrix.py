"""Measurement matrix construction."""

from __future__ import annotations

import numpy as np

from ._models import WAVELENGTH_GRID, FluorophoreRecord, HardwareConfig, RoundConfig


def _excitation_factor(
    fp: FluorophoreRecord,
    laser_wavelengths: np.ndarray,
    laser_powers: np.ndarray,
) -> float:
    """Compute E(r,f) = sum_l p_l * ex_spectrum(lambda_l).

    Parameters
    ----------
    fp : FluorophoreRecord
        The fluorophore.
    laser_wavelengths : np.ndarray
        Active laser wavelengths (nm).
    laser_powers : np.ndarray
        Relative laser powers (0-1).

    Returns
    -------
    float
        Total excitation factor.
    """
    total = 0.0
    for lam, power in zip(laser_wavelengths, laser_powers):
        # Find closest wavelength index on the grid
        idx = int(np.argmin(np.abs(WAVELENGTH_GRID - lam)))
        total += power * fp.ex_spectrum[idx]
    return total


def _detection_factor(
    fp: FluorophoreRecord,
    lo: float,
    hi: float,
) -> float:
    """Compute D(r,d,f) = fraction of emission in detector band [lo, hi].

    Parameters
    ----------
    fp : FluorophoreRecord
        The fluorophore.
    lo, hi : float
        Band edges in nm.

    Returns
    -------
    float
        Fraction of total emission captured in the band.
    """
    mask = (WAVELENGTH_GRID >= lo) & (WAVELENGTH_GRID <= hi)
    if not mask.any():
        return 0.0
    band_integral = float(np.trapezoid(fp.em_spectrum[mask], WAVELENGTH_GRID[mask]))
    total_integral = float(np.trapezoid(fp.em_spectrum, WAVELENGTH_GRID))
    if total_integral <= 0:
        return 0.0
    return band_integral / total_integral


def build_measurement_matrix(
    fluorophores: list[FluorophoreRecord],
    round_configs: list[RoundConfig],
    hardware: HardwareConfig,
) -> np.ndarray:
    """Build the measurement matrix M.

    Shape: (R*D, N) where R=rounds, D=detectors per round, N=fluorophores.

    M[(r,d), f] = E(r,f) * D(r,d,f) * brightness_f

    Parameters
    ----------
    fluorophores : list of FluorophoreRecord
        Selected fluorophores.
    round_configs : list of RoundConfig
        Per-round laser/band configurations.
    hardware : HardwareConfig
        Hardware constraints.

    Returns
    -------
    np.ndarray
        Measurement matrix of shape (R*D, N).
    """
    N = len(fluorophores)
    D = hardware.detectors_per_round
    R = len(round_configs)
    M = np.zeros((R * D, N))

    for r, rc in enumerate(round_configs):
        for f_idx, fp in enumerate(fluorophores):
            E_rf = _excitation_factor(fp, rc.laser_wavelengths, rc.laser_powers)
            brightness = fp.brightness

            if hardware.detector_mode == "discrete":
                if rc.band_edges is None:
                    continue
                for d, (lo, hi) in enumerate(rc.band_edges):
                    D_rdf = _detection_factor(fp, lo, hi)
                    M[r * D + d, f_idx] = E_rf * D_rdf * brightness
            else:
                # Spectral detector: fixed hardware bins
                for d in range(hardware.spectral_n_bins):
                    lo = hardware.spectral_start + d * hardware.spectral_bin_width
                    hi = lo + hardware.spectral_bin_width
                    D_rdf = _detection_factor(fp, lo, hi)
                    M[r * D + d, f_idx] = E_rf * D_rdf * brightness

    return M


def compute_sigma_min(M: np.ndarray) -> float:
    """Compute the minimum singular value of M.

    Returns 0.0 if M is degenerate (all zeros, or fewer rows than columns).
    """
    if M.size == 0 or M.shape[0] < M.shape[1]:
        return 0.0
    # Normalize columns by max brightness to make sigma_min scale-invariant
    col_norms = np.linalg.norm(M, axis=0)
    if np.any(col_norms == 0):
        return 0.0
    sv = np.linalg.svd(M, compute_uv=False)
    return float(sv[-1])


def compute_condition_number(M: np.ndarray) -> float:
    """Compute the condition number kappa(M) = sigma_max / sigma_min."""
    if M.size == 0 or M.shape[0] < M.shape[1]:
        return float("inf")
    sv = np.linalg.svd(M, compute_uv=False)
    if sv[-1] == 0:
        return float("inf")
    return float(sv[0] / sv[-1])


def compute_diagnostics(
    M: np.ndarray, fluorophores: list[FluorophoreRecord]
) -> dict:
    """Compute diagnostic metrics for a measurement matrix.

    Returns
    -------
    dict with keys:
        sigma_min, sigma_max, condition_number,
        per_fluorophore_signal, per_fluorophore_crosstalk_ratio,
        pairwise_crosstalk
    """
    sv = np.linalg.svd(M, compute_uv=False)
    sigma_min = float(sv[-1]) if len(sv) > 0 else 0.0
    sigma_max = float(sv[0]) if len(sv) > 0 else 0.0
    cond = sigma_max / sigma_min if sigma_min > 0 else float("inf")

    N = M.shape[1]

    # Per-fluorophore signal: L2 norm of column
    signals = np.linalg.norm(M, axis=0)

    # Pairwise crosstalk: cosine similarity between column pairs
    pairwise = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            ni = np.linalg.norm(M[:, i])
            nj = np.linalg.norm(M[:, j])
            if ni > 0 and nj > 0:
                pairwise[i, j] = float(
                    np.dot(M[:, i], M[:, j]) / (ni * nj)
                )

    # Per-fluorophore signal-to-crosstalk ratio
    # max off-diagonal cosine similarity for each fluorophore
    crosstalk_ratios = np.zeros(N)
    for i in range(N):
        off_diag = [pairwise[i, j] for j in range(N) if j != i]
        max_crosstalk = max(off_diag) if off_diag else 0.0
        crosstalk_ratios[i] = (
            (1.0 - max_crosstalk) / max_crosstalk if max_crosstalk > 0 else float("inf")
        )

    return {
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "condition_number": cond,
        "per_fluorophore_signal": signals,
        "per_fluorophore_crosstalk_ratio": crosstalk_ratios,
        "pairwise_crosstalk": pairwise,
        "fluorophore_names": [fp.name for fp in fluorophores],
    }
