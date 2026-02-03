"""Inner optimizer: band-edge optimization via scipy."""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize

from ._matrix import build_measurement_matrix, compute_sigma_min
from ._models import FluorophoreRecord, HardwareConfig, RoundConfig


def _pack_band_edges(round_configs: list[RoundConfig], D: int) -> np.ndarray:
    """Pack band edges from round configs into a flat vector."""
    x = []
    for rc in round_configs:
        if rc.band_edges is None:
            # Initialize with evenly spaced bands
            for d in range(D):
                x.extend([400.0 + d * 80.0, 400.0 + d * 80.0 + 60.0])
        else:
            for lo, hi in rc.band_edges:
                x.extend([lo, hi])
    return np.array(x)


def _unpack_band_edges(
    x: np.ndarray, R: int, D: int
) -> list[list[tuple[float, float]]]:
    """Unpack flat vector into per-round band edges."""
    bands_per_round = []
    idx = 0
    for _r in range(R):
        bands = []
        for _d in range(D):
            lo, hi = x[idx], x[idx + 1]
            bands.append((float(lo), float(hi)))
            idx += 2
        bands_per_round.append(bands)
    return bands_per_round


def _initial_band_edges(
    fluorophores: list[FluorophoreRecord],
    round_assignments: list[int],
    hardware: HardwareConfig,
) -> list[list[tuple[float, float]]]:
    """Generate initial band edges based on fluorophore emission peaks.

    Places bands centered on emission peaks of the assigned fluorophores,
    then fills remaining bands evenly in gaps.
    """
    R = hardware.n_rounds
    D = hardware.n_detectors
    lo_bound, hi_bound = hardware.band_range
    min_w = hardware.min_band_width

    all_bands: list[list[tuple[float, float]]] = []

    for r in range(R):
        # Get emission peaks of fluorophores assigned to this round
        peaks = sorted(
            fp.em_peak
            for fp, assign in zip(fluorophores, round_assignments)
            if assign == r
        )

        # Also include emission peaks of ALL fluorophores (permanent labels)
        all_peaks = sorted(fp.em_peak for fp in fluorophores)

        # Use all peaks, prioritizing assigned ones
        target_peaks = peaks + [p for p in all_peaks if p not in peaks]
        target_peaks = target_peaks[:D]

        # If we still don't have enough, fill evenly
        while len(target_peaks) < D:
            target_peaks.append(
                lo_bound
                + (hi_bound - lo_bound) * (len(target_peaks) + 0.5) / D
            )
        target_peaks = sorted(target_peaks)

        # Create bands centered on peaks, width ~40nm
        bands = []
        default_half_width = 20.0
        for peak in target_peaks:
            lo = max(lo_bound, peak - default_half_width)
            hi = min(hi_bound, peak + default_half_width)
            if hi - lo < min_w:
                hi = lo + min_w
            bands.append((lo, hi))

        # Resolve overlaps by pushing bands apart
        bands = _resolve_overlaps(bands, lo_bound, hi_bound, min_w)
        all_bands.append(bands)

    return all_bands


def _resolve_overlaps(
    bands: list[tuple[float, float]],
    lo_bound: float,
    hi_bound: float,
    min_width: float,
) -> list[tuple[float, float]]:
    """Resolve overlapping bands by pushing them apart."""
    bands = sorted(bands)
    resolved = list(bands)

    for _ in range(20):  # iterate to convergence
        changed = False
        for i in range(len(resolved) - 1):
            lo_i, hi_i = resolved[i]
            lo_j, hi_j = resolved[i + 1]
            if hi_i > lo_j:
                # Overlap: split the difference
                mid = (hi_i + lo_j) / 2
                resolved[i] = (lo_i, mid - 0.5)
                resolved[i + 1] = (mid + 0.5, hi_j)
                changed = True
        if not changed:
            break

    # Enforce min width and bounds
    final = []
    for lo, hi in resolved:
        lo = max(lo_bound, lo)
        hi = min(hi_bound, hi)
        if hi - lo < min_width:
            hi = lo + min_width
        if hi > hi_bound:
            hi = hi_bound
            lo = max(lo_bound, hi - min_width)
        final.append((lo, hi))

    return final


def optimize_band_edges(
    fluorophores: list[FluorophoreRecord],
    round_assignments: list[int],
    laser_configs: list[tuple[np.ndarray, np.ndarray]],
    hardware: HardwareConfig,
    max_iter: int = 100,
) -> tuple[list[RoundConfig], float]:
    """Optimize detector band edges to maximize sigma_min(M).

    Given fixed fluorophores, round assignments, and laser lines, this
    optimizes the continuous band-edge variables using scipy SLSQP.

    Parameters
    ----------
    fluorophores : list of FluorophoreRecord
        Selected fluorophores.
    round_assignments : list of int
        Round assignment for each fluorophore (0..R-1).
    laser_configs : list of (wavelengths, powers)
        Per-round laser configurations.
    hardware : HardwareConfig
        Hardware constraints.
    max_iter : int
        Maximum iterations for scipy optimizer.

    Returns
    -------
    round_configs : list of RoundConfig
        Optimized per-round configurations.
    sigma_min : float
        Achieved minimum singular value.
    """
    if hardware.detector_mode == "spectral":
        # No band optimization needed for spectral detector
        round_configs = []
        for r, (lasers, powers) in enumerate(laser_configs):
            targets = [i for i, a in enumerate(round_assignments) if a == r]
            round_configs.append(
                RoundConfig(
                    laser_wavelengths=lasers,
                    laser_powers=powers,
                    band_edges=None,
                    target_fluorophores=targets,
                )
            )
        M = build_measurement_matrix(fluorophores, round_configs, hardware)
        return round_configs, compute_sigma_min(M)

    R = hardware.n_rounds
    D = hardware.n_detectors
    lo_bound, hi_bound = hardware.band_range
    min_w = hardware.min_band_width

    # Initialize band edges
    init_bands = _initial_band_edges(fluorophores, round_assignments, hardware)

    # Build initial round configs to get initial x0
    init_round_configs = []
    for r in range(R):
        targets = [i for i, a in enumerate(round_assignments) if a == r]
        lasers, powers = laser_configs[r]
        init_round_configs.append(
            RoundConfig(
                laser_wavelengths=lasers,
                laser_powers=powers,
                band_edges=init_bands[r],
                target_fluorophores=targets,
            )
        )

    x0 = _pack_band_edges(init_round_configs, D)

    # Bounds
    bounds = []
    for _r in range(R):
        for _d in range(D):
            bounds.append((lo_bound, hi_bound - min_w))  # lo
            bounds.append((lo_bound + min_w, hi_bound))  # hi

    def objective(x: np.ndarray) -> float:
        bands_per_round = _unpack_band_edges(x, R, D)
        configs = []
        for r in range(R):
            targets = [i for i, a in enumerate(round_assignments) if a == r]
            lasers, powers = laser_configs[r]
            configs.append(
                RoundConfig(
                    laser_wavelengths=lasers,
                    laser_powers=powers,
                    band_edges=bands_per_round[r],
                    target_fluorophores=targets,
                )
            )
        M = build_measurement_matrix(fluorophores, configs, hardware)
        return -compute_sigma_min(M)  # minimize negative sigma_min

    # Constraints: lo + min_width <= hi, and non-overlapping within round
    constraints = []
    for r in range(R):
        for d in range(D):
            base = r * D * 2 + d * 2
            # hi - lo >= min_width
            constraints.append(
                {
                    "type": "ineq",
                    "fun": lambda x, b=base, mw=min_w: x[b + 1] - x[b] - mw,
                }
            )
            # Non-overlapping: hi_d <= lo_{d+1}
            if d < D - 1:
                next_base = base + 2
                constraints.append(
                    {
                        "type": "ineq",
                        "fun": lambda x, b=base, nb=next_base: x[nb] - x[b + 1],
                    }
                )

    result = minimize(
        objective,
        x0,
        method="SLSQP",
        bounds=bounds,
        constraints=constraints,
        options={"maxiter": max_iter, "ftol": 1e-8},
    )

    # Build final configs from optimized x
    bands_per_round = _unpack_band_edges(result.x, R, D)
    final_configs = []
    for r in range(R):
        targets = [i for i, a in enumerate(round_assignments) if a == r]
        lasers, powers = laser_configs[r]
        final_configs.append(
            RoundConfig(
                laser_wavelengths=lasers,
                laser_powers=powers,
                band_edges=bands_per_round[r],
                target_fluorophores=targets,
            )
        )

    M = build_measurement_matrix(fluorophores, final_configs, hardware)
    sigma_min = compute_sigma_min(M)

    return final_configs, sigma_min
