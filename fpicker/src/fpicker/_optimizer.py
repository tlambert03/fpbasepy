"""Outer evolutionary optimizer using pymoo NSGA-II."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.repair import Repair
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.optimize import minimize as pymoo_minimize
from pymoo.termination import get_termination

from ._inner import optimize_band_edges
from ._matrix import build_measurement_matrix, compute_sigma_min
from ._models import (
    FluorophoreRecord,
    HardwareConfig,
    PanelResult,
    RoundConfig,
)

if TYPE_CHECKING:
    pass


class IntegerRoundingRepair(Repair):
    """Repair operator that rounds variables to integers and fixes duplicates."""

    def __init__(self, n_select: int, n_library: int, n_rounds: int) -> None:
        super().__init__()
        self.n_select = n_select
        self.n_library = n_library
        self.n_rounds = n_rounds

    def _do(self, problem, X, **kwargs):  # type: ignore[no-untyped-def]
        N = self.n_select

        for i in range(len(X)):
            # Round to integers
            X[i] = np.round(X[i]).astype(int)

            # Clip to bounds
            X[i, :N] = np.clip(X[i, :N], 0, self.n_library - 1)
            X[i, N:] = np.clip(X[i, N:], 0, self.n_rounds - 1)

            # Fix duplicate fluorophore selections
            fp_indices = X[i, :N].copy()
            seen = set()
            available = set(range(self.n_library)) - set(fp_indices)
            for j in range(N):
                if fp_indices[j] in seen:
                    if available:
                        new_idx = available.pop()
                        fp_indices[j] = new_idx
                    else:
                        # Fallback: pick random
                        fp_indices[j] = np.random.randint(0, self.n_library)
                seen.add(fp_indices[j])
            X[i, :N] = fp_indices

        return X


class PanelOptProblem(ElementwiseProblem):
    """pymoo problem for fluorophore panel optimization.

    Decision variables (all integers):
    - x[0:N]: fluorophore indices into library
    - x[N:2N]: round assignments (0..R-1)

    Objectives (to minimize):
    - -sigma_min(M): negative of minimum singular value
    - -min_brightness: negative of worst-case fluorophore brightness
    """

    def __init__(
        self,
        library: list[FluorophoreRecord],
        hardware: HardwareConfig,
        n_select: int,
        autofluorescence: FluorophoreRecord | None = None,
        optimize_bands: bool = True,
        band_opt_max_iter: int = 50,
    ) -> None:
        self.library = library
        self.hardware = hardware
        self.n_select = n_select
        self.autofluorescence = autofluorescence
        self.optimize_bands = optimize_bands
        self.band_opt_max_iter = band_opt_max_iter

        n_vars = 2 * n_select
        xl = np.zeros(n_vars)
        xu = np.array(
            [len(library) - 1] * n_select
            + [max(0, hardware.n_rounds - 1)] * n_select
        )

        super().__init__(
            n_var=n_vars,
            n_obj=2,
            xl=xl,
            xu=xu,
        )

    def _evaluate(self, x, out, *args, **kwargs):  # type: ignore[no-untyped-def]
        N = self.n_select
        fp_indices = np.round(x[:N]).astype(int)
        round_assignments = np.round(x[N:]).astype(int)

        # Clip
        fp_indices = np.clip(fp_indices, 0, len(self.library) - 1)
        round_assignments = np.clip(round_assignments, 0, self.hardware.n_rounds - 1)

        # Get selected fluorophores
        selected = [self.library[i] for i in fp_indices]

        # Add autofluorescence column if provided
        all_fps = list(selected)
        if self.autofluorescence is not None:
            all_fps.append(self.autofluorescence)

        # Build laser configs heuristically
        laser_configs = self._build_laser_configs(selected, round_assignments.tolist())

        if self.optimize_bands and self.hardware.detector_mode == "discrete":
            # Inner optimization of band edges
            round_configs, sigma_min = optimize_band_edges(
                all_fps,
                list(round_assignments)
                + (
                    [0] if self.autofluorescence is not None else []
                ),  # AF in round 0
                laser_configs,
                self.hardware,
                max_iter=self.band_opt_max_iter,
            )
        else:
            # Build configs directly (spectral mode or no band opt)
            round_configs = []
            for r, (lasers, powers) in enumerate(laser_configs):
                targets = [
                    i for i, a in enumerate(round_assignments) if a == r
                ]
                round_configs.append(
                    RoundConfig(
                        laser_wavelengths=lasers,
                        laser_powers=powers,
                        band_edges=None,
                        target_fluorophores=targets,
                    )
                )
            M = build_measurement_matrix(all_fps, round_configs, self.hardware)
            sigma_min = compute_sigma_min(M)

        # Compute brightness objectives on selected (not AF)
        brightness_values = np.array([fp.brightness for fp in selected])
        min_brightness = float(brightness_values.min()) if len(brightness_values) > 0 else 0.0

        out["F"] = np.array([-sigma_min, -min_brightness])

    def _build_laser_configs(
        self,
        fluorophores: list[FluorophoreRecord],
        round_assignments: list[int],
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        """Build laser configurations heuristically.

        For each round, place one laser at the excitation peak of each
        fluorophore assigned to that round. If extra lines are available,
        add lines at excitation peaks of non-assigned fluorophores.
        """
        R = self.hardware.n_rounds
        max_lines = self.hardware.max_lasers_per_round
        lo, hi = self.hardware.laser_range

        configs = []
        for r in range(R):
            # Primary: peaks of assigned fluorophores
            assigned_peaks = [
                fp.ex_peak
                for fp, a in zip(fluorophores, round_assignments)
                if a == r
            ]
            # Secondary: peaks of non-assigned fluorophores (permanent labels)
            other_peaks = [
                fp.ex_peak
                for fp, a in zip(fluorophores, round_assignments)
                if a != r
            ]

            all_peaks = assigned_peaks + other_peaks
            # Deduplicate (within 2nm) and limit
            unique_peaks = _deduplicate_peaks(all_peaks, min_separation=2.0)
            # Clip to laser range
            unique_peaks = [p for p in unique_peaks if lo <= p <= hi]
            unique_peaks = unique_peaks[:max_lines]

            if not unique_peaks:
                unique_peaks = [550.0]  # fallback

            wavelengths = np.array(unique_peaks)
            powers = np.ones(len(wavelengths))  # uniform power

            configs.append((wavelengths, powers))

        return configs


def _deduplicate_peaks(
    peaks: list[float], min_separation: float = 2.0
) -> list[float]:
    """Remove near-duplicate peaks, keeping the first occurrence."""
    if not peaks:
        return []
    result = [peaks[0]]
    for p in peaks[1:]:
        if all(abs(p - r) >= min_separation for r in result):
            result.append(p)
    return result


def optimize_panel(
    library: list[FluorophoreRecord],
    n_select: int,
    hardware: HardwareConfig | None = None,
    autofluorescence: FluorophoreRecord | None = None,
    pop_size: int = 100,
    n_generations: int = 200,
    optimize_bands: bool = True,
    band_opt_max_iter: int = 50,
    seed: int | None = None,
    verbose: bool = True,
) -> list[PanelResult]:
    """Run the fluorophore panel optimizer.

    Parameters
    ----------
    library : list of FluorophoreRecord
        Candidate fluorophore library.
    n_select : int
        Number of fluorophores to select for the panel.
    hardware : HardwareConfig, optional
        Hardware constraints. Defaults to a Leica Stellaris-like config.
    autofluorescence : FluorophoreRecord, optional
        Autofluorescence spectrum to include as extra column in M.
    pop_size : int
        Population size for NSGA-II.
    n_generations : int
        Number of generations.
    optimize_bands : bool
        Whether to run inner band-edge optimization (discrete mode).
    band_opt_max_iter : int
        Max iterations for inner scipy optimizer.
    seed : int, optional
        Random seed for reproducibility.
    verbose : bool
        Print progress information.

    Returns
    -------
    list of PanelResult
        Pareto-optimal solutions, sorted by sigma_min (best first).
    """
    if hardware is None:
        hardware = HardwareConfig()

    if n_select > len(library):
        raise ValueError(
            f"Cannot select {n_select} fluorophores from library of {len(library)}"
        )

    if n_select < 2:
        raise ValueError("Must select at least 2 fluorophores")

    problem = PanelOptProblem(
        library=library,
        hardware=hardware,
        n_select=n_select,
        autofluorescence=autofluorescence,
        optimize_bands=optimize_bands,
        band_opt_max_iter=band_opt_max_iter,
    )

    repair = IntegerRoundingRepair(n_select, len(library), hardware.n_rounds)

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=IntegerRandomSampling(),
        crossover=SBX(prob=0.9, eta=3.0, vtype=float, repair=repair),
        mutation=PM(eta=3.0, vtype=float, repair=repair),
        repair=repair,
    )

    termination = get_termination("n_gen", n_generations)

    if verbose:
        print(
            f"Optimizing: selecting {n_select} from {len(library)} fluorophores, "
            f"{hardware.n_rounds} round(s), {hardware.detector_mode} detector"
        )

    result = pymoo_minimize(
        problem,
        algorithm,
        termination,
        seed=seed,
        verbose=verbose,
    )

    if result.F is None or result.X is None:
        return []

    # Build PanelResult objects from Pareto front
    results = []
    X = result.X if result.X.ndim == 2 else result.X.reshape(1, -1)
    F = result.F if result.F.ndim == 2 else result.F.reshape(1, -1)

    for x_i, f_i in zip(X, F):
        fp_indices = np.round(x_i[: n_select]).astype(int)
        round_assignments = np.round(x_i[n_select:]).astype(int)
        fp_indices = np.clip(fp_indices, 0, len(library) - 1)
        round_assignments = np.clip(round_assignments, 0, hardware.n_rounds - 1)

        selected = [library[i] for i in fp_indices]
        all_fps = list(selected)
        if autofluorescence is not None:
            all_fps.append(autofluorescence)

        # Rebuild configs
        laser_configs = problem._build_laser_configs(
            selected, round_assignments.tolist()
        )

        if optimize_bands and hardware.detector_mode == "discrete":
            af_assignments = list(round_assignments) + (
                [0] if autofluorescence is not None else []
            )
            round_configs, sigma_min = optimize_band_edges(
                all_fps,
                af_assignments,
                laser_configs,
                hardware,
                max_iter=band_opt_max_iter,
            )
        else:
            round_configs = []
            for r, (lasers, powers) in enumerate(laser_configs):
                targets = [
                    i for i, a in enumerate(round_assignments) if a == r
                ]
                round_configs.append(
                    RoundConfig(
                        laser_wavelengths=lasers,
                        laser_powers=powers,
                        band_edges=None,
                        target_fluorophores=targets,
                    )
                )

        M = build_measurement_matrix(all_fps, round_configs, hardware)
        sigma_min = compute_sigma_min(M)
        sv = np.linalg.svd(M, compute_uv=False)
        cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf")

        brightness_values = np.array([fp.brightness for fp in selected])

        results.append(
            PanelResult(
                fluorophore_indices=fp_indices.tolist(),
                fluorophores=selected,
                round_configs=round_configs,
                measurement_matrix=M,
                sigma_min=sigma_min,
                condition_number=cond,
                brightness_values=brightness_values,
                min_brightness=float(brightness_values.min()),
                total_brightness=float(brightness_values.sum()),
            )
        )

    # Sort by sigma_min descending (best first)
    results.sort(key=lambda r: r.sigma_min, reverse=True)
    return results
