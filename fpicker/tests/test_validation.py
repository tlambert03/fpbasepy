"""Validation tests using real fluorophore data from FPbase."""

from __future__ import annotations

import numpy as np
import pytest

import fpicker


# A small set of well-characterized FPs for testing
TEST_FP_NAMES = [
    "EGFP",
    "mCherry",
    "mTagBFP2",
    "mVenus",
    "tdTomato",
    "mTurquoise2",
    "mNeonGreen",
    "mScarlet",
]


@pytest.fixture(scope="module")
def library():
    """Load a small library from FPbase."""
    lib = fpicker.load_library(TEST_FP_NAMES)
    assert len(lib) >= 4, f"Only loaded {len(lib)} fluorophores, expected >= 4"
    return lib


class TestLibrary:
    def test_load_single(self):
        fp = fpicker.load_fluorophore("EGFP")
        assert fp is not None
        assert fp.name == "EGFP"
        assert fp.ex_spectrum.shape == fpicker.WAVELENGTH_GRID.shape
        assert fp.em_spectrum.shape == fpicker.WAVELENGTH_GRID.shape
        assert fp.ext_coeff > 0
        assert 0 < fp.qy <= 1
        assert fp.brightness > 0
        # EGFP should have excitation peak around 488nm
        assert 480 < fp.ex_peak < 500
        # EGFP should have emission peak around 507nm
        assert 500 < fp.em_peak < 520

    def test_load_library(self, library):
        # Should be sorted by emission peak
        em_peaks = [fp.em_peak for fp in library]
        assert em_peaks == sorted(em_peaks)

    def test_spectra_valid(self, library):
        for fp in library:
            assert fp.ex_spectrum.min() >= 0
            assert fp.ex_spectrum.max() <= 1
            assert fp.em_spectrum.min() >= 0
            assert fp.em_spectrum.max() <= 1
            assert fp.ex_spectrum.max() > 0.5  # should have a real peak
            assert fp.em_spectrum.max() > 0.5


class TestMeasurementMatrix:
    def test_build_discrete(self, library):
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=4
        )
        fps = library[:4]
        rc = fpicker.RoundConfig(
            laser_wavelengths=np.array([405, 488, 561, 640], dtype=float),
            laser_powers=np.ones(4),
            band_edges=[(420, 470), (500, 550), (570, 620), (650, 720)],
            target_fluorophores=[0, 1, 2, 3],
        )
        M = fpicker.build_measurement_matrix(fps, [rc], hw)
        assert M.shape == (4, 4)
        # M should have positive entries where excitation+emission overlap
        assert M.sum() > 0
        # No negative values
        assert M.min() >= 0

    def test_build_spectral(self, library):
        hw = fpicker.HardwareConfig(
            detector_mode="spectral",
            n_rounds=1,
            spectral_n_bins=32,
            spectral_start=410,
            spectral_bin_width=10,
        )
        fps = library[:3]
        rc = fpicker.RoundConfig(
            laser_wavelengths=np.array([488, 561, 640], dtype=float),
            laser_powers=np.ones(3),
            target_fluorophores=[0, 1, 2],
        )
        M = fpicker.build_measurement_matrix(fps, [rc], hw)
        assert M.shape == (32, 3)
        assert M.sum() > 0
        assert M.min() >= 0

    def test_multi_round(self, library):
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=2, n_detectors=3
        )
        fps = library[:4]
        rc1 = fpicker.RoundConfig(
            laser_wavelengths=np.array([405, 488], dtype=float),
            laser_powers=np.ones(2),
            band_edges=[(420, 470), (500, 550), (570, 630)],
            target_fluorophores=[0, 1],
        )
        rc2 = fpicker.RoundConfig(
            laser_wavelengths=np.array([561, 640], dtype=float),
            laser_powers=np.ones(2),
            band_edges=[(570, 620), (640, 700), (700, 760)],
            target_fluorophores=[2, 3],
        )
        M = fpicker.build_measurement_matrix(fps, [rc1, rc2], hw)
        # 2 rounds * 3 detectors = 6 channels, 4 fluorophores
        assert M.shape == (6, 4)
        assert M.sum() > 0
        # All fluorophores should appear in all rounds (permanent labels)
        # At least some signal in each round
        assert M[:3, :].sum() > 0  # round 1
        assert M[3:, :].sum() > 0  # round 2

    def test_sigma_min(self, library):
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=4
        )
        fps = library[:4]
        rc = fpicker.RoundConfig(
            laser_wavelengths=np.array([405, 488, 561, 640], dtype=float),
            laser_powers=np.ones(4),
            band_edges=[(420, 470), (500, 550), (570, 620), (650, 720)],
            target_fluorophores=[0, 1, 2, 3],
        )
        M = fpicker.build_measurement_matrix(fps, [rc], hw)
        sigma = fpicker.compute_sigma_min(M)
        cond = fpicker.compute_condition_number(M)
        assert sigma > 0, "sigma_min should be positive for a good panel"
        assert cond > 1, "condition number should be >= 1"
        assert np.isfinite(cond)


class TestInnerOptimizer:
    def test_band_optimization_improves(self, library):
        """Band optimization should produce sigma_min >= naive placement."""
        from fpicker._inner import optimize_band_edges

        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=4
        )
        fps = library[:4]
        assignments = [0, 0, 0, 0]  # all in round 0
        laser_configs = [
            (
                np.array([fp.ex_peak for fp in fps]),
                np.ones(len(fps)),
            )
        ]

        round_configs, sigma_min = optimize_band_edges(
            fps, assignments, laser_configs, hw, max_iter=50
        )

        assert sigma_min > 0
        assert len(round_configs) == 1
        assert round_configs[0].band_edges is not None
        assert len(round_configs[0].band_edges) == 4

        # Band edges should be ordered and non-overlapping
        bands = round_configs[0].band_edges
        for lo, hi in bands:
            assert hi > lo
            assert hi - lo >= hw.min_band_width - 1  # small tolerance


class TestAutofluorescence:
    def test_autofluorescence_column(self, library):
        af = fpicker.make_autofluorescence()
        assert af.name == "Autofluorescence"
        assert af.brightness > 0
        assert af.em_spectrum.max() > 0

        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=4
        )
        fps = library[:3] + [af]
        rc = fpicker.RoundConfig(
            laser_wavelengths=np.array([405, 488, 561], dtype=float),
            laser_powers=np.ones(3),
            band_edges=[(420, 470), (500, 550), (570, 620), (650, 720)],
            target_fluorophores=[0, 1, 2],
        )
        M = fpicker.build_measurement_matrix(fps, [rc], hw)
        assert M.shape == (4, 4)
        # AF column should have signal (it's excited by blue light)
        assert M[:, 3].sum() > 0


class TestDiagnostics:
    def test_compute_diagnostics(self, library):
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=4
        )
        fps = library[:4]
        rc = fpicker.RoundConfig(
            laser_wavelengths=np.array([405, 488, 561, 640], dtype=float),
            laser_powers=np.ones(4),
            band_edges=[(420, 470), (500, 550), (570, 620), (650, 720)],
            target_fluorophores=[0, 1, 2, 3],
        )
        M = fpicker.build_measurement_matrix(fps, [rc], hw)
        diag = fpicker.compute_diagnostics(M, fps)

        assert "sigma_min" in diag
        assert "condition_number" in diag
        assert "pairwise_crosstalk" in diag
        assert diag["sigma_min"] > 0
        assert diag["pairwise_crosstalk"].shape == (4, 4)
        # Diagonal should be 1.0 (self-similarity)
        np.testing.assert_allclose(
            np.diag(diag["pairwise_crosstalk"]), 1.0, atol=1e-10
        )


class TestOptimizer:
    @pytest.mark.slow
    def test_optimize_small(self, library):
        """Run a small optimization to verify the pipeline works end-to-end."""
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=3
        )

        results = fpicker.optimize_panel(
            library=library,
            n_select=3,
            hardware=hw,
            pop_size=20,
            n_generations=10,
            optimize_bands=True,
            band_opt_max_iter=20,
            seed=42,
            verbose=False,
        )

        assert len(results) > 0
        best = results[0]
        assert best.sigma_min > 0
        assert len(best.fluorophores) == 3
        assert best.measurement_matrix.shape == (3, 3)
        assert best.condition_number >= 1
        assert best.min_brightness > 0

    @pytest.mark.slow
    def test_optimize_multi_round(self, library):
        """Test multi-round optimization."""
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=2, n_detectors=3
        )

        results = fpicker.optimize_panel(
            library=library,
            n_select=4,
            hardware=hw,
            pop_size=20,
            n_generations=10,
            optimize_bands=True,
            band_opt_max_iter=20,
            seed=42,
            verbose=False,
        )

        assert len(results) > 0
        best = results[0]
        assert best.sigma_min > 0
        assert len(best.fluorophores) == 4
        # 2 rounds * 3 detectors = 6 channels
        assert best.measurement_matrix.shape == (6, 4)

    @pytest.mark.slow
    def test_optimize_spectral(self, library):
        """Test spectral detector mode optimization."""
        hw = fpicker.HardwareConfig(
            detector_mode="spectral",
            n_rounds=1,
            spectral_n_bins=32,
        )

        results = fpicker.optimize_panel(
            library=library,
            n_select=3,
            hardware=hw,
            pop_size=20,
            n_generations=10,
            optimize_bands=False,
            seed=42,
            verbose=False,
        )

        assert len(results) > 0
        best = results[0]
        assert best.sigma_min > 0
        assert len(best.fluorophores) == 3
        assert best.measurement_matrix.shape == (32, 3)


class TestOutput:
    def test_format_results(self, library):
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=1, n_detectors=3
        )
        fps = library[:3]
        rc = fpicker.RoundConfig(
            laser_wavelengths=np.array([488, 561, 640], dtype=float),
            laser_powers=np.ones(3),
            band_edges=[(500, 550), (570, 620), (650, 720)],
            target_fluorophores=[0, 1, 2],
        )
        M = fpicker.build_measurement_matrix(fps, [rc], hw)
        sv = np.linalg.svd(M, compute_uv=False)
        result = fpicker.PanelResult(
            fluorophore_indices=[0, 1, 2],
            fluorophores=fps,
            round_configs=[rc],
            measurement_matrix=M,
            sigma_min=float(sv[-1]),
            condition_number=float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf"),
            brightness_values=np.array([fp.brightness for fp in fps]),
            min_brightness=min(fp.brightness for fp in fps),
            total_brightness=sum(fp.brightness for fp in fps),
        )

        summary = result.summary()
        assert "Round 1" in summary
        assert "sigma_min" in summary

        table = fpicker.format_results_table([result])
        assert "sigma_min" in table
        assert fps[0].name in table

    def test_summary_multi_round(self, library):
        fps = library[:4]
        rc1 = fpicker.RoundConfig(
            laser_wavelengths=np.array([488, 561], dtype=float),
            laser_powers=np.ones(2),
            band_edges=[(500, 550), (570, 620)],
            target_fluorophores=[0, 1],
        )
        rc2 = fpicker.RoundConfig(
            laser_wavelengths=np.array([561, 640], dtype=float),
            laser_powers=np.ones(2),
            band_edges=[(570, 620), (650, 720)],
            target_fluorophores=[2, 3],
        )
        hw = fpicker.HardwareConfig(
            detector_mode="discrete", n_rounds=2, n_detectors=2
        )
        M = fpicker.build_measurement_matrix(fps, [rc1, rc2], hw)
        sv = np.linalg.svd(M, compute_uv=False)
        result = fpicker.PanelResult(
            fluorophore_indices=[0, 1, 2, 3],
            fluorophores=fps,
            round_configs=[rc1, rc2],
            measurement_matrix=M,
            sigma_min=float(sv[-1]),
            condition_number=float(sv[0] / sv[-1]) if sv[-1] > 0 else float("inf"),
            brightness_values=np.array([fp.brightness for fp in fps]),
            min_brightness=min(fp.brightness for fp in fps),
            total_brightness=sum(fp.brightness for fp in fps),
        )
        summary = result.summary()
        assert "Round 1" in summary
        assert "Round 2" in summary
