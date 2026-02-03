"""Output formatting and visualization."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

    from ._models import FluorophoreRecord, PanelResult, RoundConfig


def plot_panel(
    result: PanelResult,
    figsize: tuple[float, float] = (14, 5),
) -> Figure:
    """Plot spectra with laser lines and detection bands overlaid.

    Parameters
    ----------
    result : PanelResult
        A panel optimization result.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    from ._models import WAVELENGTH_GRID

    n_rounds = len(result.round_configs)
    fig, axes = plt.subplots(1, n_rounds, figsize=figsize, squeeze=False)

    # Color cycle
    colors = plt.cm.tab10(np.linspace(0, 1, len(result.fluorophores)))

    for r, (rc, ax) in enumerate(zip(result.round_configs, axes[0])):
        ax.set_title(f"Round {r + 1}")
        ax.set_xlabel("Wavelength (nm)")
        ax.set_ylabel("Normalized intensity")
        ax.set_xlim(380, 780)
        ax.set_ylim(0, 1.15)

        # Plot detection bands as shaded regions
        if rc.band_edges is not None:
            for lo, hi in rc.band_edges:
                rect = Rectangle(
                    (lo, 0), hi - lo, 1.1,
                    alpha=0.08, color="gray", zorder=0
                )
                ax.add_patch(rect)

        # Plot laser lines
        for lam in rc.laser_wavelengths:
            ax.axvline(lam, color="black", linestyle="--", alpha=0.4, linewidth=0.8)

        # Plot fluorophore spectra
        for i, (fp, color) in enumerate(zip(result.fluorophores, colors)):
            is_target = i in rc.target_fluorophores
            lw = 2.0 if is_target else 0.8
            alpha = 1.0 if is_target else 0.4

            # Excitation (dashed)
            ax.plot(
                WAVELENGTH_GRID,
                fp.ex_spectrum,
                color=color,
                linestyle="--",
                linewidth=lw * 0.7,
                alpha=alpha * 0.7,
            )
            # Emission (solid)
            ax.plot(
                WAVELENGTH_GRID,
                fp.em_spectrum,
                color=color,
                linestyle="-",
                linewidth=lw,
                alpha=alpha,
                label=fp.name if r == 0 else None,
            )

        if r == 0:
            ax.legend(fontsize=7, loc="upper right")

    fig.tight_layout()
    return fig


def plot_measurement_matrix(
    result: PanelResult,
    figsize: tuple[float, float] = (8, 6),
) -> Figure:
    """Plot the measurement matrix as a heatmap.

    Parameters
    ----------
    result : PanelResult
        A panel optimization result.
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    M = result.measurement_matrix
    fig, ax = plt.subplots(figsize=figsize)

    im = ax.imshow(M, aspect="auto", cmap="viridis")
    ax.set_xlabel("Fluorophore")
    ax.set_ylabel("Channel (round x detector)")

    fp_names = [fp.name for fp in result.fluorophores]
    ax.set_xticks(range(len(fp_names)))
    ax.set_xticklabels(fp_names, rotation=45, ha="right", fontsize=8)

    n_rounds = len(result.round_configs)
    D = M.shape[0] // n_rounds
    ytick_labels = [f"R{r + 1}D{d + 1}" for r in range(n_rounds) for d in range(D)]
    ax.set_yticks(range(len(ytick_labels)))
    ax.set_yticklabels(ytick_labels, fontsize=8)

    plt.colorbar(im, ax=ax, label="M[(r,d), f]")
    ax.set_title(
        f"Measurement Matrix  |  "
        f"sigma_min={result.sigma_min:.4f}  cond={result.condition_number:.2f}"
    )
    fig.tight_layout()
    return fig


def format_results_table(results: list[PanelResult]) -> str:
    """Format a list of Pareto-optimal results as a text table.

    Parameters
    ----------
    results : list of PanelResult

    Returns
    -------
    str
        Formatted table.
    """
    if not results:
        return "No results found."

    lines = []
    header = (
        f"{'#':>3}  {'sigma_min':>10}  {'cond':>8}  "
        f"{'min_bright':>11}  {'total_bright':>13}  {'Panel'}"
    )
    lines.append(header)
    lines.append("-" * len(header))

    for i, r in enumerate(results):
        names = ", ".join(fp.name for fp in r.fluorophores)
        lines.append(
            f"{i + 1:>3}  {r.sigma_min:>10.4f}  {r.condition_number:>8.2f}  "
            f"{r.min_brightness:>11.0f}  {r.total_brightness:>13.0f}  {names}"
        )

    return "\n".join(lines)
