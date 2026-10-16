#!/usr/bin/env python3
#*******************************************************************************
#  test_block_average.py
#
#  Synthetic .prp files for the equilibration detector.
#  A drifting prefix must be trimmed.  A stationary series must not be.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/16/26
#*******************************************************************************
"""Check block_average equilibration detection on made-up .prp files."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.prp import read_prp
from utilities.block_average import (
    BlockingResult,
    PropertyAverage,
    analyze_prp,
    blocking_analysis,
    plot_band_sem,
    statistical_inefficiency,
    write_plot,
)


def _write_prp(path: Path, density: np.ndarray) -> None:
    """Write a two-column .prp (density and a matching cubic volume)."""
    # 800 molecules of a 120 g/mol liquid are not needed here.  Volume is
    # just an inverse stand-in so the Volume column exists for the restart
    # length formula.  N_fake * mass is absorbed into a constant of 1e8.
    volume = 1.0e8 / density
    lines = [
        "# Instantaneous properties",
        "         # MC_SWEEP      Mass_Density            Volume",
        "                  #          (kg/m^3)             (A^3)",
    ]
    for i, (rho, vol) in enumerate(zip(density, volume), start=1):
        sweep = i * 1000
        lines.append(f"{sweep:18d} {rho:16.6E} {vol:16.6E}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    rng = np.random.default_rng(16)
    # Uncorrelated constant series: g must stay near 1.
    g_flat_noise = statistical_inefficiency(rng.normal(size=400))
    if not 0.5 <= g_flat_noise <= 3.0:
        raise SystemExit(f"expected g near 1 for white noise, got {g_flat_noise}")

    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)

        drift = np.empty(200)
        drift[:50] = np.linspace(1000.0, 1400.0, 50)
        drift[50:] = 1400.0 + rng.normal(0.0, 1.5, 150)
        drift_path = root / "drift.prp"
        _write_prp(drift_path, drift)
        drift_report = analyze_prp(read_prp(drift_path))
        if drift_report.discard_index < 40:
            raise SystemExit(
                f"drifting prefix was not trimmed (discard index {drift_report.discard_index})"
            )
        if abs(drift_report.properties[0].mean - 1400.0) > 5.0:
            raise SystemExit(
                f"kept density mean {drift_report.properties[0].mean} is far from 1400"
            )

        flat = 1400.0 + rng.normal(0.0, 1.5, 200)
        flat_path = root / "flat.prp"
        _write_prp(flat_path, flat)
        flat_report = analyze_prp(read_prp(flat_path))
        if flat_report.discard_index > 20:
            raise SystemExit(
                f"stationary series discarded too much (index {flat_report.discard_index})"
            )
        if not flat_report.flat:
            raise SystemExit(f"stationary series was not flagged flat: {flat_report.flat_notes}")
        if flat_report.cubic_length is None or flat_report.volume_average is None:
            raise SystemExit("volume average or cubic length missing")
        expected_l = flat_report.volume_average ** (1.0 / 3.0)
        if abs(flat_report.cubic_length - expected_l) > 1e-8:
            raise SystemExit("cubic length is not V^(1/3)")

        plot_path = root / "flat.png"
        flat_data = read_prp(flat_path)
        write_plot(flat_data, flat_report, plot_path, property_name="Volume", sem_band=2)
        if not plot_path.is_file() or plot_path.stat().st_size < 1000:
            raise SystemExit("standard-error band figure was not written")
        try:
            write_plot(flat_data, flat_report, plot_path, sem_band=4)
        except ValueError:
            pass
        else:
            raise SystemExit("sem_band outside 1, 2, 3 was accepted")

    wider = plot_band_sem(
        PropertyAverage(
            name="Energy_Total",
            unit="(kJ/mol)-Ext",
            mean=10.0,
            sem=1.5,
            g=1.0,
            n_used=10,
            n_blocks=5,
            block_length=2,
            sem_method="autocorr-blocks",
            blocking=BlockingResult(
                points=[],
                plateau_length=1,
                plateau_sem=2.0,
                plateau_n_blocks=10,
                plateau_found=True,
                note="",
            ),
        )
    )
    if wider != (2.0, "blocking"):
        raise SystemExit(f"band should use the larger SEM, got {wider}")

    white = rng.normal(size=1024)
    white_blocks = blocking_analysis(white)
    naive = float(white.std(ddof=1) / np.sqrt(white.size))
    if not white_blocks.plateau_found or white_blocks.plateau_length > 2:
        raise SystemExit(
            f"white noise should plateau at block length 1 or 2, got {white_blocks.plateau_length}"
        )
    if abs(white_blocks.plateau_sem - naive) / naive > 0.35:
        raise SystemExit(
            f"white-noise blocking SEM {white_blocks.plateau_sem} is far from {naive}"
        )

    # An AR(1) process with phi = 0.95 stays correlated for many samples.
    # The SEM at B = 1 treats those samples as independent and is too small.
    phi = 0.95
    correlated = np.empty(4096)
    correlated[0] = rng.normal()
    for i in range(1, correlated.size):
        correlated[i] = phi * correlated[i - 1] + np.sqrt(1.0 - phi**2) * rng.normal()
    correlated_blocks = blocking_analysis(correlated)
    sem_one = correlated_blocks.points[0].sem
    if correlated_blocks.plateau_sem < 1.5 * sem_one:
        raise SystemExit(
            f"correlated series did not raise the blocking SEM "
            f"({sem_one} at B=1, plateau {correlated_blocks.plateau_sem})"
        )
    if correlated_blocks.plateau_length < 4:
        raise SystemExit(
            f"correlated series plateaued too early (B={correlated_blocks.plateau_length})"
        )

    print("block_average synthetic tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
