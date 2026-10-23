#!/usr/bin/env python3
#*******************************************************************************
#  block_average.py
#
#  Average a Cassandra ``.prp`` file after discarding the part of the run that
#  has not settled down, and estimate the uncertainty with block averages.
#
#  Why this exists
#  ---------------
#  An NPT simulation writes an instantaneous density, energy, and pressure
#  every ``prop_freq`` sweeps.  The first part of that list can still be
#  relaxing.  Averaging the whole file then mixes the relaxed liquid with the
#  leftover drift.  This tool chooses a starting sweep, averages only what
#  comes after it, and reports an error bar.
#
#  Which series decides the starting sweep?
#  ----------------------------------------
#  By default the choice is made from ``Mass_Density`` (else ``Volume``, else
#  the first property column).  Every other column is then averaged over that
#  same window.  Density is the right guide for an NPT liquid: a drifting
#  density is the sign that the box has not settled.  The NVT restart box
#  itself is set from the mean volume, not from the mean of the density
#  column.  See "Two densities" below.
#
#  Equilibration rule (Chodera, J. Chem. Theory Comput. 2016, 12, 1799)
#  -------------------------------------------------------------------
#  Let the samples be y_0 ... y_{N-1}.  For a candidate discard point t0,
#  keep y_{t0} ... y_{N-1} and estimate how many of those samples are
#  independent.  The autocorrelation C(t) measures whether sample i still
#  looks like sample i+t.  The statistical inefficiency g is
#
#      g = 1 + 2 * sum_t C(t)
#
#  summed only while C(t) stays positive (Geyer's initial positive sequence,
#  Statistical Science 1992, 7, 473).
#  g = 1 means neighboring samples are independent.  g = 10 means you need
#  about 10 stored samples to get one independent sample.
#
#  The number of independent samples left after the discard is
#
#      N_eff = (N - t0) / g
#
#  The tool tries discard points in the first half of the series and keeps
#  the t0 with the largest N_eff.  A flat series keeps t0 = 0, because
#  throwing data away only reduces N_eff.  A drifting prefix is discarded
#  when doing so produces a tighter mean.
#
#  ``--discard-sweeps`` skips this search and uses the sweep you name.
#
#  Uncertainty: two methods, selected with --error-method
#  ----------------------------------------------------
#  autocorr
#      Estimate g from C(t) (Geyer, Statistical Science 1992, 7, 473), then
#      form blocks of length ceil(3 g) stored samples and take the standard
#      error of those block means.  If the series is too short for two such
#      blocks, use SEM = s * sqrt(g / n), which is Geyer's variance of the
#      mean.  The block length ceil(3 g) is a heuristic: it trusts the
#      autocorrelation sum.
#
#  blocking
#      Do not use g.  Repeat the standard-error calculation at block lengths
#      B = 1, 2, 4, 8, ... stored samples (Flyvbjerg and Petersen, J. Chem.
#      Phys. 1989, 91, 461).  For uncorrelated data the SEM is about the same
#      at every B.  For correlated data it starts too small and rises until
#      a block is long enough to be independent, then it plateaus.  The tool
#      takes the start of that plateau: the smallest B whose SEM agrees with
#      the SEM at 2B, within the uncertainty of those two SEM estimates.
#      If the curve is still rising at the largest B that still leaves four
#      blocks, that largest SEM is reported and the output says no plateau
#      was found.
#
#  both
#      Print the two error bars side by side, plus the SEM-versus-B table.
#      The mean is the same either way.  ``both`` is the default.
#
#  The discard point is still chosen from g (the autocorr method).  Blocking
#  changes the error bar, not which samples enter the average.
#
#  Figure (--plot)
#  --------------
#  ``--plot file.png`` draws one column: the instantaneous trace, a line at
#  the average, and a shaded band around that line.  ``--sem-band`` sets the
#  half-width of the band in multiples of the standard error of the mean.
#  The choices are 1 (the default), 2, and 3.  When both error bars were
#  computed, the band uses the larger one.  The printed table is always one
#  standard error; the multiplier is only on the figure.
#
#  The column drawn is the one that chose the discard point, unless
#  ``--plot-property`` names another column.  That option does not move the
#  discard point.  A prefix left out of the average is shaded on the figure.
#
#  Two densities
#  ------------
#  For a fixed number of molecules the instantaneous density is proportional
#  to 1/V.  The mean of the Mass_Density column is therefore proportional to
#  <1/V>.  The density of the same molecules placed in the mean volume <V>
#  is proportional to 1/<V>.  These are not the same number:
#
#      <1/V>  >  1/<V>
#
#  whenever the volume fluctuates.  The difference is small, usually well
#  inside the uncertainty on <rho>, but an NVT run built from <V> will sit
#  on the second number.  The report prints both when Volume and
#  Mass_Density are present.  The total mass used for the second number is
#  the average of (instantaneous density) * (instantaneous volume).
#
#  A separate "flat?" flag is printed for the series that chose t0.  The
#  kept half-means must agree within two combined standard errors, and a
#  straight-line slope must be consistent with zero at the same level.
#  This flag does not change the average.  It tells you whether to trust it.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/16/26
#
#*******************************************************************************
"""Block-average a Cassandra ``.prp`` file after equilibration detection."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.prp import PrpData, read_prp, resolve_property


@dataclass
class BlockingPoint:
    """Standard error at one block length."""

    block_length: int
    n_blocks: int
    sem: float
    sem_uncertainty: float


@dataclass
class BlockingResult:
    """SEM versus block length, and the plateau value used as an error bar."""

    points: list[BlockingPoint]
    plateau_length: int
    plateau_sem: float
    plateau_n_blocks: int
    plateau_found: bool
    note: str


@dataclass
class PropertyAverage:
    """Mean and uncertainty of one ``.prp`` column on the kept window."""

    name: str
    unit: str
    mean: float
    sem: float
    g: float
    n_used: int
    n_blocks: int
    block_length: int
    sem_method: str
    blocking: BlockingResult | None = None


@dataclass
class AverageReport:
    """Everything the restart tool and the printed summary need."""

    source: Path
    step_label: str
    equilibration_property: str
    discard_index: int
    discard_sweep: float
    n_total: int
    n_used: int
    flat: bool
    flat_notes: list[str]
    properties: list[PropertyAverage]
    volume_average: float | None
    cubic_length: float | None
    density_mean: float | None = None
    density_at_mean_volume: float | None = None
    error_method: str = "both"

    def as_dict(self) -> dict:
        """JSON-ready copy (paths and NumPy values become plain Python)."""
        return {
            "source": str(self.source),
            "step_label": self.step_label,
            "equilibration_property": self.equilibration_property,
            "discard_index": self.discard_index,
            "discard_sweep": self.discard_sweep,
            "n_total": self.n_total,
            "n_used": self.n_used,
            "flat": self.flat,
            "flat_notes": list(self.flat_notes),
            "volume_average": self.volume_average,
            "cubic_length_A": self.cubic_length,
            "density_mean": self.density_mean,
            "density_at_mean_volume": self.density_at_mean_volume,
            "error_method": self.error_method,
            "properties": [
                {
                    "name": p.name,
                    "unit": p.unit,
                    "mean": p.mean,
                    "sem": p.sem,
                    "g": p.g,
                    "n_used": p.n_used,
                    "n_blocks": p.n_blocks,
                    "block_length": p.block_length,
                    "sem_method": p.sem_method,
                    "blocking": None
                    if p.blocking is None
                    else {
                        "plateau_length": p.blocking.plateau_length,
                        "plateau_sem": p.blocking.plateau_sem,
                        "plateau_n_blocks": p.blocking.plateau_n_blocks,
                        "plateau_found": p.blocking.plateau_found,
                        "note": p.blocking.note,
                        "curve": [
                            {
                                "block_length": pt.block_length,
                                "n_blocks": pt.n_blocks,
                                "sem": pt.sem,
                                "sem_uncertainty": pt.sem_uncertainty,
                            }
                            for pt in p.blocking.points
                        ],
                    },
                }
                for p in self.properties
            ],
        }


def statistical_inefficiency(series: np.ndarray) -> float:
    """Return g, the number of stored samples per independent sample.

    g = 1 for an uncorrelated series.  The sum over the autocorrelation
    stops at the first pair of lags whose autocorrelations add to a
    negative number, which is Geyer's initial positive sequence
    (Statistical Science 1992, 7, 473).
    """
    x = np.asarray(series, dtype=float)
    n = int(x.size)
    if n < 4:
        return 1.0
    x = x - float(x.mean())
    c0 = float(np.dot(x, x) / n)
    if c0 <= 0.0:
        return 1.0

    nfft = 1 << int(np.ceil(np.log2(2 * n)))
    fft = np.fft.fft(x, n=nfft)
    acf = np.fft.ifft(fft * np.conjugate(fft))[:n].real / (n * c0)
    acf = acf / acf[0]

    g = 1.0
    t = 1
    while t < n - 1:
        pair = float(acf[t])
        if t + 1 < n:
            pair += float(acf[t + 1])
        if pair <= 0.0:
            break
        g += 2.0 * float(acf[t])
        if t + 1 < n:
            g += 2.0 * float(acf[t + 1])
        t += 2
    return max(g, 1.0)


def equilibration_index(series: np.ndarray) -> tuple[int, float]:
    """Index t0 in the first half of ``series`` that maximizes N_eff.

    Returns ``(t0, g)`` for the kept suffix ``series[t0:]``.
    """
    x = np.asarray(series, dtype=float)
    n = int(x.size)
    if n < 8:
        return 0, statistical_inefficiency(x)

    step = max(1, n // 100)
    best_i = 0
    best_neff = -1.0
    best_g = 1.0
    for i in range(0, n // 2 + 1, step):
        y = x[i:]
        if y.size < 8:
            continue
        g = statistical_inefficiency(y)
        neff = float(y.size) / g
        if neff > best_neff:
            best_neff = neff
            best_i = i
            best_g = g
    return best_i, best_g


def _block_sem(series: np.ndarray, g: float) -> tuple[float, int, int, str]:
    """Standard error from block means, or the g-based formula if needed.

    Returns ``(sem, n_blocks, block_length, method)``.
    """
    y = np.asarray(series, dtype=float)
    n = int(y.size)
    if n < 2:
        return float("nan"), 0, n, "too-short"

    block_length = max(1, int(np.ceil(3.0 * g)))
    n_blocks = n // block_length
    if n_blocks < 2:
        sem = float(y.std(ddof=1) * np.sqrt(g / n))
        return sem, 1, n, "g-formula"

    trimmed = y[: n_blocks * block_length].reshape(n_blocks, block_length)
    means = trimmed.mean(axis=1)
    sem = float(means.std(ddof=1) / np.sqrt(n_blocks))
    return sem, n_blocks, block_length, "autocorr-blocks"


def blocking_analysis(series: np.ndarray, min_blocks: int = 4) -> BlockingResult:
    """Standard error at block lengths 1, 2, 4, ... and the plateau value.

    ``sem_uncertainty`` is the estimated noise on the SEM itself,
    SEM / sqrt(2 (N_blocks - 1)).  A plateau is the smallest block length
    whose SEM agrees with the SEM at twice that length, within the sum of
    those two uncertainties.  Agreement means further blocking no longer
    moves the error bar by more than its own noise.
    """
    y = np.asarray(series, dtype=float)
    n = int(y.size)
    points: list[BlockingPoint] = []
    block_length = 1
    while n // block_length >= min_blocks:
        n_blocks = n // block_length
        means = y[: n_blocks * block_length].reshape(n_blocks, block_length).mean(axis=1)
        sem = float(means.std(ddof=1) / np.sqrt(n_blocks))
        sem_uncertainty = float(sem / np.sqrt(2.0 * (n_blocks - 1)))
        points.append(
            BlockingPoint(
                block_length=block_length,
                n_blocks=n_blocks,
                sem=sem,
                sem_uncertainty=sem_uncertainty,
            )
        )
        block_length *= 2

    if not points:
        naive = float(y.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
        return BlockingResult(
            points=[],
            plateau_length=n,
            plateau_sem=naive,
            plateau_n_blocks=1,
            plateau_found=False,
            note="too few samples for a blocking curve; SEM is the ordinary standard error",
        )

    plateau_at: int | None = None
    for i in range(len(points) - 1):
        here, nxt = points[i], points[i + 1]
        tolerance = here.sem_uncertainty + nxt.sem_uncertainty
        if nxt.sem <= here.sem + tolerance:
            plateau_at = i
            break

    if plateau_at is None:
        chosen = max(points, key=lambda pt: pt.sem)
        return BlockingResult(
            points=points,
            plateau_length=chosen.block_length,
            plateau_sem=chosen.sem,
            plateau_n_blocks=chosen.n_blocks,
            plateau_found=False,
            note=(
                "SEM was still changing at the largest block length that leaves "
                f"at least {min_blocks} blocks; reported that largest SEM"
            ),
        )

    chosen = points[plateau_at]
    return BlockingResult(
        points=points,
        plateau_length=chosen.block_length,
        plateau_sem=chosen.sem,
        plateau_n_blocks=chosen.n_blocks,
        plateau_found=True,
        note=(
            f"plateau begins at block length {chosen.block_length} "
            f"({chosen.n_blocks} blocks); the SEM there agrees with the next doubling"
        ),
    )


def _flatness(steps: np.ndarray, series: np.ndarray, g: float) -> tuple[bool, list[str]]:
    """True when the kept series has no obvious leftover drift.

    Both tests use a two-standard-error band so ordinary noise in a short
    liquid run is not called drift.
    """
    notes: list[str] = []
    y = np.asarray(series, dtype=float)
    t = np.asarray(steps, dtype=float)
    n = int(y.size)
    if n < 8:
        return False, ["fewer than 8 points remain; flatness was not tested"]

    mid = n // 2
    left, right = y[:mid], y[mid:]
    g_left = statistical_inefficiency(left)
    g_right = statistical_inefficiency(right)
    sem_left = float(left.std(ddof=1) * np.sqrt(g_left / left.size))
    sem_right = float(right.std(ddof=1) * np.sqrt(g_right / right.size))
    combined = float(np.sqrt(sem_left**2 + sem_right**2))
    half_gap = abs(float(left.mean()) - float(right.mean()))
    halves_ok = half_gap <= 2.0 * combined
    if not halves_ok:
        notes.append(
            f"first and second half differ by {half_gap:.6g}, "
            f"which is outside 2 combined SEM ({2.0 * combined:.6g})"
        )

    # Slope uncertainty for correlated samples: inflate the residual
    # variance by g, then use the usual straight-line formula.
    slope = float(np.polyfit(t, y, 1)[0])
    centered = t - t.mean()
    ss = float(np.dot(centered, centered))
    sigma2 = float(y.var(ddof=1) * g)
    sem_slope = float(np.sqrt(sigma2 / ss)) if ss > 0.0 else float("inf")
    slope_ok = abs(slope) <= 2.0 * sem_slope
    if not slope_ok:
        notes.append(
            f"slope {slope:.6g} per sweep is outside 2 SEM ({2.0 * sem_slope:.6g})"
        )
    if halves_ok and slope_ok:
        notes.append("kept window looks flat (halves and slope within 2 SEM)")
    return halves_ok and slope_ok, notes


def _driver_name(data: PrpData, requested: str | None) -> str:
    """Property whose equilibration window is applied to every column."""
    if requested:
        return resolve_property(data.property_names, requested)
    for candidate in ("Mass_Density", "Volume"):
        if candidate in data.property_names:
            return candidate
    return data.property_names[0]


def analyze_prp(
    data: PrpData,
    *,
    driver: str | None = None,
    discard_sweeps: float | None = None,
    error_method: str = "both",
) -> AverageReport:
    """Average every column of ``data`` on one shared production window.

    ``error_method`` is ``autocorr``, ``blocking``, or ``both``.  It selects
    the error bar.  The samples that enter the mean are chosen first, from
    the autocorrelation rule or from ``discard_sweeps``.
    """
    if error_method not in {"autocorr", "blocking", "both"}:
        raise ValueError("error_method must be autocorr, blocking, or both")

    driver_name = _driver_name(data, driver)
    driver_series = data.column(driver_name)

    if discard_sweeps is None:
        discard_index, _ = equilibration_index(driver_series)
    else:
        hits = np.where(data.steps + 1e-9 >= float(discard_sweeps))[0]
        if hits.size == 0:
            raise ValueError(
                f"--discard-sweeps {discard_sweeps} is past the last "
                f"{data.step_label} ({data.steps[-1]:.0f}) in {data.path}"
            )
        discard_index = int(hits[0])

    kept_steps = data.steps[discard_index:]
    want_blocking = error_method in {"blocking", "both"}
    averages: list[PropertyAverage] = []
    for name, unit in zip(data.property_names, data.property_units):
        kept = data.column(name)[discard_index:]
        g = statistical_inefficiency(kept)
        sem_g, n_blocks, block_length, method = _block_sem(kept, g)
        blocking = blocking_analysis(kept) if want_blocking else None
        if error_method == "blocking" and blocking is not None:
            reported_sem = blocking.plateau_sem
            reported_blocks = blocking.plateau_n_blocks
            reported_length = blocking.plateau_length
            reported_method = "blocking-plateau" if blocking.plateau_found else "blocking-no-plateau"
        else:
            reported_sem = sem_g
            reported_blocks = n_blocks
            reported_length = block_length
            reported_method = method
        averages.append(
            PropertyAverage(
                name=name,
                unit=unit,
                mean=float(kept.mean()),
                sem=reported_sem,
                g=g,
                n_used=int(kept.size),
                n_blocks=reported_blocks,
                block_length=reported_length,
                sem_method=reported_method,
                blocking=blocking,
            )
        )

    flat, notes = _flatness(kept_steps, driver_series[discard_index:], statistical_inefficiency(driver_series[discard_index:]))
    by_name = {item.name: item for item in averages}
    volume = by_name["Volume"].mean if "Volume" in by_name else None
    length = float(volume ** (1.0 / 3.0)) if volume is not None and volume > 0.0 else None
    density_mean = by_name["Mass_Density"].mean if "Mass_Density" in by_name else None
    density_from_volume = None
    if density_mean is not None and volume is not None and volume > 0.0:
        kept_density = data.column("Mass_Density")[discard_index:]
        kept_volume = data.column("Volume")[discard_index:]
        if np.all(kept_volume > 0.0):
            # Instantaneous density times volume is the fixed total mass,
            # in the units Cassandra used for those two columns.
            total_mass = float(np.mean(kept_density * kept_volume))
            density_from_volume = total_mass / volume

    return AverageReport(
        source=data.path,
        step_label=data.step_label,
        equilibration_property=driver_name,
        discard_index=discard_index,
        discard_sweep=float(data.steps[discard_index]),
        n_total=data.n_rows,
        n_used=int(kept_steps.size),
        flat=flat,
        flat_notes=notes,
        properties=averages,
        volume_average=volume,
        cubic_length=length,
        density_mean=density_mean,
        density_at_mean_volume=density_from_volume,
        error_method=error_method,
    )


def _format_blocking_curve(prop: PropertyAverage) -> list[str]:
    """SEM versus block length for one property."""
    blocking = prop.blocking
    if blocking is None:
        return []
    lines = [
        f"Blocking curve for {prop.name}" + (f" [{prop.unit}]" if prop.unit else ""),
        f"  {blocking.note}",
        f"  {'B':>8} {'blocks':>8} {'SEM':>14} {'SEM noise':>14}",
    ]
    for pt in blocking.points:
        marker = "  <-- plateau" if pt.block_length == blocking.plateau_length else ""
        lines.append(
            f"  {pt.block_length:8d} {pt.n_blocks:8d} {pt.sem:14.6g} "
            f"{pt.sem_uncertainty:14.6g}{marker}"
        )
    return lines


def format_report(report: AverageReport) -> str:
    """Plain-text summary a student can read in the terminal."""
    lines = [
        f"File: {report.source}",
        f"Equilibration property: {report.equilibration_property}",
        f"Error method: {report.error_method}",
        (
            f"Discard before {report.step_label} {report.discard_sweep:.0f} "
            f"(index {report.discard_index}); use {report.n_used} of {report.n_total} samples"
        ),
        f"Flat?: {'yes' if report.flat else 'no'}",
    ]
    lines.extend(f"  {note}" for note in report.flat_notes)
    if report.volume_average is not None and report.cubic_length is not None:
        lines.append(
            f"Average volume <V>: {report.volume_average:.8g} A^3"
            f"    cubic length L = <V>^(1/3) = {report.cubic_length:.8g} A"
        )
    if report.density_mean is not None and report.density_at_mean_volume is not None:
        unit = next(item.unit for item in report.properties if item.name == "Mass_Density")
        lines.append(
            f"Mean instantaneous density <rho>: {report.density_mean:.8g} {unit}"
        )
        lines.append(
            "  <rho> is the average of the Mass_Density column. For fixed N it is proportional to <1/V>."
        )
        lines.append(
            f"Density at the mean volume: {report.density_at_mean_volume:.8g} {unit}"
        )
        lines.append(
            "  This is (total mass)/<V>, proportional to 1/<V>. "
            "An NVT box built from <V> has this density."
        )
    lines.append("")
    if report.error_method == "both":
        lines.append(
            f"{'Property':<16} {'mean':>14} {'SEM autocorr':>14} "
            f"{'SEM blocking':>14}  unit"
        )
        for prop in report.properties:
            block_sem = float("nan") if prop.blocking is None else prop.blocking.plateau_sem
            lines.append(
                f"{prop.name:<16} {prop.mean:14.6g} {prop.sem:14.6g} "
                f"{block_sem:14.6g}  {prop.unit}  "
                f"[g={prop.g:.2f}, B={prop.blocking.plateau_length if prop.blocking else 0}]"
            )
    else:
        lines.append(f"{'Property':<16} {'mean':>14} {'SEM':>14} {'blocks':>8}  unit")
        for prop in report.properties:
            lines.append(
                f"{prop.name:<16} {prop.mean:14.6g} {prop.sem:14.6g} "
                f"{prop.n_blocks:8d}  {prop.unit}  "
                f"[{prop.sem_method}, g={prop.g:.2f}, B={prop.block_length}]"
            )
    if report.error_method in {"blocking", "both"}:
        for prop in report.properties:
            lines.append("")
            lines.extend(_format_blocking_curve(prop))
    return "\n".join(lines)


def plot_band_sem(prop: PropertyAverage) -> tuple[float, str]:
    """One-standard-error half-width for the figure, and which estimate it is.

    When both the autocorrelation and blocking estimates exist, the larger
    one is used.  A blocking-only report has no separate autocorrelation
    value stored, so its reported SEM is the band.
    """
    blocking = prop.blocking
    if blocking is not None and not prop.sem_method.startswith("blocking"):
        if blocking.plateau_sem > prop.sem:
            return blocking.plateau_sem, "blocking"
        return prop.sem, "autocorrelation"
    if prop.sem_method.startswith("blocking"):
        return prop.sem, "blocking"
    return prop.sem, "autocorrelation"


def write_plot(
    data: PrpData,
    report: AverageReport,
    path: Path,
    *,
    property_name: str | None = None,
    sem_band: int = 1,
) -> None:
    """Save one series with its average and a shaded standard-error band.

    ``property_name`` selects the curve.  The discard window still comes
    from ``report``.  ``sem_band`` is 1, 2, or 3, and multiplies the
    standard error of the mean.
    """
    if sem_band not in (1, 2, 3):
        raise ValueError("sem_band must be 1, 2, or 3")

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    name = report.equilibration_property if property_name is None else resolve_property(
        data.property_names, property_name
    )
    prop = next(item for item in report.properties if item.name == name)
    sem, which = plot_band_sem(prop)
    half = sem_band * sem
    y = data.column(name)
    x = data.steps
    x0 = float(x[report.discard_index])
    x1 = float(x[-1])

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    ax.plot(x, y, color="0.55", lw=0.9, label="instantaneous", zorder=2)
    if report.discard_index > 0:
        ax.axvspan(
            float(x[0]),
            x0,
            color="#e69f00",
            alpha=0.18,
            label="not in the average",
            zorder=1,
        )
    if np.isfinite(half):
        count = "1 SEM" if sem_band == 1 else f"{sem_band} SEM"
        ax.fill_between(
            [x0, x1],
            [prop.mean - half, prop.mean - half],
            [prop.mean + half, prop.mean + half],
            color="#0072b2",
            alpha=0.35,
            label=f"± {_format_band(half)} ({count}, {which})",
            zorder=3,
        )
    ax.plot(
        [x0, x1],
        [prop.mean, prop.mean],
        color="#0072b2",
        lw=1.8,
        label=f"average {_format_band(prop.mean)}",
        zorder=4,
    )
    ax.set_xlabel(data.step_label.replace("_", " "))
    unit = data.unit_for(name)
    ax.set_ylabel(name + (f" [{unit}]" if unit else ""))
    ax.set_title(data.path.name)
    ax.legend(frameon=True, framealpha=0.92, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _format_band(value: float) -> str:
    """Compact legend number: one decimal for large values, more for small ones."""
    magnitude = abs(value)
    if magnitude >= 10.0:
        return f"{value:.1f}"
    if magnitude >= 1.0:
        return f"{value:.2f}"
    return f"{value:.3g}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Block-average a Cassandra .prp file after dropping the unequilibrated prefix."
    )
    parser.add_argument("prp", type=Path, help="Cassandra .prp property file")
    parser.add_argument(
        "--property",
        help="Series used to choose the discard point (default: Mass_Density, else Volume)",
    )
    parser.add_argument(
        "--discard-sweeps",
        type=float,
        help="Use this sweep as the start of the average instead of the automatic choice",
    )
    parser.add_argument(
        "--error-method",
        choices=("autocorr", "blocking", "both"),
        default="both",
        help="Uncertainty estimate: autocorrelation blocks, the SEM-vs-block-size plateau, or both (default)",
    )
    parser.add_argument("--json", type=Path, help="Write a JSON sidecar for npt_restart.py")
    parser.add_argument(
        "--plot",
        type=Path,
        help="Write a PNG of one series with its average and a shaded standard-error band",
    )
    parser.add_argument(
        "--plot-property",
        help="Column drawn by --plot (default: the series that chose the discard point). Does not move that cut.",
    )
    parser.add_argument(
        "--sem-band",
        type=int,
        choices=(1, 2, 3),
        default=1,
        help="Half-width of the --plot band, in multiples of the standard error of the mean (default: 1)",
    )
    args = parser.parse_args(argv)
    if args.plot is None and (args.plot_property or args.sem_band != 1):
        parser.error("--plot-property and --sem-band apply to the figure; pass --plot as well")

    data = read_prp(args.prp)
    report = analyze_prp(
        data,
        driver=args.property,
        discard_sweeps=args.discard_sweeps,
        error_method=args.error_method,
    )
    print(format_report(report))
    if args.json is not None:
        args.json.write_text(json.dumps(report.as_dict(), indent=2) + "\n", encoding="utf-8")
        print(f"\nWrote {args.json}")
    if args.plot is not None:
        write_plot(
            data,
            report,
            args.plot,
            property_name=args.plot_property,
            sem_band=args.sem_band,
        )
        print(f"Wrote {args.plot}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
