#!/usr/bin/env python3
#*******************************************************************************
#  property_plotter.py
#
#  Cassandra V2 diagnostic tool — plot a thermodynamic property from a .prp
#  file versus Monte Carlo sweep (or step).
#
#  ----------------------------
#  After an MC run, Cassandra writes instantaneous properties to a ``.prp``
#  file at intervals set by ``prop_freq`` in the input.  Looking at a single
#  final average is not enough to judge whether the run is behaving well.
#  Plotting the series vs sweep/step shows:
#
#    * drift or equilibration (running average still changing)
#    * large fluctuations (noisy instantaneous curve)
#    * whether the production portion looks stationary
#
#  This script is intentionally simple: read → select a column → plot
#  instantaneous values and a cumulative running average → print the mean.
#
#  Algorithm notes
#  ---------------
#  Running (cumulative) average at point i:
#
#      <y>_i  =  (1 / (i+1)) * sum_{k=0}^{i} y_k
#
#  Implemented with ``numpy.cumsum(y) / (1, 2, ..., N)`` so it is O(N) and
#  does not re-sum the series at every point.  The last running-average
#  value equals the ordinary sample mean of all plotted points.
#
#  Sample standard deviation uses ``ddof=1`` (divide by N−1), the usual
#  unbiased estimator for a sample — not a block-averaged MC error bar.
#
#  Related modules
#  ---------------
#  ``diagnostics/io/prp.py``  — shared .prp parser (read this next)
#  ``docs/output-formats.md`` — .prp / .xyz layout assumptions
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#  Usage examples
#  --------------
#    python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp
#    python diagnostics/property_plotter.py nvt.out.prp -p Energy_Total
#    python diagnostics/property_plotter.py nvt.out.prp -p enthalpy --save H.png
#    python diagnostics/property_plotter.py nvt.out.prp --list
#
#*******************************************************************************
"""Plot a Cassandra ``.prp`` property vs MC sweep/step.

Shows the instantaneous series and a running (cumulative) average, and prints
the final mean to the terminal.  See the file header above for teaching notes.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Import path
#
# Students often run this as a script from the repo root:
#   python diagnostics/property_plotter.py path/to/file.prp
# In that case Python does not automatically treat the repo root as a package
# root.  We insert the repository root onto ``sys.path`` so
# ``import diagnostics.io.prp`` works without installing Cassandra as a pip
# package.
# ---------------------------------------------------------------------------
_DIAG_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _DIAG_ROOT.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.prp import PrpData, read_prp, resolve_property


def _list_properties(data: PrpData) -> None:
    """Print the columns available in the .prp file (names + units)."""
    print(f"File: {data.path}")
    print(f"Title: {data.title}")
    print(f"Step column: {data.step_label}  ({data.n_rows} points)")
    print("Properties:")
    for i, (name, unit) in enumerate(zip(data.property_names, data.property_units), start=1):
        unit_str = f"  [{unit}]" if unit else ""
        print(f"  {i:2d}. {name}{unit_str}")


def _prompt_property(data: PrpData) -> str:
    """Ask the user which property to plot when ``-p`` was not given.

    Accepts either a 1-based index from the printed list or a property name
    (same matching rules as ``resolve_property``: case-insensitive, unique
    substring).  Empty input cancels.
    """
    _list_properties(data)
    print()
    while True:
        raw = input("Property to plot (name or number, empty to cancel): ").strip()
        if not raw:
            raise SystemExit("Cancelled.")
        # Numeric choice → index into property_names
        if raw.isdigit():
            idx = int(raw)
            if 1 <= idx <= data.n_properties:
                return data.property_names[idx - 1]
            print(f"Enter a number between 1 and {data.n_properties}.")
            continue
        # Name / substring choice
        try:
            return resolve_property(data.property_names, raw)
        except ValueError as exc:
            print(exc)


def _running_average(y: np.ndarray) -> np.ndarray:
    """Return the cumulative mean of ``y``.

    For a series y_0, y_1, ..., y_{N-1}::

        avg[i] = mean(y[0], y[1], ..., y[i])

    Implementation:
      * ``np.cumsum(y)`` builds S[i] = y[0] + ... + y[i]
      * divide by the counts 1, 2, ..., N

    This is the simplest equilibration diagnostic: if ``avg`` is still
    trending at the end of the file, the run may not be equilibrated (or
    you may need to discard early points — a future enhancement).
    """
    n = y.size
    if n == 0:
        return y.copy()
    return np.cumsum(y, dtype=float) / np.arange(1, n + 1, dtype=float)


def plot_property(
    data: PrpData,
    property_name: str,
    *,
    save_path: Path | None = None,
    show: bool = True,
) -> float:
    """Plot instantaneous values and running average; print summary stats.

    Parameters
    ----------
    data :
        Parsed ``.prp`` contents from :func:`read_prp`.
    property_name :
        Exact column name as stored in ``data.property_names`` (resolve
        aliases with :func:`resolve_property` before calling).
    save_path :
        Optional path to write the figure (PNG, PDF, …).
    show :
        If True, open Matplotlib's interactive window (``plt.show()``).
        Use False with ``--no-show`` on headless machines.

    Returns
    -------
    float
        Final cumulative mean (same as the ordinary mean of all points).
    """
    # Import here so ``--list`` still works if matplotlib is missing/broken.
    import matplotlib.pyplot as plt

    # y(t): one thermodynamic observable sampled along the MC trajectory
    y = data.column(property_name)
    run_avg = _running_average(y)

    # Summary statistics printed below the plot
    final_mean = float(run_avg[-1])  # == np.mean(y)
    # ddof=1 → sample stdev (N−1); with one point there is no sample stdev
    final_std = float(np.std(y, ddof=1)) if y.size > 1 else 0.0
    unit = data.unit_for(property_name)

    # --- figure ----------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.5))
    # Instantaneous: the raw .prp values at each write interval
    ax.plot(data.steps, y, "o-", markersize=3, linewidth=1.0, label="Instantaneous")
    # Running average: smooths short-term noise; watch for long-term drift
    ax.plot(data.steps, run_avg, "-", linewidth=2.0, label="Running average")
    ax.set_xlabel(data.step_label)
    ylabel = property_name if not unit else f"{property_name} {unit}"
    ax.set_ylabel(ylabel)
    ax.set_title(f"{property_name} vs {data.step_label}\n{data.path.name}")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    # --- text report (always, even if the window is not shown) -----------
    print()
    print(f"Property:     {property_name}" + (f"  {unit}" if unit else ""))
    print(f"Points:       {data.n_rows}")
    print(f"Final mean:   {final_mean:.6e}")
    if y.size > 1:
        print(f"Sample stdev: {final_std:.6e}")
    print(f"Last value:   {float(y[-1]):.6e}")

    if save_path is not None:
        save_path = Path(save_path)
        fig.savefig(save_path, dpi=150)
        print(f"Saved plot:   {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)

    return final_mean


def main(argv: list[str] | None = None) -> int:
    """CLI entry point: parse arguments, load .prp, plot one property.

    Flow
    ----
    1. Parse command-line options (``argparse``).
    2. Read and parse the ``.prp`` file (``read_prp``).
    3. Either list columns and exit, or choose a property
       (``-p``, interactive prompt, or error if non-interactive).
    4. Call :func:`plot_property`.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Plot a Cassandra .prp property vs MC sweep/step "
            "(instantaneous + running average)."
        )
    )
    parser.add_argument(
        "prp_file",
        type=Path,
        help="Path to a Cassandra .prp property file",
    )
    parser.add_argument(
        "-p",
        "--property",
        dest="property_name",
        default=None,
        help="Property name to plot (case-insensitive; unique substring OK)",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List properties in the file and exit",
    )
    parser.add_argument(
        "--save",
        type=Path,
        default=None,
        metavar="PATH",
        help="Save the figure to PATH (PNG/PDF/…)",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Do not open an interactive plot window",
    )
    args = parser.parse_args(argv)

    # Step 2 — file I/O and header parsing live in diagnostics.io.prp
    data = read_prp(args.prp_file)

    if args.list:
        _list_properties(data)
        return 0

    # Step 3 — resolve which column to plot
    if args.property_name:
        prop = resolve_property(data.property_names, args.property_name)
    else:
        # Interactive mode only when a human is at the terminal
        if sys.stdin.isatty():
            prop = _prompt_property(data)
        else:
            parser.error(
                "No --property given and stdin is not a TTY. "
                "Pass -p/--property or use --list."
            )

    # Step 4 — analysis + plot
    plot_property(
        data,
        prop,
        save_path=args.save,
        show=not args.no_show,
    )
    return 0


if __name__ == "__main__":
    # ``SystemExit`` lets the shell see a non-zero status on failure.
    raise SystemExit(main())
