#!/usr/bin/env python3
"""Plot a Cassandra ``.prp`` property vs MC sweep/step.

Shows the instantaneous series and a running (cumulative) average, and prints
the final mean to the terminal.

Examples
--------
Interactive property choice::

    python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp

Named property::

    python diagnostics/property_plotter.py nvt.out.prp -p Energy_Total
    python diagnostics/property_plotter.py nvt.out.prp -p enthalpy --save H.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# Allow running as ``python diagnostics/property_plotter.py`` from repo root
# or from elsewhere without installing a package.
_DIAG_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _DIAG_ROOT.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.prp import PrpData, read_prp, resolve_property


def _list_properties(data: PrpData) -> None:
    print(f"File: {data.path}")
    print(f"Title: {data.title}")
    print(f"Step column: {data.step_label}  ({data.n_rows} points)")
    print("Properties:")
    for i, (name, unit) in enumerate(zip(data.property_names, data.property_units), start=1):
        unit_str = f"  [{unit}]" if unit else ""
        print(f"  {i:2d}. {name}{unit_str}")


def _prompt_property(data: PrpData) -> str:
    _list_properties(data)
    print()
    while True:
        raw = input("Property to plot (name or number, empty to cancel): ").strip()
        if not raw:
            raise SystemExit("Cancelled.")
        if raw.isdigit():
            idx = int(raw)
            if 1 <= idx <= data.n_properties:
                return data.property_names[idx - 1]
            print(f"Enter a number between 1 and {data.n_properties}.")
            continue
        try:
            return resolve_property(data.property_names, raw)
        except ValueError as exc:
            print(exc)


def _running_average(y: np.ndarray) -> np.ndarray:
    """Cumulative mean: avg[i] = mean(y[0:i+1])."""
    return np.cumsum(y, dtype=float) / np.arange(1, y.size + 1, dtype=float)


def plot_property(
    data: PrpData,
    property_name: str,
    *,
    save_path: Path | None = None,
    show: bool = True,
) -> float:
    """Plot instantaneous + running average; return final mean."""
    import matplotlib.pyplot as plt

    y = data.column(property_name)
    run_avg = _running_average(y)
    final_mean = float(run_avg[-1])
    final_std = float(np.std(y, ddof=1)) if y.size > 1 else 0.0
    unit = data.unit_for(property_name)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(data.steps, y, "o-", markersize=3, linewidth=1.0, label="Instantaneous")
    ax.plot(data.steps, run_avg, "-", linewidth=2.0, label="Running average")
    ax.set_xlabel(data.step_label)
    ylabel = property_name if not unit else f"{property_name} {unit}"
    ax.set_ylabel(ylabel)
    ax.set_title(f"{property_name} vs {data.step_label}\n{data.path.name}")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

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

    data = read_prp(args.prp_file)

    if args.list:
        _list_properties(data)
        return 0

    if args.property_name:
        prop = resolve_property(data.property_names, args.property_name)
    else:
        if sys.stdin.isatty():
            prop = _prompt_property(data)
        else:
            parser.error(
                "No --property given and stdin is not a TTY. "
                "Pass -p/--property or use --list."
            )

    plot_property(
        data,
        prop,
        save_path=args.save,
        show=not args.no_show,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
