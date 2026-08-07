#!/usr/bin/env python3
#*******************************************************************************
#  angle_distribution.py
#
#  Cassandra V2 diagnostic tool — bond-angle probability distribution from a
#  ``.xyz`` movie, using topology and harmonic parameters from the ``.mcf``.
#
#  ----------------------------
#  Why look at angle distributions?
#  ----------------------------
#  Cassandra samples configurations with a Metropolis Monte Carlo algorithm.
#  For a flexible (harmonic) bond angle, the ideal-gas reference distribution
#  at temperature T is the Boltzmann distribution of the angle potential,
#  including the geometric ``sin(θ)`` weight that comes from uniform sampling
#  in cos(θ).  Comparing the histogram from the trajectory to that ideal PDF
#  is a quick sanity check:
#
#    * Is the mean near θ₀?
#    * Is the width roughly consistent with T and K_θ?
#    * Did I use the right MCF / temperature for this run?
#
#  In a dense fluid the distribution can shift slightly relative to ideal gas
#  (intermolecular packing).  Large disagreements usually mean wrong inputs,
#  too few samples, or a bug in the run setup — not a small liquid shift.
#
#  Two analysis modes
#  ------------------
#  1. ALL molecules (default): pool angle #N from every molecule in every
#     frame.  This is the usual research mode — good statistics.
#  2. ONE molecule (``--molecule M``): only molecule ID M (1-based within
#     the selected species).  Useful for teaching and for spotting odd
#     individuals.  Needs many frames; short example movies will look noisy.
#
#  Multi-species XYZ
#  -----------------
#  Pass ``--inp`` so ``# Molecule_Files`` supplies species order and counts.
#  Select the species with ``--species N`` (default 1).  If the movie's
#  molecule counts differ from the ``.inp`` (e.g. GCMC / zero occupancy),
#  override with ``--nmols n1,n2,...``.  The positional MCF must be the MCF
#  for the selected species.
#
#  Two-panel figure
#  ----------------
#  Left:  sampled histogram (PDF) with mean and ± sample stdev marked.
#  Right: ideal-gas Boltzmann PDF at the simulation temperature (harmonic only).
#
#  Normalization of the plotted distributions (important for teaching)
#  -------------------------------------------------------------------
#  Both panels show a **probability density function** (PDF) in θ, with θ in
#  **degrees**, not a probability mass function (PMF).
#
#  Matplotlib ``hist(..., density=True)`` sets each bar height to
#
#      h_i  =  n_i / (N * Δθ)
#
#  so that the **area** of the bars integrates to one::
#
#      ∫ P(θ) dθ  =  Σ_i  h_i * Δθ  =  1
#
#  with dθ in degrees.  Consequently:
#
#    * Bar *heights* do **not** sum to 1 (unless Δθ = 1° by chance).
#    * With few samples (e.g. one molecule and a short movie) you see a
#      handful of tall, thin bars — that is normal for a density plot.
#    * Y-axis units are 1/degree.  The ideal Boltzmann curve on the right
#      is normalized the same way (∫ P(θ) dθ = 1 over degrees) so the
#      panels are directly comparable.
#
#  Ideal-gas Boltzmann PDF (harmonic angle)
#  ----------------------------------------
#  MCF stores K in K/rad² and θ₀ in degrees.  Energy form (no ½)::
#
#      U / k_B  =  K * (θ − θ₀)²     [θ in radians]
#
#  Unnormalized weight (θ in radians)::
#
#      w(θ)  ∝  sin(θ) * exp( −K (θ − θ₀)² / T )
#
#  The ``sin(θ)`` factor is the Jacobian for transforming a distribution that
#  is uniform in cos(θ) (see Src/angle_dist_pick.f90).  We evaluate w on a
#  fine θ-grid in **degrees**, then normalize with a trapezoidal integral
#  so ∫ P(θ) dθ = 1 in those degree units — matching ``density=True``.
#
#  Related modules
#  ---------------
#  ``diagnostics/io/xyz.py``  — multi-frame XYZ reader
#  ``diagnostics/io/mcf.py``  — Atom_Info / Angle_Info / T from .inp
#  ``docs/output-formats.md`` — XYZ atom order assumptions
#  ``docs/mcf-setup-workflow.md`` — MCF angle unit conventions
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#  Usage examples
#  --------------
#    python diagnostics/angle_distribution.py \\
#      Examples/NVT/dimethylether/nvt.out.xyz \\
#      Examples/NVT/dimethylether/dme.mcf \\
#      --inp Examples/NVT/dimethylether/nvt.inp --angle 1 --list
#
#    python diagnostics/angle_distribution.py \\
#      Examples/NPT/pentane/equil.out.xyz \\
#      Examples/NPT/pentane/pentane.mcf \\
#      --inp Examples/NPT/pentane/equil.inp --angle 3 --save ang3.png --no-show
#
#    python diagnostics/angle_distribution.py \\
#      equil.out.xyz pentane.mcf --temperature 336 --angle 3 --molecule 25
#
#*******************************************************************************
"""Bond-angle PDF from Cassandra ``.xyz`` + ``.mcf`` vs ideal-gas Boltzmann.

See the file header above for teaching notes on modes, units, and the PDF.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Import path (same pattern as property_plotter.py)
# ---------------------------------------------------------------------------
_DIAG_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _DIAG_ROOT.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.inp import read_temperature_from_inp
from diagnostics.io.layout import parse_nmols_arg, resolve_layout
from diagnostics.io.mcf import (
    McfAngle,
    McfData,
    read_mcf,
)
from diagnostics.io.xyz import XyzTrajectory, read_xyz


# ---------------------------------------------------------------------------
# Geometry and Boltzmann helpers
# ---------------------------------------------------------------------------

def bond_angle_degrees(coords: np.ndarray, i: int, j: int, k: int) -> float:
    """Return the angle i−j−k in degrees for one molecule's atom coordinates.

    Parameters
    ----------
    coords :
        Shape ``(n_atoms, 3)`` Cartesian coordinates of one molecule (Å).
    i, j, k :
        0-based atom indices; ``j`` is the vertex.
    """
    rij = coords[i] - coords[j]
    rkj = coords[k] - coords[j]
    nij = np.linalg.norm(rij)
    nkj = np.linalg.norm(rkj)
    if nij == 0.0 or nkj == 0.0:
        return float("nan")
    cos_th = float(np.dot(rij, rkj) / (nij * nkj))
    # Numerical clamp — floating-point can put cos slightly outside [-1, 1]
    cos_th = max(-1.0, min(1.0, cos_th))
    return float(np.degrees(np.arccos(cos_th)))


def collect_angles(
    traj: XyzTrajectory,
    layout,
    angle: McfAngle,
    *,
    species_index: int = 1,
    molecule: int | None = None,
) -> np.ndarray:
    """Collect θ (degrees) for one MCF angle over frames / molecules.

    Parameters
    ----------
    layout :
        :class:`~diagnostics.io.layout.SpeciesLayout` describing the XYZ atom
        order (one or more species).
    species_index :
        1-based species whose MCF defines ``angle``.
    molecule :
        If ``None``, use all molecules of that species.  If set, 1-based
        molecule ID within the species only.
    """
    sp = layout.block(species_index)
    if molecule is not None:
        mol_ids = [molecule]
        # Validate range early
        layout.molecule_coords(traj.coords[0], species_index, molecule)
    else:
        mol_ids = list(range(1, sp.nmol + 1))

    values: list[float] = []
    for iframe in range(traj.n_frames):
        frame = traj.coords[iframe]
        for mid in mol_ids:
            mol_xyz = layout.molecule_coords(frame, species_index, mid)
            values.append(bond_angle_degrees(mol_xyz, angle.i, angle.j, angle.k))

    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def boltzmann_angle_pdf(
    theta_deg: np.ndarray,
    *,
    k_theta: float,
    theta0_deg: float,
    temperature: float,
) -> np.ndarray:
    """Ideal-gas harmonic angle PDF on a grid of θ in **degrees**.

    Returns a probability **density** P(θ) with

        ∫ P(θ) dθ ≈ 1

    when θ is integrated in **degrees** (same convention as
    ``matplotlib.hist(..., density=True)``).  Heights therefore have units
    of 1/degree; do not expect Σ P(θ_i) = 1 on the grid.

    Internally the Boltzmann factor uses radians; only the normalization
    integral is performed in degree space so the curve matches the histogram.
    """
    if temperature <= 0.0:
        raise ValueError(f"Temperature must be positive, got {temperature}")
    if k_theta < 0.0:
        raise ValueError(f"K_theta must be non-negative, got {k_theta}")

    theta_rad = np.radians(theta_deg)
    theta0_rad = np.radians(theta0_deg)
    # P ∝ sin(θ) exp(−K (θ−θ₀)² / T); sin(θ)≥0 on (0, π)
    sin_th = np.sin(theta_rad)
    sin_th = np.clip(sin_th, 0.0, None)
    delta = theta_rad - theta0_rad
    log_w = -k_theta * delta * delta / temperature
    # Stabilize exp for stiff springs
    log_w -= np.max(log_w)
    w = sin_th * np.exp(log_w)
    # Integrate in degree space so the plotted PDF matches histogram density
    # from numpy.histogram(..., density=True) with edges in degrees.
    area = np.trapezoid(w, theta_deg)
    if area <= 0.0:
        raise ValueError("Boltzmann weight integrated to zero — check K, θ₀, T")
    return w / area


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def _list_angles(mcf: McfData) -> None:
    """Print MCF angles with indices students will pass to ``--angle``."""
    print(f"File: {mcf.path}")
    print(f"Atoms per molecule: {mcf.n_atoms}")
    print("Angles (use --angle <index>):")
    if not mcf.angles:
        print("  (none)")
        return
    for ang in mcf.angles:
        label = ang.label(mcf.atoms)
        if ang.is_harmonic:
            params = f"harmonic  K={ang.k_theta:.4g} K/rad²  θ₀={ang.theta0_deg:.2f}°"
        else:
            params = f"fixed  θ_eq={ang.theta0_deg:.2f}°"
        print(f"  {ang.index:2d}.  {label}   {params}")


def _prompt_angle(mcf: McfData) -> int:
    """Interactive choice of MCF angle index."""
    _list_angles(mcf)
    print()
    while True:
        raw = input("Angle index (empty to cancel): ").strip()
        if not raw:
            raise SystemExit("Cancelled.")
        if not raw.isdigit():
            print("Enter a positive integer angle index.")
            continue
        idx = int(raw)
        try:
            mcf.angle_by_index(idx)
            return idx
        except ValueError as exc:
            print(exc)


def plot_angle_distribution(
    samples_deg: np.ndarray,
    angle: McfAngle,
    mcf: McfData,
    *,
    temperature: float,
    n_molecules: int,
    molecule: int | None,
    species_index: int,
    n_frames: int,
    bins: int,
    title_path: str,
    save_path: Path | None = None,
    show: bool = True,
) -> None:
    """Two-panel figure: sampled PDF (left) and ideal Boltzmann PDF (right).

    Both panels use probability **densities** in θ (degrees): bar/curve
    *areas* integrate to 1, not the sum of bar heights.  See the module
    header section "Normalization of the plotted distributions".
    """
    import matplotlib.pyplot as plt

    mean = float(np.mean(samples_deg))
    std = float(np.std(samples_deg, ddof=1)) if samples_deg.size > 1 else 0.0
    label = angle.label(mcf.atoms)
    if molecule is not None:
        pop = f"species {species_index}, molecule {molecule}"
    else:
        pop = f"species {species_index}, all {n_molecules} molecules"

    # Shared x-limits: data range padded, or a window around θ₀ for empty-ish sets
    pad = max(3.0 * std, 5.0)
    x_lo = max(0.0, mean - pad)
    x_hi = min(180.0, mean + pad)
    if angle.theta0_deg < x_lo or angle.theta0_deg > x_hi:
        x_lo = max(0.0, min(x_lo, angle.theta0_deg - pad))
        x_hi = min(180.0, max(x_hi, angle.theta0_deg + pad))

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)

    # --- Left: sampled histogram ----------------------------------------
    # density=True → PDF, not probability mass:
    #   h_i = n_i / (N * Δθ)   so   Σ h_i * Δθ = 1
    # Bar heights alone do not sum to 1; y-units are 1/degree.
    # Sparse single-molecule histograms (few tall bars) are expected when N
    # is small — the *areas* still normalize.
    ax_l.hist(
        samples_deg,
        bins=bins,
        range=(x_lo, x_hi),
        density=True,
        color="C0",
        alpha=0.75,
        edgecolor="white",
        linewidth=0.5,
        label="Sampled PDF",
    )
    ax_l.axvline(mean, color="C3", linewidth=1.8, label=f"mean = {mean:.2f}°")
    if samples_deg.size > 1:
        ax_l.axvline(mean - std, color="C3", linestyle="--", linewidth=1.2, alpha=0.8)
        ax_l.axvline(
            mean + std,
            color="C3",
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
            label=f"± stdev = {std:.2f}°",
        )
    ax_l.axvline(
        angle.theta0_deg,
        color="0.35",
        linestyle=":",
        linewidth=1.5,
        label=f"θ₀ = {angle.theta0_deg:.2f}°",
    )
    ax_l.set_xlim(x_lo, x_hi)
    ax_l.set_xlabel("Angle θ (degrees)")
    # Explicit 1/deg reminds students this is a density (area = 1), not a PMF
    ax_l.set_ylabel(r"Probability density (1/deg)")
    ax_l.set_title(f"Sampled distribution\n{pop}")
    ax_l.legend(loc="best", fontsize=8)
    ax_l.grid(True, alpha=0.3)

    # --- Right: ideal Boltzmann -----------------------------------------
    # Same degree-based density normalization as the left panel (see
    # boltzmann_angle_pdf docstring and module-header notes).
    if angle.is_harmonic and angle.k_theta is not None:
        grid = np.linspace(x_lo, x_hi, 500)
        pdf = boltzmann_angle_pdf(
            grid,
            k_theta=angle.k_theta,
            theta0_deg=angle.theta0_deg,
            temperature=temperature,
        )
        ax_r.plot(grid, pdf, color="C1", linewidth=2.0, label="Ideal-gas Boltzmann")
        ax_r.fill_between(grid, pdf, alpha=0.25, color="C1")
        ax_r.axvline(
            angle.theta0_deg,
            color="0.35",
            linestyle=":",
            linewidth=1.5,
            label=f"θ₀ = {angle.theta0_deg:.2f}°",
        )
        ax_r.set_title(f"Ideal-gas Boltzmann\nT = {temperature:.2f} K")
        ax_r.legend(loc="best", fontsize=8)
    else:
        ax_r.text(
            0.5,
            0.5,
            "No Boltzmann PDF\n(fixed / non-harmonic angle)",
            ha="center",
            va="center",
            transform=ax_r.transAxes,
        )
        ax_r.set_title("Ideal-gas Boltzmann")

    ax_r.set_xlim(x_lo, x_hi)
    ax_r.set_xlabel("Angle θ (degrees)")
    ax_r.set_ylabel(r"Probability density (1/deg)")
    ax_r.grid(True, alpha=0.3)

    fig.suptitle(f"MCF angle #{angle.index}: {label}\n{title_path}", fontsize=11)
    fig.tight_layout()

    # --- terminal report ------------------------------------------------
    print()
    print(f"Angle:        #{angle.index}  {label}  ({angle.potential_type})")
    if angle.is_harmonic and angle.k_theta is not None:
        print(f"MCF params:   K = {angle.k_theta:.6g} K/rad²   θ₀ = {angle.theta0_deg:.4g}°")
    else:
        print(f"MCF params:   θ_eq = {angle.theta0_deg:.4g}°  (fixed)")
    print(f"Temperature:  {temperature:.4g} K")
    print(f"Population:   {pop}")
    print(f"Frames:       {n_frames}")
    print(f"Samples:      {samples_deg.size}")
    print(f"Mean θ:       {mean:.6g}°")
    if samples_deg.size > 1:
        print(f"Sample stdev: {std:.6g}°")
    print(f"Min / max:    {float(np.min(samples_deg)):.6g}° / {float(np.max(samples_deg)):.6g}°")

    if save_path is not None:
        save_path = Path(save_path)
        fig.savefig(save_path, dpi=150)
        print(f"Saved plot:   {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    """CLI: parse args, collect angles, two-panel plot."""
    parser = argparse.ArgumentParser(
        description=(
            "Histogram one MCF bond angle from a Cassandra .xyz movie and "
            "compare to the ideal-gas Boltzmann PDF at temperature T."
        )
    )
    parser.add_argument("xyz_file", type=Path, help="Cassandra .xyz movie file")
    parser.add_argument(
        "mcf_file",
        type=Path,
        help="MCF for the selected species (Atom_Info + Angle_Info)",
    )
    parser.add_argument(
        "-a",
        "--angle",
        type=int,
        default=None,
        metavar="N",
        help="1-based MCF angle index (from # Angle_Info)",
    )
    parser.add_argument(
        "-m",
        "--molecule",
        type=int,
        default=None,
        metavar="ID",
        help="1-based molecule ID within the species (default: all molecules)",
    )
    parser.add_argument(
        "--species",
        type=int,
        default=1,
        metavar="N",
        help="1-based species index in # Molecule_Files / XYZ order (default: 1)",
    )
    parser.add_argument(
        "--nmols",
        type=str,
        default=None,
        metavar="N1,N2,...",
        help="Override molecule counts per species (needed if .inp ≠ movie)",
    )
    parser.add_argument(
        "-T",
        "--temperature",
        type=float,
        default=None,
        metavar="K",
        help="Simulation temperature in Kelvin",
    )
    parser.add_argument(
        "--inp",
        type=Path,
        default=None,
        metavar="PATH",
        help="Cassandra .inp — T and (for multi-species) # Molecule_Files",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List MCF angles and exit",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=40,
        help="Number of histogram bins (default: 40)",
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

    if args.bins < 2:
        parser.error("--bins must be at least 2")
    if args.species < 1:
        parser.error("--species must be >= 1")

    mcf = read_mcf(args.mcf_file)

    if args.list:
        _list_angles(mcf)
        return 0

    temperature: float | None = args.temperature
    if temperature is None and args.inp is not None:
        temperature = read_temperature_from_inp(args.inp)
    if temperature is None:
        parser.error("Provide --temperature T or --inp path.inp (needed for Boltzmann panel)")

    if args.angle is not None:
        angle_idx = args.angle
    elif sys.stdin.isatty():
        angle_idx = _prompt_angle(mcf)
    else:
        parser.error("No --angle given and stdin is not a TTY. Pass -a/--angle or use --list.")

    angle = mcf.angle_by_index(angle_idx)
    traj = read_xyz(args.xyz_file)
    nmols_override = parse_nmols_arg(args.nmols) if args.nmols else None
    layout = resolve_layout(
        traj_n_atoms=traj.n_atoms,
        topology_mcf=mcf,
        inp_path=args.inp,
        species_index=args.species,
        nmols_override=nmols_override,
    )
    n_mol = layout.block(args.species).nmol

    samples = collect_angles(
        traj,
        layout,
        angle,
        species_index=args.species,
        molecule=args.molecule,
    )
    if samples.size == 0:
        raise SystemExit("No finite angle samples collected.")

    if not angle.is_harmonic:
        print(
            "Warning: selected angle is not harmonic — "
            "sampled histogram will be shown, but the Boltzmann panel is empty."
        )

    plot_angle_distribution(
        samples,
        angle,
        mcf,
        temperature=temperature,
        n_molecules=n_mol,
        molecule=args.molecule,
        species_index=args.species,
        n_frames=traj.n_frames,
        bins=args.bins,
        title_path=traj.path.name,
        save_path=args.save,
        show=not args.no_show,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
