#!/usr/bin/env python3
#*******************************************************************************
#  bond_distribution.py
#
#  Cassandra V2 **diagnostic** — check that observed bond lengths in a ``.xyz``
#  movie match the fixed lengths in the ``.mcf``.
#
#  ----------------------------
#  Why this diagnostic matters
#  ----------------------------
#  Cassandra normally uses **fixed** (constrained) bond lengths.  In a correct
#  run, every occurrence of bond #N should sit at the MCF length r₀ (Å), aside
#  from tiny scatter from coordinate output precision (e.g. ``F12.3`` →
#  ~0.001 Å).  A broad distribution or a mean shifted from r₀ is a red flag:
#
#    * wrong MCF / wrong bond index
#    * coordinates not written in MCF atom order
#    * multi-species layout mismatch (wrong ``--species`` / ``--nmols``)
#    * a bug in constraints / setup
#
#  This is a **reliability check**, not a production analysis tool.  Broader
#  structural analysis (g(r), block averages, etc.) belongs in post-processing.
#
#  Two analysis modes
#  ------------------
#  1. ALL molecules of the selected species (default).
#  2. ONE molecule (``--molecule M``).
#
#  Multi-species XYZ: ``--inp``, ``--species``, optional ``--nmols`` (same as
#  angle / dihedral diagnostics).
#
#  Two-panel figure
#  ----------------
#  Left:  sampled bond-length PDF (density in Å; areas integrate to 1).
#  Right: MCF reference — a vertical line at the fixed length r₀ (a δ-function
#         idealization of a perfect constraint).
#
#  Pass / warn criterion
#  ---------------------
#  By default we warn if max |r − r₀| exceeds ``--tol`` (default 0.01 Å).
#  With ``F12.3`` XYZ output, deviations of a few ×10⁻³ Å are normal.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#  Usage
#  -----
#    python diagnostics/bond_distribution.py \\
#      ~/CassandraV2/projects/hfc125/R125equil.out.xyz \\
#      ~/CassandraV2/projects/hfc125/r125.mcf --list
#
#    python diagnostics/bond_distribution.py \\
#      R125equil.out.xyz r125.mcf --bond 1 --save r125_bond1.png --no-show
#
#*******************************************************************************
"""Bond-length diagnostic: observed distances vs MCF fixed lengths."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_DIAG_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _DIAG_ROOT.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.layout import parse_nmols_arg, resolve_layout
from diagnostics.io.mcf import McfBond, McfData, read_mcf
from diagnostics.io.xyz import XyzTrajectory, read_xyz


def bond_length_angstrom(coords: np.ndarray, i: int, j: int) -> float:
    """Euclidean distance |r_i − r_j| in Å for one molecule."""
    return float(np.linalg.norm(coords[i] - coords[j]))


def collect_bonds(
    traj: XyzTrajectory,
    layout,
    bond: McfBond,
    *,
    species_index: int = 1,
    molecule: int | None = None,
) -> np.ndarray:
    """Collect bond lengths (Å) over frames / molecules."""
    sp = layout.block(species_index)
    if molecule is not None:
        mol_ids = [molecule]
        layout.molecule_coords(traj.coords[0], species_index, molecule)
    else:
        mol_ids = list(range(1, sp.nmol + 1))

    values: list[float] = []
    for iframe in range(traj.n_frames):
        frame = traj.coords[iframe]
        for mid in mol_ids:
            mol_xyz = layout.molecule_coords(frame, species_index, mid)
            values.append(bond_length_angstrom(mol_xyz, bond.i, bond.j))
    return np.asarray(values, dtype=float)


def _list_bonds(mcf: McfData) -> None:
    print(f"File: {mcf.path}")
    print(f"Atoms per molecule: {mcf.n_atoms}")
    print("Bonds (use --bond <index>):")
    if not mcf.bonds:
        print("  (none)")
        return
    for b in mcf.bonds:
        print(
            f"  {b.index:2d}.  {b.label(mcf.atoms)}   "
            f"{b.potential_type}  r₀ = {b.length_ang:.6g} Å"
        )


def _prompt_bond(mcf: McfData) -> int:
    _list_bonds(mcf)
    print()
    while True:
        raw = input("Bond index (empty to cancel): ").strip()
        if not raw:
            raise SystemExit("Cancelled.")
        if not raw.isdigit():
            print("Enter a positive integer bond index.")
            continue
        idx = int(raw)
        try:
            mcf.bond_by_index(idx)
            return idx
        except ValueError as exc:
            print(exc)


def plot_bond_distribution(
    samples: np.ndarray,
    bond: McfBond,
    mcf: McfData,
    *,
    n_molecules: int,
    molecule: int | None,
    species_index: int,
    n_frames: int,
    bins: int,
    tol: float,
    title_path: str,
    save_path: Path | None = None,
    show: bool = True,
) -> int:
    """Two-panel figure; return 0 if within tol, 1 if warning-level deviation."""
    import matplotlib.pyplot as plt

    r0 = bond.length_ang
    mean = float(np.mean(samples))
    std = float(np.std(samples, ddof=1)) if samples.size > 1 else 0.0
    abs_dev = np.abs(samples - r0)
    max_dev = float(np.max(abs_dev))
    rms_dev = float(np.sqrt(np.mean(abs_dev * abs_dev)))
    label = bond.label(mcf.atoms)
    if molecule is not None:
        pop = f"species {species_index}, molecule {molecule}"
    else:
        pop = f"species {species_index}, all {n_molecules} molecules"

    # Zoom around r0; still show tiny scatter clearly
    span = max(5.0 * max(std, max_dev, 1e-4), 0.02)
    x_lo = r0 - span
    x_hi = r0 + span

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(11, 4.5))

    # density=True → PDF in 1/Å; areas integrate to 1
    ax_l.hist(
        samples,
        bins=bins,
        range=(x_lo, x_hi),
        density=True,
        color="C0",
        alpha=0.75,
        edgecolor="white",
        linewidth=0.5,
        label="Sampled PDF",
    )
    ax_l.axvline(mean, color="C3", linewidth=1.8, label=f"mean = {mean:.6f} Å")
    if samples.size > 1 and std > 0.0:
        ax_l.axvline(mean - std, color="C3", linestyle="--", linewidth=1.2, alpha=0.8)
        ax_l.axvline(
            mean + std,
            color="C3",
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
            label=f"± stdev = {std:.3e} Å",
        )
    ax_l.axvline(
        r0, color="0.35", linestyle=":", linewidth=1.5, label=f"r₀ = {r0:.6g} Å"
    )
    ax_l.set_xlim(x_lo, x_hi)
    ax_l.set_xlabel("Bond length r (Å)")
    ax_l.set_ylabel(r"Probability density (1/Å)")
    ax_l.set_title(f"Sampled distribution\n{pop}")
    ax_l.legend(loc="best", fontsize=8)
    ax_l.grid(True, alpha=0.3)

    # Right: ideal fixed-bond δ at r₀
    ax_r.axvline(r0, color="C1", linewidth=2.5, label=f"MCF fixed r₀ = {r0:.6g} Å")
    ax_r.set_xlim(x_lo, x_hi)
    # Match left y-scale roughly so the spike is visible as a reference marker
    ymin, ymax = ax_l.get_ylim()
    ax_r.set_ylim(0.0, max(ymax, 1.0))
    ax_r.set_xlabel("Bond length r (Å)")
    ax_r.set_ylabel(r"Probability density (1/Å)")
    ax_r.set_title("MCF reference (fixed constraint)\nideal: δ(r − r₀)")
    ax_r.legend(loc="best", fontsize=8)
    ax_r.grid(True, alpha=0.3)
    ax_r.text(
        0.02,
        0.98,
        "Fixed bonds should sit at r₀\n"
        "(~0.001 Å scatter OK for F12.3 XYZ)",
        transform=ax_r.transAxes,
        va="top",
        fontsize=8,
        color="0.35",
    )

    fig.suptitle(f"MCF bond #{bond.index}: {label}\n{title_path}", fontsize=11)
    fig.tight_layout()

    within = max_dev <= tol
    print()
    print(f"Bond:         #{bond.index}  {label}  ({bond.potential_type})")
    print(f"MCF r₀:       {r0:.8g} Å")
    print(f"Population:   {pop}")
    print(f"Frames:       {n_frames}")
    print(f"Samples:      {samples.size}")
    print(f"Mean r:       {mean:.8g} Å")
    if samples.size > 1:
        print(f"Sample stdev: {std:.6e} Å")
    print(f"RMS |r−r₀|:   {rms_dev:.6e} Å")
    print(f"Max |r−r₀|:   {max_dev:.6e} Å")
    print(f"Min / max r:  {float(np.min(samples)):.8g} / {float(np.max(samples)):.8g} Å")
    if within:
        print(f"Check:        OK (max deviation ≤ tol = {tol:g} Å)")
        status = 0
    else:
        print(
            f"Check:        WARNING — max deviation {max_dev:.6e} Å exceeds "
            f"tol = {tol:g} Å"
        )
        status = 1

    if save_path is not None:
        save_path = Path(save_path)
        fig.savefig(save_path, dpi=150)
        print(f"Saved plot:   {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)
    return status


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic: histogram one MCF bond length from a Cassandra .xyz "
            "movie and compare to the fixed MCF length r₀."
        )
    )
    parser.add_argument("xyz_file", type=Path, help="Cassandra .xyz movie file")
    parser.add_argument(
        "mcf_file",
        type=Path,
        help="MCF for the selected species (Atom_Info + Bond_Info)",
    )
    parser.add_argument(
        "-b",
        "--bond",
        type=int,
        default=None,
        metavar="N",
        help="1-based MCF bond index (from # Bond_Info)",
    )
    parser.add_argument(
        "-m",
        "--molecule",
        type=int,
        default=None,
        metavar="ID",
        help="1-based molecule ID within the species (default: all)",
    )
    parser.add_argument(
        "--species",
        type=int,
        default=1,
        metavar="N",
        help="1-based species index (default: 1)",
    )
    parser.add_argument(
        "--nmols",
        type=str,
        default=None,
        metavar="N1,N2,...",
        help="Override molecule counts per species",
    )
    parser.add_argument(
        "--inp",
        type=Path,
        default=None,
        metavar="PATH",
        help="Cassandra .inp — multi-species # Molecule_Files if needed",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.01,
        metavar="ANG",
        help="Warn if max |r−r₀| exceeds this (Å; default: 0.01)",
    )
    parser.add_argument("--list", action="store_true", help="List MCF bonds and exit")
    parser.add_argument("--bins", type=int, default=40, help="Histogram bins (default: 40)")
    parser.add_argument("--save", type=Path, default=None, metavar="PATH", help="Save figure")
    parser.add_argument("--no-show", action="store_true", help="Do not open a plot window")
    args = parser.parse_args(argv)

    if args.bins < 2:
        parser.error("--bins must be at least 2")
    if args.species < 1:
        parser.error("--species must be >= 1")
    if args.tol < 0.0:
        parser.error("--tol must be >= 0")

    mcf = read_mcf(args.mcf_file)
    if args.list:
        _list_bonds(mcf)
        return 0

    if args.bond is not None:
        bond_idx = args.bond
    elif sys.stdin.isatty():
        bond_idx = _prompt_bond(mcf)
    else:
        parser.error("No --bond given and stdin is not a TTY.")

    bond = mcf.bond_by_index(bond_idx)
    if not bond.is_fixed:
        print(
            f"Warning: bond #{bond.index} type is '{bond.potential_type}' "
            "(reference panel still marks the MCF length parameter)."
        )

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

    samples = collect_bonds(
        traj,
        layout,
        bond,
        species_index=args.species,
        molecule=args.molecule,
    )
    if samples.size == 0:
        raise SystemExit("No bond-length samples collected.")

    return plot_bond_distribution(
        samples,
        bond,
        mcf,
        n_molecules=n_mol,
        molecule=args.molecule,
        species_index=args.species,
        n_frames=traj.n_frames,
        bins=args.bins,
        tol=args.tol,
        title_path=traj.path.name,
        save_path=args.save,
        show=not args.no_show,
    )


if __name__ == "__main__":
    raise SystemExit(main())
