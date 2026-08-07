#!/usr/bin/env python3
#*******************************************************************************
#  msd_com.py
#
#  Cassandra V2 **diagnostic** — mean-squared displacement of molecular centers
#  of mass vs MC sweep/step (from a ``.xyz`` movie).
#
#  ----------------------------
#  IMPORTANT teaching caveats
#  ----------------------------
#  This is **not** molecular-dynamics diffusion analysis.
#
#  * The independent variable is an MC sweep or step index, **not** physical
#    time.  There is no Einstein relation ``MSD = 6 D t`` unless you invent a
#    fictitious mapping from sweeps → time (not done here).
#  * Monte Carlo moves (translate, rotate, CBMC regrowth, volume, …) are
#    **not** Newtonian trajectories.  Regrowth can jump a COM by a large
#    distance in one accepted move.
#  * The plot is an **exploration / sampling** diagnostic: are molecules
#    moving through the box, or stuck after a bad start?
#
#  Supported ensembles for this tool: **NVT and NPT only** (fixed N).
#  GCMC / GEMC (changing N or identity across boxes) are refused — molecule
#  tracking is ambiguous there.
#
#  Algorithm
#  ---------
#  1. For each frame, compute each molecule's COM using MCF atomic masses.
#  2. Unwrap COMs across frames with the minimum-image convention (MIC) so
#     that a molecule crossing a periodic boundary does not look like a jump
#     of nearly one box length.
#  3. MSD(frame k) = (1/N_mol) Σ_m |r_m(k) − r_m(0)|^2   (Å²)
#     or, with ``--molecule M``, just |r_M(k) − r_M(0)|^2 for that molecule.
#
#  Box sizes
#  ---------
#  * NPT: companion ``.H`` file required (cell matrix per frame).
#  * NVT: ``.H`` if present; else cubic edge from ``# Box_Info`` in ``--inp``.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#  Usage
#  -----
#    python diagnostics/msd_com.py \\
#      Examples/NPT/pentane/equil.out.xyz \\
#      Examples/NPT/pentane/pentane.mcf \\
#      --inp Examples/NPT/pentane/equil.inp --save msd.png --no-show
#
#*******************************************************************************
"""COM mean-squared displacement vs MC sweep (NVT/NPT exploration diagnostic)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_DIAG_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _DIAG_ROOT.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.hfile import cubic_cell, read_h
from diagnostics.io.inp import (
    read_cubic_box_length_from_inp,
    read_sim_type_from_inp,
)
from diagnostics.io.layout import parse_nmols_arg, resolve_layout
from diagnostics.io.mcf import McfData, read_mcf
from diagnostics.io.xyz import XyzTrajectory, read_xyz

# Ensembles allowed for COM tracking with stable molecule IDs
_ALLOWED_SIM = {
    "nvt",
    "nvt_mc",
    "nvtmc",
    "npt",
    "npt_mc",
    "nptmc",
}
_FORBIDDEN_HINT = {
    "gcmc",
    "gemc",
    "gemc_nvt",
    "gemc_npt",
}


def _classify_ensemble(sim_type: str) -> str:
    """Return ``'nvt'``, ``'npt'``, or raise for unsupported types."""
    s = sim_type.lower().strip()
    if s in _FORBIDDEN_HINT or "gcmc" in s or "gemc" in s:
        raise ValueError(
            f"Sim_Type '{sim_type}' is not supported by msd_com "
            "(GCMC/GEMC change N or move molecules between boxes). "
            "Use NVT or NPT only."
        )
    if s in {"nvt", "nvt_mc", "nvtmc"} or s.startswith("nvt"):
        return "nvt"
    if s in {"npt", "npt_mc", "nptmc"} or s.startswith("npt"):
        return "npt"
    if s in _ALLOWED_SIM:
        return "npt" if "npt" in s else "nvt"
    raise ValueError(
        f"Unrecognized Sim_Type '{sim_type}'. "
        "msd_com supports nvt_mc / npt_mc only."
    )


def molecule_com(coords: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Mass-weighted COM of one molecule; ``coords`` shape ``(n_atoms, 3)``."""
    m = masses.reshape(-1, 1)
    return np.sum(coords * m, axis=0) / float(np.sum(masses))


def mic_delta(delta: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Minimum-image ``delta`` for an orthorhombic / cubic diagonal cell.

    For a general cell matrix H (columns = cell vectors), fractional
    coordinates are ``s = H^{-1} r``, then MIC in s, then ``r = H s``.
    Cassandra cubic/ortho boxes are diagonal in practice.
    """
    h_inv = np.linalg.inv(cell)
    s = h_inv @ delta
    s -= np.rint(s)
    return cell @ s


def collect_coms(
    traj: XyzTrajectory,
    layout,
    masses: np.ndarray,
    *,
    species_index: int,
    molecule: int | None,
) -> np.ndarray:
    """Return COM array shape ``(n_frames, n_mol, 3)`` for the selected set."""
    sp = layout.block(species_index)
    if molecule is not None:
        mol_ids = [molecule]
        layout.molecule_coords(traj.coords[0], species_index, molecule)
    else:
        mol_ids = list(range(1, sp.nmol + 1))

    n_frames = traj.n_frames
    n_mol = len(mol_ids)
    coms = np.empty((n_frames, n_mol, 3), dtype=float)
    for iframe in range(n_frames):
        frame = traj.coords[iframe]
        for im, mid in enumerate(mol_ids):
            mol_xyz = layout.molecule_coords(frame, species_index, mid)
            coms[iframe, im] = molecule_com(mol_xyz, masses)
    return coms


def unwrap_coms(coms: np.ndarray, cells: np.ndarray) -> np.ndarray:
    """Unwrap COM trajectories with MIC frame-to-frame displacements.

    ``coms``: ``(n_frames, n_mol, 3)`` wrapped coordinates.
    ``cells``: ``(n_frames, 3, 3)`` cell matrices (Å).
    """
    n_frames, n_mol, _ = coms.shape
    if cells.shape[0] != n_frames:
        raise ValueError(
            f"Need one cell per XYZ frame ({n_frames}), got {cells.shape[0]}"
        )
    unwrapped = np.empty_like(coms)
    unwrapped[0] = coms[0]
    for f in range(1, n_frames):
        # Use the later frame's cell for MIC (NPT: box may have changed)
        cell = cells[f]
        for m in range(n_mol):
            d = mic_delta(coms[f, m] - coms[f - 1, m], cell)
            unwrapped[f, m] = unwrapped[f - 1, m] + d
    return unwrapped


def msd_from_unwrapped(unwrapped: np.ndarray) -> np.ndarray:
    """MSD vs frame relative to frame 0; average over molecules.

    Returns shape ``(n_frames,)`` in Å².
    """
    dr = unwrapped - unwrapped[0:1]
    sq = np.sum(dr * dr, axis=2)  # (n_frames, n_mol)
    return np.mean(sq, axis=1)


def resolve_cells(
    *,
    ensemble: str,
    n_frames: int,
    xyz_path: Path,
    h_path: Path | None,
    inp_path: Path | None,
) -> np.ndarray:
    """Return ``(n_frames, 3, 3)`` cell matrices for MIC unwrapping."""
    # Prefer explicit --h, else same basename as xyz
    candidates: list[Path] = []
    if h_path is not None:
        candidates.append(h_path)
    else:
        auto = xyz_path.with_suffix(".H")
        # Cassandra often uses run.out.xyz → run.out.H
        if xyz_path.name.endswith(".out.xyz"):
            candidates.append(xyz_path.with_name(xyz_path.name.replace(".xyz", ".H")))
        candidates.append(auto)

    h_file = next((p for p in candidates if p.is_file()), None)

    if ensemble == "npt":
        if h_file is None:
            raise ValueError(
                "NPT MSD requires the companion .H file (cell matrix per frame). "
                f"Looked for: {', '.join(str(p) for p in candidates)}. "
                "Pass --h PATH."
            )
        htraj = read_h(h_file)
        if htraj.n_frames != n_frames:
            raise ValueError(
                f".H has {htraj.n_frames} frames but XYZ has {n_frames}. "
                "They must share coord_freq / the same movie."
            )
        print(f"Box source:   {h_file} (per-frame cells)")
        return htraj.cells

    # NVT
    if h_file is not None:
        htraj = read_h(h_file)
        if htraj.n_frames != n_frames:
            raise ValueError(
                f".H has {htraj.n_frames} frames but XYZ has {n_frames}."
            )
        print(f"Box source:   {h_file} (per-frame cells)")
        return htraj.cells

    if inp_path is None:
        raise ValueError(
            "NVT without a .H file needs --inp to read cubic # Box_Info."
        )
    L = read_cubic_box_length_from_inp(inp_path)
    print(f"Box source:   # Box_Info cubic L = {L:g} Å (constant)")
    cell = cubic_cell(L)
    return np.broadcast_to(cell, (n_frames, 3, 3)).copy()


def plot_msd(
    steps: np.ndarray,
    msd: np.ndarray,
    *,
    pop: str,
    title_path: str,
    ensemble: str,
    save_path: Path | None,
    show: bool,
) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(steps, msd, "o-", markersize=3, linewidth=1.2, label="MSD")
    ax.set_xlabel("MC sweep / step (from XYZ)")
    ax.set_ylabel(r"MSD ($\mathrm{\AA}^2$)")
    ax.set_title(
        f"COM MSD vs MC index ({ensemble.upper()})\n{pop}\n{title_path}"
    )
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    ax.text(
        0.02,
        0.98,
        "Not MD diffusion — exploration vs MC index",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout()

    print()
    print(f"Ensemble:     {ensemble.upper()} (fixed N)")
    print(f"Population:   {pop}")
    print(f"Frames:       {msd.size}")
    print(f"Final MSD:    {float(msd[-1]):.6g} Å²")
    print(f"Final RMS:    {float(np.sqrt(msd[-1])):.6g} Å")
    if msd.size > 1 and float(steps[-1] - steps[0]) != 0.0:
        # Crude slope of last half — teaching only, not a diffusivity
        mid = msd.size // 2
        slope = float(
            (msd[-1] - msd[mid]) / (steps[-1] - steps[mid])
        )
        print(f"Slope (2nd half, Å² / MC unit): {slope:.6g}  [not a D]")

    if save_path is not None:
        save_path = Path(save_path)
        fig.savefig(save_path, dpi=150)
        print(f"Saved plot:   {save_path}")

    if show:
        plt.show()
    else:
        plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnostic: molecular COM mean-squared displacement vs MC "
            "sweep/step (NVT/NPT only). Not an MD diffusion coefficient."
        )
    )
    parser.add_argument("xyz_file", type=Path, help="Cassandra .xyz movie")
    parser.add_argument(
        "mcf_file",
        type=Path,
        help="MCF for the selected species (masses for COM)",
    )
    parser.add_argument(
        "--inp",
        type=Path,
        required=True,
        metavar="PATH",
        help="Cassandra .inp (Sim_Type; Box_Info for NVT without .H)",
    )
    parser.add_argument(
        "--h",
        type=Path,
        default=None,
        metavar="PATH",
        help="Companion .H file (required for NPT; optional for NVT)",
    )
    parser.add_argument(
        "-m",
        "--molecule",
        type=int,
        default=None,
        metavar="ID",
        help="1-based molecule ID only (default: average over all)",
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
    parser.add_argument("--save", type=Path, default=None, metavar="PATH")
    parser.add_argument("--no-show", action="store_true")
    args = parser.parse_args(argv)

    if args.species < 1:
        parser.error("--species must be >= 1")

    sim_type = read_sim_type_from_inp(args.inp)
    ensemble = _classify_ensemble(sim_type)

    mcf = read_mcf(args.mcf_file)
    masses = np.asarray([a.mass for a in mcf.atoms], dtype=float)
    if np.any(masses <= 0.0):
        raise SystemExit("MCF has non-positive atomic masses — cannot compute COM.")

    traj = read_xyz(args.xyz_file)
    if traj.n_frames < 2:
        raise SystemExit("Need at least 2 XYZ frames for MSD.")

    nmols_override = parse_nmols_arg(args.nmols) if args.nmols else None
    layout = resolve_layout(
        traj_n_atoms=traj.n_atoms,
        topology_mcf=mcf,
        inp_path=args.inp,
        species_index=args.species,
        nmols_override=nmols_override,
    )
    n_mol = layout.block(args.species).nmol

    cells = resolve_cells(
        ensemble=ensemble,
        n_frames=traj.n_frames,
        xyz_path=traj.path,
        h_path=args.h,
        inp_path=args.inp,
    )

    coms = collect_coms(
        traj,
        layout,
        masses,
        species_index=args.species,
        molecule=args.molecule,
    )
    unwrapped = unwrap_coms(coms, cells)
    msd = msd_from_unwrapped(unwrapped)

    if args.molecule is not None:
        pop = f"species {args.species}, molecule {args.molecule}"
    else:
        pop = f"species {args.species}, average over {n_mol} molecules"

    # Prefer XYZ step labels; fall back to frame index
    steps = traj.steps.astype(float)
    if not np.all(np.isfinite(steps)):
        steps = np.arange(traj.n_frames, dtype=float)

    plot_msd(
        steps,
        msd,
        pop=pop,
        title_path=traj.path.name,
        ensemble=ensemble,
        save_path=args.save,
        show=not args.no_show,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
