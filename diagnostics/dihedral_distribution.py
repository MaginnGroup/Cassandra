#!/usr/bin/env python3
#*******************************************************************************
#  dihedral_distribution.py
#
#  Cassandra V2 diagnostic tool — dihedral (torsion) probability density from a
#  ``.xyz`` movie, using topology and dihedral parameters from the ``.mcf``.
#
#  ----------------------------
#  Why look at dihedral distributions?
#  ----------------------------
#  The MCF dihedral potential is only *part* of what shapes the torsion
#  probability density in a condensed-phase simulation.  Cassandra also has
#  intramolecular 1–4 nonbonded (LJ / Coulomb) interactions, and the torsion
#  is coupled to other degrees of freedom.  Comparing the sampled PDF to the
#  **ideal-gas Boltzmann curve of the bare dihedral potential** is therefore
#  especially instructive:
#
#    * Agreement near the minima → the MCF dihedral term dominates there.
#    * Systematic shifts / suppressed barriers → 1–4 and packing effects
#      (exactly what students should expect — not a bug by itself).
#    * Completely wrong shape → wrong MCF, wrong dihedral index, or T.
#
#  Two analysis modes
#  ------------------
#  1. ALL molecules of the selected species (default).
#  2. ONE molecule (``--molecule M``, 1-based within the species).
#
#  Multi-species XYZ
#  -----------------
#  Pass ``--inp`` so ``# Molecule_Files`` supplies species order and counts.
#  Select the species with ``--species N``.  Override counts with
#  ``--nmols n1,n2,...`` when the movie differs from the ``.inp``.
#
#  Two-panel figure
#  ----------------
#  Left:  sampled histogram (PDF) with mean and ± sample stdev.
#  Right: ideal-gas Boltzmann PDF from the MCF dihedral potential only.
#
#  Normalization (same teaching point as angle_distribution.py)
#  ------------------------------------------------------------
#  Both panels are probability **densities** in φ (degrees)::
#
#      ∫ P(φ) dφ  =  1     (dφ in degrees; y-units = 1/deg)
#
#  Matplotlib ``hist(..., density=True)`` uses bar heights
#  h_i = n_i / (N Δφ), so **areas** sum to 1, not the heights.
#
#  Ideal-gas Boltzmann (bare dihedral potential)
#  ---------------------------------------------
#  P(φ) ∝ exp( −U_dihedral(φ) / (R T) ) for OPLS/CHARMM/RB (U in kJ/mol),
#  or exp( −K (φ−φ₀)² / T ) for harmonic (K in K/rad²).  No geometric
#  Jacobian (unlike bond angles).  φ ∈ (−180°, 180°] matches Cassandra's
#  SIGN(ACOS(cos φ), …) convention (Src/internal_coordinate_routines.f90).
#
#  Energy forms (MCF units; see docs/mcf-setup-workflow.md)
#  -------------------------------------------------------
#  OPLS:   U = a0 + a1(1+cos φ) + a2(1−cos 2φ) + a3(1+cos 3φ)   [kJ/mol]
#  CHARMM: U = a0 (1 + cos(n φ − δ))                            [kJ/mol]
#  RB:     U = Σ_{k=0}^{5} c_k [cos φ]^k                        [kJ/mol]
#  harmonic: U/k_B = K (φ − φ₀)²                                [K]
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#  Usage examples
#  --------------
#    python diagnostics/dihedral_distribution.py \\
#      Examples/NPT/pentane/equil.out.xyz \\
#      Examples/NPT/pentane/pentane.mcf \\
#      --inp Examples/NPT/pentane/equil.inp --dihedral 1 --list
#
#    python diagnostics/dihedral_distribution.py \\
#      equil.out.xyz pentane.mcf --inp equil.inp --dihedral 2 \\
#      --save dih2.png --no-show
#
#    # R125 (CHARMM dihedrals), all molecules
#    python diagnostics/dihedral_distribution.py \\
#      R125equil.out.xyz r125.mcf --inp r125-equil.inp --dihedral 3
#
#*******************************************************************************
"""Dihedral PDF from Cassandra ``.xyz`` + ``.mcf`` vs bare-potential Boltzmann.

See the file header for teaching notes (1–4 effects, normalization, units).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

_DIAG_ROOT = Path(__file__).resolve().parent
_REPO_ROOT = _DIAG_ROOT.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.inp import read_temperature_from_inp
from diagnostics.io.layout import parse_nmols_arg, resolve_layout
from diagnostics.io.mcf import McfData, McfDihedral, read_mcf
from diagnostics.io.xyz import XyzTrajectory, read_xyz

# Gas constant for OPLS / CHARMM / RB energies stored in kJ/mol
_R_KJ_MOL_K = 8.314462618e-3


# ---------------------------------------------------------------------------
# Geometry (Cassandra convention)
# ---------------------------------------------------------------------------

def dihedral_degrees(
    coords: np.ndarray,
    i: int,
    j: int,
    k: int,
    l: int,
) -> float:
    """Return the torsion φ for atoms i−j−k−l in degrees, ∈ (−180, 180].

    Matches ``Get_Dihedral_Angle`` / ``Get_Dihedral_Angle_COS`` in
    ``Src/internal_coordinate_routines.f90``:

      r12 = r_i − r_j
      r32 = r_k − r_j
      r34 = r_k − r_l
      m = r12 × r32   (normal to plane i-j-k)
      n = r32 × r34   (normal to plane j-k-l)
      cos φ = (m·n) / (|m||n|)
      φ = sign(acos(cos φ), r12·n)
    """
    r12 = coords[i] - coords[j]
    r32 = coords[k] - coords[j]
    r34 = coords[k] - coords[l]
    m = np.cross(r12, r32)
    n = np.cross(r32, r34)
    msq = float(np.dot(m, m))
    nsq = float(np.dot(n, n))
    if msq == 0.0 or nsq == 0.0:
        return float("nan")
    cosphi = float(np.dot(m, n) / np.sqrt(msq * nsq))
    cosphi = max(-1.0, min(1.0, cosphi))
    r12dn = float(np.dot(r12, n))
    phi = float(np.arccos(cosphi))
    if r12dn < 0.0:
        phi = -phi
    return float(np.degrees(phi))


def collect_dihedrals(
    traj: XyzTrajectory,
    layout,
    dihedral: McfDihedral,
    *,
    species_index: int = 1,
    molecule: int | None = None,
) -> np.ndarray:
    """Collect φ (degrees) for one MCF dihedral over frames / molecules."""
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
            values.append(
                dihedral_degrees(
                    mol_xyz, dihedral.i, dihedral.j, dihedral.k, dihedral.l
                )
            )
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def dihedral_energy_kj_per_mol(dihedral: McfDihedral, phi_rad: np.ndarray) -> np.ndarray:
    """Bare MCF dihedral energy U(φ) in kJ/mol (vectorized over φ).

    Harmonic dihedrals are returned as U = R * K * (φ−φ₀)² so that
    U/(R T) = K (φ−φ₀)² / T with K in K/rad².
    """
    p = dihedral.params
    t = dihedral.potential_type
    phi = np.asarray(phi_rad, dtype=float)

    if t == "opls":
        a0, a1, a2, a3 = p[0], p[1], p[2], p[3]
        return (
            a0
            + a1 * (1.0 + np.cos(phi))
            + a2 * (1.0 - np.cos(2.0 * phi))
            + a3 * (1.0 + np.cos(3.0 * phi))
        )
    if t == "charmm":
        a0, n, delta_deg = p[0], p[1], p[2]
        delta = np.radians(delta_deg)
        return a0 * (1.0 + np.cos(n * phi - delta))
    if t in {"rb", "ryckaert-bellemans"}:
        c = list(p[:6]) + [0.0] * max(0, 6 - len(p))
        cosphi = np.cos(phi)
        u = np.zeros_like(phi)
        pow_c = np.ones_like(phi)
        for ck in c[:6]:
            u = u + ck * pow_c
            pow_c = pow_c * cosphi
        return u
    if t == "harmonic":
        k_phi, phi0_deg = p[0], p[1]
        # Convert to kJ/mol via R so the Boltzmann factor is consistent
        dphi = phi - np.radians(phi0_deg)
        return _R_KJ_MOL_K * k_phi * dphi * dphi
    raise ValueError(f"No energy function for dihedral type '{t}'")


def boltzmann_dihedral_pdf(
    phi_deg: np.ndarray,
    dihedral: McfDihedral,
    *,
    temperature: float,
) -> np.ndarray:
    """Ideal-gas PDF for the bare dihedral potential; φ in degrees.

    Normalized so ∫ P(φ) dφ = 1 with φ in degrees (matches density=True).
    """
    if temperature <= 0.0:
        raise ValueError(f"Temperature must be positive, got {temperature}")
    if not dihedral.has_boltzmann:
        raise ValueError(
            f"Dihedral type '{dihedral.potential_type}' has no Boltzmann PDF"
        )

    phi_rad = np.radians(phi_deg)
    u = dihedral_energy_kj_per_mol(dihedral, phi_rad)
    log_w = -u / (_R_KJ_MOL_K * temperature)
    log_w -= np.max(log_w)
    w = np.exp(log_w)
    area = np.trapezoid(w, phi_deg)
    if area <= 0.0:
        raise ValueError("Boltzmann weight integrated to zero — check params / T")
    return w / area


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------

def _list_dihedrals(mcf: McfData) -> None:
    print(f"File: {mcf.path}")
    print(f"Atoms per molecule: {mcf.n_atoms}")
    print("Dihedrals (use --dihedral <index>):")
    if not mcf.dihedrals:
        print("  (none)")
        return
    for dih in mcf.dihedrals:
        print(f"  {dih.index:2d}.  {dih.label(mcf.atoms)}   {dih.param_summary()}")


def _prompt_dihedral(mcf: McfData) -> int:
    _list_dihedrals(mcf)
    print()
    while True:
        raw = input("Dihedral index (empty to cancel): ").strip()
        if not raw:
            raise SystemExit("Cancelled.")
        if not raw.isdigit():
            print("Enter a positive integer dihedral index.")
            continue
        idx = int(raw)
        try:
            mcf.dihedral_by_index(idx)
            return idx
        except ValueError as exc:
            print(exc)


def plot_dihedral_distribution(
    samples_deg: np.ndarray,
    dihedral: McfDihedral,
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
    """Two-panel figure: sampled PDF (left) and bare-potential Boltzmann (right)."""
    import matplotlib.pyplot as plt

    mean = float(np.mean(samples_deg))
    std = float(np.std(samples_deg, ddof=1)) if samples_deg.size > 1 else 0.0
    label = dihedral.label(mcf.atoms)
    if molecule is not None:
        pop = f"species {species_index}, molecule {molecule}"
    else:
        pop = f"species {species_index}, all {n_molecules} molecules"

    # Full torsion range (−180, 180] — natural domain for Cassandra φ
    x_lo, x_hi = -180.0, 180.0

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)

    # density=True → PDF: h_i = n_i/(N Δφ); areas integrate to 1 (1/deg)
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
    ax_l.set_xlim(x_lo, x_hi)
    ax_l.set_xlabel("Dihedral φ (degrees)")
    ax_l.set_ylabel(r"Probability density (1/deg)")
    ax_l.set_title(f"Sampled distribution\n{pop}")
    ax_l.legend(loc="best", fontsize=8)
    ax_l.grid(True, alpha=0.3)

    if dihedral.has_boltzmann:
        grid = np.linspace(x_lo, x_hi, 721)
        pdf = boltzmann_dihedral_pdf(grid, dihedral, temperature=temperature)
        ax_r.plot(grid, pdf, color="C1", linewidth=2.0, label="Bare dihedral Boltzmann")
        ax_r.fill_between(grid, pdf, alpha=0.25, color="C1")
        ax_r.set_title(
            f"Ideal-gas Boltzmann (MCF dihedral only)\nT = {temperature:.2f} K"
        )
        ax_r.legend(loc="best", fontsize=8)
        ax_r.text(
            0.02,
            0.98,
            "Expect shifts from 1–4 / packing",
            transform=ax_r.transAxes,
            va="top",
            fontsize=8,
            color="0.35",
        )
    else:
        ax_r.text(
            0.5,
            0.5,
            f"No Boltzmann PDF\n(type: {dihedral.potential_type})",
            ha="center",
            va="center",
            transform=ax_r.transAxes,
        )
        ax_r.set_title("Ideal-gas Boltzmann")

    ax_r.set_xlim(x_lo, x_hi)
    ax_r.set_xlabel("Dihedral φ (degrees)")
    ax_r.set_ylabel(r"Probability density (1/deg)")
    ax_r.grid(True, alpha=0.3)

    fig.suptitle(
        f"MCF dihedral #{dihedral.index}: {label}\n{title_path}", fontsize=11
    )
    fig.tight_layout()

    print()
    print(f"Dihedral:     #{dihedral.index}  {label}  ({dihedral.potential_type})")
    print(f"MCF params:   {dihedral.param_summary()}")
    print(f"Temperature:  {temperature:.4g} K")
    print(f"Population:   {pop}")
    print(f"Frames:       {n_frames}")
    print(f"Samples:      {samples_deg.size}")
    print(f"Mean φ:       {mean:.6g}°")
    if samples_deg.size > 1:
        print(f"Sample stdev: {std:.6g}°")
    print(
        f"Min / max:    {float(np.min(samples_deg)):.6g}° / "
        f"{float(np.max(samples_deg)):.6g}°"
    )
    print(
        "Note: sampled PDF includes 1–4 nonbonded / packing; "
        "right panel is the bare MCF dihedral potential only."
    )

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
            "Histogram one MCF dihedral from a Cassandra .xyz movie and "
            "compare to the ideal-gas Boltzmann PDF of the bare dihedral potential."
        )
    )
    parser.add_argument("xyz_file", type=Path, help="Cassandra .xyz movie file")
    parser.add_argument(
        "mcf_file",
        type=Path,
        help="MCF for the selected species (Atom_Info + Dihedral_Info)",
    )
    parser.add_argument(
        "-d",
        "--dihedral",
        type=int,
        default=None,
        metavar="N",
        help="1-based MCF dihedral index (from # Dihedral_Info)",
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
        help="1-based species index in # Molecule_Files / XYZ order (default: 1)",
    )
    parser.add_argument(
        "--nmols",
        type=str,
        default=None,
        metavar="N1,N2,...",
        help="Override molecule counts per species",
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
    parser.add_argument("--list", action="store_true", help="List MCF dihedrals and exit")
    parser.add_argument("--bins", type=int, default=72, help="Histogram bins (default: 72)")
    parser.add_argument("--save", type=Path, default=None, metavar="PATH", help="Save figure")
    parser.add_argument("--no-show", action="store_true", help="Do not open a plot window")
    args = parser.parse_args(argv)

    if args.bins < 2:
        parser.error("--bins must be at least 2")
    if args.species < 1:
        parser.error("--species must be >= 1")

    mcf = read_mcf(args.mcf_file)
    if args.list:
        _list_dihedrals(mcf)
        return 0

    temperature: float | None = args.temperature
    if temperature is None and args.inp is not None:
        temperature = read_temperature_from_inp(args.inp)
    if temperature is None:
        parser.error("Provide --temperature T or --inp path.inp")

    if args.dihedral is not None:
        dih_idx = args.dihedral
    elif sys.stdin.isatty():
        dih_idx = _prompt_dihedral(mcf)
    else:
        parser.error("No --dihedral given and stdin is not a TTY.")

    dihedral = mcf.dihedral_by_index(dih_idx)
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

    samples = collect_dihedrals(
        traj,
        layout,
        dihedral,
        species_index=args.species,
        molecule=args.molecule,
    )
    if samples.size == 0:
        raise SystemExit("No finite dihedral samples collected.")

    if not dihedral.has_boltzmann:
        print(
            "Warning: selected dihedral has no Boltzmann overlay "
            f"(type={dihedral.potential_type})."
        )

    plot_dihedral_distribution(
        samples,
        dihedral,
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
