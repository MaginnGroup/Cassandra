#!/usr/bin/env python3
#*******************************************************************************
#  npt_restart.py
#
#  Build an NVT (or a GCMC template) from a finished NPT run.
#
#  What problem this solves
#  ------------------------
#  The last frame of an NPT run has *some* volume, not the average volume.
#  An NVT simulation needs one fixed box.  Using the last frame unchanged
#  would freeze in a random fluctuation.  This script:
#
#    1. Reads the equilibrated average volume, <V>, from block_average.py
#       (or from a JSON file that script already wrote).
#    2. Reads the single-frame ``*.restart.xyz`` and ``*.restart.H``.
#       It refuses the multi-frame movie, because ``read_config`` would
#       load the first frame of that movie.
#    3. Maps each molecular center of mass from its old fractional coordinate
#       into the new cubic cell, then translates the whole molecule.  It
#       reconstructs molecules across PBC before finding each COM, so bonds,
#       angles, and dihedrals are not scaled or otherwise changed.
#    4. Writes a new ``.inp``, ``.xyz``, and ``.H`` in a directory you
#       choose.  The original NPT directory is not modified.
#
#  Which density the new box has
#  ---------------------------
#  block_average.py prints two densities, and they are not interchangeable.
#
#  Mean instantaneous density, <rho>
#      The average of the Mass_Density column.  Density is proportional to
#      1/V, so this average is proportional to <1/V>.
#
#  Density at the mean volume
#      (total mass) / <V>, proportional to 1/<V>.  Because <1/V> is larger
#      than 1/<V> whenever V fluctuates, <rho> is slightly larger than the
#      density of the mean-volume box.  The gap is second order in the
#      volume fluctuations and is usually smaller than the uncertainty on
#      <rho>.  At 243.15 K for the HFC-125 production run it was
#      1433.71 versus 1433.57 kg/m^3.
#
#  This script fixes the NVT cell at <V>.  The NVT density is therefore the
#  second number, the density at the mean volume.  It does not adjust the
#  box until the NVT density equals <rho>.
#
#  NVT versus muVT
#  ---------------
#  ``--ensemble nvt`` removes the volume move and ``# Pressure_Info``,
#  sets ``nvt_mc``, and rescales the remaining move probabilities so they
#  sum to 1.
#
#  ``--ensemble gcmc`` does the same scaling of the box, then writes a
#  template: chemical potential 0.0 and insertion/deletion probabilities
#  of 0.15 each.  Those numbers are placeholders.  Cassandra will accept
#  the file, but the chemical potential is not a physical result of the
#  NPT run.  Edit it before submitting the job.  Also raise the maximum
#  molecule count on ``# Molecule_Files`` if insertions must be allowed
#  to increase N above the NPT value.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/16/26
#
#*******************************************************************************
"""Create an NVT input, or a GCMC template, at the NPT average volume."""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.hfile import read_h
from diagnostics.io.inp import read_molecule_files_from_inp
from diagnostics.io.mcf import McfData, read_mcf
from diagnostics.io.prp import read_prp
from diagnostics.io.xyz import read_xyz
from utilities.block_average import analyze_prp

# Move sections whose first data line is a probability.  Volume is removed
# rather than rescaled.  Insertion and deletion are added for GCMC and are
# not taken from the NPT input.
_SCALABLE_PROBS = (
    "# Prob_Translation",
    "# Prob_Rotation",
    "# Prob_Regrowth",
    "# Prob_Angle",
    "# Prob_Dihedral",
    "# Prob_Atom_Displacement",
    "# Prob_Ring",
    "# Prob_Swap",
)
_GCMC_INSERT = 0.15


def _strip(line: str) -> str:
    return line.strip()


def _ignorable(line: str) -> bool:
    text = _strip(line)
    return text == "" or text.startswith("!")


def _header(line: str) -> str | None:
    """Return ``# Section_Name`` for a Cassandra section line.

    Input files write ``# Run_Name`` with a space after the hash.  Extra
    tokens, such as the ``1`` in ``# Property_Info 1``, are ignored.
    """
    parts = line.strip().split()
    if not parts or not parts[0].startswith("#"):
        return None
    if parts[0] == "#":
        if len(parts) < 2:
            return "#"
        return "# " + parts[1]
    return parts[0]


def _first_token_replaced(line: str, new_token: str) -> str:
    parts = line.split()
    if not parts:
        return new_token
    parts[0] = new_token
    return " ".join(parts)


def _retarget(token: str, inp_dir: Path, out_dir: Path) -> str:
    """Point a filename token at the original input directory.

    The new input is written in ``out_dir``, so a relative name such as
    ``r125.mcf`` has to be rewritten.  Forward slashes are used so the file
    is valid on the CRC as well as on a Mac.
    """
    source = Path(token)
    target = source if source.is_absolute() else inp_dir / source
    # Resolve both sides: macOS exposes /tmp as /private/tmp, and a relative
    # path calculated from the alias is wrong when Cassandra opens it later.
    return Path(os.path.relpath(target.resolve(), out_dir.resolve())).as_posix()


def parse_move_probabilities(lines: list[str]) -> dict[str, float]:
    """First numeric line under each scalable ``# Prob_`` section."""
    found: dict[str, float] = {}
    pending: str | None = None
    seen_data = False
    skipping_volume = False
    for line in lines:
        header = _header(line)
        if header is not None:
            skipping_volume = header == "# Prob_Volume"
            pending = header if header in _SCALABLE_PROBS else None
            seen_data = False
            continue
        if skipping_volume or pending is None or seen_data or _ignorable(line):
            continue
        seen_data = True
        found[pending] = float(line.split()[0])
    return found


def scaled_probabilities(original: dict[str, float], ensemble: str) -> dict[str, float]:
    """Renormalize NPT move probabilities after the volume move is removed."""
    if not original:
        raise ValueError("No move probabilities found to rescale")
    total = sum(original.values())
    if total <= 0.0:
        raise ValueError("Move probabilities sum to zero")
    reserved = 0.0 if ensemble == "nvt" else 2.0 * _GCMC_INSERT
    if reserved >= 1.0:
        raise ValueError("GCMC insertion and deletion leave no probability for other moves")
    factor = (1.0 - reserved) / total
    return {name: value * factor for name, value in original.items()}


def transform_inp(
    lines: list[str],
    *,
    ensemble: str,
    run_name: str,
    box_length: float,
    start_line: str,
    nspecies: int,
    inp_dir: Path,
    out_dir: Path,
) -> list[str]:
    """Return a new input-file body for ``nvt`` or ``gcmc``.

    ``ensemble`` is ``nvt`` or ``gcmc``.  File paths in ``# Molecule_Files``
    and ``# Fragment_Files`` are rewritten relative to ``out_dir``.
    """
    if ensemble not in {"nvt", "gcmc"}:
        raise ValueError(f"ensemble must be nvt or gcmc, not {ensemble!r}")

    scaled = scaled_probabilities(parse_move_probabilities(lines), ensemble)
    out: list[str] = [
        f"! Generated by utilities/npt_restart.py from an NPT average volume.",
        f"! Ensemble: {ensemble}. Cubic box length {box_length:.10f} A.",
    ]
    if ensemble == "gcmc":
        out.append("! Chemical potential and insertion/deletion probabilities are placeholders.")
        out.append("! Set mu (kJ/mol) before running. 0.0 is not a result from the NPT simulation.")

    pending: str | None = None
    seen_data = False
    skip_section = False
    box_stage = 0
    path_section = False

    for line in lines:
        header = _header(line)
        if header is not None:
            pending = None
            seen_data = False
            box_stage = 0
            path_section = header in {"# Molecule_Files", "# Fragment_Files"}
            skip_section = header in {"# Pressure_Info", "# Prob_Volume"}
            if header == "# Pressure_Info" and ensemble == "gcmc":
                out.append("# Chemical_Potential_Info")
                out.append("! Placeholder kJ/mol for each species. Replace before running.")
                for _ in range(nspecies):
                    out.append("0.0")
                continue
            if skip_section:
                continue
            if header == "# Done_Probability_Info" and ensemble == "gcmc":
                method = "cbmc" if any(_header(item) == "# Fragment_Files" for item in lines) else "none"
                out.append("# Prob_Insertion")
                out.append(f"{_GCMC_INSERT:.6f}")
                out.append(" ".join([method] * nspecies))
                out.append("")
                out.append("# Prob_Deletion")
                out.append(f"{_GCMC_INSERT:.6f}")
                out.append("")
            out.append(line.rstrip("\n"))
            if header in {"# Run_Name", "# Sim_Type", "# Start_Type", "# Box_Info"} or header in scaled:
                pending = header
            continue

        if skip_section:
            continue

        if pending is not None and not seen_data and not _ignorable(line):
            seen_data = True
            if pending == "# Run_Name":
                out.append(run_name)
                continue
            if pending == "# Sim_Type":
                out.append("nvt_mc" if ensemble == "nvt" else "gcmc")
                continue
            if pending == "# Start_Type":
                out.append(start_line)
                continue
            if pending == "# Box_Info":
                nboxes = int(line.split()[0])
                if nboxes != 1:
                    raise ValueError("npt_restart.py supports one box only")
                box_stage = 1
                out.append(line.rstrip("\n"))
                continue
            if pending in scaled:
                out.append(_first_token_replaced(line, f"{scaled[pending]:.6f}"))
                continue

        if pending == "# Box_Info" and box_stage > 0 and not _ignorable(line):
            box_stage += 1
            if box_stage == 2:
                shape = line.split()[0].lower()
                if shape != "cubic":
                    raise ValueError(f"npt_restart.py supports cubic boxes only, found {shape}")
                out.append(line.rstrip("\n"))
                continue
            if box_stage == 3:
                out.append(f"{box_length:.10f}")
                pending = None
                continue

        if path_section and not _ignorable(line):
            parts = line.split()
            parts[0] = _retarget(parts[0], inp_dir, out_dir)
            out.append(" ".join(parts))
            continue

        out.append(line.rstrip("\n"))

    return out


def _load_average(prp: Path | None, averages: Path | None, discard_sweeps: float | None) -> dict:
    if averages is not None:
        report = json.loads(averages.read_text(encoding="utf-8"))
    elif prp is not None:
        report = analyze_prp(read_prp(prp), discard_sweeps=discard_sweeps).as_dict()
    else:
        raise ValueError("Provide --prp or --averages")
    if report.get("cubic_length_A") is None or report.get("volume_average") is None:
        raise ValueError("Average report has no Volume column, so no cubic length is available")
    return report


def _mic_delta(delta: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Minimum-image displacement in a cell stored as a 3-by-3 H matrix."""
    fractional = np.linalg.solve(cell, delta)
    fractional -= np.rint(fractional)
    return cell @ fractional


def _unwrapped_molecule(
    coords: np.ndarray,
    mcf: McfData,
    cell: np.ndarray,
) -> np.ndarray:
    """Make one molecular image contiguous by traversing its MCF bond graph."""
    n_atoms = mcf.n_atoms
    if coords.shape != (n_atoms, 3):
        raise ValueError(
            f"{mcf.path}: expected {n_atoms} coordinates per molecule, got {coords.shape}"
        )
    if n_atoms == 1:
        return coords.copy()

    neighbors: list[list[int]] = [[] for _ in range(n_atoms)]
    for bond in mcf.bonds:
        neighbors[bond.i].append(bond.j)
        neighbors[bond.j].append(bond.i)
    if not mcf.bonds:
        raise ValueError(
            f"{mcf.path}: a multi-atom molecule needs bonds to reconstruct it across PBC"
        )

    unwrapped = np.empty_like(coords)
    visited = np.zeros(n_atoms, dtype=bool)
    visited[0] = True
    unwrapped[0] = coords[0]
    queue = [0]
    while queue:
        parent = queue.pop(0)
        for child in neighbors[parent]:
            candidate = unwrapped[parent] + _mic_delta(coords[child] - coords[parent], cell)
            if not visited[child]:
                unwrapped[child] = candidate
                visited[child] = True
                queue.append(child)
            elif not np.allclose(
                unwrapped[child],
                candidate,
                rtol=0.0,
                atol=1.0e-6,
            ):
                # A closed path may reach the same atom through a different
                # periodic image. It is consistent when the difference is a
                # lattice vector, i.e. its minimum-image value is zero.
                if np.allclose(
                    _mic_delta(candidate - unwrapped[child], cell),
                    0.0,
                    rtol=0.0,
                    atol=1.0e-6,
                ):
                    continue
                raise ValueError(
                    f"{mcf.path}: periodic bond graph is inconsistent for a ring; "
                    "cannot reconstruct one contiguous molecular image"
                )
    if not np.all(visited):
        missing = ", ".join(str(i + 1) for i in np.where(~visited)[0])
        raise ValueError(
            f"{mcf.path}: disconnected multi-atom topology; atoms {missing} "
            "cannot share one molecule-wide COM"
        )
    return unwrapped


def _move_molecule_com(
    coords: np.ndarray,
    mcf: McfData,
    old_cell: np.ndarray,
    new_cell: np.ndarray,
) -> np.ndarray:
    """Map a molecule COM between cells without altering its internal vectors."""
    contiguous = _unwrapped_molecule(coords, mcf, old_cell)
    masses = np.asarray([atom.mass for atom in mcf.atoms], dtype=float)
    total_mass = float(masses.sum())
    if total_mass <= 0.0:
        raise ValueError(f"{mcf.path}: molecule has nonpositive total mass")

    com_old = np.average(contiguous, axis=0, weights=masses)
    fractional_com = np.linalg.solve(old_cell, com_old)
    com_new = new_cell @ fractional_com
    moved = contiguous + (com_new - com_old)

    # Keep one molecular image together. Read_Config will apply the same
    # whole-molecule PBC shift when it loads this configuration.
    new_fractional_com = np.linalg.solve(new_cell, com_new)
    lattice_shift = new_cell @ (-np.floor(new_fractional_com))
    return moved + lattice_shift


def _move_all_molecule_coms(
    coords: np.ndarray,
    *,
    inp_path: Path,
    nmols: list[int],
    old_cell: np.ndarray,
    new_cell: np.ndarray,
) -> np.ndarray:
    """Move every molecular COM while preserving coordinates within molecules."""
    entries = read_molecule_files_from_inp(inp_path)
    if len(entries) != len(nmols):
        raise ValueError(
            f"{inp_path}: # Molecule_Files has {len(entries)} species but {len(nmols)} "
            f"counts appear in the restart .H"
        )

    moved = np.empty_like(coords)
    offset = 0
    for entry, count in zip(entries, nmols):
        mcf = read_mcf(entry.mcf_path)
        if count < 0:
            raise ValueError(f"Restart .H has a negative molecule count for species {entry.species_index}")
        for _ in range(count):
            end = offset + mcf.n_atoms
            if end > coords.shape[0]:
                raise ValueError(
                    f"Restart XYZ ends inside species {entry.species_index}; "
                    "atom count disagrees with the .H and MCF files"
                )
            moved[offset:end] = _move_molecule_com(
                coords[offset:end], mcf, old_cell, new_cell
            )
            offset = end
    if offset != coords.shape[0]:
        raise ValueError(
            f"Restart XYZ has {coords.shape[0]} atoms but the .H and MCF files "
            f"describe {offset}"
        )
    return moved


def scale_configuration(
    xyz_path: Path,
    h_path: Path,
    length_avg: float,
    inp_path: Path,
) -> tuple[list[str], np.ndarray, np.ndarray, float, list[int]]:
    """Move molecular COMs into a cubic cell of side ``length_avg``.

    Returns elements are ``(xyz_lines, cell, coords, last_length, nmols)``.
    ``coords`` preserve every molecule's internal Cartesian geometry.
    """
    xyz = read_xyz(xyz_path)
    if xyz.n_frames != 1:
        raise ValueError(
            f"{xyz_path} has {xyz.n_frames} frames. Pass the single-frame "
            "*.restart.xyz file, not the movie."
        )
    htraj = read_h(h_path)
    if htraj.n_frames != 1:
        raise ValueError(
            f"{h_path} has {htraj.n_frames} frames. Pass the single-frame *.restart.H file."
        )
    frame = htraj.frames[0]
    cell = np.array(frame.cell, dtype=float)
    volume_cell = abs(float(np.linalg.det(cell)))
    if volume_cell <= 0.0:
        raise ValueError(f"Cell in {h_path} has non-positive volume")
    if abs(frame.volume - volume_cell) / volume_cell > 1.0e-4:
        raise ValueError(
            f"Volume in {h_path} ({frame.volume:.8g}) disagrees with det(H) ({volume_cell:.8g})"
        )
    length_last = volume_cell ** (1.0 / 3.0)
    off = cell.copy()
    np.fill_diagonal(off, 0.0)
    if np.max(np.abs(off)) > 1.0e-4 * length_last:
        raise ValueError("Cell is not cubic (off-diagonal H elements are nonzero)")
    diagonal = np.diag(cell)
    if np.max(np.abs(diagonal - diagonal[0])) > 1.0e-4 * length_last:
        raise ValueError("Cell is not cubic (the three edge lengths differ)")

    new_cell = cell * (length_avg / length_last)
    new_coords = _move_all_molecule_coms(
        xyz.coords[0],
        inp_path=inp_path,
        nmols=list(frame.nmols),
        old_cell=cell,
        new_cell=new_cell,
    )
    lines = [
        f"{xyz.n_atoms:8d}",
        f" RESTART COM-mapped L={length_avg:.8f} from {xyz_path.name}",
    ]
    for element, xyz_row in zip(xyz.elements, new_coords):
        lines.append(
            f"{element:<2s} {xyz_row[0]:18.8f} {xyz_row[1]:18.8f} {xyz_row[2]:18.8f}"
        )
    return lines, new_cell, new_coords, length_last, list(frame.nmols)


def write_h(path: Path, cell: np.ndarray, nmols: list[int]) -> None:
    """Write one Cassandra ``.H`` block for a cubic (or general) cell."""
    volume = abs(float(np.linalg.det(cell)))
    rows = [
        f"   {volume:.14f}",
        "   " + "  ".join(f"{value:.14f}" for value in cell[0]),
        "   " + "  ".join(f"{value:.14f}" for value in cell[1]),
        "   " + "  ".join(f"{value:.14f}" for value in cell[2]),
        "",
        f"   {len(nmols):d}",
    ]
    for index, count in enumerate(nmols, start=1):
        rows.append(f"   {index:d}   {count:d}")
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def default_restart_paths(prp: Path) -> tuple[Path, Path]:
    xyz = prp.with_name(prp.stem + ".restart.xyz")
    hfile = prp.with_name(prp.stem + ".restart.H")
    return xyz, hfile


def original_run_name(lines: list[str]) -> str:
    pending = False
    for line in lines:
        if _header(line) == "# Run_Name":
            pending = True
            continue
        if pending and not _ignorable(line):
            return line.split()[0]
    raise ValueError("Input file has no # Run_Name value")


def build_restart(
    *,
    inp: Path,
    out_dir: Path,
    ensemble: str,
    prp: Path | None,
    averages: Path | None,
    restart_xyz: Path | None,
    restart_h: Path | None,
    discard_sweeps: float | None,
    run_name: str | None,
) -> Path:
    """Write the scaled configuration and the new input. Return the ``.inp`` path."""
    report = _load_average(prp, averages, discard_sweeps)
    if restart_xyz is None or restart_h is None:
        if prp is None:
            raise ValueError("Pass --restart-xyz and --restart-h, or pass --prp so they can be inferred")
        guessed_xyz, guessed_h = default_restart_paths(prp)
        restart_xyz = restart_xyz or guessed_xyz
        restart_h = restart_h or guessed_h

    length = float(report["cubic_length_A"])
    xyz_lines, cell, _coords, length_last, nmols = scale_configuration(
        restart_xyz, restart_h, length, inp
    )
    inp_lines = inp.read_text(encoding="utf-8").splitlines()
    name = run_name or (original_run_name(inp_lines) + ("-nvt" if ensemble == "nvt" else "-gcmc"))
    start_line = "read_config " + " ".join(str(count) for count in nmols) + f" {name}.xyz"

    out_dir.mkdir(parents=True, exist_ok=True)
    written_inp = transform_inp(
        inp_lines,
        ensemble=ensemble,
        run_name=name,
        box_length=length,
        start_line=start_line,
        nspecies=len(nmols),
        inp_dir=inp.parent,
        out_dir=out_dir,
    )
    inp_path = out_dir / f"{name}.inp"
    xyz_path = out_dir / f"{name}.xyz"
    h_path = out_dir / f"{name}.H"
    json_path = out_dir / f"{name}.averages.json"
    inp_path.write_text("\n".join(written_inp) + "\n", encoding="utf-8")
    xyz_path.write_text("\n".join(xyz_lines) + "\n", encoding="utf-8")
    write_h(h_path, cell, nmols)
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(f"Last-frame length: {length_last:.8f} A")
    print(f"Average length:    {length:.8f} A")
    print(f"Cell scale factor: {length / length_last:.8f}")
    print(f"Wrote {inp_path}")
    print(f"Wrote {xyz_path}")
    print(f"Wrote {h_path}")
    if ensemble == "gcmc":
        print("GCMC template: edit # Chemical_Potential_Info before running.")
        print("GCMC template: insertion and deletion are each 0.15. Change them if you want a different mix.")
        print("GCMC template: raise the molecule maximum on # Molecule_Files if N must be allowed to grow.")
    return inp_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Map last-frame molecular COMs into the NPT average-volume cell "
            "and write an NVT input or a GCMC template."
        )
    )
    parser.add_argument("--inp", type=Path, required=True, help="NPT Cassandra input file")
    parser.add_argument("--out-dir", type=Path, required=True, help="New directory for the restart files")
    parser.add_argument("--ensemble", choices=("nvt", "gcmc"), default="nvt")
    parser.add_argument("--prp", type=Path, help="NPT .prp file. Restart files are inferred beside it.")
    parser.add_argument("--averages", type=Path, help="JSON written by block_average.py")
    parser.add_argument("--restart-xyz", type=Path, help="Single-frame *.restart.xyz")
    parser.add_argument("--restart-h", type=Path, help="Single-frame *.restart.H")
    parser.add_argument("--discard-sweeps", type=float, help="Passed through when --prp is analyzed here")
    parser.add_argument("--run-name", help="Run_Name for the new input (default: old name plus -nvt or -gcmc)")
    args = parser.parse_args(argv)
    build_restart(
        inp=args.inp,
        out_dir=args.out_dir,
        ensemble=args.ensemble,
        prp=args.prp,
        averages=args.averages,
        restart_xyz=args.restart_xyz,
        restart_h=args.restart_h,
        discard_sweeps=args.discard_sweeps,
        run_name=args.run_name,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
