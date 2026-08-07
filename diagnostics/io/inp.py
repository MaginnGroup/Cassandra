#*******************************************************************************
#  diagnostics/io/inp.py
#
#  Shared helpers for Cassandra ``.inp`` sections used by diagnostics.
#
#  Currently:
#    * ``# Temperature_Info`` → T (K)
#    * ``# Molecule_Files``  → ordered list of (mcf path, nmol) per species
#    * ``# Sim_Type``        → ensemble keyword (nvt_mc, npt_mc, …)
#    * ``# Box_Info``        → cubic edge length for box 1 (when cubic)
#
#  Multi-species XYZ layout (see docs/output-formats.md)
#  -----------------------------------------------------
#  Cassandra writes atoms species-outermost, then molecule, then MCF atom
#  order.  ``# Molecule_Files`` gives that species order and the intended
#  molecule counts for fixed-N runs (NVT/NPT).  Variable-N (GCMC) needs the
#  companion ``.H`` file — not handled here yet; pass ``--nmols`` if the
#  counts in the ``.inp`` do not match the movie.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#*******************************************************************************
"""Parse temperature and Molecule_Files from Cassandra ``.inp`` files."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass
class MoleculeFileEntry:
    """One ``# Molecule_Files`` line: species MCF path and molecule count."""

    mcf_path: Path
    nmol: int
    species_index: int  # 1-based species ID (order in the .inp)


def read_temperature_from_inp(path: str | Path) -> float:
    """Read the first box temperature (K) from ``# Temperature_Info``."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for i, raw in enumerate(lines):
        s = raw.strip()
        if s.startswith("#") and s.lstrip("#").strip().lower() == "temperature_info":
            for raw2 in lines[i + 1 :]:
                t = raw2.strip()
                if not t or t.startswith("!") or t.startswith("#"):
                    continue
                try:
                    return float(t.split()[0])
                except (ValueError, IndexError) as exc:
                    raise ValueError(
                        f"{path}: could not parse temperature after "
                        f"# Temperature_Info: {raw2!r}"
                    ) from exc
            raise ValueError(f"{path}: # Temperature_Info has no temperature value")
    raise ValueError(f"{path}: # Temperature_Info section not found")


def read_molecule_files_from_inp(path: str | Path) -> list[MoleculeFileEntry]:
    """Parse ``# Molecule_Files`` (one ``mcf  nmol`` line per species).

    MCF paths are resolved relative to the directory containing the ``.inp``.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    base = path.parent
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for i, raw in enumerate(lines):
        s = raw.strip()
        if s.startswith("#") and s.lstrip("#").strip().lower() == "molecule_files":
            start = i + 1
            break
    if start is None:
        raise ValueError(f"{path}: # Molecule_Files section not found")

    entries: list[MoleculeFileEntry] = []
    for raw in lines[start:]:
        s = raw.strip()
        if not s or s.startswith("!"):
            continue
        if s.startswith("#"):
            break
        tok = s.split()
        if len(tok) < 2:
            raise ValueError(
                f"{path}: bad # Molecule_Files line (need 'mcf nmol'): {raw!r}"
            )
        mcf_name, nmol_s = tok[0], tok[1]
        try:
            nmol = int(nmol_s)
        except ValueError as exc:
            raise ValueError(
                f"{path}: bad molecule count in # Molecule_Files: {raw!r}"
            ) from exc
        if nmol < 0:
            raise ValueError(f"{path}: negative nmol in # Molecule_Files: {raw!r}")
        mcf_path = Path(mcf_name)
        if not mcf_path.is_absolute():
            mcf_path = (base / mcf_path).resolve()
        entries.append(
            MoleculeFileEntry(
                mcf_path=mcf_path,
                nmol=nmol,
                species_index=len(entries) + 1,
            )
        )

    if not entries:
        raise ValueError(f"{path}: # Molecule_Files has no species lines")
    return entries


def read_sim_type_from_inp(path: str | Path) -> str:
    """Return the ``# Sim_Type`` keyword (lowercased), e.g. ``npt_mc``."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    for i, raw in enumerate(lines):
        s = raw.strip()
        if s.startswith("#") and s.lstrip("#").strip().lower() == "sim_type":
            for raw2 in lines[i + 1 :]:
                t = raw2.strip()
                if not t or t.startswith("!") or t.startswith("#"):
                    continue
                return t.split()[0].lower()
            raise ValueError(f"{path}: # Sim_Type has no value")
    raise ValueError(f"{path}: # Sim_Type section not found")


def read_cubic_box_length_from_inp(path: str | Path) -> float:
    """Read the cubic box edge (Å) for the first box from ``# Box_Info``.

    Expected layout::

        # Box_Info
        1
        cubic
        54.2
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    start = None
    for i, raw in enumerate(lines):
        s = raw.strip()
        if s.startswith("#") and s.lstrip("#").strip().lower() == "box_info":
            start = i + 1
            break
    if start is None:
        raise ValueError(f"{path}: # Box_Info section not found")

    vals: list[str] = []
    for raw in lines[start:]:
        t = raw.strip()
        if not t or t.startswith("!"):
            continue
        if t.startswith("#"):
            break
        vals.extend(t.split())
        if len(vals) >= 3:
            break

    if len(vals) < 3:
        raise ValueError(f"{path}: # Box_Info incomplete (need nboxes, type, size)")
    box_type = vals[1].lower()
    if box_type != "cubic":
        raise ValueError(
            f"{path}: Box_Info type is '{vals[1]}' — this helper only reads cubic "
            "edges; for non-cubic / NPT movies use the companion .H file."
        )
    try:
        return float(vals[2])
    except ValueError as exc:
        raise ValueError(f"{path}: bad cubic box length in # Box_Info: {vals[2]!r}") from exc
