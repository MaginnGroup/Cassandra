#*******************************************************************************
#  diagnostics/io/mcf.py
#
#  Minimal reader for Cassandra molecular connectivity (``.mcf``) files.
#
#  Sections read for diagnostics
#  -----------------------------
#    * ``# Atom_Info``
#    * ``# Bond_Info``
#    * ``# Angle_Info``
#    * ``# Dihedral_Info``
#
#  Bond format
#  -----------
#  ::
#
#      # Bond_Info
#      <n_bonds>
#      <index>  <i>  <j>  <type>  <param...>
#
#  * ``i``, ``j`` are **1-based** MCF atom indices.
#  * ``fixed``: bond length in **Å** (Cassandra constraint; expected sharp peak).
#
#  Angle format (see docs/mcf-setup-workflow.md)
#  --------------------------------------------
#  ::
#
#      # Angle_Info
#      <n_angles>
#      <index>  <i>  <j>  <k>  <type>  <param1>  <param2>
#
#  * ``i``, ``j``, ``k`` are **1-based**; ``j`` is the vertex.
#  * ``harmonic``: K_θ in **K/rad²**, θ₀ in **degrees** (no ½ in energy).
#  * ``fixed``: equilibrium angle in degrees.
#
#  Dihedral format
#  ---------------
#  ::
#
#      # Dihedral_Info
#      <n_dihedrals>
#      <index>  <i>  <j>  <k>  <l>  <type>  <parameters...>
#
#  Supported types for Boltzmann overlays:
#    * ``OPLS``   — a0 a1 a2 a3 in **kJ/mol**
#    * ``CHARMM`` — a0 (kJ/mol), n, delta (degrees)
#    * ``RB``     — c0..c5 in **kJ/mol** (Ryckaert–Bellemans)
#    * ``harmonic`` — K in **K/rad²**, φ₀ in degrees
#    * ``none``   — listed only; no Boltzmann curve
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#*******************************************************************************
"""Parse Atom_Info, Bond_Info, Angle_Info, and Dihedral_Info from a Cassandra ``.mcf``."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

# Re-export so older imports ``from diagnostics.io.mcf import read_temperature_from_inp``
# keep working after the helper moved to ``inp.py``.
from diagnostics.io.inp import read_temperature_from_inp  # noqa: F401


@dataclass
class McfAtom:
    """One atom row from ``# Atom_Info`` (subset of fields used by diagnostics)."""

    index: int  # 1-based MCF index
    atom_type: str
    element: str
    mass: float
    charge: float


@dataclass
class McfBond:
    """One bond row from ``# Bond_Info``.

    Attributes
    ----------
    index :
        1-based bond ID (user-facing ``--bond``).
    i, j :
        0-based atom indices.
    i1, j1 :
        Original 1-based MCF atom indices.
    potential_type :
        ``fixed`` (usual Cassandra constraint) or other lowercase type.
    length_ang :
        Constraint / equilibrium length in Å.
    """

    index: int
    i: int
    j: int
    i1: int
    j1: int
    potential_type: str
    length_ang: float

    @property
    def is_fixed(self) -> bool:
        """True if this is a fixed-length constraint."""
        return self.potential_type == "fixed"

    def label(self, atoms: list[McfAtom] | None = None) -> str:
        """Short human label, e.g. ``1-2 (C-C)``."""
        if atoms is None:
            return f"{self.i1}-{self.j1}"
        try:
            return (
                f"{self.i1}-{self.j1} "
                f"({atoms[self.i].element}-{atoms[self.j].element})"
            )
        except IndexError:
            return f"{self.i1}-{self.j1}"


@dataclass
class McfAngle:
    """One angle row from ``# Angle_Info``."""

    index: int
    i: int
    j: int
    k: int
    i1: int
    j1: int
    k1: int
    potential_type: str
    k_theta: float | None
    theta0_deg: float

    @property
    def is_harmonic(self) -> bool:
        """True if this angle has a harmonic Boltzmann PDF."""
        return self.potential_type == "harmonic"

    def label(self, atoms: list[McfAtom] | None = None) -> str:
        """Short human label, e.g. ``1-2-3 (C-O-C)``."""
        if atoms is None:
            return f"{self.i1}-{self.j1}-{self.k1}"
        try:
            elems = (
                atoms[self.i].element,
                atoms[self.j].element,
                atoms[self.k].element,
            )
            return f"{self.i1}-{self.j1}-{self.k1} ({elems[0]}-{elems[1]}-{elems[2]})"
        except IndexError:
            return f"{self.i1}-{self.j1}-{self.k1}"


@dataclass
class McfDihedral:
    """One dihedral row from ``# Dihedral_Info``.

    Attributes
    ----------
    index :
        1-based dihedral ID (user-facing ``--dihedral``).
    i, j, k, l :
        0-based atom indices (central bond is ``j``–``k``).
    potential_type :
        Lowercased type string (``opls``, ``charmm``, ``rb``, ``harmonic``, …).
    params :
        Raw MCF parameters (units depend on type; see module header).
    """

    index: int
    i: int
    j: int
    k: int
    l: int
    i1: int
    j1: int
    k1: int
    l1: int
    potential_type: str
    params: list[float] = field(default_factory=list)

    @property
    def has_boltzmann(self) -> bool:
        """True if an ideal-gas Boltzmann PDF can be built from MCF params."""
        return self.potential_type in {"opls", "charmm", "rb", "harmonic", "ryckaert-bellemans"}

    def label(self, atoms: list[McfAtom] | None = None) -> str:
        """Short human label, e.g. ``1-2-3-4 (C-C-C-C)``."""
        ids = f"{self.i1}-{self.j1}-{self.k1}-{self.l1}"
        if atoms is None:
            return ids
        try:
            elems = (
                atoms[self.i].element,
                atoms[self.j].element,
                atoms[self.k].element,
                atoms[self.l].element,
            )
            return f"{ids} ({elems[0]}-{elems[1]}-{elems[2]}-{elems[3]})"
        except IndexError:
            return ids

    def param_summary(self) -> str:
        """One-line parameter summary for ``--list`` output."""
        p = self.params
        t = self.potential_type
        if t == "opls" and len(p) >= 4:
            return (
                f"OPLS  a0={p[0]:.4g} a1={p[1]:.4g} a2={p[2]:.4g} "
                f"a3={p[3]:.4g} kJ/mol"
            )
        if t == "charmm" and len(p) >= 3:
            return f"CHARMM  a0={p[0]:.4g} kJ/mol  n={p[1]:.4g}  δ={p[2]:.4g}°"
        if t in {"rb", "ryckaert-bellemans"} and p:
            coeffs = " ".join(f"{c:.4g}" for c in p[:6])
            return f"RB  c={coeffs} kJ/mol"
        if t == "harmonic" and len(p) >= 2:
            return f"harmonic  K={p[0]:.4g} K/rad²  φ₀={p[1]:.4g}°"
        if t in {"none", "fixed"}:
            return t
        return f"{t}  params={p}"


@dataclass
class McfData:
    """Subset of an MCF needed for intramolecular diagnostics."""

    path: Path
    atoms: list[McfAtom]
    bonds: list[McfBond]
    angles: list[McfAngle]
    dihedrals: list[McfDihedral]

    @property
    def n_atoms(self) -> int:
        return len(self.atoms)

    @property
    def n_bonds(self) -> int:
        return len(self.bonds)

    @property
    def n_angles(self) -> int:
        return len(self.angles)

    @property
    def n_dihedrals(self) -> int:
        return len(self.dihedrals)

    def bond_by_index(self, index: int) -> McfBond:
        for bond in self.bonds:
            if bond.index == index:
                return bond
        available = ", ".join(str(b.index) for b in self.bonds) or "(none)"
        raise ValueError(
            f"No bond with index {index} in {self.path}. Available: {available}"
        )

    def angle_by_index(self, index: int) -> McfAngle:
        for ang in self.angles:
            if ang.index == index:
                return ang
        available = ", ".join(str(a.index) for a in self.angles) or "(none)"
        raise ValueError(
            f"No angle with index {index} in {self.path}. Available: {available}"
        )

    def dihedral_by_index(self, index: int) -> McfDihedral:
        for dih in self.dihedrals:
            if dih.index == index:
                return dih
        available = ", ".join(str(d.index) for d in self.dihedrals) or "(none)"
        raise ValueError(
            f"No dihedral with index {index} in {self.path}. Available: {available}"
        )


def _section_lines(lines: list[str], header: str) -> list[str]:
    """Return non-comment body lines after ``# header`` until the next ``#`` / END."""
    header_l = header.lower()
    start = None
    for i, raw in enumerate(lines):
        s = raw.strip()
        if s.startswith("#") and s.lstrip("#").strip().lower() == header_l:
            start = i + 1
            break
    if start is None:
        raise ValueError(f"MCF section '# {header}' not found")

    body: list[str] = []
    for raw in lines[start:]:
        s = raw.strip()
        if not s or s.startswith("!"):
            continue
        if s.upper() == "END":
            break
        if s.startswith("#"):
            break
        body.append(s)
    return body


def _parse_atoms(path: Path, lines: list[str]) -> list[McfAtom]:
    atom_body = _section_lines(lines, "Atom_Info")
    if not atom_body:
        raise ValueError(f"{path}: empty # Atom_Info section")
    try:
        n_atoms = int(atom_body[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: bad atom count in # Atom_Info: {atom_body[0]!r}") from exc
    if len(atom_body) < 1 + n_atoms:
        raise ValueError(f"{path}: expected {n_atoms} atom rows, found {len(atom_body) - 1}")

    atoms: list[McfAtom] = []
    for row in atom_body[1 : 1 + n_atoms]:
        tok = row.split()
        if len(tok) < 5:
            raise ValueError(f"{path}: bad Atom_Info row: {row!r}")
        atoms.append(
            McfAtom(
                index=int(tok[0]),
                atom_type=tok[1],
                element=tok[2],
                mass=float(tok[3]),
                charge=float(tok[4]),
            )
        )
    return atoms


def _parse_bonds(path: Path, lines: list[str], n_atoms: int) -> list[McfBond]:
    bond_body = _section_lines(lines, "Bond_Info")
    if not bond_body:
        raise ValueError(f"{path}: empty # Bond_Info section")
    try:
        n_bonds = int(bond_body[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: bad bond count in # Bond_Info: {bond_body[0]!r}") from exc

    bonds: list[McfBond] = []
    if n_bonds == 0:
        return bonds
    if len(bond_body) < 1 + n_bonds:
        raise ValueError(f"{path}: expected {n_bonds} bond rows, found {len(bond_body) - 1}")

    for row in bond_body[1 : 1 + n_bonds]:
        tok = row.split()
        if len(tok) < 5:
            raise ValueError(f"{path}: bad Bond_Info row: {row!r}")
        idx = int(tok[0])
        i1, j1 = int(tok[1]), int(tok[2])
        pot = tok[3].lower()
        if pot == "fixed":
            if len(tok) < 5:
                raise ValueError(f"{path}: fixed bond {idx} needs length: {row!r}")
            length = float(tok[4])
        else:
            # Rare / future types: take first numeric parameter as length if present
            if len(tok) < 5:
                raise ValueError(f"{path}: bond {idx} type '{tok[3]}' needs a length: {row!r}")
            length = float(tok[4])

        for a1, name in ((i1, "i"), (j1, "j")):
            if a1 < 1 or a1 > n_atoms:
                raise ValueError(
                    f"{path}: bond {idx} atom {name}={a1} out of range [1, {n_atoms}]"
                )

        bonds.append(
            McfBond(
                index=idx,
                i=i1 - 1,
                j=j1 - 1,
                i1=i1,
                j1=j1,
                potential_type=pot,
                length_ang=length,
            )
        )
    return bonds


def _parse_angles(path: Path, lines: list[str], n_atoms: int) -> list[McfAngle]:
    angle_body = _section_lines(lines, "Angle_Info")
    if not angle_body:
        raise ValueError(f"{path}: empty # Angle_Info section")
    try:
        n_angles = int(angle_body[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(f"{path}: bad angle count in # Angle_Info: {angle_body[0]!r}") from exc

    angles: list[McfAngle] = []
    if n_angles == 0:
        return angles
    if len(angle_body) < 1 + n_angles:
        raise ValueError(f"{path}: expected {n_angles} angle rows, found {len(angle_body) - 1}")

    for row in angle_body[1 : 1 + n_angles]:
        tok = row.split()
        if len(tok) < 5:
            raise ValueError(f"{path}: bad Angle_Info row: {row!r}")
        idx = int(tok[0])
        i1, j1, k1 = int(tok[1]), int(tok[2]), int(tok[3])
        pot = tok[4].lower()
        if pot == "harmonic":
            if len(tok) < 7:
                raise ValueError(f"{path}: harmonic angle {idx} needs K and theta0: {row!r}")
            k_theta = float(tok[5])
            theta0 = float(tok[6])
        elif pot == "fixed":
            if len(tok) < 6:
                raise ValueError(
                    f"{path}: fixed angle {idx} needs equilibrium angle: {row!r}"
                )
            k_theta = None
            theta0 = float(tok[5])
        else:
            raise ValueError(f"{path}: unsupported angle type '{tok[4]}' on angle {idx}")

        for a1, name in ((i1, "i"), (j1, "j"), (k1, "k")):
            if a1 < 1 or a1 > n_atoms:
                raise ValueError(
                    f"{path}: angle {idx} atom {name}={a1} out of range [1, {n_atoms}]"
                )

        angles.append(
            McfAngle(
                index=idx,
                i=i1 - 1,
                j=j1 - 1,
                k=k1 - 1,
                i1=i1,
                j1=j1,
                k1=k1,
                potential_type=pot,
                k_theta=k_theta,
                theta0_deg=theta0,
            )
        )
    return angles


def _parse_dihedrals(path: Path, lines: list[str], n_atoms: int) -> list[McfDihedral]:
    try:
        dih_body = _section_lines(lines, "Dihedral_Info")
    except ValueError:
        # Some tiny MCFs may omit the section; treat as zero dihedrals
        return []

    if not dih_body:
        return []
    try:
        n_dih = int(dih_body[0].split()[0])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            f"{path}: bad dihedral count in # Dihedral_Info: {dih_body[0]!r}"
        ) from exc

    dihedrals: list[McfDihedral] = []
    if n_dih == 0:
        return dihedrals
    if len(dih_body) < 1 + n_dih:
        raise ValueError(
            f"{path}: expected {n_dih} dihedral rows, found {len(dih_body) - 1}"
        )

    for row in dih_body[1 : 1 + n_dih]:
        tok = row.split()
        if len(tok) < 6:
            raise ValueError(f"{path}: bad Dihedral_Info row: {row!r}")
        idx = int(tok[0])
        i1, j1, k1, l1 = int(tok[1]), int(tok[2]), int(tok[3]), int(tok[4])
        pot_raw = tok[5]
        pot = pot_raw.lower()
        if pot == "ryckaert-bellemans":
            pot = "rb"

        params = [float(x) for x in tok[6:]]

        if pot == "opls" and len(params) < 4:
            raise ValueError(f"{path}: OPLS dihedral {idx} needs a0..a3: {row!r}")
        if pot == "charmm" and len(params) < 3:
            raise ValueError(f"{path}: CHARMM dihedral {idx} needs a0, n, delta: {row!r}")
        if pot == "harmonic" and len(params) < 2:
            raise ValueError(
                f"{path}: harmonic dihedral {idx} needs K and phi0: {row!r}"
            )
        if pot == "rb" and len(params) < 1:
            raise ValueError(f"{path}: RB dihedral {idx} needs at least c0: {row!r}")

        for a1, name in ((i1, "i"), (j1, "j"), (k1, "k"), (l1, "l")):
            if a1 < 1 or a1 > n_atoms:
                raise ValueError(
                    f"{path}: dihedral {idx} atom {name}={a1} out of range "
                    f"[1, {n_atoms}]"
                )

        dihedrals.append(
            McfDihedral(
                index=idx,
                i=i1 - 1,
                j=j1 - 1,
                k=k1 - 1,
                l=l1 - 1,
                i1=i1,
                j1=j1,
                k1=k1,
                l1=l1,
                potential_type=pot,
                params=params,
            )
        )
    return dihedrals


def read_mcf(path: str | Path) -> McfData:
    """Read Atom_Info, Bond_Info, Angle_Info, and Dihedral_Info from a ``.mcf``."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"MCF file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    atoms = _parse_atoms(path, lines)
    bonds = _parse_bonds(path, lines, len(atoms))
    angles = _parse_angles(path, lines, len(atoms))
    dihedrals = _parse_dihedrals(path, lines, len(atoms))
    return McfData(
        path=path,
        atoms=atoms,
        bonds=bonds,
        angles=angles,
        dihedrals=dihedrals,
    )
