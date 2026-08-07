#*******************************************************************************
#  diagnostics/io/layout.py
#
#  Map multi-species Cassandra XYZ frames onto per-molecule coordinate blocks.
#
#  Atom order in ``.xyz`` (Write_Coords_XYZ)
#  ----------------------------------------
#      for species = 1 .. nspecies:
#          for molecule = 1 .. nmol(species):
#              for atom = 1 .. natoms(species):   # MCF order
#                  write element x y z
#
#  Fixed-N diagnostics need ``nmol`` per species (from ``# Molecule_Files``
#  and/or ``--nmols``).  Total atoms per frame must equal
#  Σ_s natoms(s) * nmol(s).
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#*******************************************************************************
"""Species / molecule indexing helpers for Cassandra XYZ trajectories."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from diagnostics.io.inp import read_molecule_files_from_inp
from diagnostics.io.mcf import McfData, read_mcf


@dataclass
class SpeciesBlock:
    """One species in the XYZ atom stream."""

    species_index: int  # 1-based
    mcf: McfData
    nmol: int
    atom_offset: int  # first atom index of this species in a frame (0-based)

    @property
    def n_atoms(self) -> int:
        return self.mcf.n_atoms

    @property
    def n_frame_atoms(self) -> int:
        """Atoms contributed by this species in one frame."""
        return self.n_atoms * self.nmol


@dataclass
class SpeciesLayout:
    """Full fixed-N layout for one XYZ movie."""

    species: list[SpeciesBlock]

    @property
    def n_species(self) -> int:
        return len(self.species)

    @property
    def n_atoms_frame(self) -> int:
        return sum(s.n_frame_atoms for s in self.species)

    def block(self, species_index: int) -> SpeciesBlock:
        """Return the block for 1-based ``species_index``."""
        for sp in self.species:
            if sp.species_index == species_index:
                return sp
        avail = ", ".join(str(s.species_index) for s in self.species)
        raise ValueError(f"Species {species_index} not in layout. Available: {avail}")

    def molecule_coords(
        self,
        frame: np.ndarray,
        species_index: int,
        molecule: int,
    ) -> np.ndarray:
        """Return coords ``(n_atoms, 3)`` for one molecule (1-based IDs)."""
        sp = self.block(species_index)
        if molecule < 1 or molecule > sp.nmol:
            raise ValueError(
                f"Molecule {molecule} out of range [1, {sp.nmol}] "
                f"for species {species_index}"
            )
        start = sp.atom_offset + (molecule - 1) * sp.n_atoms
        return frame[start : start + sp.n_atoms]


def single_species_layout(mcf: McfData, nmol: int) -> SpeciesLayout:
    """Layout when the XYZ contains only one species."""
    if nmol < 1:
        raise ValueError(f"nmol must be >= 1, got {nmol}")
    return SpeciesLayout(
        species=[
            SpeciesBlock(species_index=1, mcf=mcf, nmol=nmol, atom_offset=0)
        ]
    )


def build_layout_from_inp(
    inp_path: Path,
    *,
    topology_mcf: McfData,
    species_index: int = 1,
    nmols_override: list[int] | None = None,
) -> SpeciesLayout:
    """Build layout from ``# Molecule_Files`` (+ optional ``--nmols`` override).

    ``topology_mcf`` is the MCF used for angle/dihedral analysis of
    ``species_index``.  Other species are loaded from Molecule_Files only to
    compute atom offsets in the XYZ stream.
    """
    entries = read_molecule_files_from_inp(inp_path)
    if species_index < 1 or species_index > len(entries):
        raise ValueError(
            f"--species {species_index} out of range [1, {len(entries)}]"
        )

    nmols = [e.nmol for e in entries]
    if nmols_override is not None:
        if len(nmols_override) != len(entries):
            raise ValueError(
                f"--nmols has {len(nmols_override)} values but "
                f"# Molecule_Files has {len(entries)} species"
            )
        nmols = list(nmols_override)

    # Selected species: user topology MCF must match atom count of Molecule_Files MCF
    sel = entries[species_index - 1]
    if sel.mcf_path.is_file():
        sel_n = read_mcf(sel.mcf_path).n_atoms
        if sel_n != topology_mcf.n_atoms:
            raise ValueError(
                f"Species {species_index} Molecule_Files MCF ({sel.mcf_path.name}) "
                f"has {sel_n} atoms but topology MCF ({topology_mcf.path.name}) "
                f"has {topology_mcf.n_atoms}. Pass the MCF for the selected species."
            )

    blocks: list[SpeciesBlock] = []
    offset = 0
    for entry, nmol in zip(entries, nmols):
        if entry.species_index == species_index:
            mcf = topology_mcf
        else:
            if not entry.mcf_path.is_file():
                raise FileNotFoundError(
                    f"Species {entry.species_index} MCF not found: {entry.mcf_path} "
                    "(needed to locate atoms of later species in the XYZ)"
                )
            mcf = read_mcf(entry.mcf_path)
        blocks.append(
            SpeciesBlock(
                species_index=entry.species_index,
                mcf=mcf,
                nmol=int(nmol),
                atom_offset=offset,
            )
        )
        offset += blocks[-1].n_frame_atoms

    return SpeciesLayout(species=blocks)


def resolve_layout(
    *,
    traj_n_atoms: int,
    topology_mcf: McfData,
    inp_path: Path | None,
    species_index: int,
    nmols_override: list[int] | None,
) -> SpeciesLayout:
    """Resolve single- or multi-species layout and verify XYZ atom count.

    * Multi-species path: ``--inp`` with more than one ``# Molecule_Files``
      entry, or ``--species`` ≠ 1, or ``--nmols`` given.
    * Otherwise: single-species layout with
      ``nmol = n_atoms / mcf.n_atoms`` (``--inp`` may still supply T only).
    """
    use_molecule_files = False
    if nmols_override is not None or species_index != 1:
        use_molecule_files = True
    if inp_path is not None:
        entries = read_molecule_files_from_inp(inp_path)
        if len(entries) > 1:
            use_molecule_files = True

    if use_molecule_files:
        if inp_path is None:
            raise ValueError(
                "Multi-species layout needs --inp with # Molecule_Files "
                "(and optional --nmols)."
            )
        layout = build_layout_from_inp(
            inp_path,
            topology_mcf=topology_mcf,
            species_index=species_index,
            nmols_override=nmols_override,
        )
    else:
        if traj_n_atoms % topology_mcf.n_atoms != 0:
            raise ValueError(
                f"XYZ has {traj_n_atoms} atoms/frame but MCF has "
                f"{topology_mcf.n_atoms} atoms/molecule — not divisible. "
                "For multi-species movies pass --inp (and --nmols if needed)."
            )
        layout = single_species_layout(
            topology_mcf, traj_n_atoms // topology_mcf.n_atoms
        )

    if layout.n_atoms_frame != traj_n_atoms:
        detail = ", ".join(
            f"sp{s.species_index}:{s.mcf.n_atoms}×{s.nmol}" for s in layout.species
        )
        raise ValueError(
            f"XYZ has {traj_n_atoms} atoms/frame but layout implies "
            f"{layout.n_atoms_frame} ({detail}). "
            "For fixed-N multi-species, check # Molecule_Files or pass "
            "--nmols n1,n2,... matching the movie. Variable-N (GCMC) needs "
            "per-frame counts from the .H file (not yet supported)."
        )
    return layout


def parse_nmols_arg(text: str) -> list[int]:
    """Parse ``--nmols 300,300,0`` into a list of ints."""
    parts = [p.strip() for p in text.replace(" ", ",").split(",") if p.strip()]
    if not parts:
        raise ValueError("--nmols is empty")
    try:
        return [int(p) for p in parts]
    except ValueError as exc:
        raise ValueError(
            f"Bad --nmols '{text}' (expected comma-separated ints)"
        ) from exc
