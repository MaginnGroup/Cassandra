#!/usr/bin/env python3
#*******************************************************************************
#  rdf.py
#
#  Radial distribution function and coordination number from a Cassandra
#  trajectory.
#
#  What is computed
#  ---------------
#  g(r) is the pair density at distance r, divided by the density of an
#  ideal gas at the same number of neighbors.  The coordination number is
#  how many neighbors sit inside a cutoff:
#
#      n(R) = integral_0^R  4 pi r^2 rho g(r) dr
#
#  The integral is a sum over spherical shells.  A shell from r to r+dr has
#  volume (4/3) pi [(r+dr)^3 - r^3], which is the exact integral of 4 pi r^2
#  over that shell.  For a thin shell this is 4 pi r^2 dr.
#
#  Two pair choices
#  ----------------
#  --mode com
#      Centers of mass of one species.  Each molecule is rebuilt across the
#      periodic boundary along its MCF bonds before the center is taken.
#      Averaging the stored coordinates instead would put a straddling
#      molecule near the middle of the box.
#  --mode atom
#      Atom–atom pairs inside one species.  Each side is chosen by element
#      (C, F, H), by an exact MCF atom type such as c3_1_s1, or by a 1-based
#      MCF atom index.  Atom type is the preferred selection because it
#      distinguishes sites that share an element.  Omit the second selector
#      to pair a selection with itself.
#
#  Intramolecular pairs are left out unless --include-intramolecular is set.
#  The report states the selections and their sites per molecule.  Atom
#  distances use the same per-frame triclinic minimum-image calculation as
#  COM distances.  They do not require molecule reconstruction: an atom is
#  already a point, whereas a wrapped molecule must be reconstructed before
#  calculating its COM.
#
#  Normalization
#  -------------
#  All-site selections use rho_j = N_j / V; a selection paired with itself
#  uses N-1 rather than N because a particle is not a neighbor of itself,
#  and each pair is stored once.  With that choice g(r) approaches 1 at
#  large r in an ideal gas.  N and V are taken from the frame, so an NPT
#  movie can change box size from frame to frame.  When intramolecular pairs
#  are excluded, the atom calculation uses the exact finite-system count of
#  sites on different molecules in both g(r) and n(r), rather than including
#  omitted same-molecule pairs in the ideal denominator.
#
#  Coordination cutoff
#  -------------------
#  --r-cut R integrates through the last bin whose center is at or below R.
#  Omit --r-cut and the script estimates R from g(r).  It finds the second
#  and third crossings of g(r)=1: down from the first peak, then up from the
#  first-shell valley.  The lowest bin between those crossings is the first
#  minimum and the coordination cutoff.  The report names this choice.  Pass
#  --r-cut when the printed minimum is not the shell you want.  If the three
#  crossings cannot be found, the script stops and asks for --r-cut.
#
#  Uncertainty
#  -----------
#  The standard error is the frame-to-frame scatter divided by sqrt(Nframes).
#  Neighboring frames are correlated, so this error bar checks that the
#  histogram is stable.  It is not a production uncertainty.
#
#  The script refuses a changing molecule count (GCMC).  An NPT movie is
#  accepted when the .H file has one cell per XYZ frame.
#
#  Output files
#  ------------
#  --r-max is how far g(r) is tabulated.  The default is half the shortest
#  box edge.  --r-cut is only the upper limit of the coordination integral.
#  Pass --plot to also write a figure: rdf, cn, or both on opposite axes.
#  The RDF figure does not mark the coordination cutoff.  Pass --mark-r-cut
#  to draw that line.  The CN and combined figures always mark it.  Default
#  titles name the pair, for example "COM-COM RDF, species 1" and
#  "COM-COM CN, species 1".  Replace the figure you are writing with
#  --title-rdf, --title-cn, or --title-both.
#  --legend is none, curves, or with-cutoff.  The default is none for an
#  RDF or CN figure, and with-cutoff for the combined figure, which names
#  g(r), n(r), and the cutoff.  The uncertainty band is drawn but never
#  listed.  --x-unit is nm by default, or A for angstroms, drawn as Å.
#  --r-cut, --r-max, and the table stay in angstroms.  The CN figure axis
#  runs from zero to the coordination cutoff.  The RDF and combined figures
#  run to --r-max.
#  The file names use the movie stem and the pair.  For a COM run of
#  [run_name].xyz they are:
#
#      [run_name]_com.rdf
#      [run_name]_com_rdf.png
#      [run_name]_com_cn.png
#      [run_name]_com_rdf_cn.png
#
#  An atom RDF names the pair with the MCF atom types.  Pass --atom-type-i,
#  and --atom-type-j when the second type is different.  c3_1_s1 with itself
#  uses c3_1_s1-c3_1_s1.  c3_1_s1 with f_2_s1 uses c3_1_s1-f_2_s1.  The files
#  are written next to the movie.  --out replaces the table path.
#
#  Example
#  -------
#      python utilities/rdf.py movie.xyz species.mcf --h movie.H --plot both
#      python utilities/rdf.py movie.xyz species.mcf --h movie.H --mode atom \
#          --atom-type-i c3_1_s1 --atom-type-j f_2_s1 --plot rdf
#
#  The first command is a center-of-mass run.  The second uses atom types
#  from the MCF.  The table and figure are written next to the movie.
#  --r-max is half the shortest box edge, and the coordination cutoff is the
#  lowest bin between the second and third crossings of g(r)=1.
#  The sections above name the other choices: the pair, the cutoffs, the
#  figure, the title, the legend, and the x-axis unit.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/24/26
#
#*******************************************************************************
"""COM and atom–atom radial distribution functions, and their coordination numbers."""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from diagnostics.io.hfile import read_h
from diagnostics.io.layout import SpeciesLayout, resolve_layout
from diagnostics.io.mcf import McfData, read_mcf
from diagnostics.io.xyz import read_xyz
from utilities.pbc import unwrapped_molecule


@dataclass(frozen=True)
class RdfResult:
    """One averaged g(r) and the coordination number accumulated from it."""

    r_center: np.ndarray
    g: np.ndarray
    g_sem: np.ndarray
    mean_shell_count: np.ndarray
    n_running: np.ndarray
    r_cut: float
    cutoff_note: str
    coordination: float
    n_frames: int
    exclusion_note: str
    normalization_note: str
    pair_note: str

    def table(self) -> str:
        """Comment header plus one row per bin."""
        lines = [
            f"# {self.pair_note}",
            f"# {self.exclusion_note}",
            f"# {self.normalization_note}",
            f"# frames: {self.n_frames}",
            (
                "# Frame SEM treats frames as independent. Neighboring frames are "
                "correlated, so this is not a production uncertainty."
            ),
            f"# coordination cutoff: {self.r_cut:.6g} A ({self.cutoff_note})",
            f"# coordination number: {self.coordination:.8g}",
            f"# {'r_A':>14} {'g':>14} {'g_sem':>14} {'mean_shell_count':>18} {'n_running':>14}",
        ]
        for r, g, sem, count, n_run in zip(
            self.r_center, self.g, self.g_sem, self.mean_shell_count, self.n_running
        ):
            lines.append(
                f"  {r:14.6f} {g:14.6f} {sem:14.6f} {count:18.6f} {n_run:14.6f}"
            )
        return "\n".join(lines) + "\n"


def shell_volumes(edges: np.ndarray) -> np.ndarray:
    """Exact spherical-shell volumes for bins bounded by ``edges`` (Å)."""
    outer = edges[1:]
    inner = edges[:-1]
    return (4.0 / 3.0) * np.pi * (outer**3 - inner**3)


def bin_edges(r_max: float, bin_width: float) -> np.ndarray:
    """Bin edges from 0 through the last edge that does not pass ``r_max``."""
    if r_max <= 0.0 or bin_width <= 0.0:
        raise ValueError("r_max and bin width must be positive")
    n_bins = int(np.floor(r_max / bin_width))
    if n_bins < 1:
        raise ValueError("bin width is larger than r_max, so there are no bins")
    return np.arange(n_bins + 1, dtype=float) * float(bin_width)


def running_coordination(g: np.ndarray, edges: np.ndarray, rho: float) -> np.ndarray:
    """Cumulative n(R) at each bin's outer edge for number density ``rho``."""
    return np.cumsum(np.asarray(g, dtype=float) * float(rho) * shell_volumes(edges))


def first_minimum_after_peak(r_center: np.ndarray, g: np.ndarray) -> float | None:
    """Lowest-bin center between the second and third crossings of ``g(r)=1``.

    The first crossing rises into the first peak, the second falls out of it,
    and the third rises into the second shell.  Using the lowest bin in the
    interval bracketed by the second and third crossings avoids accepting a
    small downward fluctuation on the shoulder of the first peak.
    """
    values = np.asarray(g, dtype=float)
    if values.size < 2:
        return None
    rising = np.flatnonzero((values[:-1] <= 1.0) & (values[1:] > 1.0)) + 1
    if rising.size == 0:
        return None
    first_rising = int(rising[0])
    falling = np.flatnonzero(
        (values[:-1] > 1.0) & (values[1:] <= 1.0)
    ) + 1
    falling = falling[falling > first_rising]
    if falling.size == 0:
        return None
    second_crossing = int(falling[0])
    third_candidates = rising[rising > second_crossing]
    if third_candidates.size == 0:
        return None
    third_crossing = int(third_candidates[0])
    valley = values[second_crossing : third_crossing + 1]
    minimum = second_crossing + int(np.argmin(valley))
    return float(r_center[minimum])


def _shortest_edge(cell: np.ndarray) -> float:
    return float(np.min(np.linalg.norm(cell, axis=0)))


def _pair_histogram(
    points_a: np.ndarray,
    mol_a: np.ndarray,
    points_b: np.ndarray,
    mol_b: np.ndarray,
    cell: np.ndarray,
    edges: np.ndarray,
    *,
    like: bool,
    exclude_intramolecular: bool,
) -> np.ndarray:
    """Histogram minimum-image pair distances.

    Like selections count each unordered pair once.  Unlike selections count
    every A–B pair once.  Same-molecule pairs are dropped when requested.
    """
    n_a = points_a.shape[0]
    n_b = points_b.shape[0]
    counts = np.zeros(len(edges) - 1, dtype=np.int64)
    if n_a == 0 or n_b == 0:
        return counts
    block = 128
    r_max = float(edges[-1])
    for i0 in range(0, n_a, block):
        i1 = min(n_a, i0 + block)
        delta = points_b[None, :, :] - points_a[i0:i1, None, :]
        flat = np.moveaxis(delta, -1, 0).reshape(3, -1)
        fractional = np.linalg.solve(cell, flat)
        fractional -= np.rint(fractional)
        mic = np.moveaxis((cell @ fractional).reshape(3, i1 - i0, n_b), 0, -1)
        dist = np.linalg.norm(mic, axis=-1)
        mask = dist < r_max
        if exclude_intramolecular:
            mask &= mol_a[i0:i1, None] != mol_b[None, :]
        if like:
            # Identical lists are passed twice. Keep the upper triangle only.
            row = np.arange(i0, i1)[:, None]
            col = np.arange(n_b)[None, :]
            mask &= col > row
        if np.any(mask):
            counts += np.histogram(dist[mask], bins=edges)[0]
    return counts


def _frame_g(
    counts: np.ndarray,
    shells: np.ndarray,
    volume: float,
    *,
    nmol: int,
    sites_i: int,
    sites_j: int,
    like: bool,
    exclude_intramolecular: bool,
) -> tuple[np.ndarray, float]:
    """g(r) and its neighbor density, with an exact finite-system denominator."""
    n_i = nmol * sites_i
    n_j = nmol * sites_j
    if exclude_intramolecular:
        if like:
            pair_count = nmol * (nmol - 1) / 2.0 * sites_i**2
            rho = (nmol - 1) * sites_i / volume
        else:
            pair_count = nmol * (nmol - 1) * sites_i * sites_j
            rho = (nmol - 1) * sites_j / volume
    elif like:
        pair_count = n_i * (n_i - 1) / 2.0
        rho = (n_i - 1) / volume
    else:
        pair_count = n_i * n_j
        rho = n_j / volume
    ideal = pair_count * shells / volume
    g = np.zeros(counts.shape, dtype=float)
    positive = ideal > 0.0
    g[positive] = counts[positive] / ideal[positive]
    return g, float(rho)


def molecule_coms(frame: np.ndarray, layout: SpeciesLayout, species: int, cell: np.ndarray) -> np.ndarray:
    """COM of every molecule of ``species`` after PBC reconstruction."""
    block = layout.block(species)
    masses = np.asarray([atom.mass for atom in block.mcf.atoms], dtype=float)
    if float(masses.sum()) <= 0.0:
        raise ValueError(f"{block.mcf.path}: molecule has nonpositive total mass")
    coms = np.empty((block.nmol, 3), dtype=float)
    for molecule in range(1, block.nmol + 1):
        coords = layout.molecule_coords(frame, species, molecule)
        contiguous = unwrapped_molecule(coords, block.mcf, cell)
        coms[molecule - 1] = np.average(contiguous, axis=0, weights=masses)
    return coms


def select_atom_indices(
    mcf: McfData,
    *,
    element: str | None,
    atom_type: str | None,
    atom_index: int | None,
    label: str,
) -> np.ndarray:
    """0-based MCF indices matching exactly one of the three selectors."""
    chosen = [element is not None, atom_type is not None, atom_index is not None]
    if sum(chosen) != 1:
        raise ValueError(f"{label}: pass exactly one of element, atom type, or atom index")
    if atom_index is not None:
        if atom_index < 1 or atom_index > mcf.n_atoms:
            raise ValueError(
                f"{label}: atom index {atom_index} is outside 1..{mcf.n_atoms}"
            )
        return np.asarray([atom_index - 1], dtype=int)
    if element is not None:
        key = element.strip().lower()
        hits = [i for i, atom in enumerate(mcf.atoms) if atom.element.strip().lower() == key]
        if not hits:
            known = ", ".join(sorted({atom.element.strip() for atom in mcf.atoms}))
            raise ValueError(f"{label}: element {element!r} is not in the MCF ({known})")
        return np.asarray(hits, dtype=int)
    key = atom_type.strip()
    hits = [i for i, atom in enumerate(mcf.atoms) if atom.atom_type.strip() == key]
    if not hits:
        known = ", ".join(atom.atom_type.strip() for atom in mcf.atoms)
        raise ValueError(f"{label}: atom type {atom_type!r} is not in the MCF ({known})")
    return np.asarray(hits, dtype=int)


def selector_description(
    *,
    element: str | None,
    atom_type: str | None,
    atom_index: int | None,
) -> str:
    """Human-readable atom selector for the output-table header."""
    if atom_type is not None:
        return f"atom type {atom_type.strip()!r}"
    if element is not None:
        return f"element {element.strip()!r}"
    assert atom_index is not None
    return f"MCF atom index {atom_index}"


def _gather_atoms(
    frame: np.ndarray,
    layout: SpeciesLayout,
    species: int,
    local_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    block = layout.block(species)
    n_pick = len(local_indices)
    points = np.empty((block.nmol * n_pick, 3), dtype=float)
    mols = np.empty(block.nmol * n_pick, dtype=int)
    for molecule in range(block.nmol):
        coords = layout.molecule_coords(frame, species, molecule + 1)
        start = molecule * n_pick
        points[start : start + n_pick] = coords[local_indices]
        mols[start : start + n_pick] = molecule
    return points, mols


def _cells_for_frames(h_path: Path, n_frames: int) -> tuple[list[np.ndarray], list[int]]:
    trajectory = read_h(h_path)
    if trajectory.n_frames != n_frames:
        raise ValueError(
            f"{h_path} has {trajectory.n_frames} frames but the XYZ has {n_frames}"
        )
    nmols = [tuple(frame.nmols) for frame in trajectory.frames]
    if len(set(nmols)) != 1:
        raise ValueError(
            f"{h_path} changes molecule count between frames. "
            "GCMC trajectories are not supported."
        )
    cells = [np.asarray(frame.cell, dtype=float) for frame in trajectory.frames]
    return cells, list(nmols[0])


def compute_rdf(
    xyz_path: Path,
    mcf_path: Path,
    h_path: Path,
    *,
    mode: str,
    species: int = 1,
    inp_path: Path | None = None,
    element_i: str | None = None,
    element_j: str | None = None,
    atom_type_i: str | None = None,
    atom_type_j: str | None = None,
    atom_index_i: int | None = None,
    atom_index_j: int | None = None,
    bin_width: float = 0.1,
    r_max: float | None = None,
    r_cut: float | None = None,
    exclude_intramolecular: bool = True,
) -> RdfResult:
    """Average g(r) over the movie and integrate it to a coordination number."""
    if mode not in {"com", "atom"}:
        raise ValueError("mode must be 'com' or 'atom'")
    traj = read_xyz(xyz_path)
    mcf = read_mcf(mcf_path)
    layout = resolve_layout(
        traj_n_atoms=traj.n_atoms,
        topology_mcf=mcf,
        inp_path=inp_path,
        species_index=species,
        nmols_override=None,
    )
    cells, h_nmols = _cells_for_frames(h_path, traj.n_frames)
    block = layout.block(species)
    if species > len(h_nmols) or h_nmols[species - 1] != block.nmol:
        raise ValueError(
            f"{h_path} molecule count for species {species} does not match the XYZ layout"
        )
    if r_max is None:
        r_max = 0.5 * min(_shortest_edge(cell) for cell in cells)
    edges = bin_edges(r_max, bin_width)
    shells = shell_volumes(edges)
    like = True
    local_i = local_j = None
    selection_i_note = selection_j_note = ""
    if mode == "atom":
        local_i = select_atom_indices(
            layout.block(species).mcf,
            element=element_i,
            atom_type=atom_type_i,
            atom_index=atom_index_i,
            label="selection i",
        )
        selection_i_note = selector_description(
            element=element_i, atom_type=atom_type_i, atom_index=atom_index_i
        )
        j_specified = any(value is not None for value in (element_j, atom_type_j, atom_index_j))
        if j_specified:
            local_j = select_atom_indices(
                layout.block(species).mcf,
                element=element_j,
                atom_type=atom_type_j,
                atom_index=atom_index_j,
                label="selection j",
            )
            selection_j_note = selector_description(
                element=element_j, atom_type=atom_type_j, atom_index=atom_index_j
            )
        else:
            local_j = local_i
            selection_j_note = selection_i_note
        like = np.array_equal(local_i, local_j)

    g_frames = []
    count_frames = []
    n_frames_running = []
    for iframe, cell in enumerate(cells):
        frame = traj.coords[iframe]
        volume = abs(float(np.linalg.det(cell)))
        if volume <= 0.0:
            raise ValueError(f"Frame {iframe} has a non-positive cell volume")
        if mode == "com":
            points = molecule_coms(frame, layout, species, cell)
            mols = np.arange(points.shape[0])
            counts = _pair_histogram(
                points, mols, points, mols, cell, edges,
                like=True, exclude_intramolecular=True,
            )
            sites_i = sites_j = 1
            like_frame = True
            exclude_frame = True
        else:
            assert local_i is not None and local_j is not None
            points_i, mol_i = _gather_atoms(frame, layout, species, local_i)
            points_j, mol_j = _gather_atoms(frame, layout, species, local_j)
            counts = _pair_histogram(
                points_i, mol_i, points_j, mol_j, cell, edges,
                like=like, exclude_intramolecular=exclude_intramolecular,
            )
            sites_i, sites_j = len(local_i), len(local_j)
            like_frame = like
            exclude_frame = exclude_intramolecular
        g_frame, rho = _frame_g(
            counts,
            shells,
            volume,
            nmol=block.nmol,
            sites_i=sites_i,
            sites_j=sites_j,
            like=like_frame,
            exclude_intramolecular=exclude_frame,
        )
        g_frames.append(g_frame)
        count_frames.append(counts.astype(float))
        n_frames_running.append(running_coordination(g_frame, edges, rho))

    g_stack = np.vstack(g_frames)
    g_mean = np.mean(g_stack, axis=0)
    if traj.n_frames > 1:
        g_sem = np.std(g_stack, axis=0, ddof=1) / np.sqrt(traj.n_frames)
    else:
        g_sem = np.zeros_like(g_mean)
    r_center = 0.5 * (edges[:-1] + edges[1:])
    n_running = np.mean(np.vstack(n_frames_running), axis=0)
    if r_cut is None:
        found = first_minimum_after_peak(r_center, g_mean)
        if found is None:
            raise ValueError(
                "No first-shell valley bracketed by the second and third g(r)=1 "
                "crossings. Pass --r-cut to choose the coordination cutoff."
            )
        cutoff = found
        cutoff_note = "minimum between the second and third g(r)=1 crossings"
    else:
        cutoff = float(r_cut)
        cutoff_note = "user --r-cut"
    included = r_center <= cutoff
    if not np.any(included):
        raise ValueError(f"r-cut {cutoff} A is below the first bin center")
    coordination = float(n_running[included][-1])
    if mode == "com":
        pair_note = f"COM-COM RDF, species {species}"
        exclusion_note = "Pairs are different molecules only."
        normalization_note = (
            "Intermolecular normalization: 1 COM per molecule; "
            f"{block.nmol * (block.nmol - 1) / 2:.0f} candidate pairs per frame."
        )
    else:
        pair_note = (
            f"Atom-atom RDF, species {species}: "
            f"{selection_i_note} - {selection_j_note}"
        )
        if like:
            pair_note += " (same selection)"
        else:
            pair_note += " (two selections)"
        exclusion_note = (
            "Pairs inside one molecule are excluded."
            if exclude_intramolecular
            else "Pairs inside one molecule are included."
        )
        if exclude_intramolecular:
            if like:
                candidate_pairs = block.nmol * (block.nmol - 1) / 2 * len(local_i) ** 2
            else:
                candidate_pairs = block.nmol * (block.nmol - 1) * len(local_i) * len(local_j)
            normalization_note = (
                "Intermolecular normalization: "
                f"{len(local_i)} and {len(local_j)} selected site(s) per molecule; "
                f"{candidate_pairs:.0f} candidate pairs per frame."
            )
        else:
            normalization_note = (
                "All-site normalization: "
                f"{len(local_i)} and {len(local_j)} selected site(s) per molecule."
            )
    return RdfResult(
        r_center=r_center,
        g=g_mean,
        g_sem=g_sem,
        mean_shell_count=np.mean(np.vstack(count_frames), axis=0),
        n_running=n_running,
        r_cut=cutoff,
        cutoff_note=cutoff_note,
        coordination=coordination,
        n_frames=traj.n_frames,
        exclusion_note=exclusion_note,
        normalization_note=normalization_note,
        pair_note=pair_note,
    )


def selection_tag(
    mode: str,
    *,
    element_i: str | None,
    element_j: str | None,
    atom_type_i: str | None,
    atom_type_j: str | None,
    atom_index_i: int | None,
    atom_index_j: int | None,
) -> str:
    """Short file-name piece: ``com``, ``c3_1_s1-c3_1_s1``, ``c3_1_s1-f_2_s1``."""
    if mode == "com":
        return "com"

    def one(element: str | None, atom_type: str | None, atom_index: int | None) -> str | None:
        if element:
            return element
        if atom_type:
            return atom_type
        if atom_index is not None:
            return f"a{atom_index}"
        return None

    left = one(element_i, atom_type_i, atom_index_i) or "atom"
    right = one(element_j, atom_type_j, atom_index_j) or left
    raw = f"{left}-{right}"
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in raw)


def output_paths(xyz: Path, tag: str, plot: str | None) -> tuple[Path, Path | None]:
    """Table and figure paths beside the movie.

    ``plot`` is ``rdf``, ``cn``, ``both``, or None.
    """
    stem = xyz.resolve().parent / f"{xyz.stem}_{tag}"
    table = stem.with_suffix(".rdf")
    suffix = {"rdf": "_rdf.png", "cn": "_cn.png", "both": "_rdf_cn.png"}.get(plot or "")
    figure = Path(str(stem) + suffix) if suffix else None
    return table, figure


def normalize_x_unit(unit: str) -> str:
    """Return ``nm`` or ``A``. ``A`` is drawn as Å."""
    key = unit.strip()
    if key == "nm":
        return "nm"
    if key in {"A", "angstrom"}:
        return "A"
    raise ValueError("x unit must be nm or A")


def distance_for_axis(r_angstrom: float | np.ndarray, unit: str) -> np.ndarray:
    """Convert an angstrom distance onto the plot axis."""
    values = np.asarray(r_angstrom, dtype=float)
    if normalize_x_unit(unit) == "nm":
        return values / 10.0
    return values


def distance_axis_label(unit: str) -> str:
    """X-axis label. Angstroms use Å."""
    if normalize_x_unit(unit) == "nm":
        return "r (nm)"
    return "r (Å)"


def cutoff_label(r_cut_angstrom: float, unit: str) -> str:
    """Legend text for the coordination cutoff in the plot unit."""
    value = float(distance_for_axis(r_cut_angstrom, unit))
    symbol = "nm" if normalize_x_unit(unit) == "nm" else "Å"
    return f"r_cut {value:.3g} {symbol}"


def resolve_legend(kind: str, choice: str | None) -> str:
    """Legend choice: none, curves, or with-cutoff.

    RDF and CN figures have no legend unless one is requested.  The combined
    figure names g(r), n(r), and the cutoff.
    """
    if kind not in {"rdf", "cn", "both"}:
        raise ValueError("plot kind must be rdf, cn, or both")
    if choice is None:
        return "with-cutoff" if kind == "both" else "none"
    if choice not in {"none", "curves", "with-cutoff"}:
        raise ValueError("legend must be none, curves, or with-cutoff")
    return choice


def plot_title(pair_note: str, kind: str, override: str | None = None) -> str:
    """Title for one figure. ``override`` replaces the default for that figure."""
    if override:
        return override
    title = pair_note
    for suffix in (" (same selection)", " (two selections)"):
        if title.endswith(suffix):
            title = title[: -len(suffix)]
    if kind == "cn":
        return title.replace(" RDF", " CN", 1)
    return title


def write_plot(
    result: RdfResult,
    path: Path,
    kind: str,
    mark_cutoff: bool = False,
    title: str | None = None,
    x_unit: str = "nm",
    legend: str | None = None,
) -> None:
    """Save g(r), n(r), or both on opposite axes. ``kind`` is rdf, cn, or both.

    The RDF figure omits the coordination cutoff unless ``mark_cutoff`` is set.
    The CN and combined figures always draw it.  Distances on the axis use
    ``x_unit`` (nm or A).  The uncertainty band is never named in the legend.
    """
    if kind not in {"rdf", "cn", "both"}:
        raise ValueError("plot kind must be rdf, cn, or both")
    unit = normalize_x_unit(x_unit)
    legend_choice = resolve_legend(kind, legend)
    name_curves = legend_choice in {"curves", "with-cutoff"}
    name_cutoff = legend_choice == "with-cutoff"

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = distance_for_axis(result.r_center, unit)
    fig, ax_g = plt.subplots(figsize=(8.2, 4.6))
    ax_n = None
    if kind != "rdf" or mark_cutoff:
        cutoff_name = cutoff_label(result.r_cut, unit) if name_cutoff else None
        ax_g.axvline(float(distance_for_axis(result.r_cut, unit)), color="0.45", lw=1.0, ls="--", label=cutoff_name)
    if kind in {"rdf", "both"}:
        ax_g.fill_between(
            r,
            result.g - result.g_sem,
            result.g + result.g_sem,
            color="#0072b2",
            alpha=0.25,
        )
        ax_g.plot(r, result.g, color="#0072b2", lw=1.6, label="g(r)" if name_curves else None)
        ax_g.axhline(1.0, color="#0072b2", lw=0.7, alpha=0.5)
        ax_g.set_ylabel("g(r)")
    if kind == "cn":
        ax_g.plot(r, result.n_running, color="#d55e00", lw=1.6, label="n(r)" if name_curves else None)
        ax_g.set_ylabel("n(r)")
    if kind == "both":
        ax_n = ax_g.twinx()
        ax_n.plot(r, result.n_running, color="#d55e00", lw=1.6, label="n(r)" if name_curves else None)
        ax_n.set_ylabel("n(r)")
    if legend_choice != "none":
        lines, labels = ax_g.get_legend_handles_labels()
        if ax_n is not None:
            extra_lines, extra_labels = ax_n.get_legend_handles_labels()
            lines += extra_lines
            labels += extra_labels
        if labels:
            ax_g.legend(lines, labels, frameon=True, framealpha=0.92, fontsize=9)
    ax_g.set_xlabel(distance_axis_label(unit))
    if kind == "cn":
        x_right = float(distance_for_axis(result.r_cut, unit))
        visible = r <= x_right + 1e-9
        y_top = float(np.max(result.n_running[visible])) if np.any(visible) else 1.0
        ax_g.set_ylim(0.0, y_top * 1.05 if y_top > 0.0 else 1.0)
    else:
        x_right = float(r[-1]) if len(r) else 1.0
    ax_g.set_xlim(0.0, x_right if x_right > 0.0 else 1.0)
    ax_g.set_title(plot_title(result.pair_note, kind, title))
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def _optional_text(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Radial distribution function and coordination number "
            "from a Cassandra XYZ movie and its .H file."
        )
    )
    parser.add_argument("xyz", type=Path, help="Coordinate movie")
    parser.add_argument("mcf", type=Path, help="MCF for the selected species")
    parser.add_argument("--h", type=Path, required=True, dest="h_path", help="Companion .H movie")
    parser.add_argument("--inp", type=Path, help="Input file, required for a multi-species movie")
    parser.add_argument("--mode", choices=("com", "atom"), default="com")
    parser.add_argument("--species", type=int, default=1, help="1-based species index")
    parser.add_argument("--element-i")
    parser.add_argument("--element-j")
    parser.add_argument("--atom-type-i")
    parser.add_argument("--atom-type-j")
    parser.add_argument("--atom-index-i", type=int)
    parser.add_argument("--atom-index-j", type=int)
    parser.add_argument("--bin-width", type=float, default=0.1, help="Bin width in Angstroms")
    parser.add_argument("--r-max", type=float, help="Histogram limit in Angstroms (default: half the shortest edge)")
    parser.add_argument("--r-cut", type=float, help="Coordination cutoff in Angstroms")
    parser.add_argument(
        "--include-intramolecular",
        action="store_true",
        help="Count atom pairs inside the same molecule. COM pairs stay intermolecular.",
    )
    parser.add_argument("--out", type=Path, help="Table path. Default when plotting: <movie>_<pair>.rdf beside the movie")
    parser.add_argument(
        "--plot",
        choices=("rdf", "cn", "both"),
        help="Write a PNG: g(r), n(r), or both on opposite axes",
    )
    parser.add_argument(
        "--mark-r-cut",
        action="store_true",
        help="Draw the coordination cutoff on an RDF figure. CN and combined figures always mark it.",
    )
    parser.add_argument(
        "--title-rdf",
        help="Title for --plot rdf. Default names the pair, for example 'COM-COM RDF, species 1'.",
    )
    parser.add_argument(
        "--title-cn",
        help="Title for --plot cn. Default names the pair, for example 'COM-COM CN, species 1'.",
    )
    parser.add_argument(
        "--title-both",
        help="Title for --plot both. Default matches the RDF title.",
    )
    parser.add_argument(
        "--legend",
        choices=("none", "curves", "with-cutoff"),
        help=(
            "Legend. Default: none for RDF and CN; with-cutoff for the combined "
            "figure, which names g(r), n(r), and the cutoff. The uncertainty band is never listed."
        ),
    )
    parser.add_argument(
        "--x-unit",
        choices=("nm", "A"),
        default="nm",
        help="X-axis unit. Default nm. A is drawn as Å. Cutoffs and the table stay in angstroms.",
    )
    args = parser.parse_args(argv)
    result = compute_rdf(
        args.xyz,
        args.mcf,
        args.h_path,
        mode=args.mode,
        species=args.species,
        inp_path=args.inp,
        element_i=_optional_text(args.element_i),
        element_j=_optional_text(args.element_j),
        atom_type_i=_optional_text(args.atom_type_i),
        atom_type_j=_optional_text(args.atom_type_j),
        atom_index_i=args.atom_index_i,
        atom_index_j=args.atom_index_j,
        bin_width=args.bin_width,
        r_max=args.r_max,
        r_cut=args.r_cut,
        exclude_intramolecular=not args.include_intramolecular,
    )
    text = result.table()
    sys.stdout.write(text)
    tag = selection_tag(
        args.mode,
        element_i=_optional_text(args.element_i),
        element_j=_optional_text(args.element_j),
        atom_type_i=_optional_text(args.atom_type_i),
        atom_type_j=_optional_text(args.atom_type_j),
        atom_index_i=args.atom_index_i,
        atom_index_j=args.atom_index_j,
    )
    table_path, figure_path = output_paths(args.xyz, tag, args.plot)
    if args.out is not None:
        table_path = args.out
    if args.out is not None or args.plot is not None:
        table_path.write_text(text, encoding="utf-8")
        print(f"Wrote {table_path}", file=sys.stderr)
    if figure_path is not None:
        title_for = {
            "rdf": _optional_text(args.title_rdf),
            "cn": _optional_text(args.title_cn),
            "both": _optional_text(args.title_both),
        }[args.plot]
        write_plot(
            result,
            figure_path,
            args.plot,
            mark_cutoff=args.mark_r_cut,
            title=title_for,
            x_unit=args.x_unit,
            legend=args.legend,
        )
        print(f"Wrote {figure_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
