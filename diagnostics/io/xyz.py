#*******************************************************************************
#  diagnostics/io/xyz.py
#
#  Shared reader for Cassandra multi-frame ``.xyz`` movie files.
#
#  Why this module exists
#  ----------------------
#  Angle / bond / displacement diagnostics all need the same understanding of
#  Cassandra's coordinate output.  Keep the parser here so each CLI stays short
#  and students have one place to study the format.
#
#  File format (see also docs/output-formats.md)
#  ---------------------------------------------
#  Standard multi-frame XYZ::
#
#      <n_atoms>
#       MC_STEP: <step>
#      element  x  y  z
#      ...
#
#  Coordinates are in Å.  Newer Cassandra writes ``F12.3``; older example
#  files may have more digits.  Always parse with whitespace ``float()`` —
#  do not assume fixed column widths.
#
#  Atom order (critical for topology tools)
#  ----------------------------------------
#  Species loop outermost, then molecule index, then atoms in **MCF order**.
#  For NVT/NPT with fixed N, molecule index is stable across frames.
#  Phase 1 of the diagnostics assumes fixed N (no ``.H`` / GCMC yet).
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#*******************************************************************************
"""Parse Cassandra multi-frame ``.xyz`` trajectory files.

Expected layout (see ``docs/output-formats.md``)::

    300
     MC_STEP:         1000
     C    1.234    5.678   -2.345
     ...

Whitespace-split parsing; fixed column widths are not required.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class XyzTrajectory:
    """In-memory representation of one Cassandra ``.xyz`` movie file.

    Attributes
    ----------
    path :
        Source file path (for titles / error messages).
    steps :
        MC step/sweep index from each frame comment, shape ``(n_frames,)``.
    elements :
        Element symbols for frame 0 (Cassandra keeps atom order fixed),
        length ``n_atoms``.
    coords :
        Cartesian coordinates in Å, shape ``(n_frames, n_atoms, 3)``.
    """

    path: Path
    steps: np.ndarray
    elements: list[str]
    coords: np.ndarray

    @property
    def n_frames(self) -> int:
        """Number of XYZ frames in the file."""
        return int(self.coords.shape[0])

    @property
    def n_atoms(self) -> int:
        """Number of atoms per frame."""
        return int(self.coords.shape[1])


def _parse_mc_step(comment: str) -> float:
    """Extract the numeric step from a comment like ``MC_STEP: 1000``.

    Falls back to ``nan`` if the comment does not contain a parseable number
    so odd third-party XYZ files still load (frames are still indexed 0..N-1).
    """
    tokens = comment.replace(":", " ").split()
    for tok in tokens:
        try:
            return float(tok)
        except ValueError:
            continue
    return float("nan")


def read_xyz(path: str | Path) -> XyzTrajectory:
    """Read a Cassandra multi-frame ``.xyz`` file into an :class:`XyzTrajectory`.

    Steps
    -----
    1. Read the file as text.
    2. For each frame: parse atom count, comment (MC step), then atom lines.
    3. Require every frame to have the same atom count (fixed-N assumption).
    4. Stack coordinates into a 3-D NumPy array for vectorized analysis.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"XYZ file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise ValueError(f"XYZ file is empty: {path}")

    steps: list[float] = []
    elements: list[str] | None = None
    frames: list[np.ndarray] = []
    i = 0
    n_lines = len(lines)
    frame_no = 0

    while i < n_lines:
        # Skip blank lines between frames (tolerant of extra whitespace)
        while i < n_lines and not lines[i].strip():
            i += 1
        if i >= n_lines:
            break

        # --- atom count --------------------------------------------------
        try:
            n_atoms = int(lines[i].split()[0])
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f"{path}: frame {frame_no}: expected atom count, got: {lines[i]!r}"
            ) from exc
        i += 1
        if i >= n_lines:
            raise ValueError(f"{path}: frame {frame_no}: missing comment line")

        # --- comment / MC_STEP -------------------------------------------
        comment = lines[i].strip()
        steps.append(_parse_mc_step(comment))
        i += 1

        if i + n_atoms > n_lines:
            raise ValueError(
                f"{path}: frame {frame_no}: need {n_atoms} atom lines, file ended early"
            )

        # --- atom coordinates --------------------------------------------
        frame_elems: list[str] = []
        xyz = np.empty((n_atoms, 3), dtype=float)
        for ia in range(n_atoms):
            tokens = lines[i + ia].split()
            if len(tokens) < 4:
                raise ValueError(
                    f"{path}: frame {frame_no}, atom {ia + 1}: "
                    f"expected 'element x y z', got: {lines[i + ia]!r}"
                )
            frame_elems.append(tokens[0])
            try:
                xyz[ia, 0] = float(tokens[1])
                xyz[ia, 1] = float(tokens[2])
                xyz[ia, 2] = float(tokens[3])
            except ValueError as exc:
                raise ValueError(
                    f"{path}: frame {frame_no}, atom {ia + 1}: "
                    f"could not parse coordinates: {lines[i + ia]!r}"
                ) from exc
        i += n_atoms

        if elements is None:
            elements = frame_elems
        elif len(frame_elems) != len(elements):
            raise ValueError(
                f"{path}: frame {frame_no}: atom count changed "
                f"({len(elements)} → {len(frame_elems)}); "
                "GCMC / variable-N trajectories need the companion .H file "
                "(not supported in Phase 1)."
            )

        frames.append(xyz)
        frame_no += 1

    if not frames or elements is None:
        raise ValueError(f"No XYZ frames found in {path}")

    return XyzTrajectory(
        path=path,
        steps=np.asarray(steps, dtype=float),
        elements=elements,
        coords=np.stack(frames, axis=0),
    )
