#*******************************************************************************
#  diagnostics/io/hfile.py
#
#  Reader for Cassandra companion ``.H`` (movie header) files written at the
#  same ``coord_freq`` as the ``.xyz`` movie.
#
#  Per-frame layout (Write_Coords_XYZ in Src/write_properties.f90)
#  --------------------------------------------------------------
#  ::
#
#      <volume>
#      H_11 H_12 H_13
#      H_21 H_22 H_23
#      H_31 H_32 H_33
#      <blank>
#      <nspecies>
#      <species_index>  <nmol>
#      ...
#
#  ``H`` is the cell matrix (Å); for cubic boxes it is diagonal.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#*******************************************************************************
"""Parse Cassandra ``.H`` movie-header files (volume, cell matrix, nmols)."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class HFrame:
    """One frame of a ``.H`` file."""

    volume: float
    cell: np.ndarray  # (3, 3) Å
    nmols: list[int]  # nmol per species in file order (species 1..)


@dataclass
class HTrajectory:
    """All frames from one ``.H`` file."""

    path: Path
    frames: list[HFrame]

    @property
    def n_frames(self) -> int:
        return len(self.frames)

    @property
    def volumes(self) -> np.ndarray:
        return np.asarray([f.volume for f in self.frames], dtype=float)

    @property
    def cells(self) -> np.ndarray:
        """Shape ``(n_frames, 3, 3)``."""
        return np.stack([f.cell for f in self.frames], axis=0)


def read_h(path: str | Path) -> HTrajectory:
    """Read a Cassandra ``.H`` file into an :class:`HTrajectory`."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"H file not found: {path}")

    tokens: list[str] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        s = raw.strip()
        if not s:
            continue
        tokens.extend(s.split())

    frames: list[HFrame] = []
    i = 0
    ntok = len(tokens)
    while i < ntok:
        try:
            volume = float(tokens[i])
            i += 1
            cell_vals = [float(tokens[i + k]) for k in range(9)]
            i += 9
            nspecies = int(float(tokens[i]))
            i += 1
            nmols: list[int] = []
            for _ in range(nspecies):
                # species index then nmol
                _is = int(float(tokens[i]))
                nmol = int(float(tokens[i + 1]))
                i += 2
                nmols.append(nmol)
        except (IndexError, ValueError) as exc:
            raise ValueError(
                f"{path}: failed parsing .H frame {len(frames)} near token {i}: {exc}"
            ) from exc
        cell = np.asarray(cell_vals, dtype=float).reshape(3, 3)
        frames.append(HFrame(volume=volume, cell=cell, nmols=nmols))

    if not frames:
        raise ValueError(f"No frames found in {path}")
    return HTrajectory(path=path, frames=frames)


def cubic_cell(box_length: float) -> np.ndarray:
    """Return a diagonal cubic cell matrix with edge ``box_length`` (Å)."""
    L = float(box_length)
    return np.diag([L, L, L])
