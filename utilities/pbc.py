#!/usr/bin/env python3
#*******************************************************************************
#  pbc.py
#
#  Periodic-boundary helpers shared by the restart writer and the RDF.
#
#  A molecule that straddles the box is stored as pieces on opposite faces.
#  The mass-weighted average of those stored coordinates falls near the
#  center of the box, which is not the molecule's center.  Walking the MCF
#  bond graph with minimum-image steps rebuilds one contiguous image first.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 10/24/26
#
#*******************************************************************************
"""Minimum-image displacements and contiguous molecular images."""

from __future__ import annotations

import numpy as np

from diagnostics.io.mcf import McfData


def mic_delta(delta: np.ndarray, cell: np.ndarray) -> np.ndarray:
    """Minimum-image displacement in a cell stored as a 3-by-3 H matrix."""
    fractional = np.linalg.solve(cell, np.asarray(delta, dtype=float))
    fractional -= np.rint(fractional)
    return cell @ fractional


def unwrapped_molecule(
    coords: np.ndarray,
    mcf: McfData,
    cell: np.ndarray,
) -> np.ndarray:
    """Make one molecular image contiguous by traversing its MCF bond graph."""
    n_atoms = mcf.n_atoms
    coords = np.asarray(coords, dtype=float)
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
            candidate = unwrapped[parent] + mic_delta(coords[child] - coords[parent], cell)
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
                    mic_delta(candidate - unwrapped[child], cell),
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
