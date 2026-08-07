#*******************************************************************************
#  diagnostics/io/prp.py
#
#  Shared reader for Cassandra ``.prp`` (property) files.
#
#  Why this module exists
#  ----------------------
#  Diagnostic tools (property plots, future block-average checks, etc.) all
#  need the same understanding of Cassandra's property output.  Keep that
#  logic here so each CLI script stays short and students have one place to
#  study how the file format is interpreted.
#
#  File format (see also docs/output-formats.md)
#  ---------------------------------------------
#  Cassandra writes three header lines, then one data row per property write:
#
#    # Instantaneous properties          ← title (or "# Block averages")
#             # MC_SWEEP  Energy_Total … ← column names (right-aligned)
#                      #  (kJ/mol)-Ext … ← units
#                     10  -0.362783E+04 … ← data: step, then floats
#
#  Parsing strategy: split on whitespace and convert with ``float()``.
#  Do **not** assume fixed character columns — Fortran formats (A16, E16.6)
#  control alignment for humans, but Python should treat the file as
#  whitespace-delimited.  That is robust if column widths change later.
#
#  At the end of some runs Cassandra may append ``mean`` / ``stdev`` lines
#  (from ``Write_Mean_Error``).  Those start with a non-numeric token and
#  are skipped so they are not mistaken for MC steps.
#
#  Author:  Edward J. Maginn (EJM), University of Notre Dame
#  Written: 08/07/26
#
#*******************************************************************************
"""Parse Cassandra instantaneous (or block-average) ``.prp`` property files.

Expected layout (see ``docs/output-formats.md``)::

    # Instantaneous properties
             # MC_SWEEP      Energy_Total          Pressure  ...
                      #      (kJ/mol)-Ext             (bar)  ...
                     10     -0.362783E+04      0.652585E+03  ...

Trailing ``mean`` / ``stdev`` summary lines (if present) are ignored.
Whitespace-split parsing is used; fixed column widths are not required.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


@dataclass
class PrpData:
    """In-memory representation of one Cassandra ``.prp`` file.

    Attributes
    ----------
    path :
        Source file path (for titles / error messages).
    title :
        First header line without the leading ``#``
        (e.g. ``Instantaneous properties``).
    step_label :
        Name of the independent-variable column: ``MC_SWEEP`` or ``MC_STEP``
        (depends on ``# Simulation_Length_Info`` units in the ``.inp``).
    property_names :
        Thermodynamic columns in file order (e.g. ``Energy_Total``).
    property_units :
        Parallel list of unit strings (may be empty strings if missing).
    steps :
        1-D array of sweep/step indices, shape ``(n_rows,)``.
    values :
        2-D array of property values, shape ``(n_rows, n_properties)``.
        ``values[i, j]`` is property ``j`` at row ``i``.
    """

    path: Path
    title: str
    step_label: str
    property_names: list[str]
    property_units: list[str]
    steps: np.ndarray
    values: np.ndarray

    @property
    def n_rows(self) -> int:
        """Number of property-write samples in the file."""
        return int(self.steps.shape[0])

    @property
    def n_properties(self) -> int:
        """Number of thermodynamic columns (excludes the step column)."""
        return len(self.property_names)

    def column(self, name: str) -> np.ndarray:
        """Return the 1-D series for property ``name`` (exact name)."""
        idx = self.property_names.index(name)
        return self.values[:, idx]

    def unit_for(self, name: str) -> str:
        """Return the unit string for property ``name`` (exact name)."""
        idx = self.property_names.index(name)
        return self.property_units[idx]


def _is_data_token(token: str) -> bool:
    """True if ``token`` can be parsed as a float (MC step or property value)."""
    try:
        float(token)
        return True
    except ValueError:
        return False


def read_prp(path: str | Path) -> PrpData:
    """Read a Cassandra ``.prp`` file into a :class:`PrpData` object.

    Steps
    -----
    1. Read all lines as text.
    2. Parse the three header lines (title, names, units).
    3. Convert each subsequent numeric row into a step + property vector.
    4. Skip blank lines, comment lines, and non-numeric footers.
    5. Stack rows into NumPy arrays for plotting / analysis.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Property file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Property file too short (need header + data): {path}")

    # --- header line 1: title --------------------------------------------
    title = lines[0].strip().lstrip("#").strip()

    # --- header line 2: step label + property names ----------------------
    # Example tokens after split: ['#', 'MC_SWEEP', 'Energy_Total', 'Pressure']
    name_tokens = lines[1].split()
    if not name_tokens:
        raise ValueError(f"Missing property-name header in {path}")

    if name_tokens[0] == "#":
        name_tokens = name_tokens[1:]
    if not name_tokens:
        raise ValueError(f"No step label in property-name header: {path}")

    step_label = name_tokens[0]
    property_names = name_tokens[1:]

    # --- header line 3: units (one token per property in current format) -
    # Example: ['#', '(kJ/mol)-Ext', '(bar)', '(kJ/mol)-Ext']
    unit_tokens = lines[2].split()
    if unit_tokens and unit_tokens[0] == "#":
        unit_tokens = unit_tokens[1:]
    property_units = list(unit_tokens)
    # Tolerate missing/extra unit tokens so a slightly broken file still loads
    if len(property_units) < len(property_names):
        property_units.extend([""] * (len(property_names) - len(property_units)))
    elif len(property_units) > len(property_names):
        property_units = property_units[: len(property_names)]

    # --- data rows -------------------------------------------------------
    steps: list[float] = []
    rows: list[list[float]] = []
    n_expected = 1 + len(property_names)  # step + each property

    for line_no, raw in enumerate(lines[3:], start=4):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if not tokens:
            continue
        # Footer lines like "mean ..." / "stdev ..." start with words, not numbers
        if not _is_data_token(tokens[0]):
            continue
        if len(tokens) < n_expected:
            raise ValueError(
                f"{path}:{line_no}: expected {n_expected} columns, got {len(tokens)}"
            )
        try:
            step = float(tokens[0])
            vals = [float(t) for t in tokens[1:n_expected]]
        except ValueError as exc:
            raise ValueError(f"{path}:{line_no}: could not parse numbers: {line}") from exc
        steps.append(step)
        rows.append(vals)

    if not steps:
        raise ValueError(f"No data rows found in {path}")

    return PrpData(
        path=path,
        title=title,
        step_label=step_label,
        property_names=property_names,
        property_units=property_units,
        steps=np.asarray(steps, dtype=float),
        values=np.asarray(rows, dtype=float),
    )


def resolve_property(names: Sequence[str], query: str) -> str:
    """Map a user-typed query to an exact property name in ``names``.

    Matching order
    --------------
    1. Case-insensitive exact match (``enthalpy`` → ``Enthalpy``).
    2. Otherwise, unique case-insensitive *substring* match
       (``press`` → ``Pressure`` if that is the only hit).
    3. Error if zero matches or more than one substring match
       (forces the student to be more specific).
    """
    if not query or not query.strip():
        raise ValueError("Property name is empty")

    q = query.strip()
    lower_map = {n.lower(): n for n in names}

    if q.lower() in lower_map:
        return lower_map[q.lower()]

    substr = [n for n in names if q.lower() in n.lower()]
    if len(substr) == 1:
        return substr[0]
    if len(substr) > 1:
        opts = ", ".join(substr)
        raise ValueError(f"Ambiguous property '{query}' matches: {opts}")
    available = ", ".join(names)
    raise ValueError(f"Unknown property '{query}'. Available: {available}")
