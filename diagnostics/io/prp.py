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
    """Contents of one Cassandra ``.prp`` file."""

    path: Path
    title: str
    step_label: str
    property_names: list[str]
    property_units: list[str]
    steps: np.ndarray
    values: np.ndarray  # shape (n_rows, n_properties)

    @property
    def n_rows(self) -> int:
        return int(self.steps.shape[0])

    @property
    def n_properties(self) -> int:
        return len(self.property_names)

    def column(self, name: str) -> np.ndarray:
        """Return values for ``name`` (exact match after resolve)."""
        idx = self.property_names.index(name)
        return self.values[:, idx]

    def unit_for(self, name: str) -> str:
        idx = self.property_names.index(name)
        return self.property_units[idx]


def _is_data_token(token: str) -> bool:
    try:
        float(token)
        return True
    except ValueError:
        return False


def read_prp(path: str | Path) -> PrpData:
    """Read a Cassandra ``.prp`` file into a :class:`PrpData` object."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Property file not found: {path}")

    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 3:
        raise ValueError(f"Property file too short (need header + data): {path}")

    title = lines[0].strip().lstrip("#").strip()
    name_tokens = lines[1].split()
    unit_tokens = lines[2].split()

    if not name_tokens:
        raise ValueError(f"Missing property-name header in {path}")

    # First header token is "#"; second is MC_SWEEP / MC_STEP.
    if name_tokens[0] == "#":
        name_tokens = name_tokens[1:]
    if not name_tokens:
        raise ValueError(f"No step label in property-name header: {path}")

    step_label = name_tokens[0]
    property_names = name_tokens[1:]

    # Units line: leading "#" then one unit token (or token group) per property.
    # Units may contain no spaces in current Cassandra output, e.g. (kJ/mol)-Ext.
    if unit_tokens and unit_tokens[0] == "#":
        unit_tokens = unit_tokens[1:]
    property_units = list(unit_tokens)
    if len(property_units) < len(property_names):
        property_units.extend([""] * (len(property_names) - len(property_units)))
    elif len(property_units) > len(property_names):
        property_units = property_units[: len(property_names)]

    steps: list[float] = []
    rows: list[list[float]] = []
    n_expected = 1 + len(property_names)

    for line_no, raw in enumerate(lines[3:], start=4):
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        tokens = line.split()
        if not tokens:
            continue
        # Skip mean/stdev footer lines written by Write_Mean_Error.
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
    """Match ``query`` to a property name (case-insensitive).

    Exact (case-insensitive) match preferred; otherwise unique substring match.
    Raises ``ValueError`` if zero or multiple matches.
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
