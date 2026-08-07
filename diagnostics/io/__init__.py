"""Shared I/O helpers for Cassandra diagnostic tools."""

from .prp import PrpData, read_prp, resolve_property

__all__ = ["PrpData", "read_prp", "resolve_property"]
