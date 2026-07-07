"""Python wrapper for running Cassandra simulations."""

from .api import CassandraRunError, run_cassandra

__all__ = ["CassandraRunError", "run_cassandra"]
