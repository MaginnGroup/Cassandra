"""Python entry point for running Cassandra Monte Carlo simulations.

Cassandra is a Fortran executable, not a Python library. This module launches it
as a subprocess with the same interface you would use from a shell::

    cassandra.exe path/to/run.inp

The Fortran driver (``Src/main.f90``) reads the input file path from the first
command-line argument and writes output files relative to the *current working
directory*, using the ``Run_Name`` field in the ``.inp`` file (e.g.
``water_spc_nvt.out.log``, ``water_spc_nvt.out.prp``).

Because of that I/O model, ``run_cassandra`` always sets ``cwd`` to the
directory containing the input file so outputs land alongside the example data.

Typical workflow
----------------
1. Compile Cassandra: ``make`` in ``Src/`` (produces ``cassandra_gfortran.exe``).
2. Run the smoke test from the repo root: ``python run_test.py``.
3. Call ``run_cassandra()`` from your own scripts (see ``python/README.md``).

This module is separate from the legacy test harness in ``Scripts/testSuite/``,
which invokes Cassandra with its own ``subprocess`` calls.
"""

from __future__ import annotations

import os
import subprocess


class CassandraRunError(RuntimeError):
    """Raised when Cassandra exits with a non-zero status or cannot be started."""


def run_cassandra(
    input_file_path: str,
    cassandra_exe_path: str = "cassandra.exe",
    *,
    raise_on_error: bool = True,
    verbose: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a Cassandra simulation as a subprocess.

    Args:
        input_file_path: Path to the ``.inp`` file
            (e.g. ``"Examples/NVT/water_spc/nvt.inp"``).
        cassandra_exe_path: Path to the compiled Cassandra executable.
            After ``make`` in ``Src/``, this is usually
            ``Src/cassandra_gfortran.exe``.
        raise_on_error: If ``True`` (default), raise :class:`CassandraRunError`
            when the executable is missing or exits with a non-zero status.
        verbose: If ``True`` (default), print the working directory, command,
            and captured stdout/stderr.

    Returns:
        The :class:`subprocess.CompletedProcess` from the run. Check
        ``result.returncode`` when ``raise_on_error=False``.

    Raises:
        FileNotFoundError: If the input file or executable does not exist.
        CassandraRunError: If ``raise_on_error`` is ``True`` and Cassandra
            exits with a non-zero status.

    Output files:
        Written in the same directory as the ``.inp`` file. Names are derived
        from ``Run_Name`` in the input file, with extensions such as
        ``.log``, ``.prp``, ``.xyz``, and ``.chk``.
    """
    inp_file_abs = os.path.abspath(input_file_path)
    exe_path_abs = os.path.abspath(cassandra_exe_path)

    if not os.path.isfile(inp_file_abs):
        raise FileNotFoundError(f"Input file not found: {inp_file_abs}")
    if not os.path.isfile(exe_path_abs):
        raise FileNotFoundError(
            f"Cassandra executable not found: {exe_path_abs}\n"
            "Compile Cassandra first by running 'make' in the 'Src' directory."
        )

    run_directory = os.path.dirname(inp_file_abs)
    command = [exe_path_abs, inp_file_abs]

    if verbose:
        print(f"--- Running Cassandra for {input_file_path} ---")
        print(f"Running in directory: {run_directory}")
        print(f"Running command: {command}")

    try:
        result = subprocess.run(
            command,
            cwd=run_directory,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as exc:
        raise CassandraRunError(
            f"Failed to start Cassandra executable: {exe_path_abs}"
        ) from exc

    if verbose:
        print("--- CASSANDRA STDOUT ---")
        print(result.stdout)
        if result.stderr:
            print("--- CASSANDRA STDERR ---")
            print(result.stderr)
        print("--- Cassandra run complete ---")

    if raise_on_error and result.returncode != 0:
        raise CassandraRunError(
            f"Cassandra exited with status {result.returncode}.\n"
            f"Command: {command}\n"
            f"Working directory: {run_directory}\n"
            f"stderr:\n{result.stderr or '(empty)'}"
        )

    return result
