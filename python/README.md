# Cassandra Python Wrapper

Thin Python layer for launching the Cassandra Fortran executable. This is the
starting point for Cassandra V2 modernization — a clean subprocess API separate
from the legacy scripts in `Scripts/testSuite/`.

For git branch strategy, commit/push workflow, and authentication troubleshooting,
see [`docs/V2_GIT_WORKFLOW.md`](../docs/V2_GIT_WORKFLOW.md).

For conda environments, compiling, and the day-to-day dev workflow, see
[`docs/V2_DEVELOPMENT_ENVIRONMENT.md`](../docs/V2_DEVELOPMENT_ENVIRONMENT.md).

**All V2 reminders:** [`docs/README.md`](../docs/README.md)

## How it works

Cassandra is **not** imported as a Python library. `run_cassandra()` spawns the
compiled executable with one argument: the path to a `.inp` file.

```
your_script.py  →  run_cassandra()  →  cassandra_gfortran.exe  →  nvt.inp
```

The Fortran driver (`Src/main.f90`) requires:

```bash
cassandra.exe path/to/run.inp
```

Output files (`.log`, `.prp`, `.xyz`, `.chk`, etc.) are written relative to the
**current working directory**, using the `Run_Name` field from the input file.
`run_cassandra()` therefore sets `cwd` to the directory containing the `.inp`
file so outputs appear next to the example data.

## Quick start

### 1. Compile Cassandra

```bash
cd Src
make -f Makefile.gfortran
```

This produces `Src/cassandra_gfortran.exe` (name may vary by compiler).

### 2. Run the smoke test

From the repository root:

```bash
python run_test.py
```

This runs the SPC water NVT example at `Examples/NVT/water_spc/nvt.inp`.

### 3. Use from your own script

```python
import sys
import os

# Add the python/ directory to the import path
sys.path.append(os.path.join(os.path.dirname(__file__), "python"))

from api import run_cassandra, CassandraRunError

try:
    result = run_cassandra(
        input_file_path="Examples/NVT/water_spc/nvt.inp",
        cassandra_exe_path="Src/cassandra_gfortran.exe",
    )
    print(f"Finished with exit code {result.returncode}")
except FileNotFoundError as exc:
    print(f"Missing file: {exc}")
except CassandraRunError as exc:
    print(f"Simulation failed: {exc}")
```

## API reference

### `run_cassandra(input_file_path, cassandra_exe_path="cassandra.exe", *, raise_on_error=True, verbose=True)`

| Parameter | Description |
|-----------|-------------|
| `input_file_path` | Path to the Cassandra `.inp` file |
| `cassandra_exe_path` | Path to the compiled executable |
| `raise_on_error` | Raise `CassandraRunError` on non-zero exit (default: `True`) |
| `verbose` | Print command and captured output (default: `True`) |

**Returns:** `subprocess.CompletedProcess` with `.returncode`, `.stdout`, `.stderr`

**Raises:**

- `FileNotFoundError` — input file or executable missing
- `CassandraRunError` — executable failed to start or exited non-zero (when `raise_on_error=True`)

## Relationship to other Python code

| Location | Purpose |
|----------|---------|
| `python/api.py` | New V2 wrapper (this module) |
| `diagnostics/` | Post-run analysis tools (property plots, …) |
| `run_test.py` | Smoke test that calls `run_cassandra()` |
| `Scripts/testSuite/` | Legacy regression tests with their own `subprocess` calls |
| [MoSDeF-Cassandra](https://mosdef-cassandra.readthedocs.io/) | Separate higher-level Python interface |

The test suite does **not** use this API yet. Consolidating on `run_cassandra()`
is a future modernization step.
