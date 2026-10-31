# Cassandra V2 Development Environment

A practical guide to the laptop Pixi environment, the compiler, and the
day-to-day workflow for building and running Cassandra on your Mac.

The development branch is `pixi-env`. It continues `modernization` and replaces
the `cassandra-dev` Conda environment with the Pixi environment in
[`pixi.toml`](../pixi.toml). CRC compilation is separate: on maginnfe, use the
`gcc` module, not this Pixi environment. See
[V2_CRC_WORKFLOW.md](V2_CRC_WORKFLOW.md).

**Related docs:** [docs index](README.md) · [Git workflow](V2_GIT_WORKFLOW.md) ·
[CRC / maginnfe](V2_CRC_WORKFLOW.md) ·
[Python API](../python/README.md) · [MCF setup workflow](mcf-setup-workflow.md)

---

## Returning after time away? Start here.

From the repo root, outside any activated Conda environment:

```bash
cd ~/CassandraV2/Cassandra
git branch --show-current          # expect: pixi-env
pixi install                       # first time, or after pixi.toml changes
pixi run smoke
```

`pixi run smoke` compiles `Src/cassandra_gfortran.exe` if needed, then runs
`run_test.py`. If that finishes, the laptop environment is ready.

Do not `conda activate cassandra-dev` in the same shell. An active Conda
environment exports `FC`, `LDFLAGS`, and `CMAKE_ARGS`, and those override
Pixi's compiler.

---

## What Pixi is doing

Pixi installs the packages named in `pixi.toml` into `.pixi/` inside this
repo. `.pixi/` is gitignored. `pixi.lock` records the exact versions. The
laptop platform is `osx-arm64`.

The pins match the old `cassandra-dev` Conda environment as of 2026-10-30:
Python 3.10.19, gfortran 15.1.0, CMake, NumPy, pandas, matplotlib, and RDKit.
`MACOSX_DEPLOYMENT_TARGET=26.0` is set whenever a Pixi task or `pixi shell`
runs. That avoids `invalid version number in '-mmacosx-version-min=18.0'`
on macOS 27.

| Package | Purpose |
|---------|---------|
| Python 3.10 | Python wrapper (`python/`), tests, analysis scripts |
| gfortran 15.1 | Compile `Src/cassandra_gfortran.exe` |
| cmake | Build bundled libraries under `Libraries/` |
| numpy, pandas, matplotlib | Test suite and analysis scripts |
| RDKit | `Scripts/molecule_to_pdb` |

`thermo2026` remains the thermodynamics class environment. Leave it for class
work. It is not used for Cassandra.

---

## Essential Pixi commands

Run these from `~/CassandraV2/Cassandra`:

```bash
pixi install                              # create or update .pixi/
pixi shell                                # interactive shell with Python and gfortran
pixi run build                            # compile Src/cassandra_gfortran.exe
pixi run clean                            # make clean
pixi run smoke                            # build if needed, then python run_test.py
pixi run test                             # legacy regression suite; builds first
pixi run test-utilities                   # block-average, restart, and RDF tests
pixi run python utilities/rdf.py --help   # any command, inside the environment
pixi add package-name                     # add a conda-forge package to pixi.toml
```

`pixi shell` changes the prompt. Leave it with `exit`. One-shot commands are
clearer as `pixi run ...`, because they do not depend on the shell already
being activated.

---

## Verify the environment

```bash
cd ~/CassandraV2/Cassandra
pixi run which python
pixi run which gfortran
pixi run python --version          # Python 3.10.19
pixi run gfortran --version        # GNU Fortran (GCC) 15.1.0
pixi run printenv MACOSX_DEPLOYMENT_TARGET   # 26.0
```

The Python and gfortran paths should be under
`~/CassandraV2/Cassandra/.pixi/envs/default/bin/`.

Homebrew's GCC stays installed because Homebrew Open MPI depends on it, but it
is unlinked so `/opt/homebrew/bin/gfortran` does not override Pixi's compiler.
Do not `brew uninstall gcc` while that dependency remains.

---

## Daily development workflow

### 1. Go to the repo

```bash
cd ~/CassandraV2/Cassandra
git branch --show-current          # expect: pixi-env
```

### 2. Compile

```bash
pixi run build
```

This produces `Src/cassandra_gfortran.exe`. The first build also compiles
libxdrfile and libgmxfort into `Libraries/locals`.

A clean rebuild:

```bash
pixi run clean
pixi run build
```

### 3. Smoke test

```bash
pixi run smoke
```

This calls `python/api.py` to run the SPC water NVT example at
`Examples/NVT/water_spc/nvt.inp`. See [smoke-test.md](smoke-test.md).

### 4. Utility tests and the legacy suite

```bash
pixi run test-utilities
pixi run test
```

`pixi run test` is the older regression suite. It is slower than the smoke test.

---

## Common problems

### `gfortran: command not found`

The command ran outside Pixi. Use `pixi run gfortran` or `pixi shell`, from
the repo root, after `pixi install`.

### Wrong compiler or deployment target

`pixi.toml` sets `MACOSX_DEPLOYMENT_TARGET=26.0`. If a Conda environment is
also active, leave it (`conda deactivate`) and run the Pixi command again.

### `make` fails after a long absence

```bash
pixi run clean
pixi run build
```

### `python run_test.py` — executable not found

`pixi run smoke` compiles first. A bare `python run_test.py` uses whatever
Python is on your path and expects `Src/cassandra_gfortran.exe` to exist already.

### Wrong Python packages

```bash
pixi run which python
pixi run python -c "import numpy, matplotlib, rdkit"
```

If the path is not under `.pixi/envs/default/`, the shell is not using this
environment.

### Smoke test started in the wrong directory

`run_test.py` is run from the repo root. `pixi run smoke` does that. Do not
launch it from `Src/`.

---

## Quick reference card

```bash
cd ~/CassandraV2/Cassandra
git branch --show-current          # pixi-env
pixi run build
pixi run smoke
pixi run test-utilities
```

---

## Related docs

- [docs/README.md](README.md) — index of all V2 developer reminders
- [V2_GIT_WORKFLOW.md](V2_GIT_WORKFLOW.md) — git branch, commit, push, SSH
- [mcf-setup-workflow.md](mcf-setup-workflow.md) — new molecule PDB → MCF → inp checklist
- [python/README.md](../python/README.md) — `run_cassandra()` API
