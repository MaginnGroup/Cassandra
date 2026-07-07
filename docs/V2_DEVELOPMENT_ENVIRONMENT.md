# Cassandra V2 Development Environment

A practical guide to conda environments, compilers, and the day-to-day workflow
for building and running Cassandra on your Mac.

**Related docs:** [docs index](README.md) · [Git workflow](V2_GIT_WORKFLOW.md) ·
[Python API](../python/README.md)

---

## Returning after time away? Start here.

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
which python gfortran
python run_test.py
```

If all three commands succeed, your environment is ready.

---

## What is a conda environment?

Conda environments are isolated sets of installed packages. Each project (or class)
can have its own Python version, compilers, and libraries without conflicts.

Your terminal prompt shows which environment is active:

| Prompt | Meaning |
|--------|---------|
| `(thermo2026) edwardmaginn@... %` | In your **thermodynamics class** environment |
| `(cassandra-dev) edwardmaginn@... %` | In the **Cassandra development** environment |
| `edwardmaginn@... %` (no prefix) | In `base` (default conda) — avoid doing project work here |

**Rule:** activate the right environment *before* compiling or running code.

---

## Essential conda commands

```bash
# List all environments (* marks the active one)
conda env list

# Switch to an environment
conda activate cassandra-dev

# Leave the current environment
conda deactivate

# See what is installed in the active environment
conda list

# Install a package (while environment is active)
conda install -c conda-forge package-name
```

### Creating a new environment (rarely needed)

```bash
conda create -n my-env-name python=3.10 -y
conda activate my-env-name
conda install -c conda-forge some-package -y
```

You already have `cassandra-dev` set up — you usually just `conda activate` it.

---

## Recommended environment: `cassandra-dev`

Use **`cassandra-dev`** for all Cassandra V2 work on this machine.

| Package | Purpose |
|---------|---------|
| Python 3.10 | Python wrapper (`python/`), test scripts |
| gfortran | Compile `Src/cassandra_gfortran.exe` |
| cmake | Build bundled libraries under `Libraries/` |
| numpy, pandas | Test suite and analysis scripts |

### Do not use `thermo2026` for Cassandra

`thermo2026` is for your thermodynamics class. It has a different Python version
(3.11) and lacks the Fortran compiler Cassandra needs. Switching environments is
normal — you are not doing anything wrong by having both.

```bash
# Class homework
conda activate thermo2026

# Cassandra modernization
conda activate cassandra-dev
```

---

## Verify your environment

Run these after `conda activate cassandra-dev`:

```bash
which python
# expect: /opt/anaconda3/envs/cassandra-dev/bin/python

which gfortran
# expect: /opt/anaconda3/envs/cassandra-dev/bin/gfortran

python --version
# expect: Python 3.10.x

gfortran --version
# expect: GNU Fortran version ...
```

If `gfortran` is not found, you are in the wrong environment.

---

## Daily development workflow

### 1. Activate and go to the repo

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
```

### 2. Compile (after Fortran or library changes)

```bash
cd Src
make -f Makefile.gfortran
```

This produces `Src/cassandra_gfortran.exe`.

A clean rebuild:

```bash
make -f Makefile.gfortran clean
make -f Makefile.gfortran
```

### 3. Run the Python smoke test

```bash
cd ~/CassandraV2/Cassandra
python run_test.py
```

This calls `python/api.py` to run the SPC water NVT example at
`Examples/NVT/water_spc/nvt.inp`. See [smoke-test.md](../docs/smoke-test.md) for
input/output file locations and what to expect.

### 4. Run the legacy test suite (optional, slower)

```bash
cd Scripts/testSuite
python testSuite.py ../../Src/cassandra_gfortran.exe
```

---

## Other environments on this machine

You may see these in `conda env list`:

| Environment | Use for |
|-------------|---------|
| `cassandra-dev` | **Cassandra V2 modernization** (recommended) |
| `thermo2026` | Thermodynamics class |
| `mosdef_cassandra125` | MoSDeF-Cassandra interface (separate project) |
| `cassandra`, `cassandra-tutorial` | Older Cassandra-related setups |

When in doubt, use `cassandra-dev`.

---

## Recreate `cassandra-dev` from scratch

Only if the environment is corrupted or you are setting up a new machine:

```bash
conda create -n cassandra-dev python=3.10 -y
conda activate cassandra-dev
conda install -c conda-forge cmake gfortran numpy pandas -y
```

Then recompile:

```bash
cd ~/CassandraV2/Cassandra/Src
make -f Makefile.gfortran clean
make -f Makefile.gfortran
```

---

## Common problems

### `gfortran: command not found`

You are not in `cassandra-dev`. Run:

```bash
conda activate cassandra-dev
```

### `make` fails after a long absence

Try a clean rebuild:

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra/Src
make -f Makefile.gfortran clean
make -f Makefile.gfortran
```

### `python run_test.py` — executable not found

Compile first (see above). The smoke test expects
`Src/cassandra_gfortran.exe`.

### Wrong Python packages / import errors

Confirm which Python is running:

```bash
which python
conda list numpy
```

If it points outside `cassandra-dev`, activate the correct environment.

### Accidentally ran smoke test in wrong directory

`run_test.py` must be run from the repo root (`~/CassandraV2/Cassandra`), not
from `Src/`.

---

## Quick reference card

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
cd Src && make -f Makefile.gfortran && cd ..
python run_test.py
```

---

## Related docs

- [docs/README.md](README.md) — index of all V2 developer reminders
- [V2_GIT_WORKFLOW.md](V2_GIT_WORKFLOW.md) — git branch, commit, push, SSH
- [python/README.md](../python/README.md) — `run_cassandra()` API
