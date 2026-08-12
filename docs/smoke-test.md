# Smoke Test (`run_test.py`)

A quick reminder for the Cassandra V2 smoke test: how to run it, where files live,
and what a successful run looks like.

**Related docs:** [docs index](README.md) ·
[Development environment](V2_DEVELOPMENT_ENVIRONMENT.md) ·
[Python API](../python/README.md)

---

## What this test does

`run_test.py` is a minimal end-to-end check that the V2 Python wrapper can:

1. Find the compiled Cassandra executable
2. Launch it via `python/api.py` → `run_cassandra()`
3. Run a short **NVT Monte Carlo** simulation of **SPC water**
4. Write output files in the correct location

It does **not** parse results or compare against reference values — it only
confirms the toolchain works.

---

## How to run

From the **repo root** (not `Src/`):

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
python run_test.py
```

### Prerequisites

| Requirement | Location |
|-------------|----------|
| Conda environment | `cassandra-dev` |
| Compiled executable | `Src/cassandra_gfortran.exe` |
| Input file | `Examples/NVT/water_spc/nvt.inp` |

If the executable is missing, compile first:

```bash
cd Src
make -f Makefile.gfortran
cd ..
```

---

## What happens under the hood

```
run_test.py
    └── api.run_cassandra(
            input_file_path = Examples/NVT/water_spc/nvt.inp
            cassandra_exe_path = Src/cassandra_gfortran.exe
        )
            └── subprocess: cassandra_gfortran.exe /full/path/to/nvt.inp
                └── cwd = Examples/NVT/water_spc/
```

Cassandra is run **from the input file's directory** so it can find supporting
files (`.mcf`, `.ff`, fragments) and write outputs alongside them.

---

## Input files

### Main input

```
~/CassandraV2/Cassandra/Examples/NVT/water_spc/nvt.inp
```

Key settings in the current smoke test input:

| Setting | Value | Meaning |
|---------|-------|---------|
| `Run_Name` | `nvt.out` | Base name for all output files |
| `Sim_Type` | `nvt_mc` | NVT (canonical) Monte Carlo |
| `Molecule_Files` | `spc.mcf 90` | 90 SPC water molecules |
| `Temperature_Info` | `300.0` | 300 K |
| `Start_Type` | `checkpoint nvt.inp.chk` | Restart from checkpoint |
| `Run_Type` | `equilibration 900` | Equilibration mode |
| `run` | `100` | 100 sweeps (short test) |

### Supporting files (same directory)

```
Examples/NVT/water_spc/
├── nvt.inp              ← main input (what run_test.py passes)
├── nvt.inp.chk          ← checkpoint used by Start_Type
├── spc.mcf               ← molecule connectivity
├── spc.ff                ← force field parameters
├── spc.pdb               ← molecular structure
└── species1/fragments/   ← fragment data for CBMC
```

Cassandra reads paths in `nvt.inp` relative to `Examples/NVT/water_spc/`.

---

## Output files

Outputs are written to the **same directory as the input file**:

```
~/CassandraV2/Cassandra/Examples/NVT/water_spc/
```

The `Run_Name` field in `nvt.inp` controls the filename prefix. With
`Run_Name = nvt.out`, expect:

| File | Purpose |
|------|---------|
| `nvt.out.log` | Simulation log (setup, MCF echo, progress snapshots, final stats) |
| `nvt.out.prp` | Thermodynamic properties (energy, pressure, etc.) |
| `nvt.out.xyz` | Coordinate trajectory |
| `nvt.out.chk` | Updated checkpoint for restarts |
| `nvt.out.H` | Energy histogram |

List outputs after a run, newest first:

```bash
ls -lt ~/CassandraV2/Cassandra/Examples/NVT/water_spc/nvt.out.*
```

### Older outputs you may see

Previous test runs with different `Run_Name` values may have left files like:

- `water_spc_nvt.out.*`
- `water_spc_nvt-2.out.*`

Those are **not** from the current `nvt.inp` unless you change `Run_Name` back.
They are safe to ignore or delete locally — do not commit them to git.

---

## What success looks like

### Terminal output

You should see something like:

```
Found example input: .../Examples/NVT/water_spc/nvt.inp
Using Cassandra executable: .../Src/cassandra_gfortran.exe

=== CALLING api.run_cassandra() ===
--- Running Cassandra for .../nvt.inp ---
Running in directory: .../Examples/NVT/water_spc
Running command: ['.../cassandra_gfortran.exe', '.../nvt.inp']
--- CASSANDRA STDOUT ---
 Begin Cassandra simulation
 ...
 Cassandra simulation complete
--- Cassandra run complete ---
=== TEST SCRIPT FINISHED ===
```

### Log file confirmation

Open `Examples/NVT/water_spc/nvt.out.log` and check the end:

```
************************ Cassandra simulation complete *************************
```

The log also records total wall time (typically a few seconds for this short run).

### Property output

`nvt.out.prp` contains columns such as `Energy_Total` and `Pressure` sampled
during the run. Values will vary slightly between runs due to Monte Carlo
statistics, but the file should exist and contain numeric data.

---

## What failure looks like

| Symptom | Likely cause | Fix |
|---------|--------------|-----|
| `Cassandra executable not found` | Not compiled | `make -f Makefile.gfortran` in `Src/` |
| `Input file not found` | Wrong working directory | Run from repo root, not `Src/` |
| `Cassandra exited with status 1` | Simulation error | Read `nvt.out.log` or stderr in terminal |
| `gfortran: command not found` | Wrong conda env | `conda activate cassandra-dev` |
| `Could not find the 'api.py' module` | Missing `python/` folder | Check repo is intact |

On failure, `run_cassandra()` raises `CassandraRunError` with the exit code and
stderr — read that message first, then check `nvt.out.log`.

---

## Git reminder

Running the smoke test **modifies output files** under `Examples/NVT/water_spc/`.
Before committing:

```bash
git status
```

Stage source changes only — not simulation outputs:

```bash
git add python/ docs/ run_test.py    # good
git add .                              # risky — may pick up .out.* files
```

---

## Quick reference

```bash
# Run
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
python run_test.py

# Check outputs
ls -lt Examples/NVT/water_spc/nvt.out.*
tail -5 Examples/NVT/water_spc/nvt.out.log
```

**Input:** `Examples/NVT/water_spc/nvt.inp`  
**Outputs:** `Examples/NVT/water_spc/nvt.out.*`  
**Executable:** `Src/cassandra_gfortran.exe`
