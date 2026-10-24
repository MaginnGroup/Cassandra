# Output file formats (`.xyz`, `.prp`, and `.log`)

Brief notes for V2 modernization: how Cassandra names its output files, how
it writes trajectory, property, and logfile outputs, and what future analysis
tools should assume.

**Related:** [docs index](README.md) · [smoke test](smoke-test.md) ·
[MCF setup](mcf-setup-workflow.md) · [logfile format](logfile-format.md) ·
[NVT restart utility](../utilities/README.md)

---

## File names

`# Run_Name` is the prefix for every file Cassandra creates. Call that
prefix `[run_name]`. One example is `R125-T243-NVT`; the names below use
`[run_name]` so the same pattern applies to any run. The files are written
in the directory where you launched the job. One box produces:

| File | Written | What a killed job leaves |
|------|---------|--------------------------|
| `[run_name].log` | Throughout the run | The log up to the last flush |
| `[run_name].prp` | Every `prop_freq` | Properties through the last sample |
| `[run_name].chk` | Every `coord_freq` | A checkpoint only if a dump finished. There is no checkpoint at the end of a run unless that step lands on `coord_freq` |
| `[run_name].xyz` | Every `coord_freq`, appended | The coordinate movie. This file stays open, so a kill can leave the last frame short |
| `[run_name].H` | Every `coord_freq`, appended | Volume, cell matrix, and molecule counts for the movie. Same risk as the `.xyz` |
| `[run_name].restart.xyz` | Every `coord_freq`, replaced | One complete frame at `F18.8`, if that dump finished. Cassandra closes this file before returning |
| `[run_name].restart.H` | Every `coord_freq`, replaced | The cell and molecule counts for that frame. Closed with the restart `.xyz` |

Two or more boxes insert `.box1`, `.box2`, and so on before the coordinate
suffix: `[run_name].box1.xyz`, `[run_name].box1.restart.H`. A property file for
box 2 is `[run_name].box2.prp`. Code: `Name_Files` in `Src/io_utilities.f90`
and the assignments in `Get_Simulation_Length_Info` (`Src/input_routines.f90`).

`read_config` reads exactly the number of atoms implied by the molecule
counts. It does not stop at the atom count on the first XYZ line. If the
chosen file ends early, the run aborts with an end-of-file error while
reading coordinates. When that file's name ends in `.xyz`, Cassandra also
looks for the companion formed by replacing `.xyz` with `.H` and checks the
cell and molecule counts against `# Box_Info`.

### Naming practice

Give the starting configuration a name the run will not open for writing.

```text
# Run_Name
[run_name]

# Start_Type
read_config 800 [run_name].start.xyz
```

Keep `[run_name].start.H` beside that XYZ. Cassandra will check it
against `# Box_Info`.

Leave these names for Cassandra:

- `[run_name].xyz` and `[run_name].H` are the movie. Pointing
  `read_config` at the movie stores the input and the output in one file.
  The first coordinate dump opens it, and deleting the job can leave it
  short. The next submission then reads the damaged file.
- `[run_name].restart.xyz` and `[run_name].restart.H` are the latest
  finished frame. They are the right source for a later run. Copy them to
  a new `.start` pair, or start a new `# Run_Name`, before submitting.
  The next dump replaces the `.restart` files, so a kill during that
  replacement can damage the file the job is also reading.

The multi-frame movie is also the wrong start for another reason: `read_config`
uses its first frame, which is the beginning of the trajectory rather than
the latest configuration.

---

## Precision (updated 2026-08)

| File | Format | Notes |
|------|--------|--------|
| `.xyz` coordinates | `F12.3` (Å) | Three decimal places (~0.001 Å) |
| `.prp` properties | `E16.6` | Scientific notation; width matches `A16` headers |

Previously, `.xyz` used list-directed output (~16 digits) and `.prp` used
`E16.8`. Headers for the step column and property names/units are
**right-aligned** so they line up with `I19` / `E16.6` data. Every listed
property should carry a units string (including `Enthalpy` and `Nmols`).

The new formats shrink files without affecting typical analysis or
float-based unit-test comparisons.

Code: `Src/write_properties.f90` (`Write_Coords_XYZ`, `Write_Header`,
`Write_Properties_Buffer`).

---

## `.xyz` movie file

Standard XYZ multi-frame layout:

1. Number of atoms
2. Comment line: `MC_STEP: <step>`
3. One line per atom: `element  x  y  z` (Å)

**Atom / molecule order** (important for analysis):

- Species loop outermost, then molecule index, then atoms in **MCF order**
- Atom \(k\) of molecule \(m\) of species \(s\) matches MCF atom \(k\) for that species
- For NVT/NPT with fixed \(N\), molecule index is stable across frames
- For GCMC (changing \(N\)), use the companion `.H` file for nmols per frame

**Topology is not in the XYZ.** Bonds, angles, \(K_\theta\), charges, etc. come
from the MCF (and temperature from the `.inp` / `.log`).

Companion **`.H`** file (same `coord_freq`): volume, H-matrix, species counts —
needed for box size and GCMC analysis.

---

## Restart files (`*.restart.xyz` / `*.restart.H`)

Written on every coordinate dump and again at a successful end of run
(**overwrite**, single frame — not a movie):

| File | Contents |
|------|----------|
| `run_name.restart.xyz` | One XYZ frame (`F18.8` coords); for `# Start_Type read_config` |
| `run_name.restart.H` | One `.H` block: volume, cell matrix, nspecies, `(is, nmols)` |

Multi-box: `run_name.boxN.restart.xyz` / `.H`.

Pass `*.restart.xyz` to `read_config` only after copying it to a name the
next run will not overwrite. The movie `.xyz` is the wrong start: `read_config`
uses its first frame, and that file is also the movie Cassandra keeps open.
Set `# Box_Info` and the molecule counts from the matching `.restart.H`.
See [File names](#file-names) above and [V2_CRC_WORKFLOW.md](V2_CRC_WORKFLOW.md) §5.

---

## `.log` simulation logfile

Human-readable archive of the run (input echo, MCF echo, energies, 10% progress
snapshots, final acceptance and timings). Mid-run Step/Move/Success spam is
omitted by default. Full details: [logfile-format.md](logfile-format.md).

---

## `.prp` property file

Header (comment lines):

1. Instantaneous vs block averages
2. Property names (e.g. `Energy_Total`, `Pressure`, `Volume`)
3. Units

Data columns:

1. MC step or sweep (per `# Simulation_Length_Info` units)
2. One column per requested `# Property_Info` entry

Parse with whitespace split and `float()`; do not rely on fixed column
character positions. Header names identify columns for plotting tools.

Typical units (see also the user guide):

| Property | Units |
|----------|--------|
| Energies, enthalpy | kJ/mol (extensive) |
| Pressure | bar |
| Volume | Å³ |
| Density | molec/Å³ |
| Mass density | kg/m³ |
| Nmols | molecules |

---

## Interactive stdout (updated 2026-08)

When you run Cassandra in the foreground, the terminal shows:

1. A short **setup** message while the initial configuration is built  
2. A **run summary** (ensemble, T, N, box, write frequencies, OpenMP on/off)  
3. **Property lines** at the same frequency as `.prp` (`prop_freq`), mirroring
   `# Property_Info` with the same `E16.6` formatting  
4. `Cassandra simulation complete` at the end  

The cryptic `openmp_flag = F` message has been removed; OpenMP status appears
in the run summary as `enabled` or `disabled`.

Screen output is for interactive confidence only. Analysis tools should still
read `.prp` and `.xyz` files.

Code: `Write_Stdout_Run_Banner` and stdout echo in `Write_Properties`
(`Src/write_properties.f90`, `Src/main.f90`).

---

## Implications for future Python tools

Diagnostic tools live in [`diagnostics/`](../diagnostics/README.md) (not
`Scripts/` or `python/`). Shared parsers: `diagnostics/io/prp.py`,
`diagnostics/io/xyz.py`, `diagnostics/io/mcf.py`, `diagnostics/io/inp.py`,
`diagnostics/io/layout.py`.

| Tool idea | Primary inputs | Status |
|-----------|----------------|--------|
| Property vs step (instantaneous + running average) | `.prp` | `property_plotter.py` (**diagnostic**) |
| Fixed bond length vs MCF \(r_0\) | `.xyz` + MCF | `bond_distribution.py` (**diagnostic**) |
| Bond-angle PDF + ideal-gas Boltzmann | `.xyz` + MCF; \(T\) / species from `.inp` | `angle_distribution.py` (**diagnostic**) |
| Dihedral PDF + bare-potential Boltzmann | `.xyz` + MCF; \(T\) / species from `.inp` | `dihedral_distribution.py` (**diagnostic**) |
| Multi-species fixed-\(N\) layout | `# Molecule_Files` (+ `--nmols`) | supported in bond/angle/dihedral tools |
| COM MSD vs sweep (NVT/NPT exploration; not MD \(D\)) | `.xyz` + MCF + `.inp` (+ `.H` for NPT) | `msd_com.py` (**diagnostic**) |
| Block averages | `.prp` | `utilities/block_average.py` |
| \(g(r)\) and coordination number | `.xyz` + `.H` + MCF | `utilities/rdf.py` |
| Variable-\(N\) (GCMC) from `.H` | `.xyz` + `.H` + MCF | planned |

Keep XYZ as plain XYZ (VMD-compatible). Do not require molecule IDs in the
coordinate file for Phase 1 analysis designs.
