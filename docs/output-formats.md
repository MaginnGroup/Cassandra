# Output file formats (`.xyz`, `.prp`, and `.log`)

Brief notes for V2 modernization: how Cassandra writes trajectory, property,
and logfile outputs, and what future analysis tools should assume.

**Related:** [docs index](README.md) · [smoke test](smoke-test.md) ·
[MCF setup](mcf-setup-workflow.md) · [logfile format](logfile-format.md)

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
| Block averages, \(g(r)\), … | `.prp` / `.xyz` | planned **post-processing** (separate) |
| Variable-\(N\) (GCMC) from `.H` | `.xyz` + `.H` + MCF | planned |

Keep XYZ as plain XYZ (VMD-compatible). Do not require molecule IDs in the
coordinate file for Phase 1 analysis designs.
