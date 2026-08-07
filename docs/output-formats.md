# Output file formats (`.xyz` and `.prp`)

Brief notes for V2 modernization: how Cassandra writes trajectory and property
files, and what future analysis tools should assume.

**Related:** [docs index](README.md) · [smoke test](smoke-test.md) ·
[MCF setup](mcf-setup-workflow.md)

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
`Scripts/` or `python/`). Shared `.prp` parsing is in
`diagnostics/io/prp.py`; the first CLI is `diagnostics/property_plotter.py`.

| Tool idea | Primary inputs | Status |
|-----------|----------------|--------|
| Property vs step (instantaneous + running average) | `.prp` | `property_plotter.py` |
| Bond-angle PDF (one molecule or all) | `.xyz` + MCF (+ `.H` if \(N\) changes) | planned |
| Compare to ideal-gas Boltzmann angle distribution | MCF \(K_\theta\), \(\theta_0\); \(T\) from `.inp`/`.log` | planned |

Keep XYZ as plain XYZ (VMD-compatible). Do not require molecule IDs in the
coordinate file for Phase 1 analysis designs.
