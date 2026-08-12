# Simulation logfile format (`.out.log`)

What Cassandra writes to `{run_name}.log` (Fortran unit `logunit` = 25).
Updated for the V2 logfile redesign (2026-08).

**Related:** [docs index](README.md) · [output formats](output-formats.md) ·
[smoke test](smoke-test.md)

---

## Purpose

The logfile is the **human-readable archive** of a run: what was asked for,
which force field was used, whether the run looks healthy, and how it finished.
It is **not** the place for high-frequency move diagnostics (those bloated older
logs). Instantaneous properties belong in `.prp`; coordinates in `.xyz` / `.H`.

---

## Structure (top to bottom)

| Section | Contents |
|---------|----------|
| Banner | ASCII art + one-line subtitle (`Cassandra V2 — simulation logfile`) |
| Citation | Journal reference |
| **Run info** | Aligned fields: version, input, run name, date/time, machine, verbosity, simulation type (underlined heading; no `***` box) |
| **Copy input file** | Input echoed (comments starting with `!` stripped) |
| Setup summary | Species count, box, pair style, molecule files, max molecules, … |
| **Copy MCF files** | Each species MCF echoed (force-field / topology record) |
| Remaining setup | Frequencies, moves, cutoffs, etc. through initialization |
| Begin Cassandra simulation | Initial configuration energy table |
| **Run simulation** | Progress snapshots only (see below) |
| Final configuration | Full energy table + acceptance + fragment stats + subroutine times |
| Cassandra simulation complete | Wall-clock execution time breakdown |

### Preamble style (V2)

Ordinary early sections use a short underlined title instead of an 80-character
`***` fence above and below:

```text
 Run info
 --------
  version:      Cassandra Version 1.3.0
  inputfile:    npt.inp
  ...
  verbosity:    normal
  simulation:   NPT_MC
```

Heavy `****` rules are reserved for the ASCII banner and major milestones
(Begin / Progress / Complete). Everything else — setup (temperature through
run type), initial configuration, charge check, energy headers, progress
sub-blocks (acceptance / widths / timings), and program execution time —
uses underlined headings and aligned `label: value` lines where practical.

---

## Energy tables (extensive + intensive)

Full component tables appear at **initial** configuration, **final**
configuration, and whenever `# Energy_Check` / `echeck` fires.

| Column | Meaning |
|--------|---------|
| Extensive | Total for the box (kJ/mol) |
| Intensive | Extensive ÷ \(N_{\mathrm{mol}}\) in that box (kJ/mol per molecule) |

If the box has zero molecules, intensive is printed as `n/a`. Multi-box
ensembles (GEMC, etc.) print one table (or progress line) **per box**.

Code: `Check_System_Energy` in `Src/energy_routines.f90`.

---

## Mid-run logging

### Removed by default

The old stream of lines like

```text
Step       Move   Mol Spc Box  Success  MaxWidth
 185  translate         1   1  0.18000   0.72000
```

is **no longer written**. Equilibration still **updates** translation /
rotation / volume widths in memory; current widths appear in progress
snapshots instead.

### Kept for debugging

`# Verbose_Logfile TRUE` still writes per-move accept/reject detail.

### Progress snapshots (10% … 90%)

For step-based runs with at least 10 remaining steps
(`n_mcsteps - initial_mcstep >= 10`), the drivers write a block when the
step counter first reaches each tenth of the remaining span:

```text
 Progress  40%   step     4000 /    10000
 Box 1: N = …   V = … Ang^3   density = … kg/m^3
 Total system energy:  … kJ/mol   (intensive … kJ/mol/molecule)
 Move acceptance (cumulative)   ← same tables as end-of-run
 Current move widths            ← max_disp, max_rot, dv_max
 Subroutine times so far
```

Timed runs (`Minutes`) skip percentage snapshots for now.

Code: `Init_Log_Progress`, `Maybe_Write_Log_Progress`, `Write_Log_Progress`
in `Src/read_write_checkpoint.f90`; called from NVT / NPT / GCMC / GEMC drivers.

---

## MCF echo

After molecule files are known, `Copy_MCF_Files` appends each species MCF to
the log under `Copy MCF files`. Use this when checking that a cluster run
really used the expected force-field parameters.

---

## Where to look when something fails

1. Stderr / terminal message from `Clean_Abort`
2. End of `.log` for the last completed section
3. `.prp` for thermodynamic drift (use `diagnostics/property_plotter.py`)
