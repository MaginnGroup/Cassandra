# Cassandra utilities

Tools that turn a finished simulation into a number or a new input file.
Diagnostics (`diagnostics/`) answer “does this run look healthy?” These
utilities answer “what should I report, and how do I start the next run?”

A first-year student can read the header at the top of each script. That
header is the explanation; this page is only the command reminder.

## Block averages

`block_average.py` reads a Cassandra `.prp` file.

- It chooses where the production average starts. The default series for that
  choice is `Mass_Density`. A drifting beginning, such as the early part of a
  warm NPT run, is left out. `--discard-sweeps` overrides the choice.
- It reports the mean of every column.  The error bar can come from either
  of two methods, selected with `--error-method`:
  - `autocorr` estimates a correlation time (Geyer, Statistical Science
    1992, 7, 473) and then uses blocks a few correlation times long.
  - `blocking` repeats the standard error at block lengths 1, 2, 4, 8, ...
    and keeps the value where that error bar stops rising (the plateau).
  - `both` (the default) prints the two error bars and the full blocking table.
    The average itself does not depend on this choice.
- When the file has a `Volume` column, it also reports the cubic box length
  `L = V^(1/3)` that matches the average volume.

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
python utilities/block_average.py path/to/run.prp
python utilities/block_average.py path/to/run.prp --error-method autocorr
python utilities/block_average.py path/to/run.prp --error-method blocking
python utilities/block_average.py path/to/run.prp --json run.averages.json --plot density.png
python utilities/block_average.py path/to/run.prp --plot energy.png --plot-property Energy_Total --sem-band 2
python utilities/block_average.py path/to/run.prp --discard-sweeps 100000
```

`--plot` writes one PNG: the instantaneous trace, a line at the average, and
a shaded band. `--sem-band` is 1, 2, or 3. The default is 1, so the band is
the mean plus or minus one standard error of the mean. When both error bars
were computed, the band uses the larger one. The numbers printed in the
table stay at one standard error. `--plot-property` chooses the curve.
`Mass_Density` still chooses where the average starts.

The method is John Chodera’s equilibration detector (J. Chem. Theory Comput.
2016, 12, 1799): discard the prefix that leaves the largest number of
independent samples.  That count uses the integrated autocorrelation time
(Geyer, Statistical Science 1992, 7, 473).  The blocking plateau (Flyvbjerg
and Petersen, J. Chem. Phys. 1989, 91, 461) is a second opinion on the error
bar only.  The script header explains both.

Check the synthetic cases with:

```bash
python utilities/test_block_average.py
```

## NVT restart from an NPT average

`npt_restart.py` does not modify the NPT directory. It writes a new folder.

The coordinates come from the single-frame `*.restart.xyz` and `*.restart.H`
files (the last configuration). The box volume is the **average** volume from
the block-average report, not the volume of that last frame. Coordinates and
the cell are scaled by `L_average / L_last`, so the starting density is the
average density.

```bash
python utilities/npt_restart.py \
  --inp path/to/r125-prod.inp \
  --prp path/to/R125-T243-P2-prod.prp \
  --out-dir path/to/R125-T243-nvt
```

Restart files are inferred as `<prp stem>.restart.xyz` and `.restart.H` next
to the `.prp`. Pass `--restart-xyz` and `--restart-h` if they live elsewhere.
`--averages run.averages.json` reuses a report you already inspected.

The new input is `nvt_mc`, with `# Box_Info` set to `L_average`, `# Start_Type
read_config` pointed at the scaled XYZ, and the volume move removed. Molecule
and fragment filenames are rewritten relative to `--out-dir`.

### μVT (GCMC) template

```bash
python utilities/npt_restart.py \
  --inp path/to/r125-prod.inp \
  --prp path/to/R125-T243-P2-prod.prp \
  --ensemble gcmc \
  --out-dir path/to/R125-T243-gcmc
```

The box and coordinates are scaled in the same way. The script then writes
placeholders rather than pretending the NPT run produced a chemical potential:

- `# Chemical_Potential_Info` is `0.0` kJ/mol. Replace it before running.
- Insertion and deletion each get probability 0.15. Change that mix if you want.
- The other move probabilities are rescaled so the set still sums to 1.
- `# Molecule_Files` still lists the NPT molecule count as the maximum. Raise
  it if insertions must be allowed to increase N.

Check the scaling logic with:

```bash
python utilities/test_npt_restart.py
```
