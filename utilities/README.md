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
  `L = <V>^(1/3)`.
- When both `Mass_Density` and `Volume` are present it reports two densities.
  The mean of the density column, `<rho>`, is proportional to `<1/V>`.
  The density at the mean volume, `(total mass)/<V>`, is proportional to
  `1/<V>`.  Volume fluctuations make `<rho>` slightly larger.  The NVT
  restart uses `<V>`, so its density is the second number.

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
files (the last configuration). The box volume is the **mean volume** `<V>`
from the block-average report, not the volume of that last frame. The cell
changes from `L_last` to `<V>^(1/3)`. The script reconstructs each molecule
across periodic boundaries, maps its mass-weighted center of mass to the same
fractional coordinate in the new cell, and translates the whole molecule.
It does not scale internal coordinates, so bond lengths, angles, and
dihedrals remain unchanged.

The NVT density is the density at that mean volume, `(total mass)/<V>`.
It is not the mean of the NPT `Mass_Density` column. That column average is
proportional to `<1/V>`, and `<1/V>` is slightly larger than `1/<V>` when the
NPT volume fluctuates. `block_average.py` prints both numbers. The difference
is normally smaller than the uncertainty on `<rho>`.

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
read_config` pointed at the COM-mapped XYZ, and the volume move removed.
Molecule and fragment filenames are rewritten relative to `--out-dir`.
The tool supports one cubic box, matching Cassandra’s NPT volume-move
restriction; it stops rather than guessing how to transform a non-cubic cell.

You will likely need to edit the new `.inp` before submitting. Below,
`[run_name]` means its `# Run_Name` (for example `R125-T243-NVT`). The
generated start file is `[run_name].xyz`, which is also the name Cassandra
uses for the coordinate movie. Rename it and its `.H` companion to
`[run_name].start.xyz` and `[run_name].start.H`, and point `read_config`
at the `.start.xyz`. How Cassandra names every output file, and how to
resume a killed run, is in
[docs/output-formats.md](../docs/output-formats.md).

### μVT (GCMC) template

```bash
python utilities/npt_restart.py \
  --inp path/to/r125-prod.inp \
  --prp path/to/R125-T243-P2-prod.prp \
  --ensemble gcmc \
  --out-dir path/to/R125-T243-gcmc
```

The box changes and molecular centers of mass are mapped in the same way. The
script then writes
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

## Radial distribution and coordination number

`rdf.py` reads a coordinate movie, the matching `.H` file, and the species
MCF.

- `--mode com` is the center-of-mass RDF. Each molecule is rebuilt across
  the periodic boundary before its center is calculated.
- `--mode atom` is an atom–atom RDF. Choose each side with `--element-i`,
  `--atom-type-i`, or `--atom-index-i`, and the matching `--element-j`,
  `--atom-type-j`, or `--atom-index-j`. Prefer exact MCF atom types—for
  example `c3_1_s1`—when sites with the same element need separate RDFs.
  Omitting the `j` selector pairs the selection with itself.
- Pairs inside one molecule are left out unless `--include-intramolecular`
  is set. COM pairs are always between different molecules. Atom pairs use
  the per-frame minimum-image distance directly; only COMs need molecular
  reconstruction across a periodic boundary.
- `g(r)` divides the pair count by the ideal-gas shell. A selection paired
  with itself uses `N-1` neighbors so `g(r)` approaches 1 at large `r`.
  With intramolecular atom pairs excluded, the ideal count and coordination
  density include only sites on other molecules. This is exact even when an
  atom type occurs more than once in a molecule.
- The coordination number integrates `4 pi r^2 rho g(r)` out to `--r-cut`.
  Omit `--r-cut` and the cutoff is the lowest bin between the second and
  third crossings of `g(r)=1`: out of the first peak and into the second
  shell. The CN figure axis ends at that cutoff. The RDF and combined figures
  run to `--r-max`.
- The standard error is the scatter across frames. Those frames are
  correlated, so the error bar is not a production uncertainty.

`--r-max` is how far the table goes. The default is half the shortest box
edge. `--r-cut` is only the coordination integral.

With `--plot`, the script writes a table and a figure next to the movie.
The pair name is `com`, or the MCF atom types, such as `c3_1_s1-f_2_s1`.
Pass `--atom-type-i` and, for a different second type, `--atom-type-j`. For
`[run_name].xyz` the COM files are `[run_name]_com.rdf`,
`[run_name]_com_rdf.png`, `[run_name]_com_cn.png`, and
`[run_name]_com_rdf_cn.png`. `--plot` chooses `rdf`, `cn`, or `both`
(g(r) on the left axis and n(r) on the right). The RDF figure does not mark
the coordination cutoff unless `--mark-r-cut` is set. Default titles name
the pair, such as `COM-COM RDF, species 1` and `COM-COM CN, species 1`.
Replace the figure you are writing with `--title-rdf`, `--title-cn`, or
`--title-both`. `--legend` is `none`, `curves`, or `with-cutoff`. RDF and CN
figures have no legend unless you set one. The combined figure names
`g(r)`, `n(r)`, and the cutoff. The uncertainty band is never listed.
`--x-unit` is `nm` by default, or `A` for angstroms, drawn as Å. `--r-cut`,
`--r-max`, and the table stay in angstroms. `--out` replaces the table path.

```bash
python utilities/rdf.py movie.xyz species.mcf --h movie.H --mode com --r-cut 6 --plot both
python utilities/rdf.py movie.xyz species.mcf --h movie.H --mode atom --atom-type-i c3_1_s1 --plot both
python utilities/rdf.py movie.xyz species.mcf --h movie.H --mode atom --atom-type-i c3_1_s1 --atom-type-j c3_2_s1 --plot both
```

Check the synthetic cases with:

```bash
python utilities/test_rdf.py
```
