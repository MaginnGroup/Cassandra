# Cassandra diagnostics

Post-run tools that check whether a simulation looks **reliable**. These are
**not** full scientific post-processing (block averages, \(g(r)\), production
averages, etc.) — that will live separately later.

| Location | Role |
|----------|------|
| `python/` | Launch Cassandra (`run_cassandra`) |
| `Scripts/` | Setup helpers and legacy regression tests |
| `diagnostics/` | Reliability / teaching checks on output files |
| *(planned)* post-processing | Research analysis: block averages, \(g(r)\), … |

Use after a run produces `.prp`, `.xyz`, etc. Prefer reading those files over
parsing interactive stdout.

## Diagnostic tools

### `property_plotter.py`

Instantaneous property vs MC sweep/step + running average (equilibration eye-check).

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp -p Energy_Total
```

### `bond_distribution.py`

**Fixed-bond integrity check.** Histogram one MCF bond length from `.xyz` vs
MCF \(r_0\) (near-δ expected).

```bash
python diagnostics/bond_distribution.py \
  ~/CassandraV2/projects/hfc125/R125equil.out.xyz \
  ~/CassandraV2/projects/hfc125/r125.mcf --bond 1 --save r125_bond1.png --no-show
```

### `angle_distribution.py` / `dihedral_distribution.py`

Intramolecular PDFs vs ideal / bare Boltzmann (see scripts for teaching notes).

### `msd_com.py`

**COM mean-squared displacement vs MC sweep/step** — exploration / sampling
check for **NVT and NPT only** (fixed \(N\)). **Not** an MD diffusion
coefficient (MC index ≠ physical time; regrowth jumps are allowed).

- Default: average MSD over all molecules of `--species`
- `--molecule M`: that molecule’s \(|\Delta\mathbf{r}|^2\) only
- NPT: companion `.H` required (auto-detected as `*.out.H` next to the XYZ)
- NVT: `.H` if present, else cubic `# Box_Info` from `--inp`

```bash
# NPT example (uses equil.out.H automatically)
python diagnostics/msd_com.py \
  Examples/NPT/pentane/equil.out.xyz \
  Examples/NPT/pentane/pentane.mcf \
  --inp Examples/NPT/pentane/equil.inp \
  --save pentane_msd.png --no-show

# R125 NPT
python diagnostics/msd_com.py \
  ~/CassandraV2/projects/hfc125/R125equil.out.xyz \
  ~/CassandraV2/projects/hfc125/r125.mcf \
  --inp ~/CassandraV2/projects/hfc125/r125-equil.inp \
  --save r125_msd.png --no-show
```

## Shared I/O

```
diagnostics/io/
  prp.py      # .prp
  xyz.py      # multi-frame .xyz
  hfile.py    # companion .H (volume, cell, nmols)
  mcf.py      # Atom / Bond / Angle / Dihedral_Info
  inp.py      # Temperature / Molecule_Files / Sim_Type / Box_Info
  layout.py   # multi-species XYZ atom offsets
```

## Planned diagnostics

| Idea | Notes |
|------|--------|
| Acceptance rates from `.log` | Moves working? |

## Planned post-processing (separate later)

Block averages / error bars, \(g(r)\), production averages, advanced statistics.

## Dependencies

`numpy` and `matplotlib` (`cassandra-dev` conda env).

## Related docs

- [Output formats](../docs/output-formats.md)
- [Docs index](../docs/README.md)
