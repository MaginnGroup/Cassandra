# Cassandra diagnostics

Post-run tools for checking whether a simulation looks reasonable. These are
**analysis** utilities, separate from:

| Location | Role |
|----------|------|
| `python/` | Launch Cassandra (`run_cassandra`) |
| `Scripts/` | Setup helpers and legacy regression tests |
| `diagnostics/` | Scientific checks and plots on output files |

Use after a run produces `.prp`, `.xyz`, etc. Prefer reading those files over
parsing interactive stdout.

## Tools

### `property_plotter.py`

Plot one column from a `.prp` file vs MC sweep/step: instantaneous values plus
a running (cumulative) average. Prints the final mean to the terminal.

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra

# List columns
python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp --list

# Interactive choice (TTY)
python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp

# Named property (case-insensitive; unique substring OK)
python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp -p Energy_Total
python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp -p enthalpy --save enthalpy.png

# Headless (CI / SSH without display)
python diagnostics/property_plotter.py Examples/NVT/water_spc/nvt.out.prp \
  -p Pressure --save pressure.png --no-show
```

Shared `.prp` parser: `diagnostics/io/prp.py` (reuse this in future tools).

## Layout

```
diagnostics/
  README.md                 ← this file
  property_plotter.py       ← CLI: property vs step
  io/
    prp.py                  ← shared .prp reader
    # xyz.py                ← planned
  # angle_distribution.py   ← planned
```

## Planned diagnostics

| Tool idea | Inputs |
|-----------|--------|
| Bond / angle PDF vs Boltzmann | `.xyz` + MCF (+ `.inp` / `.log` for \(T\)) |
| Molecular displacement over time | `.xyz` |
| More property statistics (burn-in, block averages) | `.prp` |

## Dependencies

`numpy` and `matplotlib` (available in the `cassandra-dev` conda environment).

## Related docs

- [Output formats](../docs/output-formats.md) — `.prp` / `.xyz` layout
- [Docs index](../docs/README.md)
