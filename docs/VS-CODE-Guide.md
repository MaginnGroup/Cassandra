# VS Code Guide for Cassandra Modernization

## Recommended setup

- Open `CassandraV2` as the workspace root.
- Install VS Code extensions:
  - Python
  - Makefile Tools
  - Modern Fortran or Fortran
  - CodeLLDB (for debugging native builds)
  - GitLens (optional)
  - Emacs keymap (optional)

## Workspace settings

The workspace file `.vscode/settings.json` is configured to:

- Use `/opt/anaconda3/envs/cassandra-dev/bin/python`
- Activate Conda in the integrated terminal
- Find `conda` at `/opt/anaconda3/bin/conda`
- Hide `.git`, `build`, and `dist` folders in Explorer

## Tasks

The workspace has three tasks:

1. **Build Cassandra (gfortran)**
   - Runs `make -f Makefile.gfortran` inside `Cassandra/Src`
   - Activates the `cassandra-dev` Conda environment first

2. **Run NVT water_spc**
   - Runs the built executable `Cassandra/Src/cassandra_gfortran.exe`
   - Passes `nvt.inp` as the input file
   - Runs from `Cassandra/Examples/NVT/water_spc`

3. **Run NVT water_spc and save log**
   - Same as above, but also writes output to `CassandraV2/run_logs/nvt_water_spc.log`

4. **Run any example**
   - Prompts for a relative example folder (for example `NVT/water_spc`)
   - Prompts for the input filename inside that folder
   - Prompts for a log filename to save under `CassandraV2/run_logs`
   - Runs the chosen example with the built executable and saves output

## Running tasks

- Build: `Shift+Cmd+B`
- Run task: `Command Palette → Tasks: Run Task → Run NVT water_spc`
- Run with log: `Command Palette → Tasks: Run Task → Run NVT water_spc and save log`

## Running examples manually

If you prefer terminal commands:

```bash
conda activate cassandra-dev
cd Cassandra/Src
make -f Makefile.gfortran
cd ../Examples/NVT/water_spc
../../Src/cassandra_gfortran.exe nvt.inp
```

If you want a specific example, replace the folder and input file accordingly.

## Notes capture

Use `run_logs/` to keep run outputs, and add notes in `Cassandra/docs/VS-CODE-Guide.md` or a new notes file.

Example note format:

- Example: `Examples/NVT/water_spc`
- Command: `./cassandra_gfortran.exe nvt.inp`
- Output log: `run_logs/nvt_water_spc.log`
- Observations: ...
