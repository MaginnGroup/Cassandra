# VS Code Guide for Cassandra

## Recommended setup

- Open `CassandraV2` as the workspace root.
- Install VS Code extensions:
  - Python
  - Makefile Tools
  - Modern Fortran or Fortran
  - CodeLLDB (for debugging native builds)
  - GitLens (optional)
  - Emacs keymap (optional)

The workspace files under `CassandraV2/.vscode/` are outside the Cassandra git
repo. They are local editor settings, not part of `pixi-env`.

## Workspace settings

`.vscode/settings.json` should point Python at the Pixi interpreter, after
`pixi install` has created it:

```text
~/CassandraV2/Cassandra/.pixi/envs/default/bin/python
```

Do not set `python.terminal.activateEnv` to a Conda environment. An active
`cassandra-dev` shell overrides Pixi's compiler.

## Tasks

The workspace tasks run through Pixi, from the Cassandra repo:

1. **Build Cassandra (gfortran)**
   - `pixi run build` (`make -f Makefile.gfortran` in `Cassandra/Src`)

2. **Run NVT water_spc**
   - Runs `Cassandra/Src/cassandra_gfortran.exe` on `nvt.inp`
   - Working directory: `Cassandra/Examples/NVT/water_spc`

3. **Run NVT water_spc and save log**
   - Same as above, and writes `CassandraV2/run_logs/nvt_water_spc.log`

4. **Run any example**
   - Prompts for a relative example folder (for example `NVT/water_spc`)
   - Prompts for the input filename inside that folder
   - Prompts for a log filename to save under `CassandraV2/run_logs`

## Running tasks

- Build: `Shift+Cmd+B`
- Run task: `Command Palette → Tasks: Run Task → Run NVT water_spc`
- Run with log: `Command Palette → Tasks: Run Task → Run NVT water_spc and save log`

## Running examples manually

From the Cassandra repo, outside any activated Conda environment:

```bash
cd ~/CassandraV2/Cassandra
pixi run build
cd Examples/NVT/water_spc
pixi run -- ../../Src/cassandra_gfortran.exe nvt.inp
```

Pixi finds `pixi.toml` in the parent repo. Replace the folder and input file
for a different example. The smoke test is `pixi run smoke` from the repo root.

## Notes capture

Use `run_logs/` to keep run outputs, and add notes in `Cassandra/docs/`.

Example note format:

- Example: `Examples/NVT/water_spc`
- Command: `pixi run -- ../../Src/cassandra_gfortran.exe nvt.inp`
- Output log: `run_logs/nvt_water_spc.log`
- Observations: ...
