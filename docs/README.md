# Cassandra V2 Developer Docs

**Start here when you return to the modernization project after time away.**

These notes live on the `modernization` branch and are meant as practical
reminders — not a replacement for the public
[Cassandra documentation](https://cassandra-mc.readthedocs.io/).

## Quick start (60 seconds)

```bash
conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
git branch --show-current          # expect: modernization
git remote -v                      # expect: git@github.com:MaginnGroup/Cassandra.git
ssh-add --apple-use-keychain ~/.ssh/id_ed25519
python run_test.py                 # smoke test
```

## Session start protocol

Use this checklist at the start of every modernization session (laptop). It
avoids the common traps: wrong conda env, wrong branch, and stale local code.

```bash
conda activate cassandra-dev       # not thermo2026 / base
cd ~/CassandraV2/Cassandra
git branch --show-current          # expect: modernization
git status -sb                     # see untracked plots / Examples outputs
git pull                           # when you want latest from origin
git log --oneline -8               # what changed since last time
```

Optional after a `git pull` that touched Fortran, or when something feels off:

```bash
python run_test.py                 # smoke test; see smoke-test.md
```

**Laptop vs CRC reminder**

| Where | Typical work |
|-------|----------------|
| **Laptop** (`cassandra-dev`) | Edit inputs, Git commits, Python diagnostics |
| **maginnfe** | Compile OpenMP binary, `git pull`, `qsub` long runs |

Bring `.prp` / `.xyz` / `.H` home with `scp`/`rsync` and plot on the laptop for
now. On the front end, CRC’s Python module provides `python3` (not always
`python`) — see [V2_CRC_WORKFLOW.md](V2_CRC_WORKFLOW.md) § Python on maginnfe.

Working style for this project: plan → approve → implement; document durable
workflow changes under `docs/`; stage → commit → push on `modernization` (no
PR unless asked). Keep `diagnostics/` for reliability checks only — research
post-processing (block averages, RDF, …) stays separate later.

## Guides

| Doc | What it covers |
|-----|----------------|
| [mcf-setup-workflow.md](mcf-setup-workflow.md) | New molecule checklist: PDB → FF → MCF → fragments → inp |
| [output-formats.md](output-formats.md) | `.xyz` / `.prp` precision, layout, and analysis assumptions |
| [smoke-test.md](smoke-test.md) | How to run `run_test.py`, input/output files, what to expect |
| [V2_DEVELOPMENT_ENVIRONMENT.md](V2_DEVELOPMENT_ENVIRONMENT.md) | Conda environments, compilers, compile & run workflow |
| [V2_CRC_WORKFLOW.md](V2_CRC_WORKFLOW.md) | maginnfe / CRC: Git, rsync, OpenMP, Python module / `.bashrc`, `qsub` |
| [V2_GIT_WORKFLOW.md](V2_GIT_WORKFLOW.md) | Branch strategy, commit/push, SSH authentication |
| [VS-CODE-Guide.md](VS-CODE-Guide.md) | VS Code / Cursor workspace, tasks, and run setup |
| [../python/README.md](../python/README.md) | Python `run_cassandra()` API |
| [../diagnostics/README.md](../diagnostics/README.md) | Post-run diagnostics (property plots, …) |

## Project layout (reminder)

```
~/CassandraV2/
├── Cassandra/              ← git repo (modernization branch)
│   ├── docs/             ← you are here
│   ├── python/           ← V2 Python wrapper
│   ├── diagnostics/      ← post-run checks and plots
│   ├── Scripts/          ← helpers (molecule_to_pdb, mcfgen, library_setup, …)
│   ├── Src/              ← Fortran source; compile here
│   ├── Examples/         ← input files for test runs
│   └── run_test.py       ← smoke test script
├── projects/             ← research runs (e.g. hfc125); not in git
└── Notes on Cassandra modernization   ← personal notes (outside git)
```

## Key facts easy to forget

- V2 work uses branch **`modernization`** on **MaginnGroup/Cassandra** (not a separate repo).
- Use conda env **`cassandra-dev`** for Cassandra — not your class env (`thermo2026`).
- Git push needs **SSH** (`git@github.com:...`), not your GitHub password.
- CRC front end for this group: **`maginnfe.crc.nd.edu`** — see [V2_CRC_WORKFLOW.md](V2_CRC_WORKFLOW.md).
- OpenMP binary name: **`cassandra_gfortran_openMP.exe`** (not `cassandra_gfortran.exe`).
- `python run_test.py` writes output files under `Examples/` — don't `git add .` blindly.
- New molecule setup: follow [mcf-setup-workflow.md](mcf-setup-workflow.md).
- Fluorinated molecules: `# Rcutoff_Low` of 2.0 is often too large (use ~1.0).
