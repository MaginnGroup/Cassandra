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

## Guides

| Doc | What it covers |
|-----|----------------|
| [smoke-test.md](smoke-test.md) | How to run `run_test.py`, input/output files, what to expect |
| [V2_DEVELOPMENT_ENVIRONMENT.md](V2_DEVELOPMENT_ENVIRONMENT.md) | Conda environments, compilers, compile & run workflow |
| [V2_GIT_WORKFLOW.md](V2_GIT_WORKFLOW.md) | Branch strategy, commit/push, SSH authentication |
| [../python/README.md](../python/README.md) | Python `run_cassandra()` API |

## Project layout (reminder)

```
~/CassandraV2/
├── Cassandra/              ← git repo (modernization branch)
│   ├── docs/             ← you are here
│   ├── python/           ← V2 Python wrapper
│   ├── Src/              ← Fortran source; compile here
│   ├── Examples/         ← input files for test runs
│   └── run_test.py       ← smoke test script
└── Notes on Cassandra modernization   ← personal notes (outside git)
```

## Key facts easy to forget

- V2 work uses branch **`modernization`** on **MaginnGroup/Cassandra** (not a separate repo).
- Use conda env **`cassandra-dev`** for Cassandra — not your class env (`thermo2026`).
- Git push needs **SSH** (`git@github.com:...`), not your GitHub password.
- `python run_test.py` writes output files under `Examples/` — don't `git add .` blindly.
