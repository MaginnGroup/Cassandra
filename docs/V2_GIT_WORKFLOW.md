# Cassandra V2 Git and GitHub Workflow

This document describes how the Cassandra V2 modernization project is organized in
Git and GitHub, and how to commit and push work safely.

## Repository strategy

**There is no separate CassandraV2 GitHub repository.** V2 development happens on
a long-lived branch of the public Cassandra repo:

| Item | Value |
|------|-------|
| GitHub repo | [github.com/MaginnGroup/Cassandra](https://github.com/MaginnGroup/Cassandra) |
| Public default branch | `master` (stable releases) |
| V2 development branch | `modernization` |
| Local clone location | `~/CassandraV2/Cassandra` (your machine) |

### Why a branch, not a new repo?

- `master` stays stable for users and CI while V2 work proceeds.
- All history is preserved in one place.
- When V2 is ready, you merge `modernization` into `master` (or open a PR) to
  replace the public code — no repo migration required.

The parent folder `~/CassandraV2/` (notes, workspace files, etc.) is **outside**
the git repo and is for your personal project notes only.

## One-time setup

If you already have the clone and are on `modernization`, skip to
[Daily workflow](#daily-workflow).

```bash
# Clone (once)
git clone https://github.com/MaginnGroup/Cassandra.git
cd Cassandra

# Create and switch to the V2 branch (once)
git checkout -b modernization

# Verify you are on the right branch
git branch
# Expected: * modernization
```

### Conda build environment (optional but recommended)

```bash
conda create --name cassandra-dev python=3.10
conda activate cassandra-dev
conda install -c conda-forge cmake pandas gfortran
```

Compile from `Src/`:

```bash
cd Src
make -f Makefile.gfortran
```

## Daily workflow

Always confirm your branch before editing:

```bash
git branch --show-current
# Must print: modernization
```

### 1. Edit code

V2 Python wrapper files live in:

- `python/api.py` — subprocess launcher
- `python/README.md` — API usage
- `run_test.py` — smoke test from repo root

See `python/README.md` for how to run simulations.

### 2. Review changes

```bash
git status
git diff                  # unstaged changes
git diff --staged         # staged changes
```

**Important:** Running `python run_test.py` regenerates simulation output files
under `Examples/`. Review `git status` carefully before staging.

### 3. Stage selectively

Prefer staging specific paths instead of `git add .`:

```bash
# Good — only source changes
git add python/ run_test.py docs/

# Risky — may include simulation outputs and local test edits
git add .
```

#### What to commit

| Commit | Do not commit |
|--------|---------------|
| `python/` source and docs | `*.exe` (compiled binary; already gitignored) |
| `run_test.py` | Local simulation outputs you generated while testing |
| Fortran / build fixes in `Src/` | Accidental edits to example `.inp` files from test runs |
| Intentional example updates | Files under `~/CassandraV2/` outside the repo |

The repo **does** track some reference output under `Examples/` and
`Scripts/testSuite/Resources/exampleResults/` for regression tests. Do not add
new output files from your local smoke tests unless you deliberately intend to
update a reference example.

### 4. Commit

```bash
git commit -m "Short summary of why the change was made"
```

Write commit messages in plain language focused on **why**, not just what changed.

### 5. Push to GitHub

```bash
git push origin modernization
```

First push of a new branch:

```bash
git push -u origin modernization
```

### 6. Verify on GitHub

Open:

[github.com/MaginnGroup/Cassandra/tree/modernization](https://github.com/MaginnGroup/Cassandra/tree/modernization)

Confirm your latest commit appears at the top.

## Authentication (push failures)

If push fails with:

```
fatal: could not read Username for 'https://github.com': Device not configured
```

the remote is using HTTPS and your terminal cannot prompt for credentials.
This often happens in IDE-integrated terminals. Fix options:

### Option A: Push from a regular Terminal.app / iTerm window

```bash
cd ~/CassandraV2/Cassandra
git push origin modernization
```

macOS Keychain (`credential.helper=osxkeychain`) should supply stored credentials.

### Option B: Use SSH instead of HTTPS

```bash
git remote set-url origin git@github.com:MaginnGroup/Cassandra.git
ssh -T git@github.com   # verify GitHub recognizes your key
git push origin modernization
```

### Option C: GitHub CLI

```bash
brew install gh
gh auth login
git push origin modernization
```

### Option D: Personal access token (HTTPS)

Create a token at GitHub → Settings → Developer settings → Personal access tokens,
then use it as the password when `git push` prompts for credentials.

## Endgame: replacing public Cassandra

When V2 is ready for release:

1. Ensure `modernization` passes the test suite and smoke tests.
2. Open a pull request: `modernization` → `master` on GitHub.
3. Review, merge, and tag a release.
4. Update conda-forge / documentation as needed.

Until that merge, users on `master` are unaffected.

## Current status (as of July 2026)

Run these commands locally to see where you stand:

```bash
git branch --show-current
git status -sb
git log origin/modernization..HEAD --oneline   # local commits not yet pushed
```

### Known issues to clean up before pushing

At the time this document was written, the local `modernization` branch had
**two unpushed commits** ahead of `origin/modernization`:

1. `92ed5b1` — Python API documentation and error handling (good to push)
2. `4f195f7` — Same topic, but also accidentally includes:
   - ~58,000 lines of local simulation output (`water_spc_nvt*.out.*`)
   - Unintended edits to `Examples/NVT/water_spc/nvt.inp` (longer test run)
   - Deletion of `Notes on Cassandra modernization` from the repo

**Do not push commit `4f195f7` as-is.** Clean up first using the steps below.

### Cleanup: keep the good commit, drop the bad one

From the repo root, with a clean understanding that this rewrites local history
(not yet pushed, so it is safe):

```bash
# 1. Move branch pointer back to the last good pushed commit,
#    keeping all file changes in your working tree
git reset --soft origin/modernization

# 2. Restore example input to the last pushed version
git restore --staged --worktree Examples/NVT/water_spc/nvt.inp

# 3. Unstage simulation outputs (leave files on disk, just don't commit them)
git restore --staged Examples/NVT/water_spc/water_spc_nvt*.out.*

# 4. Stage only the V2 source changes
git add python/ run_test.py docs/V2_GIT_WORKFLOW.md

# 5. Review what will be committed
git status
git diff --staged --stat

# 6. Commit and push
git commit -m "Add Python API docs, error handling, and V2 git workflow guide"
git push origin modernization
```

If `git restore` does not unstage the output files, use:

```bash
git reset HEAD Examples/NVT/water_spc/water_spc_nvt*.out.*
```

## Quick reference

```bash
# Where am I?
git branch --show-current
git status -sb

# What is not on GitHub yet?
git log origin/modernization..HEAD --oneline

# Safe commit cycle
git add python/ run_test.py docs/
git commit -m "Describe why"
git push origin modernization

# Run smoke test (generates local output — do not commit blindly)
python run_test.py
```

## Related docs

- `python/README.md` — Python wrapper usage
- `CONTRIBUTING.md` — upstream contribution guidelines for `master`
- `~/CassandraV2/Notes on Cassandra modernization` — personal notes (outside repo)
