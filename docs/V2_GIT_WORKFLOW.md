# Cassandra V2 Git and GitHub Workflow

This document describes how the Cassandra V2 modernization project is organized in
Git and GitHub, and how to commit and push work safely.

**Other V2 reminders:** [docs index](README.md) ·
[Development environment](V2_DEVELOPMENT_ENVIRONMENT.md) ·
[CRC / maginnfe](V2_CRC_WORKFLOW.md) ·
[Python API](../python/README.md)

## Repository strategy

**There is no separate CassandraV2 GitHub repository.** V2 development happens on
a long-lived branch of the public Cassandra repo:

| Item | Value |
|------|-------|
| GitHub repo | [github.com/MaginnGroup/Cassandra](https://github.com/MaginnGroup/Cassandra) |
| Public default branch | `master` (stable releases) |
| V2 development branch | `modernization` |
| Local clone location | `~/CassandraV2/Cassandra` (your machine) |
| Git remote (SSH) | `git@github.com:MaginnGroup/Cassandra.git` |

### Why a branch, not a new repo?

- `master` stays stable for users and CI while V2 work proceeds.
- All history is preserved in one place.
- When V2 is ready, you merge `modernization` into `master` (or open a PR) to
  replace the public code — no repo migration required.

The parent folder `~/CassandraV2/` (notes, workspace files, etc.) is **outside**
the git repo and is for your personal project notes only.

---

## Returning after time away? Start here.

Run this 30-second checklist before you edit or push:

```bash
cd ~/CassandraV2/Cassandra

# 1. Correct branch?
git branch --show-current          # expect: modernization

# 2. Remote using SSH (not HTTPS)?
git remote -v
# GOOD:  git@github.com:MaginnGroup/Cassandra.git
# BAD:   https://github.com/MaginnGroup/Cassandra.git  ← will ask for password

# 3. SSH key loaded? (may be needed after a reboot)
ssh-add --apple-use-keychain ~/.ssh/id_ed25519

# 4. GitHub recognizes you?
ssh -T git@github.com
# expect: Hi ejmaginn! You've successfully authenticated...
```

If step 2 shows `https://`, fix it once (see [One-time SSH setup](#one-time-ssh-setup)).

If step 3 asks `Enter passphrase`, type the passphrase you set when the key was
created (characters will not appear as you type — that is normal). After one
successful entry, macOS Keychain usually remembers it.

---

## One-time setup

### Clone and branch (first time only)

```bash
git clone git@github.com:MaginnGroup/Cassandra.git   # use SSH from the start
cd Cassandra
git checkout modernization    # branch already exists on GitHub
git branch                  # expect: * modernization
```

If you cloned with HTTPS earlier, switch to SSH:

```bash
git remote set-url origin git@github.com:MaginnGroup/Cassandra.git
```

### One-time SSH setup

GitHub does **not** accept your account password for git operations. Use SSH.

**1. Confirm you have a key**

```bash
ls ~/.ssh/id_ed25519.pub
```

**2. Add the public key to GitHub** (skip if already done)

```bash
pbcopy < ~/.ssh/id_ed25519.pub
```

Open [github.com/settings/keys](https://github.com/settings/keys) → **New SSH key**
→ paste → save.

**3. Load the key and store passphrase in Keychain**

```bash
ssh-add --apple-use-keychain ~/.ssh/id_ed25519
```

You may be prompted: `Enter passphrase for /Users/.../id_ed25519:` — enter the
passphrase you chose when the key was created. A wrong passphrase shows
`Bad passphrase, try again`; keep trying until you see:

```
Identity added: /Users/.../id_ed25519 (ed@nd.edu)
```

**4. Point this repo at SSH** (critical — easy to miss)

```bash
cd ~/CassandraV2/Cassandra
git remote set-url origin git@github.com:MaginnGroup/Cassandra.git
git remote -v    # confirm both fetch and push show git@github.com
```

**5. Verify**

```bash
ssh -T git@github.com
git push origin modernization   # should NOT ask for username/password
```

> **Common mistake:** Steps 3 and 5 can succeed (`Identity added`, `Hi ejmaginn!`)
> but `git push` still asks for username/password if step 4 was skipped.
> Loading the SSH key and switching the remote URL are **both** required.

### Conda build environment

Use the **`cassandra-dev`** conda environment. See
[V2_DEVELOPMENT_ENVIRONMENT.md](V2_DEVELOPMENT_ENVIRONMENT.md) for full setup,
switching between environments, and compile/run workflow.

```bash
conda activate cassandra-dev
cd Src
make -f Makefile.gfortran
```

---

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

#### Fortran revision history

When you change a Fortran routine, add a **brief dated note** in that routine’s
header comment block (create a short header if none exists). If the file also
has a module-level `Revision history` section, add a matching one-line entry
there.

Use this style (match existing `(EJM)` entries in `Src/`):

```text
! Revision history
!   08/13/26 (EJM) : Brief Relative_Error explanation above the energy table
```

When adding an explanatory routine header (not a behavior change), prefer:

```text
!   08/13/26 (EJM) : Added explanatory header for <short topic>
```

Avoid calling these “teaching” headers in the revision line. Keep entries to one
line when possible. This is part of the edit, not a separate commit.

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

Before pushing, confirm the remote is SSH:

```bash
git remote get-url origin
# expect: git@github.com:MaginnGroup/Cassandra.git
```

Then push:

```bash
git push origin modernization
```

First push of a new branch:

```bash
git push -u origin modernization
```

No username or password prompt should appear.

### 6. Verify on GitHub

Open:

[github.com/MaginnGroup/Cassandra/tree/modernization](https://github.com/MaginnGroup/Cassandra/tree/modernization)

Confirm your latest commit appears at the top.

---

## Troubleshooting push and authentication

### `git push` asks for username and password

**Cause:** The repo remote is still HTTPS.

**Fix:**

```bash
git remote set-url origin git@github.com:MaginnGroup/Cassandra.git
git remote -v
git push origin modernization
```

### `Invalid username or token. Password authentication is not supported`

**Cause:** You entered your GitHub account password over HTTPS.

**Fix:** Switch to SSH (above), or use a [Personal Access Token](#alternative-personal-access-token) as the password instead.

### `Enter passphrase for ... id_ed25519`

**Cause:** Normal. Your SSH private key is encrypted.

**Fix:** Enter the passphrase you set when the key was created. Use
`ssh-add --apple-use-keychain` so macOS remembers it.

### `Bad passphrase, try again`

**Cause:** Wrong passphrase entered.

**Fix:** Try again. If you cannot recover the passphrase, create a new key (see
[One-time SSH setup](#one-time-ssh-setup)) and add the new `.pub` file to GitHub.

### `Permission denied (publickey)` from `ssh -T git@github.com`

**Cause:** Key not loaded and/or not registered on GitHub.

**Fix:**

```bash
ssh-add --apple-use-keychain ~/.ssh/id_ed25519
pbcopy < ~/.ssh/id_ed25519.pub   # add at github.com/settings/keys if missing
ssh -T git@github.com
```

### `could not read Username for 'https://github.com': Device not configured`

**Cause:** HTTPS remote in a terminal that cannot prompt (e.g. Cursor integrated terminal).

**Fix:** Switch to SSH, or push from Terminal.app.

### Alternative: Personal Access Token

If you prefer HTTPS, create a token at
[github.com/settings/tokens](https://github.com/settings/tokens) with `repo` scope.
Use it as the **password** (not your GitHub account password) when prompted.
macOS Keychain will store it.

### Alternative: GitHub CLI

```bash
brew install gh
gh auth login
git push origin modernization
```

---

## Endgame: replacing public Cassandra

When V2 is ready for release:

1. Ensure `modernization` passes the test suite and smoke tests.
2. Open a pull request: `modernization` → `master` on GitHub.
3. Review, merge, and tag a release.
4. Update conda-forge / documentation as needed.

Until that merge, users on `master` are unaffected.

---

## Quick reference

```bash
# Returning after time away
cd ~/CassandraV2/Cassandra
git branch --show-current
git remote -v
ssh-add --apple-use-keychain ~/.ssh/id_ed25519
ssh -T git@github.com

# Where am I?
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

---

## Related docs

- [docs/README.md](README.md) — index of all V2 developer reminders
- [V2_DEVELOPMENT_ENVIRONMENT.md](V2_DEVELOPMENT_ENVIRONMENT.md) — conda, compile, run
- `python/README.md` — Python wrapper usage
- `CONTRIBUTING.md` — upstream contribution guidelines for `master`
- `~/CassandraV2/Notes on Cassandra modernization` — personal notes (outside repo)
