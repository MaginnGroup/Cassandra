# Cassandra V2 — CRC / maginnfe Workflow

Practical reminder for running Cassandra on Notre Dame’s Center for Research
Computing (CRC), especially from the Maginn group front end **maginnfe**.

**Related:** [docs index](README.md) · [Dev environment (laptop)](V2_DEVELOPMENT_ENVIRONMENT.md) ·
[Git / GitHub](V2_GIT_WORKFLOW.md) · [Diagnostics](../diagnostics/README.md) ·
CRC docs: [Quick Start](https://docs.crc.nd.edu/new_user/quick_start.html) ·
[Submitting Batch Jobs](https://docs.crc.nd.edu/new_user/submitting_batch_jobs.html)

First worked through for R125 NPT equilibration on **2026-08-07** (EJM).

---

## Roles: laptop vs CRC

| Where | Use for |
|-------|---------|
| **Laptop** | Edit inputs, short tests, Git commits, Python diagnostics (`property_plotter`, `angle_distribution`, `dihedral_distribution`) |
| **maginnfe** | Compile Cassandra, `git pull`, submit jobs with `qsub` |
| **Compute nodes** | Long MC runs (via the UGE queue — still often called “SGE”) |

Do **not** treat the front end as a place for multi-day interactive production runs.
Prefer `qsub` even though maginnfe allows longer interactive jobs than public front ends.

```
Laptop  --git push-->  GitHub  --git clone/pull-->  maginnfe
Laptop  --rsync/scp-->  project files (inp, mcf, chk, fragments)
maginnfe --qsub-->  long queue  -->  compute node  -->  .prp / .chk / .xyz
maginnfe --scp-->  laptop  -->  diagnostics/{property_plotter,angle_distribution,dihedral_distribution}.py
```

---

## Account facts (Maginn / ed)

| Item | Value |
|------|--------|
| NetID | `ed` |
| Front end | `maginnfe.crc.nd.edu` (use this, not `crcfe01` / `crcfe02`) |
| Home | `/users/ed` (100 GB quota) |
| Login | `ssh ed@maginnfe.crc.nd.edu` |
| First production queue | Public `long` (`#$ -q long`) |
| Later | Maginn-owned / priority host groups (document when you switch) |

```bash
ssh ed@maginnfe.crc.nd.edu
hostname    # maginnfe.crc.nd.edu
quota
pwd         # /users/ed
```

Off campus: ND VPN (or CRC’s off-campus instructions) before SSH.

---

## 1. GitHub access from maginnfe

The SSH key on your **laptop** does **not** unlock GitHub from maginnfe. The
cluster needs its **own** key registered with GitHub.

### 1a. Returning session — key already exists?

```bash
ls ~/.ssh/id_ed25519_github.pub
```

If that prints a path (e.g. `/users/ed/.ssh/id_ed25519_github.pub`), **skip
keygen**. Do this every new login (or put the `ssh-add` in `~/.bashrc`):

```bash
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519_github

# Ensure GitHub always uses this key
mkdir -p ~/.ssh
chmod 700 ~/.ssh
cat > ~/.ssh/config << 'EOF'
Host github.com
  HostName github.com
  User git
  IdentityFile ~/.ssh/id_ed25519_github
  IdentitiesOnly yes
EOF
chmod 600 ~/.ssh/config

ssh -T git@github.com
# expect: Hi ejmaginn! You've successfully authenticated...
```

| `ssh -T` result | What to do |
|-----------------|------------|
| `Hi … You've successfully authenticated` | Done — go to §2 |
| `Permission denied (publickey)` | Public key is missing from GitHub, or wrong `IdentityFile` — see §1c |
| Asks for passphrase | Type the passphrase you set at `ssh-keygen` (empty is fine if you used none) |

### 1b. First time only — create the key

```bash
ssh-keygen -t ed25519 -C "ed@maginnfe-github" -f ~/.ssh/id_ed25519_github
# empty passphrase is fine for convenience on the cluster
eval "$(ssh-agent -s)"
ssh-add ~/.ssh/id_ed25519_github
cat ~/.ssh/id_ed25519_github.pub
```

Then do §1c and finish with the `~/.ssh/config` + `ssh -T` steps in §1a.

### 1c. Register the public key on GitHub (browser / laptop)

```bash
# On maginnfe — copy this entire line
cat ~/.ssh/id_ed25519_github.pub
```

On your laptop: GitHub → **Settings → SSH and GPG keys → New SSH key**  
Title: `maginnfe` → paste the `.pub` line → Save.

If a key titled `maginnfe` already exists, confirm its fingerprint matches this
`.pub` file; if not, add a new key or replace the old one.

Never share the **private** key (`id_ed25519_github` without `.pub`).

See also [V2_GIT_WORKFLOW.md](V2_GIT_WORKFLOW.md) for laptop Git habits.

---

## 2. Clone or update Cassandra (`modernization`) on maginnfe

**Source code** travels via GitHub. Do not SFTP the whole `Src/` tree every time.

### First clone

```bash
mkdir -p /users/ed/CassandraV2
cd /users/ed/CassandraV2
git clone git@github.com:MaginnGroup/Cassandra.git
cd Cassandra
git checkout modernization
git pull
```

### Later updates (usual case after a laptop push)

```bash
cd /users/ed/CassandraV2/Cassandra   # adjust if your clone path differs

git status
# If Example smoke outputs are dirty, discard them (do not commit):
git restore Examples/NVT/water_spc/nvt.out.* 2>/dev/null || true

git checkout modernization
git pull origin modernization
# then recompile if Fortran changed (see §3)
```

`Your branch is up to date with 'origin/modernization'` only means **local
tracking info** is current — it does **not** fetch new commits until `git pull`
succeeds. If pull fails with `Permission denied (publickey)`, fix §1 first.

### After a successful pull with Fortran changes

```bash
cd Src
module purge
module load gcc/15.2.0
make -f Makefile.gfortran.openMP clean    # or Makefile.gfortran if that is what you use
make -f Makefile.gfortran.openMP
```

(Use the same compiler module you used for the previous build.)

Layout on CRC mirrors the laptop idea:

```
/users/ed/CassandraV2/
├── Cassandra/          ← git clone (modernization)
└── projects/
    └── hfc125/         ← R125 inputs / outputs (NOT in git)
```

---

## 3. Copy project files (rsync / scp)

Simulation directories (`.inp`, `.mcf`, fragments, `.chk`) live **outside** the
Cassandra repo. Transfer them explicitly.

`rsync` will **not** create missing parent directories. Create them first:

```bash
# on maginnfe
mkdir -p /users/ed/CassandraV2/projects/hfc125
```

From the **laptop**:

```bash
rsync -avz --progress \
  --exclude '*.xyz' --exclude '*.xlsx' --exclude '~$*' \
  ~/CassandraV2/projects/hfc125/ \
  ed@maginnfe.crc.nd.edu:/users/ed/CassandraV2/projects/hfc125/
```

Keep at least: input `.inp`, `.mcf`, force-field files, `species1/` fragments,
and the `.chk` you need to restart. Large `.xyz` movies can wait until you need them.

GUI alternative: [Cyberduck](https://cyberduck.io/) (CRC recommends it for Mac).

---

## 4. Compile with OpenMP (modules, not conda)

On CRC, use **modules** for the Fortran toolchain. Do **not** compile Cassandra
with a conda `gfortran` on maginnfe (mismatched compilers break `.mod` files).

Confirmed modules on maginnfe (2026-08):

| Need | Module |
|------|--------|
| gfortran | Comes with `gcc/15.2.0` (no separate `gfortran` module) |
| cmake | `cmake` (default `cmake/4.1.0`) |

```bash
module load gcc/15.2.0
module load cmake
which gfortran && gfortran --version

cd /users/ed/CassandraV2/Cassandra/Src
make -f Makefile.gfortran.openMP clean
make -f Makefile.gfortran.openMP
ls -la cassandra_gfortran_openMP.exe
```

### Pitfall: executable name

| Makefile | Executable |
|----------|------------|
| `Makefile.gfortran` | `cassandra_gfortran.exe` |
| `Makefile.gfortran.openMP` | **`cassandra_gfortran_openMP.exe`** |

Job scripts and interactive tests must use the **`_openMP`** name if that is
what you built. A first `qsub` failed for this reason (2026-08-07).

Load the **same** `gcc` module in the batch script that you used to compile.

Optional later: try Intel (`Makefile.intel.openMP`) for speed — only after a
working gfortran OpenMP path exists.

---

## 5. Checkpoint restart (important semantics)

### `# Start_Type`

```text
# Start_Type
checkpoint R125equil.out.chk
```

Seeds in the `.inp` are overwritten by those in the checkpoint (Cassandra docs).

### Equilibration → **new** production run (not checkpoint)

Do **not** use the multi-frame movie `run_name.xyz` with `read_config` — that
loads the **first** frame. Use the overwrite restart pair written every
`coord_freq` (and again at a clean finish):

| File | Role |
|------|------|
| `run_name.restart.xyz` | Latest single-frame coordinates |
| `run_name.restart.H` | Latest volume, H-matrix, species molecule counts |

GEMC (2 boxes): `run_name.box1.restart.xyz` / `.H` and `run_name.box2.restart.*`.

Production example:

```text
# Start_Type
read_config  200  equil.out.restart.xyz
```

Set `# Box_Info` from `equil.out.restart.H` (especially after NPT). Molecule
counts on the `read_config` line must match that `.H`. If the companion
`*.restart.H` is present, Cassandra **aborts** when volume/cell or \(N\)
disagree with the input (legacy XYZ without a sibling `.H` skips the check).

Checkpoint continues the **same** run (RNG, step counter). Restart XYZ/H
starts a **fresh** production run with new seeds.

### `run` is a **total** end point, not “how many more”

If the checkpoint is at sweep **1000**:

| `run` in `.inp` | Behavior |
|-----------------|----------|
| `1100` | ~100 more sweeps (stops at 1100) |
| `10000` | Continues until total sweep 10000 |
| `1000` or less | Exits immediately |

### When is `.chk` written?

The checkpoint is written on the **same schedule as coordinates** (`coord_freq`).
There is **no** automatic `.chk` at the end of a run.

Example: restart from 1000 with `coord_freq 500` and stop at 1100 → **no new
`.chk`** (next write would be at 1500). For long jobs, keep `coord_freq` modest
(e.g. 100 or 500) so you can restart if the queue kills the job mid-run.

### Restart input hygiene

- Use a **new** `# Run_Name` (e.g. `R125equil_cont1.out`) so you do not overwrite
  the original equilibration outputs unless you intend to.
- Prefer a **copy** of the `.inp` (`r125-equil-restart.inp`) rather than editing
  the only copy of the initial setup.

---

## 6. UGE job script and `qsub`

Example script (paths and NetID as used 2026-08):

```bash
#!/bin/bash
#$ -M ed@nd.edu
#$ -m abe
#$ -pe smp 8
#$ -q long
#$ -N r125_equil
#$ -cwd

module load gcc/15.2.0
export OMP_NUM_THREADS=$NSLOTS

cd /users/ed/CassandraV2/projects/hfc125
/users/ed/CassandraV2/Cassandra/Src/cassandra_gfortran_openMP.exe r125-equil-restart.inp
```

Notes:

- `#$ -pe smp N` reserves **N cores on one node** (correct for OpenMP).
- `OMP_NUM_THREADS=$NSLOTS` must match that request (do not oversubscribe).
- Start with 4–8 threads for NPT; scaling may be sublinear — time a short test
  before asking for many cores.
- Do **not** use MPI parallel environments for this OpenMP binary.

```bash
cd /users/ed/CassandraV2/projects/hfc125
qsub r125_equil.job
qstat -u ed
qstat -j <job_id>     # details
qdel <job_id>         # cancel
```

UGE emails (`-m abe`) notify on abort / begin / end if `-M` is set.

---

## 7. Recommended smoke tests before a long `qsub`

**1. Binary + OpenMP (water example)**

```bash
module load gcc/15.2.0
export OMP_NUM_THREADS=4
cd /users/ed/CassandraV2/Cassandra/Examples/NVT/water_spc
/users/ed/CassandraV2/Cassandra/Src/cassandra_gfortran_openMP.exe nvt.inp
```

**2. Tiny checkpoint continuation** (e.g. `run 1100` if chk is at 1000)

Interactive on maginnfe is fine for a short test; then set `run` to the long
total and `qsub`.

---

## 8. Python on maginnfe (interactive shells)

CRC’s `python` module provides **`python3`**, not always a `python` command.
If you see `python: command not found` after `module load python`, call
`python3` (or add the alias below).

Available versions (2026-08): `python/3.12.13`, `python/3.14.5`. Pin a
version so a default change later does not surprise you.

**Recommended interactive `~/.bashrc` on maginnfe** (keep CRC’s Modules init
block at the top; put this under the existing `#Additional aliases` /
`#Additional modules` comments):

```bash
#Additional aliases
alias python=python3

#Additional modules
# Interactive shells only (avoid surprising non-interactive scripts)
if [[ $- == *i* ]]; then
        module load python/3.12.13
fi

export PATH=$HOME/packmol:$PATH
```

After editing:

```bash
source ~/.bashrc
module list          # expect python/3.12.13
which python3
python --version     # expect Python 3.12.13
```

Do **not** put `module load gcc/15.2.0` in `.bashrc` unless you always want
that compiler on every login — keep it in compile steps and job scripts.
Job scripts must still `module load` what they need; do not rely on `.bashrc`
for `qsub`.

Diagnostics need **numpy** (and **matplotlib** for plots). Prefer bringing
results home and using `cassandra-dev` on the laptop (next section). A conda
env on maginnfe for in-cluster plotting is still deferred.

---

## 9. Bring results home for diagnostics

```bash
# on laptop
scp ed@maginnfe.crc.nd.edu:/users/ed/CassandraV2/projects/hfc125/R125equil_cont1.out.prp \
  ~/CassandraV2/projects/hfc125/

conda activate cassandra-dev
cd ~/CassandraV2/Cassandra
python diagnostics/property_plotter.py \
  ~/CassandraV2/projects/hfc125/R125equil_cont1.out.prp -p Mass_Density
```

For bond / angle / dihedral diagnostics, also copy the `.xyz`, species `.mcf`,
and `.inp` when needed:

```bash
python diagnostics/bond_distribution.py \
  path/to/run.out.xyz path/to/species.mcf --bond 1 --save bond1.png --no-show

python diagnostics/angle_distribution.py \
  path/to/run.out.xyz path/to/species.mcf \
  --inp path/to/run.inp --angle 1 --save angle1.png --no-show

python diagnostics/dihedral_distribution.py \
  path/to/run.out.xyz path/to/species.mcf \
  --inp path/to/run.inp --dihedral 1 --save dih1.png --no-show
```

A conda env on maginnfe for plotting is optional later; laptop plotting is enough
for now.

---

## Pitfalls checklist

| Symptom | Likely cause |
|---------|----------------|
| `Permission denied (publickey)` on `git clone` / `git pull` | Key not loaded, missing `~/.ssh/config`, or `.pub` not on GitHub — see §1 |
| `Your branch is up to date` but code looks old | Tracking is stale until a successful `git pull`; fix SSH then pull |
| Dirty `Examples/.../nvt.out.*` blocking a clean tree | Smoke-test leftovers — `git restore Examples/NVT/water_spc/nvt.out.*` |
| `rsync: mkdir ... failed` | Parent `projects/hfc125` missing — `mkdir -p` first |
| Job fails / “No such file” for exe | Used `cassandra_gfortran.exe` instead of **`cassandra_gfortran_openMP.exe`** |
| Run ends instantly after checkpoint | `run` ≤ current sweep in `.chk` |
| No new `.chk` after short restart | Did not land on a `coord_freq` boundary |
| `.mod` / compiler errors on rebuild | Mixed compilers — stick to `module load gcc/15.2.0` |
| Windows-edited job script rejected | Run `dos2unix script.job` |
| `python: command not found` after `module load python` | CRC provides **`python3`** — use it or `alias python=python3` |

---

## What we did for R125 (reference)

1. Equilibrated on laptop → `R125equil.out.chk` at sweep 1000 (density still drifting).
2. Cloned `modernization` on maginnfe; rsync’d `projects/hfc125`.
3. Built OpenMP binary with `gcc/15.2.0`.
4. Smoke-tested water + short restart (`run 1100`).
5. Submitted longer continuation from **original** `.chk` with `run 10000` on
   public `long` (after fixing the OpenMP executable name in the job script).

---

## Deferred / later

- Maginn priority queues and owned-machine host groups  
- `/groups/...` storage if home fills up  
- Conda on maginnfe for in-cluster diagnostics  
- Intel compiler comparison for wall-time tuning  
