# Future fixes / deferred work

Parked ideas that are **not** being done in the current session.
Add new items at the top of the list (newest first). When you pick one up,
move it out (fix or delete) and implement on `modernization` as usual.

**Agent / session rule:** If Maginn says to do something later (defer, park,
“not now,” add to the list, etc.), **do not implement** — add a proper entry
here. See also `.cursor/rules/future-fixes.mdc` and the note in [README.md](README.md).

**How to use**

- One item = one concrete change (or a small cluster that must ship together).
- Prefer enough context that a future you (or an agent) can act without re-reading
  the whole chat: *where*, *symptom*, *proposed fix*, *why deferred*, *risks*.
- Do not put active sprint work here — use commits / PR description for that.

---

## Open items

### Stale comments in the R125 NPT input files

| | |
|--|--|
| **Filed** | 2026-10-16 |
| **Area** | Project inputs, not the Cassandra parser |
| **Files** | `projects/hfc125/NPT-Runs/NPT-{243,258,273,298}/r125-prod.inp` (outside the git repo) |

**What happens today**

- Those files still say “short simulation of pentane” and comments such as
  “run 1000 sweeps”, left over from an example that was copied forward.
- Cassandra does not use that text. `Parse_String` in `Src/io_utilities.f90`
  skips a line whose first character is `!`, and an inline `!` ends a data
  line. `Copy_Inputfile` strips `!` lines from the logfile echo.
- The production runs completed with the numeric settings (500,000 sweeps,
  HFC-125, the pressures in `# Pressure_Info`).

**Why it matters**

- The comments will mislead the next person who opens the input. They are not
  a parser bug.

**Proposed fix (when ready)**

- Edit those four project `.inp` files so the comments match HFC-125 and the
  real run length. No Cassandra source change.

**Status:** deferred — no code change required.

---

### Sphinx / Read the Docs — restart.xyz/H equil→prod workflow

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Public user documentation (Read the Docs) |
| **Files** | `Documentation/source/guides/input_files.rst` (`# Start_Type` / `read_config`); `Documentation/source/getting_started/cassandra_basics.rst` (output files) |
| **Related** | Local notes already on `modernization`: `docs/V2_CRC_WORKFLOW.md` §5, `docs/output-formats.md`, `docs/mcf-setup-workflow.md` |

**What happens today**

- Cassandra writes `*.restart.xyz` + companion `*.restart.H` (latest single
  frame) for equilibration → production handoff; `read_config` aborts if
  companion `.H` disagrees with `# Box_Info` / molecule counts.
- Movie `*.xyz` is multi-frame; `read_config` only uses the **first** frame —
  wrong for production start after equilibration.
- V2 `docs/` describe this; Sphinx sources that feed
  [cassandra-mc.readthedocs.io](https://cassandra-mc.readthedocs.io/) still
  document only generic `.xyz` / `.H` / `.chk` and old `read_config` examples.

**Why it matters**

- Without Sphinx updates, RTD users will not see the recommended workflow after
  a docs rebuild/publish.

**Proposed fix (when ready)** — treat as a dedicated docs pass

1. **Output files** (`cassandra_basics.rst`): document overwrite
   `run_name.restart.xyz` / `.H` (and GEMC `boxN` variants); contrast with
   movie `.xyz` / `.H` and checkpoint `.chk`.
2. **`read_config`** (`input_files.rst`): prefer `*.restart.xyz`; warn not to
   pass the movie; set `# Box_Info` and \(N\) from companion `*.restart.H`;
   document abort-on-mismatch when companion `.H` exists.
3. Short equil → production example; clarify checkpoint = same-run continue,
   restart XYZ/H = fresh production start.

**Status:** deferred — local `docs/` coverage is enough for now; Sphinx/RTD
later as a larger task.

---

### `intra_overlap` API — never cleared to false

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Intramolecular nonbond / CBMC overlap flag |
| **Files** | `Src/energy_routines.f90` — `Compute_Molecule_Nonbond_Intra_Energy` |
| **Related** | Header note (caller must initialize); ring fragment uses a local unused flag |

**What happens today**

- `intra_overlap` is set `.true.` on hard-core hit (early return) but is **never**
  set `.false.` at routine entry.
- Correct use requires every caller to initialize `.false.` first. Some do
  (`Compute_System_Total_Energy`, some moves); some historically may not.

**Proposed fix (when ready)**

1. Set `intra_overlap = .FALSE.` at the start of
   `Compute_Molecule_Nonbond_Intra_Energy` (true `INTENT(OUT)` contract).
2. Audit callers; drop redundant initializations if desired.
3. Update header accordingly.

**Status:** deferred — documented in header; do not change API until picked up.

---

### Improper energy ignores `atom_list%exist` (unlike dihedrals)

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Intramolecular / partial CBMC |
| **Files** | `Src/energy_routines.f90` — `Compute_Molecule_Improper_Energy` vs `Compute_Molecule_Dihedral_Energy` |

**What happens today**

- Dihedral energy skips terms when any of the four atoms has `exist=.FALSE.`
  (partial fragment growth).
- Improper energy has **no** `%exist` guard — assumes a complete molecule.

**Proposed fix (when ready)**

- Mirror dihedral’s exist checks in the improper loop (or document a deliberate
  policy if impropers are never evaluated mid-growth).

**Status:** deferred — noted in improper header.

---

### Species-template self energy (`im` unused; no `%exist`)

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Ewald / DSF self terms |
| **Files** | `Src/energy_routines.f90` — `Compute_Molecule_Self_Energy` |

**What happens today**

- Sums `q^2` over **all** `nonbond_list(1:natoms(is),is)%charge` (MCF template).
- Argument `im` is unused; coordinates / `%exist` ignored.
- Partial ghosts still get full-species self energy in ΔE bookkeeping.

**Proposed fix (when ready)**

- Optionally sum only existing atoms for molecule `im` (or Widom atoms), if
  partial-molecule electrostatics become first-class.
- Until then, keep template behavior; headers now describe it accurately.

**Status:** deferred — clarified in header only.

---

### Ewald reciprocal updates ignore `%exist` (full-species S_mol)

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Ewald reciprocal / CBMC / Widom |
| **Files** | `Update_System_Ewald_Reciprocal_Energy`, `Update_System_Ewald_Reciprocal_Energy_Widom` in `Src/energy_routines.f90` |

**What happens today**

- Real-move update (translation / rotation / intra): loops `ia = 1..natoms(is)`
  with **no** `%exist` filter when building S_mol.
- Widom ghost: uses all species charges × `widom_atoms` coords, no `%exist`.

**Proposed fix (when ready)**

- Filter on `%exist` (and Widom exist flags) when constructing S_mol / S_ghost,
  consistent with pair-energy paths — only if partial-molecule Ewald is required.

**Status:** deferred — noted in both routine headers.

---

### Flexible (harmonic) bonds for multi-fragment / CBMC molecules

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Intramolecular degrees of freedom / CBMC molecule construction |
| **Files** | `Src/input_routines.f90` (`Get_Bond_Info`); `Src/energy_routines.f90` (`Compute_Molecule_Bond_Energy`); fragment / CBMC growth (`Src/participation.f90`, `Src/fragment_growth.f90`, regrow / insert moves); fragment libraries |
| **Related** | Header note under “Practical scope” in `Compute_Molecule_Bond_Energy` |

**What happens today**

- MCF `fixed` bonds are the production path: energy contribution is zero;
  lengths are enforced (integrity check in `Compute_Molecule_Bond_Energy`).
- MCF `harmonic` bonds **are** evaluated in `Compute_Molecule_Bond_Energy`
  (\(E = k(r-r_0)^2\)), but `Get_Bond_Info` **rejects** harmonic bonds unless
  `nfragments(is) == 1` (“Harmonic bonds are only supported for single-fragment
  species”).
- Fragment growth / CBMC place atoms assuming **fixed** bond lengths (e.g.
  points on a sphere of radius \(l_0\) in `participation.f90`). Geometry
  libraries and regrowth inherit that assumption.

**Why it matters**

- Even “vanilla” NVT in Cassandra still uses CBMC-style moves (insertion,
  regrowth, angle/dihedral sampling, etc.) to sample intramolecular degrees
  of freedom. Those paths assume fixed bond lengths when building molecules.
  Without flexible-bond CBMC, **fixed bonds are effectively the only feasible
  production setup** — not merely a preference for multi-fragment species.
- Many force fields use flexible bonds. Cassandra cannot use them for normal
  MC without breaking construction / biasing.
- Enabling harmonic (or other flexible) bonds for CBMC species is **not** a
  local change to the bond-energy routine; the energy evaluator is already
  largely ready for the harmonic form.

**Proposed direction (when ready — large project)**

1. Bond-length sampling in CBMC / fragment growth (not only angles/dihedrals).
2. Correct Rosenbluth weights / Jacobians for variable bond lengths.
3. Regrowth, insertion, and fragment-library generation consistent with
   flexible bonds.
4. Relax or replace the single-fragment restriction in `Get_Bond_Info` once
   construction supports it.
5. Tests: energy bookkeeping, acceptance rates, comparison to fixed-bond
   limits where \(k \to \infty\).

**Risks / decisions**

- Touches core CBMC paths; high regression risk for existing fixed-bond
  molecules (almost all current use cases).
- Scope should be planned as a dedicated modernization epic, not a drive-by
  patch in `energy_routines.f90`.

**Status:** deferred — documented in bond-energy header + this list; do not
implement until explicitly picked up.

---

## Done / cancelled

- **2026-10-16:** The end-of-run logfile block now labels CPU time and prints wall-clock time from `DATE_AND_TIME` (`Main` in `Src/main.f90`). `docs/logfile-format.md` matches. Per-move subroutine timers were left as they were.
- **2026-10-16:** Fixed-angle integrity checks now skip `sim_pregen`, matching fixed bonds (`Compute_Molecule_Angle_Energy`). Normal Monte Carlo still aborts.
