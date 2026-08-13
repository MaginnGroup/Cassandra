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

### Align fixed-angle integrity check with fixed bonds for `sim_pregen`

| | |
|--|--|
| **Filed** | 2026-08-13 |
| **Area** | Intramolecular energy / pregenerated trajectory analysis |
| **Files** | `Src/energy_routines.f90` — `Compute_Molecule_Angle_Energy`, `Compute_Molecule_Bond_Energy` |
| **Related** | `Src/pregen_driver.f90` (calls `Compute_System_Total_Energy` when `need_energy`); trajectory/Widom-from-LAMMPS path (`sim_pregen`) |

**What happens today**

- Both routines treat MCF type `fixed` (`int_none`) as **zero energy** plus a
  geometry integrity check against MCF value ± tolerance.
- Fixed **bonds** skip that abort when `int_sim_type == sim_pregen`:

  ```fortran
  IF ((int_sim_type /= sim_pregen) .AND. (abs(l0 - length) > ltol)) THEN
  ```

- Fixed **angles** always abort if `|theta0 - theta| > theta_tol` (degrees).
  There is no `sim_pregen` gate.

**Why it matters**

- Pregen is meant to score / analyze existing frames (e.g. LAMMPS → XYZ/H →
  Cassandra Widom). Those frames often do not match MCF fixed lengths/angles
  within tolerance.
- Bonds were gated for that workflow in PR #106 (commit `81acea75`, 2022:
  “Add trajectory reader and enhance Widom insertion support”). Angles were
  left ungated — almost certainly an oversight, not a second check elsewhere.
- Fixed bonds/angles are **not** re-validated in another routine; fragment
  growth uses fixed values for *placement*, not for checking loaded coords.

**Proposed fix (when ready)**

1. Mirror the bond condition in `Compute_Molecule_Angle_Energy` for the fixed
   (`int_none`) branch.
2. Update the angle routine header note that currently says the check is
   *not* skipped for `sim_pregen`.
3. Optional: one-line comments on both bond and angle checks explaining that
   pregen skips rigid-geometry aborts so foreign trajectories can be analyzed.
4. Smoke / mental test: pregen + `need_energy` with a molecule that has fixed
   angles slightly off MCF should complete; normal MC should still abort.

**Risks / decisions**

- Aligning weakens pregen as a rigidity sanity check for angles (bonds already
  weak). If someone relied on angle aborts during pregen, they lose that signal.
- Alternative (probably worse): remove the bond skip — would break the intended
  LAMMPS-trajectory use case.

**Status:** deferred — documented only; do not implement until explicitly picked up.

---

## Done / cancelled

_(Move completed items here with a one-line note and date, or delete them.)_
