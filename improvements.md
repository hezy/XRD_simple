# Improvement backlog

Open suggestions for XRD_simple, collected during the April 2026 code review
and refactor. Grouped by category and ordered roughly by value-per-effort
within each section.

## Current state (2026-10-05)

- **Done:** the full structure factor as a second reflection method (commits
  808351d, 5ea05af, 9fa8210, d09a036), and the public `data.toml` updated with
  the ring, Debye–Waller, `[model]` and example-cell sections (b97aaa0).
- **Not done:** non-cubic lattices (item 1 below); the private analysis repo
  cannot read cell samples (item 2 below).
- **Next step:** write a plan for non-cubic lattices, starting with which
  crystal systems the course needs.

---

## Already completed

Listed for context; no further action needed.

- Removed dead `sinθ_cleaned` in `bragg_angles`
- Removed misleading broadcast dots in scalar `Voigt_peak` / `pseudo_Voigt_peak` validation
- Fixed "Scaler" → "Scalar" typos
- Removed the orphan `abstract_peak` docstring; its argument list is now in
  the `Voigt_peak` docstring
- Refactored `Miller_indices` to return `(indices, multiplicities)` with canonical `h ≥ k ≥ l ≥ 0` enumeration
- Added `cubic_multiplicity` helper with cubic point-group counts
- Added `bragg_max_hkl_sq(a, λ)` — physics-driven cutoff from `sin(θ) ≤ 1`
- Removed hardcoded `MILLER_INDEX_MIN/MAX = ±5` constants
- Output matches pre-refactor to floating-point epsilon
- `read_xrd_config` now returns every uncommented `[lattice.*]` entry as a
  `(structure, element, a)` triple — no more silent overwrite when multiple
  elements are uncommented in one block. Any N samples, including 0, is valid.
- Merged `main_VScode.jl` into `main.jl`; VS Code is auto-detected via
  `isdefined(Main, :VSCodeServer)` and the between-plot pause / `closeall()`
  are skipped there.
- Archived five early-stage / legacy files (`functions_simple.jl`,
  `simple_XRD.txt`, `example_peaks_width.jl`, `example_use_Voigt.jl`,
  `width.jl`) into `archive/`. The three example scripts and the review notes
  on `functions.jl` were later deleted (September 2026); they no longer ran
  after the refactor.
- September 2026 refactor (`REFACTOR_PLAN.md`), which also closed the former
  code-cleanup items: `using Distributions: Normal`; the vector methods of the
  peak functions (and with them the `peak_fwhm` scalar+vector MethodError) were
  removed; `intensity_vs_angle` and `do_it_zero` were removed; the repeated
  checks in `reflection_table` were removed.
- Fixed the Voigt vs. pseudo-Voigt width discrepancy: `pseudo_Voigt_peak` now
  uses the combined FWHM for both components (refactor Phase 1).
- Added the Lorentz–polarization factor to X-ray peak areas
  (`Lorentz_polarization`, normalized to 1 at 2θ = 90°), through the mode step
  `peak_weights`; electron weights are 1 (September 2026).
- Added the Debye–Waller factor exp(−2B (sin θ/λ)²) in both modes, with B per
  element from the optional `[debye_waller]` section of `data.toml`
  (September 2026).
- Added the X-ray atomic form factor: Waasmaier–Kirfel coefficients for H–Cf
  in `src/form_factors.jl`; peak areas carry (f/Z)² (September 2026).
- Added the electron scattering factor: `electron_form_factor` by the
  Mott–Bethe relation on the Waasmaier–Kirfel coefficients; electron peak
  weights carry (f_e/f_e(0))² (October 2026).
- Added geometric distortion of the electron ring image: ellipticity and its
  axis, pattern-centre offset, and barrel/pincushion distortion, as optional
  `[instrument]` keys with default 0 (October 2026).
- Added the full structure factor as a second reflection method, selected by
  `reflections` in `[model]` (`AbsenceRules` or `StructureFactor`): F(hkl)
  over the atoms of the unit cell, with multi-atom cells given as
  `[cell.NAME]` sections (lattice, a, basis). For monatomic cells both methods
  give the same pattern (October 2026).
- B is keyed by structure, `[debye_waller.STRUCTURE]`, since phases of one
  element differ (BCC and FCC Fe). `data.toml` holds the 293 K values of
  Peng, Ren, Dudarev & Whelan (1996), supplement SUP82472, Table 1, for the
  tabulated elements of the lattice menu; the others use `default` (October
  2026).

---

## Architecture / extensibility (open)

### 1. Non-cubic crystal systems

Current code is cubic-only. Generalizing to tetragonal, hexagonal,
orthorhombic, etc. needs:

- Per-system d-spacing formula (today: `1/d² = (h²+k²+l²)/a²`)
- Per-system multiplicity (today: `cubic_multiplicity`)
- Canonical enumeration order (today: `h ≥ k ≥ l ≥ 0`; tetragonal would
  fix the c-axis separately, etc.)
- Config schema accepting multiple lattice parameters (a, c for tetragonal;
  a, b, c, α, β, γ for triclinic)
- Extended systematic-absence rules (screw axes, glide planes)

The current refactor was designed to generalize here cleanly: the pipeline
shape (per-family multiplicity weighting, physics-driven cutoff) is unchanged;
only the inner helpers become per-system.

The structure-factor method adds two cubic-only helpers that must also become
per-system: `family_members` (the m-3m sign and permutation variants; a
non-cubic family has the members of its Laue group) and
`CENTERING_TRANSLATIONS` in `unit_cell` (only P, I, F today; C and R
centerings are needed).

### 2. Analysis tool cannot read cell samples (private repo)

`analysis/analyze_results.jl` (private repo `hezy/XRD-analysis`) takes the
sample names in `results/XRD_results.csv` to be `element-STRUCTURE`. A
`[cell.NAME]` sample is named `NAME` (e.g. `NaCl`), so the tool cannot parse
it, and its SC/BCC/FCC identification does not cover multi-atom cells.

---

## Suggested order

1. **Non-cubic lattices:** the main remaining item; the structure factor, the
   form factors and the Debye–Waller factor are in place, and need only the
   per-system d-spacing, multiplicity, enumeration and family members.
2. **Analysis tool and cell samples:** only when multi-atom cells are given
   to students as unknowns.
