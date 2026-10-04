# Improvement backlog

Open suggestions for XRD_simple, collected during the April 2026 code review
and refactor. Grouped by category and ordered roughly by value-per-effort
within each section.

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
- B is keyed by structure, `[debye_waller.STRUCTURE]`, since phases of one
  element differ (BCC and FCC Fe). `data.toml` holds the 293 K values of
  Peng, Ren, Dudarev & Whelan (1996), supplement SUP82472, Table 1, for the
  tabulated elements of the lattice menu; the others use `default` (October
  2026).

---

## Physics model (open)

All enhancements operate on the intensity of each reflection. Current code:
`I ∝ multiplicity × LP(θ) × (f/Z)² × DW` for X-rays, `I ∝ multiplicity × DW`
for electrons.
Missing weights:

### 1. Electron scattering factor f_e(s) (Mott–Bethe)

```
f_e(s) = 0.023934 Å · (Z − f(s)) / s²,   s = sin θ / λ = g / 2
```

Electron-mode heights carry no atomic scattering factor, so they are not
quantitative (see `problems.md`). The Mott–Bethe relation gives f_e from the
X-ray `atomic_form_factor` already in `src/form_factors.jl`; no new data table
is needed. Multiply (f_e(s)/f_e(s_ref))², or another normalization of order 1,
into `peak_weights(::Electron, g₀, element, B)`. The relation is singular at
s = 0, but every reflection has s > 0. Optional: the relativistic factor γ at
`voltage_kV`, which scales all f_e equally and does not change relative
heights.

### 2. Full structure factor |F|²

```
F(hkl) = Σⱼ fⱼ · exp(2πi · (h xⱼ + k yⱼ + l zⱼ))
```

Replaces the hardcoded "BCC means h+k+l even" rule with a general sum over
atom positions in the unit cell. Collapses to the current centering filters
for single-element cubic, generalizes to multi-element cells (NaCl, diamond,
perovskites, alloys). Requires a data model change: unit cell = list of
(element, fractional position).

---

## Architecture / extensibility (open)

### 3. Non-cubic crystal systems

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

---

## Electron ring image (open)

### 4. Geometric distortion of the ring image

`ring_image` now draws ideal, exactly circular rings centred in the frame.
Real SAED patterns are distorted, and students must measure through that:

- Elliptical distortion (projector-lens astigmatism): the radius depends on
  azimuth, r(φ) = r₀ · (1 + η cos 2(φ − φ₀)), with ellipticity η of about
  0.5–2 % and axis angle φ₀
- Pattern centre offset from the image centre (beam-stop position)
- Optional: barrel or pincushion distortion, a radius error that grows with r

Implementation: in the radial lookup of `ring_image`, replace r by the
corrected radius before mapping to g = r / `camera_constant`. New `[instrument]`
keys (e.g. `ring_ellipticity`, `ring_axis_deg`, `ring_centre_mm`), default 0.

Tests: the ring-radius check of the `ring_image` test set assumes circular
rings (the (100) ring within one pixel in four directions). With distortion,
check the azimuthally averaged radius instead, or allow a tolerance equal to
the distortion.

---

## Suggested order

1. **Electron scattering factor:** a few lines on top of the X-ray form factor;
   removes the open item of `problems.md`.
2. **Full structure factor:** changes the data model and opens the door to
   multi-element cells; the atomic form factor it needs is in place.
3. **Non-cubic lattices:** largest scope; best tackled after the physics model
   is richer, so the per-system modules have complete equations to implement.
4. **Ring-image distortion:** independent of the rest; any time.
