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

---

## Physics model (open)

All enhancements operate on the intensity of each reflection. Current code:
`I ∝ multiplicity` only. Missing weights:

### 1. Lorentz–polarization factor (highest value, lowest cost)

```
LP(θ) = (1 + cos²(2θ)) / (sin²(θ) · cos(θ))
```

Combined geometric (Lorentz) + polarization terms for unpolarized lab X-rays.
**Changes relative peak heights by 5–10× across a typical scan** — boosts low
angles, suppresses high angles. One-line multiplication inside the peak-sum
loop. Biggest realism gain per line of code.

### 2. Debye–Waller (thermal) factor

```
exp(−2M) = exp(−B · (sin θ / λ)²)
```

Atomic thermal vibration smears scattering; damps high-angle peaks. `B` is
per-element (typical 0.3–1.5 Å² at room temperature). Trivial once the
multiplicative-weight pipeline exists.

### 3. Atomic form factor f(θ)

```
f(sin θ / λ) = Σᵢ aᵢ · exp(−bᵢ · (sin θ / λ)²) + c      [Cromer–Mann, 9 params]
```

X-ray scattering amplitude of a single atom; falls off with angle because the
electron cloud is not a point. Needs a table of Cromer–Mann coefficients per
element. Weight becomes `|F|² ∝ (m · f(θ))²`.

Sources: International Tables Vol. C, §6.1; Waasmaier–Kirfel (1995).

### 4. Full structure factor |F|²

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

### 5. Non-cubic crystal systems

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

### 6. Geometric distortion of the ring image

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

1. **Lorentz–polarization:** biggest realism gain for ~5 lines of code.
2. **Debye–Waller:** natural follow-on; shares the per-family weight hook with LP.
3. **Atomic form factor + full structure factor:** larger project, best done together —
   changes the data model and opens the door to multi-element cells.
4. **Non-cubic lattices:** largest scope; best tackled after the physics model
   is richer, so the per-system modules have complete equations to implement.
5. **Ring-image distortion:** independent of the rest; any time.
