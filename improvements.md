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
- Corrected `abstract_peak` docstring type to `Union{Float64, Vector{Float64}}`
- Refactored `Miller_indices` to return `(indices, multiplicities)` with canonical `h ≥ k ≥ l ≥ 0` enumeration
- Added `cubic_multiplicity` helper with cubic point-group counts
- Added `bragg_max_hkl_sq(a, λ)` — physics-driven cutoff from `sin(θ) ≤ 1`
- Removed hardcoded `MILLER_INDEX_MIN/MAX = ±5` constants
- Output matches pre-refactor to floating-point epsilon

---

## Code cleanup (open)

| # | Item | Effort | Impact |
|---|---|---|---|
| 1 | Scope `using Distributions` → `using Distributions: Normal` | trivial | minor (faster load, clearer deps) |
| 2 | `peak_fwhm` on mixed scalar+vector args gives a bare MethodError | trivial | low (diagnostic only) |
| 3 | Redundant λ/a validation in `intensity_vs_angle` (callees re-validate) | trivial | low (code clarity) |
| 4 | `do_it_zero` reads full config just to grab instrument params | trivial | low (minor waste) |
| 5 | `read_xrd_config` silently keeps only the last uncommented element per lattice block | small | **medium** — silent data loss if user forgets to comment one out |

**5** is the one with real bite. Current code (functions.jl:657-663):

```julia
for (structure, elements) in config["lattice"]
    if !isempty(elements)
        for (element, value) in elements
            lattice[structure] = (element, value)  # overwrites silently
        end
    end
end
```

Fix options: error on duplicates, document the single-element invariant, or
accept a list and let the caller choose.

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

## Known issues

### Voigt vs. pseudo-Voigt width discrepancy

Documented in `problems.md`: Voigt peaks are ~2× broader than pseudo-Voigt
with identical `w_L` and `w_G`. Not a cleanup — a real physics/implementation
bug. Likely causes:

- FWHM vs. HWHM mismatch in the erfcx argument
- Sigma derivation: `σ = w_G / (2√(2 ln 2))` is Gaussian-standard but needs
  verifying against the complex-error-function normalization
- Cutoff region using `w_eff` from `peak_fwhm` may mask the true profile

Worth a separate investigation session — not a quick fix.

---

## Suggested order

1. **Quick cleanup commit:** items 1, 3, and 5 (the silent-overwrite) together.
   Small surface, meaningful safety improvement.
2. **Lorentz–polarization:** biggest realism gain for ~5 lines of code.
3. **Debye–Waller:** natural follow-on; shares the per-family weight hook with LP.
4. **Atomic form factor + full structure factor:** larger project, best done together —
   changes the data model and opens the door to multi-element cells.
5. **Non-cubic lattices:** largest scope; best tackled after the physics model
   is richer, so the per-system modules have complete equations to implement.
6. **Voigt width bug:** separate track, whenever the user wants to diagnose.
