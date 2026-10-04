# CLAUDE.md - AI Assistant Guide

This file provides context for AI assistants working on the XRD_simple project.

## Project Overview

**XRD_simple** is a Julia-based powder diffraction simulation tool for cubic crystal structures (SC, BCC, FCC). It generates realistic diffraction patterns with physics-based modeling of instrumental broadening, crystallite size effects, and microstrain. It runs in two radiation modes, selected by `radiation` in `data.toml`:
- **X-ray** (default) — intensity vs 2θ (degrees), Bragg geometry.
- **Electron** — 1D powder profile, intensity vs scattering vector g = 1/d (1/Å), reciprocal-space geometry (kinematical approximation).

**Primary Use Case:** Educational and research tool for understanding how crystal structure affects diffraction patterns.

**Technology:** Julia ≥ 1.8, uses Plots.jl, TOML.jl, SpecialFunctions.jl, Distributions.jl

## Project History & Design Decisions

### Configuration Evolution
- **Early 2023:** Used `simple_XRD.txt` (plain text parser)
- **Late 2023-2024:** Migrated to `data.toml` (TOML format) - **THIS IS THE CURRENT STANDARD**
- All active scripts now use `data.toml`

### File Structure Evolution
- **Original:** Multiple main*.jl files for different environments
- **Current:** One unified entry point — `main.jl` — that auto-detects VS Code
  (via `isdefined(Main, :VSCodeServer)`) and skips the between-plots pause +
  `closeall()` when running there. Terminal REPL still pauses for each plot.
- **Rationale:** VS Code's plot pane retains all figures; a terminal REPL
  overwrites each plot and needs a pause to view them.

### Function Library Evolution
- **archive/functions_simple.jl** (April 2023) - Legacy educational version
  - Uses text file parsing
  - Fixed mixing factors
  - Works in degrees
  - Basic error handling
  - **STATUS:** Kept as simplified reference, NOT used by active scripts

- **functions.jl** (2023–2026) - production version until the September 2026
  refactor (`REFACTOR_PLAN.md`), which split it into files under `src/`
  - Uses TOML parsing
  - Advanced physics (Scherrer, Caglioti, Voigt profiles)
  - Works in radians internally
  - Comprehensive error handling
  - Performance optimizations (cutoff regions)

## Critical Files

### Active Scripts (Use These)
- `main.jl` - Unified entry point (auto-detects VS Code vs terminal REPL);
  includes `src/XRDSim.jl`, then `src/plotting.jl`
- `src/XRDSim.jl` - Entry file of the physics (plain includes, not a module):
  loads SpecialFunctions, `Distributions: Normal`, TOML, and includes in order
  - `src/form_factors.jl` - `FORM_FACTOR_COEFFICIENTS` (Waasmaier–Kirfel,
    H–Cf, generated from the DABAX file; do not edit by hand),
    `atomic_form_factor`
  - `src/config.jl` - `Radiation`, `XRay`, `Electron`, `ReflectionModel`,
    `AbsenceRules`, `StructureFactor`, `Atom`, `Sample`, `lattice_sample`,
    `XRDConfig`, `read_xrd_config`
  - `src/crystal.jl` - `cubic_multiplicity`, `family_members`, `Miller_indices`,
    `d_list`, `g_list`, `unit_cell`
  - `src/profiles.jl` - `Voigt_peak`, `pseudo_Voigt_peak`, `peak_fwhm`,
    `Debye_Waller`, `sum_peaks`
  - `src/xray.jl` - Caglioti and Scherrer widths, `Lorentz_polarization`, `bragg_angles`,
    `bragg_max_hkl_sq`, X-ray background, `XRay` methods
  - `src/electron.jl` - `electron_wavelength`, `ed_max_hkl_sq`,
    `Lorentzian_peaks_width_g`, electron background, `Electron` methods,
    `reflection_table`, `ring_image`
  - `src/reflections.jl` - the two reflection methods: `reflections`,
    `scattering_weights` for `AbsenceRules` and `StructureFactor`
  - `src/simulate.jl` - `simulate(cfg, sample)`, and
    `simulate(cfg, structure, element, a, B)` for a monatomic sample
- `src/plotting.jl` - Every Plots.jl call: `plot_title`, `plot_pattern`,
  `plot_ring_image`. Not included by the tests.
- `data.toml` - Configuration file - **THE STANDARD CONFIG FORMAT**

### Analysis / Answer Key (separate private repo)
- The blind indexing & identification tool (`analysis/analyze_results.jl` plus
  `PLAN.md`) was moved **out of this public repo** into a separate **private**
  repo (`git@github.com:hezy/XRD-analysis.git`) so students don't receive the
  answer key. Its history was purged from this repo.
- The `analysis/` directory is **gitignored here** and is its own git repo (the
  private one). It physically stays at `analysis/` so it still reads
  `../results/XRD_results.csv` and `../data.toml` and runs in place. Don't
  re-add it to this repo's tracking.

### Reference/Legacy (Do Not Modify)
- `archive/functions_simple.jl` - Simplified educational reference
- `archive/simple_XRD.txt` - Legacy config (no longer used but kept)

### Documentation
- `README.md` - User-facing documentation
- `xrd-peak-broadening.md` - Physics and equations
- `xrd-broadening-references.md` - Academic citations
- `problems.md` - Known issues

## Key Architecture Patterns

### Radiation Modes
- `read_xrd_config` returns one `XRDConfig`; its field `mode` is an `XRay` or an
  `Electron` (subtypes of `abstract type Radiation`), chosen by `radiation`
  (`"xray"` default, or `"electron"`) and holding the instrument parameters of
  that mode. No other function reads the file or supplies a default.
- The mode is selected by dispatch. One generic `simulate(cfg, sample)`
  returns `(x, y)`, x in display units; its steps are methods on the mode:
  `grid`, `max_hkl_sq`, `peak_centres`, `angular_factor`, `scattering_s`,
  `form_factor`, `peak_widths`, `background`, `display_axis`, `axis_label`
  (and `plot_title` in `plotting.jl`), and on the reflection method
  `cfg.model`: `reflections`, `scattering_weights`.
- **X-ray path:** Bragg geometry, x-axis 2θ (degrees). Uses `lambda`,
  `two_theta_min/max`, Caglioti U/V/W. Peak areas are multiplicity ×
  `Lorentz_polarization(θ)` (unpolarized beam, normalized to 1 at 2θ = 90°) ×
  the scattering weight (see Peak Intensity Weights). Element keys in
  `[lattice.*]` must be chemical symbols (checked in `read_xrd_config`).
- **Electron path:** reciprocal-space geometry, x-axis g = 1/d (1/Å). Positions
  are `g = √(h²+k²+l²)/a` (no Bragg's law); reflection cutoff is `ed_max_hkl_sq`
  (g ≤ `g_max`), not `bragg_max_hkl_sq`. Uses `voltage_kV`, `g_min`/`g_max`,
  `G_inst`. Heights are multiplicity × the scattering weight with f_e
  (`electron_form_factor`, from the X-ray form factor by Mott–Bethe), s = g/2.
- The crystallography (`Miller_indices`, `cubic_multiplicity`, the reflection
  methods) and the peak profiles (`Voigt_peak`, `pseudo_Voigt_peak`,
  `peak_fwhm`, `sum_peaks`) are shared by both paths.

### Angle Convention (X-ray path)
- **Grid and peak centres:** 2θ in **radians**. Widths from Scherrer,
  Stokes–Wilson and Caglioti are FWHM in radians of 2θ, evaluated at θ_B.
- **User input/output:** Degrees (2θ)
- **Conversion:** Done at I/O boundaries (`deg2rad`, `rad2deg`)

### Peak Profile Functions
Two implementations with identical interfaces:
- `Voigt_peak()` - Accurate convolution using complex error function (erfcx)
- `pseudo_Voigt_peak()` - Thompson–Cox–Hastings approximation: Lorentzian and
  Gaussian both with the combined FWHM `peak_fwhm(w_L, w_G)`

Both support:
- Scalar widths only: `sum_peaks` evaluates w_L, w_G once per reflection, at
  its centre, and passes them to the profile
- Cutoff optimization (only calculate near peak center)
- Normalization option
- Error validation

### Peak Width Modeling
- **Gaussian component** (instrumental): `Gaussian_peaks_width()` - Caglioti formula
- **Lorentzian component** (sample): `Lorentzian_peaks_width()` - Scherrer + Stokes-Wilson
- **Effective FWHM:** `peak_fwhm()` combines both
- `peak_widths(mode, x₀, cfg)` returns `(w_L, w_G)` at one peak centre
- **Electron (g-space):** `peak_widths(::Electron, …)` — Lorentzian via
  `Lorentzian_peaks_width_g()` (constant Scherrer K/D + strain 2εg), Gaussian a
  constant `G_inst`. `peak_fwhm()` and the profiles are reused unchanged.

### Peak Intensity Weights
- Each reflection is one pseudo-Voigt peak of area multiplicity ×
  `angular_factor(mode, x₀)` × `scattering_weights(model, mode, sample, indices, s)`,
  s = `scattering_s(mode, x₀)` (sin θ/λ for X-rays, g/2 for electrons).
- **Angular factor:** X-ray `Lorentz_polarization(θ)`, normalized to 1 at
  2θ = 90°; electron 1 (the Lorentz factor is constant in g).
- **Form factor:** `form_factor(mode, element, s)` is `atomic_form_factor` (X-ray)
  or `electron_form_factor` (electron; Mott–Bethe on the Waasmaier–Kirfel
  coefficients, with Z = f₀(0) of the fit, so f_e(0) is finite).
- **Scattering weight, `AbsenceRules`:** (f(s)/f(0))² × `Debye_Waller(s, B)` of
  the single atom. Dividing by f(0) keeps the weights of order 1, so the peaks
  keep their scale relative to the background.
- **Scattering weight, `StructureFactor`:** |F(hkl)|²/F(000)², averaged over
  the family members (`family_members`), F = Σⱼ fⱼ(s) exp(−Bⱼs²)
  exp(2πi h·rⱼ) over `unit_cell(sample)`. Equal to the rules weight for
  monatomic cells (tested in `test/test_structure_factor.jl`).
- B comes from `[debye_waller]` in `data.toml` (per structure and element
  under `[debye_waller.STRUCTURE]`, else `default`, else 0), or from the `B`
  of a basis atom in `[cell.*]`; it is stored in each `Atom` of the sample.

### Reflection Methods
`reflections` in `[model]` selects, once per run, `AbsenceRules` (`"rules"`,
default) or `StructureFactor` (`"structure_factor"`), stored in `cfg.model`.
A `Sample` is a centering (SC/BCC/FCC), `a`, and a basis of `Atom`s;
`unit_cell(sample)` repeats the basis by the centering translations and
rejects two atoms on one site. `[lattice.*]` entries are one-atom samples
(`lattice_sample`), named "element-STRUCTURE"; `[cell.NAME]` sections (keys
`lattice`, `a`, `basis`) are multi-atom samples named NAME and require
`StructureFactor`.
- **`AbsenceRules`:** `Miller_indices(sample.centering, …)`; one-atom basis only.
- **`StructureFactor`:** `Miller_indices("SC", …)`, then a family is dropped
  when the phase sum Σ exp(2πi h·rⱼ) over the sites of each kind of atom
  (element, B) vanishes for every member: a systematic absence, independent
  of s. Weak reflections (KCl odd hkl) stay.

### Miller Index Generation
`Miller_indices(cell_type::String, max_hkl_sq::Int)` enumerates the canonical
`h ≥ k ≥ l ≥ 0` wedge and returns `(indices, multiplicities)`. Systematic
absences:
- **SC:** All indices allowed (except [0,0,0])
- **BCC:** Only h+k+l = even
- **FCC:** All odd or all even

The cutoff `max_hkl_sq` is derived from Bragg physics via `bragg_max_hkl_sq(a, λ)`,
not a hard-coded range. Per-reflection multiplicity comes from `cubic_multiplicity`.

## Configuration Schema (data.toml)

```toml
[instrument]
radiation = "xray"           # "xray" | "electron" (selects the physics path)
two_theta_min = 10.0         # X-ray: degrees (auto-converted to radians)
two_theta_max = 120.0        # X-ray: degrees (auto-converted to radians)
lambda = 1.5418              # X-ray: wavelength in Angstroms (Cu Kα)
voltage_kV = 200.0           # electron: accelerating voltage (sets λ ≈ 0.025 Å)
g_min = 0.0                  # electron: min scattering vector (1/Å)
g_max = 1.2                  # electron: max scattering vector (1/Å)
N = 1000                     # number of points
noise_level = 0.15           # multiplicative noise 0–1 (optional)
camera_constant = 50.0       # electron rings: λL (mm·Å); ring radius r = camera_constant·g
image_px = 800               # electron rings: ring image side length (px)
beam_stop_mm = 2.5           # electron rings: central beam-stop radius (mm)
ring_phosphor = true         # electron rings: phosphor-green colormap (false = grayscale)
ring_gamma = 0.5             # electron rings: display gamma (<1 lifts faint outer rings)
ring_noise = 0.0             # electron rings: per-pixel multiplicative noise (0–1)
ring_ellipticity = 0.0       # electron rings: ellipticity η, semi-axes r(1 ± η)
ring_axis_deg = 0.0          # electron rings: long-axis direction φ₀ (degrees from +x)
ring_centre_x_mm = 0.0       # electron rings: pattern-centre offset (mm)
ring_centre_y_mm = 0.0
ring_radial_distortion = 0.0 # electron rings: barrel (< 0) / pincushion (> 0) κ

[peak_width]
U = 0.0001                   # X-ray: Caglioti parameter (instrumental)
V = 0.00005                  # X-ray: Caglioti parameter (instrumental)
W = 0.00001                  # X-ray: Caglioti parameter (instrumental)
G_inst = 0.005               # electron: instrumental Gaussian FWHM (1/Å)
K = 0.9                      # Scherrer constant (both)
Epsilon = 0.001              # Microstrain (both)
D = 500.0                    # Crystallite size (nm, both)

# Each [lattice.*] block holds one or more element = a (Å) entries.
# Every uncommented line becomes one simulated pattern.
[lattice.SC]
Po = 3.352

[lattice.BCC]
V  = 3.0399
# Fe = 2.866

[lattice.FCC]
Ag = 4.079
# Cu = 3.594

[model]
reflections = "rules"        # "rules" | "structure_factor"

# [cell.NaCl]                # multi-atom cell; needs "structure_factor"
# lattice = "FCC"
# a = 5.640
# basis = [ { element = "Na", xyz = [0, 0, 0] },
#           { element = "Cl", xyz = [0.5, 0, 0], B = 1.2 } ]   # B optional

[debye_waller]               # optional; Debye–Waller B (Å²), both modes
default = 0.0                # elements not listed (0 = no thermal damping)

[debye_waller.BCC]           # B per structure: phases of one element differ
# Fe = 0.3250
```

**Important:** Angular parameters in config are in degrees and automatically
converted to radians by `read_xrd_config()`. Its `XRDConfig` holds, in the
field `samples`, a vector of `Sample`s sorted by centering, then name — one per
uncommented lattice entry or `[cell.*]` section, any N (including 0) supported.
Sample names must be unique. Keys for the unused radiation mode are
ignored, so both X-ray and electron parameters can coexist in one file — flip
`radiation` to switch.

## Known Issues

### Compatibility
- JSON.jl v1.3.0 had compatibility issues with LanguageServer (documented in JSON_compatibility_fix.md)
- Resolved by updating dependencies

## Common Tasks

### Adding a New Crystal Structure
1. Extend `Miller_indices()` with the new `cell_type` branch and its systematic
   absence rule.
2. If its point group differs from cubic, add a new multiplicity helper
   alongside `cubic_multiplicity` and call it from `Miller_indices`.
3. Add its lattice translations to `CENTERING_TRANSLATIONS` (`unit_cell`), so
   that `StructureFactor` gives the same absences.
4. Add a `[lattice.NEWTYPE]` block to `data.toml` with one or more
   `element = a` entries. The main loop picks it up automatically — no
   changes needed in `main.jl`.

A cubic structure with several atoms (NaCl, diamond, ...) needs no code: write
it as a `[cell.NAME]` section and run with `reflections = "structure_factor"`.

### Modifying Peak Profiles
- Primary functions: `Voigt_peak()` and `pseudo_Voigt_peak()`
- Both take scalar widths; keep their FWHM equal to `peak_fwhm` (tested in
  `test/test_peak_profiles.jl`).
- Maintain the cutoff optimization (`cutoff_sigma * w_eff`) for performance.

### Changing Background Model
- Function: `background(::XRay, two_θ)` — air scattering (exponential at low
  angles) + fluorescence (constant) + amorphous Gaussian hump.
- Function: `background(::Electron, g)` — exponential central-beam tail +
  constant inelastic floor, over the g axis.
- Keep the non-negative intensity constraint in both. Noise is applied once,
  in `simulate`.

### Switching / Tuning Radiation Mode
- Set `radiation` in `data.toml` to `"xray"` or `"electron"`.
- Electron path: the `Electron` methods of `simulate`'s steps (in
  `src/electron.jl`); helpers `electron_wavelength()`, `g_list()`,
  `ed_max_hkl_sq()`, `Lorentzian_peaks_width_g()`.
- **Ring output (electron only):** `ring_image(g, y, mode::Electron)` maps the
  1D g-profile to a 2D Debye–Scherrer ring image by radial lookup
  (r = camera_constant·g, so r² ∝ N) and returns `(coords, img)`;
  `plot_ring_image(coords, img, mode)` draws it. `reflection_table(model, mode, sample)`
  returns the discrete answer key (hkl, N, g, multiplicity, scattering weight)
  of the chosen reflection method. `main.jl` calls
  `write_ring_outputs(…)` per electron sample → `results/rings/{title}.png` +
  `{title}_reflections.csv`. Sanity check: `test/ring_sanity.jl` (ring radii vs
  analytic g=√N/a). Knobs: `camera_constant`, `image_px`, `beam_stop_mm`,
  `ring_phosphor`, `ring_gamma`, `ring_noise`.
- **Ring distortion:** `ring_image` looks up each pixel at the undistorted
  radius `ring_true_radius(ρ, φ, mode)`, the inverse of
  ρ = r (1 + η cos 2(φ − φ₀)) (1 + κ (ρ/R)²), R = camera_constant·g_max, with ρ
  measured from the offset pattern centre (where the beam stop sits). Keys
  `ring_ellipticity`, `ring_axis_deg`, `ring_centre_x_mm`, `ring_centre_y_mm`,
  `ring_radial_distortion`, all 0 by default. The answer key keeps the
  undistorted radii. `img[i, j]` is at x = `coords[j]`, y = `coords[i]`
  (`heatmap` convention).
- Electron knobs: `voltage_kV`, `g_min`/`g_max` (detector range), `G_inst`
  (instrumental Gaussian FWHM), plus the shared `K`, `Epsilon`, `D`.

## Testing Approach

A `test/` directory with a runtests.jl harness exists (crystal functions, peak
profiles, widths, pattern computation, background, config, errors, reference
output). Run it with `julia --project=. test/runtests.jl`. The tests include
`src/XRDSim.jl` only, not Plots. `test/test_reference.jl` compares `simulate`
with the saved patterns in `test/reference/`; regenerate them with
`julia --project=. test/reference/generate.jl` only when a change of the
numbers is intended. Manual testing via:
1. Run `main.jl` and verify plots look reasonable
2. Verify CSV output in `results/XRD_results.csv`

**Visual checks:**
- SC: All peaks present
- BCC: Missing peaks follow h+k+l=odd rule
- FCC: Only unmixed parity peaks present
- Peak widths increase with angle (for typical U,V,W values)

## Git Workflow

- **Main branch:** `main`
- **Recent commits:** See `git log --oneline -5` for style
- **Commit style:** Lowercase, descriptive, brief
- Uses Claude Code attribution footer

## Dependencies Management

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()  # Install from Project.toml/Manifest.toml
```

All dependencies tracked in Project.toml. Update with:
```julia
Pkg.add("PackageName")
```

## Physics References

See `xrd-broadening-references.md` for foundational papers:
- Scherrer (1918) - Crystallite size
- Williamson & Hall (1953) - Size-strain separation
- Caglioti et al. (1958) - Instrumental resolution

Equations documented in `xrd-peak-broadening.md`.

## Performance Considerations

### Cutoff Optimization
Peak functions only calculate values within `cutoff_sigma * w_eff` of peak center, where `w_eff` is the combined FWHM. Default: `cutoff_sigma = 5`.

**Trade-off:** Accuracy vs speed. Adjust `cutoff_sigma` parameter if needed.

### Pre-allocation
Use pre-allocated arrays where possible (see `d_list()` for a representative pattern).

### Vectorization
Use broadcasting (`@.` macro) for element-wise operations.

## Important Notes for AI Assistants

1. **Always modify the files in `src/`, never archive/functions_simple.jl**
2. **data.toml is the standard config** - archive/simple_XRD.txt is legacy
3. **Angles:** Internally radians, externally degrees
4. **Widths per reflection:** Evaluate widths at each peak centre, not per grid point
5. **Error handling:** Validate inputs with descriptive ArgumentError messages
6. **Documentation:** Follow existing docstring format (Arguments, Returns, Throws, Examples)
7. **Don't over-engineer:** Keep solutions focused and simple (per project philosophy)
8. **No emojis** in code/documentation unless explicitly requested

## Future Enhancements (Ideas)

The planned work, with its order, is in `improvements.md` (full structure
factor, non-cubic lattices) and
the known defects are in `problems.md`. Further ideas, not planned:

- Implement Rietveld refinement
- Add preferred orientation modeling
- Interactive parameter fitting
- Export to common XRD data formats (XRDML, UXD)

---

**Last Updated:** 2026-10-05 (second reflection method: full structure factor
over the atoms of the unit cell, selected by `[model] reflections`; samples are
`Sample` unit cells, multi-atom cells in `[cell.*]`; 2026-09-28: peak intensities: Lorentz–polarization factor,
Debye–Waller factor with per-element B, X-ray atomic form factor; earlier in
2026-09: refactor of `REFACTOR_PLAN.md`: `functions.jl`
split into `src/`, config parsed once into `XRDConfig`, radiation mode selected
by dispatch, plotting separated from physics; 2026-06: moved the analysis answer key to a separate private
repo and gitignored `analysis/` here; earlier: added electron-diffraction mode —
1D powder profile in g-space, selected by `radiation` in data.toml; main.jl /
main_VScode.jl merge, archive move, multi-lattice config support)
**Maintainer:** Hezy Amiel
**AI Assistant Notes:** Created to provide context for future development sessions
