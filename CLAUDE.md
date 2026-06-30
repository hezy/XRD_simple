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

- **functions.jl** (2023–2026) - **CURRENT PRODUCTION VERSION**
  - Uses TOML parsing
  - Advanced physics (Scherrer, Caglioti, Voigt profiles)
  - Works in radians internally
  - Comprehensive error handling
  - Multiple dispatch (scalar/vector versions)
  - Performance optimizations (cutoff regions)

## Critical Files

### Active Scripts (Use These)
- `main.jl` - Unified entry point (auto-detects VS Code vs terminal REPL)
- `functions.jl` - Core physics engine (both X-ray and electron paths) - **THE AUTHORITATIVE VERSION**
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

### Examples (archived)
- `archive/example_peaks_width.jl` - Demonstrates peak width calculations
- `archive/example_use_Voigt.jl` - Compares Voigt vs pseudo-Voigt
- `archive/width.jl` - Peak width analysis utility

## Key Architecture Patterns

### Radiation Modes
- `read_xrd_config` returns `instrument["radiation"]` = `"xray"` (default) or
  `"electron"`. `do_it` and `do_it_zero` branch on it; `do_it_electron` is the
  g-space counterpart of the X-ray `do_it` body.
- **X-ray path:** Bragg geometry, x-axis 2θ (degrees). Uses `lambda`,
  `two_theta_min/max`, Caglioti U/V/W.
- **Electron path:** reciprocal-space geometry, x-axis g = 1/d (1/Å). Positions
  are `g = √(h²+k²+l²)/a` (no Bragg's law); reflection cutoff is `ed_max_hkl_sq`
  (g ≤ `g_max`), not `bragg_max_hkl_sq`. Uses `voltage_kV`, `g_min`/`g_max`,
  `G_inst`. Heights are multiplicity-only (no f_e(s) yet).
- The crystallography (`Miller_indices`, `cubic_multiplicity`, absences) and the
  peak profiles (`Voigt_peak`, `pseudo_Voigt_peak`, `peak_fwhm`, `sum_peaks`) are
  shared by both paths.

### Angle Convention (X-ray path)
- **Internal calculations:** Work in **radians** (θ, not 2θ)
- **User input/output:** Degrees (2θ)
- **Conversion:** Done at I/O boundaries (`deg2rad`, `rad2deg`)

### Peak Profile Functions
Two implementations with identical interfaces:
- `Voigt_peak()` - Accurate convolution using complex error function (erfcx)
- `pseudo_Voigt_peak()` - Fast linear approximation

Both support:
- Scalar and vector width parameters (multiple dispatch)
- Cutoff optimization (only calculate near peak center)
- Normalization option
- Error validation

### Peak Width Modeling
- **Gaussian component** (instrumental): `Gaussian_peaks_width()` - Caglioti formula
- **Lorentzian component** (sample): `Lorentzian_peaks_width()` - Scherrer + Stokes-Wilson
- **Effective FWHM:** `peak_fwhm()` combines both
- **Electron (g-space):** `compute_peak_widths_g()` — Lorentzian via
  `Lorentzian_peaks_width_g()` (constant Scherrer K/D + strain 2εg), Gaussian a
  constant `G_inst`. `peak_fwhm()` and the profiles are reused unchanged.

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
```

**Important:** Angular parameters in config are in degrees and automatically
converted to radians by `read_xrd_config()`. That function returns a flat
vector of `(structure, element, a)` triples — one per uncommented lattice
entry, any N (including 0) supported. Keys for the unused radiation mode are
ignored, so both X-ray and electron parameters can coexist in one file — flip
`radiation` to switch.

## Known Issues

### Voigt Width Discrepancy (see problems.md)
- Voigt peaks are ~2× broader than pseudo-Voigt with identical parameters
- Under investigation
- May be related to FWHM calculation vs parameter interpretation

### Electron Intensities (multiplicity-only)
- Electron-mode peak heights are weighted by multiplicity only; the electron
  scattering factor f_e(s) is not modelled, so relative intensities are
  geometric, not quantitative. Adding f_e(s) (Doyle–Turner/Kirkland, or
  Mott–Bethe on X-ray f_x) would make low-g reflections correctly dominant —
  and would also upgrade the X-ray heights, which are likewise multiplicity-only.

### Compatibility
- JSON.jl v1.3.0 had compatibility issues with LanguageServer (documented in JSON_compatibility_fix.md)
- Resolved by updating dependencies

## Common Tasks

### Adding a New Crystal Structure
1. Extend `Miller_indices()` with the new `cell_type` branch and its systematic
   absence rule.
2. If its point group differs from cubic, add a new multiplicity helper
   alongside `cubic_multiplicity` and call it from `Miller_indices`.
3. Add a `[lattice.NEWTYPE]` block to `data.toml` with one or more
   `element = a` entries. The main loop picks it up automatically — no
   changes needed in `main.jl`.

### Modifying Peak Profiles
- Primary functions: `Voigt_peak()` and `pseudo_Voigt_peak()`
- Each has both a scalar and a vector method — update both when changing
  behavior.
- Maintain the cutoff optimization (`cutoff_sigma * w_eff`) for performance.

### Changing Background Model
- Function: `background()` (X-ray) — air scattering (exponential at low angles)
  + fluorescence (constant) + optional amorphous Gaussian hump.
- Function: `background_electron()` (electron) — exponential central-beam tail +
  constant inelastic floor, over the g axis.
- Keep the non-negative intensity constraint in both.

### Switching / Tuning Radiation Mode
- Set `radiation` in `data.toml` to `"xray"` or `"electron"`.
- Electron path entry point: `do_it_electron()`; helpers `electron_wavelength()`,
  `g_list()`, `ed_max_hkl_sq()`, `Lorentzian_peaks_width_g()`,
  `compute_peak_widths_g()`, `intensity_vs_g()`, `compute_ed_pattern()`.
- **Ring output (electron only):** `render_ring_image(g, y, camera_constant; …)`
  maps the 1D g-profile to a 2D Debye–Scherrer ring image by radial lookup
  (r = camera_constant·g, so r² ∝ N); `reflection_table(structure, a, g_max)`
  returns the discrete answer key (hkl, N, g, multiplicity). `main.jl` calls
  `write_ring_outputs(…)` per electron sample → `results/rings/{title}.png` +
  `{title}_reflections.csv`. Sanity check: `test/ring_sanity.jl` (ring radii vs
  analytic g=√N/a). Knobs: `camera_constant`, `image_px`, `beam_stop_mm`,
  `ring_phosphor`, `ring_gamma`, `ring_noise`.
- Electron knobs: `voltage_kV`, `g_min`/`g_max` (detector range), `G_inst`
  (instrumental Gaussian FWHM), plus the shared `K`, `Epsilon`, `D`.
- To add the electron scattering factor f_e(s), weight each reflection in
  `intensity_vs_g`/`sum_peaks` by m·|F|² instead of m (see Known Issues).

## Testing Approach

A `test/` directory with a runtests.jl harness exists (crystal functions, peak
profiles, widths, pattern computation, background, config, errors). Manual
testing via:
1. Run `main.jl` and verify plots look reasonable
2. Optionally check archived demo scripts: `archive/example_peaks_width.jl`,
   `archive/example_use_Voigt.jl`
3. Verify CSV output in `results/XRD_results.csv`

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
Peak functions only calculate values within `cutoff_sigma * w_eff` of peak center. Default: 5σ.

**Trade-off:** Accuracy vs speed. Adjust `cutoff_sigma` parameter if needed.

### Pre-allocation
Use pre-allocated arrays where possible (see `d_list()` for a representative pattern).

### Vectorization
Use broadcasting (`@.` macro) for element-wise operations.

## Important Notes for AI Assistants

1. **Always use functions.jl, never archive/functions_simple.jl** for modifications
2. **data.toml is the standard config** - archive/simple_XRD.txt is legacy
3. **Angles:** Internally radians, externally degrees
4. **Multiple dispatch:** Maintain both scalar and vector versions of width functions
5. **Error handling:** Validate inputs with descriptive ArgumentError messages
6. **Documentation:** Follow existing docstring format (Arguments, Returns, Throws, Examples)
7. **Don't over-engineer:** Keep solutions focused and simple (per project philosophy)
8. **No emojis** in code/documentation unless explicitly requested

## Future Enhancements (Ideas)

- Add hexagonal crystal structures
- Implement Rietveld refinement
- Add preferred orientation modeling
- Create formal test suite
- Interactive parameter fitting
- Export to common XRD data formats (XRDML, UXD)

---

**Last Updated:** 2026-06 (moved the analysis answer key to a separate private
repo and gitignored `analysis/` here; earlier: added electron-diffraction mode —
1D powder profile in g-space, selected by `radiation` in data.toml; main.jl /
main_VScode.jl merge, archive move, multi-lattice config support)
**Maintainer:** Hezy Amiel
**AI Assistant Notes:** Created to provide context for future development sessions
