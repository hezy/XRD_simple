# CLAUDE.md - AI Assistant Guide

This file provides context for AI assistants working on the XRD_simple project.

## Project Overview

**XRD_simple** is a Julia-based powder X-ray diffraction (XRD) simulation tool for cubic crystal structures (SC, BCC, FCC). It generates realistic diffraction patterns with physics-based modeling of instrumental broadening, crystallite size effects, and microstrain.

**Primary Use Case:** Educational and research tool for understanding how crystal structure affects XRD patterns.

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
- `functions.jl` - Core physics engine - **THE AUTHORITATIVE VERSION**
- `data.toml` - Configuration file - **THE STANDARD CONFIG FORMAT**

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

### Angle Convention
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
two_theta_min = 10.0         # degrees (auto-converted to radians)
two_theta_max = 120.0        # degrees (auto-converted to radians)
N = 1000                     # number of points
lambda = 1.5418              # wavelength in Angstroms (Cu Kα)
noise_level = 0.15           # multiplicative noise 0–1 (optional)

[peak_width]
U = 0.0001                   # Caglioti parameter (instrumental)
V = 0.00005                  # Caglioti parameter (instrumental)
W = 0.00001                  # Caglioti parameter (instrumental)
K = 0.9                      # Scherrer constant
Epsilon = 0.001              # Microstrain
D = 500.0                    # Crystallite size (nm)

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
entry, any N (including 0) supported.

## Known Issues

### Voigt Width Discrepancy (see problems.md)
- Voigt peaks are ~2× broader than pseudo-Voigt with identical parameters
- Under investigation
- May be related to FWHM calculation vs parameter interpretation

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
- Function: `background()`
- Current model: air scattering (exponential at low angles) + fluorescence
  (constant) + optional amorphous Gaussian hump
- Keep the non-negative intensity constraint.

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

**Last Updated:** 2026-04 (after main.jl / main_VScode.jl merge, archive move,
and multi-lattice config support)
**Maintainer:** Hezy Amiel
**AI Assistant Notes:** Created to provide context for future development sessions
