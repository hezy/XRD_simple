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
- **functions_simple.jl** (256 lines, April 2023) - Legacy educational version
  - Uses text file parsing
  - Fixed mixing factors
  - Works in degrees
  - Basic error handling
  - **STATUS:** Kept as simplified reference, NOT used by active scripts

- **functions.jl** (695 lines, 2023-2024) - **CURRENT PRODUCTION VERSION**
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
- `functions_simple.jl` - Simplified educational reference
- `simple_XRD.txt` - Legacy config (no longer used but kept)

### Documentation
- `README.md` - User-facing documentation
- `xrd-peak-broadening.md` - Physics and equations
- `xrd-broadening-references.md` - Academic citations
- `problems.md` - Known issues

### Examples
- `example_peaks_width.jl` - Demonstrates peak width calculations
- `example_use_Voigt.jl` - Compares Voigt vs pseudo-Voigt
- `width.jl` - Peak width analysis utility

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
`Miller_indices(cell_type, min, max)` implements systematic absences:
- **SC:** All indices allowed (except [0,0,0])
- **BCC:** Only h+k+l = even
- **FCC:** All odd or all even

## Configuration Schema (data.toml)

```toml
[instrument]
two_theta_min = 2.0          # degrees (auto-converted to radians)
two_theta_max = 170.0        # degrees (auto-converted to radians)
N = 8000                     # number of points
lambda = 1.5418              # wavelength in Angstroms (Cu Kα)

[peak_width]
U = 0.01                     # Caglioti parameter (instrumental)
V = -0.02                    # Caglioti parameter (instrumental)
W = 0.01                     # Caglioti parameter (instrumental)
K = 0.9                      # Scherrer constant
Epsilon = 0.001              # Microstrain
D = 100.0                    # Crystallite size (nm)

[lattice.SC]
element = "Generic"          # Material name
a = 3.0                      # Lattice parameter (Angstroms)

[lattice.BCC]
element = "Generic"
a = 3.0

[lattice.FCC]
element = "Generic"
a = 3.0
```

**Important:** Angular parameters in config are in degrees and automatically converted to radians by `read_xrd_config()`.

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
1. Update `Miller_indices()` function in functions.jl (line 520-558)
2. Add systematic absence rules
3. Add lattice parameters to data.toml
4. Update main loop to include new structure

### Modifying Peak Profiles
- Primary location: `Voigt_peak()` and `pseudo_Voigt_peak()` (lines 62-268)
- Remember to update both scalar and vector versions
- Maintain cutoff optimization for performance

### Changing Background Model
- Function: `background()` (lines 577-595)
- Current model: Air scattering (exponential) + fluorescence (constant)
- Keep non-negative intensity constraint

## Testing Approach

No formal test suite currently. Manual testing via:
1. Run `main.jl` and verify plots look reasonable
2. Check example scripts: `example_peaks_width.jl`, `example_use_Voigt.jl`
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
Use pre-allocated arrays where possible (see `d_list()` implementation, line 451).

### Vectorization
Use broadcasting (`@.` macro) for element-wise operations.

## Important Notes for AI Assistants

1. **Always use functions.jl, never functions_simple.jl** for modifications
2. **data.toml is the standard config** - simple_XRD.txt is legacy
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

**Last Updated:** 2024 (following cleanup of duplicate files)
**Maintainer:** Hezy Amiel
**AI Assistant Notes:** Created to provide context for future development sessions
