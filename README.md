# XRD_simple

A Julia-based simulation tool for powder diffraction patterns of cubic crystal structures. This project generates realistic diffraction patterns for Simple Cubic (SC), Body-Centered Cubic (BCC), and Face-Centered Cubic (FCC) lattices with physics-based modeling of instrumental broadening, crystallite size effects, and microstrain. It supports two radiation modes: **X-ray** (intensity vs 2θ) and **electron** (1D powder profile, intensity vs scattering vector g = 1/d).

## Features

- **Accurate Physics Modeling**
  - Bragg's law for diffraction angles
  - Scherrer equation for crystallite size broadening
  - Caglioti formula for instrumental broadening
  - Voigt and pseudo-Voigt peak profiles
  - Systematic absences for BCC and FCC structures

- **Two Radiation Modes**
  - **X-ray** — intensity vs 2θ (degrees), Cu Kα by default
  - **Electron** — 1D powder profile vs scattering vector g = 1/d (1/Å),
    relativistic wavelength from accelerating voltage, broadening in
    reciprocal-space units (kinematical approximation)

- **Realistic Simulations**
  - Angle-dependent peak broadening
  - Background signal generation
  - Experimental noise simulation
  - Williamson-Hall analysis support

- **Multiple Output Formats**
  - Interactive plots (PNG export)
  - CSV data export
  - Excel spreadsheet export

## Requirements

- Julia version ≥ 1.8
- Dependencies (automatically installed via Project.toml):
  - Plots.jl
  - DataFrames.jl
  - CSV.jl
  - JSON.jl
  - TOML.jl
  - Distributions.jl
  - SpecialFunctions.jl

## Installation

1. Clone this repository:
```bash
git clone <repository-url>
cd XRD_simple
```

2. Install dependencies:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Quick Start

Run the main simulation script:

```julia
include("main.jl")
```

This will generate XRD patterns for all three cubic lattice types (SC, BCC, FCC) and save results to the `results/` directory.

**Note:** `main.jl` auto-detects VS Code (via the `VSCodeServer` module) and
skips interactive pauses there so every figure stays in the plot pane. In a
terminal REPL it pauses between plots so each can be viewed before the next
overwrites it. Use `--no-interactive` for batch or CI runs.

## Configuration

Edit `data.toml` to customize simulation parameters:

```toml
[instrument]
radiation = "xray"           # "xray" | "electron"
two_theta_min = 10.0         # X-ray: minimum 2θ angle (degrees)
two_theta_max = 120.0        # X-ray: maximum 2θ angle (degrees)
lambda = 1.5418              # X-ray: wavelength (Å) — Cu Kα
voltage_kV = 200.0           # electron: accelerating voltage (sets λ ≈ 0.025 Å)
g_min = 0.0                  # electron: min scattering vector (1/Å)
g_max = 1.2                  # electron: max scattering vector (1/Å); d_min ≈ 0.83 Å
N = 1000                     # Number of data points
noise_level = 0.15           # Multiplicative noise level (0–1)

[peak_width]
U = 0.0001                   # X-ray: Caglioti instrumental parameter
V = 0.00005                  # X-ray: Caglioti instrumental parameter
W = 0.00001                  # X-ray: Caglioti instrumental parameter
G_inst = 0.005               # electron: instrumental Gaussian FWHM (1/Å)
K = 0.9                      # Scherrer constant (both)
Epsilon = 0.001              # Microstrain (both)
D = 500.0                    # Crystallite size (nm, both)

[lattice.SC]
Po = 3.352                   # Element = lattice parameter (Å)

[lattice.BCC]
Fe = 2.866
# V  = 3.0399
# W  = 3.155

[lattice.FCC]
Pd = 3.859
# Ag = 4.079
# Cu = 3.594
```

Each uncommented entry under a `[lattice.*]` block produces one pattern.
Leave entries commented out to skip them; add more to run several at once.

**Radiation mode.** `radiation` selects the physics path. With `"xray"` the
X-ray-only keys are used (2θ window, `lambda`, Caglioti U/V/W); with
`"electron"` the electron keys are used (`voltage_kV`, `g_min`/`g_max`,
`G_inst`). The unused keys for the other mode are simply ignored, so both sets
can coexist in one file — flip `radiation` to switch.

## Usage Examples

### Basic Simulation

From a shell:

```bash
julia --project=. main.jl
```

Flags: `--config PATH` (default `data.toml`), `--theme NAME`, `--seed N`,
`--no-interactive`, `--no-plots`.

Or from the REPL / VS Code:

```julia
include("main.jl")
```

### Peak Width Analysis

Explore how peak widths vary with angle:

```julia
include("archive/example_peaks_width.jl")
```

### Voigt vs Pseudo-Voigt Comparison

Compare different peak profile models:

```julia
include("archive/example_use_Voigt.jl")
```

## Output

Running the simulation generates:

- **PNG files**: `results/{element}-{structure}.png` — one per uncommented
  lattice entry (e.g. `Fe-BCC.png`, `Pd-FCC.png`).
- **CSV file**: `results/XRD_results.csv` — an x-axis column plus one intensity
  column per sample, named `{element}-{structure}`. The x column is
  `2θ (deg)` in X-ray mode and `g (1/Å)` in electron mode.

In **electron mode** each sample additionally produces a 2D Debye–Scherrer ring
image (`results/rings/{element}-{structure}.png`) and a reflection answer-key CSV
(`results/rings/{element}-{structure}_reflections.csv`). The ring image is a pure
radial map of the 1D profile: a pixel at radius `r` (mm) takes the intensity at
`g = r / camera_constant`, so ring radius `r = camera_constant · g` and `r² ∝ N`
(`N = h²+k²+l²`). The answer key lists every allowed reflection — `h k l`, `N`,
`g`, ring radius (mm), multiplicity — sorted by `g`; it is the hidden key for the
lab's ring-identification exercise (measure radii → `r²` ratios → `N`-sequence →
SC/BCC/FCC selection rule → lattice constant `a`). Ring cosmetics are tunable in
`[instrument]`: `camera_constant` (λL, mm·Å), `image_px`, `beam_stop_mm`,
`ring_phosphor` (phosphor-green vs grayscale), `ring_gamma`, `ring_noise`.

The final line printed on every run reports how many samples were produced.

**Note:** the `results/` directory is regenerable and is gitignored, so it doesn't
show up in `git status`. Regenerate it any time by re-running the scripts; force-add
a file (`git add -f <path>`) only if you want to snapshot a specific result.

## Physics Background

The simulation implements several key concepts in powder diffraction:

### Bragg's Law
```
nλ = 2d sin(θ)
```
where n is the diffraction order, λ is the wavelength, d is the d-spacing, and θ is the Bragg angle.

### Peak Broadening

**Instrumental Broadening (Gaussian):**
```
β_inst = √(U tan²θ + V tanθ + W)
```

**Size Broadening (Lorentzian):**
```
β_L = Kλ / (L cos θ)
```
where L is the crystallite size and K is the Scherrer constant.

**Strain Broadening (Lorentzian):**
```
β_ε = 4ε tan θ
```
where ε is the microstrain.

For detailed equations and derivations, see [`xrd-peak-broadening.md`](xrd-peak-broadening.md).

### Electron Diffraction (1D powder)

At electron wavelengths (≈ 0.025 Å at 200 kV, from the relativistic de Broglie
relation) every Bragg angle is a fraction of a degree, so a 2θ axis is not
useful. The electron path instead works in **scattering vector** g = 1/d (1/Å),
where reflections sit at purely geometric positions and Bragg's law drops out:
```
g = |G| = √(h² + k² + l²) / a
```
The reflection set is capped by the plotted range (g ≤ `g_max`) rather than the
Bragg `sin θ ≤ 1` bound, which is never binding at electron wavelengths.
Broadening is expressed in reciprocal-space units — size broadening is the
constant Scherrer width Δg = K/D, and strain broadening is Δg = 2ε·g — while the
instrumental term is a single constant Gaussian FWHM `G_inst` (the Caglioti
U/V/W terms are degenerate at θ ≈ 0). The crystallography (Miller indices,
multiplicities, systematic absences) and the Voigt / pseudo-Voigt peak profiles
are shared with the X-ray path.

Peak **heights** are multiplicity-weighted only (the same fidelity as the X-ray
path); the electron atomic scattering factor f_e(s) is not yet modelled, so the
relative intensities are geometric rather than quantitative. The model is
kinematical — valid for thin specimens; real selected-area electron diffraction
is dynamical.

## Project Structure

```
XRD_simple/
├── main.jl                      # Unified entry point (auto-detects VS Code)
├── functions.jl                 # Core physics engine
├── data.toml                    # Configuration file
├── archive/                     # Legacy files and early-stage demo scripts
│   ├── functions_simple.jl      # Simplified reference version (256 lines)
│   ├── simple_XRD.txt           # Legacy text config
│   ├── example_peaks_width.jl   # Peak width demonstration
│   ├── example_use_Voigt.jl     # Voigt profile comparison
│   └── width.jl                 # Peak width analysis utility
├── xrd-peak-broadening.md       # Detailed physics documentation
├── xrd-broadening-references.md # Academic references
└── results/                     # Output directory (gitignored)
```

**Note:** `archive/functions_simple.jl` is a simplified legacy version kept for educational reference. All current scripts use `functions.jl`.

## Documentation

- **[xrd-peak-broadening.md](xrd-peak-broadening.md)** - Comprehensive mathematical background and equations
- **[xrd-broadening-references.md](xrd-broadening-references.md)** - Academic references and foundational papers
- **[problems.md](problems.md)** - Known issues and ongoing investigations

## Known Issues

- Voigt peak widths are approximately 2× broader than pseudo-Voigt profiles with identical input parameters (under investigation)
- Electron mode: peak heights are multiplicity-only — the electron scattering factor f_e(s) is not yet modelled, so relative intensities are geometric rather than quantitative

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## References

Key foundational papers:
- Scherrer, P. (1918) - Crystallite size determination
- Williamson, G. K., & Hall, W. H. (1953) - Size-strain separation
- Caglioti, G., Paoletti, A., & Ricci, F. P. (1958) - Instrumental resolution function

See [`xrd-broadening-references.md`](xrd-broadening-references.md) for the complete reference list.

## License
**[GPL-3.0](LICENSE)**
