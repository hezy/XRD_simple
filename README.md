# XRD_simple

A Julia-based simulation tool for powder X-ray diffraction (XRD) patterns of cubic crystal structures. This project generates realistic diffraction patterns for Simple Cubic (SC), Body-Centered Cubic (BCC), and Face-Centered Cubic (FCC) lattices with physics-based modeling of instrumental broadening, crystallite size effects, and microstrain.

## Features

- **Accurate Physics Modeling**
  - Bragg's law for diffraction angles
  - Scherrer equation for crystallite size broadening
  - Caglioti formula for instrumental broadening
  - Voigt and pseudo-Voigt peak profiles
  - Systematic absences for BCC and FCC structures

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
two_theta_min = 10.0         # Minimum 2θ angle (degrees)
two_theta_max = 120.0        # Maximum 2θ angle (degrees)
N = 1000                     # Number of data points
lambda = 1.5418              # X-ray wavelength (Å) — Cu Kα
noise_level = 0.15           # Multiplicative noise level (0–1)

[peak_width]
U = 0.0001                   # Caglioti instrumental parameter
V = 0.00005                  # Caglioti instrumental parameter
W = 0.00001                  # Caglioti instrumental parameter
K = 0.9                      # Scherrer constant
Epsilon = 0.001              # Microstrain
D = 500.0                    # Crystallite size (nm)

[lattice.SC]
Po = 3.352                   # Element = lattice parameter (Å)

[lattice.BCC]
# Fe = 2.866
V  = 3.0399
# W  = 3.155

[lattice.FCC]
Ag = 4.079
# Cu = 3.594
# Au = 4.065
```

Each uncommented entry under a `[lattice.*]` block produces one XRD pattern.
Leave entries commented out to skip them; add more to run several at once.

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
  lattice entry (e.g. `V-BCC.png`, `Ag-FCC.png`).
- **CSV file**: `results/XRD_results.csv` — a `θ` column plus one intensity
  column per sample, named `{element}-{structure}`.

The final line printed on every run reports how many samples were produced.

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
└── results/                     # Output directory
```

**Note:** `archive/functions_simple.jl` is a simplified legacy version kept for educational reference. All current scripts use `functions.jl`.

## Documentation

- **[xrd-peak-broadening.md](xrd-peak-broadening.md)** - Comprehensive mathematical background and equations
- **[xrd-broadening-references.md](xrd-broadening-references.md)** - Academic references and foundational papers
- **[problems.md](problems.md)** - Known issues and ongoing investigations

## Known Issues

- Voigt peak widths are approximately 2× broader than pseudo-Voigt profiles with identical input parameters (under investigation)

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
