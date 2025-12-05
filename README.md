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

**Note:**
- `main.jl` - Interactive version (pauses between plots, press Enter to continue) - recommended for REPL
- `main_VScode.jl` - Non-interactive version (no pauses) - recommended for VS Code or Jupyter

## Configuration

Edit `data.toml` to customize simulation parameters:

```toml
[parameters]
two_theta_min = 2.0          # Minimum 2θ angle (degrees)
two_theta_max = 170.0        # Maximum 2θ angle (degrees)
n_points = 8000              # Number of data points
wavelength = 1.5418          # X-ray wavelength (Å) - Cu Kα

[peak_widths]
U = 0.01                     # Instrumental broadening coefficient
V = -0.02                    # Instrumental broadening coefficient
W = 0.01                     # Instrumental broadening coefficient
K = 0.9                      # Scherrer constant
epsilon = 0.001              # Microstrain
D = 100.0                    # Crystallite size (nm)

[lattice_constants]
SC = 3.0                     # Simple cubic lattice parameter (Å)
BCC = 3.0                    # BCC lattice parameter (Å)
FCC = 3.0                    # FCC lattice parameter (Å)
```

## Usage Examples

### Basic Simulation

```julia
include("functions.jl")

# Read configuration
config = read_xrd_config("data.toml")

# Generate pattern for BCC structure
lattice_type = "BCC"
do_it(lattice_type, config)
```

### Peak Width Analysis

Explore how peak widths vary with angle:

```julia
include("example_peaks_width.jl")
```

### Voigt vs Pseudo-Voigt Comparison

Compare different peak profile models:

```julia
include("example_use_Voigt.jl")
```

## Output

Running the simulation generates:

- **PNG files**: `results/XRD - SC.png`, `results/XRD - BCC.png`, `results/XRD - FCC.png`
- **CSV file**: `results/XRD_results.csv` (tabular data with θ and intensity values)
- **Excel file**: `results/XRD_results.xlsx` (optional, same data)

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
├── main.jl                      # Primary entry point
├── functions.jl                 # Core physics engine (695 lines)
├── data.toml                    # Configuration file
├── example_peaks_width.jl       # Peak width demonstration
├── example_use_Voigt.jl         # Voigt profile comparison
├── width.jl                     # Peak width analysis utility
├── xrd-peak-broadening.md       # Detailed physics documentation
├── xrd-broadening-references.md # Academic references
└── results/                     # Output directory
```

## Documentation

- **[xrd-peak-broadening.md](xrd-peak-broadening.md)** - Comprehensive mathematical background and equations
- **[xrd-broadening-references.md](xrd-broadening-references.md)** - Academic references and foundational papers
- **[problems.md](problems.md)** - Known issues and ongoing investigations

## Known Issues

- Voigt peak widths are approximately 2× broader than pseudo-Voigt profiles with identical input parameters (under investigation)

## Alternative Entry Points

- `main_VScode.jl` - Non-interactive version for VS Code (displays all plots without pausing)

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## References

Key foundational papers:
- Scherrer, P. (1918) - Crystallite size determination
- Williamson, G. K., & Hall, W. H. (1953) - Size-strain separation
- Caglioti, G., Paoletti, A., & Ricci, F. P. (1958) - Instrumental resolution function

See [`xrd-broadening-references.md`](xrd-broadening-references.md) for the complete reference list.

## License

[Add your license information here]
