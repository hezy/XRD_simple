# Functions Comparison: `functions.jl` vs `functions_simple.jl`

This document analyzes the differences between the two function libraries in the XRD_simple project.

## Overview

| Aspect | `functions_simple.jl` | `functions.jl` |
|--------|----------------------|----------------|
| Lines | 256 | 695 |
| Date | April 2023 | 2023-2024 |
| Status | Legacy/educational reference | **Current production** |
| Config format | Plain text parser | TOML |

---

## Key Differences

### 1. Angle Convention

- **functions_simple.jl**: Works entirely in degrees (2θ)
- **functions.jl**: Works in radians internally, converts at I/O boundaries

```julia
# functions_simple.jl:73 - returns degrees
return 2 * (180 / π) * asin.(sinθ_cleaned)

# functions.jl:415 - returns radians
angles = asin.(sinθ[valid_idx])
```

---

### 2. Peak Profile Functions

#### Simple Version
Fixed mixing factor passed as parameter:

```julia
# functions_simple.jl:40-52
function pseudo_Voigt_peak(θ, θ₀, A, w, n)  # n = mixing factor
    return @. A * (n * pdf.(Cauchy...) + (1-n) * pdf.(Normal...))
end
```

#### Full Version
Mixing factor calculated from widths, multiple dispatch, cutoff optimization:

```julia
# functions.jl:168-216 (scalar) and 219-268 (vector)
function pseudo_Voigt_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)
    # Calculates η from Humps2 approximation
    η = 1.36603 * (w_L/w_eff) - 0.47719 * (w_L/w_eff)^2 + 0.11116 * (w_L/w_eff)^3
    # Only calculates within cutoff region for performance
end
```

---

### 3. Peak Width Modeling

#### Simple Version
Single function using Caglioti formula only:

```julia
# functions_simple.jl:56-63
function peaks_width(two_θ_deg, U, V, W)
    return @. √(U * tan(two_θ_deg * π / 360)^2 + V * tan(...) + W)
end
```

#### Full Version
Separates instrumental (Gaussian) from sample (Lorentzian) broadening:

```julia
# functions.jl:344-351 - Gaussian (Caglioti formula for instrumental broadening)
function Gaussian_peaks_width(θ, U, V, W)
    return @. √(U * tan(θ)^2 + V * tan(θ) + W)
end

# functions.jl:358-377 - Lorentzian (Scherrer + Stokes-Wilson for sample broadening)
function Lorentzian_peaks_width(θ, K, E, λ, D)
    w_L_strain = @. 4 * E * tan(θ)      # Microstrain contribution
    w_L_size = @. K * λ / (D * cos(θ))  # Crystallite size contribution
    return @. w_L_strain + w_L_size
end
```

---

### 4. Miller Indices Types

```julia
# functions_simple.jl:79 - uses Int8 (memory efficient but limited)
function d_list(indices::Vector{Vector{Int8}}, a::Float64)

# functions.jl:440 - uses Int (more flexible)
function d_list(indices::Vector{Vector{Int}}, a::Float64)
```

---

### 5. Background Model

#### Simple Version
Parabolic background:

```julia
# functions_simple.jl:161-165
function background(θ)
    return @. 2 + θ * (360 - θ) / 15000
end
```

#### Full Version
Physics-based model with optional noise:

```julia
# functions.jl:577-595
function background(θ; add_noise=false)
    air_scatter = @. 50 * exp(-5θ)   # Exponential decay at low angles
    fluorescence = 10.0               # Constant background
    bg = @. max(air_scatter + fluorescence, 0.0)
    # Optional Gaussian noise can be added
    return bg
end
```

---

### 6. Configuration Parsing

#### Simple Version
Custom text file parser for `simple_XRD.txt`:

```julia
# functions_simple.jl:178-211
function read_file(filename)
    # Parses key-value pairs from plain text format
end
```

#### Full Version
TOML parser with automatic unit conversion for `data.toml`:

```julia
# functions.jl:625-643
function read_xrd_config(filename)
    # Parses TOML format
    # Auto-converts two_theta_min/max from degrees to radians
end
```

---

### 7. Error Handling

#### Simple Version
Minimal error handling:

```julia
# functions_simple.jl:126
error("Invalid cell_type: $cell_type...")
```

#### Full Version
Comprehensive validation with descriptive ArgumentError messages:

```julia
# functions.jl:71-75
A > 0 || throw(ArgumentError("Amplitude A must be positive"))
w_L .> 0 || throw(ArgumentError("Lorentzian width w_L must be positive"))
w_G .> 0 || throw(ArgumentError("Gaussian width w_G must be positive"))
```

---

### 8. Additional Functions in `functions.jl`

| Function | Purpose | Lines |
|----------|---------|-------|
| `estimate_peak_bounds` | Calculate cutoff distance for peak calculation optimization | 271-306 |
| `peak_fwhm` | Combine Gaussian and Lorentzian widths into effective FWHM | 383-410 |
| `bragg_angles` | Returns both angles AND valid indices (not just angles) | 413-436 |

---

## Architectural Differences

### Multiple Dispatch
`functions.jl` provides both scalar and vector versions of peak/width functions:
- Scalar `Voigt_peak`: lines 62-105
- Vector `Voigt_peak`: lines 108-151
- Scalar `pseudo_Voigt_peak`: lines 168-216
- Vector `pseudo_Voigt_peak`: lines 219-268

### Performance Optimization
`functions.jl` uses cutoff regions to skip calculations far from peak centers:

```julia
# Only calculate intensity within cutoff_sigma * w_eff of peak center
cutoff = cutoff_sigma * w_eff
mask = abs.(θ .- θ₀) .< cutoff
```

### Normalization Option
Full version has optional `normalize=true` to scale peak height to 1.0:

```julia
function Voigt_peak(θ, θ₀, A, w_L, w_G; cutoff_sigma=5.0, normalize=false)
    # ...
    if normalize
        I_peak = I_peak ./ maximum(I_peak)
    end
end
```

### Enhanced Return Values
`bragg_angles` in full version returns `Tuple{Vector{Float64}, Vector{Int}}` (angles + valid indices) vs simple version returning just angles.

---

## Function Correspondence Table

| Purpose | `functions_simple.jl` | `functions.jl` |
|---------|----------------------|----------------|
| True Voigt peak | `Voigt_peak(θ, θ₀, A, w, n)` | `Voigt_peak(θ, θ₀, A, w_L, w_G)` |
| Pseudo-Voigt peak | `pseudo_Voigt_peak(θ, θ₀, A, w, n)` | `pseudo_Voigt_peak(θ, θ₀, A, w_L, w_G)` |
| Peak width | `peaks_width(two_θ_deg, U, V, W)` | `Gaussian_peaks_width` + `Lorentzian_peaks_width` + `peak_fwhm` |
| Bragg angles | `bragg_angles(d_list, λ)` | `bragg_angles(d_list, λ)` (enhanced return) |
| d-spacing list | `d_list(indices, a)` | `d_list(indices, a)` |
| Miller indices | `Miller_indices(cell_type, min, max)` | `Miller_indices(cell_type, min, max)` |
| Peak intensities | `peak_intensity(indices)` | `peak_intensity(indices)` |
| Background | `background(θ)` | `background(θ; add_noise=false)` |
| Config parsing | `read_file(filename)` | `read_xrd_config(filename)` |

---

## Summary

| Feature | Simple | Full |
|---------|--------|------|
| Educational clarity | High | Medium |
| Physics accuracy | Basic | Advanced |
| Performance | Baseline | Optimized |
| Error handling | Minimal | Comprehensive |
| Configurability | Low | High |
| Maintenance status | Frozen | Active |

**Bottom Line**: `functions_simple.jl` is a teaching-friendly 256-line implementation suitable for understanding the basics. `functions.jl` is a production-quality 695-line version with proper physics modeling (Scherrer crystallite size, Stokes-Wilson microstrain), performance optimizations, and robust error handling.

---

*Document generated: December 2024*
