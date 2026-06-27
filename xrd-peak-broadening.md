# X-Ray Diffraction Peak Broadening Equations

## 1. Instrumental Broadening (βinst)
- Generally approximated using a standard reference material (like LaB₆)
- Instrumental Resolution Function (IRF):
```
βinst = U tan²θ + V tanθ + W
```
where U, V, and W are refinable parameters determined from standard measurements

## 2. Crystallite Size Broadening (βL)
- Scherrer Equation:
```
βL = Kλ/(L cosθ)
```
where:
- βL is the peak width due to size effects (in radians)
- K is the Scherrer constant (typically 0.9-1.0)
- λ is the X-ray wavelength
- L is the volume-weighted crystallite size
- θ is the Bragg angle

## 3. Strain Broadening (βε)
- Uniform Strain:
```
βε = 4ε tanθ
```
where ε is the strain parameter

- Wilson Formula for microstrain:
```
βε = 4(⟨ε²⟩)½ tanθ
```
where ⟨ε²⟩ is the mean square strain

## 4. Peak Shape Functions
Common profile functions used to model peak shapes:

### Gaussian Profile:
```
G(x) = (1/(σ√(2π))) exp(-(x-x₀)²/(2σ²))
FWHM = 2√(2ln2)σ ≈ 2.355σ
```

### Lorentzian Profile:
```
L(x) = (1/π) * (γ/2)/((x-x₀)² + (γ/2)²)
FWHM = γ
```

### Pseudo-Voigt Profile (combination of Gaussian and Lorentzian):
```
pV(x) = ηL(x) + (1-η)G(x)
```
where η is the mixing parameter (0 ≤ η ≤ 1)

## 5. Total Peak Broadening

### For Gaussian components:
```
β²total = β²inst + β²size + β²strain
```

### For Lorentzian components:
```
βtotal = βinst + βsize + βstrain
```

### Williamson-Hall Plot Equation:
```
βtotal cosθ = Kλ/L + 4ε sinθ
```
This equation allows separation of size and strain effects by plotting βcosθ vs 4sinθ:
- Slope = strain (ε)
- Y-intercept = Kλ/L (related to crystallite size)

## 6. Integral Breadth Methods
For more accurate analysis using integral breadth (β):

### Warren-Averbach Method:
```
ln A(L) = ln As(L) + ln Ad(L)
```
where:
- A(L) is the Fourier coefficient
- As(L) is the size coefficient
- Ad(L) is the distortion coefficient

### Double-Voigt Method:
```
βL = 1/⟨D⟩v
βG = 4ε tanθ
```
where:
- ⟨D⟩v is the volume-weighted crystallite size
- ε is the upper limit of strain distribution

## 7. Electron Diffraction (reciprocal-space form)

For 1D powder electron diffraction the natural coordinate is the scattering
vector g = 1/d (1/Å) rather than 2θ, because at electron wavelengths
(λ ≈ 0.025 Å) all Bragg angles are sub-degree. The broadening above is then
re-expressed in reciprocal-space units.

### Electron wavelength (relativistic de Broglie)
```
λ = h / √(2 m₀ eV (1 + eV / 2 m₀c²))   ≈   12.2643 / √(V (1 + 0.978476×10⁻⁶ V))   [Å, V in volts]
```
e.g. 200 kV → λ ≈ 0.0251 Å.

### Peak positions
```
g = |G| = √(h² + k² + l²) / a
```
Bragg's law is not needed; positions are purely geometric.

### Size broadening (Lorentzian, in g)
The Scherrer width in reciprocal units is constant in g:
```
Δg_size = K / D
```
(the angle-space β_L = Kλ/(D cosθ) maps to a constant Δ(1/d) under
g = 2 sinθ/λ).

### Strain broadening (Lorentzian, in g)
From Δd/d = ε:
```
Δg_strain = 2 ε g
```

### Instrumental broadening (Gaussian, in g)
A single constant FWHM `G_inst` replaces the Caglioti U tan²θ + V tanθ + W,
whose tanθ terms are degenerate near θ = 0:
```
β_inst = G_inst   (constant)
```

The peak-shape functions (Gaussian, Lorentzian, pseudo-Voigt) and the
size/strain separation are otherwise identical to the X-ray case. Heights are
multiplicity-weighted; the electron atomic scattering factor f_e(s) — obtainable
from the X-ray factor via the Mott–Bethe relation — is not yet included.
