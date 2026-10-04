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

## 7. Peak Intensity

The integrated intensity of reflection hkl is its multiplicity m times the
Lorentz–polarization factor (X-ray, unpolarized beam without monochromator),
the squared atomic form factor (X-ray) and the Debye–Waller factor:
```
I_hkl ∝ m · LP(θ) · (f(s)/Z)² · exp(−2B s²)
LP(θ) = (1 + cos²2θ) / (sin²θ · cos θ),   s = sin θ / λ = 1/(2d)
f(s)  = c + Σᵢ aᵢ exp(−bᵢ s²),   i = 1…5
```
LP is normalized to 1 at 2θ = 90°. It is large at low angles and has its
minimum near 2θ ≈ 100°–120°. B = 8π²⟨u²⟩ (Å²) is the parameter of the atomic
temperature factor exp(−B s²) on the amplitude, so the intensity carries
exp(−2B s²); typical room-temperature values are 0.2–2 Å².

f(s) is the X-ray scattering amplitude of one neutral atom, in electrons. It
equals Z at s = 0 and falls with s because the electron cloud has a finite
size. The coefficients are those of Waasmaier & Kirfel (1995), Acta Cryst.
A51, 416–431, fitted to International Tables Vol. C, Table 6.1.1.1, and valid
for 0 ≤ s ≤ 6 1/Å. Dividing by Z keeps the weights of order 1; for a
single-element cubic cell the structure factor is n·f (n atoms per cell), so
the relative intensities within one pattern are unchanged by this scaling.
Anomalous dispersion (f′, f″) is omitted.

In electron mode s = g/2, the Lorentz–polarization factor is omitted, and f
is replaced by the electron scattering factor f_e (Å), obtained from f by the
Mott–Bethe relation:
```
I_hkl ∝ m · (f_e(s)/f_e(0))² · exp(−2B s²)
f_e(s) = C (Z − f(s)) / s² = C Σᵢ aᵢ (1 − exp(−bᵢ s²)) / s²,   C = 0.023934 1/Å
f_e(0) = C Σᵢ aᵢ bᵢ
```
The second form takes Z = f(0) of the fit, which removes the 0/0 at s = 0.
f_e is the scattering of the electrostatic potential: the nucleus minus the
electron cloud. It falls much faster with s than f, so low-g rings dominate.
The relativistic factor γ scales every f_e equally and is omitted.

### Structure factor

With `reflections = "structure_factor"` the factor (f/Z)² (or (f_e/f_e(0))²)
and the Debye–Waller factor are replaced by the structure factor of the unit
cell:
```
I_hkl ∝ m · LP(θ) · |F(hkl)|² / F(000)²
F(hkl) = Σⱼ fⱼ(s) exp(−Bⱼ s²) exp(2πi (h xⱼ + k yⱼ + l zⱼ))
F(000) = Σⱼ fⱼ(0)
```
The sum runs over every atom j of the cubic cell, at fractional position
(xⱼ, yⱼ, zⱼ): the basis atoms and their copies by the centering translations,
(½,½,½) for BCC and (0,½,½), (½,0,½), (½,½,0) for FCC. Each atom carries its own
B. |F|² is averaged over the m members of the family, which changes nothing
for a cell with the full cubic symmetry m-3m. A reflection is systematically
absent when the phase sum over the sites of each kind of atom vanishes; it is
then zero at every s.

For a monatomic cell of n atoms F = n f exp(−B s²) on the allowed reflections
and zero on the others, and F(000) = n Z, so the weight is (f/Z)² exp(−2B s²)
and the BCC and FCC absence rules follow. With several kinds of atom the
pattern carries more information:

- NaCl (FCC, Na at 0, Cl at (½,0,0)): F = 4(f_Na + f_Cl) for h, k, l all even,
  4(f_Na − f_Cl) for all odd. In KCl, K⁺ and Cl⁻ have nearly the same f, so the
  odd reflections almost vanish and the pattern looks like SC with a/2.
- CsCl (SC, Cs at 0, Cl at (½,½,½)): F = f_Cs ± f_Cl; the h+k+l odd reflections
  are weak but present, so the cell is SC, not BCC.
- Diamond (FCC, atoms at 0 and (¼,¼,¼)): F = 4f (1 + i^(h+k+l)) for the FCC
  reflections, which is zero for h+k+l = 4n+2 (200, 222).
- Ordered Cu₃Au (SC, Au at 0, Cu at the face centres): F = f_Au + 3f_Cu on the
  FCC reflections and f_Au − f_Cu on the others, the superlattice lines.

## 8. Electron Diffraction (reciprocal-space form)

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
multiplicity × (f_e(s)/f_e(0))² × Debye–Waller factor, with f_e from the
Mott–Bethe relation (§7).
