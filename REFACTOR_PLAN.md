# Refactor plan: physics fixes and structure

Plan from the September 2026 review of `main.jl` and `functions.jl`. Work
through the phases in order. Each phase ends with passing tests and one commit
(or one PR). Mark items done as they are finished, and update **Status**.

Physics *extensions* (Lorentz–polarization factor, structure factors, f_e(s))
are not part of this plan; they stay in `improvements.md`.

**Status:** Phases 1–5 done. Phase 6 not started.

---

## Phase 1: Fix the X-ray widths and the pseudo-Voigt profile

These errors change the output. Fix them first, so that the later phases have
correct reference values.

- [x] **1.1 Unit mismatch in Scherrer size broadening.**
  `Lorentzian_peaks_width` computes Kλ/(D cos θ) with λ in Å and D in nm, so
  the size term is 10× too large. Convert D to Å, as
  `Lorentzian_peaks_width_g` already does.
- [x] **1.2 Widths in 2θ, profile on θ.** Scherrer, Stokes–Wilson and Caglioti
  give FWHM in 2θ, but the grid and peak centres are in θ, so every X-ray
  peak is 2× too wide. Decision: compute the X-ray pattern on a 2θ grid
  (the axis the user sees and the axis the widths refer to). Keep radians
  internally.
- [x] **1.3 Caglioti docstring.** The formula gives FWHM², not HWHM². Also
  state the unit of the result (radians of 2θ).
- [x] **1.4 Pseudo-Voigt uses two widths instead of one.** A
  Thompson–Cox–Hastings pseudo-Voigt mixes a Lorentzian and a Gaussian that
  both have the combined FWHM f = `peak_fwhm(w_L, w_G)`. Rewrite
  `pseudo_Voigt_peak` accordingly. Then check whether the "Voigt is 2× broader"
  entry in `problems.md` is resolved, and remove it if so.
- [x] **1.5 Width evaluated at the peak centre.** The vector-width methods use
  w(θ[i]) at each grid point, so the width changes across one peak. Evaluate
  w_L, w_G once per reflection at its centre, and pass scalars to the profile.
  Then delete the vector methods of `Voigt_peak` and `pseudo_Voigt_peak`
  (the vector method of `peak_fwhm` may still be useful; keep it if used).
- [x] **1.6 Retune `data.toml` if needed.** After 1.1 and 1.2 the X-ray peaks
  are about 20× narrower for the size term. Check that the default D, ε and
  U/V/W still give readable patterns. `data.toml` is skip-worktree; change the
  committed defaults deliberately, not per-run toggles.
  Result: no change needed. FWHM is 0.3° to 1.4° of 2θ, grid step 0.11°.

**Tests to add:**
- Scherrer: one reflection, known K, λ, D, θ; compare FWHM with a hand calculation.
- Profile FWHM: measure the numerical FWHM of `Voigt_peak` and
  `pseudo_Voigt_peak` with the same w_L, w_G; both must agree with
  `peak_fwhm` within about 1 %.
- Symmetry: a single peak on a grid is symmetric about its centre.
- Area: with `normalize=false`, the numerical area equals A within about 1 %.

**Done when:** all tests pass, and `main.jl` produces plots that look correct
for SC, BCC and FCC in both radiation modes.

---

## Phase 2: Record a reference output

Phases 3–6 must not change the numbers. A reference output makes this
checkable.

- [x] **2.1** With a fixed seed and `noise_level = 0`, save the pattern for
  one sample per structure and per radiation mode to `test/reference/`.
- [x] **2.2** Add a test that recomputes these patterns and compares them with
  `isapprox` (relative tolerance about 1e-12).

  Result: `test/reference/{xray,electron}.toml` (fixed configs; Po-SC,
  Fe-BCC, Cu-FCC), saved patterns in `{xray,electron}.csv`,
  `test/test_reference.jl`. The patterns are computed through `simulate`
  (formerly `do_it`), in `test/reference/reference.jl`.
  Regenerate with `julia --project=. test/reference/generate.jl`, only when a
  change of the numbers is intended.

**Done when:** the reference test passes on the Phase 1 code.

---

## Phase 3: Parse the configuration once, into one typed value

Now: the file is parsed 2 + N times (`main`, `do_it_zero`, `do_it` per
sample), and defaults are applied by scattered `get(...)` calls that disagree
(`image_px` is 700 in `render_ring_image`, 800 in `main`).

- [x] **3.1** `read_xrd_config` returns one struct (or NamedTuple) holding
  instrument, peak-width and sample data. It applies every default and
  validates every value in one place, with descriptive `ArgumentError`s
  (missing key, negative width, unknown `radiation`).
- [x] **3.2** All other functions receive this value, not a file name. Remove
  every `get(instrument, ..., default)` outside `read_xrd_config`.
- [x] **3.3** Delete `do_it_zero`. It only builds the first DataFrame column,
  which the loop overwrites. (With zero samples it writes θ in radians under
  the label "2θ (deg)".) Build the x column from the first computed pattern,
  or from the grid function of Phase 4.
- [x] **3.4** Update `test/test_config.jl`.

  Result: `XRDConfig` struct; `read_xrd_config` has a file and a `Dict`
  method (the tests use the `Dict` one). Parameters of the unused mode are
  `NaN` when absent. Functions that took `instrument`/`peak_width` now take
  `cfg`; `intensity_vs_angle`, `compute_xrd_pattern` and `compute_peak_widths`
  take λ from it. `render_ring_image(g, y, cfg)` has no keyword defaults.

**Done when:** the reference test and all other tests pass.

---

## Phase 4: Select the radiation mode by dispatch

Now: the mode is checked in `main`, `do_it_zero` and `do_it`;
`do_it`/`do_it_electron` and `compute_xrd_pattern`/`compute_ed_pattern` are
near copies.

- [x] **4.1** Define `abstract type Radiation end` with `struct XRay` and
  `struct Electron`. Each holds its own parameters (λ and U/V/W; voltage and
  G_inst). `read_xrd_config` constructs the right one.
- [x] **4.2** Give each type a small set of methods:
  - `grid(mode)`: x axis (2θ or g)
  - `max_hkl_sq(mode, a)`: reflection cutoff
  - `peak_centres(mode, indices, a)`: centre and multiplicity of each visible reflection
  - `peak_widths(mode, x₀, sample)`: (w_L, w_G) at one centre
  - `background(mode, x)`
  - `axis_label(mode)`
- [x] **4.3** Write one generic `simulate(mode, cfg, structure, a)` that
  returns `(x, y)`. It replaces `do_it`, `do_it_electron`,
  `compute_xrd_pattern`, `compute_ed_pattern`, `intensity_vs_angle` and
  `intensity_vs_g`.
- [x] **4.4** Remove the unused `noise_level` keyword of `background` and
  `background_electron`; noise is applied once, in `simulate`.

  Result: `XRDConfig` holds `mode::Radiation` (`XRay` or `Electron`, with the
  instrument parameters of that mode, including the ring settings) in place of
  `radiation` and the per-mode fields. Keys of the unused mode are no longer
  read or checked. `simulate(cfg, structure, a)` takes the mode from `cfg`,
  and returns x in display units (2θ in degrees, or g). Mode methods:
  `grid(mode, N)`, `max_hkl_sq`, `peak_centres(mode, indices, multiplicities,
  a)`, `peak_widths(mode, x₀, cfg)`, `background`, and in addition
  `display_axis`, `axis_label`, `plot_title`. The width functions take a
  scalar angle or g. `plot_pattern(mode, x, y, title, theme)` holds the plot
  code of the old `do_it` until Phase 5; `render_ring_image` takes the
  `Electron` mode; `write_ring_outputs` in `main.jl` has an empty `XRay`
  method. The electron plot x label is now "g (1/Å)", the same as the CSV
  column (was "g = 1/d (1/Å)").

**Done when:** no `== "electron"` test remains outside `read_xrd_config`, and
all tests pass.

---

## Phase 5: Separate plotting from physics

Now: `functions.jl` loads Plots and calls `theme()`, so the tests load Plots.

- [x] **5.1** Physics functions return numbers only. `simulate` returns
  `(x, y)`; the title is built by the caller.
- [x] **5.2** `render_ring_image` is split: `ring_image(...)` returns the
  matrix (physics, testable); a plotting function draws it with the colormap.
- [x] **5.3** All `Plots` calls (`theme`, `plot`, `heatmap`, `savefig`) live
  in the plotting file and in `main.jl`.
- [x] **5.4** The tests no longer load Plots, except `test/ring_sanity.jl` if
  it needs to.

  Result: new file `plotting.jl` (included by `main.jl` after `functions.jl`)
  holds `using Plots`, `PHOSPHOR_RAMP`, `plot_title`, `plot_pattern` and
  `plot_ring_image(coords, img, mode)`. `ring_image(g, y, mode)` in
  `functions.jl` returns `(coords, img)`, the normalised, gamma-compressed
  matrix; it replaces `render_ring_image`. Neither the tests nor
  `test/ring_sanity.jl` load Plots. The CSV, the pattern plots and the ring
  images of the reference configs are byte-identical to the Phase 4 output.

**Done when:** `functions.jl` (or its successors) contains no `using Plots`.

---

## Phase 6: Split `functions.jl` into files

- [ ] **6.1** Split by subject, included in this order by one entry file:
  - `config.jl`: config struct, `read_xrd_config`
  - `crystal.jl`: `Miller_indices`, `cubic_multiplicity`, `d_list`, `g_list`
  - `profiles.jl`: `Voigt_peak`, `pseudo_Voigt_peak`, `peak_fwhm`, `sum_peaks`
  - `xray.jl`: `XRay` methods, Bragg angles, Caglioti, Scherrer
  - `electron.jl`: `Electron` methods, `electron_wavelength`, ring image matrix
  - `plotting.jl`: pattern plot, ring plot
- [ ] **6.2** Decision to take at this point: plain `include` files in `src/`,
  or a module `XRDSim` (then the tests use `using XRDSim`). A module is
  cleaner for the tests; plain includes are simpler for students. Choose
  before starting 6.1.

**Done when:** `main.jl` and all tests run unchanged in behaviour.

---

## Phase 7: Small cleanups and documentation

- [ ] **7.1** Change the banner strings ("Constants", "Functions",
  "Electron diffraction (1D)") and the orphan `abstract_peak` docstring into
  comments. As strings, Julia attaches them to the next expression.
- [ ] **7.2** Relax argument types from `Vector{Float64}` to `AbstractVector`
  and `Float64` to `Real` where nothing depends on the concrete type.
- [ ] **7.3** Remove validation that is repeated in both caller and callee.
- [ ] **7.4** `using Distributions: Normal` (from `improvements.md`, item 1).
- [ ] **7.5** Update `CLAUDE.md` (architecture, file list, function names),
  `README.md`, `problems.md` and `improvements.md` to match the new code.

**Done when:** the documentation names only functions that exist.
