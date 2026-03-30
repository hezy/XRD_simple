# Code Review — 10 commits (2026-03-31)

## Commit Summary

| # | Hash | Type | Description |
|---|------|------|-------------|
| 1 | e4b5ae3 | perf | optimize sum_peaks with in-place broadcasting |
| 2 | 6895600 | perf | remove dead pre-allocation in do_it and use in-place multiplication |
| 3 | 0db8ec6 | refactor | use typed exceptions in Miller_indices |
| 4 | 0a6928c | refactor | extract magic numbers into named constants |
| 5 | f9bae64 | refactor | add explicit Plots import and clean up dead imports |
| 6 | e4ca701 | fix | ensure results directory exists before saving |
| 7 | 10526c1 | docs | add full docstrings to all public functions |
| 8 | 3d5fe47 | feat | read lattice types dynamically from config |
| 9 | b185e03, ee4d1f8 | chore | remove unused make_noisy and estimate_peak_bounds functions |
| 10 | 72be414 | feat | add CLI argument parsing with ArgParse.jl |
| 11 | 2d7540a | refactor | rename the_title to title |
| 12 | f51cff7 | test | add comprehensive test suite with 90 tests |

---

## Detailed Review

### 1. e4b5ae3 — perf: optimize sum_peaks with in-place broadcasting

Replaces `y = y + ...` with `y .+= ...` and `zeros(size(θ))` with `zeros(length(θ))`. Good, correct micro-optimization.

### 2. 6895600 — perf: remove dead pre-allocation in do_it and use in-place multiplication

Removes unused `y = zeros(N)`, switches to `y .*= rand(...)`, and uses `.+` for background addition. Clean.

**Note:** This commit also changes `noise_level` in data.toml from 0.1 to 0.15 without mentioning it in the commit message.

### 3. 0db8ec6 — refactor: use typed exceptions in Miller_indices

Replaces generic `error()` calls with `ArgumentError`. Also removes the redundant `isa(min, Int)` check (already enforced by the type signature). Good improvement.

### 4. 0a6928c — refactor: extract magic numbers into named constants

Extracts background parameters and Miller index range into module-level constants. Clean and improves readability. Constants added:

- `AIR_SCATTER_AMPLITUDE`, `AIR_SCATTER_DECAY`
- `FLUORESCENCE_LEVEL`
- `AMORPHOUS_AMPLITUDE`, `AMORPHOUS_CENTER`, `AMORPHOUS_WIDTH`
- `MILLER_INDEX_MIN`, `MILLER_INDEX_MAX`

### 5. f9bae64 — refactor: add explicit Plots import and clean up dead imports

Replaces commented-out imports with explicit `using Plots`. Straightforward cleanup.

### 6. e4ca701 — fix: ensure results directory exists before saving

Adds `isdir("results") || mkdir("results")`. Simple and correct fix for fresh clones.

### 7. 10526c1 — docs: add full docstrings to all public functions

121 lines of docstrings added. Well-structured with Arguments/Returns/Throws format. Fixes typos ("colecting" -> "Collecting", "Lorntzian" -> "Lorentzian").

**Note:** Docstrings were added to `make_noisy`, which gets deleted two commits later. Wasted effort but no harm.

### 8. 3d5fe47 — feat: read lattice types dynamically from config

Replaces hard-coded ("SC", "BCC", "FCC") loop in main.jl with lattice types read from data.toml. The `isempty(elements)` check in `read_xrd_config` skips structures with no materials, but doesn't handle malformed entries (e.g., a structure key that maps to a non-dict value). Minor risk given the TOML schema.

### 9. b185e03, ee4d1f8 — chore: remove unused functions

Clean removal of `make_noisy` and `estimate_peak_bounds`. No issues.

### 10. 72be414 — feat: add CLI argument parsing with ArgParse.jl

Good addition with flags for `--config`, `--theme`, `--seed`, `--no-interactive`, and `--no-plots`. Code is clean and well-structured.

One minor observation: the DataFrame initialization `DataFrame(Dict(lt => θ₀ for lt in lattice_types))` initializes all lattice columns with θ₀ values that are immediately overwritten in the loop.

### 11. 2d7540a — refactor: rename the_title to title

Clean rename applied consistently in functions.jl and main_VScode.jl.

### 12. f51cff7 — test: add comprehensive test suite with 90 tests

8 test files covering crystal functions, peak profiles, width calculations, pattern computation, background, config parsing, and error handling.

---

## Issues Found

### Bugs / Correctness

**~~Tautological tests in test/test_config.jl (lines 22-24)~~ FIXED**

The assertions `haskey(x) || !haskey(x)` (always true) have been replaced with proper `haskey(lattice, "SC")` etc.

### Commit Hygiene

- **Undocumented config change:** Commit 6895600 silently changes `noise_level` in data.toml from 0.1 to 0.15 without mentioning it in the commit message.
- **Wasted work:** Commit 10526c1 adds docstrings to `make_noisy`, which commit b185e03 deletes. Better ordering would have avoided this.

### Test Structure

- **~~Repeated includes~~ FIXED:** functions.jl is now included once in runtests.jl instead of in every test file.
- **~~Relative config path~~ FIXED:** Tests now use `joinpath(@__DIR__, "..", "data.toml")` via a `DATA_TOML` constant.

---

## Overall Assessment

Solid set of changes. Well-structured, incremental, and following conventional commit conventions. The codebase is cleaner, better documented, and now has a proper test suite and CLI interface. All actionable issues from this review have been resolved.
