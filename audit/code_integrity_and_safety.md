# Code Integrity, Numerical Stability, Concurrency & Safety Audit

**Scope**: All 6 Rust crates in the workspace (`astrosat-fits`, `astrosat-core`, `astrosat-czti`, `astrosat-laxpc`, `astrosat-cli`, `astrosat-bench`).  
**Audit Verdict**: **PASS — All reported `HIGH` and `MEDIUM` integrity, safety, numerical stability, determinism, and CLI/config findings have been fixed and re-verified.**

---

## Executive Summary

| Severity | Found | Fixed & Verified | Open |
|---|---:|---:|---:|
| **CRITICAL** | 0 | 0 | 0 |
| **HIGH** | 9 | 9 | 0 |
| **MEDIUM** | 13 | 13 | 0 |
| **LOW** | 2 | 1 | 1 (documented `mmap` external-truncation caveat) |
| **INFO** | 4 | — | — |

### Post-Fix & Post-Neon Verification Summary
- **Workspace Clippy**: `cargo clippy --workspace --all-targets -- -D warnings` — **0 warnings, 0 errors**.
- **Crate Test Suites (all executed post-Neon fixes, 0 failed)**:
  - `astrosat-fits`: `cargo test -p astrosat-fits` — **7 passed, 0 failed** (including `malformed_fits_headers_and_heap_rejected_without_panic` with `TZERO = 1e30` & `2^63` overflow regressions).
  - `astrosat-core`: `cargo test -p astrosat-core` — **7 passed, 0 failed** (including `safety_edge_cases_nan_double_bary_and_primary_dates`).
  - `astrosat-czti`: `cargo test -p astrosat-czti` — **8 passed, 0 failed** (including `audit_integrity_and_legacy_parity_regression_suite` and `real_obsid_9000000406_end_to_end_pipeline_and_reference_parity`).
  - `astrosat-laxpc`: `cargo test -p astrosat-laxpc --test unit_tests` — **8 passed, 0 failed** (`audit_regressions_safety_math_and_edge_cases` + 7 stage unit tests; full `--release` suite including 2/2 Crab `9000000252` integration tests passed 10/10 in `#114`).
  - `astrosat-cli`: `cargo test -p astrosat-cli` — **10 passed, 0 failed** (`3` unit tests including `z2_fap_never_produces_nan_on_large_z2_and_high_harmonics` + `7` integration tests including `cli_defaults_symmetrically_override_config_and_energy_range_requires_caldb` and `neon_review_regressions_negative_args_mkf_precedence_and_timing_guards`).
  - `astrosat-bench`: `cargo test -p astrosat-bench` — **3 passed, 0 failed**.
- **Empirical Thread Determinism (post-Neon `target/release/astrosat`, `--threads 1` vs `--threads 16`)**:
  - `astrosat czti run --gen-evt` (ObsID `9000000406`): **31/31** output FITS files bit-for-bit identical (SHA-256).
  - `astrosat laxpc run --format-a` (ObsID `9000000252`): **14/14** output FITS files bit-for-bit identical (SHA-256).
  - `astrosat bary` (using `AS1T01_112T01_9000000406czt_level1.orb`): bit-for-bit identical (`cmp -s`).

---

## 1. Safety & Panic Surface (`unsafe`, `.unwrap()`, Indexing, Integer Casts, Memory Bounds)

### `unsafe` Inventory
- **[INFO / LOW] Single `unsafe` block across entire workspace** (`crates/astrosat-fits/src/io.rs:68`):
  `unsafe { Mmap::map(&file) }` opens FITS files read-only and falls back to `std::fs::read` if mmap fails (`io.rs:68-74`).
  *Residual caveat (LOW)*: Standard POSIX/OS `mmap` semantics apply—external truncation of an open FITS file on disk while mapped can trigger `SIGBUS`; in-memory `FitsFile::from_bytes` avoids `mmap`.

### Findings & Verified Fixes

1. **[HIGH — FIXED & VERIFIED] Negative FITS header integers, invalid `BITPIX`, and unchecked `usize`/heap arithmetic** (`@fits-core`, `#46` -> `#56`)
   - `crates/astrosat-fits/src/io.rs:126-152,254-294`: `compute_data_len` now validates `Bitpix::from_i64(bitpix)` (rejecting non-standard `BITPIX` such as `9`), checks `0 <= NAXIS <= 999`, converts `NAXISi`, `PCOUNT`, and `GCOUNT` via `usize::try_from`, and uses `checked_mul`/`checked_add`/`checked_next_multiple_of` in both `compute_data_len` and `from_store` (`data_start + data_len <= buf.len()`).
   - `crates/astrosat-fits/src/bintable.rs:141-210,774-830`: Replaced unchecked `as usize` casts on `NAXIS1`, `NAXIS2`, `PCOUNT`, `TFIELDS`, and `THEAP` with `usize::try_from`, guarded `tfields <= header.cards().len()` before `Vec::with_capacity`, and added checked arithmetic (`checked_mul`/`checked_add`) for `parse_tform`, `NAXIS1 * NAXIS2`, `THEAP + PCOUNT`, and `P32`/`Q64` heap slices in `cell_array_bytes`.

2. **[HIGH — FIXED & VERIFIED] Negative `NAXIS` and pixel-byte overflow in `ImageRef`** (`@fits-core`, `#46` -> `#56`)
   - `crates/astrosat-fits/src/image.rs:51-98,120-136`: Rejects negative/`> 999` `NAXIS` and negative `NAXISi` via `usize::try_from`; computes `num_pixels` via `try_fold(..., checked_mul)` and byte count via `checked_mul(bitpix.byte_size())` with EOF verification in both `from_header_and_data` and `checked_data_bytes`.

3. **[HIGH — FIXED & VERIFIED] Unconditional `xs[0]` indexing before length validation in `OrbitTable::from_fits`** (`@fits-core`, `#46` -> `#56`)
   - `crates/astrosat-core/src/orbit.rs:133-155`: Validates `time.len() < 2` and matching column lengths before indexing `xs[0]` for unit inference.

4. **[HIGH — FIXED & VERIFIED] Unchecked JPL DE ephemeris descriptors and `span == 0.0` in `JplDeEphemeris`** (`@fits-core`, `#46` -> `#56`)
   - `crates/astrosat-core/src/bary.rs:1069-1097,1150-1155,1184-1220`: Validates required body `ipt` entries (`2`, `9`, `10` with `start_1 > 0`, `ncf > 0`, `nsub > 0`), positive finite `span_days`/`emrat`, per-record `span > 0.0`, and checked coefficient slice bounds in `eval_body` returning `Result`.

5. **[HIGH — FIXED & VERIFIED] Cross-frame event column desynchronization in CZTI `science2event`** (`@czti-dev`, `#47` -> `#60`)
   - `crates/astrosat-czti/src/science2event.rs:359-421,566-696`: 3-record onboard bunch reports (`nframe_bunch == 0/1/2`) are now staged in `PendingOnboardBunch` and committed atomically across `bunches` and all `QuadEvents` columns on record 2, eliminating mid-bunch cross-frame column skew.

6. **[HIGH — FIXED & VERIFIED] Unchecked `lld_pi` indexing, unbounded `nbins`, and duplicate `if/else` branch in `noiseclean`** (`@czti-dev`, `#47` -> `#60`)
   - `crates/astrosat-czti/src/noiseclean.rs:523-524`: `heavybunchclean_quad_with_exposure` uses bounds-checked `lld_pi.get(dpix).copied().unwrap_or(0)`.
   - `crates/astrosat-czti/src/noiseclean.rs:606-690`: `flickpixclean_quad` validates positive finite `flicker_time_bin_sec` and `t1 >= t0`, and steps directly from populated bin to populated bin (`O(N)` regardless of timestamp gaps) instead of allocating `Vec` over `(t1 - t0) / bin_sec`.
   - `crates/astrosat-czti/src/noiseclean.rs:810-845`: `eventsep_quad` classifies contiguous clusters into suppressed (`0` for `cluster_len >= 3`), singles (`1`), and doubles (`2`), preventing 3rd+ burst events from leaking into singles.

7. **[HIGH — FIXED & VERIFIED] Unchecked 64x64 indexing in `image.rs`, `.expect()` panic & unchecked response bounds in `bindata.rs`** (`@czti-dev`, `#47` -> `#60`)
   - `crates/astrosat-czti/src/image.rs:250-255,337-348`: `balanced_cross_correlation_from_mask` and `fft_cross_correlation_with_shift` use `.get(y).and_then(|r| r.get(x))` instead of direct indexing.
   - `crates/astrosat-czti/src/bindata.rs:866-880,1036-1126`: Removed `.expect("mask weights required...")` from `bindata_quad`; guarded `n_ebins == 0 || n_pibins == 0` in `generate_quad_response` alongside checked `grp_base` multiplication and slice bounds checks on `resp_id`, `lld`, `uld`, `pixresp_matrix`, `eff_area`, `avg_open_frac`, `pix_exp_frac`, and `weights_512`.

8. **[HIGH — FIXED & VERIFIED] Signed `i32` casts, `count * 8` overflow, and `assert!` panic in LAXPC 2D background parser** (`@laxpc-dev`, `#50` -> `#61`)
   - `crates/astrosat-laxpc/src/background.rs:780-803,867-872,979-993`: Enforces `no >= 2`, `na >= 2`, `iw > 0`, `mbk > 0`, checked `(no + 2) * (na + 2)`, and upper bounds before converting to `usize`; `read_f64_vec` uses `count.checked_mul(8)` and validates `rec.len() >= needed` prior to `Vec::with_capacity(count)`; `eval_bspline_basis` returns safely on `nx < 2 || k == 0` instead of panicking on `assert!`.

9. **[MEDIUM — FIXED & VERIFIED] Card UTF-8/i64 handling, `N_EBIN`/`N_PIBIN` casts, unbounded MKF/BBC allocations, `--seed` overflow, and empty-row bench indexing**
   - `crates/astrosat-fits/src/card.rs:533-545,639-661`: `parse_quoted_string` iterates `char_indices()` (preserving valid UTF-8); `get_i64` validates finite integral `f64` values within `[i64::MIN, 2^63)`.
   - `crates/astrosat-czti/src/fits_io.rs:669-677`: Rejects `N_EBIN <= 0 || N_PIBIN <= 0` before casting to `usize`.
   - `crates/astrosat-laxpc/src/pipeline.rs:999-1015` & `bb.rs:151-182`: `write_orbit_from_mkf` rejects non-finite/inverted spans and bounds duration (`<= MAX_ORBIT_SPAN_SEC = 40_000_000.0`); `BbcSeries::from_records` guards non-positive/non-finite `bin_sec` and caps bins at `MAX_BBC_BINS = 5_000_000`.
   - `crates/astrosat-cli/src/czti_cmd.rs:917-921`: Clamps `--seed` to `1..=i32::MAX - 12` before `+ 6` / `+ 12`.
   - `crates/astrosat-bench/src/bench_fits.rs:48-57`: Uses `slice.first()` / `slice.last()` so `NAXIS1 == 0` BINTABLE rows cannot panic.

---

## 2. Numerical Stability & Domain Guards

### Timestamp Precision (`f64` vs `f32`)
- **[INFO] Verified `f64` precision across all timestamps**: All event times, `TSTART`/`TSTOP`, GTI bounds, orbit tables, TCT lookups, and barycentric corrections use IEEE 754 `f64`. Where `f32` is computed (`evtencor.rs:92` days since gain-drift epoch; `pha2energy.rs:144` normalized time fraction), epoch offsets are subtracted in `f64` prior to `f32` conversion.

### Findings & Verified Fixes

1. **[MEDIUM — FIXED & VERIFIED] `NaN` range-check bypasses & extreme-timestamp panics in `astrosat-core`** (`@fits-core`, `#46` -> `#56`)
   - `crates/astrosat-core/src/orbit.rs:197-205` & `bary.rs:1114-1117`: Explicit `!t_sec.is_finite()` and `!jd.is_finite()` checks prevent `NaN` from bypassing `< min || > max` guards.
   - `crates/astrosat-core/src/tct.rs:254-256`: `CztiTct::interpolate` returns `Err(CoreError::TimeOutOfRange(cztseccnt))` for non-finite `cztseccnt` or `feb_lag` instead of `Ok(0.0)`.
   - `crates/astrosat-core/src/time.rs:312-339`: `astrosat_utc_to_iso8601` guards non-finite and out-of-range timestamps via `(-6.0e12..=6.0e12).contains(&floored)` and `Duration::try_seconds` + `checked_add_signed`.

2. **[HIGH — FIXED & VERIFIED] `NaN` FOV bypass and non-positive bin / zero-`FRACEXP` hazards in CZTI `bindata`** (`@czti-dev`, `#47` -> `#60`)
   - `crates/astrosat-czti/src/bindata.rs:363-375`: `validate_fov` rejects `!tx_deg.is_finite() || !ty_deg.is_finite()`.
   - `crates/astrosat-czti/src/bindata.rs:789-810,967-975`: Sanitizes `time_binsize <= 0.0` / non-finite, caps `n_lcbins <= MAX_LC_BINS = 10_000_000`, and requires `fe.is_finite() && fe > min_fe && fe > 0.0` (`min_fe = cfg.min_fracexp.max(0.0)`) before dividing by `(fe as f64) * time_binsize`.

3. **[MEDIUM — FIXED & VERIFIED] Non-positive `tbinsize` div-by-zero & empty-event CALDB badpix bypass in `pixclean` / `bunchclean`** (`@czti-dev`, `#47` -> `#60`)
   - `crates/astrosat-czti/src/pixclean.rs:150-185`: Guards `det_tbinsize` / `pix_tbinsize` against non-finite/`<= 0.0` inputs and masks CALDB bad pixels (`>= NOISY`) to `0.0` exposure even when `nrows == 0`.
   - `crates/astrosat-czti/src/bunchclean.rs:258-271`: Sanitizes `livetime_binsize <= 0.0` / non-finite and bounds `num_timebins <= MAX_LV_BINS = 10_000_000`.

4. **[MEDIUM — FIXED & VERIFIED] `0.0 * Inf = NaN` in timing `z2_fap`** (`@cli-dev`, `#48` -> `#64`)
   - `crates/astrosat-cli/src/timing.rs:429-451`: Rewrote `z2_fap` term-by-term in log space (`(-x + sum ln(x/j)).exp()`), eliminating underflow/overflow `0.0 * Inf = NaN` for large $Z_n^2$ and high `--harmonics`, with finite clamping in `fap_to_sigma`.

5. **[INFO] Verified Well-Guarded Domains**
   - `crates/astrosat-laxpc/src/products.rs`: Non-paralyzable deadtime denominator `1.0 - R_all * tau` and `RATE` zero-exposure paths explicitly guard against non-positive exposure and denominator saturation.
   - `crates/astrosat-czti/src/image.rs` (`balanced_cross_correlation`) and `evtencor.rs`: Zero-variance / non-positive stdev paths are guarded cleanly.

---

## 3. Parallelism & Determinism (`rayon`)

- **[INFO] CZTI `WichmannHillRng` & quadrant parallelism**: `WichmannHillRng` is instantiated deterministically per quadrant (`evtencor.rs`, `bindata.rs:784`) rather than shared across worker threads. Verified bit-for-bit identical across `--threads 1` vs `--threads 16` (**31/31 files**).
- **[INFO] Barycentric correction (`astrosat bary`)**: Event-parallel ephemeris and Roemer delay evaluation preserves row order deterministically. Verified bit-for-bit identical across `--threads 1` vs `--threads 16`.
- **[MEDIUM — FIXED & VERIFIED] Unstable parallel sort on tied timestamps in LAXPC Format A** (`@laxpc-dev`, `#50` -> `#61`)
  - `crates/astrosat-laxpc/src/pipeline.rs:1126-1133`: Replaced `par_sort_unstable_by(|a, b| a.time_utc.total_cmp(&b.time_utc))` with deterministic `par_sort_by` and total tie-breakers `(time_utc, detector_id, anode_id, channel, layer)`. Verified bit-for-bit identical across `--threads 1` vs `--threads 16` (**14/14 files**).

---

## 4. CLI & Config Integrity (`astrosat-cli` & `astrosat-laxpc`)

1. **[MEDIUM — FIXED & VERIFIED] Symmetric CLI vs `--config` precedence & `CztiFileConfig::output_type`** (`@cli-dev`, `#48` -> `#64`)
   - `crates/astrosat-cli/src/config.rs:46`: Added `output_type: Option<String>` to `CztiFileConfig`.
   - `crates/astrosat-cli/src/cli.rs:8-38`, `czti_cmd.rs:307-383`, `laxpc_cmd.rs:220-338`, `bary_cmd.rs:57-89`, `main.rs:125-150`: Replaced default-value comparisons (`!= default` and asymmetric `||` bools) with `clap::ArgMatches::value_source(ValueSource::CommandLine)` (`leaf_matches`, `is_cli_arg`, `resolve_flag`) and paired positive/negative boolean flags (`--checksum`/`--no-checksum`, `--mask-weight`/`--no-mask-weight`, `--resp`/`--no-resp`, `--gen-evt`/`--no-gen-evt`, `--format-a`/`--no-format-a`, `--format-b`/`--no-format-b`, `--orbit`/`--no-orbit`, `--ascii`/`--no-ascii`, `--single-only`/`--no-single-only`, `--layer1-only`/`--no-layer1-only`, `--user-gti-only`/`--no-user-gti-only`), allowing explicit CLI defaults to override `--config`.

2. **[MEDIUM — FIXED & VERIFIED] Explicit error when `--energy-range` lacks CALDB `.channel` map** (`@cli-dev` `#48` -> `#64`, `@laxpc-dev` `#50` -> `#61`)
   - `crates/astrosat-laxpc/src/pipeline.rs:353-360` & `crates/astrosat-cli/src/laxpc_cmd.rs:236-242,547-560`: Both `astrosat laxpc run` and `astrosat laxpc lightcurve` validate `lxp{10,20,30}.channel` upfront and `run_unit_internal` returns `LaxpcError::Config` if `energy_range_kev` is requested without an available channel map, preventing silent unfiltered fallback.

3. **[LOW — FIXED & VERIFIED] `(0.0, 0.0)` RA/DEC sentinel ambiguity in CZTI `bindata`** (`@czti-dev`, `#47` -> `#60`, `@cli-dev`, `#64`)
   - `crates/astrosat-czti/src/pipeline.rs:1022` & `crates/astrosat-cli/src/czti_cmd.rs:346-348`: Added `BinDataConfig::source_coords: Option<(f64, f64)>` + `resolve_source_radec` so explicitly targeting `(0.0, 0.0)` is distinguished from omitting `--ra`/`--dec`.

---

## 5. Neon Fresh-Eyes Self-Review Addendum (Post-Oracle Codebase)

Following the real-data oracle wave (`#74`–`#94`), a complete fresh-eyes static and manual re-audit was executed across all 6 crates, covering both pre-existing code and all post-Wave-1 additions (`JplDeEphemeris::from_heasarc_fits_bytes`, 16-sig-fig `format_fits_float`, CZTI `compute_spectrum_quality_and_grouping`/`read_response_file`/FFT & `SHADOWLIB` fixes, LAXPC `BackgroundCatalog`/`ChannelMap::from_caldb_for_day`, and CLI `laxpc spectrum --energy-range`/`czti pha2energy --ss`).

### Workspace-Wide Static Invariants Re-Verified
- **`unsafe` blocks**: Still exactly **1** (`crates/astrosat-fits/src/io.rs:68`).
- **Sorting & Panic Hygiene**: `0` `par_sort_unstable` calls, `0` float `partial_cmp().unwrap()` calls, `1` `.expect()` (`checksum.rs:96`, proven ASCII bytes), and all `asin`/`acos`/`ln`/`log10` call sites domain-guarded (`bary.rs:220` clamps `cos_theta` to `[-1.0 + 1e-12, 1.0]`, `bindata.rs:656-672` restricts `absco_ta` to `e >= 1.0`, `gaas.rs:96` and `bindata.rs:238,255,295` clamp to `[-1.0, 1.0]`).

### Additional Findings Caught & Verified Fixed in Neon Review

1. **[MEDIUM — FIXED & VERIFIED] Debug-mode integer overflow panics on untrusted FITS `TZERO` values** (`@fits-core`, `#102` -> `#107`)
   - `crates/astrosat-fits/src/bintable.rs:438-493` (`read_col_i64`): Previously saturated huge/`2^63` `TZERO` to `i64::MAX` and executed unchecked `raw + z` on `I32`/`I16`/`Byte` arms (panicking on `raw > 0` in debug builds). Now guards `tzero.is_finite() && tzero.fract() == 0.0`, maps standard unsigned 64-bit `TZERO = 2^63` (`9_223_372_036_854_775_808.0`) to `i64::MIN`, bounds `[i64::MIN, 2^63)`, and uses `raw.wrapping_add(z)` across all 4 integer arms.
   - `crates/astrosat-fits/src/bintable.rs:660-668` (`read_vec_u16`): Replaced unchecked `(raw + tzero) as u16` with `raw.wrapping_add(tzero) as u16`.

2. **[MEDIUM — FIXED & VERIFIED] Unchecked `DETCHANS` allocation overflow & unbounded FFT `oversampling`** (`@czti-dev`, `#103` -> `#117`)
   - `crates/astrosat-czti/src/fits_io.rs:1903-1922` (`read_response_file`): Now rejects `DETCHANS <= 0`, computes `n_ebins.checked_mul(n_pibins)`, and rejects matrices exceeding `MAX_RSP_MATRIX_ELEMS = 64_000_000` with `Err(CztiError::InvalidConfig(...))`.
   - `crates/astrosat-czti/src/image.rs:346` (`fft_cross_correlation_with_shift`): Clamps `ov = oversampling.clamp(1, 64)` to prevent `usize` overflow / OOM on huge oversampling factors.

3. **[MEDIUM — FIXED & VERIFIED] `NaN`/`Inf` timing parameter bypass & unbounded `--n-trials` OOM** (`@cli-dev`, `#104` -> `#116`)
   - `crates/astrosat-cli/src/timing.rs:74-103,203-275`: Added `.is_finite()` validation on `freq_min`, `freq_max`, `freq_hz`, `period_sec`, `fdot`, `fddot`, `search_span_hz`, `epoch`, and `df_hz`; enforced explicit bounds `--n-trials` in `1..=2_000_000`, `--harmonics` in `1..=1024`, and `--nbins` in `2..=65536`.
   - `crates/astrosat-cli/src/czti_cmd.rs:399-404`: Validated `--oversampling` in `1..=64` on CLI entrypoints.
