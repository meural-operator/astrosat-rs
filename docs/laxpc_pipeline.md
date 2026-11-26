# AstroSat LAXPC Pipeline (`astrosat-laxpc`)

`astrosat-laxpc` is a pure-Rust, `rayon`-parallelized implementation of the AstroSat LAXPC Level-1 to Level-2 reduction pipeline for all three proportional counter units (**LAXPC10** / `lxp1`, **LAXPC20** / `lxp2`, **LAXPC30** / `lxp3`). It preserves bit-level compatibility with TIFR `LaxpcSoft v3.4.6` (`laxpcl1.f` and `backshiftv3.f`) while removing fixed-size Fortran array limits, `call system('rm ...')` shell calls, byte-swapping fragility, and interactive `.inp` prompts.

---

## 1. Pipeline Architecture & Module Reference

### 1.1 Frame Parsing, Sorting & Deduplication (`frames.rs`)
Implemented in [`frames.rs`](../crates/astrosat-laxpc/src/frames.rs):
* **2242-Byte Telemetry Frames (`L1_FRAME_SIZE = 2242`)**: Parses `LAXPC_MODEEA` and `LAXPC_MODEBB` `BINTABLE` rows directly from memory-mapped FITS bytes.
* **Header Integrity (`FrameHeader`)**: Extracts `l1_time` (`0..8`), 64-bit sequence `raw_seq` (`36..44`), `FRAME_SYNC_BYTE = 0xDE` (`44`), `unit_id` (`45`, `1..=3`), `mode_byte` (`46`: `MODE_EA_BYTE = 0xEA` or `MODE_BB_BYTE = 0xBB`), 56-bit and 40-bit `100 kHz` STBG clock (`47..54`, `49..54`), `status_flag == 0` (`57`), 16-bit onboard `frame_counter` (`58..60`), and error flags `err_flags = [ie1, ie2, ie3]` (`2092..2098`, requiring `ie1 <= 1`).
* **Unwrapping, Stray Rollback & Deduplication (`sort_and_dedup_frames`)**:
  * Detects and rolls back stray sequence jumps in `(6000, 60000)`.
  * Unwraps 16-bit (`65536`) frame-counter rollovers and long-gap (`> 2000 s`) wraps into monotonic `i64` counters stored in a `BTreeMap<i64, OrderedFrame>`.
  * Replaces/ignores duplicate frames and records deterministic `FrameDiagnostics` (`single_frame_gaps`, `multi_frame_gaps`, `total_missing_frames`).
* **Multi-Orbit Merging (`merge_multi_orbit_frames`)**: Eliminates boundary duplicates sharing `(frame_counter, stbg_40)` across consecutive orbits.

### 1.2 Event Analysis (`modeEA`) Unpacking (`ea.rs`)
Implemented in [`ea.rs`](../crates/astrosat-laxpc/src/ea.rs):
* **STBG Clock & 5th-Byte Rollover Fix (`decode_ea_frame`)**: Reconstructs the 40-bit `100 kHz` (`STBG_CLOCK_HZ = 100_000.0`) header time and applies the Dec 2024 5th-byte rollover correction using `STBG_EPOCH_OFFSET_SEC = 1.8121823219164e8` and `POW256_4 = 4_294_967_296.0`.
* **5-Byte Word Stream (`pos = 60..=2091`)**:
  * Stops at `END_OF_FRAME_MARKER = 0xEE` (or truncates non-monotonic frame-end tails when `na0 > 4` or `0xEE` is hit).
  * Processes `COARSE_TIME_MARKER = 0xEF` coarse time updates with carry propagation, single-marker glitch filtering, and retroactive alignment of pre-marker events.
  * Unpacks 5-byte event words (`RawEaWord::unpack`): `fine_ticks = i4.wrapping_sub(2)` (`10 us` resolution), 4-bit anodes `ia1`, `ia2`, 10-bit pulse heights `ic1`, `ic2`, and Xe K X-ray flags `k1`, `k2`.
* **Single & Xe K-Escape Double Events (`RawEaWord::resolve`)**:
  * **Single**: `ic2 == 0` $\implies$ `(ic1, ia1, EventKind::Single)`.
  * **Xe K-escape double**: `(k1 == 1 || k2 == 1) && ia1 < 8 && ia2 < 8` $\implies$ `pha = min(ic1 + ic2, 1023)`, assigning the dominant non-escape anode (`ia2` if `k1 == 1 && (ic1 <= ic2 || (k2 == 0 && ic2 > 1))`, else `ia1`).
  * **Validity & Filtering (`decode_ea_orbits`, `passes_event_filter`)**: Keeps `time_ins > 1.0`, `1 <= pha <= 999`, `1 <= anode <= 10`, detects internal frame gaps (`>= 50 ms`) and large gaps (`>= 10.0 s`), and filters by `EventFilter` (`single_only`, `AnodeSelection` `A1..A10`, `pha_min..=pha_max` or `energy_range_kev` — explicitly returning an error if `energy_range_kev` is requested without `--caldb` / `lxp{10,20,30}.channel`).

### 1.3 Broad Band Counting (`modeBB`) Decoding (`bb.rs`)
Implemented in [`bb.rs`](../crates/astrosat-laxpc/src/bb.rs):
* Decodes `31 x 64`-byte subframes per 2242-byte BB frame (`decode_bb_frame`), accumulating `lld` (12-bit), `jan` (12-bit), `ldu` / ULD (11-bit), and `ige` (11-bit) into 16 `BbcRecord`s per frame (emitting on odd subframes `j = 1, 3, ..., 31`, with `> 2.06 s` timestamp jump repair).
* Accumulates `BbcSeries` into `BBC_BIN_SEC = 32.768 s` (`0.128 * 256 s`) bins, providing deadtime-corrected linear (`NUSE=2`) interpolation for `uld_rate_cps`, `rejected_rate_cps`, and 4-quartile ULD binning (`uld_quartile_bounds`: `(750, 1050)` for Lxp10, `(780, 1100)` for Lxp20, `(740, 1090)` for Lxp30).

### 1.4 TCT Conversion & Standard GTI Screening (`gti.rs`)
Implemented in [`gti.rs`](../crates/astrosat-laxpc/src/gti.rs):
* Merges per-orbit `.tct` and `.mkf` files (`load_laxpc_tct_files`, `load_merged_mkf_files`), rejecting TCT clock glitches outside `(1.000, 1.001]` via `astrosat_core::LaxpcTct`.
* Screens each time bin against `GtiConfig` (`generate_standard_gti`):
  * **Earth Elevation**: `elv * je > min_elv_deg` (default `3.5°`; `je = -1` when `earth_occultation = true`).
  * **SAA Exclusion**: Uses wrap-safe longitude interpolation (`interp_lon_deg`, guarding `359° -> 1°` wraparound) with `lon < saa_entry_lon_deg` (`249.0°`) and `lon > saa_exit_lon_deg` (`3.0°`), plus optional diurnal SAA exit extension to `40.0°` (`diurnal_period_sec = 84_494.9 s`).
  * **Active Detector Threshold**: All-main-anode rate `>= 80.0` cps (`Lxp10`/`Lxp20`) or `50.0` / `20.0` cps (`Lxp30`).
  * **LAXPC10 ULD Anomaly**: Time-dependent `lxp10_uld_anomaly_threshold` post-anomaly cut.
  * **Post-processing**: Merges gaps `< merge_gap_sec` (`0.3 s`), drops intervals `< min_duration_sec` (`5.0 s`), optionally intersects `user_gti_utc` via strict set intersection (`gti_utc.intersect(user_gti)`), and computes `frame_gap_loss_sec` (`net_exposure_sec = gross_exposure_sec - frame_gap_loss_sec`).

### 1.5 Non-Paralyzable Deadtime & OGIP Products (`products.rs`)
Implemented in [`products.rs`](../crates/astrosat-laxpc/src/products.rs):
* **Deadtime Formula (`deadtime_correct_rate`, `deadtime_factor`)**:
  $$R_{\text{corr}} = \frac{R_{\text{obs}}}{1 - \tau\,R_{\text{all}}}, \qquad \tau = \text{DEFAULT\_DEADTIME\_SEC} = 43.0\,\mu\text{s}.$$
  * **Singularity Guard**: When $R_{\text{all}} \ge 23\,000\text{ cps}$ (`DEADTIME_RATE_GUARD_CPS`) or $1 - \tau R_{\text{all}} \le 0$, returns $R_{\text{obs}}$ unscaled.
  * **Sub-0.2 s Light Curve Smoothing (`smooth_all_counts_for_deadtime`)**: When `tbin < 0.2 s`, $R_{\text{all}}$ is evaluated from an $O(N)$ prefix-sum 1-second running mean of all-main-anode counts (matching `laxpcl1.f`).
  * **Frame-Gap Compensation**: Per-bin `fracexp = 1.0 - gap_fracs[i]` scales raw bin counts when `gap_fracs[i] < 0.90`.
* **Default Spectral Channels (`default_spectrum_channels`)**: `512` for `Lxp10`/`Lxp30` and `256` for `Lxp20` (avoiding even-odd ADC binning artifacts; `1024` also supported).
* **Output Formats**:
  * **Format B (per-unit, TIFR/LaxpcSoft style)**: `{prefix}level2.evn` (`LAXPC_MODE_EA` + `GTI`), `{prefix}level2.spec` (`SPECTRUM` v1.2.0), `{prefix}level2.lcurv` (`RATE` + `GTI`), `{prefix}level2.gtis` (`GTI`), plus optional ASCII `{prefix}level2.pha` / `.lc`.
  * **Format A (combined multi-unit)**: `level2.event.fits` (`LAXPC_EVENTS` with `TIME`, `Channel`, `Layer`, `LAXPC_No`, plus `GTI`). Uses a **uniform channel scale** across all units (`config.spectrum_channels.unwrap_or(1024)` with `TLMIN2 = 0`, `TLMAX2 = fa_channels - 1`), deterministic tie-broken sorting across threads, and complete `MJDREFI=55197`, `MJDREFF=0.0`, `MJDREF=55197.0`, `RA_NOM`/`DEC_NOM`, and `TIMESYS`/`CLOCKCOR` on both HDU 0 and HDU 1 so `astrosat bary` runs directly without `preplaxbary`.
  * **Orbit File**: `AstroSat.orb` (`ORBIT` HDU at 1-s cadence interpolated from MKF `POSX..VELZ`).

### 1.6 Background Modeling, Gain Alignment & RMF Recommendation (`background.rs`)
Implemented in [`background.rs`](../crates/astrosat-laxpc/src/background.rs):
* **2D B-Spline + 5-Harmonic Diurnal Model (`BackgroundModel2D`)**: Parses Fortran sequential unformatted `backfit.lxp?0*` files directly with checked dimensions (`read_fortran_records`), evaluates the tensor-product cubic (`k=4`) B-spline `eval_bspline_2d` in `(EarthLon, EarthLat)` using wrap-safe `interp_lon_deg`, adds 5 diurnal harmonics (`DEFAULT_DIURNAL_PERIOD_SEC = 84_496.0 s`, or file header `per`), and clamps to `mean_rate ± 40.0` cps (returning `(mean_rate - 40.0).max(0.0)` for `lon > 252°`). Automatically substitutes `lxp10`/`lxp1` $\to$ `lxp20`/`lxp30` per unit for both `backfit.*` and `backlxp*.pha`, and writes `{prefix}level2.lcbk` from the 2D model even if `background_pha_path` is omitted.
* **Gain Tracking & Background Catalog (`GainTable`, `HvJumpTable`, `BackgroundCatalog`)**:
  * Parses `gaina8.dat` (A8 30/60 keV peaks vs days since launch) and `lxhv.dat`/`lxhv.data` (HV jump epochs), applying one-sided interpolation within `±1.0` day of HV jumps, LAXPC30 exponential extrapolation (`ndays > 883`), and paired LAXPC20 linear drift (`source_day > 550 && background_day > 550` / `> 2977`) inside `compute_gain_shift` (matching `backshiftv3.f:298-307`).
  * `BackgroundCatalog` resolves the exact background 30 keV calibration peak channel (`chb`) from `backshiftv3.inp` (with a built-in TIFR background catalog fallback) so epochs near HV jumps match `backshiftv3.e` exactly; custom backgrounds can also pass `--bkg-peak-ch` / `background_peak_ch: Option<f64>`.
* **Spectral & Light-Curve Shifting (`shift_background_spectrum`, `apply_lightcurve_gain_shift`)**: Rebins background rates across channel boundaries using unit quadratic non-linearity `gain_shift_e2` (`-1.98e-3`, `-2.67e-3`, `-1.97e-3`) and `lxp{10,20,30}.channel` (plus epoch-shifted `lxp10cshm{141,175,197,212}v2.2.channel` maps for post-day-917 LAXPC10) energy edges (`ChannelMap`), scales by secular polynomials (`secular_background_scale`), writes `{prefix}level2back_shifted.spec` and `{prefix}level2.lcbk`, and warns on `crossed_hv_jumps`.
* **Automatic RMF Recommendation (`recommend_rmf`)**: Selects the nominal/gain-shifted TIFR `.rmf` (e.g. `lx10cshm01v1.1.rmf`, `lx20v1.0.rmf`, `lx30rho90cshp04v1.0.rmf`), including `L1` (top layer), `SE` (single events), and LAXPC30 gas-density `rho{irho:02}` (`estimate_lxp30_density`).

---

## 2. Configuration & Programmatic Entrypoints

`PipelineConfig` ([`pipeline.rs`](../crates/astrosat-laxpc/src/pipeline.rs)) orchestrates automatic L1 orbit discovery (`discover_l1_orbits`) and parallel execution (`run_pipeline`):

```rust
let cfg = PipelineConfig {
    input_dir: Some(l1_root.into()),
    output_dir: out_dir.into(),
    units: LaxpcUnit::ALL.to_vec(), // [Lxp10, Lxp20, Lxp30]
    tbin: 1.0,
    deadtime_sec: 43.0e-6,
    spectrum_channels: None,        // defaults: 512 (Lxp10/30), 256 (Lxp20); Format A uses 1024
    caldb_dir: Some(caldb_dir.into()),
    write_format_b_events: true,
    write_format_a_combined: true,
    write_orbit_file: true,
    ..Default::default()
};
let result = astrosat_laxpc::run_pipeline(&cfg)?;
```

---

## 3. Verification Notes: Crab ObsID `9000000252` (Orbit `01509`)

The public AstroSat Crab pulsar dataset (`20160107_G02_010T01_9000000252_level1`, 2016-01-07 / MJD `57394`, target ICRS `RA = 83.6331°`, `DEC = 22.0145°`) is exercised end-to-end in [`crates/astrosat-laxpc/tests/crab_integration.rs`](../crates/astrosat-laxpc/tests/crab_integration.rs):

* **Official L2 Reference Parity (`crab_9000000252_parity_against_official_l2_reference`)**: Running all 3 units with `spectrum_channels = Some(1024)` and the official `AS1G02_010T01_9000000252lxp{1,2,3}EA_level2.gti` verifies:
  * Gross GTI exposure matches the official `.pha` `EXPOSURE` within `< 0.2%`.
  * Total main-anode raw counts match within `~0.3%` across all 3 units (`Lxp10`/`Lxp20`/`Lxp30`).
  * 1024-channel spectral shape achieves Pearson correlation $r = 1.000000$ (`0.99999997`) against the official L2 `.pha` across the Crab continuum and Xe K edge.
* **Auto-GTI, Format B, Format A & Barycentric Check (`crab_9000000252_full_multi_unit_auto_gti_format_b_and_format_a`)**:
  * Verifies `net_exposure_sec > 3000 s` and Crab on-source deadtime-corrected count rates across all 3 units (`3379` / `3314` / `3005` cps), plus `25,824,880` combined events in `level2.event.fits`.
  * Validates every Format B (`.spec`, `.lcurv`, `.gtis`, `.evn`), Format A (`level2.event.fits`), `AstroSat.orb`, and barycentered `.evn` / `level2.event.fits` product with both `astrosat validate --strict` (`astrosat_fits::validate_ogip`) and NASA HEASARC `fitsverify 4.22` (`0 warnings, 0 errors`).
