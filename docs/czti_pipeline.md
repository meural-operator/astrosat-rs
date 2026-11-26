# AstroSat CZTI Pipeline (`astrosat-czti`)

`astrosat-czti` is a pure-Rust, `rayon`-parallelized implementation of all three stages of the ISRO AstroSat CZTI Level-1 to Level-2 scientific pipeline (matching CZTI Pipeline v3.0.1). It replaces the 16 standalone C++/Fortran/Python executables (`cztscience2event`, `cztbunchclean`, `cztpha2energy`, `cztgtigen`, `cztgaas`, `cztdatasel`, `cztpixclean`, `cztnoiseclean`, `cztevtencor`, `cztevtclean`, `cztflagbadpix`, `cztdpigen`, `cztimage`, `cztbindata`, `cztrspgen`, `cztaddspec`) with an in-memory, quadrant-parallel Rust pipeline orchestrated by [`run_pipeline`](../crates/astrosat-czti/src/pipeline.rs) and standalone file-level stage entrypoints in [`stages.rs`](../crates/astrosat-czti/src/stages.rs).

---

## 1. Detector Geometry & Physical Constants

Defined in [`crates/astrosat-czti/src/types.rs`](../crates/astrosat-czti/src/types.rs):

| Constant | Value | Meaning |
| :--- | :--- | :--- |
| `NUM_QUADS` | `4` | Quadrants `Q0`, `Q1`, `Q2`, `Q3` |
| `NUM_DET_PER_QUAD` | `16` | 4x4 CZT detector modules per quadrant |
| `NUM_PIX_PER_DET` | `256` | 16x16 pixels per detector module |
| `PIXELS_PER_QUAD` | `4096` | Total pixels per quadrant (`64 x 64`, `QUAD_DIM = 64`) |
| `FULL_DIM` | `128` | Full-plane pixel edge (`128 x 128`) |
| `PI_CHANNELS` | `512` | Output PI energy channels (`5.0..=261.0` keV, `0.5` keV/ch) |
| `VETOSPEC_SIZE` / `CZTSPEC_SIZE` | `256` / `512` | Veto and SSM CZT spectrum channel counts |
| `NUM_TEMP_BINS` | `5` | CALDB temperature bins (`0, 5, 10, 15, 20` °C) |
| `MASK_HEIGHT_MM` | `478.0` mm | Coded-mask height above detector plane |
| `XLDET_MM` | `39.0` mm | Detector module active side length for WCS scale |
| `MAX_FOV_DEG` | `5.8`° | Coded-mask FOV half-angle limit (`|tx|, |ty| <= 5.8°`) |
| `GAIN_DRIFT_EPOCH_SEC` | `181_094_400.0` s | Secular gain drift reference epoch (AstroSat UTC s) |
| `CZT_TICK_HZ` | `50_000.0` Hz | Onboard event fractional clock (`20 us` per tick) |

### Coordinate Mappings
* **`detid_pixid_to_detxy(quad_id, det_id, pix_id)`**:
  `dx = det_id % 4`, `dy = det_id / 4`, `px = pix_id % 16`, `py = pix_id / 16`, `detx = dx*16 + px`, `dety = dy*16 + py`; then `dety = 63 - dety` for `Q0`/`Q3` or `detx = 63 - detx` for `Q1`/`Q2`.
* **`detxy_to_locxy(quad_id, detx, dety)`** (128×128 full-plane coordinates):
  `Q0 -> (detx, 64 + dety)`, `Q1 -> (64 + detx, 64 + dety)`, `Q2 -> (64 + detx, dety)`, `Q3 -> (detx, dety)`.
* **Pixel Quality Flags (`pix_flag`)**: `GOOD = 0`, `SPEC_BAD = 1`, `FLICKERING = 2`, `NOISY = 3`, `DEAD = 4`.

---

## 2. Pipeline Stages

### Stage 1: Packet Decoding, Bunch Cleaning, and Energy Calibration

1. **`science2event` ([`science2event.rs`](../crates/astrosat-czti/src/science2event.rs))**
   * Reads `CZT_QUAD1..4` `DataArray` columns from Level-1 `modeM0` (and optional `modeSS`) FITS files in parallel across all 4 quadrants.
   * Decodes 2048-byte (`PACKET_SIZE_BYTES`, `1024` big-endian `u16` words) packets matching `SYNC1 = 0xF9A4`, `SYNC2 = 0x2BB1`.
   * Extracts 12-word `PacketFrameHeader` and 8/24-word `FebHeader` from `pkt_no == 0`, assembles contiguous frames (`assemble_frames`), and extracts:
     * **3-word normal events** (`mode_id in {0, 1, 4, 5}`, `eventsize = 3`): `cztntick` (`20 us`), 10-bit `PHA`, `DetID`, `pixID`, 7-bit `veto`, 1-bit `alpha`.
     * **2-word reduced events** (`mode_id in {2, 3, 6, 7}`, `eventsize = 2`).
     * **Onboard bunch reports** (`bunch_sync == 7` and `boot_page_no == 2`): decodes 3-word triplets into `BunchTable` (`Time`, `Time_dfs`, `Time_dsl`, `NumEvent`, `DetId1..4`, `RowNum`, `DetId/PixId_{fevt,sevt,levt}`), keeping event columns strictly synchronized even if a frame ends mid-bunch.
     * **`VETOSPECTRUM`** (`256`-ch), **`SSM Data`** (`data_id in 8..=11`, `CZTSpec`, `CZTSpecVeto`, `CZTSpecAlpha`), and **`TEMP`** (16 module temperatures per quadrant, cleaned via `clean_junk_temperatures`).
   * Converts `cztseccnt` $\to$ AstroSat UTC `Time` via `CztiTct` (`-2.0 s` FEB lag, `mod 1048576.0 s`).
   * Writes `{prefix}.hdr` (`HEADER` HDU).

2. **`bunchclean` ([`bunchclean.rs`](../crates/astrosat-czti/src/bunchclean.rs))**
   * Excludes the first `400 ms` of every `100 s` CZT clock tick.
   * Identifies cosmic-ray bunches (`dt <= bunch_def_time_us`, default `30 us`, `length >= 3`) when onboard bunches are not already flagged, retaining single (`evt_index == 1`) and double (`evt_index == 2`) events.
   * Applies post-bunch deadtime (`BunchCleanConfig`):
     * Default mode: `skip_t1 = 0.0 s`, `skip_t2 = 1e-3 s` (dominant module when `length < bunch_length_thresh = 20`), `skip_t3 = 1e-3 s` (quadrant-wide when `length >= 20`).
     * Generalized v3.0 preset (`BunchCleanConfig::generalized_compat()`): `skip_t1 = skip_t2 = skip_t3 = 0.0`, setting `NCOMPAT = 1` in output headers.
   * Computes 1-s binned `QuadLivetime` (`livetime_binsize = 1.0`) and per-pixel fractional exposure (`Vec<f64>`, 4096 pixels).
   * Writes `{prefix}_bc_bunch.fits` and `{prefix}_bc_livetime.fits`.

3. **`pha2energy` ([`pha2energy.rs`](../crates/astrosat-czti/src/pha2energy.rs))**
   * Linearly interpolates module temperatures from `modeSS` `TEMP` (`QuadTempSeries`, defaulting to `10.0 °C` if absent; standalone `astrosat czti pha2energy` supports both L2 `TEMP` extensions and L1 `modeSS` decoding when `--tct` is provided).
   * Selects the nearest of the 5 CALDB `GAIN` temperature bins (`temp_to_bin`: `<2.5`, `[2.5, 7.5)`, `[7.5, 12.5)`, `[12.5, 17.5)`, `>=17.5 °C`).
   * Evaluates $E = \text{gain} \times \text{PHA} + \text{offset}$ (`0.0` if `PHA == 0`) in parallel via `rayon`, and maps $E$ to 512-channel `PI` via binary search over CALDB `EBOUNDS`.
   * Writes `{prefix}_bc.evt`.

---

### Stage 2: GTI, Aspect, Noise Cleaning, and Energy Correction

4. **`gtigen` ([`gtigen.rs`](../crates/astrosat-czti/src/gtigen.rs))**
   * Evaluates `MkfTable` rows against `mkfThresholds.txt` rules (or `default_mkf_thresholds()`: `ELV >= 5.0`, `ANG_OFFSET <= 0.03`, `Q1..Q4_MODEID == 0`, `Q1_VetoCounter <= 450`, `Q2_VetoCounter <= 500`, `Q3_VetoCounter <= 500`, `Q4_VetoCounter <= 480`).
   * Intersects MKF GTIs with M0 `GTI` and `Q0_GTI..Q3_GTI`, subtracts optional `BTI`, intersects optional `user_gti_file`, and drops intervals `< min_exposure` (`1.05 s`).
   * Writes `{prefix}.gti` containing `GTI` and `Q0_GTI..Q3_GTI`.

5. **`gaas` ([`gaas.rs`](../crates/astrosat-czti/src/gaas.rs))**
   * Interpolates MKF Roll/Pitch/Yaw `(RA, DEC)` at 1-s cadence across GTI (`interpolate_rpy`), transforms detector pointing `(0, 0, 1)` and twist reference `(0, 1, 0)` through CALDB `TELDEF` (`ALIGNM11..33`) into ICS (`det_to_inertial`), and computes average `(RA, DEC, TWIST)` per quadrant (`calculate_ra_dec_twist_deg`).
   * Writes `{prefix}.aspect_Q0` through `{prefix}.aspect_Q3`.

6. **`datasel` ([`datasel.rs`](../crates/astrosat-czti/src/datasel.rs))**
   * Filters each quadrant's events in parallel against `COMMON` (`GTI`) or `QUAD` (`Q{q}_GTI`) intervals (`GtiType`), updating `TSTART`, `TSTOP`, and `EXPOSURE`.

7. **Pixel & Noise Cleaning (`EventSelectionAlgo` — parallel across quadrants)**
   * **`EventSelectionAlgo::Default` ([`pixclean.rs`](../crates/astrosat-czti/src/pixclean.rs))**:
     * Separates single vs double (`dt <= 25 us`) events with 30 us companion verification.
     * Accumulates 4096-pixel quadrant DPH and runs iterative `nsigma = 5.0` sigma-clipping (`find_bad_pixels`), marking `dph == 0` as `DEAD (4)` and outliers as `NOISY (3)` on top of CALDB `BADPIX` (retaining `GOOD (0)` and `SPEC_BAD (1)`).
     * Applies 1-s (`det_tbinsize = 1.0`, `pix_tbinsize = 1.0`) burst cuts: module counts `> 25` (`det_count_thresh`) and pixel counts `> 2` (`pix_count_thresh`), updating livetime (`re_update_livetime`) and per-pixel exposure.
   * **`EventSelectionAlgo::Generalized` ([`noiseclean.rs`](../crates/astrosat-czti/src/noiseclean.rs))**:
     * `noisypixclean_quad`: Applies pixel-area normalization (`2.28^2` edge vs `2.46^2` interior), marks `dph < dead_pix_threshold` (`10`) as `DEAD`, runs two-sided `5.0`-sigma clipping separately for `GOOD (0)` and `SPEC_BAD (1)` pixels, and resets `pixel_exposure = 1.0` for `flag < 2`.
     * `datasel_quad` $\to$ `superbunchclean_quad`: Merges consecutive saturated onboard bunches (`merge_split_bunches`, `num_event == 63` within `30 us`), flags super-bunches (`size >= 100` or `>= 4` heavy bunches of `size >= 15` within `1 ms` or `10 ms`), evaluates weighted spatial hotness $m_a m_b / d$ with a `3x3` neighbor double-event guard in `dph_time = 0.1 s` windows, and subtracts resulting BTI intervals from both livetime and `Q*_GTI`.
     * `heavybunchclean_quad` / `heavybunchclean_quad_with_exposure`: Merges saturated bunches, drops events in the `25x26` detector box around pixel `119` within `heavy_bunch_time_det = 0.005 s`, rejects low-energy (`< LLD_keV + 10.0 keV`) pairs within `<= 2 ms` over `heavy_bunch_time = 0.25 s`, removes small-bunch tails (`small_bunch_time_det = 0.001 s`), and updates per-pixel exposure.
     * `flickpixclean_quad` + `pix_threshold_fix`: Evaluates `< 50 keV` and `>= 50 keV` normalized light curves (`pix_lc / lc`) across `100.0 s` bins with adaptive `(sigma, num_times)` thresholds (`(5, 1)` / `(6, 2)` / `(8, 2)`) plus 1-s/10-s burst suppression, flagging flickering pixels as `FLICKERING (2)`.
     * `eventsep_quad`: Separates singles from Compton doubles within `double_time_sec = 30e-6 s` (`30 us`, matching `czteventsep.c`), optionally emitting `*_double.evt`.

8. **`evtencor` ([`evtencor.rs`](../crates/astrosat-czti/src/evtencor.rs)) & `evtclean` ([`evtclean.rs`](../crates/astrosat-czti/src/evtclean.rs))**
   * Applies per-pixel CALDB `GAINCORPAR` (where `GAINCORFLAG == 1`) and secular gain drift:
     $$\text{tgaincorpar} = 0.994 + \frac{0.00014925}{30.8955} \cdot \frac{\text{TSTART} - 181094400.0}{86400.0}$$
     with deterministic Wichmann-Hill PRNG sub-bin PI randomization (`rng_seed = (5, 11, 17)`, moduli `30269, 30307, 30319` matching `genurandt_c`), `[LLD, ULD]` cuts, and an idempotent `GAINCOR = 'YES'` guard.
   * `evtclean_quad` filters by `allowed_alpha = [0]` and `veto_range = (0, 0)` (and optional `--badpix` keeping `flag < pix_flag::NOISY`).
   * Writes `{prefix}_{quad|common}_badpix.fits`, `{prefix}_{quad|common}_livetime.fits`, and `{prefix}_{quad|common}_clean.evt`.

---

### Stage 3: High-Level Scientific Products & Co-Addition

9. **`dpigen` ([`dpigen.rs`](../crates/astrosat-czti/src/dpigen.rs))**
   * Builds 64×64 per-quadrant and 128×128 full-plane integer counts (`DPH`) for pixels with `flag <= badpix_threshold` (`0`), plus normalized CALDB `EFF_AREA`-corrected `DPI`.
   * Writes `{prefix}_{quad|common}_clean.dph` and `{prefix}_{quad|common}_clean.dpi`.

10. **`image` ([`image.rs`](../crates/astrosat-czti/src/image.rs))**
    * **`ImageAlgorithm::BalCC`** (default): Balanced Cross-Correlation using CALDB `SHADOWLIB` (mapping `(tx_i, ty_i)` to `shadows[tx_i * grid_dim + ty_i]` with WCS scale `6.0 / (grid_dim - 1)` = `0.10 deg/pix` for 61×61 `[-3°, +3°]` grids; falls back to shifted 64×64 `MASK_PATTERN`), parallelized across sky pixels:
      $$I_{jk} = \frac{\sum_{\text{DPI}\neq 0} s_p\,\text{DPI}_p}{\sum_{\text{DPI}\neq 0} s_p} - \frac{\sum_{\text{DPI}\neq 0} (1 - s_p)\,\text{DPI}_p}{\sum_{\text{DPI}\neq 0} (1 - s_p)}.$$
    * **`ImageAlgorithm::Fft`**: 2D `rustfft` cross-correlation of oversampled (`oversampling_factor = 16`) `DPI` and balanced mask $2M - 1$, applying `yshift = 1.0` Fourier phase ramp for `Q2`/`Q3` and outputting complex magnitude.
    * Computes both ICRS (`RA---TAN`/`DEC--TAN`, `CD*_*`) and local camera (`TX---TAN`/`TY---TAN`, `CD*_*B`) WCS via `CztiWcs`.
    * Writes `{prefix}_{quad|common}_clean.img_Q{0..3}` (and combined `{prefix}_common_clean.img` via `combine_sky_images` in `GtiType::Common` mode).

11. **`bindata`, `cztrspgen` & `cztaddspec` ([`bindata.rs`](../crates/astrosat-czti/src/bindata.rs), [`stages.rs`](../crates/astrosat-czti/src/stages.rs))**
    * Projects source `(RA, DEC)` (`source_coords: Option<(f64, f64)>`, unambiguously supporting `(0.0, 0.0)`, or falling back to MKF `RA_PNT`/`DEC_PNT` when unset) into CZTI camera coordinates `(tx, ty)` (`radec_to_czti_tx_ty`) and enforces `|tx|, |ty| <= 5.8°` (`validate_fov`).
    * Applies measured per-module boresight offsets (`DETWISE_TX_SHIFT`, `DETWISE_TY_SHIFT`) plus `delta_tx`/`delta_ty` and per-module `ADDITIONAL_LLD` floors (`apply_additional_lld`), and ray-traces the 20-micron (`MASK_CELL_MM = 0.02`) compressed CALDB `MASK_OVERSAMPLED` at height `478.0 mm` (`2.46 mm` interior / `2.28 mm` edge pixel pitch, `2.5 mm` inter-module gap) to compute `open_frac` (including 5×5 sub-pixel averaging `compute_5x5_avg_open_frac`) and `raw_weights = 2 * open_frac - 1`.
    * Computes module-wise `BKGSCAL`-renormalized energy-dependent weights (`calculate_renormalized_weights`, requiring `> min_active_pix_per_module = 100` active pixels per 16×16 module) guaranteeing exact background cancellation $\sum w_i B_i T_i = 0$.
    * Accumulates weighted spectra ($\sum w$, $\sigma = \sqrt{\sum w^2}$, CALDB `SYSERR`, exact `QUALITY = 5` for empty/`< 34`/`> 390` channels and SNR-based logarithmic `GROUPING` via `compute_spectrum_quality_and_grouping`) and livetime-corrected light curves (`FRACEXP > min_fracexp = 0.05`).
    * Generates per-quadrant OGIP 1.3.0 `.rsp` (`generate_quad_response`, `write_response_file` with `lo_thresh = 1e-10`) using CALDB `PIXRESP`, `RESP_PAR` (nominal temperature bin 2 = 10 °C), `EFF_AREA`, `EFF_AREACOR` (`AREACOR1`, `AREACOR2`), `[LLD, ULD]`, and Tantalum 0.5 mm mask transmission (`absco_ta`, `ta_mask_transmission`), parallelized across all 16 modules.
    * In `GtiType::Common` mode, combines quadrants via `combine_spectra`, `combine_responses`, and `combine_lightcurves`.
    * **`cztaddspec` (`stages::run_addspec` / `astrosat czti addspec`)**: Co-adds arbitrary lists of `.pha` spectra (exact counts sum, quadrature `STAT_ERR`, counts-weighted `SYS_ERR`, averaged or summed exposure via `--add-exposure`) and optional exposure-weighted `.rsp` matrices.

---

## 3. Configuration & Programmatic Entrypoints

`CztiPipelineConfig` ([`pipeline.rs`](../crates/astrosat-czti/src/pipeline.rs)) drives `run_pipeline(&config) -> Result<CztiPipelineSummary>`:

```rust
let mut cfg = CztiPipelineConfig::new(
    l1_m0_file,   // *cztM0_level1.fits
    Some(l1_ss),  // optional *cztSS_level1.fits for TEMP
    tct_file,     // *_level1.tct
    mkf_file,     // *.mkf
    caldb_dir,    // root containing data/as1/czti/caldb.indx
    output_dir,
);
cfg.entry_stage = 1; // 1..=3
cfg.exit_stage = 3;  // 1..=3
cfg.gti_type = GtiType::Quad; // Quad or Common
cfg.evt_sel_algo = EventSelectionAlgo::Default; // Default or Generalized
let summary = astrosat_czti::run_pipeline(&cfg)?;
```

---

## 4. CALDB Requirements (`CztiCaldb`)

`CztiCaldb::open` resolves all calibration files automatically from `caldb.indx` (`data/as1/czti/caldb.indx`), transparently decompressing gzipped `.fits` files (`open_fits_maybe_gz`, e.g., `PIXRESP` and `SHADOWLIB`) via `flate2`:

| CIF Codename | Example File in CZTI CALDB `20221209` | Pipeline Usage |
| :--- | :--- | :--- |
| `EBOUNDS` | `AS1cztebounds20150707v01.fits` | 512-ch energy $\leftrightarrow$ PI grid (`pha2energy`, `.rsp`) |
| `GAIN` | `AS1cztgain20150707v01.fits` | 5-temperature gain/offset tables (`GAIN_00..GAIN_20`) |
| `TELDEF` | `AS1czt{0..3}20151005v01.teldef` | 3×3 detector-to-spacecraft alignment (`ALIGNM11..33`) |
| `BADPIX` | `AS1cztbadpix20220720v03.fits` | Baseline bad-pixel map (`Q0..Q3`) |
| `GAINCOR` | `AS1cztgaincor20220720v01.fits` | In-flight pixel gain factors (`GAINCORPAR`, `GAINCORFLAG`) |
| `LLD` / `ULD` | `AS1cztlld20220720v04.fits`, `AS1cztuld20220720v01.fits` | Per-pixel 5-temp (`5I`) PI thresholds (`evtencor`, `bindata`, `.rsp`) |
| `EFF_AREA` | `AS1czteff_area20220805v04.fits` | Pixel effective area (`dpigen`, `.rsp`) |
| `EFF_AREACOR` | `AS1czteff_areacor20220914v01.fits` | Module-PI (`AREACOR1`) & energy (`AREACOR2`) corrections |
| `BKGSCAL` | `AS1cztbkgscal20210618v01.fits` | 4096×512 background normalization weights (`bindata`) |
| `SYSERR` | `AS1cztsyserr20220921v01.fits` | Per-channel fractional systematic error (`.pha`) |
| `MASK_PATTERN` | `AS1cztmask_pattern20150707v01.fits` | 64×64 mask pattern (`fft` / fallback `balcc`) |
| `MASK_OVERSAMPLED` | `AS1cztmask_oversampled20150707v01.fits` | Compressed 20 $\mu$m bitmask (`bindata` ray-tracing) |
| `SHADOWLIB` | `AS1cztshadowlib20221004v01.fits` | Precomputed 61×61 sky shadow templates (`balcc`, gzipped FITS supported) |
| `PIXRESP` / `RESP_PAR` | `AS1cztpixresp20150226v03.fits`, `AS1cztresp_par20150922v01.fits` | Pixel response groups & redistribution matrices (`.rsp`, gzipped FITS supported) |

---

## 5. Exact FITS Output Products

Assuming `prefix` is derived from `l1_m0_file` (or `output_prefix`) and `gti_word` is `quad` or `common`:

| Product | Extensions / Format | Produced By |
| :--- | :--- | :--- |
| `{prefix}.hdr` | `PRIMARY`, `HEADER` (30 frame/FEB telemetry columns) | Stage 1 (`science2event`) |
| `{prefix}_bc_bunch.fits` | `PRIMARY`, `Q0..Q3` (bunch records) | Stage 1 (`bunchclean`) |
| `{prefix}_bc_livetime.fits` | `PRIMARY`, `Q0..Q3` (`TIME`, `FRACEXP`) | Stage 1 (`bunchclean`) |
| `{prefix}_bc.evt` | `PRIMARY`, `Q0..Q3`, `VETOSPECTRUM`, `SSM Data`, `TEMP`, `GTI`, `Q0_GTI..Q3_GTI`, `EXPOSURE` | Stage 1 (`pha2energy`) |
| `{prefix}.gti` | `PRIMARY`, `GTI`, `Q0_GTI..Q3_GTI` | Stage 2 (`gtigen`) |
| `{prefix}.aspect_Q0..3` | `PRIMARY`, `ASPECT` (`TIME`, `Nx..Nzt`, `RA`, `DEC`, `TWIST`) | Stage 2 (`gaas`) |
| `{prefix}_{gti_word}_badpix.fits` | `PRIMARY`, `Q0..Q3` (`PIX_FLAG`), `BADPIX` (128×128 `U8` IMAGE) | Stage 2 (`pixclean`/`noiseclean`) |
| `{prefix}_{gti_word}_livetime.fits` | `PRIMARY`, `Q0..Q3` (`TIME`, `FRACEXP`) | Stage 2 (`pixclean`/`noiseclean`) |
| `{prefix}_{gti_word}_clean.evt` | Same 13 HDUs as `_bc.evt` with `GAINCOR='YES'`, `LLDCUT`, `ULDCUT` | Stage 2 (`evtencor`/`evtclean`) |
| `{prefix}_{gti_word}_clean.dph` | `PRIMARY`, `Q0..Q3` (64×64 `I32`), `FULL_DPH` (128×128 `I32`) | Stage 3 (`dpigen`) |
| `{prefix}_{gti_word}_clean.dpi` | `PRIMARY`, `Q0..Q3` (64×64 `F32`), `FULL_DPI` (128×128 `F32`) | Stage 3 (`dpigen`) |
| `{prefix}_{gti_word}_clean.img_Q0..3` | `PRIMARY`, `Q{q}` (`F32` sky IMAGE + ICRS & `TX---TAN`/`TY---TAN` WCS) | Stage 3 (`image`) |
| `{prefix}_{gti_word}_clean_Q{0..3}.pha` | `PRIMARY`, `SPECTRUM` (OGIP 1.2.0, 512 PI channels, `QUALITY` & `GROUPING`) | Stage 3 (`bindata`) |
| `{prefix}_{gti_word}_clean_Q{0..3}.lc` | `PRIMARY`, `LIGHTCURVE` (`TIME`, `RATE`, `ERROR`, `FRACEXP`) | Stage 3 (`bindata`) |
| `{prefix}_{gti_word}_clean_Q{0..3}.rsp` | `PRIMARY`, `EBOUNDS`, `MATRIX` (OGIP 1.3.0, sparse `lo_thresh=1e-10`) | Stage 3 (`bindata`) |
| `{prefix}_{gti_word}_clean_weight.evt` | Optional (`generate_eventfile=true`): adds `OPEN_FRACTION`, `WEIGHT` | Stage 3 (`bindata`) |
