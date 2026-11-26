# `astrosat-rs` Architecture Overview

`astrosat-rs` is a pure-Rust, `rayon`-parallelized redesign of the ISRO AstroSat **CZTI** (Cadmium Zinc Telluride Imager) and **LAXPC** (Large Area X-ray Proportional Counter) Level-1 to Level-2 data processing pipelines and Solar System Barycentric correction toolchain.

Unlike the legacy ISRO software stack—which requires compiling C++ against bundled `cfitsio`, `glog_modified-0.3.2`, and `pil` (`InstallLibs`), compiling Fortran-77/90 monoliths with platform-specific `recl=2880` and byte-swapping assumptions (`laxpcl1.f`, `backshiftv3.f`), and initializing HEASoft (`$HEADAS/headas-init.sh`) for `as1bary`—`astrosat-rs` builds into a self-contained pure-Rust workspace with **zero external C, Fortran, or HEASoft runtime dependencies**.

---

## 1. Workspace Crate Graph

Defined in [`Cargo.toml`](../Cargo.toml) (`edition = "2021"`, `resolver = "2"`, release profile `opt-level = 3`, `lto = "thin"`, `codegen-units = 1`):

```text
                       ┌────────────────────┐
                       │    astrosat-cli    │  Unified `astrosat` CLI + TOML/JSON config
                       └────┬───┬───────┬───┘
                            │   │       │
        ┌───────────────────┘   │       └───────────────────┐
        ▼                       ▼                           ▼
┌───────────────┐       ┌───────────────┐           ┌───────────────┐
│ astrosat-czti │       │ astrosat-laxpc│           │ astrosat-bench│
│ (3-stage CZTI │       │ (LaxpcSoft    │           │ (throughput & │
│  L1 -> L2)    │       │  L1 -> L2)    │           │  RSS suite)   │
└───────┬───────┘       └───────┬───────┘           └───────┬───────┘
        │                       │                           │
        └───────────┬───────────┴───────────────────────────┘
                    ▼
            ┌───────────────┐
            │ astrosat-core │  Timing, TCT, GTI, MKF, Orbit, Bary, CALDB CIF
            └───────┬───────┘
                    ▼
            ┌───────────────┐
            │ astrosat-fits │  Zero-copy memmap2 FITS 4.0 + CHECKSUM + OGIP
            └───────────────┘
```

| Crate | Path | Responsibility | Key External Crates |
| :--- | :--- | :--- | :--- |
| `astrosat-fits` | [`crates/astrosat-fits/`](../crates/astrosat-fits/Cargo.toml) | Memory-mapped FITS 4.0 parser/writer, Seaman/Pence checksums, `BINTABLE`/`IMAGE`, OGIP product builders & validator | `memmap2`, `thiserror` |
| `astrosat-core` | [`crates/astrosat-core/`](../crates/astrosat-core/Cargo.toml) | AstroSat `MJDREFI=55197` timing, leap seconds, CZTI/LAXPC TCT, `GtiSet`, `MkfTable`, `OrbitTable`, `bary`, `CaldbIndex` | `astrosat-fits`, `rayon`, `chrono`, `serde`, `thiserror` |
| `astrosat-czti` | [`crates/astrosat-czti/`](../crates/astrosat-czti/Cargo.toml) | Full 3-stage CZTI Level-1 -> Level-2 pipeline (`science2event` through `bindata` / `.rsp` + `addspec`) | `astrosat-fits`, `astrosat-core`, `rayon`, `rustfft`, `num-complex`, `flate2`, `serde`, `thiserror` |
| `astrosat-laxpc` | [`crates/astrosat-laxpc/`](../crates/astrosat-laxpc/Cargo.toml) | Full LAXPC10/20/30 Level-1 -> Level-2 pipeline (`frames`, `ea`, `bb`, `gti`, `products`, `background`) | `astrosat-fits`, `astrosat-core`, `rayon`, `serde`, `serde_json`, `thiserror`, `tracing` |
| `astrosat-cli` | [`crates/astrosat-cli/`](../crates/astrosat-cli/Cargo.toml) | Unified standalone `astrosat` binary CLI (`czti`, `laxpc`, `bary`, `inspect`, `validate`, `timing`), TOML/JSON config loader | All library crates, `clap`, `anyhow`, `rayon`, `serde_json`, `tempfile`, `tracing`, `tracing-subscriber` |
| `astrosat-bench` | [`crates/astrosat-bench/`](../crates/astrosat-bench/Cargo.toml) | Reproducible latency, throughput, `rayon` scaling, and `/proc/self/status` `VmHWM` peak RSS benchmarking harness | All library crates, `clap`, `anyhow`, `rayon`, `serde_json`, `tempfile` |

---

## 2. Zero-Copy FITS 4.0 & OGIP Engine (`astrosat-fits`)

Implemented in [`crates/astrosat-fits/src/`](../crates/astrosat-fits/src/lib.rs), `astrosat-fits` replaces CFITSIO with direct big-endian byte-slice indexing over memory-mapped FITS files.

### 2.1 Memory-Mapped Reader & Streaming Writer (`io.rs`)
* **`FitsFile::open`**: Maps the FITS file read-only via `memmap2::Mmap` (`BackingStore::Mmap`), automatically falling back to `std::fs::read` (`BackingStore::Owned`) if `mmap` is unsupported by the underlying filesystem.
* **Block Alignment & HDU Indexing**: Scans 2880-byte (`FITS_BLOCK_SIZE`) FITS blocks once at open time, validating `SIMPLE = T` on Primary HDUs, valid `BITPIX` values (`8, 16, 32, 64, -32, -64`), non-negative `NAXIS`/`PCOUNT`/`GCOUNT` with checked arithmetic against file length, and `XTENSION` (`BINTABLE`, `IMAGE`, `TABLE`) on extensions.
* **Zero-Copy Views**: `HduRef::raw_data()`, `BinTableRef::row_bytes(row)`, `BinTableRef::cell_bytes(row, col_idx)`, and `ImageRef::raw_data()` borrow `&[u8]` slices directly from the OS page cache without heap allocation—critical when scanning 2048-byte CZTI packets or 2242-byte LAXPC frames.
* **Writer (`FitsWriter` / `write_fits`)**: Uses a 256 KiB `BufWriter<File>`, enforces 2880-byte block padding (`0u8` for binary/image HDUs, ASCII spaces `b' '` for `TABLE`), and optionally updates `DATASUM` and `CHECKSUM` on each HDU prior to serialization.

### 2.2 FITS 4.0 Header Cards & `CONTINUE` Long Strings (`card.rs`)
* Parses and formats 80-byte FITS cards (`Card`, `CardValue::{Logical, Integer, Float, Str, Empty}`), including legacy Fortran `D` exponent notation (e.g., `1.23D-04`) on input and escaped single quotes (`''`) on output.
* Automatically emits FITS 4.0 / OGIP 1.0 `CONTINUE` cards (splitting with trailing `&` without breaking `''` escape pairs) and injects `LONGSTRN = 'OGIP 1.0'` when any string value exceeds 68 bytes.

### 2.3 Seaman / Pence / Rots `CHECKSUM` & `DATASUM` (`checksum.rs`)
Implements the FITS 4.0 ones'-complement 32-bit big-endian checksum specification:
* `accumulate_checksum` / `datasum`: Sums 4-byte big-endian words with end-around carry folding (`(sum & 0xFFFF_FFFF) + (sum >> 32)`).
* `compute_hdu_checksums`: Writes `DATASUM`, temporarily sets `CHECKSUM = '0000000000000000'`, sums the 2880-byte-aligned header + `datasum`, and encodes `!total` via `encode_checksum` (16-char ASCII encoding avoiding punctuation ranges `0x3A..=0x40` and `0x5B..=0x60`, rotated right by 1 byte for column-12 alignment).
* `verify_hdu_checksum`: Independently validates `DATASUM` and checks that `header_sum + datasum` in ones'-complement equals `0xFFFF_FFFF` (or `0`).

### 2.4 `BINTABLE`, `IMAGE`, and OGIP Builders (`bintable.rs`, `image.rs`, `ogip.rs`)
* **`BinTableRef` / `BinTableBuilder`**: Supports fixed scalar and vector columns (`L`, `X`, `B`, `I`, `J`, `K`, `E`, `D`, `A`, `C`, `M`, e.g., `2242B`, `2048B`, `512E`, `16E`), `TZERO`/`TSCAL` scaling (preserving exact `i64` arithmetic for unsigned 16/32-bit `TZERO = 32768` / `2147483648` columns), and checked variable-length `P` (32-bit) / `Q` (64-bit) heap descriptors (`read_var_i16`, `read_var_f32` for sparse RMF `F_CHAN`, `N_CHAN`, `MATRIX`).
* **`ImageRef` / `ImageBuilder`**: Supports `Bitpix::{U8, I16, I32, I64, F32, F64}` (`8, 16, 32, 64, -32, -64`), `BZERO`/`BSCALE`, and `WcsParams` (`RA---TAN`, `DEC--TAN`, `CRPIX*`, `CRVAL*`, `CDELT*`, `CROTA2`, `RADESYS`, `EQUINOX`).
* **OGIP Science Builders (`ogip.rs`)**: `apply_ogip_meta`, `build_gti_hdu` (`GTI`), `build_spectrum_hdu` (OGIP 1.2.0 `SPECTRUM`, `COUNTS` or `RATE`), `build_lightcurve_hdu` (`RATE` with `TIMEDEL`, `TIMEPIXR = 0.5`, `FRACEXP`), `build_response_hdus` (OGIP 1.3.0 `MATRIX` + `EBOUNDS`), `build_orbit_hdu` (`ORBIT` with `TIME`, `X..Z`, `VX..VZ`), and `validate_ogip`.

---

## 3. Shared Science & Timing Engine (`astrosat-core`)

Implemented in [`crates/astrosat-core/src/`](../crates/astrosat-core/src/lib.rs), `astrosat-core` provides mission-wide AstroSat timing, calibration indexing, orbit interpolation, and barycentric physics.

### 3.1 AstroSat Time System & Leap Seconds (`time.rs`)
* **Mission Reference Epoch**: `ASTROSAT_MJDREFI = 55197`, `ASTROSAT_MJDREFF = 0.0` (`2010-01-01T00:00:00.000 UTC`), `TT_MINUS_TAI = 32.184` s, `SECS_PER_DAY = 86_400.0`, `JD_MINUS_MJD = 2_400_000.5`, `JD_J2000 = 2_451_545.0`.
* **Embedded `tai-utc.dat` Table (`LEAP_SECONDS_TABLE`)**: Covers every integer leap second from 1972-01-01 (`10.0` s) through 2026+ (`34.0` s at AstroSat epoch MJD `55197`, `35.0` s at `56109` / 2012-07-01, `36.0` s at `57204` / 2015-07-01 active at AstroSat's 2015-09-28 launch, and `37.0` s at `57754` / 2017-01-01).
* **Precision Time Scale Chain**:
  $$\text{UTC (AstroSat s)} \xrightarrow{+\,(\text{TAI}-\text{UTC})(\text{MJD}) + 32.184\text{ s}} \text{TT} \xrightarrow{\text{split }(\text{JD}_{\text{hi}},\,\text{JD}_{\text{lo}})} \text{TDB}$$
  Split Julian Date representation `(jd_hi, jd_lo)` preserves sub-nanosecond precision across multi-decade baselines. `tdb_minus_tt` evaluates all **127 exact terms** (`93 T0 + 27 T1 + 6 T2 + 1 T3`) of the Fairhead & Bretagnon (1990) analytical series matching `as1bary/ctatv.c` to `< 1 ns` across 1900–2100.

### 3.2 Time Correlation Table (`.tct`) Interpolation (`tct.rs`)
* **`CztiTct`**: Subtracts `CZTI_DEFAULT_FEB_LAG = 2.0` s from `cztseccnt`, unwraps the 20-bit counter across `CZTI_CLOCK_WRAP_SEC = 1_048_576.0` s ($2^{20}$ s), splits non-rollover onboard clock resets into monotonic `CztiSegment`s, and rejects non-finite inputs. Supports `rayon`-parallel `interpolate_batch` for `>= 16,384` timestamps.
* **`LaxpcTct`**: Filters out telemetry/clock glitch rows where `dt_sps / dt_ins` lies outside `(LAXPC_TCT_RATIO_MIN = 1.000, LAXPC_TCT_RATIO_MAX = 1.001]`. Provides $O(\log N)$ binary-search linear interpolation (`interpolate`) and `order`-point Neville polynomial interpolation (`interpolate_neville`, up to order 8), plus `rayon`-parallel `interpolate_batch`.

### 3.3 Good Time Interval Algebra (`gti.rs`) & MKF (`mkf.rs`)
* **`GtiSet`**: Automatically drops invalid intervals (`stop <= start`), sorts, and merges overlaps. Provides $O(N + M)$ `intersect`, `union`, `subtract` (for BTI or frame-gap removal), `merge_gaps(max_gap)`, `drop_short(min_exposure)`, `exposure()`, `overlap_duration(start, stop)` (for light-curve `FRACEXP`), $O(\log N)$ `contains`, and `rayon`-parallel `filter_mask` / `filter_indices` for `>= 16,384` events.
* **`MkfTable`**: Loads all scalar numeric columns from `.mkf` FITS tables, supports case-insensitive and multi-alias lookups (`column_any`), linear interpolation (`interpolate`), and predicate-based GTI extraction (`gti_where`).

### 3.4 Spacecraft Orbit & Solar System Barycentric Correction (`orbit.rs`, `bary.rs`)
* **`OrbitTable`**: Reads AstroSat `.orb` (`ORBIT` HDU, `X/Y/Z` and `VX/VY/VZ` or `POS_*/VEL_*`), auto-converting `km` / `km/s` to SI `m` / `m/s` (via `TUNIT` or orbital radius magnitude). `interpolate(t_sec)` uses the exact 3-point centered-difference 2nd-order Taylor expansion from `as1bary/scorbit.c` (with 120 s endpoint tolerance and strict NaN/bounds guards):
  $$\mathbf{a}_i = \frac{\mathbf{v}_{i+1} - \mathbf{v}_{i-1}}{t_{i+1} - t_{i-1}}, \qquad \mathbf{r}(t) = \mathbf{r}_i + \mathbf{v}_i\,\Delta t + \frac{1}{2}\mathbf{a}_i\,\Delta t^2, \qquad \mathbf{v}(t) = \mathbf{v}_i + \mathbf{a}_i\,\Delta t.$$
* **`BaryCalculator` (`bary.rs`)**: Implements exact `as1bary` relativistic and geometric corrections (`C_LIGHT = 299_792_458.0` m/s, `AU_M = 149_597_870_691.0` m, `GM_SUN_OVER_C3 = 4.925_490_947e-6` s, `EMRAT_DE405 = 81.300_56`):
  1. **Coordinate Frame**: `source_unit_vector_icrs` converts `(ra_deg, dec_deg)` in `CoordFrame::Icrs` or applies Standish's DE200/FK5 $\to$ DE405/ICRS rotation (`fk5_to_icrs`).
  2. **Einstein term**: `einstein_s = tdb_minus_tt(jd_tt_hi, jd_tt_lo)` (127-term Fairhead & Bretagnon series).
  3. **Earth & Spacecraft Roemer**: `roemer_earth_s = (r_earth . n) / c`, `roemer_sc_s = (r_sc . n) / c`.
  4. **Spacecraft Lorentz**: `lorentz_sc_s = (r_sc . v_earth) / c^2`.
  5. **Solar Shapiro delay**: `shapiro_s = 2 * GM_SUN_OVER_C3 * ln(1 + n . r_hat_sun_sc)`.
* **Dual Ephemeris Backends (`BaryEphemeris`)**:
  * `BaryEphemeris::Analytical`: Zero-external-file built-in analytical ephemeris combining Standish (1992) Table 2a + 34-term frequency-domain Hill EMB perturbations + ELP2000/Meeus lunar reflex (`~4–9 ms` / `< 10 ms` accuracy vs JPL DE405, `PLEPHEM = 'ASTROSAT-ANALYTICAL'`).
  * `BaryEphemeris::JplDe(Box<JplDeEphemeris>)`: Pure-Rust reader for **BOTH** official raw JPL binary ephemerides (`lnxp*.405`, DE200/DE405/DE421/DE430, auto-detecting big/little endianness) **AND** NASA HEASARC FITS-wrapped ephemerides (`JPLEPH.405` / `JPLEPH.200`, reading `DE1`/`DE2`/`DE3` HDUs), achieving **`<= 29.8 ns` (1 ULP of `f64`)** parity with compiled NASA/IUCAA `as1bary` (`PLEPHEM = 'JPL-DE<numde>'`).
* **Universal `barycenter_fits`**: Works seamlessly on CZTI, LAXPC Format B (`.evn`/`.lcurv`), and LAXPC Format A (`level2.event.fits`). Resolves `MJDREFI`/`MJDREFF` from HDU 1 first then HDU 0 (avoiding the `as1bary` `preplaxbary` failure mode), folds any pre-existing non-zero `TIMEZERO` into timestamps before setting `TIMEZERO = 0.0`, rejects already-barycentered files (`TIMESYS = 'TDB'` / `'TCB'`), updates `TIME` and GTI `START`/`STOP` columns in parallel, updates `TSTART`/`TSTOP` and `DATE-OBS`/`DATE-END` consistently across event/GTI/Primary HDUs, and recomputes `DATASUM`/`CHECKSUM`.

### 3.5 CALDB CIF Index (`caldb.rs`)
`CaldbIndex::open` reads OGIP `caldb.indx` (`CIF` HDU), probing standard roots (`caldb.indx`, `data/as1/{czti,laxpc}/caldb.indx`, `data/astrosat/{czti,laxpc}/caldb.indx`, `bcf/caldb.indx`). `query` and `query_best` filter on `CAL_QUAL == 0`, `TELESCOP`, `INSTRUME`, `CAL_CNAM`, `DETNAM` (preferring exact matches over `NONE`/`-`), and validity start `REF_TIME` (`CAL_VSD`/`CAL_VST`), returning the newest valid calibration entry.

---

## 4. `rayon` Data-Parallel Architecture

`astrosat-rs` exploits multi-core parallelism at the natural structural boundaries of AstroSat telemetry and event lists while maintaining **bit-for-bit deterministic** outputs across thread counts (`--threads 1` vs `--threads 16`):

1. **CZTI Quadrant & Pixel/Sky Parallelism (`astrosat-czti`)**:
   * **Stage 1**: All 4 quadrants (`Q0..Q3`) are decoded (`decode_quad_packets`), bunch-cleaned (`bunchclean_quad`), and energy-calibrated (`pha2energy_quad`, parallelized per event within each quadrant) concurrently via `rayon::par_iter`.
   * **Stage 2**: `datasel`, `pixclean` / `noiseclean` (`noisypixclean`, `superbunchclean`, `heavybunchclean`, `flickpixclean`, `eventsep`), `evtencor`, and `evtclean` execute concurrently across all 4 quadrants via `(0..NUM_QUADS).into_par_iter()`.
   * **Stage 3**: `generate_dph_dpi` runs per-quadrant; `balanced_cross_correlation` parallelizes across all `grid_dim * grid_dim` sky pixels; `fft_cross_correlation` parallelizes row/column `rustfft` transforms via `par_chunks_mut`; `bindata` attitude projection (`radec_to_czti_tx_ty`) and `raytrace_quad_shadow` parallelize across attitude samples and detector rows; `generate_quad_response` parallelizes `.rsp` matrix accumulation across the 16 detector modules (`0..16`) per quadrant with all CALDB arrays resident in memory (zero inner-loop FITS I/O).
2. **LAXPC Multi-Unit & Multi-Orbit Parallelism (`astrosat-laxpc`)**:
   * `run_pipeline` processes all requested units (`LAXPC10`, `LAXPC20`, `LAXPC30`) concurrently via `unit_orbits.par_iter()`.
   * Within each unit (`run_unit_internal`), all orbit FITS files (`modeEA` and `modeBB`) are memory-mapped, parsed, and frame-sorted concurrently via `orbits.par_iter()` / `bb_paths.par_iter()`.
   * Spectrum accumulation (`build_spectrum`), raw light-curve binning (`bin_events_raw`), 2D B-spline background evaluation, event UTC calibration (`build_calibrated_events`), and deterministic Format A time-sorting all execute in parallel across event chunks.
3. **Core Vectorized Operations (`astrosat-core`)**:
   * `CztiTct::interpolate_batch` and `LaxpcTct::interpolate_batch` switch to `par_iter()` for `>= 16,384` timestamps.
   * `GtiSet::filter_mask` and `GtiSet::filter_indices` parallelize for `>= 16,384` events.
   * `BaryCalculator::correct_times` parallelizes for `>= 8,192` timestamps.

---

## 5. Zero-Dependency Design & Conformance

* **No External Shared Libraries**: No `libcfitsio.so`, `libglog.so`, `libpil.a`, `libgfortran`, or `$HEADAS` runtime environment is required (linking only standard glibc/libm/libgcc).
* **No Subprocess File Hacking**: Legacy Fortran `call system('rm ...')` workarounds and fixed-size static array bounds (`parameter (nn=1000000, ns=9740000, nl=10000000)` in `laxpcl1.f`) are eliminated in favor of dynamically sized Rust `Vec`s and deterministic atomic file writes.
* **Verified Standards Conformance**: All generated FITS products include standard FITS 4.0 `DATASUM`/`CHECKSUM` cards and pass both `astrosat validate --strict` (`astrosat_fits::validate_ogip`) and NASA HEASARC `fitsverify 4.22` (`0 warnings, 0 errors`).
