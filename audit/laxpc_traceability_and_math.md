# LAXPC Traceability & Mathematical Verification Audit Report

**Scope**: Exhaustive ground-truth code and mathematical audit comparing legacy ISRO/TIFR **LaxpcSoft v3.4.6** (`data/ref_src/laxpc/{laxpcl1.f,backshiftv3.f,readme.laxpcl1,readme.backshiftv3,README}`) and Format A conventions against the pure-Rust implementation (`crates/astrosat-laxpc/src/`, `crates/astrosat-core/src/tct.rs`, `crates/astrosat-cli/src/laxpc_cmd.rs`).

**Audit Classification Legend**:
- `EXACT_MATCH` — Bit-level, algebraic, and algorithmic parity with legacy `LaxpcSoft v3.4.6` / Format A.
- `INTENTIONAL_IMPROVEMENT` — Deliberate precision, safety, OGIP/FITS compliance, or bug-fix improvement over legacy Fortran (documented with exact differences).
- `PARTIAL_GAP` — Legacy feature or optional sub-mode omitted or partially implemented (low/moderate impact, explicitly cataloged).
- `BUG/MATH_DRIFT` — Mathematical drift, data-routing bug, or CLI mismatch discovered during audit.

---

## 1. Executive Summary

- **`crates/astrosat-laxpc`**: Core telemetry decoding (`frames.rs`, `ea.rs`, `bb.rs`), TCT glitch rejection (`astrosat-core::tct::LaxpcTct`), deadtime mathematics (`products.rs`), 2D B-spline + 5-harmonic diurnal background evaluation (`background.rs`), and Format A/B FITS product writers match `laxpcl1.f` (4057 lines) and `backshiftv3.f` (1243 lines) to IEEE-754 double precision (`~1.7e-13 – 4.0e-13` verified against `gfortran` on `BSPEV2` + harmonics).
- **Bugs Found & Fixed During Audit (`@laxpc-dev`, #51 -> #61, Verified)**:
  1. **[FIXED, HIGH]** `pipeline.rs`: `background_pha_path` did not substitute `lxp10`/`lxp1` -> `lxp20`/`lxp30` per unit, causing `Lxp20` (256 ch) to skip background spectrum/LC generation and `Lxp30` (512 ch) to silently use LAXPC10's background `.pha`; also `.lcbk` is now written from `BackgroundModel2D` even when `background_pha_path` is `None`.
  2. **[FIXED, MEDIUM]** `pipeline.rs` + `products.rs`: Format A `level2.event.fits` previously mixed `0..511` (`Lxp10`/`Lxp30`) and `0..255` (`Lxp20`) channel scales when `spectrum_channels == None`; fixed via uniform `fa_channels = config.spectrum_channels.unwrap_or(1024)` (`TLMIN2 = 0`, `TLMAX2 = fa_channels - 1`).
  3. **[FIXED, MEDIUM]** `background.rs`: LAXPC20 post-day-550 / post-day-2977 linear drift overrides were applied unconditionally inside single-epoch `GainTable::peak_30kev_at`; moved paired `(source_day > 550 && background_day > 550)` / `> 2977` logic into `compute_gain_shift` matching `backshiftv3.f:298-307`.
  4. **[FIXED, MEDIUM]** `pipeline.rs` + `gti.rs`: 2D background evaluation used raw `interp_linear` on `EARTHLON` without `360 -> 0 deg` wrap protection; replaced with wrap-safe `interp_lon_deg` + `< 0 -> +360` normalization.
  5. **[FIXED, LOW/MEDIUM]** `gti.rs`: `generate_standard_gti` previously fell back to `user_gti.clone()` when standard auto-GTI was empty; fixed to strict `gti_utc.intersect(user_gti)`.
- **CLI Findings Found & Fixed (`@cli-dev`, #52 -> #64, Verified)**:
  1. **[FIXED, MEDIUM]** `cli.rs` + `laxpc_cmd.rs:547-650`: `astrosat laxpc lightcurve` now validates `--caldb` when `--energy-range` is passed, applies `--pha-min`/`--pha-max`, uses wrap-safe `interp_lon_deg` for `--bkg-model`, and sets `meta.exposure` to net exposure (`gross - frame_gap_loss`).
  2. **[FIXED, LOW]** `laxpc_cmd.rs:862-873`: `astrosat laxpc backshift` now computes and emits `lcbk_shift_fraction` (`fshbk`) over `--pha-min..=--pha-max`.
- **Post-Wave-1 Oracle & E2E Fixes Verified During Neon Self-Review (`#79 -> #80`, `#81 -> #83/#84`)**:
  1. **[FIXED, HIGH]** `background.rs:303-547` + `pipeline.rs:753-763` + `laxpc_cmd.rs:348,878-888` (`#83/#84`): Added `BackgroundCatalog` (`backshiftv3.inp` parser + 103-entry `BUILTIN_BACKSHIFT_CATALOG` fallback) and `--bkg-peak-ch` CLI override so `compute_gain_shift` uses the measured blank-sky `chb` instead of interpolating `gaina8.dat` across `lxhv.dat` HV jumps (verified to `2.22e-16 – 8.88e-16` cps against compiled `backshiftv3.e` in `audit/laxpc_head2head_verification.md`).
  2. **[FIXED, LOW/MEDIUM]** `cli.rs:891-907` + `laxpc_cmd.rs:720-808,429-440` (`#80`): `astrosat laxpc spectrum` now supports `--caldb`, `--energy-range`, `--pha-min`, `--pha-max` (masking out-of-band channels to `0` and recomputing `total_corrected_rate_cps`); `astrosat laxpc decode -o <evn>` falls back to `RA_PNT`/`DEC_PNT` when `RA_NOM`/`RA_OBJ` are absent.
- **Remaining Intentional `PARTIAL_GAP` Items**: `nul = -2` (L1 <20 keV scaling with `ncl = 100` after `2.59736e8 s`), `nul = 0` (4-quartile `.pha` blending), `nul = -3` (deprecated 3D `BSPEVN`), post-day-917 `lxp10cshm{141,175,197,212}v2.2.channel` auto-loading, `laxpcl1.f:FFT` quicklook `.four`, and `backshiftv3.f:diur` `_corr.lcsr` / `_corr.lc` output.

---

## 2. Complete Feature & Subroutine Traceability Matrix

### 2.1 `laxpcl1.f` (4057 lines) -> Rust Mapping

| Legacy Subroutine / Block | Legacy Location | Rust Module / Symbol | Rust Location | Classification |
| :--- | :--- | :--- | :--- | :--- |
| `program laxpcl1` (CLI/inputs/orchestration) | `laxpcl1.f:64-428` | `PipelineConfig`, `run_pipeline`, `run_unit_pipeline`, `discover_l1_orbits` | `pipeline.rs:55-264,1066-1194` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine frameno` (2242-byte frame sort, dedup, stray, rollover, gaps) | `laxpcl1.f:2666-2845` | `FrameHeader::parse`, `sort_and_dedup_frames`, `merge_multi_orbit_frames` | `frames.rs:66-429` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine conver1` (40-bit STBG + Dec 2024 5th-byte fix + `0xEF`/`0xEE` + 5-byte EA) | `laxpcl1.f:2128-2340` | `RawEaWord::{unpack,resolve}`, `EaDecoderState`, `decode_ea_frame` | `ea.rs:48-404` | `EXACT_MATCH` |
| Multi-orbit EA stream & overlap/gap tracking (`tim <= t00`, `>= 50 ms`, `>= 10 s`) | `laxpcl1.f:513-653` | `decode_ea_orbits`, `passes_event_filter` | `ea.rs:427-534` | `INTENTIONAL_IMPROVEMENT` |
| TCT FITS read & glitch rejection `dt_sps/dt_ins in (1.000, 1.001]` | `laxpcl1.f:656-710` | `load_laxpc_tct_files`, `astrosat_core::tct::LaxpcTct` | `gti.rs:92-117`, `tct.rs:325-470` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine converk` & MKF merge | `laxpcl1.f:713-803,2499-2568` | `load_merged_mkf_files`, `astrosat_core::mkf::MkfTable` | `gti.rs:120-163`, `mkf.rs:14-194` | `EXACT_MATCH` |
| `subroutine converb` & `32.768 s` BBC binning | `laxpcl1.f:805-861,2573-2631` | `decode_bb_frame`, `BbcSeries::{from_ordered_orbits,from_records}` | `bb.rs:31-224` | `EXACT_MATCH` |
| `subroutine converb1` (unused single-record helper) | `laxpcl1.f:2636-2660` | Not needed (`decode_bb_frame` handles 31x64-byte frames) | — | `EXACT_MATCH` |
| `subroutine convert` (unused 2-byte/EA diagnostic) | `laxpcl1.f:2442-2496` | Covered by `RawEaWord::resolve` & unit tests | `ea.rs:82-98` | `EXACT_MATCH` |
| GTI frame-gap loss `tgap` & spectrum deadtime `fcd = 1/(1 - tau*R_all)` | `laxpcl1.f:878-957` | `compute_gti_frame_gap_loss`, `deadtime_factor`, `build_spectrum` | `products.rs:37-73,170-300` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine converc` & `.spec` / ASCII `.pha` output | `laxpcl1.f:959-1075,2413-2439` | `write_spectrum_fits`, `write_ascii_pha` | `products.rs:675-731,831-855` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine orbit` (`AstroSat.orb` 1-s interpolation) | `laxpcl1.f:2850-3016` | `write_orbit_from_mkf`, `write_orbit_fits` | `pipeline.rs:997-1063`, `products.rs:809-828` | `INTENTIONAL_IMPROVEMENT` |
| Light curve & standard GTI screening (`ELV`, `SAA`, `cthr`, `uthr`, `merge <0.3s`, `drop <5s`) | `laxpcl1.f:1187-1617` | `bin_events_raw`, `lxp10_uld_anomaly_threshold`, `generate_standard_gti` | `products.rs:339-423`, `gti.rs:198-411` | `INTENTIONAL_IMPROVEMENT` |
| Sub-`0.2 s` 1-s running mean & GTI-trimmed `.lcurv` / `.lcbk` + `subroutine converl` | `laxpcl1.f:1273-1286,1623-1737,2080-2103,2343-2364` | `smooth_all_counts_for_deadtime`, `build_lightcurve`, `write_lightcurve_fits` | `products.rs:113-143,426-504,735-789` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine convere` & `.evn` FITS + `subroutine converg` & `.gtis` FITS | `laxpcl1.f:1740-1934,2367-2410` | `build_calibrated_events`, `write_format_b_event_fits`, `write_gti_fits` | `products.rs:523-599,792-805` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine FFT` & `.four` power spectrum | `laxpcl1.f:1936-1970,3054-3124` | Not emitted in `astrosat-laxpc` (`astrosat timing` provides $Z_n^2$/epoch folding) | `crates/astrosat-cli/src/timing.rs` | `PARTIAL_GAP` |
| Background `.pha`/`.spec` modes (`nul = 0, -1, -2, -3`) | `laxpcl1.f:1136-1183,1973-2078` | `BackgroundModel2D`, `BackgroundPha`, `run_optional_background_and_gain` (`nul=-1` equivalent) | `background.rs:912-1107`, `pipeline.rs:708-925` | `PARTIAL_GAP` |
| `function swp8` (8-byte BE->LE swap) | `laxpcl1.f:3020-3035` | `f64::from_be_bytes` in `astrosat-fits` / `frames.rs` | `frames.rs:76-77` | `EXACT_MATCH` |
| `subroutine DIVDIF`, `function NEARST` | `laxpcl1.f:3148-3252` | `interp_linear` (`NUSE=2`), `interp_lon_deg`, `LaxpcTct::{interpolate,interpolate_neville}` | `bb.rs:296-329`, `gti.rs:414-440`, `tct.rs:425-490` | `EXACT_MATCH` |
| `subroutine BSPLIN`, `BSPEVL`, `BSPEV2` | `laxpcl1.f:3436-3625,3951-3985` | `eval_bspline_basis`, `eval_bspline_2d` | `background.rs:1113-1197` | `EXACT_MATCH` |
| `function BSPEVN` (3D B-spline for `nul = -3`) | `laxpcl1.f:4013-4057` | Not implemented (deprecated by TIFR `readme.laxpcl1:165-168`) | — | `PARTIAL_GAP` |
| `subroutine BSPFIT`, `SVD`, `SVDEVL` | `laxpcl1.f:3253-3408,3643-3918` | Offline fitter routines (unused at runtime in `laxpcl1.f:1107-1110`) | — | `EXACT_MATCH` |

### 2.2 `backshiftv3.f` (1243 lines) & Format A Products -> Rust Mapping

| Legacy Subroutine / Product | Legacy Location | Rust Module / Symbol | Rust Location | Classification |
| :--- | :--- | :--- | :--- | :--- |
| `program backshift`: `backshiftv3.inp` + `gaina8.dat` + `lxhv.dat` gain & HV-jump tracking | `backshiftv3.f:108-313` | `GainTable`, `HvJumpTable`, `BackgroundCatalog`, `compute_gain_shift` | `background.rs:119-547,729-785` | `EXACT_MATCH` |
| Secular background polynomial scaling (`co3`, `maxday`) | `backshiftv3.f:79-95,325-333` | `secular_background_scale` | `background.rs:574-607` | `EXACT_MATCH` |
| Post-day-917 LAXPC10 `cshm{141,175,197,212}v2.2.channel` auto-switch | `backshiftv3.f:334-347` | `resolve_channel_map` always loads `lxp{10,20,30}.channel` | `pipeline.rs:955-964`, `laxpc_cmd.rs:904-906` | `PARTIAL_GAP` |
| `subroutine gainsh` (non-linear `e2` channel rebinning + ULD fill) | `backshiftv3.f:661-726` | `shift_background_spectrum` | `background.rs:789-870` | `EXACT_MATCH` |
| `cdel` / `fshbk` & `.lcbk` gain-shift scaling | `backshiftv3.f:103-104,382-394,599-612` | `background_lightcurve_shift_fraction`, `apply_lightcurve_gain_shift` | `background.rs:874-908`, `pipeline.rs:819-896` | `INTENTIONAL_IMPROVEMENT` |
| RMF selection (`gar`, LAXPC30 `rho` density, `idel`/`id2`/`id4`/`id5`, LAXPC10 `cshm` >917) | `backshiftv3.f:417-525` | `estimate_lxp30_density`, `recommend_rmf` | `background.rs:610-723` | `INTENTIONAL_IMPROVEMENT` |
| `subroutine converc` (`.spec` FITS writer) | `backshiftv3.f:629-655` | `write_spectrum_fits` | `products.rs:675-731` | `EXACT_MATCH` |
| `subroutine diur` (`_corr.lcsr`, `_corr.lc`, MKF bad-interval warning) | `backshiftv3.f:736-823` | Not emitted (`lcurv` + `lcbk` emitted separately) | — | `PARTIAL_GAP` |
| `subroutine DIVDIF`, `NEARST`, `SVD`, `SVDEVL` | `backshiftv3.f:848-1243` | `interp_linear` (`SVD`/`SVDEVL` unused at runtime after Apr 2021) | `bb.rs:296-329` | `EXACT_MATCH` |
| Format A `laxpc_make_event` (`level2.event.fits`: `TIME`, `Channel`, `Layer`, `LAXPC_No`, `GTI`) | Format A spec (`bary.rs:10`) | `write_format_a_event_fits`, `write_format_a_event_fits_with_channels` | `products.rs:604-671`, `pipeline.rs:590-605,1146-1186` | `INTENTIONAL_IMPROVEMENT` |
| Format A `laxpc_make_stdgti`, `laxpc_make_spectra`, `laxpc_make_backspectra`, `laxpc_make_lightcurve`, `laxpc_make_backlightcurve` | Format A workflow | Unified via `run_pipeline` + CLI subcommands (`gtigen`, `spectrum`, `backshift`, `lightcurve`) | `pipeline.rs`, `laxpc_cmd.rs` | `EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT` (see §5) |

---

## 3. Bit-Level & Mathematical Verification

### 3.1 2242-Byte Frame Header, Sync, Ordering, Deduplication & Gap Accounting (`frames.rs` vs `laxpcl1.f:frameno`)

- **Frame geometry & fields (`EXACT_MATCH`)**:
  - `L1_FRAME_SIZE = 2242` (`types.rs:311`).
  - `l1_time`: bytes `0..8` BE `f64` (`frames.rs:76` vs `laxpcl1.f:2156`).
  - `raw_seq`: bytes `36..44` BE `i64` (`frames.rs:77` vs `laxpcl1.f:2157,2771`).
  - `sync = raw[44] == 0xDE` (`-34` as `i8`, `frames.rs:12,78,115` vs `laxpcl1.f:2730`).
  - `unit_id = raw[45]` in `1..=3` (`frames.rs:79,116` vs `laxpcl1.f:2730` `abs(ibuf(kj+46)-2) <= 1`).
  - `mode_byte = raw[46]`: `0xEA` (`-22`) for EA, `0xBB` (`-69`) for BB (`frames.rs:14-16,80`).
  - `stbg_56`: bytes `47..54` BE (`frames.rs:82-84` vs `laxpcl1.f:2773-2775`); `stbg_40`: bytes `49..54` BE (`frames.rs:86-88` vs `laxpcl1.f:2172`).
  - `status_flag = raw[57] == 0` (`frames.rs:90,117` vs `laxpcl1.f:2731`).
  - `frame_counter = u16::from_be_bytes([raw[58], raw[59]])` (`frames.rs:91` vs `laxpcl1.f:2720-2721`).
  - `err_flags[0..3]` (`ie1, ie2, ie3`): bytes `2092..2098` BE `i16`, requiring `ie1 <= 1` (`frames.rs:93-95,118` vs `laxpcl1.f:2724-2731`).
- **Stray rejection, 16-bit rollover, duplicate & gap logic**:
  - Stray jump rollback (`frames.rs:244-267` vs `laxpcl1.f:2737-2755`): rejects previous frame when `|ifram - prev| in (6000, 60000)` and `|prev - prev2| in (6000, 60000)`. `EXACT_MATCH`.
  - 16-bit rollover (`frames.rs:273-289` vs `laxpcl1.f:2759-2767`): increments `rollover` by 1 when `ifram < prev_nfo - 12000` and `(cur_idx - last_rollover_idx) > 5000`, and adjusts late pre-rollover frames `> 55000` by `-65536`. `EXACT_MATCH`.
  - Long-gap wrap & duplicates (`frames.rs:299-330` vs `laxpcl1.f:2790-2817`): handles `> 2000 s` wrap cleanly via `unwrapped += 65_536` (fixing Fortran's `j1 = j1 + 65535` off-by-1), rejects `stbg_56 > 5e9` and `dt_new > 2 * dt_old` using actual `existing.header.l1_time` (fixing Fortran's `timen(j1)` uninitialized-array quirk). `INTENTIONAL_IMPROVEMENT`.
  - Note (`frames.rs:408-425`): `merge_multi_orbit_frames` updates `max_time` without reading it; in practice `pipeline.rs` passes per-orbit batches to `decode_ea_orbits`, matching `laxpcl1.f`'s `tim <= t00` cross-orbit overlap filter.

### 3.2 5-Byte EA Event Decoding (`ea.rs` vs `laxpcl1.f:conver1`)

- **40-bit STBG clock & Dec 2024 5th-byte rollover fix (`EXACT_MATCH`)**:
  - `STBG_CLOCK_HZ = 100_000.0` (`10 us` resolution), `STBG_EPOCH_OFFSET_SEC = 1.8121823219164e8`, `POW256_4 = 4_294_967_296.0` (`ea.rs:16-18,216-226` vs `laxpcl1.f:2191-2199`).
  - `dt = l1_time - timf / 1e5 - 1.8121823219164e8`, `d4 = dt * 1e5 / 256^4`, `n4 = trunc(d4 ± 0.3)`, `i5 += n4`, `timf += n4 * 256^4`. Exact match.
- **`0xEF` coarse time markers, glitch filter & retroactive alignment (`EXACT_MATCH`)**:
  - Carry propagation across `m2, m3, m4, it5` (`ea.rs:282-297` vs `laxpcl1.f:2235-2247`).
  - Single-marker glitch filter against `j2_hist[nt-2]`, `j3_hist[nt-2]`, `j4_hist[nt-2]` (`ea.rs:300-326` vs `laxpcl1.f:2248-2259`).
  - Retroactive alignment of pre-marker events when `ichk == 0` (`ea.rs:332-347` vs `laxpcl1.f:2261-2282`).
- **5-byte record bit-unpacking (`ea.rs:51-77` vs `laxpcl1.f:2290-2300`, `EXACT_MATCH`)**:
  - `fine_ticks = bytes[1].wrapping_sub(2)` (`(i4 - 2) mod 256`).
  - `ia1 = bytes[0] & 0x0F`, `ia2 = (bytes[0] >> 4) & 0x0F`.
  - `k1 = bytes[4] & 0x01`, `l1 = (bytes[4] >> 2) & 0x03`, `k2 = (bytes[4] >> 4) & 0x01`, `l2 = (bytes[4] >> 6) & 0x03`.
  - `ic1 = (bytes[2] << 2) | l1`, `ic2 = (bytes[3] << 2) | l2`.
- **Single vs Xe K-escape double events (`ea.rs:82-98` vs `laxpcl1.f:2303-2314`, `EXACT_MATCH`)**:
  - Single: `ic2 == 0` -> `(ic1, ia1, EventKind::Single)`.
  - Double: `(k1 == 1 || k2 == 1) && ia1 < 8 && ia2 < 8` -> `pha = min(ic1 + ic2, 1023)`, `anode = ia2` iff `k1 == 1 && (ic1 <= ic2 || (k2 == 0 && ic2 > 1))` else `ia1`.
  - Frame-end truncation on `na0 > 4 || hit_eof_marker` (`0xEE`, `ea.rs:373-383` vs `laxpcl1.f:2330-2338`) and stream validity `time_ins > 1.0`, `1 <= pha <= 999`, `1 <= anode <= 10` (`ea.rs:450,472` vs `laxpcl1.f:534,549-550`): exact match.
  - Note (`ea.rs:531-533`): `passes_event_filter` compares `ev.pha` directly in `0..=1023` rather than `laxpcl1.f:363-364`'s `icl+1..=icu+1` shift (identical at default `0..=1023`, cleaner 0-based convention for custom sub-ranges; `INTENTIONAL_IMPROVEMENT`).

### 3.3 BBC 31x64-Byte Subframe Decoding (`bb.rs` vs `laxpcl1.f:converb`)

- **Subframe decoding (`EXACT_MATCH`)**:
  - 31 subframes `j = 1..=31`, `base = j * 64`; 56-bit STBG clock at `base+47..base+54` (`bb.rs:48` vs `laxpcl1.f:2591`).
  - 16-bit signed big-endian counters `i1` (`LLD` at `base+56`), `i2` (`JAN` at `base+58`), `i3` (`LDU` at `base+60`), `i4` (`IGE` at `base+64`), summed with Euclidean moduli `4096, 4096, 2048, 2048` (`bb.rs:53-61` vs `laxpcl1.f:2594-2612`).
  - Emits 16 records on odd `j`, applying `> 2.06 s` jump repair (`prev + (prev - prev2)` for `j > 3`, `prev + 1.28e-3` for `j == 3`) and `total_rate_counts = jan + ldu + ige` (`bb.rs:63-84` vs `laxpcl1.f:2613-2622`).
- **32.768 s binning & ULD quartiles**:
  - `BBC_BIN_SEC = 0.128 * 256.0 = 32.768 s`; bin index `floor((time_ins - t0_ins + 1e-3) / 32.768)` (`bb.rs:165` vs `laxpcl1.f:853`), `t0_ins` selected from first valid record in `(1e5, 1e9]` (`INTENTIONAL_IMPROVEMENT` over Fortran's `kb == 1` assumption at `laxpcl1.f:841`).
  - Quartile bounds `(750, 1050)`, `(780, 1100)`, `(740, 1090)` for Lxp10/20/30 (`types.rs:87-93` vs `laxpcl1.f:98`), `uld_quartile_index` clamped to `0..=5` (`bb.rs:232-237`). `EXACT_MATCH` on quartiles `1..=4`.

### 3.4 TCT Glitch Rejection & Interpolation (`tct.rs`, `gti.rs` vs `laxpcl1.f`)

- **Glitch filter (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**:
  - `LAXPC_TCT_RATIO_MIN = 1.000`, `LAXPC_TCT_RATIO_MAX = 1.001`: rejects rows with `dt_sps / dt_ins <= 1.000 || > 1.001` (`tct.rs:336-405` vs `laxpcl1.f:692-703`). Rust also scans for the first valid anchor pair cleanly rather than only patching row 1 (`INTENTIONAL_IMPROVEMENT`).
- **Interpolation accuracy (`EXACT_MATCH`)**:
  - `laxpcl1.f` uses `DIVDIF` with `NUSE = 2` (2-point piecewise linear) for all instrument-to-UTC conversions (`lines 1233, 1604, 1608, 1642, 1668, 1672, 1862, 1900, 1904`); `LaxpcTct::interpolate` uses $O(\log N)$ piecewise linear interpolation (`tct.rs:422-487`), and `interpolate_neville` is additionally provided for higher-order polynomial interpolation.
  - `utc_to_ins` (`gti.rs:168-194`) inverts `LaxpcTct` via 12 Newton/secant steps to `< 1e-10 s`.

### 3.5 Standard GTI Screening (`gti.rs` vs `laxpcl1.f:1210-1614`)

- **Elevation & SAA thresholds (`EXACT_MATCH`)**:
  - `elv * je > 3.5` deg (`je = -1.0` for Earth occultation, `gti.rs:262,337` vs `laxpcl1.f:1337`).
  - `lon < 249.0 && lon > saaex` (`saaex` default `3.0`, extended to `40.0` when `tu0 = xb mod per` falls in `(1e4, 3.5e4)` for `nday > 1104` or `(1e4, 3e4)` for `186 < nday <= 1100`, `gti.rs:307-317,338` vs `laxpcl1.f:257-258,1321-1327`).
  - Wrap-safe `interp_lon_deg` (`gti.rs:414-440`) prevents `359 -> 1 deg` interpolation glitches (`INTENTIONAL_IMPROVEMENT` over `laxpcl1.f:1318`).
- **Active detector & LAXPC10 ULD anomaly (`EXACT_MATCH`)**:
  - `count_rate_threshold`: `80.0` cps (Lxp10/20), `50.0` cps (Lxp30 `t <= 2.37e8`), `20.0` cps (Lxp30 `t > 2.37e8`) (`types.rs:98-109` vs `laxpcl1.f:1210,1224-1225`).
  - `lxp10_uld_anomaly_threshold` (`gti.rs:198-231` vs `laxpcl1.f:1237-1246`): all 9 UTC epoch thresholds (`2.600928e8..4.351461e8 s`) match verbatim.
- **Interval boundaries, merge/drop & exposure (`INTENTIONAL_IMPROVEMENT`)**:
  - Exact bin boundaries `[tc - 0.5*tbin, tc + 0.5*tbin]` avoid `laxpcl1.f:1311-1314`'s integer-second truncation (`il1, iu1`) and entry-bin skip (`tu1+1`).
  - `merge_gaps(0.3).drop_short(5.0)` (`gti.rs:366-368` vs `laxpcl1.f:1509-1564`), exact `frame_gap_loss_sec` subtraction (`products.rs:49-73` vs `laxpcl1.f:881-896`), and strict `gti_utc.intersect(user_gti)` (`gti.rs:380-389`, fixed in #61).

### 3.6 Deadtime, Running Mean & Spectral/LC Products (`products.rs` vs `laxpcl1.f`)

- **Non-paralyzable deadtime (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**:
  - Formula $R_{\text{corr}} = R_{\text{obs}} / (1 - \tau R_{\text{all}})$ with $\tau = 43\,\mu\text{s}$ (`DEFAULT_DEADTIME_SEC = 43.0e-6`) and singularity guard $R_{\text{all}} < 23{,}000\text{ cps}$ (`DEADTIME_RATE_GUARD_CPS = 23_000.0`, `products.rs:25-45` vs `laxpcl1.f:354,917,1303,1660`). Spectrum `deadtime_factor` also enforces the `< 23000 cps` and `denom > 0` guards.
  - Single-only rejected main-anode doubles are included in $R_{\text{all}}$ (`products.rs:204-208,233-235` vs `laxpcl1.f:623-625,911-914`).
  - Running 1-s mean for `tbin < 0.2 s` (`smooth_all_counts_for_deadtime`, `products.rs:113-143`): exact $O(N)$ prefix-sum equivalence to `laxpcl1.f:1226,1274-1286,1648-1659` (`ntb = floor(0.5/tbin + 0.2)`).
  - Error propagation: `selected_stat_err[ch] = sqrt(rate[ch] / exposure_sec)` (`products.rs:285` vs `laxpcl1.f:954`) and `error_cps[i] = sqrt(corr_rate[i] / tbin)` (`products.rs:473` vs `laxpcl1.f:2094`).
  - Default channels: `256` for `Lxp20`, `512` for `Lxp10`/`Lxp30` (`types.rs:78-83` vs `laxpcl1.f:357,360`).

### 3.7 `backshiftv3.f` Gain, 2D B-Spline Background & RMF Logic (`background.rs` vs `backshiftv3.f`, `laxpcl1.f`)

- **2D B-spline + 5-harmonic diurnal model (`EXACT_MATCH`)**:
  - `BackgroundModel2D::from_fortran_bytes` (`background.rs:1002-1063`) parses Fortran unformatted sequential records (`no, na, iw, mbk, cav, per`, knots `xo`, `xa`, `bcof`) with strict signed/overflow validation.
  - `eval_bspline_basis` + `eval_bspline_2d` + 5 harmonics (`per = 84494.9 s` from `backfit.*`, `lon > 252.0` cut, `[cav - 40, cav + 40]` clamp, `background.rs:1067-1197`) verified directly against `gfortran` executing `laxpcl1.f:BSPLIN`/`BSPEV2` on `backfit.lxp10mar16`: max absolute discrepancy across test points is `3.98e-13` cps (machine precision).
- **`backshiftv3.inp` / `gaina8.dat` / `lxhv.dat` & secular polynomial (`EXACT_MATCH`)**:
  - `BackgroundCatalog` (`background.rs:303-547` vs `backshiftv3.f:108-169,295`): parses `backshiftv3.inp` or falls back to `BUILTIN_BACKSHIFT_CATALOG` (all **103/103** epoch entries and all secular polynomial coefficients verified `0` diffs against `data/ref_src/laxpc/backshiftv3.inp`). Resolves measured blank-sky `chb` by filename stem tag or `|bkg_day - entry.day| < 2.0`, with `--bkg-peak-ch` manual override (`#83/#84`).
  - `GainTable::peak_30kev_at` (`background.rs:191-242` vs `backshiftv3.f:248-294`): linear `gaina8.dat` interpolation, one-sided interpolation within `±1.0` day of `lxhv.dat` HV jumps (with OOB guards at `it >= days.len()` and `it == 0` fixing Fortran's unguarded `t(ndd+1)` / `t(0)` accesses), LAXPC30 `> 892` off error and `> 883` exponential extrapolation `450 * exp((ndays - 871)/30)`.
  - `compute_gain_shift` (`background.rs:729-785` vs `backshiftv3.f:298-333`): paired LAXPC20 linear drift fits (`230.1 - 0.00571*d` when both `> 550`, `228.2 - 0.00571*(d - 2977)` when both `> 2977`, fixed in #61), XSPEC offset `0.08919 - 0.0003001*ndays`, and `secular_background_scale` (`DEFAULT_SECULAR_POLY_LXP{10,20,30}`, `maxday = 2321` for Lxp20 / `875` for Lxp10/30, `background.rs:574-607`).
- **Spectral & LC gain shifting (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**:
  - `shift_background_spectrum` (`background.rs:789-870` vs `backshiftv3.f:661-726`): `e2i = [-1.98e-3, -2.67e-3, -1.97e-3]`, `iuld = [631, 556, 629]` (with `icu = 158` for Lxp10 `ndays > 917`), fractional channel rebinning (with `.min(nch)` and `last_n1 >= 2` / `.max(2)` OOB guards), `de < 0` gap fill, unmodified `icu+1..nch` tail, followed by `secular_scale` and `A1..A7` re-sum. Verified against compiled `backshiftv3.e` to `2.22e-16 – 8.88e-16` cps across positive and negative shifts.
  - `background_lightcurve_shift_fraction` & `apply_lightcurve_gain_shift` (`background.rs:874-908` vs `backshiftv3.f:103-104,382-394,606`): uses exact 0-based `ch_min..=ch_max`, fixing Fortran's 1-based vs 0-based off-by-1 shift (`backshiftv3.f:103-104` sets `lcl=lcl/ncf, lcu=lcu/ncf` and tests 1-based `i.ge.lcl.and.i.le.lcu`, shifting sub-bands down by 1 channel and dropping the final channel `255`/`511` at default `0..1023`).
- **RMF selection (`INTENTIONAL_IMPROVEMENT`)**:
  - `estimate_lxp30_density` (`background.rs:610-619` vs `backshiftv3.f:428-430`): `1 - (ndays - 130)*1.5e-3` (`> 130`), `0.4 * exp(-(ndays - 500)/200)` (`> 700`).
  - `recommend_rmf` (`background.rs:623-723` vs `backshiftv3.f:417-524`): matches `irho` 5% quantization (`< 0.08 -> 6`, `< 0.05 -> 4`), exact `irho < 84` > `dch > 32.5` > `id2` precedence, `id2`/`id4`/`id5` rounding (fixing Fortran line 451's `id4=95` typo to `id5=95`), LAXPC10 post-917 `gar` & offset steps (`140, 175, 197, 212, 220` and `cshm{idel:03}{L1}{SE}v2.2.rmf`, adding `SE` tag matching actual CALDB files).

### 3.8 Format B & Format A FITS Schemas (`products.rs` vs `laxpcl1.f` & Format A)

| Product | Columns (`Rust` vs `Legacy`) | Key Headers & Differences | Status |
| :--- | :--- | :--- | :--- |
| Format B `.evn` (`LAXPC_MODE_EA` + `GTI`) | `TIME(1D)`, `PHA(1I)`, `ANODE_ID(1B)`, `DETECTOR_ID(1B)` (legacy had spaces `'ANODE ID'`, `'DETECTOR  ID'`) | `MJDREF=55197.0`, `MJDREFI=55197`, `MJDREFF=0.0`, `TIMEDEL=1e-5`, `TIMEPIXR=0`, `TLMIN2=0`, `TLMAX2=nch-1` (fixes legacy `TLMIN1=1, TLMAX1=nca` & duplicate 2nd `GTI` HDU) | `INTENTIONAL_IMPROVEMENT` |
| Format B `.spec` (`SPECTRUM`) | `CHANNEL(1I)`, `RATE(1D)`, `STAT_ERR(1D)` | OGIP `1.2.0`, `CHANTYPE='PHA'`, `TLMIN1=0`, `TLMAX1=nch-1`, `POISSERR=F`, `DEADTIME`, `RESPFILE`, `BACKFILE` | `EXACT_MATCH` |
| Format B `.lcurv` / `.lcbk` (`RATE` + `GTI`) | `TIME(1D)`, `RATE(1D)`, `ERROR(1D)`, `FRACEXP(1D)` (`FRACEXP` added over 3-col `laxpcl1.f`, matching Format A `.lc`) | `TIMEPIXR=0.5` (accurate bin center vs legacy `0`), `MJD-OBS`, `TIMEZERO=0.0`, `CHANRANG`, `DEADAPP=T`, `VIGNAPP=T` | `INTENTIONAL_IMPROVEMENT` |
| Format B `.gtis` (`GTI`) | `START(1D)`, `STOP(1D)` | `HDUCLASS='OGIP'`, `HDUCLAS1='GTI'`, `MJDREF=55197.0` | `EXACT_MATCH` |
| `AstroSat.orb` (`ORBIT`) | `TIME(1D)`, `X..VZ(1D)` (upgrades legacy `J` int `TIME` + `E` `f32` pos/vel) | 1-s cadence `[t0.floor()-3, t1.ceil()+5]`, `TIMEDEL=1.0`, `MJDREF=55197.0` | `INTENTIONAL_IMPROVEMENT` |
| Format A `level2.event.fits` (`LAXPC_EVENTS` + `GTI`) | `TIME(1D)`, `Channel(1I)`, `Layer(1B)`, `LAXPC_No(1B)` | Uniform `Channel` grid `0..1023` (or explicit `spectrum_channels`), `TLMIN2`/`TLMAX2`, `MJDREFI=55197`, `RA_NOM`/`DEC_NOM`, `TIMESYS='UTC'` on HDU 0 & 1 | `INTENTIONAL_IMPROVEMENT` |

---

## 4. Verification of Fixes & Test Evidence

- **Differential Numerical Verification (Fortran vs Rust formula)**:
  - `gfortran` execution of `laxpcl1.f:BSPLIN` + `BSPEV2` + 5 diurnal harmonics on `backfit.lxp10mar16` (`40x10` knots, `cav=262.631383615356`, `per=84494.9 s`):
    - `(100.0, 0.0, 189910000.0)`: Fortran `250.605414405066` vs `250.605414405066` (`diff = 3.98e-13`)
    - `(100.0, 0.0, 189931124.0)`: Fortran `241.776540684699` vs `241.776540684699` (`diff = 2.84e-13`)
    - `(240.0, -4.5, 190000000.0)`: Fortran `274.097753135343` vs `274.097753135343` (`diff = 1.71e-13`)
- **Post-Fix `astrosat-laxpc` Test & Lint Suite** (run after `@laxpc-dev` #61):
  - `cargo clippy -p astrosat-laxpc --all-targets -- -D warnings`: **0 warnings, 0 errors**.
  - `cargo test -p astrosat-laxpc`: **10/10 passed, 0 failed**:
    - `tests/unit_tests.rs` (**8/8 passed**): includes `stage1`..`stage6` and `audit_regressions_safety_math_and_edge_cases` (verifying `interp_lon_deg` wrap guard, empty auto-GTI intersection, LAXPC20 paired vs unpaired gain drift, and corrupt-input bounds).
    - `tests/crab_integration.rs` (**2/2 passed**): `crab_9000000252_parity_against_official_l2_reference` (1024-ch parity against official L2 `.pha`/`.gti`) and `crab_9000000252_full_multi_unit_auto_gti_format_b_and_format_a` (verifying Format B + Format A `TLMIN2=0, TLMAX2=1023`, all 3 units' `.spec` & `.lcbk` generation, NASA `fitsverify 4.22` 0 errors, and native `barycenter_fits` on both Format B `.evn` and Format A `level2.event.fits`).

- **Post-Fix `astrosat-cli` Test & Lint Suite** (run after `@cli-dev` #64):
  - `cargo clippy -p astrosat-cli --all-targets -- -D warnings`: **0 warnings, 0 errors**.
  - `cargo test -p astrosat-cli`: **9/9 passed, 0 failed** (`3` unit tests + `6` integration tests including `cli_defaults_symmetrically_override_config_and_energy_range_requires_caldb`).

---

## 5. Verified CLI Audit Items (`astrosat-cli`)

All items raised in `#52` (`@cli-dev` `#64`), `#79` (`@cli-dev` `#80`), and `#81` (`@cli-dev` `#84`) were verified against current source and tests:
1. **[FIXED & VERIFIED, #64]** `crates/astrosat-cli/src/laxpc_cmd.rs:554-669`: `run_laxpc_lightcurve` validates `--caldb` + channel map when `--energy-range` is passed, applies `--pha-min`/`--pha-max`, uses shortest-arc wrap-safe `interp_lon_deg` for `--bkg-model`, and sets `meta.exposure = net_exposure_sec` (`gross - frame_gap_loss`).
2. **[FIXED & VERIFIED, #64/#84]** `crates/astrosat-cli/src/laxpc_cmd.rs:872-940`: `run_laxpc_backshift` resolves `bkg_peak_ch` via `--bkg-peak-ch` or `BackgroundCatalog::load_or_builtin` (`backshiftv3.inp` + 103-entry builtin), shifts the background spectrum, and computes/emits `lcbk_shift_fraction` (`fshbk`) over `--pha-min..=--pha-max`.
3. **[FIXED & VERIFIED, #80]** `crates/astrosat-cli/src/laxpc_cmd.rs:720-808,429-440`: `run_laxpc_spectrum` supports `--caldb`, `--energy-range`, `--pha-min`, and `--pha-max` (zeroing out-of-band channels and recomputing `total_corrected_rate_cps`); `run_laxpc_decode` falls back to `RA_PNT`/`DEC_PNT` for `RA_OBJ`/`DEC_OBJ`.
4. **[FIXED & VERIFIED, #84]** `crates/astrosat-cli/src/laxpc_cmd.rs:348` + `cli.rs:787,932`: `--bkg-peak-ch` wired symmetrically in both `astrosat laxpc run` and `astrosat laxpc backshift`.

---

## 6. Neon Self-Review Verification Addendum

A fresh first-principles re-audit of every LAXPC module (`crates/astrosat-laxpc/src/{frames,ea,bb,gti,products,background,pipeline,types}.rs`, `crates/astrosat-core/src/tct.rs`, `crates/astrosat-cli/src/{laxpc_cmd,cli}.rs`) against legacy `laxpcl1.f`, `backshiftv3.f`, `backshiftv3.inp`, `lxhv.dat`, and the head-to-head oracle report (`audit/laxpc_head2head_verification.md`) verified **zero remaining `BUG/MATH_DRIFT` defects**:

1. **Programmatic `backshiftv3.inp` Catalog & Secular Polynomial Parity**: Automated comparison of all **103** entries of `BUILTIN_BACKSHIFT_CATALOG` (`background.rs:325-429`) and all 3 secular polynomial arrays (`DEFAULT_SECULAR_POLY_LXP{10,20,30}`, `maxday=2321/875`) against `data/ref_src/laxpc/backshiftv3.inp` confirmed **0 discrepancies**.
2. **Combined `laxpcl1.f` + `backshiftv3.f` Background Scaling Proof (`pipeline.rs:784-896`)**: Verified step-by-step that `scaled_orig = pha.rates * (mean_bk / orig_total)` -> `shifted = shift_background_spectrum(..., secular_scale, ..., &scaled_orig)` -> `bk_lc *= (sbi / orig_total) * (1 + fshbk)` is algebraically identical to running `laxpcl1.f` (`nul = -1`, `fbk = bkav / sbk`, `frbk = sbi / sb0`) followed by `backshiftv3.f` (`fact = ss(ndays)/ss(ndayb)`, `rb8 *= (1 + fshbk)`).
3. **Bit-Level & Index Boundary Proofs**:
   - `frames.rs:82-88,339`: Verified little-endian `ib4(9-k) = ibuf(kj+46+k)` with `ib4(8)=0` (`laxpcl1.f:2773-2775`) maps bit-for-bit to `raw[47..54]` big-endian `stbg_56`, `raw[49..54]` to `stbg_40`, and `data[36..44]` to `ibuf(kj+45-k)=ib2(k)`.
   - `ea.rs:218-224,261,300-326,493`: Verified 407-word loop (`pos = 60, 65, ..., 2090` == Fortran `do i=61,2092,5`), 0-based `j2_hist[nt-2]` indexing == Fortran 1-based `j2(nt-2)`, explicit `n4 = 0` when `d4 == 0.0` (avoiding Fortran's uninitialized `save` `n4`), and persistent `prev_valid_tim` across orbits (fixing `laxpcl1.f:504` `jp = 0` dropping cross-orbit gaps).
   - `products.rs:113-143,251-265,888-904`: Verified exact 0-based algebraic equivalence of `smooth_all_counts_for_deadtime` (`itbl`/`itbu`) to `laxpcl1.f:1274-1286` for all `nt >= 1`, exact `512`/`256`-channel rebinning (`2*ch..2*ch+1` / `4*ch..4*ch+3`), strict `(start, stop)` event membership (`in_sorted_intervals`), and non-strict `[tc - tbin, tc] ⊆ [start, stop]` bin membership (`bin_inside_gti`).
   - `background.rs:206-234,709-717,827-850,874-899`: Verified HV-jump one-sided interpolation matches `backshiftv3.f:272-292` while preventing Fortran's `t(ndd+1)`/`t(0)` OOB accesses; verified `recommend_rmf` preserves Fortran's `irho < 84` > `dch > 32.5` > `id2` precedence while fixing `backshiftv3.f:451`'s `id4=95` typo; verified `background_lightcurve_shift_fraction` fixes Fortran `backshiftv3.f:103-104,382`'s 0-based vs 1-based off-by-1 channel range shift (explaining the `0.0023` `fshbk` difference in `ep1_lxp2_mar16_band`).
   - `lxhv.dat` note: Three unit-epochs have back-to-back HV jumps 1 day apart (`Lxp30`: day `25/26`; `Lxp10`: days `913/914` and `952/953`); `nearest_jump_within` (`.find`) selects the earlier jump while Fortran's overwrite loop selects the later jump strictly between them, which is immaterial since neither side of a 1-day inter-jump gap is jump-free.
   - `products.rs:25-45`: Verified `< 23,000 cps` singularity guard protects both spectrum and light curve (preventing `laxpcl1.f:917`'s unguarded `-2.27e6 cps` sign-flip explosion at `23,500 cps`).
