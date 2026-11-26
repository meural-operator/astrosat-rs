# AstroSat-rs Public-Release End-User Scenario, Stress & Cross-Instrument Verification Report

- **Date:** 2026-11-26 (Post-Neon Self-Review & Final Verification)
- **Verifier Role:** Public-Release End-User Scenario, Stress & Cross-Instrument Verifier (`@e2e-scenarios`)
- **Executable Tested:** Standalone release binary installed via [`scripts/install.sh`](file:///d:/ISRO/scripts/install.sh) (`--prefix=/mnt/d/ISRO/.tmp/public_scenarios`) to `/mnt/d/ISRO/.tmp/public_scenarios/bin/astrosat` (incorporating all Neon-review fixes `#107`, `#110`, `#113`, `#116`, `#117`)
- **Verification Script:** [`verification/verify_public_scenarios.py`](file:///d:/ISRO/verification/verify_public_scenarios.py)
- **JSON Evidence Output:** [`verification/reports/public_scenarios.json`](file:///d:/ISRO/verification/reports/public_scenarios.json)
- **Overall Verdict:** **PASS (`8/8` scenario suites passed; `149/149` FITS files passed NASA `fitsverify 4.22` with `0` errors and `0` warnings; `0` panics across `69` adversarial inputs; `109/109` generated products passed `astrosat inspect --full` & `astrosat validate --strict`; clean `env -i` deterministic parity verified)**

---

## 1. Standalone Binary & Zero-Dependency Audit

Installed cleanly via:
```bash
PREFIX=/mnt/d/ISRO/.tmp/public_scenarios bash /mnt/d/ISRO/scripts/install.sh --prefix=/mnt/d/ISRO/.tmp/public_scenarios
```
- **Build & Install:** `cargo build --release -p astrosat-cli` succeeded cleanly; `scripts/install.sh` smoke verification passed (`--version`, `czti addspec --help`, `laxpc backshift --help`, etc.).
- **Shared Library Audit (`ldd`):** Links only `linux-vdso.so.1`, `libgcc_s.so.1`, `libm.so.6`, `libc.so.6`, `/lib64/ld-linux-x86-64.so.2`. Zero external C/Fortran scientific libraries (`cfitsio`, `heasoft`, `glog`, `pil`, `gfortran`).
- **Storage Compliance:** All temporary outputs, synthetic FITS datasets, and installed test binaries strictly stored under `/mnt/d/ISRO/.tmp/public_scenarios/`; zero writes to `/home/diat`, `/tmp`, or `C:`.

---

## 2. Pass / Fail Summary Matrix

| ID | Scenario Suite | Status | Runtime | Key Observed Metrics |
|---|---|---|---|---|
| `wf1` | Chained Per-Stage CZTI CLI Pipeline (ObsID `9000000406`) | **PASS** | 100.50 s | Exact default event parity (`Q0..Q3`: `43281, 44276, 41337, 37140`), exact `TIME`/`ENERGY`/`DPH`/`BalCC` image (`0.10°/pix` WCS)/`.pha` (`COUNTS`, `STAT_ERR`, `QUALITY`, `GROUPING`)/`.lc` (`RATE`)/`.rsp` (bit-exact `MATRIX` table+heap & `EBOUNDS`)/`bary` match vs `czti run`; generalized chain within `0.20%–0.58%` (post-`#117` temporal `flickpixclean` + pairwise `superbunchclean`); `--ra 0.0 --dec 0.0` off-axis rejected & omitted pointing reports true MKF coords; `119/119` FITS `fitsverify` (`0` err, `0` warn) |
| `wf2` | Chained Per-Stage LAXPC CLI Pipeline (Crab `9000000252`) | **PASS** | 111.82 s | Exact `.evn` (`8,939,817` events), exposure (`3029.793 s`), `.lcurv` (`3379.692 cps`), `.spec` (`3379.406 cps`), & `bary` parity vs `laxpc run`; filtered `3:20 keV --tbin 0.01 --layer1-only --single-only` LC (`2694.096 cps`) & spectrum (`2694.023 cps`); `backshift` standard (`+1.800 ch`), LXP3 HV-jump catalog (`+29.8 ch`) + `--bkg-peak-ch 220` (`-7.0 ch`), & LXP10 day-950 `cshm141` (`fshbk = 0.550028`); Format A `bary` + `timing` (`--fdot -3.693e-10`, `29.658 Hz`); `18/18` FITS `fitsverify` (`0` err, `0` warn) |
| `wf3` | Declarative Config (`TOML` + `JSON`) & CLI Default Overrides | **PASS** | 39.84 s | Exact `TOML` vs `JSON` parity for CZTI, LAXPC, and Bary; all CLI default overrides verified; HEASARC FITS `JPLEPH.405` matches binary DE405 (`< 1e-9 s`), `JPLEPH.200` + `--frame fk5` -> `JPL-DE200`; explicit `bary --mkf` outranks `cfg.orbit_file` & sibling `AstroSat.orb` (`< 1e-6 s` vs `--orbit`) |
| `wf4` | Idempotency & Guardrails | **PASS** | 2.88 s | `czti evtencor` on `GAINCOR='YES'` exact no-op (`[43281, 44276, 41337, 37140]`); double `bary` on `TIMESYS='TDB'` rejected cleanly (`rc != 0`); `laxpc run/lightcurve/spectrum --energy-range` without `--caldb` rejected cleanly with actionable errors |
| `sc2a` | Planted Spectral Power-Law ($\Gamma=1.7$) + Gaussian Line Recovery | **PASS** | 0.52 s | Exposure `49.9995 s` (exact); deadtime-corrected rate `615.8963 cps` (`rel_err = 2.95e-11`); Gaussian centroid `99.82 ch` (true `100.0`), excess rate `119.05` vs `122.10 cps` (`1.79 sigma`); reduced $\chi^2 = 1.15$; slope `-1.669` (true `-1.70`); `4/4` FITS `fitsverify` (`0` err, `0` warn) |
| `sc2b` | Planted `17.3456 Hz` Pulsar Across `2017-01-01` Leap Second + Orbit | **PASS** | 2.36 s | Dual-branch pre/post leap-second grid (zero boundary clamping); Analytical, binary DE405, and HEASARC FITS `JPLEPH.405` `bary` + `timing` all recover `17.3456102 Hz` (`err = 1.02e-5 Hz < 1e-4 Hz`), $Z_2^2 = 27,091.7$, pre/post leap-second phase jump **`0.000215` cycles (`12.4 µs`)** (analytical `0.000216` cycles); binary DE405 vs FITS `JPLEPH.405` `< 1e-9 s`; un-barycentered raw `.evn` degraded (`Z_2^2 = 14,319.6`); `8/8` FITS `fitsverify` (`0` err, `0` warn) |
| `sec3` | Adversarial & Corrupt-Input Stress Testing | **PASS** | 6.38 s | `69/69` corrupt/truncated/overflow/empty-GTI/missing-file/invalid-CLI invocations (incl. `--n-trials 99999999`, `--freq NaN`, `--fdot inf`, `--oversampling 1000000`, unpaired `--ra`) returned clean non-zero exits (or valid 0-row outputs on empty FITS GTI) with **0 panics, 0 segfaults, 0 backtraces** |
| `sec4` | Cross-Instrument Validation & Clean `env -i` Determinism | **PASS** | 35.32 s | `109/109` generated CZTI & LAXPC products passed `inspect --full` & `validate --strict` (`0` errors, `0` warnings); clean `env -i PATH=...` ran `--help`, `czti run`, `laxpc run`, `bary`, `timing`, `inspect`, `validate` with **exact deterministic parity** vs normal-env runs |

---

## 3. Defects Found During Scenario Testing & Verified Fixed

### Wave 1 (Initial Scenario Verification — Board `#79` → Fixed `#80`)
Four real CLI/pipeline defects in `crates/astrosat-cli` were uncovered, reported to `@cli-dev` (`#79`), fixed in `#80`, rebuilt via [`scripts/install.sh`](file:///d:/ISRO/scripts/install.sh), and verified fixed:
1. **`astrosat laxpc spectrum` silently ignored `--energy-range`, `--pha-min`, `--pha-max` and lacked `--caldb`** (`crates/astrosat-cli/src/cli.rs`, `crates/astrosat-cli/src/laxpc_cmd.rs`): Now rejects `--energy-range` without `--caldb` (`rc=1`, actionable error) and masks out-of-band channels (`2694.02 cps` for `3:20 keV` L1 single-anode, boundary channels `0` and `511` = `0.0`).
2. **`astrosat laxpc decode -o <evn>` omitted `RA_OBJ`/`DEC_OBJ` when L1 `modeEA` only carried `RA_PNT`/`DEC_PNT`** (`crates/astrosat-cli/src/laxpc_cmd.rs`): Now falls back to `RA_PNT`/`DEC_PNT`, enabling seamless chained `decode -> bary`.
3. **`astrosat czti evtclean --badpix` rejected `SPECTROSCOPIC` (`flag == 1`) pixels** (`crates/astrosat-cli/src/czti_cmd.rs`): Now keeps `flag < pix_flag::NOISY` (`<= 1`), matching `czti run` bit-for-bit (`[43281, 44276, 41337, 37140]`).
4. **`astrosat czti pha2energy --ss` silently ignored Level-1 `modeSS` FITS files** (`crates/astrosat-cli/src/cli.rs`, `crates/astrosat-cli/src/czti_cmd.rs`): Added `--tct` so L1 `modeSS` packets decode cleanly into `TEMP`, yielding exact `PI`, `ENERGY`, and event counts matching `czti run`.

### Wave 2 (Neon Self-Review — Board `#99` → Fixed `#116`)
Two additional real CLI edge-case defects were caught during the Neon self-review, reported to `@cli-dev` (`#99`), fixed in `#116`, rebuilt via [`scripts/install.sh`](file:///d:/ISRO/scripts/install.sh), and verified fixed:
5. **`astrosat bary --mkf` was silently shadowed by sibling `AstroSat.orb` and `cfg.orbit_file`** (`crates/astrosat-cli/src/bary_cmd.rs`):
   - *Initial bug:* Passing explicit `--mkf <custom.mkf>` on an event file sitting next to `AstroSat.orb` (always produced by `astrosat laxpc run`) or alongside a `--config` setting `bary.orbit_file` silently ignored `--mkf`.
   - *Verified fix:* Precedence is now `args.orbit` → `args.mkf` → `cfg.orbit_file` → sibling `AstroSat.orb`; verified in `wf3` (`r_mkf["orbit"] == LAXPC_MKF` even when both `--config` `orbit_file` and sibling `AstroSat.orb` exist, matching `--orbit` to `< 1e-6 s`).
6. **`astrosat czti bindata` JSON/human report printed bogus `(RA=0.0000, DEC=0.0000)` when `--ra`/`--dec` were omitted** (`crates/astrosat-cli/src/czti_cmd.rs`):
   - *Initial bug:* Even though `run_bindata` resolved `RA_PNT`/`DEC_PNT` from MKF, `CztiBindataReport` hardcoded `args.ra.unwrap_or(0.0)` / `args.dec.unwrap_or(0.0)`.
   - *Verified fix:* Now resolves actual MKF pointing coordinates into `CztiBindataReport`, verified in `wf1` alongside explicit `--ra 0.0 --dec 0.0` off-axis rejection.

### Test-Suite Improvements Made During Neon Self-Review (`verification/verify_public_scenarios.py`)
- Eliminated dead/duplicate code (`czti_toml` overwrite in `wf3` and dead `einstein_grid` line in `sc2b`).
- Replaced the coarse single-grid leap-second interpolation in `sc2b` with separate dense 2000-point pre-leap (`<= t_leap - 1 ms`) and post-leap (`>= t_leap + 1 ms`) branches excluding the 1-second leap-second gap, removing boundary clamping and tightening measured pre/post leap-second phase coherence from `0.00194` to **`0.000215` cycles (`12.4 µs`)** across all three ephemeris backends (analytical, binary DE405, HEASARC FITS `JPLEPH.405`).
- Expanded end-to-end assertions to cover: BalCC sky image parity + `0.10°/pix` WCS plate scale, `.pha` `STAT_ERR`/`QUALITY`/`GROUPING`, `.lc` `RATE`, `.rsp` bit-exact `MATRIX` table+heap & `EBOUNDS`, `--ra 0.0 --dec 0.0` off-axis rejection, LAXPC LXP3 HV-jump catalog + `--bkg-peak-ch` + LXP10 day-950 `lxp10cshm141v2.2.channel` (`fshbk = 0.550028`), space-separated negative `--fdot -3.693e-10`, HEASARC FITS `JPLEPH.405`/`JPLEPH.200`, `bary --mkf` precedence, 5 new adversarial input-bound cases (`--n-trials 99999999`, `--freq NaN`, `--fdot inf`, `--oversampling 1000000`, unpaired `--ra`), and explicit `env -i` determinism equality assertions.

---

## 4. Detailed Verification Results

### 4.1 Workflow 1 — Chained Per-Stage CZTI CLI Pipeline
- **Subcommands exercised:** `czti decode` → `czti bunchclean` (default & `--generalized`) → `czti pha2energy` (with `--ss` + `--tct`) → `czti gtigen` → `czti gaas` → `czti datasel` → `czti clean` (`--algo default` & `--algo generalized`) → `czti evtencor` → `czti evtclean` (with `--badpix`) → `czti dpigen` (`--emin 0.0 --emax 1000.0`) → `czti image` (`--algo balcc` & `--algo fft --oversampling 16`) → `czti bindata` (`--gen-evt`, omitted coords & explicit off-axis `--ra 0.0 --dec 0.0`) → `czti addspec` (both averaged-exposure & `--add-exposure`) → `bary`.
- **Parity vs `astrosat czti run`:**
  - Default event counts (`Q0..Q3`): `43281`, `44276`, `41337`, `37140` — **exact match**.
  - Event `TIME` (`atol=1e-9`), `ENERGY` (`atol=1e-5`), `DPH` arrays, `BalCC` `61×61` sky images (WCS `|CDELT1| = |CDELT2| = 0.10°/pix`), `.pha` (`COUNTS`, `STAT_ERR`, `QUALITY`, `GROUPING`), `.lc` (`RATE`), `.rsp` (bit-exact `MATRIX` table+heap and `EBOUNDS`), and `bary` `TIME` (`atol=1e-9`) — **exact match**.
  - Generalized clean event counts (`run` vs `chain`, post-`#117` temporal `flickpixclean` & pairwise `superbunchclean`): `Q0: 39365 vs 39230` (`0.34%`), `Q1: 39527 vs 39446` (`0.20%`), `Q2: 37461 vs 37251` (`0.56%`), `Q3: 34238 vs 34041` (`0.58%`) — expected minor ordering difference (`run_pipeline` runs `noisypixclean` before `datasel`, whereas chained `clean --algo generalized` runs after `datasel`).
  - `czti addspec`: co-added 4-quadrant exposure `364.481 s` (average mode) and verified `--add-exposure` sums exposure across inputs (`> 1.5x`).
  - NASA `fitsverify 4.22`: **`119/119`** FITS products passed with `0` errors and `0` warnings.

### 4.2 Workflow 2 — Chained Per-Stage LAXPC CLI Pipeline
- **Subcommands exercised:** `laxpc gtigen`, `laxpc decode`, `laxpc orbit`, `laxpc lightcurve` (default + filtered `--energy-range 3:20 --tbin 0.01 --layer1-only --single-only`), `laxpc spectrum` (default + filtered `--energy-range 3:20 --layer1-only --single-only`), `laxpc backshift` (standard `backlxp10jan16.pha`, HV-jump `backlxp30aug17.pha` + `--bkg-peak-ch 220.0`, and post-day-917 `backlxp10apr18.pha` with `lxp10cshm141v2.2.channel`), `bary` (Format B `.evn` and Format A `level2.event.fits`), and `timing` (`--period`, space-separated `--fdot -3.693e-10`, `--unit 1 --layer 1 --pi-min 10 --pi-max 300`).
- **Parity vs `astrosat laxpc run`:**
  - Event count (`lxp1`): `8,939,817` — **exact match**.
  - Net GTI exposure: `3029.7935 s` — **exact match**.
  - Default light curve mean rate: `3379.6921 cps` — **exact match**.
  - Default spectrum total rate & per-channel rates: `3379.4058 cps` (`atol=1e-9`) — **exact match**.
  - Barycentered Format B `TIME`: **exact match** (`atol=1e-9`).
  - Filtered `3:20 keV`, `0.01 s`, Layer-1, Single-only light curve: `2694.0959 cps`, `TIMEDEL = 0.01 s`.
  - Filtered `3:20 keV`, Layer-1, Single-only spectrum: `2694.0233 cps`, out-of-band channels strictly `0.0`.
  - `laxpc backshift`: LXP1 standard `delta_ch_30kev = 1.8002` (`lx10v1.1.rmf`); LXP3 day-700 `backlxp30aug17.pha` catalog `delta_ch ~ +29.8` & `--bkg-peak-ch 220.0` exact `-7.0 ch` shift; LXP10 day-950 `backlxp10apr18.pha` `cshm141` `lcbk_shift_fraction = 0.550028` (exact match to `backshiftv3.e`).
  - Format A `bary` + `timing` (`--fdot -3.693e-10 --unit 1 --layer 1 --pi-min 10 --pi-max 300`): recovered Crab frequency `29.658 Hz` (`Z_2^2 > 100`).
  - NASA `fitsverify 4.22`: **`18/18`** FITS products passed with `0` errors and `0` warnings.

### 4.3 Workflow 3 — Declarative Config (`TOML` / `JSON`), Ephemerides & CLI Overrides
- Verified identical CZTI and LAXPC outputs for equivalent `TOML` and `JSON` `--config` files across all config keys (`CztiFileConfig`, `LaxpcFileConfig`, `BaryFileConfig`).
- Verified explicit CLI flags override **every** config value, including when the CLI supplies a default value (`--entry-stage 1`, `--exit-stage 3`, `--output-type both`, `--gti-type quad`, `--algo default`, `--image-algo balcc`, `--no-gen-evt`, `--no-format-a`, `--no-orbit`, `--no-ascii`, `--no-single-only`, `--no-layer1-only`, `--no-user-gti-only`, `--tbin 1.0`, `--deadtime 0.000043`, `--elv-min 3.5`, `--frame icrs`, `--ephem analytical`, `--checksum`).
- Verified HEASARC FITS `JPLEPH.405` matches binary JPL `lnxp1600p2200.405` to `< 1e-9 s` on event timestamps, HEASARC FITS `JPLEPH.200` + `--frame fk5` emits `PLEPHEM='JPL-DE200'` / `FK5`, and explicit `bary --mkf` overrides both `cfg.orbit_file` and sibling `AstroSat.orb` (`< 1e-6 s` vs `--orbit`).
- Verified `--json` output across all top-level and instrument subcommands parses cleanly as valid JSON with expected scientific fields.

### 4.4 Workflow 4 — Idempotency & Guardrails
- `astrosat czti evtencor` re-run on `GAINCOR='YES'` `.evt`: clean exit `0`, identical event counts (`[43281, 44276, 41337, 37140]`), identical `PI` and `ENERGY` (no double gain shift).
- `astrosat bary` re-run on `TIMESYS='TDB'` `.evt`: clean non-zero exit (`rc != 0`), no output file written, clear `TDB` / already-barycentered error message.
- `astrosat laxpc run --energy-range 3-20` / `3:20`, `laxpc lightcurve --energy-range 3:20`, and `laxpc spectrum --energy-range 3:20` without `--caldb`: all fail cleanly (`rc != 0`) with actionable error messages.

### 4.5 Ground-Truth Planted Physics Scenarios (`sc2a`, `sc2b`)
- **Scenario 2A (Synthetic LAXPC L1 → L2 Spectral Recovery):**
  - Constructed synthetic 2242-byte AstroSat LAXPC `modeEA` L1 FITS + TCT + MKF with `30,000` in-GTI events over `T_ins = 49.9995 s` (`R_obs = 600.006 cps`, deadtime $\tau = 43\ \mu\text{s}$, analytical deadtime-corrected rate $R_{\text{corr}} = 615.8963\text{ cps}$), distributed as 80% power-law ($\propto \text{PHA}^{-1.7}$) + 20% Gaussian line ($\mu = 200$ PHA / `100` in 512-ch, $\sigma = 8$ PHA).
  - Recovered `EXPOSURE = 49.9995 s`, total rate `615.8963 cps` (`rel_err = 2.95e-11`), exact `STAT_ERR = sqrt(RATE / EXPOSURE)`, Gaussian centroid `99.82 ch`, line flux `119.05` vs `122.10 cps` (`1.79 sigma`), power-law slope `-1.669`, and full-band reduced $\chi^2 = 1.15$ (`4/4` FITS `fitsverify` `0` errors, `0` warnings).
- **Scenario 2B (Synthetic `17.3456 Hz` Pulsar Across `2017-01-01` Leap Second + LEO + Earth Orbit):**
  - Forward-modeled `~8,000` TDB pulsar photons at $f_0 = 17.3456\text{ Hz}$ spanning `2016-12-31T23:55:00` to `2017-01-01T00:05:00` UTC (`MJD 57754.0` leap-second boundary, `TAI-UTC`: `36 s → 37 s`), inverted on separate dense pre-leap / post-leap UTC branches through JPL DE405 Earth Roemer + `650 km` inclined LEO AstroSat orbit + Einstein (`TDB - TT`) into synthetic L1 `modeEA` frames.
  - Built-in analytical, binary JPL DE405, and HEASARC FITS `JPLEPH.405` `astrosat bary` → `astrosat timing` (`--df 2e-5`) all recovered $f = 17.3456102\text{ Hz}$ (`|Δf| = 1.02e-5 Hz`, 10x tighter than `1e-4 Hz`), $Z_2^2 = 27,091.7$, and pre- vs post-leap-second pulse phase coherence $\Delta\phi = 0.000215$ cycles (`12.4 µs`, analytical `0.000216` cycles), with binary DE405 vs FITS `JPLEPH.405` agreeing to `< 1e-9 s` (vs degraded raw un-barycentered $Z_2^2 = 14,319.6$; `8/8` FITS `fitsverify` `0` errors, `0` warnings).

### 4.6 Adversarial Stress Testing & Cross-Instrument Clean-Env Determinism (`sec3`, `sec4`)
- **Adversarial & Corrupt Inputs (`69/69` passed, `0` panics):** Tested empty `0`-byte files, `135`-byte non-FITS garbage, header-only truncated FITS, mid-data truncated FITS, `BITPIX=9`, `NAXIS1=-5`, `PCOUNT=2147483647`, `NAXIS2=100000000`, `TFIELDS=9999999`, empty FITS/ASCII GTIs, inverted `START > STOP` GTIs, missing MKF/TCT/CALDB/ORBIT paths, and invalid CLI flags/parameters (`--entry-stage 9`, `--tbin -1.0`, `--channels 333`, `--harmonics 0`, `--nbins 1`, `--n-trials 99999999`, `--freq NaN`, `--fdot inf`, `--oversampling 1000000`, unpaired `--ra`). Zero Rust panics, zero segfaults, zero backtraces.
- **Cross-Instrument Consistency (`109/109` products):** `astrosat inspect --full` verified `100%` valid `CHECKSUM`/`DATASUM`, monotonic `TIME` across all sorted event/rate/orbit HDUs, and zero non-finite timestamps; `astrosat validate --strict` reported `0` errors and `0` warnings across all `109` generated CZTI & LAXPC products.
- **Clean Environment (`env -i PATH=/mnt/d/ISRO/.tmp/public_scenarios/bin:/usr/bin:/bin`):** Verified `--help`, `czti run`, `laxpc run`, `bary`, `timing`, `inspect --full`, and `validate --strict` execute cleanly with zero external environment variable dependencies (`HEADAS`, `CALDB`, `PFILES`, `HOME`, `TMPDIR` all unset) and **exact deterministic parity** (`CZTI = [43281, 44276, 41337, 37140]`, `LAXPC = 8,939,817`, `bary = 8,939,817`, `Crab = 29.65805147 Hz`).
