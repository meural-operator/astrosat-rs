# Barycentric Timing, AstroSat Core & FITS 4.0 Mathematical Audit

**Scope**: Exhaustive line-by-line mathematical and specification audit comparing legacy **`as1bary`** (`data/ref_src/as1bary/as1bary/as1bary/{axBary.c,bary.c,scorbit.c,clock.c,dpleph.c,ctatv.c,bary.h,tai-utc.dat,README}`) and the **FITS 4.0 / OGIP 1.2–1.3** standards against **`crates/astrosat-core/src/{time.rs,tct.rs,gti.rs,mkf.rs,orbit.rs,bary.rs,caldb.rs}`** and **`crates/astrosat-fits/src/{card.rs,checksum.rs,io.rs,bintable.rs,image.rs,ogip.rs}`**.

**Post-Oracle & Neon Fresh-Eyes Review Status**:
- All 3 initial `BUG/MATH_DRIFT` findings (`#54` -> fixed `#56`), all post-oracle enhancements (`#82` -> `#85`: HEASARC FITS `JPLEPH.405`/`JPLEPH.200`, native `DE200`/`FK5` frame handling, Primary `TSTART`/`TSTOP` preservation, 16-sig-fig float formatting), and the Neon fresh-eyes review findings (`#98`/`#102` -> fixed `#107`: 0-row `EVENTS`/`GTI` `TSTART`/`TSTOP`/`DATE-OBS`/`DATE-END` barycentering, unsorted GTI `min`/`max` fold + GTI-only fallback, `TZERO = 2^63` exact two's-complement wrapping, and `127-term` docstring correction) have been independently verified:
  - `cargo test -p astrosat-fits -p astrosat-core`: **14/14 passed** (7 `astrosat-core` + 7 `astrosat-fits`, 0 failed, 0 ignored).
  - `cargo clippy -p astrosat-fits -p astrosat-core --all-targets -- -D warnings`: **0 errors, 0 warnings**.
  - Independent numerical verification (`time.rs::tdb_minus_tt` vs `ctatv.c` across 1900–2100): max residual **`1.301e-18 s` (`0.001 fs`)**.
  - Independent HEASARC FITS ephemeris constant check (`data/ephem/JPLEPH.405` & `JPLEPH.200`): `CLIGHT = 299792.458 km/s`, `EMRAT` dynamically read from `DE1` (`81.30056` / `81.300587`), `GM_SUN_OVER_C3` matches `GMS*(AU/CLIGHT)^3/86400^2` to `< 2e-15 s` (`< 10^-5 ns` Shapiro difference even at solar limb).
  - Cross-checked against live compiled `as1bary` C oracle (`#94`, `audit/heasoft_as1bary_numerical_verification.md`): DE405/DE200 binary and HEASARC FITS parity confirmed at **`0–1 ULP` (`<= 29.80 ns`)** on real CZTI `9000000406`, LAXPC Crab `9000000252` (8.94M events, `+1.26 ns` phase offset), 2016→2017 leap-second crossing, solar conjunction, `FK5 → ICRS`, and `jd == jd_end`.

---

## Summary Classification Matrix

| Area | Item | Classification | Legacy / Standard Reference | Rust Implementation | Notes |
|---|---|---|---|---|---|
| **Time & Leap Seconds** | AstroSat epoch & constants (`MJDREFI=55197`, `MJDREFF=0.0`, `TAItoTT=32.184`) | `EXACT_MATCH` | `axBary.c:275-280`, `bary.h:65,71` | `time.rs:7-17` | Identical mission constants |
| **Time & Leap Seconds** | Leap-second table (`1972-01-01` .. `2017-01-01+`, 28 entries) | `EXACT_MATCH` | `tai-utc.dat:14-41` | `time.rs:21-50` | Exact MJD & `TAI-UTC` values across 2012-07-01, 2015-07-01, 2017-01-01 |
| **Time & Leap Seconds** | Pre-1972 rubber-second rows (1961–1968) | `PARTIAL_GAP` | `tai-utc.dat:1-13`, `bary.c:151-155` | `time.rs:53-63` | Clamps pre-1972 to `10.0 s`; AstroSat launched 2015 |
| **Time & Leap Seconds** | Per-event `(leapevt - leapref)` across leap boundaries | `EXACT_MATCH` | `axBary.c:496-500,794-821` | `time.rs:92-126`, `bary.rs:198-230` | Algebraically identical `t_utc + leapevt + 32.184 + total_delay` |
| **Time & Leap Seconds** | Fairhead & Bretagnon (1990) `tdb_minus_tt` (`ctatv`) | `EXACT_MATCH` *(fixed `#56`, doc fixed `#107`)* | `ctatv.c:59-204`, `bary.c:179-203` | `time.rs:128-299` | All 127 exact `ctatv.c` terms (`93 T0 + 27 T1 + 6 T2 + 1 T3`, `1.3e-18 s` match; avoids `bary.c` daily linearization `-31.6 ns` sagitta) |
| **Barycentric Physics** | Roemer delay `(r_earth . n)/c + (r_sc . n)/c` | `EXACT_MATCH` | `bary.c:801,823` | `bary.rs:208-209` | Exact SI meter / `C_LIGHT` equivalence |
| **Barycentric Physics** | Earth-velocity Lorentz term `(r_sc . v_earth) / c^2` | `EXACT_MATCH` | `bary.c:721,798,824`, `dpleph.c:89,557` | `bary.rs:210` | `bary.c` `vce` is `ls/s = v/c`, so `DOT(scposn, vce)/c` = `(r_sc . v_earth)/c^2` |
| **Barycentric Physics** | Solar Shapiro delay `+(2 GM_sun / c^3) ln(1 + n . r_hat_{Sun->sc})` | `EXACT_MATCH` | `bary.c:802,809-810,825`, `dpleph.c:551-552` | `bary.rs:32,213-226` | Identical `Sun -> S/C` vector & sign + safe `-1+1e-12` clamp |
| **Barycentric Physics** | Standish FK5 -> ICRS rotation `c200to405` & native `DE200`/`FK5` | `EXACT_MATCH` *(updated `#85`)* | `bary.c:1236-1254`, `axBary.c:315,339-354` | `bary.rs:118-129,176-186` | Exact `(0.002, 0.012, 0.006)/206265` cross-product + renormalization; skips rotation when `DE200` + `Fk5` |
| **Orbit Interpolation** | `km -> m`, 3-pt centered acceleration, 2nd-order Taylor | `EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT` | `scorbit.c:158-235` | `orbit.rs:93-246` | Exact interior formula; improves boundary handling (`<= 120 s`), non-uniform sorting, `TUNIT` detection; avoids `scorbit.c:185` `+0.5` row offset |
| **Ephemeris** | Built-in analytical Earth/Moon/Sun ephemeris | `INTENTIONAL_IMPROVEMENT` *(fixed `#56`, docs updated `#85`)* | *(none in `as1bary`)* | `bary.rs:7-9,58,563-994` | Fixed from `35.7 ms` RMS / `73.3 ms` max to Mercury + 34-term Hill series (`~4–9 ms` vs DE405, honestly documented; use `JplDe` for `< 30 ns`) |
| **Ephemeris** | JPL DE200/405/421/430 Chebyshev reader (Binary + HEASARC FITS) | `EXACT_MATCH` *(completed `#85`)* | `dpleph.c:97-624` | `bary.rs:1004-1439` | Supports BOTH raw JPL binary (`lnxp*.405`) AND HEASARC FITS (`JPLEPH.405`/`JPLEPH.200`); exact `pc`/`vc`, `velfac`, `EMRAT`, `jd_end` clamp |
| **FITS Barycenter** | `barycenter_fits` header/table/GTI updates | `EXACT_MATCH` *(fixed `#56`, `#85`, `#107`)* | `axBary.c:207-887` | `bary.rs:294-563` | Preserves Primary `TSTART`/`TSTOP`, syncs `DATE-OBS`/`DATE-END`, handles 0-row `EVENTS`/`GTI`, guards `TDB`/`TCB`, incorporates `TIMEZERO`; `PARTIAL_GAP`: `MJDREF != 55197` / `TIMEUNIT='d'` |
| **TCT & MKF & GTI & CALDB** | `CztiTct`, `LaxpcTct`, `GtiSet`, `MkfTable`, `CaldbIndex` | `EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT` | AstroSat / OGIP CIF | `tct.rs`, `gti.rs`, `mkf.rs`, `caldb.rs` | Clean algebraic & edge-case behavior |
| **FITS 4.0 & OGIP** | 2880-byte blocking, 16-sig-fig floats, `CONTINUE`+`LONGSTRN`, `CHECKSUM`/`DATASUM`, `TZERO`/`BZERO` (incl. `2^63`), `P`/`Q` heaps, OGIP builders | `EXACT_MATCH` *(updated `#85`, `#107`; minor `PARTIAL_GAP`s)* | FITS 4.0 Std, Seaman/Pence/Rots (2012), OGIP 1.2/1.3 | `astrosat-fits/src/{card,checksum,io,bintable,image,ogip}.rs` | Verified `fitsverify 4.22` (0 warnings, 0 errors); minor gaps: `TSCAL`/`TZERO` on non-`u16` array cols, `Q64` writer, Random Groups |

---

## 1. AstroSat Time System & Leap Seconds (`time.rs` vs `axBary.c`, `bary.c`, `ctatv.c`, `tai-utc.dat`)

### 1.1 Mission Epoch, Physical Constants & Leap-Second Table (`EXACT_MATCH` / `PARTIAL_GAP`)
- **Constants**: `time.rs:7-17` defines `ASTROSAT_MJDREFI = 55197`, `ASTROSAT_MJDREFF = 0.0` (2010-01-01T00:00:00 UTC), `TT_MINUS_TAI = 32.184`, `SECS_PER_DAY = 86_400.0`, `JD_MINUS_MJD = 2_400_000.5`, `JD_J2000 = 2_451_545.0`. These match `axBary.c:275-280` (`mission = ASTROSAT; mjdrefi = 55197; mjdreff = 0.0; timesys = "UTC"`), `bary.h:65,71` (`SECDAY = 86400`, `TAItoTT = 32.184`), and `ctatv.c:64` (`2451545`).
- **Leap-Second Table**: All 28 integer leap-second epochs from `1972-01-01` (`MJD 41317.0`, `10.0 s`) through `2017-01-01` (`MJD 57754.0`, `37.0 s`) in `time.rs:21-50` match `tai-utc.dat:14-41` (`JD - 2400000.5`) entry-for-entry, including all AstroSat-relevant epochs:
  - `2009-01-01`: `MJD 54832.0`, `34.0 s` (active at `MJDREFI = 55197`) — `time.rs:46` vs `tai-utc.dat:38`.
  - `2012-07-01`: `MJD 56109.0`, `35.0 s` — `time.rs:47` vs `tai-utc.dat:39`.
  - `2015-07-01`: `MJD 57204.0`, `36.0 s` (active at AstroSat launch `2015-09-28`) — `time.rs:48` vs `tai-utc.dat:40`.
  - `2017-01-01`: `MJD 57754.0`, `37.0 s` — `time.rs:49` vs `tai-utc.dat:41`.
- **Pre-1972 Epochs (`PARTIAL_GAP`)**: `tai-utc.dat:1-13` contains 13 pre-1972 fractional/rubber-second entries (`1961-01-01` to `1968-02-01`, handled in `bary.c:151-155` via `leapcoeff[i]`). `time.rs:53-63` omits pre-1972 drift terms and returns `10.0 s` for any `mjd_utc < 41317.0`. This has zero impact on AstroSat data (`MJD >= 55197`).

### 1.2 Per-Event `(leapevt - leapref)` Handling Across Leap-Second Boundaries (`EXACT_MATCH`)
- **Legacy `as1bary` (`axBary.c:496-500, 564-609, 794-821`, `bary.c:628-639`)**:
  - At `baryinit`: `leapref = UTCtoTAI(&mtref) = 34.0` at `MJD 55197`, and `mjdRef` is shifted to `TT` (`55197 + (34.0 + 32.184)/86400.0`).
  - For each event timestamp `t_utc` (`axBary.c:795-821`):
    1. Computes `mtevt = 55197 + floor(t_utc / 86400)` and `leapevt = UTCtoTAI(&mtevt)` (`bary.c:148-152`, comparing `mjdUTC->MJDint >= leapsMJD[i]`).
    2. Sets `moreleap = leapevt - leapref`.
    3. Passes `t_utc + moreleap` to `xbarycorr(...)`, so `met2mjd` (`bary.c:77-97`) evaluates the Solar System ephemeris and `TTtoTDB` at `MJD_TT = 55197 + (t_utc + leapevt + 32.184)/86400.0`.
    4. Adds `xbarycorr(...) + (leapevt + TAItoTT)` to `t_utc`, yielding `t_TDB = t_utc + leapevt + 32.184 + (t_TDB - t_TT)` referenced to `MJDREF = 55197.0 TDB`.
- **Rust (`time.rs:53-126`, `bary.rs:198-230`)**:
  - `tai_minus_utc_astrosat(t_utc_sec)` computes `tai_minus_utc_mjd(55197.0 + t_utc_sec / 86400.0)`, stepping cleanly at `mjd_utc >= mjd_start`.
  - `astrosat_utc_to_tt_sec(t_utc_sec)` computes `t_tt_s = t_utc_sec + tai_utc + 32.184`.
  - `astrosat_utc_to_jd_tt(t_utc_sec)` splits `t_tt_s` into `(jd_hi, jd_lo)` preserving sub-nanosecond precision.
  - `bary.rs:230` adds `total_delay_s = t_TDB - t_TT` to `t_tt_s`.
  - **Equivalence**: Both produce the exact same value $t_{\text{TDB}} = t_{\text{UTC}} + (\text{TAI}-\text{UTC})(t_{\text{UTC}}) + 32.184 + \Delta_{\text{bary}}$, and `tt_sec_to_astrosat_utc` (`time.rs:98-106`) inverts `UTC <-> TT` via 3-step fixed-point iteration cleanly across all leap-second boundaries.
- **Note on `astrosat_utc_to_iso8601` (`time.rs:312-339`)**: Hardened by `@fits-core` (`#56`) against non-finite / out-of-range floats via finite checks, `[-6.0e12, 6.0e12]` bounds, `Duration::try_seconds`, and `checked_add_signed`. Also avoids `axBary.c:621`'s legacy `+66.184 s` double-counting bug in `DATE-OBS`/`DATE-END`.

### 1.3 Fairhead & Bretagnon (1990) `tdb_minus_tt` (`EXACT_MATCH` after `#56`/`#107` fixes)
- **Initial `BUG/MATH_DRIFT` (reported in `#54`)**: `time.rs` originally included only 36 of `ctatv.c`'s 127 terms (32 `T0` + 4 `T1`, omitting 61 `T0`, 23 `T1`, 6 `T2`, and 1 `T3` terms), had frequency mismatches in `TERMS_T0[0,1,2,7]`, and had erroneous amplitude/phase in `TERMS_T1[2]` (`1.454_44e-7` vs `0.269668e-6`) and `TERMS_T1[3]` (`1.008_70e-7, 5.814_773_255` vs `0.265919e-6, 5.836047367`). Numerical comparison against `ctatv.c` showed max discrepancy **`2,823.36 ns`** (`1,375.28 ns` RMS) over 2010–2030 and **`9,765.72 ns`** (`3,741.49 ns` RMS) over 1900–2100, violating the `< 100 ns` specification.
- **Fixes Verified (`time.rs:128-299`)**: Replaced in `#56` with all **127 exact terms** matching `ctatv.c:67-203` (`93 TERMS_T0` + `27 TERMS_T1` + `6 TERMS_T2` + `1` $t^3$ term), and docstring corrected to `"full 127-term"` in `#107`. Re-verified across 2010–2030 and 1900–2100 against `ctatv.c`: max discrepancy **`1.301e-18 s` (`0.001 fs`)**.
- **Note vs `bary.c:179-203`**: `bary.c` caches `ctatv` at integer day boundaries (`MJDint`) and linearly interpolates using 1-day central difference (`ctatv(l, 0.5) - ctatv(l, -0.5)`), which introduces up to `~31.6 ns` quadratic sagitta error (`#94`), whereas `time.rs:134-299` evaluates the full series directly at `(jd_tt_hi, jd_tt_lo)` (`INTENTIONAL_IMPROVEMENT`).

---

## 2. Barycentric Physics, Orbit & Ephemerides (`bary.rs`, `orbit.rs` vs `bary.c`, `scorbit.c`, `dpleph.c`)

### 2.1 Relativistic & Geometric Delay Terms (`EXACT_MATCH`)
Comparing `bary.rs:198-241` against `bary.c:742-843`:
1. **Earth Roemer Delay (`EXACT_MATCH`)**:
   - `bary.c:801,823`: `DOT(sourcedir, rce)` where `rce` is SSB->Earth in light-seconds (`dpleph.c:89`).
   - `bary.rs:208`: `dot3(ss.earth_pos_m, n) / C_LIGHT` where `earth_pos_m` is SSB->Earth in meters and `C_LIGHT = 299_792_458.0` m/s.
2. **Spacecraft Roemer Delay (`EXACT_MATCH`)**:
   - `bary.c:801,823`: `DOT(sourcedir, scposn/c)` where `scposn` is geocentric S/C in meters (`bary.c:721`) and `c = clight * 1000.0` in m/s (`dpleph.c:557`).
   - `bary.rs:209`: `dot3(sc.pos_m, n) / C_LIGHT`.
3. **Earth-Velocity Lorentz Term (`EXACT_MATCH`)**:
   - `bary.c:824`: `+ DOT(scposn, vce)/c`. In `dpleph.c:89`, `vce` is Earth velocity in **light-seconds/s** ($\mathbf{v}_\oplus / c$, dimensionless), while `scposn` is in **meters** and `c` is in **m/s**, yielding $\frac{\mathbf{r}_{\text{sc}} \cdot (\mathbf{v}_\oplus / c)}{c} = \frac{\mathbf{r}_{\text{sc}} \cdot \mathbf{v}_\oplus}{c^2}$ in seconds.
   - `bary.rs:210`: `dot3(sc.pos_m, ss.earth_vel_ms) / (C_LIGHT * C_LIGHT)`, where `pos_m` is in `m`, `earth_vel_ms` is in `m/s`, and `C_LIGHT^2` is in `m^2/s^2`, yielding seconds identically (`<= ~2.4 us` for AstroSat LEO).
4. **Solar Shapiro Delay (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**:
   - `bary.c:801-802,809-810,825`: Computes `rsa = rca - rcs` (Sun-to-Spacecraft vector: `r_sc_ssb - r_sun_ssb`), `cth = DOT(sourcedir, rsa)/|rsa|`, and adds `+ 2*msol*log(1+cth)` where `msol = GMS * (AU/CLIGHT)^3 / 86400^2` (`dpleph.c:550-552` = `4.9254909483e-6 s` in DE405, `4.9254909452e-6 s` in DE200).
   - `bary.rs:32,213-226`: Computes `r_sun_sc = ss.earth_pos_m + sc.pos_m - ss.sun_pos_m` (Sun-to-Spacecraft), `cos_theta = (dot3(n, r_sun_sc) / r_mag).clamp(-1.0 + 1e-12, 1.0)`, and adds `+ 2.0 * GM_SUN_OVER_C3 * (1.0 + cos_theta).ln()` with `GM_SUN_OVER_C3 = 4.925_490_947e-6`. Sign, vector orientation (`Sun -> S/C`), and constant match `bary.c` identically (`< 10^-5 ns` difference between DE405/DE200 `msol` constants); the `-1.0 + 1e-12` clamp prevents `ln(0) = -inf` at exact solar conjunction (`INTENTIONAL_IMPROVEMENT`).
5. **Standish FK5 -> ICRS Frame Rotation & Native DE200 Handling (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**:
   - `bary.c:1236-1254` (`c200to405`): `dir_405 = dir_200 + eps x dir_200` with `eps = (0.002, 0.012, 0.006) / 206265.0`.
   - `bary.rs:118-129` (`fk5_to_icrs`): Uses the exact same `rx, ry, rz` and cross-product signs (`x + ry*z - rz*y`, `y + rz*x - rx*z`, `z + rx*y - ry*x`), plus re-normalizes `[x/norm, y/norm, z/norm]` to preserve exact unit length (`INTENTIONAL_IMPROVEMENT`, `~5e-15` order).
   - `bary.rs:176-186` (added in `#85`): When `ephem` is `DE200` (`numde == 200`, native FK5/J2000 dynamical frame) and `frame == CoordFrame::Fk5`, skips the `FK5 -> ICRS` rotation (matching `axBary.c:339-354`, verified `0.000000 ns` vs `as1bary` in `#94`). Note: pairing `DE200` with `CoordFrame::Icrs` does not apply inverse `-eps x dir` (~25 ns, an unsupported combination in `axBary.c:339-342` as well).

### 2.2 Orbit Interpolation (`orbit.rs` vs `scorbit.c`) (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)
- **Units (`EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**: `scorbit.c:234` multiplies `km` positions by `1000.0` (`km -> m`). `orbit.rs:148-171` inspects `TUNIT` (`"m"` vs `"km"`, `"m/s"` vs `"km/s"`) with orbital-radius fallback (`r0 > 1e5`), scaling both `pos` and `vel` to SI `m` and `m/s`.
- **Taylor Interpolation (`EXACT_MATCH`)**: For interior points (`0 < i < n-1`), `orbit.rs:220-245` computes centered acceleration `acc[k] = (vel_ms[i+1][k] - vel_ms[i-1][k]) / (time[i+1] - time[i-1])`, `pos_m[k] = pos_m[i][k] + vel_ms[i][k]*dt + 0.5*dt^2*acc[k]`, matching `scorbit.c:228-235` identically.
- **Boundary & Non-Uniform Handling (`INTENTIONAL_IMPROVEMENT`)**:
  - `scorbit.c:185,190` assumes uniform `DELTAT`/`TIMEDEL` (and adds `+ 0.5` before `round(...)`, biasing sample selection for `frac < 0.5` as confirmed in `#94`), fails silently on raw L2 `.orb` files missing `DELTAT`/`TSTART`, and aborts (`*oerror = -2`) near the first/last 2 rows (`startinx <= 0 || startinx > numrows-2`).
  - `orbit.rs:61-74,198-227` sorts/deduplicates timestamps, uses exact nearest-sample binary search (`partition_point`), uses one-sided `(0, 1)` / `(n-2, n-1)` differences at endpoints, rejects non-finite timestamps, and permits up to `tol = 120.0 s` endpoint extrapolation. (Note: for `dt <= 3 s` in AstroSat LEO, 2nd-order Taylor truncation error is `< 0.1 ns`; at the extreme `120 s` extrapolation limit, LEO jerk `~0.008 m/s^3` introduces `~7.7 us`.)

### 2.3 Built-in Analytical Ephemeris vs JPL DE Reader (`bary.rs` vs `dpleph.c`)
1. **Built-in Analytical Ephemeris (`bary.rs:563-994`, `INTENTIONAL_IMPROVEMENT` after `#56`/`#85` fixes)**:
   - **Initial `BUG/MATH_DRIFT` (reported in `#54`)**: The original 9-term longitude / 6-term radius perturbation in `emb_heliocentric_ecliptic_au` had an inverted Venus `m=1` longitude sign (`cos(l_v - l_e + pi/2) = -sin(l_v - l_e)` instead of `+sin(l_v - l_e)`), erroneous phases on Venus `m=2`, a 3.8x oversized Saturn `m=1` term, and radial `d_r_au` populated with arcsec-to-radian magnitudes (`23.59e-6`, `34.96e-6`, `7.90e-6` -> `117e-6 AU` peak), producing **`35.7 ms` RMS / `73.3 ms` max** Roemer error against a 10-year N-body integration (violating the original `~50 us` / `< 100 us` claim by ~700x) and omitting Mercury from `sun_ssb_ecl`.
   - **Fix Verified (`bary.rs:7-9,58,594-927`)**: Added Mercury (`Standish 1992 Table 2a`, `M_MER = 1/6_023_600.0`, `~32 us` solar reflex) to `sun_ssb_ecl`, replaced the perturbation block with a 34-term frequency-domain Hill perturbation series (`pert_terms`), verified signs/amplitudes against linearized Newton's 2nd law perturbations (Venus `m=1` `(+4.89", -5.52e-6)` vs `(+4.96", -5.63e-6)`, Venus `m=2` `(-5.59", +15.94e-6)` vs `(-5.64", +16.08e-6)`, Jupiter `m=1` `(+7.07", +15.95e-6)` vs `(+7.18", +16.20e-6)`, Jupiter `m=2` `(-2.71", -9.23e-6)` vs `(-2.72", -9.25e-6)`, Mars `m=2` `(+2.05", +4.75e-6)` vs `(+2.08", +4.80e-6)`), and updated docstrings (`bary.rs:9,58`) to accurately report **`~4–9 ms (< 10 ms)`** accuracy vs JPL DE405 (verified at `4.15 ms` on Crab `9000000252` and `8.28 ms` on CZTI `9000000406` in `#94`; use `EphemerisKind::JplDe` for `< 30 ns` absolute phase alignment).
   - **Physical Limitations (`PARTIAL_GAP`)**: Standish Table 2a planetary orbits for `sun_ssb_ecl`, truncated ELP2000 Moon (`~4 us` Earth reflex), and J2000 ecliptic-to-equatorial rotation without ~0.02–0.04" ICRS frame bias (`~50–100 us`) remain bounded within the documented `< 10 ms` envelope.
2. **JPL DE200/405/421/430 Chebyshev Reader (`bary.rs:1004-1439` vs `dpleph.c:97-624`) (`EXACT_MATCH`)**:
   - **Dual-Format Support (`bary.rs:1020-1285`, completed in `#85`)**: Automatically detects `b"SIMPLE  ="` to load **NASA HEASARC FITS** ephemerides (`JPLEPH.405`, `JPLEPH.200`: `DE1` `Cname`/`Cvalue` for `DENUM`/`AU`/`EMRAT`, `DE2` `Pointer`/`NumCoeff`/`NumSubIntv` for `ipt[0..13]`, `DE3` `TSTART`/`TSTOP`/`TIMEDEL`/`NAXIS1` big-endian `ChebCoeffs`, matching `dpleph.c:530-605` 1:1), and otherwise loads **raw JPL binary** ephemerides (`lnxp1600p2200.405`, `linux_p1550p2650.430`, etc.: header offsets `2652/2660/2668/2680/2688/2696/2840/2844`, auto-detecting big/little endianness from `numde` at `2840`, deriving `ncoeff` with `ndim=2` for nutations `idx=11` and `ndim=3` otherwise).
   - **Chebyshev Evaluation & Scaling (`bary.rs:1287-1439`) (`EXACT_MATCH`)**: Subinterval indexing `sub_idx = floor(frac * nsub).min(nsub - 1)`, normalized time `tc = 2*(sub_f - sub_idx) - 1`, Chebyshev recurrence `pc[i] = 2*tc*pc[i-1] - pc[i-2]` & derivative `vc[i] = 2*tc*vc[i-1] + 2*pc[i-1] - vc[i-2]`, velocity factor `vfac = 2*nsub/span_days`, Earth recovery `emb - moon/(1+emrat)`, and unit scaling `km -> m` / `km/d -> m/s` match `dpleph.c:254-421,604`. Reverse-order summation `(0..ncf).rev()` reduces floating-point cancellation (`INTENTIONAL_IMPROVEMENT`).
   - **Boundary Clamping (`bary.rs:1293-1300`) (`EXACT_MATCH`)**: Clamps exact `jd == jd_end` to `max_data_rec` (both binary `data_record_offset = 2` and FITS `data_record_offset = 0`), matching `dpleph.c:182` and rejecting `jd > jd_end`.

### 2.4 Universal `barycenter_fits` (`bary.rs:294-563` vs `axBary.c:112-887`)
- **Fixes Verified (`#56`, `#85`, `#107`)**:
  1. **Primary `TSTART`/`TSTOP` & `DATE-OBS`/`DATE-END` (`bary.rs:340-350,454-551`)**: Preserves Primary HDU's own barycentered `TSTART`/`TSTOP` when present (avoiding inflation by auxiliary `VETOSPECTRUM` HDUs, `#85`), and populates `TSTART`/`TSTOP` together with `DATE-OBS`/`DATE-END` from `global_tstart`/`global_tstop` when absent (`#56`).
  2. **0-Row `EVENTS`/`GTI` Header Barycentering (`bary.rs:363,404,446,456,470`, Neon review fix `#107`)**: Initial code set `modified = true` even for `nrows == 0` while skipping `TSTART`/`TSTOP` updates, causing `if !modified` to leave header `TSTART`/`TSTOP`/`DATE-OBS`/`DATE-END` in uncorrected UTC on 0-row extensions labeled `TIMESYS='TDB'`. Fixed by tracking `updated_bounds` separately from `modified`.
  3. **GTI Bounds Folding & GTI-Only Fallback (`bary.rs:439-447,524-526`, `#107`)**: Folds `min`/`max` over `starts_tdb`/`stops_tdb` (safe for unsorted GTI rows) and propagates GTI bounds to `global_tstart`/`global_tstop` when no `TIME` extension exists.
  4. **Already-Barycentered Guard (`bary.rs:321-336,366-373`)**: Rejects input files whose principal HDU already has `TIMESYS = 'TDB'` or `'TCB'` and skips individual `TDB`/`TCB` extensions (`axBary.c:443-449`), preventing double-correction.
  5. **`TIMEZERO` Incorporation (`bary.rs:375,386-390,422-429,458-480`)**: Incorporates pre-existing `TIMEZERO` into `TIME` and `START`/`STOP` columns and uses `axBary.c:100-107`'s `timezero_heuristic_function` (`|timezero| < 0.1 * |TSTART|`) on non-table `TSTART`/`TSTOP` before zeroing `TIMEZERO = 0.0` (strictly more accurate than legacy `clock.c:340`, which dropped small `TIMEZERO` values for AstroSat!).
- **Remaining Gaps (`PARTIAL_GAP`)**:
  - `barycenter_fits` reads/writes `MJDREFI`/`MJDREFF` (`bary.rs:309-319,489-490`) but `BaryCalculator` hardcodes `ASTROSAT_MJDREFI = 55197, MJDREFF = 0.0` (`time.rs:7-9`) and `TIMEUNIT = 's'` (`axBary.c:455-468` supports `'d'` and arbitrary `MJDREF`).
  - Fine clock correction (`clock.c`, `tdc.dat`) is omitted (`clock.c:333-358` `astrosatClockCorr` itself returns `0.0` for all AstroSat instruments, so `EXACT_MATCH` for AstroSat).

---

## 3. AstroSat Core Modules (`tct.rs`, `gti.rs`, `mkf.rs`, `caldb.rs`)

- **`CztiTct` (`tct.rs:132-309`, `EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**: Subtracts `feb_lag = 2.0 s`, reduces modulo `2^20 = 1_048_576.0 s`, unwraps 20-bit counter rollovers when `(dt_inst + 1048576) - dt_sps` matches within `5%`, splits non-rollover onboard clock resets into monotonic `CztiSegment`s, and rejects non-finite inputs (`tct.rs:254-256`).
- **`LaxpcTct` (`tct.rs:322-490`, `EXACT_MATCH` / `INTENTIONAL_IMPROVEMENT`)**: Rejects telemetry glitches outside `dt_sps / dt_ins in (1.000, 1.001]` (`tct.rs:362,389`), anchors on the first valid pair, and provides `O(log N)` linear (`linear_interp`) and order-`m` Neville polynomial interpolation (`interpolate_neville`, shifted by `xs[0]` for numerical conditioning).
- **`GtiSet` (`gti.rs:36-314`, `EXACT_MATCH`)**:
  - `new`: filters `is_finite() && stop > start`, sorts by `(start, stop)`, merges overlapping/touching intervals (`iv.start <= last.stop`).
  - `intersect` / `subtract`: verified `O(N + M)` two-pointer sweeps (`subtract` keeps `j` monotonic across sorted `self` while scanning local `k >= j`), preserving sorted disjoint positive-duration invariant.
  - `merge_gaps`, `drop_short`, `interval_index` (closed `[start, stop]` via `partition_point(|iv| iv.stop < t)`), and `overlap_duration` (half-open overlap for `FRACEXP`) handle empty sets, zero/negative durations, and `NaN` queries (`NaN` returns `false`/`None`) cleanly.
- **`MkfTable` (`mkf.rs:14-194`, `EXACT_MATCH`)**: Case-insensitive column lookup, alias fallback, linear interpolation, and `gti_where` cadence-based GTI construction.
- **`CaldbIndex` (`caldb.rs:31-241`, `EXACT_MATCH` / `PARTIAL_GAP`)**:
  - Filters `cal_qual == 0`, case-insensitive `TELESCOP`, `INSTRUME`, `CAL_CNAM`, `DETNAM` (matching exact or `"NONE"`/`"-"`/`""`), and `ref_time_mjd <= obs_mjd + 1e-6`, sorting newest `ref_time_mjd` first with exact-`DETNAM` tie-break and earliest-entry fallback if observation predates all `CAL_VSD`.
  - Minor `PARTIAL_GAP`s: `tstart < 100_000.0` heuristic treats numbers `< 100,000` as MJD rather than AstroSat seconds (`caldb.rs:179`), `_tstop` (`caldb.rs:176`) is unused, and `parse_date_to_mjd` (`time.rs:342-361`) parses ISO-8601 (`YYYY-MM-DD[THH:MM:SS]`) rather than legacy `dd/mm/yy` strings.

---

## 4. FITS 4.0 & OGIP 1.2/1.3 Engine (`crates/astrosat-fits/src/`)

- **2880-Byte Block Alignment (`io.rs:83-294,406-479`, `card.rs:459-468,508-514`) — `EXACT_MATCH`**:
  - Validates `FITS_BLOCK_SIZE = 2880` alignment on read, checks `SIMPLE = T` on Primary and `XTENSION` (`BINTABLE`, `IMAGE`, `TABLE`) on extensions, rejects invalid `BITPIX`, negative/overflowing `NAXIS`/`NAXISi`/`PCOUNT`/`GCOUNT`, and verifies EOF bounds with checked arithmetic (`io.rs:140-152,254-293`).
  - Pads headers with `b' '` (`0x20`), binary/image data with `0u8`, and ASCII `TABLE` data with `b' '` (`io.rs:433-442`), with `0` data bytes for `NAXIS = 0`. (`PARTIAL_GAP`: FITS Random Groups `GROUPS=T, NAXIS1=0` not supported.)
- **Header Cards, 16-Sig-Fig Floats, Escaping & `CONTINUE` (`card.rs:90-515,644-714`) — `EXACT_MATCH`**:
  - `format_fits_float` (`card.rs:169-213`, upgraded in `#85`): preserves up to 16 significant digits in fixed-point (matching CFITSIO `%.15E` to `0.000 us` on `~2e8 s` AstroSat timestamps) and scientific notation. *(Minor precision note: for numbers in `[1e-4, 1e-1)`, `int_digits.max(1)` limits fixed-point decimals to 15, yielding 12–14 significant digits near `1e-4`.)*
  - Escapes single quotes `'` -> `''` (`card.rs:244-247`) and unescapes `''` -> `'` (`card.rs:651-653`).
  - Pads short single-card strings to `>= 8` characters inside quotes (`card.rs:256-258`).
  - Splits strings `> 68` escaped bytes into `CONTINUE` cards (`<= 67` bytes + `&`), guarding against splitting an escaped `''` pair across card boundaries (`card.rs:296-306`), automatically injects `LONGSTRN = 'OGIP 1.0'` before `DATASUM`/`CHECKSUM`/`END` (`card.rs:489-503`), parses Fortran `D` exponents (`card.rs:708`), and uses UTF-8-safe `char_indices()` (`card.rs:648-663`).
- **`DATASUM` & `CHECKSUM` (`checksum.rs:10-180`) — `EXACT_MATCH`**:
  - 32-bit big-endian ones'-complement word accumulation with end-around carry folding (`accumulate_checksum`, `add_ones_complement`, `checksum.rs:18-38,121-127`).
  - `compute_hdu_checksums` (`checksum.rs:132-157`) computes `DATASUM`, writes placeholder `CHECKSUM = '0000000000000000'` (`0xC0C0C0C0` contribution), sums 2880-padded header + `dsum`, and encodes `!total` via `encode_checksum` (`checksum.rs:50-97`), avoiding excluded ASCII punctuation (`0x3A..=0x40`, `0x5B..=0x60`) while preserving per-lane byte sums and right-rotating by 1 character for column-12 (`mod 4 = 3`) alignment per Seaman, Pence & Rots (2012).
  - Verified independently by NASA HEASARC `fitsverify 4.22` (`0 warnings, 0 errors` in `lib.rs:257-355` and across all generated L2 CZTI/LAXPC/barycentered products).
- **`BZERO`/`BSCALE`, `TZERO`/`TSCAL` & Unsigned Offsets (`bintable.rs:339-671,1149-1200,1274-1298,1543-1552`, `image.rs:154-275,380-465`) — `EXACT_MATCH` / `PARTIAL_GAP`**:
  - Scalar `u16` (`TZERO=32768, TSCAL=1`, `(v as i32 - 32768) as i16`), `u32` (`TZERO=2147483648, TSCAL=1`, `(v as i64 - 2_147_483_648) as i32`), exact unsigned 64-bit `TZERO = 2^63` wrapping (`bintable.rs:438-493`, fixed `#107`), and vector `u16` (`read_vec_u16`, `bintable.rs:649-671`, `wrapping_add` `#107`) preserve exact integer arithmetic without overflow or `f64` rounding.
  - `ImageRef::read_f64` applies `bscale * v + bzero` across all `Bitpix` variants (`8, 16, 32, 64, -32, -64`).
  - Minor `PARTIAL_GAP`: `read_vec_{i16,i32,f32,f64,u8}` and `read_var_*` do not apply `TSCAL`/`TZERO` (only scalar readers and `read_vec_u16` apply `TZERO`).
- **Variable-Length `P`/`Q` Heap Arrays (`bintable.rs:122-138,176-214,774-940,1491-1522`) — `EXACT_MATCH` / `PARTIAL_GAP`**:
  - Parses `rPt(emax)` (`P32`: 8-byte `(u32 count, u32 offset)`) and `rQt(emax)` (`Q64`: 16-byte `(u64 count, u64 offset)`), resolves offsets relative to `THEAP` (default `NAXIS1 * NAXIS2`) with overflow-checked `count * elem_size`, `heap_off + nbytes <= pcount`, and `start + nbytes <= data.len()` bounds validation.
  - `BinTableBuilder` writes big-endian `P32` `(count, heap_off)` descriptors and appends the heap buffer with `PCOUNT = heap_buf.len()`. *(Minor `PARTIAL_GAP`: writer emits `P32` descriptors only; `Q64` is supported on read.)*
- **OGIP 1.2/1.3 Science Builders & Validator (`ogip.rs:10-539`) — `EXACT_MATCH`**:
  - Produces compliant `GTI` (`HDUCLASS='OGIP'`, `HDUCLAS1='GTI'`, `HDUCLAS2='STANDARD'`), `SPECTRUM` (OGIP-92-007, `HDUVERS='1.2.0'`), `LIGHTCURVE` (`HDUCLAS1='LIGHTCURVE'`, `TIMEDEL`, `TIMEPIXR=0.5`, `FRACEXP`), `RESPONSE` (`MATRIX` `HDUVERS='1.3.0'` + `EBOUNDS` `HDUVERS='1.2.0'`), `ORBIT`, and `TAN` WCS headers.
