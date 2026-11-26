# Migration Guide: Legacy ISRO Tools to `astrosat-rs`

This guide maps legacy AstroSat **CZTI Pipeline v3.0.1**, **LaxpcSoft v3.4.6**, and **HEASoft/`as1bary` (`+ preplaxbary`)** workflows onto `astrosat-rs` (`astrosat czti run`, `astrosat laxpc run`, `astrosat bary`, and the TOML/JSON configuration schema in [`crates/astrosat-cli/src/config.rs`](../crates/astrosat-cli/src/config.rs)).

---

## 1. CZTI Pipeline v3.0.1 $\to$ `astrosat czti run`

### 1.1 Eliminating Build & Environment Fragility
The legacy `czti_pipeline` (`README`, `InstallLibs`, `czti/makefile.common`) requires:
* Running `./InstallLibs` to configure, build, and install three bundled C/C++ libraries from source (`LIBS/glog_modified-0.3.2`, `LIBS/cfitsio`, and `LIBS/pil`), injecting `CXXFLAGS="-fpermissive"` and `-std=c++14` / `-std=c++1y` workarounds for modern GCC.
* Compiling 15+ separate C++/Fortran/C binaries (`make` in `czti/`) and setting 7+ shell environment variables (`as1czt`, `PFILES`, `PATH`, `GLOG_log_dir`, `CZTI_templates`, `LD_LIBRARY_PATH`, `PYTHONPATH`, `CALDB`, plus `ulimit -s 65532`).
* Managing `PFILES` parameter state under `$HOME/pfiles` to avoid prompt collisions with HEASoft's `pil`.

In `astrosat-rs`, all FITS I/O, gzipped CALDB decompression (`flate2`), logging (`tracing`), parameter configuration, and all 3 CZTI pipeline stages (plus `addspec`) compile into a single pure-Rust binary with **no CFITSIO, glog, PIL, or `CZTI_templates` runtime dependency**.

### 1.2 Parameter Mapping: `cztpipeline` / `cztpipeline.par` $\to$ `astrosat-czti` / `[czti]`

Legacy `cztpipeline` prompts or parses `cztpipeline.par` parameters (`entrystage`, `exitstage`, `l1dir`, `l2dir`, `usergtifile`, `outputtype`, `maskWeight`, `rasrc`, `decsrc`, `timebinsize`, `emin`, `emax`, `evtselalgo`, `gtitype`, `mkfthresholdfile`). In `astrosat-rs`:

| Legacy Parameter (`cztpipeline.par` / Module) | `CztiPipelineConfig` / TOML `[czti]` Key | Default / Notes |
| :--- | :--- | :--- |
| `entrystage` / `exitstage` | `entry_stage` / `exit_stage` | `1` / `3` (`1..=3`); explicit CLI flags always override config values |
| `l1dir` (`*cztM0_level1.fits`, `*SS*`, `.tct`, `.mkf`) | `input_dir` (or explicit `l1_m0_file`, `l1_ss_file`, `tct_file`, `mkf_file`) | Explicit files or L1 directory auto-discovery |
| `l2dir` | `output_dir` | Output directory created automatically |
| `CALDB` env var | `caldb_dir` | Path containing `data/as1/czti/caldb.indx` |
| `mkfthresholdfile` (`${as1czt}/config/mkfThresholds.txt`) | `mkf_threshold_file` | Optional; built-in v3.0.1 defaults (`ELV >= 5`, `ANG_OFFSET <= 0.03`, `Q1..4_MODEID == 0`, veto limits) |
| `usergtifile` | `user_gti_file` | Optional two-column ASCII `tstart tstop` |
| `gtitype` (`QUAD` / `COMMON`) | `gti_type` (`"quad"` / `"common"`) | `GtiType::Quad` (`Common` auto-combines `.img`/`.pha`/`.lc`/`.rsp`) |
| `evtselalgo` (`default` / `generalized`) | `algo` / `evt_sel_algo` | `EventSelectionAlgo::Default` or `Generalized` |
| `rasrc` / `decsrc` | `ra_deg` / `dec_deg` (`Option<f64>`) | `None` uses MKF `RA_PNT`/`DEC_PNT`; explicit `--ra 0 --dec 0` targets `(0, 0)` cleanly |
| `maskWeight` | `mask_weight` | `true` |
| `outputtype` (`lc`/`spec`/`both`) | `output_type` / `BinDataOutputType` | `"both"` (`Spec`, `Lc`, `Both`) |
| `timebinsize` / `emin` / `emax` | `time_binsize` / `emin_kev` / `emax_kev` | `1.0 s` / `20.0 keV` / `100.0 keV` |
| `genresp` / `generate_eventfile` | `gen_resp` / `generate_eventfile` | `true` / `false` |
| `algorithm` / `oversamplingfactor` | `image_algo` / `oversampling` | `"balcc"` / `16` |
| `quadsToProcess` | `quads` | `[0, 1, 2, 3]` |
| `cztaddspec` | `astrosat czti addspec --pha ... [--rsp ...] -o <out.pha> [--out-rsp <out.rsp>] [--add-exposure]` | Co-adds spectra & exposure-weighted RMFs |

**Example TOML configuration (`czti.toml`)** matching [`CztiFileConfig`](../crates/astrosat-cli/src/config.rs):
```toml
threads = 16

[czti]
l1_m0_file = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1/modeM0/AS1T01_112T01_9000000406cztM0_level1.fits"
l1_ss_file = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1/modeSS/AS1T01_112T01_9000000406cztSS_level1.fits"
tct_file = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1/aux/AS1T01_112T01_9000000406czt_level1.tct"
mkf_file = "/mnt/d/ISRO/data/sample/czti/czti_sample_data/20160331_T01_112T01_9000000406_level1/AS1T01_112T01_9000000406czt_level2.mkf"
caldb_dir = "/mnt/d/ISRO/data/caldb/czti"
output_dir = "./czti_out"
entry_stage = 1
exit_stage = 3
gti_type = "quad"
algo = "default"
output_type = "both"
image_algo = "balcc"
```
Run with `astrosat czti run --config czti.toml` (or pass flags directly on the CLI).

---

## 2. LaxpcSoft v3.4.6 $\to$ `astrosat laxpc run`

### 2.1 Eliminating Fortran Compilation, Endianness & Fixed-Array Traps
Legacy `laxpcl1.f` and `backshiftv3.f` (`readme.laxpcl1`, `readme.backshiftv3`, `makefile.gfortran`, `makefile.ifort`) suffer from several operational hazards:
* **Compiler-Specific `recl=2880`**: Requires `gfortran -w -O3` (`makefile.gfortran`) or `ifort -O3 -assume byterecl` (`makefile.ifort`); otherwise `ifort` interprets `recl=2880` in 4-byte words (`11520` bytes) and fails.
* **Host Endianness Dependence**: `laxpcl1.f` manually byte-swaps FITS records after reading raw `direct` unformatted blocks and breaks on big-endian architectures. `astrosat-fits` explicitly decodes big-endian FITS bytes via `from_be_bytes`, making it portable across all host endiannesses.
* **Shell `call system('rm ...')` & Static Array Bounds**: `laxpcl1.f` shells out `call system('rm '//fil1//' 2> /dev/null')` before writing FITS files and hard-caps memory at `nn=1000000`, `ns=9740000`, `nl=10000000` (silently truncating high-resolution light curves). `astrosat-laxpc` uses dynamic Rust vectors and writes FITS files cleanly without spawning subshells.
* **Two-Pass Interactive `.inp` Workflow**: Legacy users had to edit `findfile` to build `laxpcl1.inp` (4 files per orbit in strict order: `.mkf`, `.tct`, `modeBB`, `modeEA` for LAXPC10 only), create a dummy `gti.inp` with `0` on line 1, run `laxpcl1.e` once with `tbin=1.0`, copy `lxp1level2.gti` to `gti.inp`, run `laxpcl1.e` a second time for each detector `lxp = 1, 2, 3`, and then run `backshiftv3.e` with fixed-column `backshiftv3.inp`.

### 2.2 Side-by-Side Workflow & Input Mapping

| Legacy LaxpcSoft Step | `astrosat-laxpc` / `astrosat laxpc run` Equivalent |
| :--- | :--- |
| `sh findfile <obsID>` + manual `laxpcl1.inp` (4 files/orbit, `lxp1` only) | `input_dir` auto-discovers all orbits & all units via `discover_l1_orbits` (or explicit `orbits_by_unit`) |
| Interactive prompts `lxp, tbin, ian`, `icl, icu`, `nca`, `nul, je, iev` | Single `PipelineConfig` / TOML `[laxpc]` (`units`, `tbin`, `spectrum_channels`, `EventFilter`) + symmetric CLI overrides (`--format-a`/`--no-format-a`, etc.) |
| Two-pass `gti.inp` (`0` -> `lxp1level2.gti` -> `gti.inp`) | Single-pass automatic `generate_standard_gti` (with optional `user_gti_file` / `user_gti_only`) |
| Running `laxpcl1.e` 3 separate times (`lxp = 1, 2, 3`) | Processes `LAXPC10`, `LAXPC20`, `LAXPC30` concurrently in one `rayon` invocation |
| `backshiftv3.e` + fixed-format `backshiftv3.inp` | Integrated Stage 6 (`caldb_dir` + `background_model_path` + `background_pha_path` + automatic `BackgroundCatalog` / `--bkg-peak-ch` override) or standalone `astrosat laxpc backshift` |
| **Format B** (`lxp?level2.{evn,spec,lcurv,gtis}`) vs **Format A** (`level2.event.fits`) | `write_format_b_events = true` and/or `write_format_a_combined = true` (uniform 1024-ch default scale) in the same run, plus `AstroSat.orb` (`write_orbit_file = true`) |

**Example TOML configuration (`laxpc.toml`)** matching [`LaxpcFileConfig`](../crates/astrosat-cli/src/config.rs):
```toml
threads = 16

[laxpc]
input_dir = "/mnt/d/ISRO/data/sample/laxpc_l1/20160107_G02_010T01_9000000252_level1"
output_dir = "./laxpc_out"
caldb_dir = "/mnt/d/ISRO/data/caldb/laxpc"
units = ["lxp1", "lxp2", "lxp3"]
tbin = 1.0
deadtime_sec = 0.000043
write_format_b_events = true
write_format_a_combined = true
write_orbit_file = true
single_only = false
layer1_only = false
```
Run with `astrosat laxpc run --config laxpc.toml`.

---

## 3. `as1bary` + `preplaxbary` + HEASoft $\to$ `astrosat bary`

### 3.1 Why Legacy `as1bary` Needed HEASoft & `preplaxbary`
Inspecting `as1bary`'s `axBary.c`, `clock.c`, and `README` reveals why legacy barycentric correction required brittle setup:
1. **HEASoft Dependency**: `as1bary` links against HEASoft (`Makefile` uses `HD_CTASK`, `HD_STD_CLIBS`, `cfitsio`) and requires `$HEADAS/headas-init.sh` plus `$TIMING_DIR` / `$ASC_DATA` containing `tai-utc.dat` and `JPLEPH.405`/`JPLEPH.200`.
2. **The `RA_NOM` / `MJDREFI` Bug (`axBary.c:295-302`)**: If `RA_NOM`/`DEC_NOM` is missing in the principal HDU, `axBary.c` falls back to `RA_PNT`/`DEC_PNT` and explicitly resets `mjdrefi = 0.0`, then aborts if `EQUINOX != 2000.0` or `RADECSYS` is missing. Legacy LAXPC Level-2 event files often lacked `MJDREFI` in HDU 0 or `RA_NOM`/`RADECSYS`, forcing users to run `preplaxbary` (or manual `fparkey`/`fthedit`) before `as1bary`.
3. **In-Place File Mutation & Edge-Case Bugs**: Omitting `-o outputfile` in `as1bary` overwrites the input event file in place; `clock.c` drops small non-zero `TIMEZERO` offsets, and `axBary.c` double-counts `+66.184 s` when formatting `DATE-OBS`/`DATE-END`.

### 3.2 `astrosat bary` Improvements
Implemented in [`astrosat_core::bary::barycenter_fits`](../crates/astrosat-core/src/bary.rs) and configured via [`BaryFileConfig`](../crates/astrosat-cli/src/config.rs):
* **Automatic HDU-0 / HDU-1 Keyword Resolution**: Checks HDU 1 first and falls back to HDU 0 for `MJDREFI`/`MJDREFF` (defaulting to `ASTROSAT_MJDREFI = 55197`, `0.0`), handling CZTI `.evt`, LAXPC Format B `.evn`/`.lcurv`, and LAXPC Format A `level2.event.fits` natively—**zero `preplaxbary` needed**.
* **Dual Ephemeris Backends**: Uses `EphemerisKind::Analytical` (`PLEPHEM = 'ASTROSAT-ANALYTICAL'`, `~4–9 ms` vs DE405) out of the box with no external `tai-utc.dat` or `JPLEPH.*` files required, OR accepts **BOTH** official raw JPL binary ephemerides (`lnxp*.405`, DE200/DE405/DE421/DE430) **AND** NASA HEASARC FITS-wrapped ephemerides (`JPLEPH.405` / `JPLEPH.200`) for **`<= 29.8 ns` (1 ULP)** parity with NASA/IUCAA `as1bary`.
* **Non-Destructive, Guarded, Parallel, Checksummed Output**: Rejects already-barycentered files (`TIMESYS = 'TDB'` / `'TCB'`), properly folds any pre-existing non-zero `TIMEZERO`, rewrites `TIME` and all `GTI` `START`/`STOP` columns in parallel via `rayon`, sets `TIMESYS = 'TDB'`, `TIMEREF = 'SOLARSYSTEM'`, `TREFPOS = 'BARYCENTER'`, `TREFDIR = 'RA_OBJ,DEC_OBJ'`, `RA_OBJ`, `DEC_OBJ`, `RADESYS`/`RADECSYS`, `PLEPHEM`, `TSTART`/`TSTOP`, `DATE-OBS`/`DATE-END`, and updates `DATASUM`/`CHECKSUM`.

```toml
# Example TOML configuration for `astrosat bary --config bary.toml`
[bary]
orbit_file = "./laxpc_out/AstroSat.orb"
ra_deg = 83.6331
dec_deg = 22.0145
frame = "ICRS"           # "ICRS" or "FK5"
ephemeris = "analytical" # "analytical", path to raw JPL binary (lnxp*.405), or HEASARC FITS (JPLEPH.405/200)
write_checksums = true
```
