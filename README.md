<p align="center">
  <img src="assets/banner.png" alt="astrosat-rs — Pure-Rust ISRO AstroSat CZTI &amp; LAXPC Pipelines" width="100%">
</p>

<h1 align="center"><code>astrosat-rs</code></h1>

<p align="center">
  <strong>Zero-patch, multi-threaded pure-Rust reduction pipelines, Solar System barycentering, and an embedded offline web workbench (AstroSat Studio) for ISRO AstroSat CZTI &amp; LAXPC.</strong>
</p>

---

## Why `astrosat-rs`?

Traditional AstroSat Level-1 to Level-2 workflows require multi-gigabyte legacy stacks (`HEASoft`, `CFITSIO`, `glog`, `pil`, `gfortran`, `preplaxbary`, `$HEADAS`) that are fragile to install and slow to batch across archives.

`astrosat-rs` replaces the entire chain with a **single ~8.4 MB standalone binary** linked solely against standard system libraries (`libc`, `libm`, `libgcc_s`):

- **End-to-End Scientific Pipelines:** Full `rayon`-parallelized **CZTI** (Stages 1–3 + `addspec`) and **LAXPC** (Units 10/20/30 + `backshift` + 2D B-spline background) Level-1 $\to$ Level-2 reduction, plus native **Solar System barycentric correction** (analytical & JPL DE200/405/421/430 binary or HEASARC FITS `JPLEPH.*`).
- **AstroSat Studio (`astrosat ui`):** An embedded, 100% offline browser workbench compiled directly into the binary — zero Node, Python, or external CDNs — with full CLI-parity pipeline execution, two-way `.toml` config round-trip, interactive light curves/spectra, WCS 2D imaging, and multi-band pulse alignment.
- **Zero-Copy FITS 4.0 & OGIP Engine:** Memory-mapped (`memmap2`) I/O with Seaman/Pence `CHECKSUM`/`DATASUM` and OGIP-compliant products verified with NASA HEASARC `fitsverify 4.22` (`0` errors, `0` warnings).
- **Publication-Grade Accuracy:** Independently cross-checked against the official **Jodrell Bank Observatory Crab Pulsar radio ephemeris** to **$3\text{–}8\text{ ppb}$** in spin frequency and **$0.00044\%$ ($4.4\text{ ppm}$)** in spin-down rate over 2.77 years.

---

## Quickstart

### Prerequisites
Only a Rust toolchain (`rustc` + `cargo` via [rustup.rs](https://rustup.rs)) — **no C/C++/Fortran scientific libraries required**.

### Install (No `sudo` Required)

```bash
git clone https://github.com/meural-operator/astrosat-rs.git
cd astrosat-rs

# Installs to $HOME/.local/bin/astrosat and verifies shared-library linkage via ldd
./scripts/install.sh

# Or build directly with Cargo
cargo build --release -p astrosat-cli
./target/release/astrosat --help
```

---

## AstroSat Studio — Embedded Web Workbench

Launch the full graphical workbench from the same binary — ideal for interactive exploration, teaching, or running over SSH port-forwarding on a remote server:

```bash
# Starts AstroSat Studio at http://127.0.0.1:8080/ using the current directory (Classic UI at /classic)
astrosat ui

# Optional: point the file browser & auto-discovery at a specific data folder or custom port
astrosat ui --data-dir /path/to/astrosat/data --port 8086 --no-open
```

> [!TIP]
> Every form in **AstroSat Studio** live-generates its exact `astrosat` shell command and supports two-way **Load / Download `.toml`** configuration files, so interactive sessions remain 100% reproducible and ready for HPC batch scripting.

- **Pipelines & Tools:** Run `czti run`, `laxpc run`, `bary`, `czti addspec` (multi-quadrant / multi-orbit spectrum & response co-adding), and `laxpc backshift` in-process, with collapsed **Advanced options** drawers for fine-grained CLI overrides.
- **Interactive Analysis & Export:** Rebinned light curves with GTI shading, channel/keV spectra, 2D CZTI DPH/DPI/sky images with live WCS coordinates, multi-band phase-folded pulse alignment + $Z^2_n$ periodogram search + circular cross-correlation phase lags, strict FITS/OGIP inspection, Light/Dark themes, and SVG/PNG export.

---

## CLI Usage

Every subcommand supports `--threads N`, `--json` (structured stdout), `--config <file.toml|json>`, and `--help`.

<details open>
<summary><strong>Core Pipeline, Barycentering &amp; Timing Examples</strong></summary>

```bash
# 1. CZTI Level-1 -> Level-2 (Stages 1-3: events, GTIs, sky images, spectra, light curves, RMFs)
astrosat czti run \
  --input <DATA_DIR>/czti/20160331_T01_112T01_9000000406_level1 \
  --caldb <CALDB_DIR>/czti \
  --output ./out_czti --threads 16

# Co-add multi-quadrant or multi-orbit CZTI spectra & exposure-weighted responses
astrosat czti addspec \
  --pha ./out_czti/*_Q0.pha ./out_czti/*_Q1.pha \
  --rsp ./out_czti/*_Q0.rsp ./out_czti/*_Q1.rsp \
  -o ./combined.pha --out-rsp ./combined.rsp

# 2. LAXPC Level-1 -> Level-2 (Units 10/20/30, Format A + B, and AstroSat.orb)
astrosat laxpc run \
  --input <DATA_DIR>/laxpc_l1/20160107_G02_010T01_9000000252_level1 \
  --caldb <CALDB_DIR>/laxpc \
  --output ./out_laxpc --format-a --format-b --orbit --threads 16

# 3. Solar System Barycentric Correction (JPL DE405 FITS/binary, zero preplaxbary hacks)
astrosat bary \
  -i ./out_laxpc/lxp1level2.evn \
  -o ./out_laxpc/lxp1level2_bary.evn \
  --orbit ./out_laxpc/AstroSat.orb \
  --ra 83.63308333 --dec 22.01450000 --frame icrs \
  --ephem <EPHEM_DIR>/JPLEPH.405

# 4. FITS Inspection, Strict OGIP Validation & Z^2_n Pulsar Periodicity Search
astrosat inspect ./out_laxpc/lxp1level2_bary.evn --full
astrosat validate ./out_czti/*.pha --strict
astrosat timing ./out_laxpc/lxp1level2_bary.evn --freq 29.658 --search-span 0.01 --harmonics 2
```

Individual micro-stages are also exposed directly (`czti decode|bunchclean|pha2energy|gtigen|gaas|datasel|clean|evtencor|evtclean|dpigen|image|bindata|addspec` and `laxpc decode|gtigen|lightcurve|spectrum|backshift|orbit`).

</details>

---

## Science Validation on Real AstroSat Data

To validate end-to-end physical fidelity — not just software regressions — `astrosat-rs` was tested on **29.34 million real photons** from public ISRO AstroSat Crab Pulsar (`PSR B0531+21`) observations spanning **2.77 years** (`2016-01-08` LAXPC ObsID `9000000252`, `2016-04-02` CZTI ObsID `9000000406`, and `2018-10-14` CZTI ObsID `9000002432`).

### 1. Crab Pulsar 2.77-Year Spin-Down vs. Jodrell Bank Radio Observatory

All three epochs were barycentered with `astrosat bary` (`JPLEPH.405`, ICRS $\alpha = 83.63308333^\circ$, $\delta = +22.01450000^\circ$) and measured with Fourier-oversampled $Z^2_2$ searches (`astrosat timing`), then compared directly against the official [Jodrell Bank Observatory Crab Pulsar Monthly Ephemeris (`crab2.txt`)](https://www.jb.man.ac.uk/~pulsar/crab/crab2.txt):

| Epoch (UTC / MJD TDB) | Instrument / ObsID | Photons | Jodrell Bank Radio $f(t)$ | `astrosat-rs` X-Ray $f(t)$ | Offset | Agreement |
| :--- | :--- | ---: | :--- | :--- | :--- | :--- |
| **2016-01-08** (`57395.0483`) | LAXPC 10/20/30 (`9000000252`) | 25.82 M | `29.658058911 Hz` | `29.658059000 ± 4.5e-7 Hz` | **`+89 nHz`** (`3.0 ppb`) | **`+0.20 σ`** |
| **2016-04-02** (`57480.8605`) | CZTI Q0–Q3, 367 s (`9000000406`) | 0.19 M | `29.655321228 Hz` | `29.654804000 ± 2.1e-4 Hz` | `-0.517 mHz` *(6-min span)* | `-2.47 σ` |
| **2018-10-14** (`58405.4327`) | CZTI Q0–Q3 (`9000002432`) | 3.33 M | `29.625857771 Hz` | `29.625858000 ± 2.1e-6 Hz` | **`+229 nHz`** (`7.7 ppb`) | **`+0.11 σ`** |

Across the $\Delta t = 1{,}010.384\text{ day}$ ($8.730 \times 10^7\text{ s}$) LAXPC $\to$ CZTI baseline, the Crab slowed by $\Delta f = -32.2010 \pm 0.0022\text{ mHz}$ ($>14{,}800\sigma$):

$$
\dot{f}_{\text{AstroSat}} = (-3.688663 \pm 0.000248) \times 10^{-10}\ \text{Hz s}^{-1}
\qquad \text{vs.} \qquad
\dot{f}_{\text{Jodrell Bank}} = -3.688679 \times 10^{-10}\ \text{Hz s}^{-1}
\quad (\mathbf{0.00044\%}\text{ diff})
$$

Derived neutron-star parameters recovered directly from AstroSat photons:
- **Period derivative:** $\dot{P} = -\dot{f}/f^2 = (+4.19812 \pm 0.00028) \times 10^{-13}\ \text{s s}^{-1}$ ($+36.27\ \text{ns day}^{-1}$)
- **Characteristic age:** $\tau_c = -f / (2\dot{f}) = 1{,}273\ \text{yr}$ *(consistent with SN 1054 and braking index $n \approx 2.51$)*
- **Spin-down luminosity:** $\dot{E} = -4\pi^2 I f \dot{f} = 4.32 \times 10^{38}\ \text{erg s}^{-1}$ *(for $I = 10^{45}\ \text{g cm}^2$)*
- **Surface dipole field:** $B_s = 3.2 \times 10^{19}\sqrt{P\dot{P}} = 3.81 \times 10^{12}\ \text{G}$

### 2. Five Additional Real-Data Astrophysics Demonstrations

| Experiment | Real Dataset | Measured Result | Literature / Physics Check |
| :--- | :--- | :--- | :--- |
| **Broadband 3–200 keV Pulse Evolution** | LAXPC10 (8.8 M) + CZTI (3.2 M) | Interpulse-to-main-pulse ratio $\text{P2}/\text{P1}$ rises from `0.900 ± 0.011` (`3–6 keV`) to `1.230 ± 0.045` (`35–55 keV`) and `1.05–1.67` (`100–200 keV`); $\text{Bridge}/\text{P1}$ triples from `0.316` to `1.063`; peak separation holds at $\Delta\phi \approx 0.406$ cycles. | Matches Kuiper et al. (2001) magnetospheric pulse-shape evolution. |
| **Phase-Resolved Spectral Hardening** | LAXPC10 (`3–80 keV`) + CZTI (`20–200 keV`) | After subtracting Off-Pulse Crab Nebula emission, LAXPC hardness $(20\text{–}80)/(3\text{–}10\ \text{keV})$ orders strictly as **P1** (`0.2115 ± 0.0029`) $<$ **P2** (`0.2582 ± 0.0035`, $+22\%$, $10.3\sigma$) $<$ **Bridge** (`0.3787 ± 0.0139`, $+79\%$, $11.8\sigma$); CZTI (`100–200`/`20–50 keV`) replicates **P1** (`0.211`) $<$ **P2** (`0.242`) $<$ **Bridge** (`0.332`). | Confirms the Bridge is the spectrally hardest region of the Crab magnetosphere. |
| **Multi-Detector Clock Synchronization** | LAXPC 10/20/30 (25.0 M) & CZTI Q0–Q3 (3.2 M) | Cross-correlation offsets: `LX20 − LX10 = +0.3 µs` ($r = 0.99915$), `LX30 − LX10 = −4.3 µs` ($r = 0.99911$); CZTI `Q1–Q3` aligned to `Q0` within `< 92 µs` ($r > 0.90$). | Sub-$5\ \mu\text{s}$ inter-unit synchronization across LAXPC. |
| **0–250 Hz High-Frequency FFT Comb & Deadtime** | LAXPC10 (3.91 M, $\Delta t = 2\ \text{ms}$, 80 segments) | Cleanly detects all harmonics $1f$–$8f$ (`29.66`–`237.26 Hz`) up to **`857.9 σ`** (`2f`), showing destructive dip at `4f` (`31 σ`) and constructive surge at `5f` (`238.5 σ`, since $5 \times 0.40 = 2.0$ cycles); `150–240 Hz` Leahy floor = **`1.551`**. | Matches Zhang et al. (1995) deadtime-suppressed white-noise floor $2(1 - R\tau_d)^2 = 1.520$ for $\tau_d = 43\ \mu\text{s}$. |
| **100–380 keV CZTI Compton Polarimetry Prototype** | CZTI Oct 2018 (3.33 M, 5.9 ks) | Isolated `806` adjacent 2-pixel coincidences ($\Delta t \le 40\ \mu\text{s}$) and `510` Compton-kinematics pairs ($1 \le E_{\text{high}}/E_{\text{low}} \le 6$); `308` Edge vs. `202` Corner pairs give **Edge/Corner = `1.52`**. | Reproduces the ~1.5–1.6 square-pixel geometry factor of Vadawale et al. (2018, *Nature Astronomy*). |

---

## Architecture (7 Crates)

| Crate | Responsibility | Docs |
| :--- | :--- | :--- |
| `astrosat-fits` | Zero-copy `memmap2` FITS 4.0 reader/writer, `CONTINUE`, `CHECKSUM`/`DATASUM`, OGIP builders | [architecture.md](docs/architecture.md) |
| `astrosat-core` | AstroSat `MJDREFI=55197` time/leap-seconds, TCT, `GtiSet`, Taylor orbit, `as1bary` SSB engine, CALDB CIF, shared TOML config | [architecture.md](docs/architecture.md) |
| `astrosat-czti` | 3-stage `rayon`-parallelized CZTI pipeline (`science2event` $\to$ `bindata`/`.rsp` + `addspec`) | [czti_pipeline.md](docs/czti_pipeline.md) |
| `astrosat-laxpc` | Multi-unit `rayon`-parallelized LAXPC 10/20/30 pipeline + 2D B-spline background & RMF selection | [laxpc_pipeline.md](docs/laxpc_pipeline.md) |
| `astrosat-ui` | Embedded **AstroSat Studio** web workbench (`/` default, `/classic` preserved) with CLI-parity jobs, TOML round-trip, WCS imaging, and pulse alignment | [`crates/astrosat-ui`](crates/astrosat-ui/src/lib.rs) |
| `astrosat-cli` | Unified `astrosat` CLI (`czti`, `laxpc`, `bary`, `inspect`, `validate`, `timing`, `ui` — 29 subcommands) | [migration_from_legacy.md](docs/migration_from_legacy.md) |
| `astrosat-bench` | Reproducible latency, throughput, `rayon` scaling, and `VmHWM` peak-RSS benchmark suite | [`crates/astrosat-bench`](crates/astrosat-bench/src/main.rs) |

---

## Legacy Parity, FITS Compliance & Benchmarks

`astrosat-rs` has been numerically audited against untouched **IUCAA CZTI Pipeline v3.0.1**, **LaxpcSoft v3.4.6**, **NASA/IUCAA `as1bary`**, and **NASA HEASARC `fitsverify 4.22`**:

- **Barycentering:** `0–29.8 ns` (`0–1` IEEE-754 `f64` ULP) vs. `as1bary` across CZTI, LAXPC Format A/B, leap-second boundaries (`0.00 ns`), and near-limb solar conjunction (`-107.81 µs` Shapiro delay, `0.00 ns` diff); **7×–15× faster** ([full report](audit/heasoft_as1bary_numerical_verification.md)).
- **CZTI L1 $\to$ L2:** `0`-event / `0`-bunch / `0.000000 s` GTI diff vs. v3.0.1; spectra match to `7.9e-7% – 2.3e-6%` ($r > 0.9999999999999$) while fixing legacy `float32` temperature cancellation and `cztimage` UB crashes ([full report](audit/czti_head2head_verification.md)).
- **LAXPC L1 $\to$ L2:** `0`-event / `0`-channel diff across `25.8M` Crab photons; bit-exact `1.0 s` light curves; `backshift` matches `backshiftv3.e` to `0 – 8.88e-16 cps` ([full report](audit/laxpc_head2head_verification.md)).
- **FITS & Adversarial Safety:** `149/149` generated FITS products pass `fitsverify 4.22` (`0` errors, `0` warnings); `0` panics across 69 corrupt/adversarial inputs ([full report](audit/public_release_scenario_report.md)).

<details>
<summary><strong>Reproducible Performance Benchmarks (52-vCPU Host, <code>./scripts/run_benchmarks.sh</code>)</strong></summary>

| Benchmark | 1T (ms) | 4T (ms) | 16T (ms) | 52T (ms) | Best Speedup | Peak RSS (MiB) | Workload |
| :--- | ---: | ---: | ---: | ---: | ---: | ---: | :--- |
| `bary_synth` | 6402.96 | 2029.18 | 601.79 | 252.24 | **25.38×** | 35.6 | 1M timestamps (`3.96M events/s` @ 52T) |
| `bary_fits_real` | 1680.83 | 810.17 | 606.76 | 565.99 | **2.97×** | 76.9 | ObsID `9000000406` `.evt` + `.orb` (11 HDUs) |
| `fits_io_real` | 763.62 | 513.07 | 532.88 | 552.77 | **1.49×** | 72.4 | CZTI + LAXPC L1 BINTABLEs (`187 MB/s` @ 4T) |
| `fits_io_synth` | 298.87 | 307.74 | 279.73 | 299.13 | **1.07×** | 72.6 | 1M rows, 17.17 MiB R+W, checksummed (`3.57M rows/s`) |
| `czti_real_9000000406` | 9457.48 | 6840.34 | 6826.38 | 6746.44 | **1.40×** | 678.9 | Stages 1–3: `2.40M` raw $\to$ `166k` clean + `.pha`/`.lc`/`.rsp` |
| `laxpc_real_9000000252` | 11899.09 | 7421.79 | 6870.03 | 6808.10 | **1.75×** | 2323.6 | LX10/20/30: `28.46M` Crab events + Format A/B + `.orb` |

</details>

---

## Running Tests

```bash
# Full unit & synthetic integration test suite (works offline without external data)
cargo test --workspace

# Include real-observation & full CLI<->UI parity integration tests
ASTROSAT_TEST_DATA=/path/to/astrosat/data cargo test --workspace
ASTROSAT_TEST_DATA=/path/to/astrosat/data cargo test -p astrosat-ui --test e2e_parity -- --include-ignored
```

---

## License

Dual-licensed under either **MIT** ([LICENSE-MIT](LICENSE-MIT)) or **Apache-2.0** ([LICENSE-APACHE](LICENSE-APACHE)), at your option.
