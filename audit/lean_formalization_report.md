# Lean 4 Formal Verification Report — `astrosat-rs`

## 1. Executive Summary

A self-contained, pure Lean 4 (`Init` + `Std.Tactic.BVDecide`, zero Mathlib dependencies) formal verification package [`AstrosatFormal`](file:///d:/ISRO/formal/lakefile.toml) was constructed in [`d:/ISRO/formal/`](file:///d:/ISRO/formal/) mirroring the exact bitmasks, coordinate reflections, interval algorithms, checksum arithmetic, deadtime formulas, mask-weight renormalization equations, leap-second table, and rollover logic of `astrosat-rs`.

- **Toolchain**: Lean `v4.36.0` (`611c117bfaf7b7ae1cfcf246ce6fba31e8f509e9`), Lake `5.0.0-src+611c117`, installed strictly under `/mnt/d/ISRO/.toolchain/elan` (`3.1 GiB` on `D:`).
- **Package Size**: `9.0 MiB` (including `.lake/build` artifacts) under `/mnt/d/ISRO/formal`.
- **Storage Audit**: WSL2 root `/dev/sdd` unchanged at `71 GiB` used (`886 GiB` free); Windows `C:` free space `11.36 GiB` (vs `10.57 GiB` before installation — zero storage consumed on `C:`); `D:` has `424.05 GiB` free.
- **Verification Status**:
  - `lake clean && lake build`: **12/12 jobs succeeded, 0 errors, 0 warnings** (`15s` full clean build).
  - `sorry` / custom `axiom` scan across `AstrosatFormal.lean` and `AstrosatFormal/*.lean`: **0 `sorry`, 0 `axiom`**.
  - **126 machine-checked theorems** across 4 modular proof files (`Telemetry.lean`: 18, `GtiAlgebra.lean`: 48, `FitsChecksum.lean`: 20, `AstroPhysics.lean`: 40).

### Clean Build Transcript (`lake clean && lake build`)

```text
Lean (version 4.36.0, x86_64-unknown-linux-gnu, commit 611c117bfaf7b7ae1cfcf246ce6fba31e8f509e9, Release)
Lake version 5.0.0-src+611c117 (Lean version 4.36.0)
✔ [2/12] Built AstrosatFormal.AstroPhysics (11s)
✔ [4/12] Built AstrosatFormal.GtiAlgebra (11s)
✔ [6/12] Built AstrosatFormal.FitsChecksum (12s)
✔ [8/12] Built AstrosatFormal.Telemetry (15s)
✔ [10/12] Built AstrosatFormal (15s)
Build completed successfully (12 jobs).
```

---

## 2. Package Structure

| File | Theorems | Scope |
| :--- | :---: | :--- |
| [`formal/lean-toolchain`](file:///d:/ISRO/formal/lean-toolchain) | — | Pinned to `leanprover/lean4:v4.36.0` |
| [`formal/lakefile.toml`](file:///d:/ISRO/formal/lakefile.toml) | — | Pure Lean 4 library target `AstrosatFormal` |
| [`formal/AstrosatFormal.lean`](file:///d:/ISRO/formal/AstrosatFormal.lean) | — | Root module importing all 4 proof suites |
| [`formal/AstrosatFormal/Telemetry.lean`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean) | 18 | CZTI 3-word/2-word unpacking, Q0–Q3 pixel & full-plane bijectivity, LAXPC 5-byte EA decoding, Xe K-escape clamping & dominant-anode bounds |
| [`formal/AstrosatFormal/GtiAlgebra.lean`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean) | 48 | GTI normalization, exposure non-negativity, intersection/subtraction soundness, exact partition identity, monotonicity, `drop_short` & `merge_gaps` |
| [`formal/AstrosatFormal/FitsChecksum.lean`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean) | 20 | FITS 4.0 2880-byte block alignment, 32-bit 1's-complement monoid laws, Seaman/Pence `CHECKSUM` complement & HDU verification soundness |
| [`formal/AstrosatFormal/AstroPhysics.lean`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean) | 40 | Exact `Q` rational field, LAXPC non-paralyzable deadtime (`tau = 43e-6`), CZTI mask-weight exact background cancellation, UTC↔TAI↔TT & leap-second continuity, CZTI 20-bit rollover |
| **Total** | **126** | **0 errors, 0 warnings, 0 `sorry`, 0 custom `axiom`** |

---

## 3. Theorem-to-Rust Traceability Matrix

### 3.1 `AstrosatFormal/Telemetry.lean` (18 Theorems)

| Lean 4 Theorem | Lean File:Line | Rust Implementation | Property Verified |
| :--- | :--- | :--- | :--- |
| `decodeCztiNormal_bounds` | [`Telemetry.lean:47`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L47) | [`science2event.rs:479-485`](file:///d:/ISRO/crates/astrosat-czti/src/science2event.rs#L479-L485) | For any `(w0, w1, w2) : BitVec 16`, `cztntick < 65536`, `pha < 1024`, `detid < 16`, `pixid < 256`, `veto < 128`, `alpha < 2` |
| `decodeCztiNormal_encode_bv` / `decodeCztiNormal_encode` | [`Telemetry.lean:69`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L69), [`L88`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L88) | [`science2event.rs:479-485`](file:///d:/ISRO/crates/astrosat-czti/src/science2event.rs#L479-L485) | Exact pack→unpack round-trip identity on both `BitVec 16` and `Nat` fields for 3-word normal events |
| `encodeCztiNormal_decode_bv` | [`Telemetry.lean:122`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L122) | [`science2event.rs:479-485`](file:///d:/ISRO/crates/astrosat-czti/src/science2event.rs#L479-L485) | Re-packing a decoded 3-word event recovers `w0`, `w2` identically and `w1` with its 2 low ADC bits |
| `decodeCztiReduced_bounds` | [`Telemetry.lean:169`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L169) | [`science2event.rs:614-622`](file:///d:/ISRO/crates/astrosat-czti/src/science2event.rs#L614-L622) | For any `(w0, w1) : BitVec 16`, `cztntick < 512`, `pha < 1024`, `pha % 2 = 0`, `detid < 16`, `pixid < 256`, `veto < 2 < 128`, `alpha < 2` |
| `decodeCztiReduced_encode_bv` / `encodeCztiReduced_decode_bv` / `decodeCztiReduced_encode` | [`Telemetry.lean:201`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L201), [`L215`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L215), [`L227`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L227) | [`science2event.rs:614-622`](file:///d:/ISRO/crates/astrosat-czti/src/science2event.rs#L614-L622) | Lossless 32-bit pack↔unpack bijection for 2-word reduced events (all 32 bits accounted for) |
| `detidPixidToDetxy_bounds` / `detxyToDetidPixid_bounds` | [`Telemetry.lean:285`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L285), [`L293`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L293) | [`czti/types.rs:97-110`](file:///d:/ISRO/crates/astrosat-czti/src/types.rs#L97-L110), [`L134-145`](file:///d:/ISRO/crates/astrosat-czti/src/types.rs#L134-L145) | Exact domain/codomain closure `[0,15] × [0,255] ↔ [0,63] × [0,63]` for all quadrants |
| `detxy_left_inv` / `detxy_right_inv` / `detidPixid_bijective` | [`Telemetry.lean:301`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L301), [`L313`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L313), [`L325`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L325) | [`czti/types.rs:97-110`](file:///d:/ISRO/crates/astrosat-czti/src/types.rs#L97-L110), [`L134-145`](file:///d:/ISRO/crates/astrosat-czti/src/types.rs#L134-L145) | Exact **bijectivity** on `Fin 16 × Fin 256 ↔ Fin 64 × Fin 64` for all 4 quadrants `Q0..Q3`, including Q0/Q3 `63 - dety` and Q1/Q2 `63 - detx` reflections |
| `detxyToLocxy_bijective` | [`Telemetry.lean:358`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L358) | [`czti/types.rs:120-130`](file:///d:/ISRO/crates/astrosat-czti/src/types.rs#L120-L130) | Full-plane 16,384-pixel bijection `Fin 4 × Fin 64 × Fin 64 ↔ Fin 128 × Fin 128` |
| `unpackEaWord_bounds` | [`Telemetry.lean:437`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L437) | [`ea.rs:51-77`](file:///d:/ISRO/crates/astrosat-laxpc/src/ea.rs#L51-L77) | Exact 5-byte EA bit-field bounds: `fineTicks < 256`, `ia1, ia2 < 16`, `ic1, ic2 < 1024`, `k1, k2 < 2` |
| `ea_double_energy_clamped` | [`Telemetry.lean:477`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L477) | [`ea.rs:82-98`](file:///d:/ISRO/crates/astrosat-laxpc/src/ea.rs#L82-L98) | `ic1 + ic2 ≤ 2046` (no `u16` overflow), `pha = min(ic1 + ic2, 1023) < 1024`, `anode < 8`, `anode ∈ {ia1, ia2}` |
| `ea_dominant_anode_total_bounds` | [`Telemetry.lean:507`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L507) | [`ea.rs:85-95`](file:///d:/ISRO/crates/astrosat-laxpc/src/ea.rs#L85-L95) | Totality & `1 ≤ anode ≤ 7` on all valid main-anode Xe K-escape double events (`1 ≤ ia1, ia2 ≤ 7`, `ic2 > 0`, `k1 = 1 ∨ k2 = 1`) |
| `accepted_double_event_bounds` | [`Telemetry.lean:536`](file:///d:/ISRO/formal/AstrosatFormal/Telemetry.lean#L536) | [`ea.rs:82-98`](file:///d:/ISRO/crates/astrosat-laxpc/src/ea.rs#L82-L98), [`L472`](file:///d:/ISRO/crates/astrosat-laxpc/src/ea.rs#L472) | Combined with `decode_ea_orbits` filter, every accepted double event satisfies `1 ≤ pha ≤ 999 ∧ 1 ≤ anode ≤ 7` |

### 3.2 `AstrosatFormal/GtiAlgebra.lean` (48 Theorems)

| Lean 4 Theorem | Lean File:Line | Rust Implementation | Property Verified |
| :--- | :--- | :--- | :--- |
| `duration_nonneg`, `exposure_nonneg`, `exposure_append` | [`GtiAlgebra.lean:49-65`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L49-L65) | [`gti.rs:24-26`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L24-L26), [`L149-151`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L149-L151) | Exposure non-negativity (`exposure L ≥ 0`) and additivity over concatenation |
| `normalize_normalized`, `memHalfOpen_normalize` (+ 14 helper lemmas) | [`GtiAlgebra.lean:248`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L248), [`L330`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L330) | [`gti.rs:42-62`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L42-L62) | `GtiSet::new` always produces strictly valid (`start < stop`), strictly sorted, strictly disjoint (`iv₁.stop < iv₂.start`) intervals preserving exact half-open membership |
| `intersect_normalized`, `memHalfOpen_intersect` (+ `intersect_normalizedAfter`) | [`GtiAlgebra.lean:362`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L362), [`L393`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L393), [`L404`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L404) | [`gti.rs:183-201`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L183-L201) | Two-pointer `intersect` preserves `Normalized` and satisfies exact set soundness `t ∈ A ∩ B ↔ t ∈ A ∧ t ∈ B` |
| `subtract_normalized`, `memHalfOpen_subtract` (+ `subtractOne_bounds`, `memHalfOpen_subtractOne`, `subtract_normalizedAfter`) | [`GtiAlgebra.lean:527`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L527), [`L597`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L597), [`L719`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L719) | [`gti.rs:212-245`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L212-L245) | `subtract` preserves `Normalized` and satisfies exact set soundness `t ∈ A \ B ↔ t ∈ A ∧ t ∉ B` |
| `exposure_intersect_add_subtract` (+ `exposure_partition_one`, `intersect_eq_flatMap_intersectOne`) | [`GtiAlgebra.lean:752`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L752), [`L827`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L827), [`L921`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L921) | [`gti.rs:149-151`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L149-L151), [`L183-245`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L183-L245) | **Exact partition identity**: `exposure(A ∩ B) + exposure(A \ B) = exposure(A)` for all normalized `A, B` |
| `exposure_intersect_le_left`, `exposure_subtract_le` | [`GtiAlgebra.lean:936`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L936), [`L944`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L944) | [`gti.rs:183-245`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L183-L245) | Monotonicity: `exposure(A ∩ B) ≤ exposure(A)` and `exposure(A \ B) ≤ exposure(A)` |
| `dropShort_normalized`, `dropShort_forall_ge`, `exposure_dropShort_le` | [`GtiAlgebra.lean:970`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L970), [`L979`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L979), [`L987`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L987) | [`gti.rs:172-180`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L172-L180) | `drop_short` preserves `Normalized`, retains only `duration ≥ min_exposure`, and never increases exposure |
| `mergeGaps_normalized`, `mergeGaps_minGap`, `exposure_le_mergeGaps` | [`GtiAlgebra.lean:1055`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L1055), [`L1070`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L1070), [`L1101`](file:///d:/ISRO/formal/AstrosatFormal/GtiAlgebra.lean#L1101) | [`gti.rs:154-169`](file:///d:/ISRO/crates/astrosat-core/src/gti.rs#L154-L169) | `merge_gaps` preserves `Normalized`, enforces gap `≥ max_gap` for `max_gap > 0`, and never decreases exposure |

### 3.3 `AstrosatFormal/FitsChecksum.lean` (20 Theorems)

| Lean 4 Theorem | Lean File:Line | Rust Implementation | Property Verified |
| :--- | :--- | :--- | :--- |
| `paddedLen_mod_zero`, `le_paddedLen`, `paddedLen_lt_add_block`, `paddedLen_zero`, `paddedLen_idempotent`, `paddedLen_minimal` | [`FitsChecksum.lean:30-62`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L30-L62) | [`io.rs:127-138`](file:///d:/ISRO/crates/astrosat-fits/src/io.rs#L127-L138), [`card.rs:503-504`](file:///d:/ISRO/crates/astrosat-fits/src/card.rs#L503-L504) | FITS 4.0 2880-byte block alignment: `paddedLen n % 2880 = 0`, `n ≤ paddedLen n < n + 2880`, idempotence & minimality |
| `writerPaddedLen_eq_paddedLen` | [`FitsChecksum.lean:47`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L47) | [`io.rs:127-138`](file:///d:/ISRO/crates/astrosat-fits/src/io.rs#L127-L138), [`L428-445`](file:///d:/ISRO/crates/astrosat-fits/src/io.rs#L428-L445) | Exact equivalence between reader `checked_next_multiple_of(2880)` / `div_ceil * 2880` and writer remainder padding |
| `add1c_single_fold_terminates`, `add1c_toNat` | [`FitsChecksum.lean:75`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L75), [`L83`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L83) | [`checksum.rs:121-127`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L121-L127) | Carry-fold loop in `add_ones_complement` always terminates in ≤ 1 iteration; exact `Nat` characterization |
| `add1c_comm`, `add1c_assoc`, `add1c_zero_left`, `add1c_zero_right`, `add1c_allOnes_right`, `add1c_mod_mersenne` | [`FitsChecksum.lean:113-161`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L113-L161) | [`checksum.rs:121-127`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L121-L127) | Exact commutative monoid laws (including strict associativity over all `BitVec 32`) and congruence mod `2^32 - 1` |
| `add1c_complement`, `add1c_complement_mod_zero` | [`FitsChecksum.lean:169`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L169), [`L175`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L175) | [`checksum.rs:147-148`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L147-L148), [`L173-174`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L173-L174) | Seaman/Pence complement theorem: `add1c x (~~~x) = 0xFFFFFFFF` (`≡ 0 mod (2^32 - 1)`) for every 32-bit word |
| `hdu_checksum_verification_sound` | [`FitsChecksum.lean:184`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L184) | [`checksum.rs:132-180`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L132-L180) | End-to-end HDU `compute_hdu_checksums` → `verify_hdu_checksum` soundness: `add1c (add1c hdrSum (~~~(add1c hdrSum dsum))) dsum = 0xFFFFFFFF` |
| `fold1c_append_zeros`, `fold1c_append` | [`FitsChecksum.lean:202`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L202), [`L213`](file:///d:/ISRO/formal/AstrosatFormal/FitsChecksum.lean#L213) | [`checksum.rs:18-45`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L18-L45) | Trailing FITS 2880-byte zero padding invariance and list-concatenation homomorphism of `accumulate_checksum` |

### 3.4 `AstrosatFormal/AstroPhysics.lean` (40 Theorems)

| Lean 4 Theorem | Lean File:Line | Rust Implementation | Property Verified |
| :--- | :--- | :--- | :--- |
| `Q.*` algebraic library (19 lemmas) | [`AstroPhysics.lean:77-285`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L77-L285) | — | Exact rational field/order laws over `Int` (`den > 0`), zero external dependencies |
| `deadtimeDenom_pos`, `deadtimeDenom_gt_eleven_thousandths`, `deadtimeDenom_num_pos` | [`AstroPhysics.lean:306`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L306), [`L320`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L320), [`L330`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L330) | [`laxpc/types.rs:305,308`](file:///d:/ISRO/crates/astrosat-laxpc/src/types.rs#L305-L308), [`products.rs:25-45`](file:///d:/ISRO/crates/astrosat-laxpc/src/products.rs#L25-L45) | For `tau = 43 / 1_000_000` and `0 ≤ R_all < 23000`, `1 - tau * R_all > 11 / 1000 > 0` (strictly positive denominator; no singularity/NaN/negative rate) |
| `deadtimeCorrectRate_ge_obs` | [`AstroPhysics.lean:349`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L349) | [`products.rs:25-33`](file:///d:/ISRO/crates/astrosat-laxpc/src/products.rs#L25-L33) | `R_obs ≥ 0 ∧ 0 ≤ R_all < 23000 → R_corr ≥ R_obs ≥ 0` |
| `deadtimeCorrectRate_strictMono_obs`, `deadtimeCorrectRate_strictMono_all`, `deadtimeFactor_ge_one_and_strictMono` | [`AstroPhysics.lean:374`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L374), [`L396`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L396), [`L436`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L436) | [`products.rs:25-45`](file:///d:/ISRO/crates/astrosat-laxpc/src/products.rs#L25-L45) | Strict monotonicity of `R_corr` in both `R_obs` and `R_all` (for `R_obs > 0`), and `f_cd ≥ 1` strictly increasing on `[0, 23000)` |
| `czti_mask_weight_background_cancellation` (+ `sumRenormWBT_decomp`, `sub_add_sub_equiv`) | [`AstroPhysics.lean:474`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L474), [`L501`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L501), [`L535`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L535) | [`bindata.rs:519-580`](file:///d:/ISRO/crates/astrosat-czti/src/bindata.rs#L519-L580) | For any finite good-pixel list with `sum_BT > 0` and `denom ≠ 0`, `w_i = (w'_i - D) / denom` with `D = (∑ w'_i B_i T_i)/(∑ B_i T_i)` satisfies `∑ w_i B_i T_i = 0` identically (both `≈ 0` and exact numerator `= 0`) |
| `ttToAstrosatUtc_astrosatUtcToTt`, `tai_utc_invertible` (+ `taiMinusUtc_congr`, `tt_step_simplify`, `ofInt_le_add_ofInt`) | [`AstroPhysics.lean:600`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L600), [`L652`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L652), [`L729`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L729) | [`time.rs:7-106`](file:///d:/ISRO/crates/astrosat-core/src/time.rs#L7-L106) | Global exact `UTC ↔ TAI ↔ TT` invertibility verifying the exact 3-step fixed-point loop `tt_sec_to_astrosat_utc` across all 4 AstroSat-era epochs (including the subtle 1-second window `[T_2012, T_2012 + 1)` requiring iteration 2 to converge!) |
| `leap_second_2012_07_01_step`, `leap_second_2015_07_01_step`, `leap_second_2017_01_01_step` (+ `leap_step_tai_tt`) | [`AstroPhysics.lean:734`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L734), [`L782`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L782), [`L802`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L802), [`L823`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L823) | [`time.rs:46-68,92-95`](file:///d:/ISRO/crates/astrosat-core/src/time.rs#L46-L95) | Exact 1-second leap-second step (`Δ(TAI - UTC) = 1`) and exact continuity `TAI(T_b) - TAI(T_b - ε) = TT(T_b) - TT(T_b - ε) = 1 + ε` across 2012-07-01, 2015-07-01, and 2017-01-01 |
| `czti_rollover_mod_invariant`, `czti_unwrapStep_exact` | [`AstroPhysics.lean:861`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L861), [`L866`](file:///d:/ISRO/formal/AstrosatFormal/AstroPhysics.lean#L866) | [`tct.rs:11,182-193,257-293`](file:///d:/ISRO/crates/astrosat-core/src/tct.rs#L182-L193) | Exact `mod 1048576` invariance under wrap offsets and exact elapsed-time recovery `unwrappedNext - unwrappedPrev = dt` + `unwrappedNext % 1048576 = nextRaw` across 20-bit counter rollovers |

---

## 4. Notable Findings & Formal Insights During Verification

No soundness bugs were found in the verified mathematical cores of `science2event.rs`, `czti/types.rs`, `ea.rs`, `gti.rs`, `checksum.rs`, `io.rs`, `products.rs`, `bindata.rs`, `time.rs`, or `tct.rs`. Formalization did surface four subtle structural properties worth documenting:

1. **Exact 3-iteration convergence of `tt_sec_to_astrosat_utc` ([`time.rs:98-106`](file:///d:/ISRO/crates/astrosat-core/src/time.rs#L98-L106))**: Because the initial guess subtracts `36.0`, for timestamps in the 1-second window `[2012-07-01T00:00:00, 2012-07-01T00:00:01)` where `TAI - UTC = 35`, `t0 = t - 1` drops into the pre-2012 epoch (`34`), iteration 1 overshoots to `t1 = t + 1` (landing back in `35`), iteration 2 hits `t2 = t` (`35`), and iteration 3 confirms the fixed point `t`. Lean proved that 3 iterations are both necessary and sufficient for exact invertibility across the entire AstroSat era (`2009-01-01` through `2026+`).
2. **CZTI `unwrapStep` requires `wrapOffset % 1048576 = 0` ([`tct.rs:174,190,205`](file:///d:/ISRO/crates/astrosat-core/src/tct.rs#L174-L205))**: Lean explicitly verified that `wrapOffset % 1048576 = 0` is an inductive loop invariant maintained by `wrap_offset = 0.0` (init/reset) and `wrap_offset += 1_048_576.0` (rollover), ensuring `unwrapped % 1048576 = raw` at every step.
3. **Strict associativity of 32-bit 1's-complement `add1c` ([`checksum.rs:121-127`](file:///d:/ISRO/crates/astrosat-fits/src/checksum.rs#L121-L127))**: Despite having two zero representations (`0x00000000` and `0xFFFFFFFF`), `add1c` is strictly associative on all of `BitVec 32` (not merely modulo `2^32 - 1`), and the carry-fold loop terminates in at most 1 iteration for any pair of 32-bit inputs.
4. **LAXPC double-event `anode == 0` guard ([`ea.rs:85-95,472`](file:///d:/ISRO/crates/astrosat-laxpc/src/ea.rs#L85-L95))**: `RawEaWord::resolve` guarantees `anode < 8` on raw 5-byte words; the strict `1 ≤ anode ≤ 7` invariant holds for all valid main-anode telemetry (`1 ≤ ia1, ia2 ≤ 7`) and unconditionally after the orbit-level validity filter at `ea.rs:472`.
