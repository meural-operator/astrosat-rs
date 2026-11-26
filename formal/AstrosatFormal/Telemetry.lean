import Std.Tactic.BVDecide

/-!
# AstrosatFormal.Telemetry

Formal verification of AstroSat CZTI and LAXPC telemetry bit-packing, unpacking,
coordinate mappings, and event resolution invariants, mirroring:
- `crates/astrosat-czti/src/science2event.rs` (3-word normal & 2-word reduced event records)
- `crates/astrosat-czti/src/types.rs` (`detid_pixid_to_detxy`, `detxy_to_detid_pixid`, `detxy_to_locxy`)
- `crates/astrosat-laxpc/src/ea.rs` (`RawEaWord::unpack`, `RawEaWord::resolve`, `decode_ea_orbits`)
-/

namespace AstrosatFormal.Telemetry

/-! ## 1. CZTI 3-Word Normal Event Unpacking & Packing (`science2event.rs:479-485`) -/

/-- Decoded fields of a CZTI 3-word normal event (`eventsize == 3`, `bunch_sync != 7`). -/
structure CztiNormalEvent where
  cztntick : Nat
  pha      : Nat
  detid    : Nat
  pixid    : Nat
  veto     : Nat
  alpha    : Nat
  deriving DecidableEq, Repr

/-- Exact mirror of `crates/astrosat-czti/src/science2event.rs:479-485`. -/
def decodeCztiNormal (w0 w1 w2 : BitVec 16) : CztiNormalEvent :=
  { cztntick := w0.toNat
    pha      := (((w1 >>> 4) &&& 0x0FFF#16) >>> 2).toNat
    detid    := (w1 &&& 0x000F#16).toNat
    pixid    := ((w2 >>> 8) &&& 0x00FF#16).toNat
    veto     := ((w2 >>> 1) &&& 0x007F#16).toNat
    alpha    := (w2 &&& 0x0001#16).toNat }

/-- Bit-level packer for CZTI 3-word normal event words (with optional 2-bit low ADC bits `low2 < 4`). -/
def encodeCztiNormal (cztntick pha low2 detid pixid veto alpha : Nat) :
    BitVec 16 × BitVec 16 × BitVec 16 :=
  let w0 : BitVec 16 := BitVec.ofNat 16 cztntick
  let w1 : BitVec 16 :=
    (BitVec.ofNat 16 pha <<< 6) ||| (BitVec.ofNat 16 low2 <<< 4) ||| BitVec.ofNat 16 detid
  let w2 : BitVec 16 :=
    (BitVec.ofNat 16 pixid <<< 8) ||| (BitVec.ofNat 16 veto <<< 1) ||| BitVec.ofNat 16 alpha
  (w0, w1, w2)

/-- Range bounds for ANY 3-word CZTI normal event record (`science2event.rs:479-485`). -/
theorem decodeCztiNormal_bounds (w0 w1 w2 : BitVec 16) :
    let e := decodeCztiNormal w0 w1 w2
    e.cztntick < 65536 ∧
    e.pha < 1024 ∧
    e.detid < 16 ∧
    e.pixid < 256 ∧
    e.veto < 128 ∧
    e.alpha < 2 := by
  dsimp [decodeCztiNormal]
  refine ⟨w0.isLt, ?_, ?_, ?_, ?_, ?_⟩
  · have h : (((w1 >>> 4) &&& 0x0FFF#16) >>> 2) < 1024#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : (w1 &&& 0x000F#16) < 16#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((w2 >>> 8) &&& 0x00FF#16) < 256#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((w2 >>> 1) &&& 0x007F#16) < 128#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : (w2 &&& 0x0001#16) < 2#16 := by bv_decide
    exact BitVec.lt_def.mp h

/-- Exact bit-level pack/unpack round-trip for CZTI 3-word normal events. -/
theorem decodeCztiNormal_encode_bv
    (tick pha low2 det pix veto alpha : BitVec 16)
    (hpha : pha < 1024#16)
    (hlow : low2 < 4#16)
    (hdet : det < 16#16)
    (hpix : pix < 256#16)
    (hveto : veto < 128#16)
    (halpha : alpha < 2#16) :
    let w1 := (pha <<< 6) ||| (low2 <<< 4) ||| det
    let w2 := (pix <<< 8) ||| (veto <<< 1) ||| alpha
    (((w1 >>> 4) &&& 0x0FFF#16) >>> 2) = pha ∧
    (w1 &&& 0x000F#16) = det ∧
    ((w2 >>> 8) &&& 0x00FF#16) = pix ∧
    ((w2 >>> 1) &&& 0x007F#16) = veto ∧
    (w2 &&& 0x0001#16) = alpha ∧
    tick = tick := by
  bv_decide

/-- Pack/unpack exactness on natural-number fields for CZTI 3-word normal events. -/
theorem decodeCztiNormal_encode
    (cztntick pha low2 detid pixid veto alpha : Nat)
    (htick : cztntick < 65536)
    (hpha : pha < 1024)
    (hlow : low2 < 4)
    (hdet : detid < 16)
    (hpix : pixid < 256)
    (hveto : veto < 128)
    (halpha : alpha < 2) :
    let (w0, w1, w2) := encodeCztiNormal cztntick pha low2 detid pixid veto alpha
    decodeCztiNormal w0 w1 w2 =
      { cztntick := cztntick, pha := pha, detid := detid,
        pixid := pixid, veto := veto, alpha := alpha } := by
  dsimp [encodeCztiNormal, decodeCztiNormal]
  have hbv := decodeCztiNormal_encode_bv
    (BitVec.ofNat 16 cztntick) (BitVec.ofNat 16 pha) (BitVec.ofNat 16 low2)
    (BitVec.ofNat 16 detid) (BitVec.ofNat 16 pixid) (BitVec.ofNat 16 veto)
    (BitVec.ofNat 16 alpha)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : pha < 65536)]; exact hpha)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : low2 < 65536)]; exact hlow)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : detid < 65536)]; exact hdet)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : pixid < 65536)]; exact hpix)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : veto < 65536)]; exact hveto)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : alpha < 65536)]; exact halpha)
  rcases hbv with ⟨h1, h2, h3, h4, h5, _⟩
  congr 1
  · simp [Nat.mod_eq_of_lt htick]
  · exact (congrArg BitVec.toNat h1).trans (by simp [Nat.mod_eq_of_lt (by omega : pha < 65536)])
  · exact (congrArg BitVec.toNat h2).trans (by simp [Nat.mod_eq_of_lt (by omega : detid < 65536)])
  · exact (congrArg BitVec.toNat h3).trans (by simp [Nat.mod_eq_of_lt (by omega : pixid < 65536)])
  · exact (congrArg BitVec.toNat h4).trans (by simp [Nat.mod_eq_of_lt (by omega : veto < 65536)])
  · exact (congrArg BitVec.toNat h5).trans (by simp [Nat.mod_eq_of_lt (by omega : alpha < 65536)])

/-- Re-packing a decoded 3-word event recovers `w0`, `w2` identically and `w1` up to the 2 dropped low ADC bits. -/
theorem encodeCztiNormal_decode_bv (w0 w1 w2 : BitVec 16) :
    let pha := ((w1 >>> 4) &&& 0x0FFF#16) >>> 2
    let low2 := (w1 >>> 4) &&& 0x0003#16
    let det := w1 &&& 0x000F#16
    let pix := (w2 >>> 8) &&& 0x00FF#16
    let veto := (w2 >>> 1) &&& 0x007F#16
    let alpha := w2 &&& 0x0001#16
    ((pha <<< 6) ||| (low2 <<< 4) ||| det) = w1 ∧
    ((pix <<< 8) ||| (veto <<< 1) ||| alpha) = w2 ∧
    w0 = w0 := by
  bv_decide

/-! ## 2. CZTI 2-Word Reduced Event Unpacking & Exact Bijection (`science2event.rs:614-622`) -/

/-- Decoded fields of a CZTI 2-word reduced event (`eventsize == 2`, `modeID ∈ {2,3,6,7}`). -/
structure CztiReducedEvent where
  cztntick : Nat
  pha      : Nat
  detid    : Nat
  pixid    : Nat
  veto     : Nat
  alpha    : Nat
  deriving DecidableEq, Repr

/-- Exact mirror of `crates/astrosat-czti/src/science2event.rs:614-622`. -/
def decodeCztiReduced (w0 w1 : BitVec 16) : CztiReducedEvent :=
  { cztntick := ((w0 &&& 0xFF80#16) >>> 7).toNat
    pha      := ((((w0 &&& 0x007F#16) <<< 2) ||| ((w1 >>> 14) &&& 0x0003#16)) <<< 1).toNat
    detid    := (w1 &&& 0x000F#16).toNat
    pixid    := ((w1 >>> 6) &&& 0x00FF#16).toNat
    veto     := ((w1 >>> 5) &&& 0x0001#16).toNat
    alpha    := ((w1 >>> 4) &&& 0x0001#16).toNat }

/-- Bit-level packer for CZTI 2-word reduced events (taking 9-bit reduced PHA `pha9 < 512`). -/
def encodeCztiReduced (cztntick pha9 detid pixid veto alpha : Nat) :
    BitVec 16 × BitVec 16 :=
  let t : BitVec 16 := BitVec.ofNat 16 cztntick
  let p : BitVec 16 := BitVec.ofNat 16 pha9
  let d : BitVec 16 := BitVec.ofNat 16 detid
  let x : BitVec 16 := BitVec.ofNat 16 pixid
  let v : BitVec 16 := BitVec.ofNat 16 veto
  let a : BitVec 16 := BitVec.ofNat 16 alpha
  let w0 := (t <<< 7) ||| (p >>> 2)
  let w1 := ((p &&& 0x0003#16) <<< 14) ||| (x <<< 6) ||| (v <<< 5) ||| (a <<< 4) ||| d
  (w0, w1)

/-- Range bounds (and evenness of reduced PHA) for ANY 2-word CZTI reduced event (`science2event.rs:614-622`). -/
theorem decodeCztiReduced_bounds (w0 w1 : BitVec 16) :
    let e := decodeCztiReduced w0 w1
    e.cztntick < 512 ∧
    e.pha < 1024 ∧
    e.pha % 2 = 0 ∧
    e.detid < 16 ∧
    e.pixid < 256 ∧
    e.veto < 2 ∧
    e.veto < 128 ∧
    e.alpha < 2 := by
  dsimp [decodeCztiReduced]
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · have h : ((w0 &&& 0xFF80#16) >>> 7) < 512#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((((w0 &&& 0x007F#16) <<< 2) ||| ((w1 >>> 14) &&& 0x0003#16)) <<< 1) < 1024#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : (((((w0 &&& 0x007F#16) <<< 2) ||| ((w1 >>> 14) &&& 0x0003#16)) <<< 1) &&& 1#16) = 0#16 := by bv_decide
    have hnat := congrArg BitVec.toNat h
    simp only [BitVec.toNat_and, BitVec.toNat_ofNat, Nat.reducePow, Nat.reduceMod, Nat.and_one_is_mod] at hnat
    exact hnat
  · have h : (w1 &&& 0x000F#16) < 16#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((w1 >>> 6) &&& 0x00FF#16) < 256#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((w1 >>> 5) &&& 0x0001#16) < 2#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((w1 >>> 5) &&& 0x0001#16) < 128#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((w1 >>> 4) &&& 0x0001#16) < 2#16 := by bv_decide
    exact BitVec.lt_def.mp h

/-- Exact left & right inverse bit-vector proof for CZTI 2-word reduced events (all 32 bits lossless). -/
theorem decodeCztiReduced_encode_bv
    (t p d x v a : BitVec 16)
    (ht : t < 512#16) (hp : p < 512#16) (hd : d < 16#16)
    (hx : x < 256#16) (hv : v < 2#16) (ha : a < 2#16) :
    let w0 := (t <<< 7) ||| (p >>> 2)
    let w1 := ((p &&& 0x0003#16) <<< 14) ||| (x <<< 6) ||| (v <<< 5) ||| (a <<< 4) ||| d
    ((w0 &&& 0xFF80#16) >>> 7) = t ∧
    ((((w0 &&& 0x007F#16) <<< 2) ||| ((w1 >>> 14) &&& 0x0003#16)) <<< 1) = (p <<< 1) ∧
    (w1 &&& 0x000F#16) = d ∧
    ((w1 >>> 6) &&& 0x00FF#16) = x ∧
    ((w1 >>> 5) &&& 0x0001#16) = v ∧
    ((w1 >>> 4) &&& 0x0001#16) = a := by
  bv_decide

theorem encodeCztiReduced_decode_bv (w0 w1 : BitVec 16) :
    let t := (w0 &&& 0xFF80#16) >>> 7
    let p := ((w0 &&& 0x007F#16) <<< 2) ||| ((w1 >>> 14) &&& 0x0003#16)
    let d := w1 &&& 0x000F#16
    let x := (w1 >>> 6) &&& 0x00FF#16
    let v := (w1 >>> 5) &&& 0x0001#16
    let a := (w1 >>> 4) &&& 0x0001#16
    ((t <<< 7) ||| (p >>> 2)) = w0 ∧
    (((p &&& 0x0003#16) <<< 14) ||| (x <<< 6) ||| (v <<< 5) ||| (a <<< 4) ||| d) = w1 := by
  bv_decide

/-- Pack/unpack exactness on natural-number fields for CZTI 2-word reduced events. -/
theorem decodeCztiReduced_encode
    (cztntick pha9 detid pixid veto alpha : Nat)
    (ht : cztntick < 512) (hp : pha9 < 512) (hd : detid < 16)
    (hx : pixid < 256) (hv : veto < 2) (ha : alpha < 2) :
    let (w0, w1) := encodeCztiReduced cztntick pha9 detid pixid veto alpha
    decodeCztiReduced w0 w1 =
      { cztntick := cztntick, pha := 2 * pha9, detid := detid,
        pixid := pixid, veto := veto, alpha := alpha } := by
  dsimp [encodeCztiReduced, decodeCztiReduced]
  have hbv := decodeCztiReduced_encode_bv
    (BitVec.ofNat 16 cztntick) (BitVec.ofNat 16 pha9) (BitVec.ofNat 16 detid)
    (BitVec.ofNat 16 pixid) (BitVec.ofNat 16 veto) (BitVec.ofNat 16 alpha)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : cztntick < 65536)]; exact ht)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : pha9 < 65536)]; exact hp)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : detid < 65536)]; exact hd)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : pixid < 65536)]; exact hx)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : veto < 65536)]; exact hv)
    (by rw [BitVec.lt_def]; simp [Nat.mod_eq_of_lt (by omega : alpha < 65536)]; exact ha)
  rcases hbv with ⟨h1, h2, h3, h4, h5, h6⟩
  congr 1
  · exact (congrArg BitVec.toNat h1).trans (by simp [Nat.mod_eq_of_lt (by omega : cztntick < 65536)])
  · exact (congrArg BitVec.toNat h2).trans (by
      simp only [BitVec.toNat_shiftLeft, BitVec.toNat_ofNat, Nat.reducePow]
      omega)
  · exact (congrArg BitVec.toNat h3).trans (by simp [Nat.mod_eq_of_lt (by omega : detid < 65536)])
  · exact (congrArg BitVec.toNat h4).trans (by simp [Nat.mod_eq_of_lt (by omega : pixid < 65536)])
  · exact (congrArg BitVec.toNat h5).trans (by simp [Nat.mod_eq_of_lt (by omega : veto < 65536)])
  · exact (congrArg BitVec.toNat h6).trans (by simp [Nat.mod_eq_of_lt (by omega : alpha < 65536)])

/-! ## 3. Exact Bijectivity of CZTI `(detid, pixid) ↔ (detx, dety)` for ALL 4 Quadrants (`types.rs:97-145`) -/

/-- Exact mirror of `crates/astrosat-czti/src/types.rs:97-110` (`detid_pixid_to_detxy`). -/
def detidPixidToDetxy (q detid pixid : Nat) : Nat × Nat :=
  let dx := detid % 4
  let dy := detid / 4
  let px := pixid % 16
  let py := pixid / 16
  let detx := dx * 16 + px
  let dety := dy * 16 + py
  if q == 0 || q == 3 then
    (detx, 63 - dety)
  else
    (63 - detx, dety)

/-- Exact mirror of `crates/astrosat-czti/src/types.rs:134-145` (`detxy_to_detid_pixid`). -/
def detxyToDetidPixid (q detx dety : Nat) : Nat × Nat :=
  let (ux, uy) :=
    if q == 0 || q == 3 then
      (detx, 63 - dety)
    else
      (63 - detx, dety)
  let dx := ux / 16
  let dy := uy / 16
  let px := ux % 16
  let py := uy % 16
  (dy * 4 + dx, py * 16 + px)

/-- Range preservation: `[0, 15] × [0, 255] → [0, 63] × [0, 63]` for every quadrant `q`. -/
theorem detidPixidToDetxy_bounds (q detid pixid : Nat)
    (hd : detid < 16) (hp : pixid < 256) :
    let (detx, dety) := detidPixidToDetxy q detid pixid
    detx < 64 ∧ dety < 64 := by
  dsimp [detidPixidToDetxy]
  split <;> omega

/-- Inverse range preservation: `[0, 63] × [0, 63] → [0, 15] × [0, 255]` for every quadrant `q`. -/
theorem detxyToDetidPixid_bounds (q detx dety : Nat)
    (hx : detx < 64) (hy : dety < 64) :
    let (detid, pixid) := detxyToDetidPixid q detx dety
    detid < 16 ∧ pixid < 256 := by
  dsimp [detxyToDetidPixid]
  split <;> omega

/-- Left inverse: `detxyToDetidPixid q (detidPixidToDetxy q detid pixid) = (detid, pixid)` for all Q0..Q3. -/
theorem detxy_left_inv (q detid pixid : Nat)
    (hd : detid < 16) (hp : pixid < 256) :
    let (detx, dety) := detidPixidToDetxy q detid pixid
    detxyToDetidPixid q detx dety = (detid, pixid) := by
  dsimp [detidPixidToDetxy, detxyToDetidPixid]
  by_cases hq : (q == 0 || q == 3) = true
  · simp only [hq, ↓reduceIte]
    ext <;> omega
  · simp only [hq, Bool.false_eq_true, ↓reduceIte]
    ext <;> omega

/-- Right inverse: `detidPixidToDetxy q (detxyToDetidPixid q detx dety) = (detx, dety)` for all Q0..Q3. -/
theorem detxy_right_inv (q detx dety : Nat)
    (hx : detx < 64) (hy : dety < 64) :
    let (detid, pixid) := detxyToDetidPixid q detx dety
    detidPixidToDetxy q detid pixid = (detx, dety) := by
  dsimp [detidPixidToDetxy, detxyToDetidPixid]
  by_cases hq : (q == 0 || q == 3) = true
  · simp only [hq, ↓reduceIte]
    ext <;> omega
  · simp only [hq, Bool.false_eq_true, ↓reduceIte]
    ext <;> omega

/-- Exact bijection on `Fin 16 × Fin 256 ↔ Fin 64 × Fin 64` for every quadrant `q : Fin 4` (Q0, Q1, Q2, Q3). -/
theorem detidPixid_bijective (q : Fin 4) :
    (∀ d : Fin 16, ∀ p : Fin 256,
      let xy := detidPixidToDetxy q.val d.val p.val
      xy.1 < 64 ∧ xy.2 < 64 ∧ detxyToDetidPixid q.val xy.1 xy.2 = (d.val, p.val)) ∧
    (∀ x : Fin 64, ∀ y : Fin 64,
      let dp := detxyToDetidPixid q.val x.val y.val
      dp.1 < 16 ∧ dp.2 < 256 ∧ detidPixidToDetxy q.val dp.1 dp.2 = (x.val, y.val)) := by
  constructor
  · intro d p
    have hb := detidPixidToDetxy_bounds q.val d.val p.val d.isLt p.isLt
    have hl := detxy_left_inv q.val d.val p.val d.isLt p.isLt
    exact ⟨hb.1, hb.2, hl⟩
  · intro x y
    have hb := detxyToDetidPixid_bounds q.val x.val y.val x.isLt y.isLt
    have hr := detxy_right_inv q.val x.val y.val x.isLt y.isLt
    exact ⟨hb.1, hb.2, hr⟩

/-- Exact mirror of `crates/astrosat-czti/src/types.rs:120-130` (`detxy_to_locxy`). -/
def detxyToLocxy (q detx dety : Nat) : Nat × Nat :=
  match q with
  | 0 => (detx, 64 + dety)
  | 1 => (64 + detx, 64 + dety)
  | 2 => (64 + detx, dety)
  | _ => (detx, dety)

/-- Inverse full-plane map `(locx, locy) ∈ [0, 127]² → (q, detx, dety) ∈ [0,3] × [0,63]²`. -/
def locxyToQuadDetxy (locx locy : Nat) : Nat × Nat × Nat :=
  if locx < 64 then
    if locy < 64 then (3, locx, locy) else (0, locx, locy - 64)
  else
    if locy < 64 then (2, locx - 64, locy) else (1, locx - 64, locy - 64)

/-- Full-plane 16,384-pixel bijectivity across Q0..Q3 (`types.rs:120-130`). -/
theorem detxyToLocxy_bijective :
    (∀ q detx dety, q < 4 → detx < 64 → dety < 64 →
      let (lx, ly) := detxyToLocxy q detx dety
      lx < 128 ∧ ly < 128 ∧ locxyToQuadDetxy lx ly = (q, detx, dety)) ∧
    (∀ lx ly, lx < 128 → ly < 128 →
      let (q, dx, dy) := locxyToQuadDetxy lx ly
      q < 4 ∧ dx < 64 ∧ dy < 64 ∧ detxyToLocxy q dx dy = (lx, ly)) := by
  constructor
  · intro q detx dety hq hx hy
    match q with
    | 0 =>
      have h1 : ¬ (64 + dety < 64) := by omega
      refine ⟨by omega, by omega, ?_⟩
      simp [locxyToQuadDetxy, hx, h1]
    | 1 =>
      have h1 : ¬ (64 + detx < 64) := by omega
      have h2 : ¬ (64 + dety < 64) := by omega
      refine ⟨by omega, by omega, ?_⟩
      simp [locxyToQuadDetxy, h1, h2]
    | 2 =>
      have h1 : ¬ (64 + detx < 64) := by omega
      refine ⟨by omega, by omega, ?_⟩
      simp [locxyToQuadDetxy, h1, hy]
    | 3 =>
      refine ⟨by omega, by omega, ?_⟩
      simp [locxyToQuadDetxy, hx, hy]
    | n + 4 => omega
  · intro lx ly hlx hly
    by_cases h1 : lx < 64 <;> by_cases h2 : ly < 64
    · simp only [locxyToQuadDetxy, detxyToLocxy, h1, h2, ↓reduceIte]
      exact ⟨by omega, trivial, trivial, trivial⟩
    · simp only [locxyToQuadDetxy, detxyToLocxy, h1, h2, ↓reduceIte]
      refine ⟨by omega, trivial, by omega, ?_⟩
      ext <;> dsimp <;> omega
    · simp only [locxyToQuadDetxy, detxyToLocxy, h1, h2, ↓reduceIte]
      refine ⟨by omega, by omega, trivial, ?_⟩
      ext <;> dsimp <;> omega
    · simp only [locxyToQuadDetxy, detxyToLocxy, h1, h2, ↓reduceIte]
      refine ⟨by omega, by omega, by omega, ?_⟩
      ext <;> dsimp <;> omega

/-! ## 4. LAXPC 5-Byte EA Event Unpacking, Xe K-Escape Clamping & Anode Dominance (`ea.rs:51-98,472`) -/

/-- Event classification mirroring `crates/astrosat-laxpc/src/ea.rs:22-27`. -/
inductive LaxpcEventKind where
  | single
  | double
  deriving DecidableEq, Repr

/-- Unpacked 5-byte LAXPC EA event word mirroring `crates/astrosat-laxpc/src/ea.rs:31-46`. -/
structure RawEaWord where
  fineTicks : Nat
  ia1       : Nat
  ia2       : Nat
  ic1       : Nat
  ic2       : Nat
  k1        : Nat
  k2        : Nat
  deriving DecidableEq, Repr

/-- Exact mirror of `RawEaWord::unpack` (`crates/astrosat-laxpc/src/ea.rs:51-77`). -/
def unpackEaWord (b0 b1 b2 b3 b4 : BitVec 8) : RawEaWord :=
  let i5 := b0
  let i4 := b1
  let i3 := b2.zeroExtend 16
  let i2 := b3.zeroExtend 16
  let i1 := b4
  let fineTicks := (i4 - 2#8).toNat
  let ia1 := (i5 &&& 0x0F#8).toNat
  let ia2 := ((i5 >>> 4) &&& 0x0F#8).toNat
  let k1 := (i1 &&& 0x01#8).toNat
  let l1 := ((i1 >>> 2) &&& 0x03#8).zeroExtend 16
  let k2 := ((i1 >>> 4) &&& 0x01#8).toNat
  let l2 := ((i1 >>> 6) &&& 0x03#8).zeroExtend 16
  let ic1 := ((i3 <<< 2) ||| l1).toNat
  let ic2 := ((i2 <<< 2) ||| l2).toNat
  { fineTicks, ia1, ia2, ic1, ic2, k1, k2 }

/-- Exact bit-field bounds for ANY 5-byte LAXPC EA record (`ea.rs:51-77`). -/
theorem unpackEaWord_bounds (b0 b1 b2 b3 b4 : BitVec 8) :
    let w := unpackEaWord b0 b1 b2 b3 b4
    w.fineTicks < 256 ∧
    w.ia1 < 16 ∧
    w.ia2 < 16 ∧
    w.ic1 < 1024 ∧
    w.ic2 < 1024 ∧
    w.k1 < 2 ∧
    w.k2 < 2 := by
  dsimp [unpackEaWord]
  refine ⟨(b1 - 2#8).isLt, ?_, ?_, ?_, ?_, ?_, ?_⟩
  · have h : (b0 &&& 0x0F#8) < 16#8 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((b0 >>> 4) &&& 0x0F#8) < 16#8 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((b2.zeroExtend 16 <<< 2) ||| ((b4 >>> 2) &&& 0x03#8).zeroExtend 16) < 1024#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((b3.zeroExtend 16 <<< 2) ||| ((b4 >>> 6) &&& 0x03#8).zeroExtend 16) < 1024#16 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : (b4 &&& 0x01#8) < 2#8 := by bv_decide
    exact BitVec.lt_def.mp h
  · have h : ((b4 >>> 4) &&& 0x01#8) < 2#8 := by bv_decide
    exact BitVec.lt_def.mp h

/-- Exact mirror of `RawEaWord::resolve` (`crates/astrosat-laxpc/src/ea.rs:82-98`). -/
def resolveEaWord (w : RawEaWord) : Option (Nat × Nat × LaxpcEventKind) :=
  if w.ic2 == 0 then
    some (w.ic1, w.ia1, LaxpcEventKind.single)
  else if (w.k1 == 1 || w.k2 == 1) && w.ia1 < 8 && w.ia2 < 8 then
    let pha := min (w.ic1 + w.ic2) 1023
    let anode :=
      if w.k1 == 1 && (w.ic1 ≤ w.ic2 || (w.k2 == 0 && w.ic2 > 1)) then
        w.ia2
      else
        w.ia1
    some (pha, anode, LaxpcEventKind.double)
  else
    none

/-- Xe K-escape double-event energy clamping: `min(ic1 + ic2, 1023) < 1024` and no `u16` overflow. -/
theorem ea_double_energy_clamped (b0 b1 b2 b3 b4 : BitVec 8) (pha anode : Nat) :
    let w := unpackEaWord b0 b1 b2 b3 b4
    w.ic1 + w.ic2 ≤ 2046 ∧
    (resolveEaWord w = some (pha, anode, LaxpcEventKind.double) →
      pha = min (w.ic1 + w.ic2) 1023 ∧ pha < 1024 ∧ anode < 8 ∧ (anode = w.ia1 ∨ anode = w.ia2)) := by
  have hb := unpackEaWord_bounds b0 b1 b2 b3 b4
  dsimp only at hb ⊢
  constructor
  · omega
  · intro hres
    unfold resolveEaWord at hres
    split at hres
    · cases hres
    · split at hres
      · rename_i _ hcond
        injection hres with htuple
        have h1 : pha = min ((unpackEaWord b0 b1 b2 b3 b4).ic1 + (unpackEaWord b0 b1 b2 b3 b4).ic2) 1023 := by
          cases htuple; rfl
        have h2 : anode = (if (unpackEaWord b0 b1 b2 b3 b4).k1 == 1 &&
            ((unpackEaWord b0 b1 b2 b3 b4).ic1 ≤ (unpackEaWord b0 b1 b2 b3 b4).ic2 ||
              ((unpackEaWord b0 b1 b2 b3 b4).k2 == 0 && (unpackEaWord b0 b1 b2 b3 b4).ic2 > 1))
            then (unpackEaWord b0 b1 b2 b3 b4).ia2 else (unpackEaWord b0 b1 b2 b3 b4).ia1) := by
          cases htuple; rfl
        simp only [Bool.and_eq_true, decide_eq_true_eq] at hcond
        refine ⟨h1, by omega, ?_, ?_⟩
        · rw [h2]; split <;> omega
        · rw [h2]; split <;> simp
      · cases hres

/-- Totality and exact `1 ≤ anode ≤ 7` bounds for valid main-anode Xe K-escape double events (`ea.rs:85-95`). -/
theorem ea_dominant_anode_total_bounds (w : RawEaWord)
    (hic2 : w.ic2 > 0)
    (hk : w.k1 = 1 ∨ w.k2 = 1)
    (ha1 : 1 ≤ w.ia1 ∧ w.ia1 ≤ 7)
    (ha2 : 1 ≤ w.ia2 ∧ w.ia2 ≤ 7) :
    ∃ pha anode,
      resolveEaWord w = some (pha, anode, LaxpcEventKind.double) ∧
      pha < 1024 ∧
      1 ≤ anode ∧ anode ≤ 7 ∧
      (anode = w.ia1 ∨ anode = w.ia2) := by
  have hne : (w.ic2 == 0) = false := by
    simp only [beq_eq_false_iff_ne]; omega
  have hcond : ((w.k1 == 1 || w.k2 == 1) && decide (w.ia1 < 8) && decide (w.ia2 < 8)) = true := by
    simp only [Bool.and_eq_true, Bool.or_eq_true, beq_iff_eq, decide_eq_true_eq]
    exact ⟨⟨hk, by omega⟩, by omega⟩
  refine ⟨min (w.ic1 + w.ic2) 1023,
    if w.k1 == 1 && (w.ic1 ≤ w.ic2 || (w.k2 == 0 && w.ic2 > 1)) then w.ia2 else w.ia1,
    ?_, by omega, ?_, ?_, ?_⟩
  · unfold resolveEaWord
    simp [hne, hcond]
  · split <;> omega
  · split <;> omega
  · split <;> simp

/-- Exact orbit-level validity filter (`ea.rs:472`): every accepted double event satisfies
    `1 ≤ pha ≤ 999` and `1 ≤ anode ≤ 7`. -/
def isValidOrbitEvent (pha anode : Nat) : Bool :=
  !(pha == 0 || pha >= 1000 || anode == 0 || anode > 10)

theorem accepted_double_event_bounds (b0 b1 b2 b3 b4 : BitVec 8) (pha anode : Nat)
    (hres : resolveEaWord (unpackEaWord b0 b1 b2 b3 b4) = some (pha, anode, LaxpcEventKind.double))
    (hacc : isValidOrbitEvent pha anode = true) :
    1 ≤ pha ∧ pha ≤ 999 ∧ 1 ≤ anode ∧ anode ≤ 7 := by
  have h := (ea_double_energy_clamped b0 b1 b2 b3 b4 pha anode).2 hres
  unfold isValidOrbitEvent at hacc
  simp only [Bool.not_eq_true', Bool.or_eq_false_iff, beq_eq_false_iff_ne, decide_eq_false_iff_not] at hacc
  omega

end AstrosatFormal.Telemetry
