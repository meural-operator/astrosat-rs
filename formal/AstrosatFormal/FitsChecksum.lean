import Std.Tactic.BVDecide

/-!
# AstrosatFormal.FitsChecksum

Formal verification of FITS 4.0 2880-byte block alignment and 32-bit 1's-complement
checksum algebra (Seaman, Pence & Rots 2012), mirroring:
- `crates/astrosat-fits/src/io.rs` (`data_padded_len`, `FitsWriter::write_hdu`)
- `crates/astrosat-fits/src/card.rs` (`Header::to_fits_bytes`)
- `crates/astrosat-fits/src/checksum.rs` (`FITS_BLOCK_SIZE`, `add_ones_complement`,
  `accumulate_checksum`, `compute_hdu_checksums`, `verify_hdu_checksum`)
-/

namespace AstrosatFormal.FitsChecksum

/-! ## 1. FITS 4.0 2880-Byte Block Alignment (`io.rs:127-131,412-414`, `card.rs:503-504`) -/

/-- FITS standard logical block size in bytes (`checksum.rs:7`). -/
def FITS_BLOCK_SIZE : Nat := 2880

/-- Exact mirror of `n.div_ceil(FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE` (`io.rs:130`, `card.rs:504`). -/
def paddedLen (n : Nat) : Nat :=
  ((n + FITS_BLOCK_SIZE - 1) / FITS_BLOCK_SIZE) * FITS_BLOCK_SIZE

/-- Exact mirror of `FitsWriter::write_hdu` remainder-based padding (`io.rs:412-414`). -/
def writerPaddedLen (n : Nat) : Nat :=
  let rem := n % FITS_BLOCK_SIZE
  if rem == 0 then n else n + (FITS_BLOCK_SIZE - rem)

theorem paddedLen_mod_zero (n : Nat) : paddedLen n % 2880 = 0 := by
  unfold paddedLen FITS_BLOCK_SIZE
  exact Nat.mul_mod_left _ 2880

theorem le_paddedLen (n : Nat) : n ≤ paddedLen n := by
  unfold paddedLen FITS_BLOCK_SIZE
  omega

theorem paddedLen_lt_add_block (n : Nat) : paddedLen n < n + 2880 := by
  unfold paddedLen FITS_BLOCK_SIZE
  omega

theorem paddedLen_zero : paddedLen 0 = 0 := by
  decide

/-- The reader's `div_ceil * 2880` (`io.rs:130`) and writer's `rem` padding (`io.rs:412-414`)
    are identically equal for all payload sizes `n : Nat`. -/
theorem writerPaddedLen_eq_paddedLen (n : Nat) : writerPaddedLen n = paddedLen n := by
  unfold writerPaddedLen paddedLen FITS_BLOCK_SIZE
  by_cases h : n % 2880 = 0
  · simp [h]; omega
  · have hne : (n % 2880 == 0) = false := beq_false_of_ne h
    simp [hne]; omega

/-- `paddedLen` is idempotent: padding an already-padded length leaves it unchanged. -/
theorem paddedLen_idempotent (n : Nat) : paddedLen (paddedLen n) = paddedLen n := by
  unfold paddedLen FITS_BLOCK_SIZE
  omega

/-- `paddedLen n` is the minimal multiple of 2880 greater than or equal to `n`. -/
theorem paddedLen_minimal (n k : Nat) (hmult : k % 2880 = 0) (hle : n ≤ k) :
    paddedLen n ≤ k := by
  unfold paddedLen FITS_BLOCK_SIZE
  omega

/-! ## 2. 32-Bit 1's-Complement Addition (`checksum.rs:121-127`) -/

/-- Exact 1-fold 32-bit 1's-complement addition (`checksum.rs:121-127`). -/
def add1c (a b : BitVec 32) : BitVec 32 :=
  let s : BitVec 64 := a.zeroExtend 64 + b.zeroExtend 64
  let folded : BitVec 64 := (s &&& 0xFFFFFFFF#64) + (s >>> 32)
  folded.truncate 32

/-- Proof that the `while (s >> 32) != 0` loop in `checksum.rs:123-125` terminates
    after AT MOST ONE fold step for ALL `a, b : BitVec 32`. -/
theorem add1c_single_fold_terminates (a b : BitVec 32) :
    let s : BitVec 64 := a.zeroExtend 64 + b.zeroExtend 64
    let folded : BitVec 64 := (s &&& 0xFFFFFFFF#64) + (s >>> 32)
    (folded >>> 32) = 0#64 := by
  bv_decide

/-- Arithmetic characterization of `add1c` on `Nat`:
    if `a + b < 2^32`, `add1c a b = a + b`; otherwise `add1c a b = a + b - (2^32 - 1)`. -/
theorem add1c_toNat (a b : BitVec 32) :
    (add1c a b).toNat =
      if a.toNat + b.toNat < 4294967296 then
        a.toNat + b.toNat
      else
        a.toNat + b.toNat - 4294967295 := by
  unfold add1c
  have ha : a.toNat < 4294967296 := a.isLt
  have hb : b.toNat < 4294967296 := b.isLt
  have h1 : (a.zeroExtend 64 + b.zeroExtend 64).toNat = a.toNat + b.toNat := by
    rw [BitVec.toNat_add, BitVec.toNat_setWidth, BitVec.toNat_setWidth]
    have h64 : (2 ^ 64 : Nat) = 18446744073709551616 := by decide
    rw [h64, Nat.mod_eq_of_lt (by omega), Nat.mod_eq_of_lt (by omega), Nat.mod_eq_of_lt (by omega)]
  have hmask : ((a.zeroExtend 64 + b.zeroExtend 64) &&& 0xFFFFFFFF#64).toNat =
      (a.toNat + b.toNat) % 4294967296 := by
    rw [BitVec.toNat_and, h1]
    have hpow : (0xFFFFFFFF#64 : BitVec 64).toNat = 2 ^ 32 - 1 := by decide
    rw [hpow, Nat.and_two_pow_sub_one_eq_mod]
  have hshift : ((a.zeroExtend 64 + b.zeroExtend 64) >>> 32).toNat =
      (a.toNat + b.toNat) / 4294967296 := by
    rw [BitVec.toNat_ushiftRight, h1]
    exact Nat.shiftRight_eq_div_pow (a.toNat + b.toNat) 32
  dsimp only
  rw [BitVec.toNat_setWidth, BitVec.toNat_add, hmask, hshift]
  have h64 : (2 ^ 64 : Nat) = 18446744073709551616 := by decide
  have h32 : (2 ^ 32 : Nat) = 4294967296 := by decide
  rw [h64, h32]
  split <;> omega

/-- Commutativity of 32-bit 1's-complement addition (`add1c a b = add1c b a`). -/
theorem add1c_comm (a b : BitVec 32) : add1c a b = add1c b a := by
  apply BitVec.eq_of_toNat_eq
  rw [add1c_toNat, add1c_toNat]
  split <;> split <;> omega

/-- Exact left identity `0` (`add1c 0 x = x`). -/
theorem add1c_zero_left (x : BitVec 32) : add1c 0#32 x = x := by
  apply BitVec.eq_of_toNat_eq
  rw [add1c_toNat]
  have hx := x.isLt
  simp [hx]

/-- Exact right identity `0` (`add1c x 0 = x`). -/
theorem add1c_zero_right (x : BitVec 32) : add1c x 0#32 = x := by
  rw [add1c_comm, add1c_zero_left]

/-- `-0` (`0xFFFFFFFF`) is also an identity for every non-zero 32-bit word. -/
theorem add1c_allOnes_right (x : BitVec 32) (hx : x ≠ 0#32) :
    add1c x 0xFFFFFFFF#32 = x := by
  apply BitVec.eq_of_toNat_eq
  rw [add1c_toNat]
  have hlt := x.isLt
  have hnz : x.toNat ≠ 0 := by
    intro h0
    apply hx
    apply BitVec.eq_of_toNat_eq
    simpa using h0
  have hff : (0xFFFFFFFF#32 : BitVec 32).toNat = 4294967295 := by decide
  rw [hff]
  split <;> omega

/-- Strict associativity of 32-bit 1's-complement addition over ALL `a, b, c : BitVec 32`! -/
theorem add1c_assoc (a b c : BitVec 32) :
    add1c (add1c a b) c = add1c a (add1c b c) := by
  apply BitVec.eq_of_toNat_eq
  have ha := a.isLt
  have hb := b.isLt
  have hc := c.isLt
  rw [add1c_toNat, add1c_toNat, add1c_toNat, add1c_toNat]
  split <;> split <;> split <;> split <;> omega

/-- Exact modular congruence: `add1c a b ≡ a + b (mod 2^32 - 1)`. -/
theorem add1c_mod_mersenne (a b : BitVec 32) :
    (add1c a b).toNat % 4294967295 = (a.toNat + b.toNat) % 4294967295 := by
  rw [add1c_toNat]
  have ha := a.isLt
  have hb := b.isLt
  split
  · rfl
  · omega

/-! ## 3. Seaman/Pence `CHECKSUM` Complement Theorem (`checksum.rs:147-148,173-174`) -/

/-- Exact Seaman/Pence complement identity: `add1c x (~~~x) = 0xFFFFFFFF` for EVERY `x : BitVec 32`.
    This guarantees that inserting the encoded `!total` into `CHECKSUM` forces the entire
    HDU 1's-complement sum to `0xFFFFFFFF` (`-0`, i.e., `0 mod (2^32 - 1)`). -/
theorem add1c_complement (x : BitVec 32) :
    add1c x (~~~x) = 0xFFFFFFFF#32 := by
  unfold add1c
  bv_decide

/-- Corollary: `add1c x (~~~x) ≡ 0 (mod 2^32 - 1)` (zero net HDU checksum). -/
theorem add1c_complement_mod_zero (x : BitVec 32) :
    (add1c x (~~~x)).toNat % 4294967295 = 0 := by
  rw [add1c_complement]
  decide

/-- Full HDU checksum verification soundness (`checksum.rs:147-148,173-174`):
    Given header sum `hdrSum` (with `CHECKSUM='0000000000000000'`), data sum `dsum`,
    `total = add1c hdrSum dsum`, replacing the placeholder with `!total` adds `!total`
    to `hdrSum`, making the final verified HDU total equal `0xFFFFFFFF`. -/
theorem hdu_checksum_verification_sound (hdrSum dsum : BitVec 32) :
    let total := add1c hdrSum dsum
    let newHdrSum := add1c hdrSum (~~~total)
    add1c newHdrSum dsum = 0xFFFFFFFF#32 := by
  dsimp only
  calc
    add1c (add1c hdrSum (~~~(add1c hdrSum dsum))) dsum
      = add1c hdrSum (add1c (~~~(add1c hdrSum dsum)) dsum) := add1c_assoc _ _ _
    _ = add1c hdrSum (add1c dsum (~~~(add1c hdrSum dsum))) := by rw [add1c_comm (~~~(add1c hdrSum dsum)) dsum]
    _ = add1c (add1c hdrSum dsum) (~~~(add1c hdrSum dsum)) := (add1c_assoc _ _ _).symm
    _ = 0xFFFFFFFF#32 := add1c_complement (add1c hdrSum dsum)

/-- Word-list fold matching `accumulate_checksum` (`checksum.rs:18-38`). -/
def fold1c (init : BitVec 32) (ws : List (BitVec 32)) : BitVec 32 :=
  ws.foldl add1c init

/-- Appending trailing zero words (e.g. FITS 2880-byte zero padding) does NOT change `DATASUM`
    (`checksum.rs:42`). -/
theorem fold1c_append_zeros (init : BitVec 32) (ws : List (BitVec 32)) (k : Nat) :
    fold1c init (ws ++ List.replicate k 0#32) = fold1c init ws := by
  unfold fold1c
  rw [List.foldl_append]
  induction k generalizing init with
  | zero => rfl
  | succ k ih =>
    simp only [List.replicate_succ, List.foldl_cons, add1c_zero_right]
    exact ih _

/-- Concatenation homomorphic property of `fold1c` (allowing independent `DATASUM` + `HDRSUM` combination). -/
theorem fold1c_append (a b : List (BitVec 32)) :
    fold1c 0#32 (a ++ b) = add1c (fold1c 0#32 a) (fold1c 0#32 b) := by
  unfold fold1c
  rw [List.foldl_append]
  suffices ∀ (acc : BitVec 32) (xs : List (BitVec 32)),
      xs.foldl add1c acc = add1c acc (xs.foldl add1c 0#32) by
    exact this (a.foldl add1c 0#32) b
  intro acc xs
  induction xs generalizing acc with
  | nil => simp [add1c_zero_right]
  | cons x xs ih =>
    simp only [List.foldl_cons, add1c_zero_left]
    rw [ih (add1c acc x), ih x, add1c_assoc]

end AstrosatFormal.FitsChecksum
