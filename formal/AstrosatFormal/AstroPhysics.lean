/-!
# AstrosatFormal.AstroPhysics

Formal verification of the core astrophysical & timing invariants of `astrosat-rs`:
1. **LAXPC Non-Paralyzable Deadtime (`tau = 43 / 1_000_000`)**
   Mirroring `crates/astrosat-laxpc/src/types.rs:305,308` & `products.rs:25-45`:
   - Strictly positive denominator `1 - tau * R_all > 0` for `0 ≤ R_all < 23000`
   - `R_corr ≥ R_obs ≥ 0`
   - Strict monotonicity in both `R_obs` and `R_all`
2. **CZTI Mask-Weight Renormalization Conservation**
   Mirroring `crates/astrosat-czti/src/bindata.rs:519-580`:
   - Exact algebraic background cancellation `∑_{i ∈ good} w_i * B_i * T_i = 0` identically
     for any finite good-pixel set with `∑ B_i * T_i > 0` and `denom ≠ 0`.
3. **AstroSat Time & Leap-Second Exactness + CZTI 20-Bit Rollover**
   Mirroring `crates/astrosat-core/src/time.rs:7-106` & `tct.rs:11,182-193,254-289`:
   - Exact `UTC ↔ TAI ↔ TT` invertibility (including the exact 3-step fixed-point loop
     `tt_sec_to_astrosat_utc`) across all AstroSat leap-second epochs
   - Exact 1-second continuity jump across 2012-07-01, 2015-07-01, and 2017-01-01 boundaries
   - Exact CZTI 20-bit rollover `mod 1048576` preservation
-/

namespace AstrosatFormal.AstroPhysics

/-! ## 0. Exact Rational Arithmetic (`Q`) over `Int` (Pure Lean 4 Core, Zero Dependencies) -/

/-- Exact rational number `num / den` with strictly positive integer denominator `den > 0`. -/
structure Q where
  num     : Int
  den     : Int
  den_pos : 0 < den
  deriving Repr

namespace Q

def ofInt (n : Int) : Q := ⟨n, 1, by omega⟩
def zero : Q := ofInt 0
def one  : Q := ofInt 1

/-- Rational equality by cross-multiplication (`a.num * b.den = b.num * a.den`). -/
def equiv (a b : Q) : Prop := a.num * b.den = b.num * a.den

infix:50 " ≈ " => Q.equiv

def lt (a b : Q) : Prop := a.num * b.den < b.num * a.den
def le (a b : Q) : Prop := a.num * b.den ≤ b.num * a.den

instance : LT Q := ⟨Q.lt⟩
instance : LE Q := ⟨Q.le⟩

instance (a b : Q) : Decidable (a < b) := inferInstanceAs (Decidable (a.num * b.den < b.num * a.den))
instance (a b : Q) : Decidable (a ≤ b) := inferInstanceAs (Decidable (a.num * b.den ≤ b.num * a.den))
instance (a b : Q) : Decidable (a ≈ b) := inferInstanceAs (Decidable (a.num * b.den = b.num * a.den))

def neg (a : Q) : Q := ⟨-a.num, a.den, a.den_pos⟩

def add (a b : Q) : Q :=
  ⟨a.num * b.den + b.num * a.den, a.den * b.den, Int.mul_pos a.den_pos b.den_pos⟩

def sub (a b : Q) : Q :=
  ⟨a.num * b.den - b.num * a.den, a.den * b.den, Int.mul_pos a.den_pos b.den_pos⟩

def mul (a b : Q) : Q :=
  ⟨a.num * b.num, a.den * b.den, Int.mul_pos a.den_pos b.den_pos⟩

/-- Safe exact division by a non-zero rational `b` (`b.num ≠ 0`), preserving `den > 0`. -/
def div (a b : Q) (hb : b.num ≠ 0) : Q :=
  if hpos : 0 < b.num then
    ⟨a.num * b.den, a.den * b.num, Int.mul_pos a.den_pos hpos⟩
  else
    ⟨-(a.num * b.den), a.den * (-b.num), Int.mul_pos a.den_pos (by omega)⟩

instance : Add Q := ⟨Q.add⟩
instance : Sub Q := ⟨Q.sub⟩
instance : Mul Q := ⟨Q.mul⟩
instance : Neg Q := ⟨Q.neg⟩

theorem equiv_refl (a : Q) : a ≈ a := rfl

theorem equiv_symm {a b : Q} (h : a ≈ b) : b ≈ a := h.symm

theorem equiv_trans {a b c : Q} (h1 : a ≈ b) (h2 : b ≈ c) : a ≈ c := by
  unfold equiv at *
  have hb : b.den ≠ 0 := by have := b.den_pos; omega
  have hmul : a.num * c.den * b.den = c.num * a.den * b.den := by
    calc a.num * c.den * b.den
        = (a.num * b.den) * c.den := by rw [Int.mul_right_comm]
      _ = (b.num * a.den) * c.den := by rw [h1]
      _ = (b.num * c.den) * a.den := by rw [Int.mul_right_comm]
      _ = (c.num * b.den) * a.den := by rw [h2]
      _ = c.num * a.den * b.den := by rw [Int.mul_right_comm]
  exact Int.eq_of_mul_eq_mul_right hb hmul

theorem equiv_zero_iff (a : Q) : a ≈ zero ↔ a.num = 0 := by
  unfold equiv zero ofInt
  simp

theorem equiv_add {a₁ a₂ b₁ b₂ : Q} (ha : a₁ ≈ a₂) (hb : b₁ ≈ b₂) :
    a₁ + b₁ ≈ a₂ + b₂ := by
  unfold equiv at *
  show (a₁.num * b₁.den + b₁.num * a₁.den) * (a₂.den * b₂.den) =
       (a₂.num * b₂.den + b₂.num * a₂.den) * (a₁.den * b₁.den)
  have t1 : a₁.num * b₁.den * (a₂.den * b₂.den) = a₂.num * b₂.den * (a₁.den * b₁.den) := by
    calc a₁.num * b₁.den * (a₂.den * b₂.den)
        = (a₁.num * a₂.den) * (b₁.den * b₂.den) := by
          rw [Int.mul_assoc, Int.mul_left_comm b₁.den a₂.den, ← Int.mul_assoc]
      _ = (a₂.num * a₁.den) * (b₁.den * b₂.den) := by rw [ha]
      _ = a₂.num * b₂.den * (a₁.den * b₁.den) := by
          rw [Int.mul_assoc, Int.mul_comm b₁.den b₂.den, Int.mul_left_comm a₁.den b₂.den,
              ← Int.mul_assoc]
  have t2 : b₁.num * a₁.den * (a₂.den * b₂.den) = b₂.num * a₂.den * (a₁.den * b₁.den) := by
    calc b₁.num * a₁.den * (a₂.den * b₂.den)
        = (b₁.num * b₂.den) * (a₁.den * a₂.den) := by
          rw [Int.mul_assoc, Int.mul_comm a₂.den b₂.den, Int.mul_left_comm a₁.den b₂.den,
              ← Int.mul_assoc]
      _ = (b₂.num * b₁.den) * (a₁.den * a₂.den) := by rw [hb]
      _ = b₂.num * a₂.den * (a₁.den * b₁.den) := by
          rw [Int.mul_assoc, Int.mul_comm a₁.den a₂.den, Int.mul_left_comm b₁.den a₂.den,
              ← Int.mul_assoc, Int.mul_comm b₁.den a₁.den]
  rw [Int.add_mul, Int.add_mul, t1, t2]

theorem equiv_sub {a₁ a₂ b₁ b₂ : Q} (ha : a₁ ≈ a₂) (hb : b₁ ≈ b₂) :
    a₁ - b₁ ≈ a₂ - b₂ := by
  unfold equiv at *
  show (a₁.num * b₁.den - b₁.num * a₁.den) * (a₂.den * b₂.den) =
       (a₂.num * b₂.den - b₂.num * a₂.den) * (a₁.den * b₁.den)
  have t1 : a₁.num * b₁.den * (a₂.den * b₂.den) = a₂.num * b₂.den * (a₁.den * b₁.den) := by
    calc a₁.num * b₁.den * (a₂.den * b₂.den)
        = (a₁.num * a₂.den) * (b₁.den * b₂.den) := by
          rw [Int.mul_assoc, Int.mul_left_comm b₁.den a₂.den, ← Int.mul_assoc]
      _ = (a₂.num * a₁.den) * (b₁.den * b₂.den) := by rw [ha]
      _ = a₂.num * b₂.den * (a₁.den * b₁.den) := by
          rw [Int.mul_assoc, Int.mul_comm b₁.den b₂.den, Int.mul_left_comm a₁.den b₂.den,
              ← Int.mul_assoc]
  have t2 : b₁.num * a₁.den * (a₂.den * b₂.den) = b₂.num * a₂.den * (a₁.den * b₁.den) := by
    calc b₁.num * a₁.den * (a₂.den * b₂.den)
        = (b₁.num * b₂.den) * (a₁.den * a₂.den) := by
          rw [Int.mul_assoc, Int.mul_comm a₂.den b₂.den, Int.mul_left_comm a₁.den b₂.den,
              ← Int.mul_assoc]
      _ = (b₂.num * b₁.den) * (a₁.den * a₂.den) := by rw [hb]
      _ = b₂.num * a₂.den * (a₁.den * b₁.den) := by
          rw [Int.mul_assoc, Int.mul_comm a₁.den a₂.den, Int.mul_left_comm b₁.den a₂.den,
              ← Int.mul_assoc, Int.mul_comm b₁.den a₁.den]
  rw [Int.sub_mul, Int.sub_mul, t1, t2]

theorem equiv_mul {a₁ a₂ b₁ b₂ : Q} (ha : a₁ ≈ a₂) (hb : b₁ ≈ b₂) :
    a₁ * b₁ ≈ a₂ * b₂ := by
  unfold equiv at *
  show (a₁.num * b₁.num) * (a₂.den * b₂.den) = (a₂.num * b₂.num) * (a₁.den * b₁.den)
  calc (a₁.num * b₁.num) * (a₂.den * b₂.den)
      = (a₁.num * a₂.den) * (b₁.num * b₂.den) := by
        rw [Int.mul_assoc, Int.mul_left_comm b₁.num a₂.den, ← Int.mul_assoc]
    _ = (a₂.num * a₁.den) * (b₂.num * b₁.den) := by rw [ha, hb]
    _ = (a₂.num * b₂.num) * (a₁.den * b₁.den) := by
        rw [Int.mul_assoc, Int.mul_left_comm a₁.den b₂.num, ← Int.mul_assoc]

theorem add_zero_left (a : Q) : zero + a ≈ a := by
  unfold equiv zero ofInt
  show (0 * a.den + a.num * 1) * a.den = a.num * (1 * a.den)
  simp

theorem sub_self_equiv_zero (a : Q) : a - a ≈ zero := by
  unfold equiv zero ofInt
  show (a.num * a.den - a.num * a.den) * 1 = 0 * (a.den * a.den)
  omega

theorem div_mul_cancel (a b : Q) (hb : b.num ≠ 0) : (a.div b hb) * b ≈ a := by
  unfold equiv div
  split
  · show (a.num * b.den * b.num) * a.den = a.num * (a.den * b.num * b.den)
    rw [Int.mul_assoc a.num b.den b.num, Int.mul_assoc a.num,
        Int.mul_comm (b.den * b.num) a.den, Int.mul_assoc a.den, Int.mul_comm b.den b.num]
  · show (-(a.num * b.den) * b.num) * a.den = a.num * (a.den * (-b.num) * b.den)
    rw [Int.neg_mul, Int.neg_mul, Int.mul_neg, Int.neg_mul, Int.mul_neg]
    congr 1
    rw [Int.mul_assoc a.num b.den b.num, Int.mul_assoc a.num,
        Int.mul_comm (b.den * b.num) a.den, Int.mul_assoc a.den, Int.mul_comm b.den b.num]

theorem div_add_dist (a b c : Q) (hc : c.num ≠ 0) :
    a.div c hc + b.div c hc ≈ (a + b).div c hc := by
  unfold equiv div
  split
  · show (a.num * c.den * (b.den * c.num) + b.num * c.den * (a.den * c.num)) *
         (a.den * b.den * c.num) =
         ((a.num * b.den + b.num * a.den) * c.den) * (a.den * c.num * (b.den * c.num))
    rw [Int.add_mul, Int.add_mul, Int.add_mul]
    congr 1 <;> (
      simp only [Int.mul_assoc]
      ac_rfl)
  · show (-(a.num * c.den) * (b.den * -c.num) + -(b.num * c.den) * (a.den * -c.num)) *
         (a.den * b.den * -c.num) =
         (-((a.num * b.den + b.num * a.den) * c.den)) * (a.den * -c.num * (b.den * -c.num))
    rw [Int.add_mul, Int.add_mul, Int.neg_add, Int.add_mul]
    congr 1 <;> (
      simp only [Int.neg_mul, Int.mul_neg, Int.neg_neg, Int.mul_assoc]
      congr 1
      ac_rfl)

theorem sub_mul_dist (a b c : Q) : (a - b) * c ≈ a * c - b * c := by
  unfold equiv
  show ((a.num * b.den - b.num * a.den) * c.num) * (a.den * c.den * (b.den * c.den)) =
       (a.num * c.num * (b.den * c.den) - b.num * c.num * (a.den * c.den)) * (a.den * b.den * c.den)
  rw [Int.sub_mul, Int.sub_mul, Int.sub_mul]
  congr 1 <;> (
    simp only [Int.mul_assoc]
    ac_rfl)

theorem div_mul_assoc (a b c : Q) (hb : b.num ≠ 0) :
    (a.div b hb) * c ≈ (a * c).div b hb := by
  unfold equiv div
  split
  · show (a.num * b.den * c.num) * (a.den * c.den * b.num) =
         (a.num * c.num * b.den) * (a.den * b.num * c.den)
    simp only [Int.mul_assoc]; ac_rfl
  · show (-(a.num * b.den) * c.num) * (a.den * c.den * -b.num) =
         (-(a.num * c.num * b.den)) * (a.den * -b.num * c.den)
    simp only [Int.neg_mul, Int.mul_neg, Int.neg_neg, Int.mul_assoc]; ac_rfl

theorem mul_add_dist (a b c : Q) : a * (b + c) ≈ a * b + a * c := by
  unfold equiv
  show (a.num * (b.num * c.den + c.num * b.den)) * (a.den * b.den * (a.den * c.den)) =
       (a.num * b.num * (a.den * c.den) + a.num * c.num * (a.den * b.den)) * (a.den * (b.den * c.den))
  rw [Int.mul_add, Int.add_mul, Int.add_mul]
  congr 1 <;> (simp only [Int.mul_assoc]; ac_rfl)

theorem div_zero_left (b : Q) (hb : b.num ≠ 0) : zero.div b hb ≈ zero := by
  unfold equiv div zero ofInt
  split <;> simp

theorem equiv_div_left {a₁ a₂ b : Q} (ha : a₁ ≈ a₂) (hb : b.num ≠ 0) :
    a₁.div b hb ≈ a₂.div b hb := by
  unfold equiv div at *
  split
  · show (a₁.num * b.den) * (a₂.den * b.num) = (a₂.num * b.den) * (a₁.den * b.num)
    calc (a₁.num * b.den) * (a₂.den * b.num)
        = (a₁.num * a₂.den) * (b.den * b.num) := by simp only [Int.mul_assoc]; ac_rfl
      _ = (a₂.num * a₁.den) * (b.den * b.num) := by rw [ha]
      _ = (a₂.num * b.den) * (a₁.den * b.num) := by simp only [Int.mul_assoc]; ac_rfl
  · show (-(a₁.num * b.den)) * (a₂.den * -b.num) = (-(a₂.num * b.den)) * (a₁.den * -b.num)
    simp only [Int.neg_mul, Int.mul_neg, Int.neg_neg]
    calc (a₁.num * b.den) * (a₂.den * b.num)
        = (a₁.num * a₂.den) * (b.den * b.num) := by simp only [Int.mul_assoc]; ac_rfl
      _ = (a₂.num * a₁.den) * (b.den * b.num) := by rw [ha]
      _ = (a₂.num * b.den) * (a₁.den * b.num) := by simp only [Int.mul_assoc]; ac_rfl

theorem equiv_sub_cancel_add (a b : Q) : (a + b) - b ≈ a := by
  unfold equiv
  show ((a.num * b.den + b.num * a.den) * b.den - b.num * (a.den * b.den)) * a.den =
       a.num * (a.den * b.den * b.den)
  rw [Int.add_mul, Int.mul_assoc b.num a.den b.den]
  have hsub : a.num * b.den * b.den + b.num * (a.den * b.den) - b.num * (a.den * b.den) =
      a.num * b.den * b.den := by omega
  rw [hsub]
  simp only [Int.mul_assoc]; ac_rfl

theorem equiv_add_sub_cancel (a b : Q) : (a - b) + b ≈ a := by
  unfold equiv
  show ((a.num * b.den - b.num * a.den) * b.den + b.num * (a.den * b.den)) * a.den =
       a.num * (a.den * b.den * b.den)
  rw [Int.sub_mul, Int.mul_assoc b.num a.den b.den]
  have hcancel : a.num * b.den * b.den - b.num * (a.den * b.den) + b.num * (a.den * b.den) =
      a.num * b.den * b.den := by omega
  rw [hcancel]
  simp only [Int.mul_assoc]; ac_rfl

theorem le_iff_of_equiv {a₁ a₂ b₁ b₂ : Q} (ha : a₁ ≈ a₂) (hb : b₁ ≈ b₂) :
    a₁ ≤ b₁ ↔ a₂ ≤ b₂ := by
  show a₁.num * b₁.den ≤ b₁.num * a₁.den ↔ a₂.num * b₂.den ≤ b₂.num * a₂.den
  unfold equiv at ha hb
  have hpos1 : 0 < a₂.den * b₂.den := Int.mul_pos a₂.den_pos b₂.den_pos
  have hpos2 : 0 < a₁.den * b₁.den := Int.mul_pos a₁.den_pos b₁.den_pos
  have hl : (a₁.num * b₁.den) * (a₂.den * b₂.den) = (a₂.num * b₂.den) * (a₁.den * b₁.den) := by
    calc (a₁.num * b₁.den) * (a₂.den * b₂.den)
        = (a₁.num * a₂.den) * (b₁.den * b₂.den) := by simp only [Int.mul_assoc]; ac_rfl
      _ = (a₂.num * a₁.den) * (b₁.den * b₂.den) := by rw [ha]
      _ = (a₂.num * b₂.den) * (a₁.den * b₁.den) := by simp only [Int.mul_assoc]; ac_rfl
  have hr : (b₁.num * a₁.den) * (a₂.den * b₂.den) = (b₂.num * a₂.den) * (a₁.den * b₁.den) := by
    calc (b₁.num * a₁.den) * (a₂.den * b₂.den)
        = (b₁.num * b₂.den) * (a₁.den * a₂.den) := by simp only [Int.mul_assoc]; ac_rfl
      _ = (b₂.num * b₁.den) * (a₁.den * a₂.den) := by rw [hb]
      _ = (b₂.num * a₂.den) * (a₁.den * b₁.den) := by simp only [Int.mul_assoc]; ac_rfl
  constructor
  · intro h
    have hm := Int.mul_le_mul_of_nonneg_right h (Int.le_of_lt hpos1)
    rw [hl, hr] at hm
    exact Int.le_of_mul_le_mul_right hm hpos2
  · intro h
    have hm := Int.mul_le_mul_of_nonneg_right h (Int.le_of_lt hpos2)
    rw [← hl, ← hr] at hm
    exact Int.le_of_mul_le_mul_right hm hpos1

end Q

/-! ## 1. LAXPC Non-Paralyzable Deadtime (`types.rs:305,308`, `products.rs:25-45`) -/

/-- Exact LAXPC deadtime constant `tau = 43 us = 43 / 1_000_000 s` (`types.rs:305`). -/
def laxpcTau : Q := ⟨43, 1000000, by omega⟩

/-- Exact LAXPC deadtime rate singularity guard `23000 cps` (`types.rs:308`). -/
def laxpcRateGuard : Q := Q.ofInt 23000

/-- Non-paralyzable deadtime denominator `1 - tau * R_all` (`products.rs:27,39`). -/
def deadtimeDenom (rAll : Q) : Q := Q.one - laxpcTau * rAll

/-- Theorem 1(a): For `0 ≤ R_all < 23000`, `1 - tau * R_all > 0` (strictly `> 11/1000`),
    guaranteeing no division by zero, NaN, or negative rate in `deadtime_correct_rate`. -/
theorem deadtimeDenom_pos (rAll : Q)
    (hlo : Q.zero ≤ rAll) (hhi : rAll < laxpcRateGuard) :
    Q.zero < deadtimeDenom rAll := by
  have h1 : 0 ≤ rAll.num := by
    change 0 * rAll.den ≤ rAll.num * 1 at hlo
    omega
  have h2 : rAll.num < 23000 * rAll.den := by
    change rAll.num * 1 < 23000 * rAll.den at hhi
    omega
  have hden := rAll.den_pos
  show 0 * (1 * (1000000 * rAll.den)) < (1 * (1000000 * rAll.den) - 43 * rAll.num * 1) * 1
  omega

/-- Stronger quantitative lower bound: `1 - tau * R_all > 11 / 1000` for `R_all < 23000`. -/
theorem deadtimeDenom_gt_eleven_thousandths (rAll : Q) (hhi : rAll < laxpcRateGuard) :
    (⟨11, 1000, by omega⟩ : Q) < deadtimeDenom rAll := by
  have h2 : rAll.num < 23000 * rAll.den := by
    change rAll.num * 1 < 23000 * rAll.den at hhi
    omega
  have hden := rAll.den_pos
  show 11 * (1 * (1000000 * rAll.den)) < (1 * (1000000 * rAll.den) - 43 * rAll.num * 1) * 1000
  omega

/-- Helper: positive `deadtimeDenom` has `num ≠ 0`. -/
theorem deadtimeDenom_num_pos (rAll : Q)
    (hlo : Q.zero ≤ rAll) (hhi : rAll < laxpcRateGuard) :
    0 < (deadtimeDenom rAll).num := by
  have h := deadtimeDenom_pos rAll hlo hhi
  change 0 * (deadtimeDenom rAll).den < (deadtimeDenom rAll).num * 1 at h
  omega

/-- Exact deadtime-corrected rate `R_corr = R_obs / (1 - tau * R_all)` (`products.rs:29`). -/
def deadtimeCorrectRate (rObs rAll : Q)
    (hlo : Q.zero ≤ rAll) (hhi : rAll < laxpcRateGuard) : Q :=
  rObs.div (deadtimeDenom rAll) (by have := deadtimeDenom_num_pos rAll hlo hhi; omega)

/-- Exact deadtime factor `f_cd = 1 / (1 - tau * R_all)` (`products.rs:41`). -/
def deadtimeFactor (rAll : Q)
    (hlo : Q.zero ≤ rAll) (hhi : rAll < laxpcRateGuard) : Q :=
  Q.one.div (deadtimeDenom rAll) (by have := deadtimeDenom_num_pos rAll hlo hhi; omega)

/-- Theorem 1(b): For `R_obs ≥ 0` and `0 ≤ R_all < 23000`, `R_corr ≥ R_obs ≥ 0`
    (deadtime correction preserves non-negativity and never underestimates observed rate). -/
theorem deadtimeCorrectRate_ge_obs (rObs rAll : Q)
    (hobs : Q.zero ≤ rObs) (hlo : Q.zero ≤ rAll) (hhi : rAll < laxpcRateGuard) :
    rObs ≤ deadtimeCorrectRate rObs rAll hlo hhi ∧ Q.zero ≤ rObs := by
  refine ⟨?_, hobs⟩
  have hdpos := deadtimeDenom_num_pos rAll hlo hhi
  unfold deadtimeCorrectRate Q.div
  simp only [hdpos, ↓reduceDIte]
  have hro : 0 ≤ rObs.num := by
    change 0 * rObs.den ≤ rObs.num * 1 at hobs; omega
  have hra : 0 ≤ rAll.num := by
    change 0 * rAll.den ≤ rAll.num * 1 at hlo; omega
  show rObs.num * (rObs.den * (deadtimeDenom rAll).num) ≤
       (rObs.num * (deadtimeDenom rAll).den) * rObs.den
  have hsub : (deadtimeDenom rAll).num ≤ (deadtimeDenom rAll).den := by
    show 1 * (1000000 * rAll.den) - 43 * rAll.num * 1 ≤ 1 * (1000000 * rAll.den)
    omega
  have hprod : 0 ≤ rObs.num * rObs.den := Int.mul_nonneg hro (Int.le_of_lt rObs.den_pos)
  have hm := Int.mul_le_mul_of_nonneg_left hsub hprod
  calc rObs.num * (rObs.den * (deadtimeDenom rAll).num)
      = (rObs.num * rObs.den) * (deadtimeDenom rAll).num := by rw [Int.mul_assoc]
    _ ≤ (rObs.num * rObs.den) * (deadtimeDenom rAll).den := hm
    _ = (rObs.num * (deadtimeDenom rAll).den) * rObs.den := by
        rw [Int.mul_assoc, Int.mul_comm rObs.den (deadtimeDenom rAll).den, ← Int.mul_assoc]

/-- Theorem 1(c-i): Strict monotonicity of deadtime correction in `R_obs` (for fixed `0 ≤ R_all < 23000`). -/
theorem deadtimeCorrectRate_strictMono_obs (rObs₁ rObs₂ rAll : Q)
    (hlt : rObs₁ < rObs₂) (hlo : Q.zero ≤ rAll) (hhi : rAll < laxpcRateGuard) :
    deadtimeCorrectRate rObs₁ rAll hlo hhi < deadtimeCorrectRate rObs₂ rAll hlo hhi := by
  have hdpos := deadtimeDenom_num_pos rAll hlo hhi
  have hdden := (deadtimeDenom rAll).den_pos
  unfold deadtimeCorrectRate Q.div
  simp only [hdpos, ↓reduceDIte]
  show (rObs₁.num * (deadtimeDenom rAll).den) * (rObs₂.den * (deadtimeDenom rAll).num) <
       (rObs₂.num * (deadtimeDenom rAll).den) * (rObs₁.den * (deadtimeDenom rAll).num)
  have h1 : rObs₁.num * rObs₂.den < rObs₂.num * rObs₁.den := hlt
  have hfac : 0 < (deadtimeDenom rAll).den * (deadtimeDenom rAll).num :=
    Int.mul_pos hdden hdpos
  have hm := Int.mul_lt_mul_of_pos_right h1 hfac
  calc (rObs₁.num * (deadtimeDenom rAll).den) * (rObs₂.den * (deadtimeDenom rAll).num)
      = (rObs₁.num * rObs₂.den) * ((deadtimeDenom rAll).den * (deadtimeDenom rAll).num) := by
        simp only [Int.mul_assoc]; ac_rfl
    _ < (rObs₂.num * rObs₁.den) * ((deadtimeDenom rAll).den * (deadtimeDenom rAll).num) := hm
    _ = (rObs₂.num * (deadtimeDenom rAll).den) * (rObs₁.den * (deadtimeDenom rAll).num) := by
        simp only [Int.mul_assoc]; ac_rfl

/-- Theorem 1(c-ii): Strict monotonicity of deadtime correction in `R_all`
    (for `R_obs > 0` and `0 ≤ R_all₁ < R_all₂ < 23000`). -/
theorem deadtimeCorrectRate_strictMono_all (rObs rAll₁ rAll₂ : Q)
    (hobs : Q.zero < rObs)
    (hlo₁ : Q.zero ≤ rAll₁) (hhi₁ : rAll₁ < laxpcRateGuard)
    (hlo₂ : Q.zero ≤ rAll₂) (hhi₂ : rAll₂ < laxpcRateGuard)
    (hlt : rAll₁ < rAll₂) :
    deadtimeCorrectRate rObs rAll₁ hlo₁ hhi₁ < deadtimeCorrectRate rObs rAll₂ hlo₂ hhi₂ := by
  have hdpos₁ := deadtimeDenom_num_pos rAll₁ hlo₁ hhi₁
  have hdpos₂ := deadtimeDenom_num_pos rAll₂ hlo₂ hhi₂
  unfold deadtimeCorrectRate Q.div
  simp only [hdpos₁, hdpos₂, ↓reduceDIte]
  have hro : 0 < rObs.num := by
    change 0 * rObs.den < rObs.num * 1 at hobs; omega
  have hlt_int : rAll₁.num * rAll₂.den < rAll₂.num * rAll₁.den := hlt
  have hdenom_lt : (deadtimeDenom rAll₂).num * (deadtimeDenom rAll₁).den <
                   (deadtimeDenom rAll₁).num * (deadtimeDenom rAll₂).den := by
    show (1 * (1000000 * rAll₂.den) - 43 * rAll₂.num * 1) * (1 * (1000000 * rAll₁.den)) <
         (1 * (1000000 * rAll₁.den) - 43 * rAll₁.num * 1) * (1 * (1000000 * rAll₂.den))
    have h43 : 43 * (rAll₁.num * rAll₂.den) < 43 * (rAll₂.num * rAll₁.den) := by omega
    have hscale := Int.mul_lt_mul_of_pos_left h43 (by omega : (0 : Int) < 1000000)
    simp only [Int.one_mul, Int.mul_one, Int.sub_mul]
    have e1 : 1000000 * rAll₂.den * (1000000 * rAll₁.den) =
              1000000 * rAll₁.den * (1000000 * rAll₂.den) := Int.mul_comm _ _
    have e2 : 43 * rAll₂.num * (1000000 * rAll₁.den) =
              1000000 * (43 * (rAll₂.num * rAll₁.den)) := by simp only [Int.mul_assoc]; ac_rfl
    have e3 : 43 * rAll₁.num * (1000000 * rAll₂.den) =
              1000000 * (43 * (rAll₁.num * rAll₂.den)) := by simp only [Int.mul_assoc]; ac_rfl
    rw [e1, e2, e3]
    omega
  show (rObs.num * (deadtimeDenom rAll₁).den) * (rObs.den * (deadtimeDenom rAll₂).num) <
       (rObs.num * (deadtimeDenom rAll₂).den) * (rObs.den * (deadtimeDenom rAll₁).num)
  have hfac : 0 < rObs.num * rObs.den := Int.mul_pos hro rObs.den_pos
  have hm := Int.mul_lt_mul_of_pos_left hdenom_lt hfac
  calc (rObs.num * (deadtimeDenom rAll₁).den) * (rObs.den * (deadtimeDenom rAll₂).num)
      = (rObs.num * rObs.den) * ((deadtimeDenom rAll₂).num * (deadtimeDenom rAll₁).den) := by
        simp only [Int.mul_assoc]; ac_rfl
    _ < (rObs.num * rObs.den) * ((deadtimeDenom rAll₁).num * (deadtimeDenom rAll₂).den) := hm
    _ = (rObs.num * (deadtimeDenom rAll₂).den) * (rObs.den * (deadtimeDenom rAll₁).num) := by
        simp only [Int.mul_assoc]; ac_rfl

/-- Corollary: `deadtimeFactor` is `≥ 1` and strictly increasing on `[0, 23000)`. -/
theorem deadtimeFactor_ge_one_and_strictMono (rAll₁ rAll₂ : Q)
    (hlo₁ : Q.zero ≤ rAll₁) (hhi₁ : rAll₁ < laxpcRateGuard)
    (hlo₂ : Q.zero ≤ rAll₂) (hhi₂ : rAll₂ < laxpcRateGuard)
    (hlt : rAll₁ < rAll₂) :
    Q.one ≤ deadtimeFactor rAll₁ hlo₁ hhi₁ ∧
    deadtimeFactor rAll₁ hlo₁ hhi₁ < deadtimeFactor rAll₂ hlo₂ hhi₂ := by
  have hge := (deadtimeCorrectRate_ge_obs Q.one rAll₁ (by decide) hlo₁ hhi₁).1
  have hmono := deadtimeCorrectRate_strictMono_all Q.one rAll₁ rAll₂ (by decide) hlo₁ hhi₁ hlo₂ hhi₂ hlt
  exact ⟨hge, hmono⟩

/-! ## 2. CZTI Mask-Weight Renormalization Conservation (`bindata.rs:519-580`) -/

/-- A single good pixel's raw mask weight `w'`, `BKGSCAL` factor `b`, and exposure fraction `t`. -/
structure GoodPixel where
  wRaw : Q
  bkg  : Q
  texp : Q

/-- Weighted background sum `det_renorm_a = ∑_{i ∈ good} B_i * T_i` (`bindata.rs:540`). -/
def sumBT : List GoodPixel → Q
  | [] => Q.zero
  | p :: ps => p.bkg * p.texp + sumBT ps

/-- Raw mask-weighted background sum `det_renorm_d = ∑_{i ∈ good} w'_i * B_i * T_i` (`bindata.rs:539`). -/
def sumWBT : List GoodPixel → Q
  | [] => Q.zero
  | p :: ps => p.wRaw * (p.bkg * p.texp) + sumWBT ps

/-- Renormalized mask weight `w_i = (w'_i - D) / denom` (`bindata.rs:577`). -/
def renormWeight (wRaw D denom : Q) (hdenom : denom.num ≠ 0) : Q :=
  (wRaw - D).div denom hdenom

/-- Weighted background sum under renormalized weights `∑_{i ∈ good} w_i * B_i * T_i`. -/
def sumRenormWBT (D denom : Q) (hdenom : denom.num ≠ 0) : List GoodPixel → Q
  | [] => Q.zero
  | p :: ps => renormWeight p.wRaw D denom hdenom * (p.bkg * p.texp) + sumRenormWBT D denom hdenom ps

/-- Algebraic helper `(a - b) + (c - d) ≈ (a + c) - (b + d)`. -/
theorem sub_add_sub_equiv (a b c d : Q) :
    (a - b) + (c - d) ≈ (a + c) - (b + d) := by
  unfold Q.equiv
  show ((a.num * b.den - b.num * a.den) * (c.den * d.den) +
        (c.num * d.den - d.num * c.den) * (a.den * b.den)) *
       (a.den * c.den * (b.den * d.den)) =
       ((a.num * c.den + c.num * a.den) * (b.den * d.den) -
        (b.num * d.den + d.num * b.den) * (a.den * c.den)) *
       (a.den * b.den * (c.den * d.den))
  rw [Int.sub_mul, Int.sub_mul, Int.add_mul, Int.add_mul, Int.add_mul, Int.sub_mul,
      Int.sub_mul, Int.sub_mul, Int.add_mul, Int.add_mul]
  have t1 : a.num * b.den * (c.den * d.den) * (a.den * c.den * (b.den * d.den)) =
            a.num * c.den * (b.den * d.den) * (a.den * b.den * (c.den * d.den)) := by
    simp only [Int.mul_assoc]; ac_rfl
  have t2 : b.num * a.den * (c.den * d.den) * (a.den * c.den * (b.den * d.den)) =
            b.num * d.den * (a.den * c.den) * (a.den * b.den * (c.den * d.den)) := by
    simp only [Int.mul_assoc]; ac_rfl
  have t3 : c.num * d.den * (a.den * b.den) * (a.den * c.den * (b.den * d.den)) =
            c.num * a.den * (b.den * d.den) * (a.den * b.den * (c.den * d.den)) := by
    simp only [Int.mul_assoc]; ac_rfl
  have t4 : d.num * c.den * (a.den * b.den) * (a.den * c.den * (b.den * d.den)) =
            d.num * b.den * (a.den * c.den) * (a.den * b.den * (c.den * d.den)) := by
    simp only [Int.mul_assoc]; ac_rfl
  rw [t1, t2, t3, t4]
  omega

/-- Exact algebraic decomposition of `sumRenormWBT` into `(sumWBT - D * sumBT) / denom`. -/
theorem sumRenormWBT_decomp (D denom : Q) (hdenom : denom.num ≠ 0) (ps : List GoodPixel) :
    sumRenormWBT D denom hdenom ps ≈ (sumWBT ps - D * sumBT ps).div denom hdenom := by
  induction ps with
  | nil =>
    simp only [sumRenormWBT, sumWBT, sumBT]
    have h0 : Q.zero - D * Q.zero ≈ Q.zero := by
      unfold Q.equiv Q.zero Q.ofInt
      show (0 * (D.den * 1) - D.num * 0 * 1) * 1 = 0 * (1 * (D.den * 1))
      simp
    exact Q.equiv_trans (Q.equiv_symm (Q.div_zero_left denom hdenom)) (Q.equiv_symm (Q.equiv_div_left h0 hdenom))
  | cons p rest ih =>
    simp only [sumRenormWBT, sumWBT, sumBT, renormWeight]
    let bt := p.bkg * p.texp
    have hhead : (p.wRaw - D).div denom hdenom * bt ≈ ((p.wRaw - D) * bt).div denom hdenom :=
      Q.div_mul_assoc (p.wRaw - D) denom bt hdenom
    have hhead2 : ((p.wRaw - D) * bt).div denom hdenom ≈ (p.wRaw * bt - D * bt).div denom hdenom :=
      Q.equiv_div_left (Q.sub_mul_dist p.wRaw D bt) hdenom
    have hstep1 :
        (p.wRaw - D).div denom hdenom * bt + sumRenormWBT D denom hdenom rest ≈
        (p.wRaw * bt - D * bt).div denom hdenom + (sumWBT rest - D * sumBT rest).div denom hdenom :=
      Q.equiv_add (Q.equiv_trans hhead hhead2) ih
    have hstep2 := Q.div_add_dist (p.wRaw * bt - D * bt) (sumWBT rest - D * sumBT rest) denom hdenom
    have hdist : D * (bt + sumBT rest) ≈ D * bt + D * sumBT rest := Q.mul_add_dist D bt (sumBT rest)
    have halg :
        (p.wRaw * bt - D * bt) + (sumWBT rest - D * sumBT rest) ≈
        (p.wRaw * bt + sumWBT rest) - D * (bt + sumBT rest) :=
      Q.equiv_trans (sub_add_sub_equiv (p.wRaw * bt) (D * bt) (sumWBT rest) (D * sumBT rest))
        (Q.equiv_sub (Q.equiv_refl _) (Q.equiv_symm hdist))
    exact Q.equiv_trans (Q.equiv_trans hstep1 hstep2) (Q.equiv_div_left halg hdenom)

/-- **CZTI Mask-Weight Renormalization Exact Background Cancellation Theorem** (`bindata.rs:539-578`):
    For ANY finite good-pixel set `ps` with `sumBT ps > 0` and any `denom ≠ 0`, setting
    `D = (∑ w'_i B_i T_i) / (∑ B_i T_i)` and `w_i = (w'_i - D) / denom` yields
    `∑_{i ∈ good} w_i * B_i * T_i = 0` IDENTICALLY (both `≈ 0` and exact numerator `= 0`). -/
theorem czti_mask_weight_background_cancellation
    (ps : List GoodPixel) (denom : Q)
    (hsumPos : Q.zero < sumBT ps)
    (hdenom : denom.num ≠ 0) :
    let hbt : (sumBT ps).num ≠ 0 := by
      change 0 * (sumBT ps).den < (sumBT ps).num * 1 at hsumPos; omega
    let D := (sumWBT ps).div (sumBT ps) hbt
    sumRenormWBT D denom hdenom ps ≈ Q.zero ∧
    (sumRenormWBT D denom hdenom ps).num = 0 := by
  have hbt : (sumBT ps).num ≠ 0 := by
    change 0 * (sumBT ps).den < (sumBT ps).num * 1 at hsumPos; omega
  dsimp only
  let D := (sumWBT ps).div (sumBT ps) hbt
  have hdecomp := sumRenormWBT_decomp D denom hdenom ps
  have hcancel : D * sumBT ps ≈ sumWBT ps := Q.div_mul_cancel (sumWBT ps) (sumBT ps) hbt
  have hzero : sumWBT ps - D * sumBT ps ≈ Q.zero :=
    Q.equiv_trans (Q.equiv_sub (Q.equiv_refl _) hcancel) (Q.sub_self_equiv_zero (sumWBT ps))
  have hdiv0 : (sumWBT ps - D * sumBT ps).div denom hdenom ≈ Q.zero :=
    Q.equiv_trans (Q.equiv_div_left hzero hdenom) (Q.div_zero_left denom hdenom)
  have hres := Q.equiv_trans hdecomp hdiv0
  exact ⟨hres, (Q.equiv_zero_iff _).mp hres⟩

/-! ## 3. AstroSat Time & Leap-Second Exactness (`time.rs:7-106`) + CZTI 20-Bit Rollover (`tct.rs:11,182-193`) -/

/-- Exact constants matching `crates/astrosat-core/src/time.rs:7-13`. -/
def ASTROSAT_MJDREFI : Int := 55197
def SECS_PER_DAY : Int := 86400
/-- Exact `TT - TAI = 32.184 s = 32184 / 1000 s` (`time.rs:11`). -/
def TT_MINUS_TAI : Q := ⟨32184, 1000, by omega⟩

/-- AstroSat-era leap-second boundary timestamps in AstroSat UTC seconds since MJD 55197.0 (`time.rs:46-49`):
    - 2009-01-01 (MJD 54832): `TAI - UTC = 34 s` (active at MJD 55197 = 2010-01-01)
    - 2012-07-01 (MJD 56109): `TAI - UTC = 35 s` (`(56109 - 55197) * 86400 = 78_796_800`)
    - 2015-07-01 (MJD 57204): `TAI - UTC = 36 s` (`(57204 - 55197) * 86400 = 173_404_800`)
    - 2017-01-01 (MJD 57754): `TAI - UTC = 37 s` (`(57754 - 55197) * 86400 = 220_924_800`) -/
def T_2009_01_01 : Q := Q.ofInt ((54832 - ASTROSAT_MJDREFI) * SECS_PER_DAY)
def T_2012_07_01 : Q := Q.ofInt ((56109 - ASTROSAT_MJDREFI) * SECS_PER_DAY)
def T_2015_07_01 : Q := Q.ofInt ((57204 - ASTROSAT_MJDREFI) * SECS_PER_DAY)
def T_2017_01_01 : Q := Q.ofInt ((57754 - ASTROSAT_MJDREFI) * SECS_PER_DAY)

/-- Exact mirror of `tai_minus_utc_astrosat` (`time.rs:53-68`) across the entire AstroSat era
    (`t ≥ 2009-01-01`, i.e., `mjd_utc ≥ 54832`). -/
def taiMinusUtcAstrosat (tUtc : Q) : Int :=
  if T_2017_01_01 ≤ tUtc then 37
  else if T_2015_07_01 ≤ tUtc then 36
  else if T_2012_07_01 ≤ tUtc then 35
  else 34

/-- `TAI = UTC + (TAI - UTC)`. -/
def astrosatUtcToTai (tUtc : Q) : Q :=
  tUtc + Q.ofInt (taiMinusUtcAstrosat tUtc)

/-- Exact mirror of `astrosat_utc_to_tt_sec` (`time.rs:92-95`):
    `TT = UTC + (TAI - UTC) + 32.184`. -/
def astrosatUtcToTt (tUtc : Q) : Q :=
  tUtc + Q.ofInt (taiMinusUtcAstrosat tUtc) + TT_MINUS_TAI

/-- Exact mirror of the 3-step fixed-point loop `tt_sec_to_astrosat_utc` (`time.rs:98-106`). -/
def ttToAstrosatUtc (tTt : Q) : Q :=
  let t0 := tTt - Q.ofInt 36 - TT_MINUS_TAI
  let t1 := tTt - Q.ofInt (taiMinusUtcAstrosat t0) - TT_MINUS_TAI
  let t2 := tTt - Q.ofInt (taiMinusUtcAstrosat t1) - TT_MINUS_TAI
  tTt - Q.ofInt (taiMinusUtcAstrosat t2) - TT_MINUS_TAI

/-- Helper: `taiMinusUtcAstrosat` is invariant under rational equivalence `a ≈ b`. -/
theorem taiMinusUtc_congr {a b : Q} (h : a ≈ b) :
    taiMinusUtcAstrosat a = taiMinusUtcAstrosat b := by
  unfold taiMinusUtcAstrosat
  have h1 := Q.le_iff_of_equiv (Q.equiv_refl T_2017_01_01) h
  have h2 := Q.le_iff_of_equiv (Q.equiv_refl T_2015_07_01) h
  have h3 := Q.le_iff_of_equiv (Q.equiv_refl T_2012_07_01) h
  by_cases c1 : T_2017_01_01 ≤ a
  · simp [c1, h1.mp c1]
  · have c1' : ¬ (T_2017_01_01 ≤ b) := fun hb => c1 (h1.mpr hb)
    simp only [c1, c1', ↓reduceIte]
    by_cases c2 : T_2015_07_01 ≤ a
    · simp [c2, h2.mp c2]
    · have c2' : ¬ (T_2015_07_01 ≤ b) := fun hb => c2 (h2.mpr hb)
      simp only [c2, c2', ↓reduceIte]
      by_cases c3 : T_2012_07_01 ≤ a
      · simp [c3, h3.mp c3]
      · have c3' : ¬ (T_2012_07_01 ≤ b) := fun hb => c3 (h3.mpr hb)
        simp [c3, c3']

/-- Helper: `(t + ofInt k + TT_MINUS_TAI) - ofInt m - TT_MINUS_TAI ≈ t + ofInt (k - m)`. -/
theorem tt_step_simplify (t : Q) (k m : Int) :
    (t + Q.ofInt k + TT_MINUS_TAI) - Q.ofInt m - TT_MINUS_TAI ≈ t + Q.ofInt (k - m) := by
  show ((((t.num * 1 + k * t.den) * 1000 + 32184 * (t.den * 1)) * 1 -
          m * (t.den * 1 * 1000)) * 1000 - 32184 * (t.den * 1 * 1000 * 1)) * (t.den * 1) =
       (t.num * 1 + (k - m) * t.den) * (t.den * 1 * 1000 * 1 * 1000)
  simp only [Int.mul_one]
  have h1 : (((t.num + k * t.den) * 1000 + 32184 * t.den) - m * (t.den * 1000)) * 1000 -
            32184 * (t.den * 1000) =
            (t.num + (k - m) * t.den) * 1000000 := by
    rw [Int.sub_mul, Int.add_mul, Int.add_mul, Int.sub_mul, Int.add_mul]
    have e1 : 32184 * t.den * 1000 = 32184 * (t.den * 1000) := Int.mul_assoc _ _ _
    have e2 : m * (t.den * 1000) * 1000 = m * t.den * 1000000 := by simp only [Int.mul_assoc]; rfl
    have e3 : t.num * 1000 * 1000 = t.num * 1000000 := by simp only [Int.mul_assoc]; rfl
    have e4 : k * t.den * 1000 * 1000 = k * t.den * 1000000 := by simp only [Int.mul_assoc]; rfl
    rw [e1, e2, e3, e4]; omega
  rw [h1]
  have h2 : (1000000 : Int) * t.den = t.den * 1000 * 1000 := by omega
  rw [Int.mul_assoc, h2]

/-- Helper: integer threshold comparison on `t + Q.ofInt d`. -/
theorem ofInt_le_add_ofInt (c d : Int) (t : Q) :
    Q.ofInt c ≤ t + Q.ofInt d ↔ (c - d) * t.den ≤ t.num := by
  show c * (t.den * 1) ≤ (t.num * 1 + d * t.den) * 1 ↔ (c - d) * t.den ≤ t.num
  simp only [Int.mul_one]
  rw [Int.sub_mul]
  omega

/-- **Global Exact `UTC ↔ TT` Invertibility Theorem** (`time.rs:92-106`):
    For EVERY timestamp `t` across all 4 AstroSat-era leap-second epochs
    (2009–2012, 2012–2015, 2015–2017, and 2017–2026+, including AT and arbitrarily close to
    every leap-second boundary), the exact 3-iteration fixed-point routine
    `ttToAstrosatUtc (astrosatUtcToTt t)` recovers `t` EXACTLY (`≈ t`)! -/
theorem ttToAstrosatUtc_astrosatUtcToTt (t : Q) :
    ttToAstrosatUtc (astrosatUtcToTt t) ≈ t := by
  unfold ttToAstrosatUtc astrosatUtcToTt
  dsimp only
  have hden : 0 < t.den := t.den_pos
  have hzero_add : t + Q.ofInt 0 ≈ t := by
    unfold Q.equiv Q.ofInt; show (t.num * 1 + 0 * t.den) * t.den = t.num * (t.den * 1); simp
  by_cases c17 : T_2017_01_01 ≤ t
  · have hk : taiMinusUtcAstrosat t = 37 := by unfold taiMinusUtcAstrosat; simp [c17]
    have ht17 : 220924800 * t.den ≤ t.num := by
      change 220924800 * t.den ≤ t.num * 1 at c17; omega
    have h0_val : taiMinusUtcAstrosat (t + Q.ofInt (37 - 36)) = 37 := by
      unfold taiMinusUtcAstrosat T_2017_01_01 ASTROSAT_MJDREFI SECS_PER_DAY
      have hge : Q.ofInt 220924800 ≤ t + Q.ofInt 1 := by rw [ofInt_le_add_ofInt]; omega
      simp [hge]
    have h1t : t + Q.ofInt 37 + TT_MINUS_TAI - Q.ofInt 37 - TT_MINUS_TAI ≈ t :=
      Q.equiv_trans (tt_step_simplify t 37 37) hzero_add
    rw [hk, taiMinusUtc_congr (tt_step_simplify t 37 36), h0_val,
        taiMinusUtc_congr h1t, hk, taiMinusUtc_congr h1t, hk]
    exact h1t
  · by_cases c15 : T_2015_07_01 ≤ t
    · have hk : taiMinusUtcAstrosat t = 36 := by unfold taiMinusUtcAstrosat; simp [c17, c15]
      have h0t : t + Q.ofInt 36 + TT_MINUS_TAI - Q.ofInt 36 - TT_MINUS_TAI ≈ t :=
        Q.equiv_trans (tt_step_simplify t 36 36) hzero_add
      rw [hk, taiMinusUtc_congr h0t, hk, taiMinusUtc_congr h0t, hk, taiMinusUtc_congr h0t, hk]
      exact h0t
    · by_cases c12 : T_2012_07_01 ≤ t
      · have hk : taiMinusUtcAstrosat t = 35 := by unfold taiMinusUtcAstrosat; simp [c17, c15, c12]
        have ht12 : 78796800 * t.den ≤ t.num := by
          change 78796800 * t.den ≤ t.num * 1 at c12; omega
        have ht15 : t.num < 173404800 * t.den := by
          change ¬ (173404800 * t.den ≤ t.num * 1) at c15; omega
        have h1t : t + Q.ofInt 35 + TT_MINUS_TAI - Q.ofInt 35 - TT_MINUS_TAI ≈ t :=
          Q.equiv_trans (tt_step_simplify t 35 35) hzero_add
        by_cases hnear : 78796801 * t.den ≤ t.num
        · have h0_val : taiMinusUtcAstrosat (t + Q.ofInt (35 - 36)) = 35 := by
            unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
            have h17 : ¬ (Q.ofInt 220924800 ≤ t + Q.ofInt (-1)) := by rw [ofInt_le_add_ofInt]; omega
            have h15 : ¬ (Q.ofInt 173404800 ≤ t + Q.ofInt (-1)) := by rw [ofInt_le_add_ofInt]; omega
            have h12 : Q.ofInt 78796800 ≤ t + Q.ofInt (-1) := by rw [ofInt_le_add_ofInt]; omega
            simp [h17, h15, h12]
          rw [hk, taiMinusUtc_congr (tt_step_simplify t 35 36), h0_val,
              taiMinusUtc_congr h1t, hk, taiMinusUtc_congr h1t, hk]
          exact h1t
        · have h0_val : taiMinusUtcAstrosat (t + Q.ofInt (35 - 36)) = 34 := by
            unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
            have h17 : ¬ (Q.ofInt 220924800 ≤ t + Q.ofInt (-1)) := by rw [ofInt_le_add_ofInt]; omega
            have h15 : ¬ (Q.ofInt 173404800 ≤ t + Q.ofInt (-1)) := by rw [ofInt_le_add_ofInt]; omega
            have h12 : ¬ (Q.ofInt 78796800 ≤ t + Q.ofInt (-1)) := by rw [ofInt_le_add_ofInt]; omega
            simp [h17, h15, h12]
          have h1_val : taiMinusUtcAstrosat (t + Q.ofInt (35 - 34)) = 35 := by
            unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
            have h17 : ¬ (Q.ofInt 220924800 ≤ t + Q.ofInt 1) := by rw [ofInt_le_add_ofInt]; omega
            have h15 : ¬ (Q.ofInt 173404800 ≤ t + Q.ofInt 1) := by rw [ofInt_le_add_ofInt]; omega
            have h12 : Q.ofInt 78796800 ≤ t + Q.ofInt 1 := by rw [ofInt_le_add_ofInt]; omega
            simp [h17, h15, h12]
          rw [hk, taiMinusUtc_congr (tt_step_simplify t 35 36), h0_val,
              taiMinusUtc_congr (tt_step_simplify t 35 34), h1_val,
              taiMinusUtc_congr h1t, hk]
          exact h1t
      · have hk : taiMinusUtcAstrosat t = 34 := by unfold taiMinusUtcAstrosat; simp [c17, c15, c12]
        have ht12 : t.num < 78796800 * t.den := by
          change ¬ (78796800 * t.den ≤ t.num * 1) at c12; omega
        have h0_val : taiMinusUtcAstrosat (t + Q.ofInt (34 - 36)) = 34 := by
          unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
          have h17 : ¬ (Q.ofInt 220924800 ≤ t + Q.ofInt (-2)) := by rw [ofInt_le_add_ofInt]; omega
          have h15 : ¬ (Q.ofInt 173404800 ≤ t + Q.ofInt (-2)) := by rw [ofInt_le_add_ofInt]; omega
          have h12 : ¬ (Q.ofInt 78796800 ≤ t + Q.ofInt (-2)) := by rw [ofInt_le_add_ofInt]; omega
          simp [h17, h15, h12]
        have h1t : t + Q.ofInt 34 + TT_MINUS_TAI - Q.ofInt 34 - TT_MINUS_TAI ≈ t :=
          Q.equiv_trans (tt_step_simplify t 34 34) hzero_add
        rw [hk, taiMinusUtc_congr (tt_step_simplify t 34 36), h0_val,
            taiMinusUtc_congr h1t, hk, taiMinusUtc_congr h1t, hk]
        exact h1t

/-- Exact `UTC ↔ TAI` invertibility within each epoch: subtracting `taiMinusUtcAstrosat t`
    recovers `t` identically. -/
theorem tai_utc_invertible (t : Q) :
    astrosatUtcToTai t - Q.ofInt (taiMinusUtcAstrosat t) ≈ t :=
  Q.equiv_sub_cancel_add t (Q.ofInt (taiMinusUtcAstrosat t))

/-- Helper: exact difference of `astrosatUtcToTai` and `astrosatUtcToTt` across any 1-second leap boundary. -/
theorem leap_step_tai_tt (tb eps : Q) (c k : Int)
    (htb : tb = Q.ofInt c)
    (ha : taiMinusUtcAstrosat tb = k + 1)
    (hb : taiMinusUtcAstrosat (tb - eps) = k) :
    astrosatUtcToTai tb - astrosatUtcToTai (tb - eps) ≈ Q.one + eps ∧
    astrosatUtcToTt tb - astrosatUtcToTt (tb - eps) ≈ Q.one + eps := by
  subst htb
  constructor
  · unfold astrosatUtcToTai; rw [ha, hb]
    unfold Q.equiv Q.one Q.ofInt
    show ((c * 1 + (k + 1) * 1) * (1 * eps.den * 1) -
          ((c * eps.den - eps.num * 1) * 1 + k * (1 * eps.den)) * (1 * 1)) * (1 * eps.den) =
         (1 * eps.den + eps.num * 1) * (1 * 1 * (1 * eps.den * 1))
    simp only [Int.mul_one, Int.one_mul]
    have hlin : (c + (k + 1)) * eps.den - (c * eps.den - eps.num + k * eps.den) =
                eps.den + eps.num := by
      rw [Int.add_mul, Int.add_mul, Int.one_mul]; omega
    rw [hlin]
  · unfold astrosatUtcToTt; rw [ha, hb]
    unfold Q.equiv Q.one Q.ofInt TT_MINUS_TAI
    show (((c * 1 + (k + 1) * 1) * 1000 + 32184 * (1 * 1)) * (1 * eps.den * 1 * 1000) -
          (((c * eps.den - eps.num * 1) * 1 + k * (1 * eps.den)) * 1000 +
            32184 * (1 * eps.den * 1)) * (1 * 1 * 1000)) * (1 * eps.den) =
         (1 * eps.den + eps.num * 1) * (1 * 1 * 1000 * (1 * eps.den * 1 * 1000))
    simp only [Int.mul_one, Int.one_mul]
    have hlin : ((c + (k + 1)) * 1000 + 32184) * (eps.den * 1000) -
                ((c * eps.den - eps.num + k * eps.den) * 1000 + 32184 * eps.den) * 1000 =
                (eps.den + eps.num) * 1000000 := by
      rw [Int.add_mul, Int.add_mul, Int.add_mul, Int.add_mul, Int.add_mul, Int.add_mul,
          Int.add_mul, Int.add_mul, Int.sub_mul, Int.sub_mul, Int.one_mul]
      have e1 : c * 1000 * (eps.den * 1000) = c * eps.den * 1000000 := by simp only [Int.mul_assoc]; ac_rfl
      have e2 : k * 1000 * (eps.den * 1000) = k * eps.den * 1000000 := by simp only [Int.mul_assoc]; ac_rfl
      have e3 : 1000 * (eps.den * 1000) = eps.den * 1000000 := by omega
      have e4 : 32184 * (eps.den * 1000) = 32184 * eps.den * 1000 := (Int.mul_assoc ..).symm
      have e5 : c * eps.den * 1000 * 1000 = c * eps.den * 1000000 := by
        rw [Int.mul_assoc (c * eps.den) 1000 1000]; rfl
      have e6 : k * eps.den * 1000 * 1000 = k * eps.den * 1000000 := by
        rw [Int.mul_assoc (k * eps.den) 1000 1000]; rfl
      have e7 : eps.num * 1000 * 1000 = eps.num * 1000000 := by
        rw [Int.mul_assoc eps.num 1000 1000]; rfl
      rw [e1, e2, e3, e4, e5, e6, e7]; omega
    rw [hlin]
    have hden1 : (1000 : Int) * (eps.den * 1000) = 1000000 * eps.den := by omega
    rw [hden1]; simp only [Int.mul_assoc]

/-- Exact 1-second continuity preservation across the 2012-07-01 leap-second boundary (MJD 56109):
    for any `0 < ε ≤ 86400`, `TAI(T_2012) - TAI(T_2012 - ε) = 1 + ε` and
    `TT(T_2012) - TT(T_2012 - ε) = 1 + ε`. -/
theorem leap_second_2012_07_01_step (eps : Q) (hpos : Q.zero < eps) (_hle : eps ≤ Q.ofInt 86400) :
    taiMinusUtcAstrosat T_2012_07_01 - taiMinusUtcAstrosat (T_2012_07_01 - eps) = 1 ∧
    astrosatUtcToTai T_2012_07_01 - astrosatUtcToTai (T_2012_07_01 - eps) ≈ Q.one + eps ∧
    astrosatUtcToTt T_2012_07_01 - astrosatUtcToTt (T_2012_07_01 - eps) ≈ Q.one + eps := by
  have hep : 0 < eps.num := by change 0 * eps.den < eps.num * 1 at hpos; omega
  have hden := eps.den_pos
  have ha : taiMinusUtcAstrosat T_2012_07_01 = 35 := by decide
  have hb : taiMinusUtcAstrosat (T_2012_07_01 - eps) = 34 := by
    unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
    have h17 : ¬ (Q.ofInt 220924800 ≤ Q.ofInt 78796800 - eps) := by
      show ¬ (220924800 * (1 * eps.den) ≤ (78796800 * eps.den - eps.num * 1) * 1); omega
    have h15 : ¬ (Q.ofInt 173404800 ≤ Q.ofInt 78796800 - eps) := by
      show ¬ (173404800 * (1 * eps.den) ≤ (78796800 * eps.den - eps.num * 1) * 1); omega
    have h12 : ¬ (Q.ofInt 78796800 ≤ Q.ofInt 78796800 - eps) := by
      show ¬ (78796800 * (1 * eps.den) ≤ (78796800 * eps.den - eps.num * 1) * 1); omega
    simp [h17, h15, h12]
  have hpair := leap_step_tai_tt T_2012_07_01 eps 78796800 34 rfl ha hb
  exact ⟨by omega, hpair.1, hpair.2⟩

/-- Exact 1-second continuity preservation across the 2015-07-01 leap-second boundary (MJD 57204). -/
theorem leap_second_2015_07_01_step (eps : Q) (hpos : Q.zero < eps) (hle : eps ≤ Q.ofInt 86400) :
    taiMinusUtcAstrosat T_2015_07_01 - taiMinusUtcAstrosat (T_2015_07_01 - eps) = 1 ∧
    astrosatUtcToTai T_2015_07_01 - astrosatUtcToTai (T_2015_07_01 - eps) ≈ Q.one + eps ∧
    astrosatUtcToTt T_2015_07_01 - astrosatUtcToTt (T_2015_07_01 - eps) ≈ Q.one + eps := by
  have hep : 0 < eps.num := by change 0 * eps.den < eps.num * 1 at hpos; omega
  have hle' : eps.num ≤ 86400 * eps.den := by change eps.num * 1 ≤ 86400 * eps.den at hle; omega
  have hden := eps.den_pos
  have ha : taiMinusUtcAstrosat T_2015_07_01 = 36 := by decide
  have hb : taiMinusUtcAstrosat (T_2015_07_01 - eps) = 35 := by
    unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
    have h17 : ¬ (Q.ofInt 220924800 ≤ Q.ofInt 173404800 - eps) := by
      show ¬ (220924800 * (1 * eps.den) ≤ (173404800 * eps.den - eps.num * 1) * 1); omega
    have h15 : ¬ (Q.ofInt 173404800 ≤ Q.ofInt 173404800 - eps) := by
      show ¬ (173404800 * (1 * eps.den) ≤ (173404800 * eps.den - eps.num * 1) * 1); omega
    have h12 : Q.ofInt 78796800 ≤ Q.ofInt 173404800 - eps := by
      show 78796800 * (1 * eps.den) ≤ (173404800 * eps.den - eps.num * 1) * 1; omega
    simp [h17, h15, h12]
  have hpair := leap_step_tai_tt T_2015_07_01 eps 173404800 35 rfl ha hb
  exact ⟨by omega, hpair.1, hpair.2⟩

/-- Exact 1-second continuity preservation across the 2017-01-01 leap-second boundary (MJD 57754). -/
theorem leap_second_2017_01_01_step (eps : Q) (hpos : Q.zero < eps) (hle : eps ≤ Q.ofInt 86400) :
    taiMinusUtcAstrosat T_2017_01_01 - taiMinusUtcAstrosat (T_2017_01_01 - eps) = 1 ∧
    astrosatUtcToTai T_2017_01_01 - astrosatUtcToTai (T_2017_01_01 - eps) ≈ Q.one + eps ∧
    astrosatUtcToTt T_2017_01_01 - astrosatUtcToTt (T_2017_01_01 - eps) ≈ Q.one + eps := by
  have hep : 0 < eps.num := by change 0 * eps.den < eps.num * 1 at hpos; omega
  have hle' : eps.num ≤ 86400 * eps.den := by change eps.num * 1 ≤ 86400 * eps.den at hle; omega
  have hden := eps.den_pos
  have ha : taiMinusUtcAstrosat T_2017_01_01 = 37 := by decide
  have hb : taiMinusUtcAstrosat (T_2017_01_01 - eps) = 36 := by
    unfold taiMinusUtcAstrosat T_2017_01_01 T_2015_07_01 T_2012_07_01 ASTROSAT_MJDREFI SECS_PER_DAY
    have h17 : ¬ (Q.ofInt 220924800 ≤ Q.ofInt 220924800 - eps) := by
      show ¬ (220924800 * (1 * eps.den) ≤ (220924800 * eps.den - eps.num * 1) * 1); omega
    have h15 : Q.ofInt 173404800 ≤ Q.ofInt 220924800 - eps := by
      show 173404800 * (1 * eps.den) ≤ (220924800 * eps.den - eps.num * 1) * 1; omega
    simp [h17, h15]
  have hpair := leap_step_tai_tt T_2017_01_01 eps 220924800 36 rfl ha hb
  exact ⟨by omega, hpair.1, hpair.2⟩

/-! ## 4. CZTI 20-Bit Counter Rollover `mod 1048576` Preservation (`tct.rs:11,182-193,254-289`) -/

/-- Exact CZTI 20-bit counter wrap modulus `2^20 = 1048576` (`tct.rs:11`). -/
def CZTI_CLOCK_WRAP : Int := 1048576

/-- Exact mirror of the CZTI wrap-around step in `CztiTct::new` (`tct.rs:182-193`). -/
def cztiUnwrapStep (prevRaw raw : Int) (wrapOffset : Int) : Int × Int :=
  let dtInst := raw - prevRaw
  if dtInst < -(CZTI_CLOCK_WRAP / 2) then
    (raw + (wrapOffset + CZTI_CLOCK_WRAP), wrapOffset + CZTI_CLOCK_WRAP)
  else
    (raw + wrapOffset, wrapOffset)

/-- **CZTI 20-Bit Rollover `mod 1048576` Exactness Theorem** (`tct.rs:11,182-193,254-289`):
    1. Adding any integer multiple `k * 1048576` preserves Euclidean `mod 1048576`.
    2. For any consecutive monotonic clock readings `tPrev` and `tNext = tPrev + dt` with
       `0 < dt < 524288` (`0.5 * 2^20`) and wrap-offset invariant `wrapOffset % 1048576 = 0`,
       wrapping modulo `1048576` and applying `cztiUnwrapStep` recovers the EXACT elapsed
       interval `dt`, preserves `unwrapped % 1048576 = nextRaw`, and preserves the wrap-offset
       multiple-of-`1048576` invariant! -/
theorem czti_rollover_mod_invariant (t k : Int) :
    (t + k * CZTI_CLOCK_WRAP) % CZTI_CLOCK_WRAP = t % CZTI_CLOCK_WRAP := by
  unfold CZTI_CLOCK_WRAP
  omega

theorem czti_unwrapStep_exact (tPrev dt wrapOffset : Int)
    (hdt_pos : 0 < dt) (hdt_lt : dt < CZTI_CLOCK_WRAP / 2)
    (hwrap : wrapOffset % CZTI_CLOCK_WRAP = 0) :
    let prevRaw := tPrev % CZTI_CLOCK_WRAP
    let nextRaw := (tPrev + dt) % CZTI_CLOCK_WRAP
    let (unwrappedNext, newOffset) := cztiUnwrapStep prevRaw nextRaw wrapOffset
    let unwrappedPrev := prevRaw + wrapOffset
    unwrappedNext - unwrappedPrev = dt ∧
    unwrappedNext % CZTI_CLOCK_WRAP = nextRaw ∧
    newOffset % CZTI_CLOCK_WRAP = 0 ∧
    (newOffset = wrapOffset ∨ newOffset = wrapOffset + CZTI_CLOCK_WRAP) := by
  unfold cztiUnwrapStep CZTI_CLOCK_WRAP at *
  dsimp only
  split <;> refine ⟨by omega, by omega, by omega, by omega⟩

end AstrosatFormal.AstroPhysics
