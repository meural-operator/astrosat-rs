/-!
# AstrosatFormal.GtiAlgebra

Formal verification of AstroSat Good Time Interval (GTI) set algebra, mirroring
`crates/astrosat-core/src/gti.rs`:
- `GtiInterval`, `GtiSet::new` (normalization: filter valid, sort, merge overlaps/touching)
- `exposure` non-negativity
- `intersect` and `subtract` soundness & normalization preservation
- Monotonicity (`exposure(A ∩ B) ≤ exposure(A)`)
- Exact partition identity (`exposure(A ∩ B) + exposure(A \ B) = exposure(A)`)
- `drop_short` and `merge_gaps` invariants
-/

namespace AstrosatFormal.GtiAlgebra

/-! ## 1. Core Definitions (`gti.rs:13-62,149-180`) -/

/-- Closed/half-open Good Time Interval `[start, stop]` over exact integer ticks. -/
structure GtiInterval where
  start : Int
  stop  : Int
  deriving DecidableEq, Repr

/-- Duration of an interval (`(stop - start).max(0)`, `gti.rs:24-26`). -/
def GtiInterval.duration (iv : GtiInterval) : Int :=
  if iv.start < iv.stop then iv.stop - iv.start else 0

/-- Half-open membership `start ≤ t < stop` (measure-exact for interval partition & exposure). -/
def GtiInterval.memHalfOpen (iv : GtiInterval) (t : Int) : Prop :=
  iv.start ≤ t ∧ t < iv.stop

/-- Closed membership `start ≤ t ≤ stop` (`gti.rs:29-31`). -/
def GtiInterval.memClosed (iv : GtiInterval) (t : Int) : Prop :=
  iv.start ≤ t ∧ t ≤ iv.stop

/-- List-level half-open membership. -/
def memHalfOpen (t : Int) (L : List GtiInterval) : Prop :=
  ∃ iv ∈ L, iv.memHalfOpen t

/-- List-level closed membership. -/
def memClosed (t : Int) (L : List GtiInterval) : Prop :=
  ∃ iv ∈ L, iv.memClosed t

/-- Total exposure (sum of interval durations, `gti.rs:149-151`). -/
def exposure : List GtiInterval → Int
  | [] => 0
  | iv :: rest => iv.duration + exposure rest

theorem exposure_append (A B : List GtiInterval) :
    exposure (A ++ B) = exposure A + exposure B := by
  induction A with
  | nil => simp [exposure]
  | cons a as ih => simp [exposure, ih]; omega

theorem duration_nonneg (iv : GtiInterval) : 0 ≤ iv.duration := by
  unfold GtiInterval.duration
  split <;> omega

theorem exposure_nonneg (L : List GtiInterval) : 0 ≤ exposure L := by
  induction L with
  | nil => simp [exposure]
  | cons a as ih =>
    simp only [exposure]
    have ha := duration_nonneg a
    omega

/-! ## 2. Normalization Invariant (`gti.rs:42-62`) -/

/-- `NormalizedAfter b L`: every interval in `L` satisfies `start < stop`,
    is strictly sorted & disjoint (`iv₁.stop < iv₂.start`), and starts strictly after `b`. -/
def NormalizedAfter (b : Int) : List GtiInterval → Prop
  | [] => True
  | iv :: rest => b < iv.start ∧ iv.start < iv.stop ∧ NormalizedAfter iv.stop rest

/-- `Normalized L`: every interval has `start < stop`, and adjacent intervals are
    strictly separated (`iv₁.stop < iv₂.start`, no overlaps or touching endpoints). -/
def Normalized : List GtiInterval → Prop
  | [] => True
  | iv :: rest => iv.start < iv.stop ∧ NormalizedAfter iv.stop rest

theorem NormalizedAfter.mono {b₁ b₂ : Int} {L : List GtiInterval}
    (h : NormalizedAfter b₂ L) (hle : b₁ ≤ b₂) : NormalizedAfter b₁ L := by
  cases L with
  | nil => trivial
  | cons iv rest =>
    rcases h with ⟨h1, h2, h3⟩
    exact ⟨by omega, h2, h3⟩

theorem NormalizedAfter.toNormalized {b : Int} {L : List GtiInterval}
    (h : NormalizedAfter b L) : Normalized L := by
  cases L with
  | nil => trivial
  | cons iv rest =>
    rcases h with ⟨_, h2, h3⟩
    exact ⟨h2, h3⟩

theorem NormalizedAfter.forall_start_gt {b : Int} {L : List GtiInterval}
    (h : NormalizedAfter b L) : ∀ iv ∈ L, b < iv.start ∧ iv.start < iv.stop := by
  induction L generalizing b with
  | nil => intro iv hmem; cases hmem
  | cons a as ih =>
    rcases h with ⟨hb, ha, hrest⟩
    intro iv hmem
    cases hmem with
    | head => exact ⟨hb, ha⟩
    | tail _ hm =>
      have hsub := ih hrest iv hm
      exact ⟨by omega, hsub.2⟩

theorem Normalized.forall_valid {L : List GtiInterval} (h : Normalized L) :
    ∀ iv ∈ L, iv.start < iv.stop := by
  cases L with
  | nil => intro iv hm; cases hm
  | cons a as =>
    rcases h with ⟨ha, hrest⟩
    intro iv hm
    cases hm with
    | head => exact ha
    | tail _ hmem => exact (hrest.forall_start_gt iv hmem).2

theorem Normalized.tail {a : GtiInterval} {as : List GtiInterval}
    (h : Normalized (a :: as)) : Normalized as :=
  h.2.toNormalized

/-- Lexicographic order on intervals `(start, stop)` (`gti.rs:47`). -/
def ivLe (a b : GtiInterval) : Bool :=
  a.start < b.start || (a.start == b.start && a.stop ≤ b.stop)

/-- Insertion sort by `(start, stop)`. -/
def insertSorted (x : GtiInterval) : List GtiInterval → List GtiInterval
  | [] => [x]
  | y :: ys => if ivLe x y then x :: y :: ys else y :: insertSorted x ys

def sortIntervals : List GtiInterval → List GtiInterval
  | [] => []
  | x :: xs => insertSorted x (sortIntervals xs)

/-- Sorted-by-start and valid predicate for pre-merge interval lists. -/
def ValidSortedFrom (s : Int) : List GtiInterval → Prop
  | [] => True
  | iv :: rest => s ≤ iv.start ∧ iv.start < iv.stop ∧ ValidSortedFrom iv.start rest

def ValidSorted : List GtiInterval → Prop
  | [] => True
  | iv :: rest => iv.start < iv.stop ∧ ValidSortedFrom iv.start rest

theorem ValidSortedFrom.mono {s₁ s₂ : Int} {L : List GtiInterval}
    (h : ValidSortedFrom s₂ L) (hle : s₁ ≤ s₂) : ValidSortedFrom s₁ L := by
  cases L with
  | nil => trivial
  | cons iv rest =>
    rcases h with ⟨h1, h2, h3⟩
    exact ⟨by omega, h2, h3⟩

theorem insertSorted_validSortedFrom {s : Int} {x : GtiInterval} {L : List GtiInterval}
    (hx : x.start < x.stop) (hxs : s ≤ x.start) (hL : ValidSortedFrom s L) :
    ValidSortedFrom s (insertSorted x L) := by
  induction L generalizing s with
  | nil => exact ⟨hxs, hx, trivial⟩
  | cons y ys ih =>
    rcases hL with ⟨hsy, hy, hrest⟩
    unfold insertSorted
    by_cases hle : ivLe x y = true
    · simp only [hle, ↓reduceIte]
      have hxy : x.start ≤ y.start := by
        unfold ivLe at hle
        simp only [Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hle
        omega
      exact ⟨hxs, hx, hxy, hy, hrest⟩
    · simp only [hle, Bool.false_eq_true, ↓reduceIte]
      have hyx : y.start ≤ x.start := by
        unfold ivLe at hle
        simp only [Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq, not_or] at hle
        omega
      exact ⟨hsy, hy, ih hyx hrest⟩

theorem insertSorted_validSorted {x : GtiInterval} {L : List GtiInterval}
    (hx : x.start < x.stop) (hL : ValidSorted L) : ValidSorted (insertSorted x L) := by
  cases L with
  | nil => exact ⟨hx, trivial⟩
  | cons y ys =>
    rcases hL with ⟨hy, hrest⟩
    unfold insertSorted
    by_cases hle : ivLe x y = true
    · simp only [hle, ↓reduceIte]
      have hxy : x.start ≤ y.start := by
        unfold ivLe at hle
        simp only [Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at hle
        omega
      exact ⟨hx, hxy, hy, hrest⟩
    · simp only [hle, Bool.false_eq_true, ↓reduceIte]
      have hyx : y.start ≤ x.start := by
        unfold ivLe at hle
        simp only [Bool.or_eq_true, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq, not_or] at hle
        omega
      exact ⟨hy, insertSorted_validSortedFrom hx hyx hrest⟩

theorem sortIntervals_validSorted {L : List GtiInterval}
    (hvalid : ∀ iv ∈ L, iv.start < iv.stop) : ValidSorted (sortIntervals L) := by
  induction L with
  | nil => trivial
  | cons x xs ih =>
    apply insertSorted_validSorted
    · exact hvalid x (List.mem_cons_self ..)
    · exact ih (fun iv hm => hvalid iv (List.mem_cons_of_mem _ hm))

/-- Merge overlapping or touching sorted intervals (`gti.rs:49-60`). -/
def mergeSortedAux (cur : GtiInterval) : List GtiInterval → List GtiInterval
  | [] => [cur]
  | iv :: rest =>
    if iv.start ≤ cur.stop then
      mergeSortedAux ⟨cur.start, max cur.stop iv.stop⟩ rest
    else
      cur :: mergeSortedAux iv rest

def mergeSorted : List GtiInterval → List GtiInterval
  | [] => []
  | iv :: rest => mergeSortedAux iv rest

/-- Exact mirror of `GtiSet::new` (`gti.rs:42-62`): filter `start < stop`, sort, merge overlaps. -/
def normalize (L : List GtiInterval) : List GtiInterval :=
  mergeSorted (sortIntervals (L.filter (fun iv => decide (iv.start < iv.stop))))

theorem mergeSortedAux_normalizedAfter {b : Int} {cur : GtiInterval} {rest : List GtiInterval}
    (hb : b < cur.start) (hcur : cur.start < cur.stop)
    (hrest : ValidSortedFrom cur.start rest) :
    NormalizedAfter b (mergeSortedAux cur rest) := by
  induction rest generalizing b cur with
  | nil => exact ⟨hb, hcur, trivial⟩
  | cons iv ys ih =>
    rcases hrest with ⟨hcs, hiv, hys⟩
    unfold mergeSortedAux
    by_cases hov : iv.start ≤ cur.stop
    · simp only [hov, ↓reduceIte]
      exact ih (cur := ⟨cur.start, max cur.stop iv.stop⟩) hb (by dsimp; omega) (hys.mono hcs)
    · simp only [hov, ↓reduceIte]
      exact ⟨hb, hcur, ih (by omega) hiv hys⟩

theorem mergeSorted_normalized {L : List GtiInterval} (h : ValidSorted L) :
    Normalized (mergeSorted L) := by
  cases L with
  | nil => trivial
  | cons cur rest =>
    rcases h with ⟨hcur, hrest⟩
    exact (mergeSortedAux_normalizedAfter (b := cur.start - 1) (by omega) hcur hrest).toNormalized

/-- Normalization invariant: `normalize L` is ALWAYS strictly valid, sorted, and disjoint (`gti.rs:42-62`). -/
theorem normalize_normalized (L : List GtiInterval) : Normalized (normalize L) := by
  unfold normalize
  apply mergeSorted_normalized
  apply sortIntervals_validSorted
  intro iv hmem
  rw [List.mem_filter] at hmem
  exact of_decide_eq_true hmem.2

/-- Membership helper for `insertSorted`. -/
theorem mem_insertSorted {x iv : GtiInterval} {L : List GtiInterval} :
    iv ∈ insertSorted x L ↔ iv = x ∨ iv ∈ L := by
  induction L with
  | nil => simp [insertSorted]
  | cons y ys ih =>
    unfold insertSorted
    split
    · simp
    · simp [ih, or_left_comm]

theorem mem_sortIntervals {iv : GtiInterval} {L : List GtiInterval} :
    iv ∈ sortIntervals L ↔ iv ∈ L := by
  induction L with
  | nil => simp [sortIntervals]
  | cons x xs ih => simp [sortIntervals, mem_insertSorted, ih]

theorem memHalfOpen_mergeSortedAux (t : Int) (cur : GtiInterval) (rest : List GtiInterval)
    (hcur : cur.start < cur.stop) (hrest : ValidSortedFrom cur.start rest) :
    memHalfOpen t (mergeSortedAux cur rest) ↔
      cur.memHalfOpen t ∨ memHalfOpen t rest := by
  induction rest generalizing cur with
  | nil =>
    simp [mergeSortedAux, memHalfOpen]
  | cons iv ys ih =>
    rcases hrest with ⟨hcs, hiv, hys⟩
    unfold mergeSortedAux
    by_cases hov : iv.start ≤ cur.stop
    · simp only [hov, ↓reduceIte]
      rw [ih ⟨cur.start, max cur.stop iv.stop⟩ (by dsimp; omega) (hys.mono hcs)]
      simp only [memHalfOpen, GtiInterval.memHalfOpen, List.mem_cons]
      constructor
      · intro h
        rcases h with ⟨h1, h2⟩ | ⟨w, hw, hwt⟩
        · by_cases hlt : t < cur.stop
          · exact Or.inl ⟨h1, hlt⟩
          · exact Or.inr ⟨iv, Or.inl rfl, by omega, by omega⟩
        · exact Or.inr ⟨w, Or.inr hw, hwt⟩
      · intro h
        rcases h with ⟨h1, h2⟩ | ⟨w, hw, hwt⟩
        · exact Or.inl ⟨h1, by omega⟩
        · rcases hw with rfl | hw'
          · exact Or.inl ⟨by omega, by omega⟩
          · exact Or.inr ⟨w, hw', hwt⟩
    · simp only [hov, ↓reduceIte]
      have hrec := ih iv hiv hys
      simp only [memHalfOpen, List.mem_cons] at hrec ⊢
      constructor
      · intro ⟨w, hw, hwt⟩
        rcases hw with rfl | hw'
        · exact Or.inl hwt
        · have hsub := hrec.mp ⟨w, hw', hwt⟩
          rcases hsub with h1 | ⟨u, hu, hut⟩
          · exact Or.inr ⟨iv, Or.inl rfl, h1⟩
          · exact Or.inr ⟨u, Or.inr hu, hut⟩
      · intro h
        rcases h with h1 | ⟨w, hw, hwt⟩
        · exact ⟨cur, Or.inl rfl, h1⟩
        · rcases hw with rfl | hw'
          · rcases hrec.mpr (Or.inl hwt) with ⟨u, hu, hut⟩
            exact ⟨u, Or.inr hu, hut⟩
          · rcases hrec.mpr (Or.inr ⟨w, hw', hwt⟩) with ⟨u, hu, hut⟩
            exact ⟨u, Or.inr hu, hut⟩

theorem memHalfOpen_mergeSorted (t : Int) {S : List GtiInterval} (hvs : ValidSorted S) :
    memHalfOpen t (mergeSorted S) ↔ memHalfOpen t S := by
  cases S with
  | nil => simp [mergeSorted, memHalfOpen]
  | cons c cs =>
    simp only [mergeSorted]
    rw [memHalfOpen_mergeSortedAux t c cs hvs.1 hvs.2]
    simp [memHalfOpen]

/-- Normalization preserves half-open membership of all valid intervals (`gti.rs:42-62`). -/
theorem memHalfOpen_normalize (t : Int) (L : List GtiInterval) :
    memHalfOpen t (normalize L) ↔ ∃ iv ∈ L, iv.start < iv.stop ∧ iv.memHalfOpen t := by
  unfold normalize
  have hvs : ValidSorted (sortIntervals (L.filter (fun iv => decide (iv.start < iv.stop)))) := by
    apply sortIntervals_validSorted
    intro iv hm
    rw [List.mem_filter] at hm
    exact of_decide_eq_true hm.2
  rw [memHalfOpen_mergeSorted t hvs]
  simp only [memHalfOpen, mem_sortIntervals, List.mem_filter, decide_eq_true_eq]
  constructor
  · intro ⟨iv, ⟨hm, hv⟩, hmem⟩
    exact ⟨iv, hm, hv, hmem⟩
  · intro ⟨iv, hm, hv, hmem⟩
    exact ⟨iv, ⟨hm, hv⟩, hmem⟩

/-! ## 3. Exact Set Intersection (`gti.rs:183-201`) -/

/-- Exact O(N + M) two-pointer intersection mirroring `GtiSet::intersect` (`gti.rs:183-201`). -/
def intersect : List GtiInterval → List GtiInterval → List GtiInterval
  | [], _ => []
  | _, [] => []
  | a :: as, b :: bs =>
    let lo := max a.start b.start
    let hi := min a.stop b.stop
    let rest :=
      if a.stop < b.stop then intersect as (b :: bs)
      else intersect (a :: as) bs
    if lo < hi then ⟨lo, hi⟩ :: rest else rest
termination_by a b => a.length + b.length

/-- Strong preservation of `NormalizedAfter` under `intersect`. -/
theorem intersect_normalizedAfter {x y : Int} (A B : List GtiInterval)
    (hA : NormalizedAfter x A) (hB : NormalizedAfter y B) :
    NormalizedAfter (max x y) (intersect A B) := by
  match A, B with
  | [], _ => simp [intersect, NormalizedAfter]
  | _ :: _, [] => simp [intersect, NormalizedAfter]
  | a :: as, b :: bs =>
    rcases hA with ⟨hxa, ha, has⟩
    rcases hB with ⟨hyb, hb, hbs⟩
    unfold intersect
    dsimp only
    by_cases hcmp : a.stop < b.stop
    · simp only [hcmp, ↓reduceIte]
      have hBcons : NormalizedAfter y (b :: bs) := ⟨hyb, hb, hbs⟩
      have hrec := intersect_normalizedAfter as (b :: bs) has hBcons
      by_cases hov : max a.start b.start < min a.stop b.stop
      · simp only [hov, ↓reduceIte]
        exact ⟨by dsimp; omega, hov, hrec.mono (by dsimp; omega)⟩
      · simp only [hov, ↓reduceIte]
        exact hrec.mono (by omega)
    · simp only [hcmp, ↓reduceIte]
      have hAcons : NormalizedAfter x (a :: as) := ⟨hxa, ha, has⟩
      have hrec := intersect_normalizedAfter (a :: as) bs hAcons hbs
      by_cases hov : max a.start b.start < min a.stop b.stop
      · simp only [hov, ↓reduceIte]
        exact ⟨by dsimp; omega, hov, hrec.mono (by dsimp; omega)⟩
      · simp only [hov, ↓reduceIte]
        exact hrec.mono (by omega)
termination_by A.length + B.length

/-- `intersect` preserves `Normalized` (`gti.rs:183-201`). -/
theorem intersect_normalized {A B : List GtiInterval}
    (hA : Normalized A) (hB : Normalized B) : Normalized (intersect A B) := by
  match A, B with
  | [], _ => simp [intersect, Normalized]
  | _ :: _, [] => simp [intersect, Normalized]
  | a :: as, b :: bs =>
    have hA' : NormalizedAfter (a.start - 1) (a :: as) := ⟨by omega, hA.1, hA.2⟩
    have hB' : NormalizedAfter (b.start - 1) (b :: bs) := ⟨by omega, hB.1, hB.2⟩
    exact (intersect_normalizedAfter (a :: as) (b :: bs) hA' hB').toNormalized

/-- Exact half-open membership soundness of `intersect`: `t ∈ A ∩ B ↔ t ∈ A ∧ t ∈ B`. -/
theorem memHalfOpen_intersect (t : Int) (A B : List GtiInterval)
    (hA : Normalized A) (hB : Normalized B) :
    memHalfOpen t (intersect A B) ↔ memHalfOpen t A ∧ memHalfOpen t B := by
  match A, B with
  | [], _ => simp [intersect, memHalfOpen]
  | _ :: _, [] => simp [intersect, memHalfOpen]
  | a :: as, b :: bs =>
    rcases hA with ⟨ha, has⟩
    rcases hB with ⟨hb, hbs⟩
    unfold intersect
    dsimp only
    by_cases hcmp : a.stop < b.stop
    · simp only [hcmp, ↓reduceIte]
      have hBcons : Normalized (b :: bs) := ⟨hb, hbs⟩
      have ih := memHalfOpen_intersect t as (b :: bs) has.toNormalized hBcons
      by_cases hov : max a.start b.start < min a.stop b.stop
      · simp only [hov, ↓reduceIte]
        simp only [memHalfOpen, List.mem_cons, exists_eq_or_imp] at ih ⊢
        rw [ih]
        constructor
        · intro h
          rcases h with hhead | ⟨⟨u, hu, hut⟩, hBmem⟩
          · unfold GtiInterval.memHalfOpen at hhead; dsimp at hhead
            refine ⟨Or.inl ⟨by omega, by omega⟩, Or.inl ⟨by omega, by omega⟩⟩
          · exact ⟨Or.inr ⟨u, hu, hut⟩, hBmem⟩
        · intro ⟨hAmem, hBmem⟩
          rcases hAmem with hat | ⟨u, hu, hut⟩
          · rcases hBmem with hbt | ⟨v, hv, hvt⟩
            · unfold GtiInterval.memHalfOpen at hat hbt ⊢; dsimp
              exact Or.inl ⟨by omega, by omega⟩
            · have hvgt := (hbs.forall_start_gt v hv).1
              exfalso; unfold GtiInterval.memHalfOpen at hat hvt; omega
          · exact Or.inr ⟨⟨u, hu, hut⟩, hBmem⟩
      · simp only [hov, ↓reduceIte]
        rw [ih]
        simp only [memHalfOpen, List.mem_cons, exists_eq_or_imp]
        constructor
        · intro ⟨hAs, hBmem⟩
          exact ⟨Or.inr hAs, hBmem⟩
        · intro ⟨hAmem, hBmem⟩
          rcases hAmem with hat | hAs
          · exfalso
            rcases hBmem with hbt | ⟨v, hv, hvt⟩
            · unfold GtiInterval.memHalfOpen at hat hbt; omega
            · have hvgt := (hbs.forall_start_gt v hv).1
              unfold GtiInterval.memHalfOpen at hat hvt; omega
          · exact ⟨hAs, hBmem⟩
    · simp only [hcmp, ↓reduceIte]
      have hAcons : Normalized (a :: as) := ⟨ha, has⟩
      have ih := memHalfOpen_intersect t (a :: as) bs hAcons hbs.toNormalized
      by_cases hov : max a.start b.start < min a.stop b.stop
      · simp only [hov, ↓reduceIte]
        simp only [memHalfOpen, List.mem_cons, exists_eq_or_imp] at ih ⊢
        rw [ih]
        constructor
        · intro h
          rcases h with hhead | ⟨hAmem, ⟨v, hv, hvt⟩⟩
          · unfold GtiInterval.memHalfOpen at hhead; dsimp at hhead
            refine ⟨Or.inl ⟨by omega, by omega⟩, Or.inl ⟨by omega, by omega⟩⟩
          · exact ⟨hAmem, Or.inr ⟨v, hv, hvt⟩⟩
        · intro ⟨hAmem, hBmem⟩
          rcases hBmem with hbt | ⟨v, hv, hvt⟩
          · rcases hAmem with hat | ⟨u, hu, hut⟩
            · unfold GtiInterval.memHalfOpen at hat hbt ⊢; dsimp
              exact Or.inl ⟨by omega, by omega⟩
            · have hugt := (has.forall_start_gt u hu).1
              exfalso; unfold GtiInterval.memHalfOpen at hbt hut; omega
          · exact Or.inr ⟨hAmem, ⟨v, hv, hvt⟩⟩
      · simp only [hov, ↓reduceIte]
        rw [ih]
        simp only [memHalfOpen, List.mem_cons, exists_eq_or_imp]
        constructor
        · intro ⟨hAmem, hBs⟩
          exact ⟨hAmem, Or.inr hBs⟩
        · intro ⟨hAmem, hBmem⟩
          rcases hBmem with hbt | hBs
          · exfalso
            rcases hAmem with hat | ⟨u, hu, hut⟩
            · unfold GtiInterval.memHalfOpen at hat hbt; omega
            · have hugt := (has.forall_start_gt u hu).1
              unfold GtiInterval.memHalfOpen at hbt hut; omega
          · exact ⟨hAmem, hBs⟩
termination_by A.length + B.length

/-! ## 4. Exact Set Subtraction (`gti.rs:212-245`) & Exact Partition Identity -/

/-- Subtract normalized `B` from a single interval `[curStart, curStop]`,
    mirroring the inner loop of `GtiSet::subtract` (`gti.rs:221-241`). -/
def subtractOne (curStart curStop : Int) : List GtiInterval → List GtiInterval
  | [] => if curStart < curStop then [⟨curStart, curStop⟩] else []
  | b :: bs =>
    if b.stop ≤ curStart then
      subtractOne curStart curStop bs
    else if curStop ≤ b.start then
      if curStart < curStop then [⟨curStart, curStop⟩] else []
    else
      let head := if curStart < b.start then [⟨curStart, b.start⟩] else []
      let nextStart := max curStart b.stop
      if nextStart < curStop then
        head ++ subtractOne nextStart curStop bs
      else
        head

/-- Exact set difference `A \ B` mirroring `GtiSet::subtract` (`gti.rs:212-245`). -/
def subtract (A B : List GtiInterval) : List GtiInterval :=
  match A with
  | [] => []
  | a :: as => subtractOne a.start a.stop B ++ subtract as B

/-- Single-interval intersection helper used to bridge `intersect` and `subtractOne`. -/
def intersectOne (curStart curStop : Int) : List GtiInterval → List GtiInterval
  | [] => []
  | b :: bs =>
    if b.stop ≤ curStart then
      intersectOne curStart curStop bs
    else if curStop ≤ b.start then
      []
    else
      let lo := max curStart b.start
      let hi := min curStop b.stop
      let head := if lo < hi then [⟨lo, hi⟩] else []
      head ++ intersectOne (max curStart b.stop) curStop bs

theorem subtractOne_bounds (lo s e : Int) (B : List GtiInterval)
    (hlo : lo < s) (hB : Normalized B) :
    NormalizedAfter lo (subtractOne s e B) ∧
    (∀ iv ∈ subtractOne s e B, iv.stop ≤ e) := by
  induction B generalizing lo s with
  | nil =>
    unfold subtractOne
    split
    · refine ⟨⟨hlo, by assumption, trivial⟩, ?_⟩
      intro iv hm; simp at hm; subst hm; dsimp; omega
    · exact ⟨trivial, fun _ hm => by cases hm⟩
  | cons b bs ih =>
    rcases hB with ⟨hb, hbs⟩
    unfold subtractOne
    by_cases h1 : b.stop ≤ s
    · simp only [h1, ↓reduceIte]
      exact ih lo s hlo hbs.toNormalized
    · simp only [h1, ↓reduceIte]
      by_cases h2 : e ≤ b.start
      · simp only [h2, ↓reduceIte]
        split
        · refine ⟨⟨hlo, by assumption, trivial⟩, ?_⟩
          intro iv hm; simp at hm; subst hm; dsimp; omega
        · exact ⟨trivial, fun _ hm => by cases hm⟩
      · simp only [h2, ↓reduceIte]
        by_cases h3 : max s b.stop < e
        · simp only [h3, ↓reduceIte]
          by_cases h4 : s < b.start
          · simp only [h4, ↓reduceIte, List.singleton_append]
            have hrec2 := ih b.start (max s b.stop) (by omega) hbs.toNormalized
            refine ⟨⟨hlo, h4, hrec2.1⟩, ?_⟩
            intro iv hm
            cases hm with
            | head => dsimp; omega
            | tail _ hmem => exact hrec2.2 iv hmem
          · simp only [h4, ↓reduceIte, List.nil_append]
            exact ih lo (max s b.stop) (by omega) hbs.toNormalized
        · simp only [h3, ↓reduceIte]
          by_cases h4 : s < b.start
          · simp only [h4, ↓reduceIte]
            refine ⟨⟨hlo, h4, trivial⟩, ?_⟩
            intro iv hm; simp at hm; subst hm; dsimp; omega
          · simp only [h4, ↓reduceIte]
            exact ⟨trivial, fun _ hm => by cases hm⟩

theorem normalizedAfter_append {b mid : Int} {X Y : List GtiInterval}
    (hX : NormalizedAfter b X) (hXub : ∀ iv ∈ X, iv.stop ≤ mid)
    (hY : NormalizedAfter mid Y) (hbmid : b ≤ mid) :
    NormalizedAfter b (X ++ Y) := by
  induction X generalizing b with
  | nil => exact hY.mono hbmid
  | cons x xs ih =>
    rcases hX with ⟨hbx, hx, hxs⟩
    have hxub : x.stop ≤ mid := hXub x (List.mem_cons_self ..)
    have hxsub : ∀ iv ∈ xs, iv.stop ≤ mid := fun iv hm => hXub iv (List.mem_cons_of_mem _ hm)
    exact ⟨hbx, hx, ih hxs hxsub hxub⟩

theorem subtract_normalizedAfter {lo : Int} {A B : List GtiInterval}
    (hA : NormalizedAfter lo A) (hB : Normalized B) :
    NormalizedAfter lo (subtract A B) := by
  induction A generalizing lo with
  | nil => trivial
  | cons a as ih =>
    rcases hA with ⟨hlo, ha, has⟩
    unfold subtract
    have h1 := subtractOne_bounds lo a.start a.stop B hlo hB
    have h2 := ih has
    exact normalizedAfter_append h1.1 h1.2 h2 (by omega)

/-- `subtract` preserves `Normalized` (`gti.rs:212-245`). -/
theorem subtract_normalized {A B : List GtiInterval}
    (hA : Normalized A) (hB : Normalized B) : Normalized (subtract A B) := by
  cases A with
  | nil => trivial
  | cons a as =>
    have hA' : NormalizedAfter (a.start - 1) (a :: as) := ⟨by omega, hA.1, hA.2⟩
    exact (subtract_normalizedAfter hA' hB).toNormalized

/-- Soundness of `subtractOne`: `t ∈ subtractOne s e B ↔ (s ≤ t < e) ∧ t ∉ B`. -/
theorem memHalfOpen_subtractOne (t s e : Int) {B : List GtiInterval} (hB : Normalized B) :
    memHalfOpen t (subtractOne s e B) ↔ (s ≤ t ∧ t < e) ∧ ¬ memHalfOpen t B := by
  induction B generalizing s with
  | nil =>
    unfold subtractOne
    by_cases hse : s < e
    · simp only [hse, ↓reduceIte, memHalfOpen, GtiInterval.memHalfOpen, List.mem_singleton]
      constructor
      · intro ⟨iv, hEq, ht⟩; subst hEq; dsimp at ht; exact ⟨ht, fun ⟨_, hm, _⟩ => by cases hm⟩
      · intro ⟨ht, _⟩; exact ⟨⟨s, e⟩, rfl, ht⟩
    · simp only [hse, ↓reduceIte, memHalfOpen, List.not_mem_nil]
      constructor
      · intro ⟨_, hm, _⟩; cases hm
      · intro ⟨ht, _⟩; omega
  | cons b bs ih =>
    rcases hB with ⟨hb, hbs⟩
    unfold subtractOne
    by_cases h1 : b.stop ≤ s
    · simp only [h1, ↓reduceIte]
      rw [ih s hbs.toNormalized]
      simp only [memHalfOpen, List.mem_cons, GtiInterval.memHalfOpen]
      constructor
      · intro ⟨ht, hnot⟩
        refine ⟨ht, fun ⟨iv, hm, hiv⟩ => ?_⟩
        rcases hm with rfl | hm'
        · omega
        · exact hnot ⟨iv, hm', hiv⟩
      · intro ⟨ht, hnot⟩
        exact ⟨ht, fun ⟨iv, hm', hiv⟩ => hnot ⟨iv, Or.inr hm', hiv⟩⟩
    · simp only [h1, ↓reduceIte]
      by_cases h2 : e ≤ b.start
      · simp only [h2, ↓reduceIte]
        have hgt : ∀ iv ∈ bs, b.stop < iv.start := fun iv hm => (hbs.forall_start_gt iv hm).1
        by_cases hse : s < e
        · simp only [hse, ↓reduceIte, memHalfOpen, List.mem_singleton, GtiInterval.memHalfOpen]
          constructor
          · intro ⟨iv, hEq, ht⟩; subst hEq; dsimp at ht
            refine ⟨ht, fun ⟨w, hw, hwt⟩ => ?_⟩
            rcases List.mem_cons.mp hw with rfl | hw'
            · omega
            · have := hgt w hw'; omega
          · intro ⟨ht, _⟩; exact ⟨⟨s, e⟩, rfl, ht⟩
        · simp only [hse, ↓reduceIte, memHalfOpen, List.not_mem_nil]
          constructor
          · intro ⟨_, hm, _⟩; cases hm
          · intro ⟨ht, _⟩; omega
      · simp only [h2, ↓reduceIte]
        have hrec := ih (max s b.stop) hbs.toNormalized
        have hgt : ∀ iv ∈ bs, b.stop < iv.start := fun iv hm => (hbs.forall_start_gt iv hm).1
        by_cases h3 : max s b.stop < e
        · simp only [h3, ↓reduceIte]
          by_cases h4 : s < b.start
          · simp only [h4, ↓reduceIte, List.singleton_append]
            simp only [memHalfOpen, List.mem_cons, exists_eq_or_imp] at hrec ⊢
            rw [hrec]
            simp only [GtiInterval.memHalfOpen]
            constructor
            · intro h
              rcases h with hhead | ⟨ht, hnot⟩
              · refine ⟨⟨by omega, by omega⟩, fun hmem => ?_⟩
                rcases hmem with hbt | ⟨w, hw', hwt⟩
                · omega
                · have := hgt w hw'; omega
              · refine ⟨⟨by omega, by omega⟩, fun hmem => ?_⟩
                rcases hmem with hbt | ⟨w, hw', hwt⟩
                · omega
                · exact hnot ⟨w, hw', hwt⟩
            · intro ⟨ht, hnot⟩
              by_cases htlt : t < b.start
              · exact Or.inl ⟨ht.1, htlt⟩
              · have hge : b.stop ≤ t := by
                  by_cases hlt : t < b.stop
                  · exact absurd (Or.inl ⟨by omega, hlt⟩) hnot
                  · omega
                exact Or.inr ⟨⟨by omega, ht.2⟩, fun ⟨w, hw', hwt⟩ => hnot (Or.inr ⟨w, hw', hwt⟩)⟩
          · simp only [h4, ↓reduceIte, List.nil_append]
            rw [hrec]
            simp only [memHalfOpen, List.mem_cons, exists_eq_or_imp, GtiInterval.memHalfOpen]
            constructor
            · intro ⟨ht, hnot⟩
              refine ⟨⟨by omega, ht.2⟩, fun hmem => ?_⟩
              rcases hmem with hbt | ⟨w, hw', hwt⟩
              · omega
              · exact hnot ⟨w, hw', hwt⟩
            · intro ⟨ht, hnot⟩
              have hge : b.stop ≤ t := by
                by_cases hlt : t < b.stop
                · exact absurd (Or.inl ⟨by omega, hlt⟩) hnot
                · omega
              exact ⟨⟨by omega, ht.2⟩, fun ⟨w, hw', hwt⟩ => hnot (Or.inr ⟨w, hw', hwt⟩)⟩
        · simp only [h3, ↓reduceIte]
          by_cases h4 : s < b.start
          · simp only [h4, ↓reduceIte, memHalfOpen, List.mem_singleton, GtiInterval.memHalfOpen]
            constructor
            · intro ⟨iv, hEq, ht⟩; subst hEq; dsimp at ht
              refine ⟨⟨ht.1, by omega⟩, fun ⟨w, hw, hwt⟩ => ?_⟩
              rcases List.mem_cons.mp hw with rfl | hw'
              · omega
              · have := hgt w hw'; omega
            · intro ⟨ht, hnot⟩
              have htlt : t < b.start := by
                by_cases hlt : t < b.start
                · exact hlt
                · exfalso
                  exact hnot ⟨b, List.mem_cons_self .., by omega, by omega⟩
              exact ⟨⟨s, b.start⟩, rfl, ht.1, htlt⟩
          · simp only [h4, ↓reduceIte, memHalfOpen, List.not_mem_nil]
            constructor
            · intro ⟨_, hm, _⟩; cases hm
            · intro ⟨ht, hnot⟩
              exact False.elim (hnot ⟨b, List.mem_cons_self .., by omega, by omega⟩)

/-- Exact soundness of `subtract`: `t ∈ A \ B ↔ t ∈ A ∧ t ∉ B` (`gti.rs:212-245`). -/
theorem memHalfOpen_subtract (t : Int) {A B : List GtiInterval}
    (hB : Normalized B) :
    memHalfOpen t (subtract A B) ↔ memHalfOpen t A ∧ ¬ memHalfOpen t B := by
  induction A with
  | nil => simp [subtract, memHalfOpen]
  | cons a as ih =>
    unfold subtract
    have happ : memHalfOpen t (subtractOne a.start a.stop B ++ subtract as B) ↔
        memHalfOpen t (subtractOne a.start a.stop B) ∨ memHalfOpen t (subtract as B) := by
      unfold memHalfOpen
      constructor
      · intro ⟨iv, hm, ht⟩
        rcases List.mem_append.mp hm with h1 | h2
        · exact Or.inl ⟨iv, h1, ht⟩
        · exact Or.inr ⟨iv, h2, ht⟩
      · intro h
        rcases h with ⟨iv, h1, ht⟩ | ⟨iv, h2, ht⟩
        · exact ⟨iv, List.mem_append_left _ h1, ht⟩
        · exact ⟨iv, List.mem_append_right _ h2, ht⟩
    rw [happ, memHalfOpen_subtractOne t a.start a.stop hB, ih]
    simp only [memHalfOpen, List.mem_cons, GtiInterval.memHalfOpen]
    constructor
    · intro h
      rcases h with ⟨ht, hnot⟩ | ⟨⟨u, hu, hut⟩, hnot⟩
      · exact ⟨⟨a, Or.inl rfl, ht⟩, hnot⟩
      · exact ⟨⟨u, Or.inr hu, hut⟩, hnot⟩
    · intro ⟨⟨u, hu, hut⟩, hnot⟩
      rcases hu with rfl | hu'
      · exact Or.inl ⟨hut, hnot⟩
      · exact Or.inr ⟨⟨u, hu', hut⟩, hnot⟩

/-- Single-interval exact partition identity:
    `exposure(intersectOne s e B) + exposure(subtractOne s e B) = max (e - s) 0`. -/
theorem exposure_partition_one (s e : Int) {B : List GtiInterval} (hB : Normalized B) :
    exposure (intersectOne s e B) + exposure (subtractOne s e B) =
      if s < e then e - s else 0 := by
  induction B generalizing s with
  | nil =>
    unfold intersectOne subtractOne
    split <;> simp [exposure, GtiInterval.duration] <;> omega
  | cons b bs ih =>
    rcases hB with ⟨hb, hbs⟩
    unfold intersectOne subtractOne
    by_cases h1 : b.stop ≤ s
    · simp only [h1, ↓reduceIte]
      exact ih s hbs.toNormalized
    · simp only [h1, ↓reduceIte]
      by_cases h2 : e ≤ b.start
      · simp only [h2, ↓reduceIte]
        split <;> simp [exposure, GtiInterval.duration] <;> omega
      · simp only [h2, ↓reduceIte]
        have hrec := ih (max s b.stop) hbs.toNormalized
        rw [exposure_append]
        by_cases h3 : max s b.stop < e
        · simp only [h3, ↓reduceIte] at hrec ⊢
          rw [exposure_append]
          have hse : s < e := by omega
          have hlo : max s b.start < min e b.stop := by omega
          simp only [hlo, hse, ↓reduceIte]
          by_cases h4 : s < b.start
          · simp only [h4, ↓reduceIte, exposure, GtiInterval.duration]
            simp only [hlo, ↓reduceIte]
            omega
          · simp only [h4, ↓reduceIte, exposure, GtiInterval.duration]
            simp only [hlo, ↓reduceIte]
            omega
        · simp only [h3, ↓reduceIte] at hrec ⊢
          have hrec0 : exposure (intersectOne (max s b.stop) e bs) = 0 := by
            have hge := exposure_nonneg (intersectOne (max s b.stop) e bs)
            have hge2 := exposure_nonneg (subtractOne (max s b.stop) e bs)
            omega
          rw [hrec0]
          by_cases hse : s < e
          · have hlo : max s b.start < min e b.stop := by omega
            simp only [hlo, hse, ↓reduceIte]
            by_cases h4 : s < b.start
            · simp only [h4, ↓reduceIte, exposure, GtiInterval.duration]
              simp only [hlo, ↓reduceIte]
              omega
            · simp only [h4, ↓reduceIte, exposure, GtiInterval.duration]
              simp only [hlo, ↓reduceIte]
              omega
          · have hlo : ¬ (max s b.start < min e b.stop) := by omega
            have h4 : ¬ (s < b.start) := by omega
            simp only [hlo, hse, h4, ↓reduceIte, exposure]
            omega

/-- Equivalence between `intersect` and `flatMap intersectOne` on normalized lists. -/
theorem intersectOne_skip_left {s e : Int} {b : GtiInterval} {bs : List GtiInterval}
    (hle : b.stop ≤ s) :
    intersectOne s e (b :: bs) = intersectOne s e bs := by
  rw [intersectOne]
  simp [hle]

theorem flatMap_intersectOne_nil (L : List GtiInterval) :
    L.flatMap (fun a => intersectOne a.start a.stop []) = [] := by
  induction L with
  | nil => rfl
  | cons _ _ ih => rw [List.flatMap_cons, ih]; rfl

theorem flatMap_congr_eq {α β : Type} {f g : α → List β} {L : List α}
    (h : ∀ x ∈ L, f x = g x) : L.flatMap f = L.flatMap g := by
  induction L with
  | nil => rfl
  | cons x xs ih =>
    simp only [List.flatMap_cons]
    rw [h x (List.mem_cons_self ..), ih (fun y hy => h y (List.mem_cons_of_mem _ hy))]

theorem intersect_eq_flatMap_intersectOne (A B : List GtiInterval)
    (hA : Normalized A) (hB : Normalized B) :
    intersect A B = A.flatMap (fun a => intersectOne a.start a.stop B) := by
  match A, B with
  | [], _ => simp [intersect]
  | a :: as, [] =>
    rw [flatMap_intersectOne_nil]
    simp [intersect]
  | a :: as, b :: bs =>
    rcases hA with ⟨ha, has⟩
    rcases hB with ⟨hb, hbs⟩
    unfold intersect
    dsimp only
    by_cases hcmp : a.stop < b.stop
    · simp only [hcmp, ↓reduceIte]
      have hBcons : Normalized (b :: bs) := ⟨hb, hbs⟩
      rw [intersect_eq_flatMap_intersectOne as (b :: bs) has.toNormalized hBcons]
      rw [List.flatMap_cons]
      have h1 : intersectOne a.start a.stop (b :: bs) =
          if max a.start b.start < min a.stop b.stop then
            [⟨max a.start b.start, min a.stop b.stop⟩]
          else [] := by
        rw [intersectOne]
        by_cases hsk : b.stop ≤ a.start
        · omega
        · simp only [hsk, ↓reduceIte]
          by_cases hdisj : a.stop ≤ b.start
          · simp only [hdisj, ↓reduceIte]
            have hnot : ¬ (max a.start b.start < min a.stop b.stop) := by omega
            simp [hnot]
          · simp only [hdisj, ↓reduceIte]
            have htail : intersectOne (max a.start b.stop) a.stop bs = [] := by
              cases bs with
              | nil => rfl
              | cons b2 bs2 =>
                have hb2 : b.stop < b2.start := hbs.1
                rw [intersectOne]
                have hsk2 : ¬ (b2.stop ≤ max a.start b.stop) := by
                  have := hbs.2.1; omega
                have hdisj2 : a.stop ≤ b2.start := by omega
                simp [hsk2, hdisj2]
            rw [htail, List.append_nil]
      rw [h1]
      split <;> simp
    · simp only [hcmp, ↓reduceIte]
      have hAcons : Normalized (a :: as) := ⟨ha, has⟩
      rw [intersect_eq_flatMap_intersectOne (a :: as) bs hAcons hbs.toNormalized]
      rw [List.flatMap_cons, List.flatMap_cons]
      have has_eq : as.flatMap (fun x => intersectOne x.start x.stop (b :: bs)) =
          as.flatMap (fun x => intersectOne x.start x.stop bs) := by
        apply flatMap_congr_eq
        intro x hx
        have hxgt := (has.forall_start_gt x hx).1
        exact intersectOne_skip_left (by omega)
      rw [has_eq]
      have h1 : intersectOne a.start a.stop (b :: bs) =
          (if max a.start b.start < min a.stop b.stop then
            [⟨max a.start b.start, min a.stop b.stop⟩]
          else []) ++ intersectOne a.start a.stop bs := by
        rw [intersectOne]
        by_cases hsk : b.stop ≤ a.start
        · simp only [hsk, ↓reduceIte]
          have hnot : ¬ (max a.start b.start < min a.stop b.stop) := by omega
          simp [hnot]
        · simp only [hsk, ↓reduceIte]
          by_cases hdisj : a.stop ≤ b.start
          · omega
          · simp only [hdisj, ↓reduceIte]
            have htail : intersectOne (max a.start b.stop) a.stop bs =
                intersectOne a.start a.stop bs := by
              cases bs with
              | nil => rfl
              | cons b2 bs2 =>
                have hb2 : b.stop < b2.start := hbs.1
                have hb2v : b2.start < b2.stop := hbs.2.1
                have h1 : ¬ (b2.stop ≤ max a.start b.stop) := by omega
                have h2 : ¬ (b2.stop ≤ a.start) := by omega
                rw [intersectOne, intersectOne]
                simp only [h1, h2, ↓reduceIte]
                by_cases h3 : a.stop ≤ b2.start
                · simp [h3]
                · simp only [h3, ↓reduceIte]
                  have hmax1 : max (max a.start b.stop) b2.start = b2.start := by omega
                  have hmax2 : max a.start b2.start = b2.start := by omega
                  have hmax3 : max (max a.start b.stop) b2.stop = b2.stop := by omega
                  have hmax4 : max a.start b2.stop = b2.stop := by omega
                  simp [hmax1, hmax2, hmax3, hmax4]
            rw [htail]
      rw [h1]
      split <;> simp
termination_by A.length + B.length

/-- **Exact Partition Identity**: for all normalized GTI sets `A, B`,
    `exposure (A ∩ B) + exposure (A \ B) = exposure A`. -/
theorem exposure_intersect_add_subtract {A B : List GtiInterval}
    (hA : Normalized A) (hB : Normalized B) :
    exposure (intersect A B) + exposure (subtract A B) = exposure A := by
  rw [intersect_eq_flatMap_intersectOne A B hA hB]
  induction A with
  | nil => simp [subtract, exposure]
  | cons a as ih =>
    rcases hA with ⟨ha, has⟩
    simp only [List.flatMap_cons, subtract, exposure_append, exposure, GtiInterval.duration]
    have hone := exposure_partition_one a.start a.stop hB
    have hrec := ih has.toNormalized
    simp only [ha, ↓reduceIte] at hone ⊢
    omega

/-- **Exposure Monotonicity (Intersection)**: `exposure (A ∩ B) ≤ exposure A`. -/
theorem exposure_intersect_le_left {A B : List GtiInterval}
    (hA : Normalized A) (hB : Normalized B) :
    exposure (intersect A B) ≤ exposure A := by
  have hpart := exposure_intersect_add_subtract hA hB
  have hnn := exposure_nonneg (subtract A B)
  omega

/-- **Exposure Monotonicity (Subtraction)**: `exposure (A \ B) ≤ exposure A`. -/
theorem exposure_subtract_le {A B : List GtiInterval}
    (hA : Normalized A) (hB : Normalized B) :
    exposure (subtract A B) ≤ exposure A := by
  have hpart := exposure_intersect_add_subtract hA hB
  have hnn := exposure_nonneg (intersect A B)
  omega

/-! ## 5. `drop_short` and `merge_gaps` Invariants (`gti.rs:154-180`) -/

/-- Exact mirror of `GtiSet::drop_short` (`gti.rs:172-180`). -/
def dropShort (L : List GtiInterval) (minExp : Int) : List GtiInterval :=
  L.filter (fun iv => decide (minExp ≤ iv.duration))

theorem dropShort_normalizedAfter {b : Int} {L : List GtiInterval} (minExp : Int)
    (h : NormalizedAfter b L) : NormalizedAfter b (dropShort L minExp) := by
  induction L generalizing b with
  | nil => trivial
  | cons a as ih =>
    rcases h with ⟨hba, ha, has⟩
    unfold dropShort
    simp only [List.filter_cons]
    split
    · exact ⟨hba, ha, ih has⟩
    · exact (ih has).mono (by omega)

/-- `dropShort` preserves `Normalized` (`gti.rs:172-180`). -/
theorem dropShort_normalized {L : List GtiInterval} (minExp : Int) (h : Normalized L) :
    Normalized (dropShort L minExp) := by
  cases L with
  | nil => trivial
  | cons a as =>
    have hA : NormalizedAfter (a.start - 1) (a :: as) := ⟨by omega, h.1, h.2⟩
    exact (dropShort_normalizedAfter minExp hA).toNormalized

/-- Every interval retained by `dropShort` has `duration ≥ minExp`. -/
theorem dropShort_forall_ge (L : List GtiInterval) (minExp : Int) :
    ∀ iv ∈ dropShort L minExp, minExp ≤ iv.duration := by
  intro iv hm
  unfold dropShort at hm
  rw [List.mem_filter] at hm
  exact of_decide_eq_true hm.2

/-- `dropShort` never increases total exposure. -/
theorem exposure_dropShort_le (L : List GtiInterval) (minExp : Int) :
    exposure (dropShort L minExp) ≤ exposure L := by
  induction L with
  | nil => simp [dropShort, exposure]
  | cons a as ih =>
    unfold dropShort at ih ⊢
    simp only [List.filter_cons, exposure]
    have ha := duration_nonneg a
    split
    · simp only [exposure]; omega
    · omega

/-- Exact mirror of `GtiSet::merge_gaps` (`gti.rs:154-169`). -/
def mergeGapsAux (maxGap : Int) (cur : GtiInterval) : List GtiInterval → List GtiInterval
  | [] => [cur]
  | iv :: rest =>
    if iv.start - cur.stop < maxGap then
      mergeGapsAux maxGap ⟨cur.start, max cur.stop iv.stop⟩ rest
    else
      cur :: mergeGapsAux maxGap iv rest

def mergeGaps (L : List GtiInterval) (maxGap : Int) : List GtiInterval :=
  if maxGap ≤ 0 then L
  else
    match L with
    | [] => []
    | iv :: rest => mergeGapsAux maxGap iv rest

/-- Minimum gap invariant: every adjacent pair `(x, y)` has `y.start - x.stop ≥ g`. -/
def MinGapAfter (g b : Int) : List GtiInterval → Prop
  | [] => True
  | iv :: rest => b + g ≤ iv.start ∧ iv.start < iv.stop ∧ MinGapAfter g iv.stop rest

def MinGap (g : Int) : List GtiInterval → Prop
  | [] => True
  | iv :: rest => iv.start < iv.stop ∧ MinGapAfter g iv.stop rest

theorem MinGapAfter.toNormalizedAfter {g b : Int} {L : List GtiInterval}
    (hg : 0 < g) (h : MinGapAfter g b L) : NormalizedAfter b L := by
  induction L generalizing b with
  | nil => trivial
  | cons a as ih =>
    rcases h with ⟨h1, h2, h3⟩
    exact ⟨by omega, h2, ih h3⟩

theorem MinGapAfter.toMinGap {g b : Int} {L : List GtiInterval}
    (h : MinGapAfter g b L) : MinGap g L := by
  cases L with
  | nil => trivial
  | cons a as => exact ⟨h.2.1, h.2.2⟩

theorem mergeGapsAux_minGapAfter {maxGap b : Int} {cur : GtiInterval} {rest : List GtiInterval}
    (_hmg : 0 < maxGap) (hb : b + maxGap ≤ cur.start) (hcur : cur.start < cur.stop)
    (hrest : NormalizedAfter cur.stop rest) :
    MinGapAfter maxGap b (mergeGapsAux maxGap cur rest) := by
  induction rest generalizing b cur with
  | nil => exact ⟨hb, hcur, trivial⟩
  | cons iv ys ih =>
    rcases hrest with ⟨hcs, hiv, hys⟩
    unfold mergeGapsAux
    by_cases hgap : iv.start - cur.stop < maxGap
    · simp only [hgap, ↓reduceIte]
      have hmax : max cur.stop iv.stop = iv.stop := by omega
      exact ih (cur := ⟨cur.start, max cur.stop iv.stop⟩) hb (by dsimp; omega) (by dsimp; rw [hmax]; exact hys)
    · simp only [hgap, ↓reduceIte]
      exact ⟨hb, hcur, ih (by omega) hiv hys⟩

/-- `mergeGaps` preserves `Normalized` for all `maxGap` (`gti.rs:154-169`). -/
theorem mergeGaps_normalized {L : List GtiInterval} (maxGap : Int) (h : Normalized L) :
    Normalized (mergeGaps L maxGap) := by
  unfold mergeGaps
  split
  · exact h
  · rename_i hpos
    cases L with
    | nil => trivial
    | cons cur rest =>
      rcases h with ⟨hcur, hrest⟩
      have hmg : 0 < maxGap := by omega
      have hmin := mergeGapsAux_minGapAfter (b := cur.start - maxGap) hmg (by omega) hcur hrest
      exact (hmin.toNormalizedAfter hmg).toNormalized

/-- When `maxGap > 0`, every remaining gap in `mergeGaps L maxGap` is at least `maxGap` (`gti.rs:161`). -/
theorem mergeGaps_minGap {L : List GtiInterval} {maxGap : Int}
    (hmg : 0 < maxGap) (h : Normalized L) : MinGap maxGap (mergeGaps L maxGap) := by
  unfold mergeGaps
  have hnot : ¬ (maxGap ≤ 0) := by omega
  simp only [hnot, ↓reduceIte]
  cases L with
  | nil => trivial
  | cons cur rest =>
    rcases h with ⟨hcur, hrest⟩
    exact (mergeGapsAux_minGapAfter (b := cur.start - maxGap) hmg (by omega) hcur hrest).toMinGap

/-- Merging gaps never decreases exposure (`exposure L ≤ exposure (mergeGaps L maxGap)`). -/
theorem exposure_le_mergeGapsAux (maxGap : Int) (cur : GtiInterval) (rest : List GtiInterval)
    (hcur : cur.start < cur.stop) (hrest : NormalizedAfter cur.stop rest) :
    cur.duration + exposure rest ≤ exposure (mergeGapsAux maxGap cur rest) := by
  induction rest generalizing cur with
  | nil => simp [mergeGapsAux, exposure]
  | cons iv ys ih =>
    rcases hrest with ⟨hcs, hiv, hys⟩
    unfold mergeGapsAux
    by_cases hgap : iv.start - cur.stop < maxGap
    · simp only [hgap, ↓reduceIte]
      have hmax : max cur.stop iv.stop = iv.stop := by omega
      have hlt' : cur.start < max cur.stop iv.stop := by omega
      have hrec := ih ⟨cur.start, max cur.stop iv.stop⟩ hlt' (by dsimp; rw [hmax]; exact hys)
      simp only [exposure, GtiInterval.duration, hlt', hcur, hiv, ↓reduceIte] at hrec ⊢
      omega
    · simp only [hgap, ↓reduceIte, exposure]
      have hrec := ih iv hiv hys
      omega

theorem exposure_le_mergeGaps {L : List GtiInterval} (maxGap : Int) (h : Normalized L) :
    exposure L ≤ exposure (mergeGaps L maxGap) := by
  unfold mergeGaps
  split
  · omega
  · cases L with
    | nil => simp [exposure]
    | cons cur rest =>
      rcases h with ⟨hcur, hrest⟩
      simp only [exposure]
      exact exposure_le_mergeGapsAux maxGap cur rest hcur hrest

end AstrosatFormal.GtiAlgebra
