/-
Copyright (c) 2026 Siddhant Sood. All rights reserved.
Authors: Siddhant Sood
-/
import KuramotoBlockLimit.Quotient
import Mathlib.Algebra.Order.Chebyshev

set_option linter.style.header false

/-!
# Clique blow-ups of a Kuramoto Hessian

This file formalizes the algebraic core of the clique-blow-up theorem.  The
base Hessian is represented by its action on a finite vertex type `V`.  A
clone block has size `m`; the blow-up operator below is exactly the Hessian of
the lexicographic product by `K_m`, with the unnormalized Kuramoto convention.

The main results prove its action on the two invariant summands and prove that
its kernel consists only of the global rotation vector when the same is true
for the base Hessian and the base diagonal is nonnegative.
-/

namespace KuramotoBlockLimit

noncomputable section

open scoped BigOperators

/-- Action of a finite matrix on a vector. -/
def finiteMatrixAction
    {V : Type*} [Fintype V]
    (H : V → V → ℝ) (x : V → ℝ) (v : V) : ℝ :=
  ∑ w : V, H v w * x w

/-- Sum of a vector over all clones of one base vertex. -/
def cloneSum
    {V : Type*} {m : ℕ}
    (x : V × Fin m → ℝ) (v : V) : ℝ :=
  ∑ a : Fin m, x (v, a)

/-- Lift a base vector to a vector constant on every clone block. -/
def liftToClones
    {V : Type*} {m : ℕ}
    (x : V → ℝ) : V × Fin m → ℝ :=
  fun p ↦ x p.1

/-- Kuramoto torque for a finite weighted coupling matrix. -/
def finiteTorque
    {V : Type*} [Fintype V]
    (A : V → V → ℝ) (phase : V → ℝ) (v : V) : ℝ :=
  ∑ w : V, A v w * Real.sin (phase w - phase v)

/-- Coupling matrix of the lexicographic clique blow-up. -/
def cliqueBlowupCoupling
    {V : Type*} [DecidableEq V] {m : ℕ}
    (A : V → V → ℝ) (p q : V × Fin m) : ℝ :=
  if p.1 = q.1 then
    if p.2 = q.2 then 0 else 1
  else A p.1 q.1

theorem cliqueBlowupCoupling_loopless
    {V : Type*} [DecidableEq V] {m : ℕ}
    (A : V → V → ℝ) (p : V × Fin m) :
    cliqueBlowupCoupling A p p = 0 := by
  simp [cliqueBlowupCoupling]

theorem cliqueBlowupCoupling_symmetric
    {V : Type*} [DecidableEq V] {m : ℕ}
    (A : V → V → ℝ) (hA : ∀ v w, A v w = A w v)
    (p q : V × Fin m) :
    cliqueBlowupCoupling A p q = cliqueBlowupCoupling A q p := by
  obtain ⟨v, a⟩ := p
  obtain ⟨w, b⟩ := q
  by_cases hvw : v = w
  · subst w
    by_cases hab : a = b
    · subst b
      rfl
    · simp [cliqueBlowupCoupling, hab, Ne.symm hab]
  · simp [cliqueBlowupCoupling, hvw, Ne.symm hvw, hA]

/-- The standard unnormalized cosine-weighted Kuramoto potential Hessian. -/
def kuramotoHessian
    {V : Type*} [Fintype V] [DecidableEq V]
    (A : V → V → ℝ) (phase : V → ℝ) : Matrix V V ℝ :=
  fun v w ↦
    if v = w then
      ∑ u : V, A v u * Real.cos (phase v - phase u)
    else -A v w * Real.cos (phase v - phase w)

/-- Kuramoto torque on the clique blow-up. -/
def cliqueBlowupTorque
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (A : V → V → ℝ) (phase : V × Fin m → ℝ) (p : V × Fin m) : ℝ :=
  ∑ q : V × Fin m,
    cliqueBlowupCoupling A p q * Real.sin (phase q - phase p)

/-- Explicit block matrix of the unnormalized clique-blow-up Hessian. -/
def cliqueBlowupHessian
    {V : Type*} [DecidableEq V] (m : ℕ) (H : V → V → ℝ) :
    Matrix (V × Fin m) (V × Fin m) ℝ :=
  fun p q ↦
    if p.1 = q.1 then
      if p.2 = q.2 then (m : ℝ) * (1 + H p.1 p.1) - 1 else -1
    else H p.1 q.1

/--
The unnormalized Hessian action after replacing every base vertex by `K_m`.

The diagonal clone block is `m(1+H_vv)I-J`; every off-diagonal clone block is
the constant matrix with entry `H_vw`.
-/
def cliqueBlowupAction
    {V : Type*} [Fintype V] [DecidableEq V] (m : ℕ)
    (H : V → V → ℝ) (x : V × Fin m → ℝ) (p : V × Fin m) : ℝ :=
  (m : ℝ) * (1 + H p.1 p.1) * x p - cloneSum x p.1 +
    ∑ w ∈ (Finset.univ.erase p.1), H p.1 w * cloneSum x w

/-- The displayed action formula is exactly matrix-vector multiplication. -/
theorem cliqueBlowupHessian_mulVec
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (x : V × Fin m → ℝ) (v : V) (a : Fin m) :
    (cliqueBlowupHessian m H).mulVec x (v, a) =
      cliqueBlowupAction m H x (v, a) := by
  change (∑ q : V × Fin m, cliqueBlowupHessian m H (v, a) q * x q) = _
  rw [Fintype.sum_prod_type]
  rw [← Finset.add_sum_erase Finset.univ
    (fun w ↦ ∑ b : Fin m, cliqueBlowupHessian m H (v, a) (w, b) * x (w, b))
    (Finset.mem_univ v)]
  have hsame :
      (∑ b : Fin m, cliqueBlowupHessian m H (v, a) (v, b) * x (v, b)) =
        (m : ℝ) * (1 + H v v) * x (v, a) - cloneSum x v := by
    rw [cloneSum]
    rw [← Finset.add_sum_erase Finset.univ (fun b ↦ x (v, b)) (Finset.mem_univ a)]
    rw [← Finset.add_sum_erase Finset.univ
      (fun b ↦ cliqueBlowupHessian m H (v, a) (v, b) * x (v, b))
      (Finset.mem_univ a)]
    simp only [cliqueBlowupHessian, reduceIte]
    have hoff :
        (∑ b ∈ Finset.univ.erase a,
          (if a = b then (m : ℝ) * (1 + H v v) - 1 else -1) * x (v, b)) =
          -∑ b ∈ Finset.univ.erase a, x (v, b) := by
      rw [← Finset.sum_neg_distrib]
      apply Finset.sum_congr rfl
      intro b hb
      have hab : a ≠ b := Ne.symm (Finset.mem_erase.mp hb).1
      simp [hab]
    rw [hoff]
    ring
  rw [hsame]
  have hoff : ∀ w ∈ Finset.univ.erase v,
      (∑ b : Fin m, cliqueBlowupHessian m H (v, a) (w, b) * x (w, b)) =
        H v w * cloneSum x w := by
    intro w hw
    have hvw : v ≠ w := Ne.symm (Finset.mem_erase.mp hw).1
    simp [cliqueBlowupHessian, hvw, cloneSum, Finset.mul_sum]
  apply congrArg (fun t : ℝ ↦
    (m : ℝ) * (1 + H v v) * x (v, a) - cloneSum x v + t)
  apply Finset.sum_congr rfl
  intro w hw
  exact hoff w hw

/-- Symmetry of the base Hessian implies symmetry of the explicit blow-up Hessian. -/
theorem cliqueBlowupHessian_isSymm
    {V : Type*} [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (hH : ∀ v w, H v w = H w v) :
    (cliqueBlowupHessian m H).IsSymm := by
  apply Matrix.IsSymm.ext
  rintro ⟨v, a⟩ ⟨w, b⟩
  by_cases hvw : v = w
  · subst w
    by_cases hab : a = b
    · subst b
      rfl
    · simp [cliqueBlowupHessian, hab, Ne.symm hab]
  · simp [cliqueBlowupHessian, hvw, Ne.symm hvw, hH]

/--
The explicit block matrix is the actual Kuramoto Hessian of the lexicographic
clique blow-up at the replicated phase assignment.
-/
theorem kuramotoHessian_cliqueBlowup
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (A : V → V → ℝ) (phase : V → ℝ) (hloop : ∀ v, A v v = 0) :
    kuramotoHessian (cliqueBlowupCoupling (m := m) A)
        (liftToClones (m := m) phase) =
      cliqueBlowupHessian m (kuramotoHessian A phase) := by
  funext p q
  obtain ⟨v, a⟩ := p
  obtain ⟨w, b⟩ := q
  by_cases hvw : v = w
  · subst w
    by_cases hab : a = b
    · subst b
      simp only [kuramotoHessian, cliqueBlowupHessian, cliqueBlowupCoupling,
        liftToClones, reduceIte]
      rw [Fintype.sum_prod_type]
      rw [← Finset.add_sum_erase Finset.univ
        (fun u ↦ ∑ c : Fin m,
          (if v = u then if a = c then 0 else 1 else A v u) *
            Real.cos (phase v - phase u)) (Finset.mem_univ v)]
      have hsame :
          (∑ c : Fin m, (if a = c then 0 else 1) *
            Real.cos (phase v - phase v)) = (m : ℝ) - 1 := by
        simp only [sub_self, Real.cos_zero, mul_one]
        calc
          (∑ c : Fin m, if a = c then (0 : ℝ) else 1) =
              ∑ c : Fin m, if a ≠ c then (1 : ℝ) else 0 := by
                apply Finset.sum_congr rfl
                intro c _
                by_cases hac : a = c <;> simp [hac]
          _ = ((Finset.univ.erase a).card : ℝ) := by
                rw [Finset.sum_boole]
                simp [Finset.filter_ne]
          _ = (m : ℝ) - 1 := by
                rw [Finset.card_erase_of_mem (Finset.mem_univ a), Finset.card_fin]
                rw [Nat.cast_sub]
                · norm_num
                · have ha := a.isLt
                  omega
      simp only [ite_true]
      rw [hsame]
      have hoff : ∀ u ∈ Finset.univ.erase v,
          (∑ c : Fin m, (if v = u then if a = c then 0 else 1 else A v u) *
            Real.cos (phase v - phase u)) =
            (m : ℝ) * (A v u * Real.cos (phase v - phase u)) := by
        intro u hu
        have hvu : v ≠ u := Ne.symm (Finset.mem_erase.mp hu).1
        simp [hvu]
      rw [Finset.sum_congr rfl hoff]
      rw [← Finset.mul_sum]
      have hdiagbase :
          (∑ u : V, A v u * Real.cos (phase v - phase u)) =
            ∑ u ∈ Finset.univ.erase v, A v u * Real.cos (phase v - phase u) := by
        symm
        apply Finset.sum_erase
        simp [hloop v]
      rw [hdiagbase]
      ring
    · simp [kuramotoHessian, cliqueBlowupHessian, cliqueBlowupCoupling,
        liftToClones, hab]
  · simp [kuramotoHessian, cliqueBlowupHessian, cliqueBlowupCoupling,
      liftToClones, hvw]

/-- Quadratic form of the base Hessian. -/
def finiteQuadratic
    {V : Type*} [Fintype V]
    (H : V → V → ℝ) (x : V → ℝ) : ℝ :=
  ∑ v : V, x v * finiteMatrixAction H x v

/-- Quadratic form of the clique-blow-up Hessian. -/
def cliqueBlowupQuadratic
    {V : Type*} [Fintype V] [DecidableEq V] (m : ℕ)
    (H : V → V → ℝ) (x : V × Fin m → ℝ) : ℝ :=
  ∑ p : V × Fin m, x p * cliqueBlowupAction m H x p

/-- The unnormalized within-block variance numerator. -/
def cloneVariance
    {V : Type*} (m : ℕ) (x : V × Fin m → ℝ) (v : V) : ℝ :=
  (m : ℝ) * ∑ a : Fin m, (x (v, a)) ^ 2 - (cloneSum x v) ^ 2

theorem cloneSum_liftToClones
    {V : Type*} {m : ℕ} (x : V → ℝ) (v : V) :
    cloneSum (m := m) (liftToClones (m := m) x) v = (m : ℝ) * x v := by
  simp [cloneSum, liftToClones]

/-- Replication preserves the existence of two distinct phase values. -/
theorem liftToClones_nonsynchronous
    {V : Type*} {m : ℕ} (hm : 0 < m) (phase : V → ℝ)
    (hphase : ∃ v w, phase v ≠ phase w) :
    ∃ p q : V × Fin m,
      liftToClones (m := m) phase p ≠ liftToClones (m := m) phase q := by
  obtain ⟨v, w, hvw⟩ := hphase
  let a : Fin m := ⟨0, hm⟩
  exact ⟨(v, a), (w, a), hvw⟩

/-- Replicating an equilibrium multiplies every torque by the clone count. -/
theorem cliqueBlowupTorque_lift
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (A : V → V → ℝ) (phase : V → ℝ) (v : V) (a : Fin m) :
    cliqueBlowupTorque A (liftToClones (m := m) phase) (v, a) =
      (m : ℝ) * finiteTorque A phase v := by
  rw [cliqueBlowupTorque, Fintype.sum_prod_type, finiteTorque]
  have hinner : ∀ w : V,
      (∑ b : Fin m,
        cliqueBlowupCoupling A (v, a) (w, b) *
          Real.sin (liftToClones phase (w, b) - liftToClones phase (v, a))) =
        (m : ℝ) * (A v w * Real.sin (phase w - phase v)) := by
    intro w
    by_cases hvw : v = w
    · subst w
      simp [cliqueBlowupCoupling, liftToClones]
    · simp [cliqueBlowupCoupling, liftToClones, hvw]
  simp_rw [hinner]
  rw [Finset.mul_sum]

/-- On block-constant vectors, the blow-up Hessian is `m` times the base Hessian. -/
theorem cliqueBlowupAction_lift
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (x : V → ℝ) (v : V) (a : Fin m) :
    cliqueBlowupAction m H (liftToClones (m := m) x) (v, a) =
      (m : ℝ) * finiteMatrixAction H x v := by
  rw [cliqueBlowupAction]
  simp only [liftToClones, cloneSum_liftToClones]
  rw [finiteMatrixAction]
  rw [← Finset.add_sum_erase Finset.univ (fun w ↦ H v w * x w) (Finset.mem_univ v)]
  have hoff :
      (∑ w ∈ Finset.univ.erase v, H v w * ((m : ℝ) * x w)) =
        (m : ℝ) * ∑ w ∈ Finset.univ.erase v, H v w * x w := by
    rw [Finset.mul_sum]
    apply Finset.sum_congr rfl
    intro w _
    ring
  rw [hoff]
  ring

/-- The clone sum of a zero-sum block vector is zero by definition. -/
theorem cliqueBlowupAction_zeroSum
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (z : V × Fin m → ℝ)
    (hz : ∀ v, cloneSum z v = 0) (v : V) (a : Fin m) :
    cliqueBlowupAction m H z (v, a) =
      (m : ℝ) * (1 + H v v) * z (v, a) := by
  simp [cliqueBlowupAction, hz]

/-- Summing the blow-up action over one clone block recovers the base action. -/
theorem cloneSum_cliqueBlowupAction
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (x : V × Fin m → ℝ) (v : V) :
    cloneSum (fun p ↦ cliqueBlowupAction m H x p) v =
      (m : ℝ) * finiteMatrixAction H (cloneSum x) v := by
  rw [cloneSum, finiteMatrixAction]
  simp_rw [cliqueBlowupAction]
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib]
  simp only [Finset.sum_const, Finset.card_fin, nsmul_eq_mul]
  have hdiag :
      (∑ b : Fin m, (m : ℝ) * (1 + H v v) * x (v, b)) =
        (m : ℝ) * (1 + H v v) * cloneSum x v := by
    simp only [cloneSum]
    rw [Finset.mul_sum]
  rw [hdiag]
  rw [← Finset.add_sum_erase Finset.univ
    (fun w ↦ H v w * cloneSum x w) (Finset.mem_univ v)]
  ring

/-- Exact quadratic-form contribution of one clone block. -/
theorem cliqueBlowupQuadratic_block
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (x : V × Fin m → ℝ) (v : V) :
    (∑ a : Fin m, x (v, a) * cliqueBlowupAction m H x (v, a)) =
      cloneSum x v * finiteMatrixAction H (cloneSum x) v +
        (1 + H v v) * cloneVariance m x v := by
  have hbase :
      finiteMatrixAction H (cloneSum x) v =
        H v v * cloneSum x v +
          ∑ w ∈ Finset.univ.erase v, H v w * cloneSum x w := by
    rw [finiteMatrixAction]
    rw [← Finset.add_sum_erase Finset.univ
      (fun w ↦ H v w * cloneSum x w) (Finset.mem_univ v)]
  rw [hbase]
  simp_rw [cliqueBlowupAction]
  have hpoint : ∀ a : Fin m,
      x (v, a) *
          ((m : ℝ) * (1 + H v v) * x (v, a) - cloneSum x v +
            ∑ w ∈ Finset.univ.erase v, H v w * cloneSum x w) =
        (m : ℝ) * (1 + H v v) * (x (v, a)) ^ 2 -
          cloneSum x v * x (v, a) +
          (∑ w ∈ Finset.univ.erase v, H v w * cloneSum x w) * x (v, a) := by
    intro a
    ring
  simp_rw [hpoint]
  rw [Finset.sum_add_distrib, Finset.sum_sub_distrib]
  rw [← Finset.mul_sum, ← Finset.mul_sum, ← Finset.mul_sum]
  simp only [cloneVariance, cloneSum]
  ring

/-- Exact orthogonal-splitting identity for the full blow-up quadratic form. -/
theorem cliqueBlowupQuadratic_decomposition
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ) (x : V × Fin m → ℝ) :
    cliqueBlowupQuadratic m H x =
      finiteQuadratic H (cloneSum x) +
        ∑ v : V, (1 + H v v) * cloneVariance m x v := by
  rw [cliqueBlowupQuadratic, Fintype.sum_prod_type, finiteQuadratic]
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro v _
  exact cliqueBlowupQuadratic_block H x v

/-- Cauchy--Schwarz makes every within-block variance numerator nonnegative. -/
theorem cloneVariance_nonnegative
    {V : Type*} {m : ℕ} (x : V × Fin m → ℝ) (v : V) :
    0 ≤ cloneVariance m x v := by
  apply sub_nonneg.mpr
  simpa [cloneVariance, cloneSum] using
    (sq_sum_le_card_mul_sum_sq
      (s := (Finset.univ : Finset (Fin m))) (f := fun a ↦ x (v, a)))

/-- A nonnegative quadratic form has nonnegative diagonal entries. -/
theorem diagonal_nonnegative_of_quadratic_nonnegative
    {V : Type*} [Fintype V]
    (H : V → V → ℝ) (hbase : ∀ y : V → ℝ, 0 ≤ finiteQuadratic H y) :
    ∀ v, 0 ≤ H v v := by
  classical
  intro v
  let e : V → ℝ := fun w ↦ if w = v then 1 else 0
  have he := hbase e
  simpa [finiteQuadratic, finiteMatrixAction, e] using he

/-- Positive semidefiniteness of the base Hessian passes to every clique blow-up. -/
theorem cliqueBlowupQuadratic_nonnegative
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ)
    (hbase : ∀ y : V → ℝ, 0 ≤ finiteQuadratic H y)
    (hdiag : ∀ v, 0 ≤ H v v)
    (x : V × Fin m → ℝ) :
    0 ≤ cliqueBlowupQuadratic m H x := by
  rw [cliqueBlowupQuadratic_decomposition]
  apply add_nonneg (hbase (cloneSum x))
  apply Finset.sum_nonneg
  intro v _
  exact mul_nonneg (by linarith [hdiag v]) (cloneVariance_nonnegative x v)

/-- Every globally constant vector is in the blow-up kernel when base rows sum to zero. -/
theorem cliqueBlowupAction_rotation
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (H : V → V → ℝ)
    (hrow : ∀ v, finiteMatrixAction H (fun _ ↦ 1) v = 0)
    (c : ℝ) (v : V) (a : Fin m) :
    cliqueBlowupAction m H (fun _ ↦ c) (v, a) = 0 := by
  have hbase : finiteMatrixAction H (fun _ ↦ c) v = 0 := by
    calc
      finiteMatrixAction H (fun _ ↦ c) v = c * ∑ w : V, H v w := by
        rw [finiteMatrixAction, Finset.mul_sum]
        apply Finset.sum_congr rfl
        intro w _
        ring
      _ = c * finiteMatrixAction H (fun _ ↦ 1) v := by
        simp [finiteMatrixAction]
      _ = 0 := by rw [hrow v, mul_zero]
  have hlift : (fun _ : V × Fin m ↦ c) = liftToClones (m := m) (fun _ : V ↦ c) := by
    funext p
    rfl
  rw [hlift, cliqueBlowupAction_lift, hbase, mul_zero]

/--
The clique blow-up has no new zero mode.  The hypotheses say that the base
matrix has zero row sums, that its kernel contains only constant vectors, and
that its diagonal is nonnegative (the latter follows from positive
semidefiniteness for an actual Hessian).
-/
theorem cliqueBlowup_kernel_eq_rotation
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (hm : 0 < m) (H : V → V → ℝ)
    (hrow : ∀ v, finiteMatrixAction H (fun _ ↦ 1) v = 0)
    (hbaseKernel : ∀ y : V → ℝ,
      (∀ v, finiteMatrixAction H y v = 0) → ∃ c, ∀ v, y v = c)
    (hdiag : ∀ v, 0 ≤ H v v)
    (x : V × Fin m → ℝ) :
    (∀ p, cliqueBlowupAction m H x p = 0) ↔ ∃ c, ∀ p, x p = c := by
  constructor
  · intro hx
    have hmreal : (m : ℝ) ≠ 0 := by exact_mod_cast (Nat.ne_of_gt hm)
    have hsumKernel : ∀ v, finiteMatrixAction H (cloneSum x) v = 0 := by
      intro v
      have hsumZero :
          cloneSum (fun p ↦ cliqueBlowupAction m H x p) v = 0 := by
        simp [cloneSum, hx]
      rw [cloneSum_cliqueBlowupAction] at hsumZero
      exact (mul_eq_zero.mp hsumZero).resolve_left hmreal
    obtain ⟨s, hs⟩ := hbaseKernel (cloneSum x) hsumKernel
    refine ⟨s / (m : ℝ), ?_⟩
    rintro ⟨v, a⟩
    have hrowsum : H v v + ∑ w ∈ Finset.univ.erase v, H v w = 0 := by
      have hv := hrow v
      rw [finiteMatrixAction] at hv
      have hv' : (∑ w : V, H v w) = 0 := by simpa using hv
      rw [← Finset.add_sum_erase Finset.univ (fun w ↦ H v w) (Finset.mem_univ v)] at hv'
      exact hv'
    have hoff :
        (∑ w ∈ Finset.univ.erase v, H v w * cloneSum x w) = -H v v * s := by
      calc
        (∑ w ∈ Finset.univ.erase v, H v w * cloneSum x w) =
            ∑ w ∈ Finset.univ.erase v, H v w * s := by
              apply Finset.sum_congr rfl
              intro w _
              rw [hs w]
        _ = s * ∑ w ∈ Finset.univ.erase v, H v w := by
              rw [Finset.mul_sum]
              apply Finset.sum_congr rfl
              intro w _
              ring
        _ = s * (-H v v) := by
              congr 1
              linarith [hrowsum]
        _ = -H v v * s := by ring
    have heq := hx (v, a)
    rw [cliqueBlowupAction, hs v, hoff] at heq
    have hfactor : (1 + H v v) * ((m : ℝ) * x (v, a) - s) = 0 := by
      nlinarith
    have hcoefficient : 1 + H v v ≠ 0 := by
      have : 0 < 1 + H v v := by linarith [hdiag v]
      exact ne_of_gt this
    have hmx : (m : ℝ) * x (v, a) = s := by
      exact sub_eq_zero.mp ((mul_eq_zero.mp hfactor).resolve_left hcoefficient)
    apply (eq_div_iff hmreal).2
    simpa [mul_comm] using hmx
  · rintro ⟨c, hc⟩ p
    obtain ⟨v, a⟩ := p
    have hxconstant : x = fun _ ↦ c := by
      funext q
      exact hc q
    rw [hxconstant]
    exact cliqueBlowupAction_rotation H hrow c v a

/--
Complete machine-checked clique-blow-up theorem for the homogeneous
first-order Kuramoto model: a replicated equilibrium remains an equilibrium,
its actual potential Hessian remains positive semidefinite, and its Hessian
kernel contains only global rotation.
-/
theorem cliqueBlowup_preserves_stable_equilibrium
    {V : Type*} [Fintype V] [DecidableEq V] {m : ℕ}
    (hm : 0 < m) (A : V → V → ℝ) (phase : V → ℝ)
    (hloop : ∀ v, A v v = 0)
    (hequilibrium : ∀ v, finiteTorque A phase v = 0)
    (hrow : ∀ v,
      finiteMatrixAction (kuramotoHessian A phase) (fun _ ↦ 1) v = 0)
    (hbaseNonnegative : ∀ y : V → ℝ,
      0 ≤ finiteQuadratic (kuramotoHessian A phase) y)
    (hbaseKernel : ∀ y : V → ℝ,
      (∀ v, finiteMatrixAction (kuramotoHessian A phase) y v = 0) →
        ∃ c, ∀ v, y v = c) :
    (∀ p, cliqueBlowupTorque A (liftToClones (m := m) phase) p = 0) ∧
      (∀ x : V × Fin m → ℝ,
        0 ≤ ∑ p : V × Fin m, x p *
          (kuramotoHessian (cliqueBlowupCoupling (m := m) A)
            (liftToClones (m := m) phase)).mulVec x p) ∧
      (∀ x : V × Fin m → ℝ,
        ((∀ p, (kuramotoHessian (cliqueBlowupCoupling (m := m) A)
            (liftToClones (m := m) phase)).mulVec x p = 0) ↔
          ∃ c, ∀ p, x p = c)) := by
  constructor
  · rintro ⟨v, a⟩
    rw [cliqueBlowupTorque_lift, hequilibrium v, mul_zero]
  constructor
  · intro x
    rw [kuramotoHessian_cliqueBlowup A phase hloop]
    simp_rw [cliqueBlowupHessian_mulVec
      (fun v w ↦ kuramotoHessian A phase v w) x]
    exact cliqueBlowupQuadratic_nonnegative
      (kuramotoHessian A phase) hbaseNonnegative
      (diagonal_nonnegative_of_quadratic_nonnegative
        (kuramotoHessian A phase) hbaseNonnegative) x
  · intro x
    rw [kuramotoHessian_cliqueBlowup A phase hloop]
    simp_rw [cliqueBlowupHessian_mulVec
      (fun v w ↦ kuramotoHessian A phase v w) x]
    exact cliqueBlowup_kernel_eq_rotation hm (kuramotoHessian A phase)
      hrow hbaseKernel
      (diagonal_nonnegative_of_quadratic_nonnegative
        (kuramotoHessian A phase) hbaseNonnegative) x

end

end KuramotoBlockLimit
