/-
Copyright (c) 2026 Siddhant Sood. All rights reserved.
Authors: Siddhant Sood
-/
import KuramotoBlockLimit.Basic

set_option linter.style.header false

/-!
# Exact equilibrium lifting from equitable neighbor counts

The theorem here is local to one vertex. Applying it to every vertex of an
equitable partition gives the exact lifting theorem in the accompanying paper
notes, without any limiting argument.
-/

namespace KuramotoBlockLimit

noncomputable section

open scoped BigOperators

/-- Regroup a finite sum of a class-dependent function by class. -/
theorem sum_by_class
    {V I : Type*} [Fintype I] [DecidableEq I]
    (neighbors : Finset V) (cls : V → I) (f : I → ℝ) :
    (∑ v ∈ neighbors, f (cls v)) =
      ∑ i : I, ((neighbors.filter fun v ↦ cls v = i).card : ℝ) * f i := by
  rw [← Finset.sum_fiberwise neighbors cls (fun v ↦ f (cls v))]
  apply Finset.sum_congr rfl
  intro i _
  calc
    (∑ v ∈ neighbors with cls v = i, f (cls v)) =
        ∑ v ∈ neighbors with cls v = i, f i := by
      apply Finset.sum_congr rfl
      intro v hv
      exact congrArg f (Finset.mem_filter.mp hv).2
    _ = ((neighbors.filter fun v ↦ cls v = i).card : ℝ) * f i := by simp

/-- Kuramoto torque contributed by a finite neighbor set at a vertex in class `i`. -/
def vertexTorque
    {V I : Type*} [Fintype I]
    (neighbors : Finset V) (cls : V → I) (phase : I → ℝ) (i : I) : ℝ :=
  ∑ v ∈ neighbors, Real.sin (phase (cls v) - phase i)

/-- Reduced class torque for prescribed neighbor counts `degree`. -/
def classTorque
    {I : Type*} [Fintype I]
    (degree : I → ℕ) (phase : I → ℝ) (i : I) : ℝ :=
  ∑ j : I, (degree j : ℝ) * Real.sin (phase j - phase i)

/-- Equitable neighbor counts make the vertex torque equal the reduced class torque exactly. -/
theorem vertexTorque_eq_classTorque
    {V I : Type*} [Fintype I] [DecidableEq I]
    (neighbors : Finset V) (cls : V → I) (degree : I → ℕ)
    (phase : I → ℝ) (i : I)
    (hcount : ∀ j : I, (neighbors.filter fun v ↦ cls v = j).card = degree j) :
    vertexTorque neighbors cls phase i = classTorque degree phase i := by
  unfold vertexTorque classTorque
  rw [sum_by_class neighbors cls (fun j ↦ Real.sin (phase j - phase i))]
  apply Finset.sum_congr rfl
  intro j _
  rw [hcount j]

/-- A zero reduced torque lifts to zero torque at the actual vertex. -/
theorem exact_equilibrium_at_vertex
    {V I : Type*} [Fintype I] [DecidableEq I]
    (neighbors : Finset V) (cls : V → I) (degree : I → ℕ)
    (phase : I → ℝ) (i : I)
    (hcount : ∀ j : I, (neighbors.filter fun v ↦ cls v = j).card = degree j)
    (hzero : classTorque degree phase i = 0) :
    vertexTorque neighbors cls phase i = 0 := by
  rw [vertexTorque_eq_classTorque neighbors cls degree phase i hcount]
  exact hzero

end


end KuramotoBlockLimit
