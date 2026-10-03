/-
Copyright (c) 2026 Siddhant Sood. All rights reserved.
Authors: Siddhant Sood
-/
import KuramotoBlockLimit.Equitable
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv

set_option linter.style.header false

/-!
# Mass-orthonormal quotient operator

The entries use square-root class masses `r i`, so the actual mass is
`(r i)^2`. This avoids introducing divisions or square-root side conditions in
the algebraic kernel proof.
-/

namespace KuramotoBlockLimit

noncomputable section

open scoped BigOperators

/-- Symmetric matrix of the class-constant Hessian in mass-orthonormal coordinates. -/
def quotientMatrix
    {I : Type*} [Fintype I] [DecidableEq I]
    (r : I → ℝ) (weight : I → I → ℝ) : Matrix I I ℝ :=
  fun i j ↦
    if i = j then
      ∑ k ∈ Finset.univ.erase i, (r k) ^ 2 * weight i k
    else
      -r i * r j * weight i j

/-- Symmetry of the block weights implies symmetry of the quotient Hessian. -/
theorem quotientMatrix_apply_symm
    {I : Type*} [Fintype I] [DecidableEq I]
    (r : I → ℝ) (weight : I → I → ℝ)
    (hweight : ∀ i j, weight i j = weight j i) (i j : I) :
    quotientMatrix r weight i j = quotientMatrix r weight j i := by
  by_cases hij : i = j
  · subst j
    rfl
  · simp only [quotientMatrix, hij, Ne.symm hij, ↓reduceIte]
    rw [hweight i j]
    ring

theorem quotientMatrix_isSymm
    {I : Type*} [Fintype I] [DecidableEq I]
    (r : I → ℝ) (weight : I → I → ℝ)
    (hweight : ∀ i j, weight i j = weight j i) :
    (quotientMatrix r weight).IsSymm := by
  apply Matrix.IsSymm.ext
  intro i j
  exact quotientMatrix_apply_symm r weight hweight j i

/-- Coordinate-free row action of the same mass-orthonormal quotient Hessian. -/
def quotientAction
    {I : Type*} [Fintype I]
    (r : I → ℝ) (weight : I → I → ℝ) (q : I → ℝ) (i : I) : ℝ :=
  ∑ j : I, r j * weight i j * (r j * q i - r i * q j)

/-- The square-root mass vector is the exact global-rotation kernel vector. -/
theorem quotientAction_rotation
    {I : Type*} [Fintype I]
    (r : I → ℝ) (weight : I → I → ℝ) (i : I) :
    quotientAction r weight r i = 0 := by
  unfold quotientAction
  apply Finset.sum_eq_zero
  intro j _
  ring

end


end KuramotoBlockLimit
