/-
Copyright (c) 2026 Siddhant Sood. All rights reserved.
Authors: Siddhant Sood
-/
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Basic
import Mathlib.Tactic

set_option linter.style.header false

/-!
# Algebraic core of the reflection-paired six-block Kuramoto model

This file starts the machine-checked layer with identities used by the direct
proof. It intentionally does not encode numerical optimization claims.
-/

namespace KuramotoBlockLimit

noncomputable section

/-- The three distinct normalized class degrees in the observed six-block topology. -/
def degreeA (a b c x y : ℝ) : ℝ := a * (1 + x) + b * y + 2 * c

def degreeB (a b c y : ℝ) : ℝ := 2 * b + a * y + c

def degreeC (a b c : ℝ) : ℝ := c + 2 * a + b

/-- The third class degree is exactly `1/2 + a` after mass normalization. -/
theorem degreeC_eq_half_add_a {a b c : ℝ} (hmass : a + b + c = 1 / 2) :
    degreeC a b c = 1 / 2 + a := by
  unfold degreeC
  linarith

/-- The three independent torques for the positive members of the reflection pairs. -/
def torqueA (a b c x y α β γ : ℝ) : ℝ :=
  a * x * Real.sin (-2 * α) +
    b * y * Real.sin (β - α) +
    c * (Real.sin (γ - α) + Real.sin (-γ - α))

def torqueB (a b c y α β γ : ℝ) : ℝ :=
  a * y * Real.sin (α - β) +
    b * Real.sin (-2 * β) +
    c * Real.sin (-γ - β)

def torqueC (a b α β γ : ℝ) : ℝ :=
  a * (Real.sin (α - γ) + Real.sin (-α - γ)) +
    b * Real.sin (-β - γ)

/-- Explicit torques of the reflected partners. -/
def torqueAReflected (a b c x y α β γ : ℝ) : ℝ :=
  a * x * Real.sin (2 * α) +
    b * y * Real.sin (α - β) +
    c * (Real.sin (γ + α) + Real.sin (α - γ))

def torqueBReflected (a b c y α β γ : ℝ) : ℝ :=
  a * y * Real.sin (β - α) +
    b * Real.sin (2 * β) +
    c * Real.sin (γ + β)

def torqueCReflected (a b α β γ : ℝ) : ℝ :=
  a * (Real.sin (γ - α) + Real.sin (γ + α)) +
    b * Real.sin (γ + β)

theorem torqueA_reflection (a b c x y α β γ : ℝ) :
    torqueAReflected a b c x y α β γ = -torqueA a b c x y α β γ := by
  simp only [torqueAReflected, torqueA]
  rw [show -2 * α = -(2 * α) by ring]
  rw [show α - β = -(β - α) by ring]
  rw [show α - γ = -(γ - α) by ring]
  rw [show -γ - α = -(γ + α) by ring]
  simp only [Real.sin_neg]
  ring

theorem torqueB_reflection (a b c y α β γ : ℝ) :
    torqueBReflected a b c y α β γ = -torqueB a b c y α β γ := by
  simp only [torqueBReflected, torqueB]
  rw [show α - β = -(β - α) by ring]
  rw [show -2 * β = -(2 * β) by ring]
  rw [show -γ - β = -(γ + β) by ring]
  simp only [Real.sin_neg]
  ring

theorem torqueC_reflection (a b α β γ : ℝ) :
    torqueCReflected a b α β γ = -torqueC a b α β γ := by
  simp only [torqueCReflected, torqueC]
  rw [show -α - γ = -(γ + α) by ring]
  rw [show -β - γ = -(γ + β) by ring]
  rw [show α - γ = -(γ - α) by ring]
  simp only [Real.sin_neg]
  ring

/-- It is enough to check the three positive-pair torques. -/
theorem all_six_torques_zero
    {a b c x y α β γ : ℝ}
    (hA : torqueA a b c x y α β γ = 0)
    (hB : torqueB a b c y α β γ = 0)
    (hC : torqueC a b α β γ = 0) :
    torqueAReflected a b c x y α β γ = 0 ∧
      torqueBReflected a b c y α β γ = 0 ∧
      torqueCReflected a b α β γ = 0 := by
  constructor
  · rw [torqueA_reflection, hA, neg_zero]
  constructor
  · rw [torqueB_reflection, hB, neg_zero]
  · rw [torqueC_reflection, hC, neg_zero]

end

end KuramotoBlockLimit
