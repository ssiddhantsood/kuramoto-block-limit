# Formalization plan

Lean, Lake, and Elan are not currently installed on this machine, so this
directory does not yet claim a checked formal proof.

The first Lean targets should be finite algebraic statements rather than the
full analytic persistence theorem:

1. reflection symmetry reduces six torque equations to three;
2. the displayed formulas for `dA`, `dB`, and `dC` follow from the block matrix;
3. `dC = 1/2 + a` when `a+b+c=1/2`;
4. the mass-orthonormal quotient Hessian is symmetric;
5. the square-root mass vector is in its kernel; and
6. biregularity plus the class torque equations gives zero torque at every
   vertex of a finite blow-up.

These statements fit Lean/mathlib's finite sums, matrices, and real
trigonometric functions. The probabilistic or expander existence theorem for
transverse stability should remain an ordinary proof until its exact statement
and dependencies are settled.

