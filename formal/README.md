# Lean formalization

Lean 4.34.1 is managed by Elan and this repository pins mathlib 4.34.1.
Build the checked layer with:

```bash
lake build
```

The first Lean targets should be finite algebraic statements rather than the
full analytic persistence theorem:

1. reflection symmetry reduces six torque equations to three;
2. the displayed formulas for `dA`, `dB`, and `dC` follow from the block matrix;
3. `dC = 1/2 + a` when `a+b+c=1/2`;
4. the mass-orthonormal quotient Hessian is symmetric;
5. the square-root mass vector is in its kernel; and
6. biregularity plus the class torque equations gives zero torque at every
   vertex of a finite blow-up.

The checked layer proves:

- the degree identity;
- the three reflection-torque identities and reduction from six torque
  equations to three; and
- exact lifting of a zero class torque to an actual vertex whenever its
  neighbor counts agree with the equitable quotient.

The next formal targets are the mass-orthonormal quotient symmetry and rotation
kernel.

The probabilistic or expander existence theorem for transverse stability
should remain an ordinary proof until its exact statement and dependencies are
settled.
