# Lean formalization

Lean 4.34.1 is managed by Elan and this repository pins mathlib 4.34.1.
Build the checked layer with:

```bash
lake build
```

The Lean layer contains finite algebraic statements rather than the full
analytic graphon-persistence theorem:

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
  neighbor counts agree with the equitable quotient;
- symmetry of the mass-orthonormal quotient matrix; and
- the exact global-rotation kernel identity;
- exact replication of a Kuramoto equilibrium under a clique blow-up;
- equality between the constructed block matrix and the actual Kuramoto
  Hessian of the blown-up coupling matrix;
- the block-constant and block-zero-sum Hessian actions;
- the exact full quadratic-form decomposition;
- preservation of positive semidefiniteness; and
- preservation of the rotation-only Hessian kernel.

The bundled theorem is
`cliqueBlowup_preserves_stable_equilibrium` in
[`KuramotoBlockLimit/CliqueBlowup.lean`](../KuramotoBlockLimit/CliqueBlowup.lean).
It checks equilibrium, Hessian nonnegativity, and the exact Hessian kernel for
every finite base type and every positive clone count. The standard dynamical
systems implication from this spectral statement to local asymptotic
stability modulo rotation is not formalized here.

The next formal target is the finite transverse comparison inequality for
partial biregular blocks.

The probabilistic or expander existence theorem for transverse stability
should remain an ordinary proof until its exact statement and dependencies are
settled.
