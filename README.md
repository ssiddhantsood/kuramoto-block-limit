# Kuramoto block-limit program

This repository studies the small limiting problem suggested by the large
finite witnesses in
[`kuramoto-11-16-witness`](https://github.com/ssiddhantsood/kuramoto-11-16-witness).
It is intentionally separate from the witness repository: numerical searches,
new conjectures, and incomplete proofs live here rather than beside certified
finite records.

## Research tracks

1. **Exact block/blow-up theorem (current priority).** Reduce the observed
   reflection-paired six-class family to a finite-dimensional problem and prove
   that a strictly stable weighted block solution produces arbitrarily large
   stable simple unweighted graphs.
2. **Finite-step graphon persistence.** Extend the persistence framework of
   Bramburger--Holzer--Williams to finite-step graphons and piecewise-constant
   equilibria, then compare that theorem with the direct blow-up theorem.

The detailed success criteria and non-goals are in [`GOALS.md`](GOALS.md).
The first exact results are in
[`docs/DIRECT_LIFTING_THEOREM.md`](docs/DIRECT_LIFTING_THEOREM.md): exact
equilibrium lifting, the symmetric quotient formula, and a rigorous transverse
comparison criterion.

## Current evidence

The five ultrasmall-gap records converge toward a six-block candidate with
minimum-degree ratio near

```text
0.691553760598
```

while one even and one odd nonrotation Hessian eigenvalue approach zero
together. This is evidence for a boundary optimizer, not yet a proof.

Run the dependency-free reconstruction with:

```bash
python3 scripts/analyze_candidate.py
```

The script reconstructs the six-class torques, degrees, mass-orthonormal
quotient Hessian, transverse multiplication terms, and a simple extrapolation
of the candidate boundary value.

Build the machine-checked algebraic layer with:

```bash
lake build
```

The Lean proofs currently cover the degree identity, reflection reduction from
six torques to three, and exact equilibrium lifting from equitable neighbor
counts. See [`formal/README.md`](formal/README.md).

## Proof status

- Exact finite-dimensional formulas: derived and numerically checked.
- Exact equilibrium lifting to equitable blow-ups: proved.
- Quotient reduction and a checkable finite transverse criterion: proved.
- Arbitrarily-large realization theorem: proof plan drafted; not proved yet.
- Sharp upper bound within the six-block topology: conjectural.
- General graphon extension: scoped; not proved here yet.
- Lean formalization: active. Lean 4.34.1/mathlib is pinned; the degree
  identity, reflection-torque reduction, and exact finite equilibrium lifting
  are machine-checked.
