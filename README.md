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

[`docs/CLIQUE_BLOWUP_THEOREM.md`](docs/CLIQUE_BLOWUP_THEOREM.md) proves an
additional unconditional result: every clique blow-up of a strictly stable
finite witness remains strictly stable. Applied to the `N=80,002` witness, it
produces infinitely many stable graphs above `11/16`. Its finite-dimensional
algebra is machine-checked in
[`KuramotoBlockLimit/CliqueBlowup.lean`](KuramotoBlockLimit/CliqueBlowup.lean).
The concrete hypothesis-by-hypothesis application to the `N=80,002` witness
is recorded in the witness repository's
[`CLIQUE_BLOWUP_APPLICATION.md`](https://github.com/ssiddhantsood/kuramoto-11-16-witness/blob/b4f5d8f4db1192c82f1984f80cfe964cb6dac991/records/n80002/CLIQUE_BLOWUP_APPLICATION.md),
backed by an independent 512-bit Arb certificate.

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
python3 scripts/clique_blowup_family.py
```

The script reconstructs the six-class torques, degrees, mass-orthonormal
quotient Hessian, transverse multiplication terms, and a simple extrapolation
of the candidate boundary value.

The actual reduced optimization over three classes plus their reflected
partners is implemented in
[`scripts/optimize_reflection_model.py`](scripts/optimize_reflection_model.py).
It reproduces the observed connectivity sequence while separately constraining
the physical even, odd, and transverse stability sectors. See
[`docs/REDUCED_OPTIMIZATION.md`](docs/REDUCED_OPTIMIZATION.md).

The density pattern is now also tested rather than assumed.  The general
nine-density search recovers the same two-interior-density face, finds two
lower stable branches, and exposes an exact degree-`3/4` marginal family with a
proved first-order stability obstruction.  A Krawczyk certificate establishes
the limiting KKT root, and interval bounds on all seven active density-bound
multipliers extend strict local optimality to the full nine-density model.
For fixed masses and phases, the complete density layer is also formulated as
a convex semidefinite program, leaving only five nonlinear outer variables.  See
[`docs/GENERAL_REFLECTION_SEARCH.md`](docs/GENERAL_REFLECTION_SEARCH.md).

The next target is a global bound within this entire three-pair family.  A
direct stability-margin search found no counterexample in 6,060 outer
evaluations (5,951 SDP calls), a fixed-dual interval calculation now bounds
generic outer boxes, and
the false degree-`3/4` marginal optimum is obstructed through second order.
These are proof ingredients, not yet a global theorem.  See
[`docs/GLOBAL_THREE_PAIR_BOUND.md`](docs/GLOBAL_THREE_PAIR_BOUND.md).

Build the machine-checked algebraic layer with:

```bash
lake build
```

The Lean proofs cover the degree identity, reflection reduction from six
torques to three, exact equilibrium lifting from equitable neighbor counts,
and the complete clique-blow-up theorem. For the latter, Lean checks the actual
Kuramoto Hessian matrix, its quadratic-form decomposition, preservation of
positive semidefiniteness, and preservation of the rotation-only kernel. See
[`formal/README.md`](formal/README.md).

## Proof status

- Exact finite-dimensional formulas: derived and numerically checked.
- Exact equilibrium lifting to equitable blow-ups: proved.
- Quotient reduction and a checkable finite transverse criterion: proved.
- Arbitrarily-large realization theorem: proof plan drafted; not proved yet.
- Strict local optimum in the full nine-density three-pair model:
  interval-certified.
- Sharp upper bound over the general nine-density six-block family:
  conjectural; global certificate machinery is in progress.
- General graphon extension: scoped; not proved here yet.
- Lean formalization: the complete uniform clique-blow-up result is checked,
  in addition to the degree, reflection, equitable-lifting, and quotient
  identities. The more general partial-biregular realization theorem remains
  open.
