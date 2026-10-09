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

A separate [random-block transfer theorem](docs/RANDOM_BLOCK_TRANSFER_THEOREM.md)
now proves that **every fixed, strictly stable finite-step block equilibrium**
in the stated sense produces stable nonsynchronous equilibria on sufficiently
large simple unweighted graphs, with normalized minimum degree approaching the
block value. This is a direct asymptotic argument, not an exact equitable lift
and not a claim that the published graphon persistence theorem applies
verbatim to discontinuous step graphons.

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

On the observed active-degree branch, the three torque equations and the
degree tie eliminate both mass ratios and both partial densities. The local
search then uses **only three phase angles**. An exact phase-defined strict
block point is interval-checked with minimum-degree ratio at least
`0.691553760567123749076110`; the transfer theorem turns it into an
asymptotic finite-simple-graph existence result. The derivation, verifier,
reproduction commands, and limitations are in
[`docs/THREE_PHASE_BRANCH.md`](docs/THREE_PHASE_BRANCH.md). An independent
[Arb-ball audit](data/three_phase_arb_certificate.json) checks the same strict
point through a separate arithmetic backend.

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
evaluations (5,951 SDP calls), and a fixed-dual interval branch-and-bound now
certifies individual boxes in the complete five-dimensional outer domain.
The first run is only a partial cover, not a global certificate.  The false
degree-`3/4` optimum has also been identified as part of an exact
six-dimensional marginal sheet: all second-flat directions are tangent to
that sheet, while the all-orders transverse exclusion remains open.  These
are proof ingredients, not yet a global theorem.  See
[`docs/GLOBAL_THREE_PAIR_BOUND.md`](docs/GLOBAL_THREE_PAIR_BOUND.md).
The exact phase triple `(-pi/4,-3pi/4,pi/4)` is also completely classified:
all weakly stable points on that slice have degree at most `3/4` and extra
zero modes. This useful exclusion is only a fixed-angle result; see
[`docs/FIXED_ANGLE_SLICE.md`](docs/FIXED_ANGLE_SLICE.md).

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
- Arbitrarily-large **random-block** realization for each fixed strict step
  equilibrium: proved, with phase corrections after sampling. The stronger
  exact equitable partial-biregular realization target remains open.
- Strict local optimum in the full nine-density three-pair model:
  interval-certified.
- Sharp upper bound over the general nine-density six-block family:
  conjectural; global certificate machinery is in progress.
- General graphon extension: scoped; not proved here yet.
- Lean formalization: the complete uniform clique-blow-up result is checked,
  in addition to the degree, reflection, equitable-lifting, and quotient
  identities. The more general partial-biregular realization theorem remains
  open.
