# A fourth-pair invasion test with a replayable interval exclusion

## What question this answers

The certified six-class endpoint has minimum block degree

```text
mu_* = 0.69155376059797087278515045074165659...
```

and is marginal: besides global rotation, its quotient Hessian has one even
and one odd zero mode. Could adding one *new reflected pair* with vanishing
mass improve degree immediately? This is a different question from optimizing
the known six-class face. It tests whether that face is locally vulnerable to
a new phase class.

Give each new class mass `epsilon`, phase `+/-delta`, and arbitrary limiting
cross-edge densities to the old three pairs. Allow all three old pair masses,
all three angles, and all nine reflection-compatible old densities to change
at arbitrary finite first-order rates. Those rates are **not capped** in the
final exclusion. The self-density of the new pair does not enter at first
order.

The result is:

> For `delta` in `[0,pi]`, outside the three explicit `+/-0.05`-radian
> intervals centered at the decimal phase values in the certificate, there
> is no differentiable reflected four-pair branch
> through this endpoint whose minimum-degree derivative is positive and whose
> Hessian is nonnegative to first order.

The three excluded centers (absolute phase values) are
`0.7296280156946933`, `0.9565524921031598`, and `2.358744770379922`.
The full decimal centers appear in the certificate, and each encloses the
actual KKT phase to within `10^-25`.
This is a local first-order statement at one endpoint, **not** a global
optimality theorem for four-pair graphons or arbitrary graphs. It does not
exclude improvements beginning at second or higher order, paths with a
non-differentiable scaling, or branches away from this endpoint. In
particular, near-duplicate pair splitting remains the obvious open case.

## Why all four weak directions matter

For class masses `p_i`, phases `theta_i`, and block densities `W_ij`, use the
mass-squared quotient form

```text
L = sum_(i<j) p_i p_j W_ij cos(theta_i-theta_j)
                 (e_i-e_j)(e_i-e_j)^T.
```

For positive masses, this is congruent to the physical quotient Hessian and
therefore has the same positivity signs. At `epsilon=0`, the extended
eight-class matrix has a five-dimensional kernel: rotation, the two old weak
modes, and the two newly massless classes. A `4 x 4` projection of `L'(0)`
onto a basis for that kernel modulo rotation must be positive semidefinite
for a stable branch. Checking the old `2 x 2` weak block alone is invalid.

There is a concrete numerical falsification of the smaller test: at
`delta=pi`, with derivatives temporarily bounded by `10` for the diagnostic
optimization, the old `2 x 2` test reports degree-gain derivative about
`+3.87` and both old weak eigenvalues above `10^-7`. But the full projected
matrix has an eigenvalue about `-0.415`. Directly assembling the eight-class
Hessian at `epsilon=10^-6` gives a negative eigenvalue about `-4.37e-7`.
Thus the apparent improvement destabilizes a *new-class mode*.

## The one-variable semidefinite reduction

For fixed `delta`, the old first-order mass, angle, and density changes and
the six new cross densities are the only variables. The first-order old
torques, the leading new torque, the active old degree derivatives, the new
class degree, the new class transverse Hessian coefficient, and `L'(0)` are
all affine in these variables. Requiring zero torque, new degree at least
`mu_*`, nonnegative new transverse coefficient, and a positive-semidefinite
`4 x 4` weak projection gives a convex SDP. The outer search is only one
phase variable `delta`.

The SDP is a **relaxation** of a genuine stable branch. For example, if the
new class degree equals `mu_*` at `epsilon=0`, its first-order degree would
also need checking for a realizable improvement; we omit that extra
restriction. An exclusion by the relaxed SDP is therefore valid.

## Why the final dual has no derivative cap

The endpoint's existing KKT certificate includes exact equations for all
old free parameters and strict signs for all seven old densities sitting at
`0` or `1`. We use those equations in the invasion SDP's dual. Let `E` and
`O` be the three-by-three even and odd sectors of the per-row Kuramoto
Hessian, and let `v_E,v_O` be their weak null vectors. Choose `v_E` also
orthogonal to rotation in the mass inner product. At the endpoint, `E` has
rank one and `O` has rank two. For a reflected perturbation,

```text
d(e_2(E)) = trace(E) / <v_E,P v_E> * <v_E,d(P E) v_E>,
d(det(O)) = e_2(O) / <v_O,P v_O> * <v_O,d(P O) v_O>.
```

Here `e_2` is the sum of principal two-by-two minors and `P` is the class
mass matrix. These identities follow by differentiating the one nonrotation
zero eigenvalue in each sector; the rotation eigenvalue stays identically
zero. They convert the endpoint KKT multipliers into *exact* coefficients of
the old weak quadratic forms. The numerical SDP supplies only the remaining
dual choices. Its PSD dual is rebuilt as a Gram matrix; the first two Gram
row lengths are rescaled to the exact KKT coefficients. Thus every free old
derivative cancels algebraically. The seven one-sided density derivatives
have nonpositive dual contribution, with signs checked by intervals. No
bound on their magnitudes enters the upper bound.

After that cancellation, the dual objective is just a constant plus six
coefficients times new cross densities in `[0,1]`. The interval verifier
upper-bounds it over an entire `delta` box by adding each positive-part
coefficient. A negative bound excludes positive degree gain throughout that
box. The base KKT root and multipliers are enclosed by the previously
certified radius `10^-25` box; all trigonometric and dual evaluations use
outward interval arithmetic.

## Reproduce and audit

The complete replayable [certificate](../data/pair_insertion_unbounded_cover_r05.json)
contains `284` certified phase boxes. Their exact-decimal endpoints partition
all four intervals outside the three `+/-0.05` neighborhoods, with no
unresolved leaves. The largest recorded interval upper bound is below
`-0.00064065`. The solver selected duals in `564` processed boxes, but the
following replay **does not call a conic solver**: it checks the partition,
rebuilds each Gram dual, checks the KKT-consistent cancellation and boundary
signs, and recomputes each negative interval bound.

```bash
.venv/bin/python scripts/cover_pair_insertion_phase.py \
  --replay data/pair_insertion_unbounded_cover_r05.json
```

For the direct finite-difference and false-positive diagnostics:

```bash
.venv/bin/python scripts/test_pair_insertion.py \
  --grid 61 --sigma 0 --new-degree-slack 0 \
  --derivative-cap 100 --audit --out /tmp/pair-insertion-scan.json

.venv/bin/python scripts/certify_pair_insertion_box.py \
  --lower 1.55 --upper 1.59 --unbounded-derivatives
```

The first command is *numerical exploration*, not the proof. The second
demonstrates one interval certificate. The replayable cover is the actual
first-order exclusion away from duplicated phases.

## What to attack next

The remaining `0.05`-radian neighborhoods are not a matter of running a finer
grid: their dual upper bound touches zero exactly when the new pair duplicates
an old one. The next mathematical question is a second-order or
blow-up analysis of pair splitting there. A successful local exclusion could
combine with this phase cover to rule out every differentiable one-pair
insertion at the endpoint. Even that would not prove the six-class endpoint
globally optimal; disconnected branches and other phase architectures would
remain possible.
