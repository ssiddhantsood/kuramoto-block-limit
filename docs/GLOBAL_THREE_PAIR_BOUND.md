# Toward a global bound in the three-pair reflected model

## Target

The proposed sharp statement is

```text
sup minimum degree over strictly stable, nonsynchronous,
three-reflected-pair systems
  = 0.69155376059797087278515045074165659...
```

This theorem is **not proved yet**.  The current interval calculation proves a
strict local maximum in the full nine-density model.  This document records
the additional reduction, falsification search, and rigorous ingredients for
a global attack.

## The correct five-dimensional outer problem

Let `p=(a,b,c)` be the three pair masses, with `a+b+c=1/2`, let
`theta=(alpha,beta,gamma)` be their reflection phases, and let `q` contain the
nine nontrivial block densities.  For fixed `(p,theta)`, all torque, degree,
Hessian-sector, and transverse expressions are affine in `q`.  The density
layer is therefore a convex semidefinite program.

Pair permutations allow the masses to be ordered `a>=b>=c>=0`.  Independent
swaps of the two classes in each reflection pair allow all three angles to lie
in `[0,pi]`.  The unit-square parameterization in
[`search_stability_margin.py`](../scripts/search_stability_margin.py) covers
the entire ordered mass simplex.

There is also a useful nonsynchrony reduction.  Minimum degree above `1/2`
forces the weighted block graph to be connected.  If all phases of a connected
Kuramoto equilibrium lie in an open semicircle, the maximum-phase torque has
one sign unless every neighboring phase agrees; connectivity then forces
synchrony.  For the reflected phase set `+/-theta_i`, containment in an open
semicircle is equivalent to all three `theta_i` lying strictly on the same
side of `pi/2`.  A nonsynchronous candidate must therefore satisfy

```text
min(theta) <= pi/2 <= max(theta).
```

This replaces the earlier arbitrary order-parameter cutoff in the global
search.

## Weak SDP value versus strict stability

Write `V0(p,theta)` for the optimal minimum degree in the fixed-outer density
SDP with positive-semidefinite, rather than positive-definite, stability
sectors.  Define the **inner Slater locus** to be the outer points for which at
least one density vector satisfies torque zero and makes every physical even,
odd, and transverse stability sector strictly positive.

At every point in the inner Slater locus,

```text
supremum over strictly stable densities = V0(p,theta).
```

Indeed, convex-combine a weak optimum with any strictly stable density vector.
Torque zero and the density box are preserved, every stability sector becomes
strict, and the degree tends to the weak optimum.  Thus the global problem is
equivalent to maximizing `V0` over the inner Slater locus.

The qualification by the Slater locus is essential.  The unrestricted weak
problem already admits the exact degree-`3/4` marginal family, which has three
additional zero modes and is not strictly stable.

## A direct counterexample oracle

For a proposed degree bound `mu0` and fixed outer point, solve

```text
maximize    sigma
subject to  torque(q) = 0
            degree(q) >= mu0
            E(q) >= sigma I_2
            O(q) >= sigma I_3
            transverse(q) >= sigma
            0 <= q <= 1.
```

This is again a convex SDP.  A positive optimum is immediately a strictly
stable counterexample to `mu0`.  A differential-evolution search made 6,060
outer evaluations, including 5,951 SDP solver calls, at a threshold `1e-8`
above the proposed bound:

```text
mu0 = 0.6915537705979709...
```

and found no positive stability margin.  Its best value was numerical zero on
the known marginal strata, not a strict counterexample.  Reproduce it with

```bash
.venv-sdp/bin/python scripts/search_stability_margin.py \
  --iterations 100 --population 12 \
  --out data/stability_margin_counterexample_search.json
```

The saved output is
[`data/stability_margin_counterexample_search.json`](../data/stability_margin_counterexample_search.json).

This is a stronger falsification test than maximizing degree at a prescribed
small stability floor, but it is still a non-exhaustive numerical search.

## Rigorous upper bounds on outer boxes

The fixed-outer SDP has a particularly useful dual.  Let `D_k,T_k,R_k,E_k,O_k`
be the degree, torque, transverse, even, and odd coefficients of density
`q_k`.  Choose

```text
lambda >= 0,  sum(lambda)=1,
tau >= 0,
Y >= 0, Z >= 0,
nu free,
```

and define

```text
c_k = lambda.D_k + nu.T_k + <Y,E_k> + <Z,O_k> + tau.R_k.
```

For stability floor `epsilon`, weak duality gives

```text
V_epsilon(p,theta)
 <= lambda.p + tau.(p-epsilon)
    - epsilon*(trace(Y)+trace(Z))
    + sum_k max(0,c_k).
```

The support terms `max(0,c_k)` handle the complete density box without having
to guess its active face.  A single fixed dual can be interval-evaluated over
an entire mass/phase box.

[`certify_density_sdp_box.py`](../scripts/certify_density_sdp_box.py) obtains a
useful dual numerically, then reconstructs its nonnegative weights as exact
decimals and its PSD matrices as exact-decimal Gram matrices.  Only the
reconstructed dual is used in the outward-rounded interval calculation.  On a
box around the known endpoint with radius `1e-10` in the two independent
masses and all three angles (the dependent third-mass interval has radius
`2e-10`) it gives

```text
V0 <= 0.691553763381935176238521624202.
```

That box bound is not as sharp as the specialized local KKT certificate.  Its
importance is that the same construction applies to generic boxes and is the
bounding primitive needed for a global branch-and-bound calculation.

Reproduce it with

```bash
.venv-sdp/bin/python scripts/certify_density_sdp_box.py \
  --mass-radius 1e-10 --angle-radius 1e-10
```

The corresponding certificate output is
[`data/density_sdp_box_certificate.json`](../data/density_sdp_box_certificate.json).

## The degree-`3/4` obstruction now extends through second order

For the exact marginal family, the first-order sum of the three weak-mode
openings is

```text
-4 a d(t_AB) - 4 (1/4-a) d(s_BC) <= 0.
```

Positive semidefiniteness forces both missing-density derivatives and every
individual first-order weak opening to vanish.  The complete second spectral
variation can also be calculated exactly.  After imposing the linearized
torque equations and the critical-tangent conditions, the fixed-kernel
quadratic term is

```text
[3(dA+4a dB)^2 + 32a^2 d(theta_A)^2] / [4a(1-4a)].
```

Coupling to the strong even and odd eigenvectors subtracts exactly the same
quantity.  Therefore

```text
F''(0) = -4a d2(t_AB) - 4(1/4-a) d2(s_BC) <= 0,
```

where `F` is the sum of the three weak eigenvalues.  Stability forces the
opposite inequality, so every linear and quadratic weak opening must vanish.
The symbolic verification, including the moving rotation-mode cancellation,
is

```bash
python scripts/analyze_marginal_second_order.py \
  --out data/marginal_second_order.json
```

This rules out ordinary linear and quadratic bifurcating branches.  It does
not exclude cubic, higher-order, Puiseux/nonanalytic, or arbitrary-sequence
approaches, so it is not yet a stable-free neighborhood theorem.

## What remains for a proof

1. Subdivide the compact five-dimensional ordered-mass/phase domain.
2. Use common-dual interval bounds to discard every box whose weak SDP value
   is below the candidate.
3. Use the existing full-model KKT certificate for the candidate neighborhood.
4. Prove an all-orders or interval neighborhood exclusion for the remaining
   higher-degree marginal strata, beginning with the `3/4` family.
5. If a box instead has positive certified stability margin above the target,
   extract it as a counterexample.

The dual-box and second-order calculations are concrete progress toward this
program.  The global theorem should not be claimed until all outer boxes and
all marginal strata have been resolved.
