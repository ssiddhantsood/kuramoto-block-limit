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

## A full-domain branch-and-bound now exists, but the cover is partial

[`certify_stability_margin_domain.py`](../scripts/certify_stability_margin_domain.py)
applies the stability-margin dual to exact rational boxes in all five outer
coordinates `(u,v,w_A,w_B,w_C)`, with `theta_i=pi*w_i`.  A numerical SDP is
used only to select multipliers.  Before a box can be discarded, the script
reconstructs nonnegative multipliers as exact rationals, reconstructs the PSD
dual matrices as exact rational Gram matrices, imposes the dual normalization
as a rational identity, and reevaluates the bound with outward interval
arithmetic.  It discards a box only when the resulting upper bound on the
common stability margin is nonpositive.

At the test threshold

```text
mu = mu_star + 10^-8,
```

the first 600-node depth-first run produced

```text
470 interval-dual boxes discarded,
116 open-semicircle boxes rejected,
15 frontier boxes unresolved,
0 dual-selector failures.
```

The exact rational volumes of the three leaf sets sum to one.  This accounting
also exposes the limitation of the run: the unresolved leaves still occupy
`0.99860858917236328125` of the original parameter-cube volume.  The small
frontier count is therefore not evidence of nearly complete coverage; several
frontier boxes are still coarse.  This is a functioning rigorous bounding
primitive and a partial cover, **not** a proof of the global bound.

Reproduce the saved run with

```bash
python scripts/certify_stability_margin_domain.py \
  --max-nodes 600 --max-depth 35 --dual-cap 1 \
  --out data/stability_margin_domain_partial.json
```

The saved JSON is a run summary.  It does not serialize every rational Gram
factor selected during the run, so it is not a standalone solver-free proof
object; independent replay currently requires rerunning the script.

## The degree-`3/4` point lies on an exact six-dimensional marginal sheet

The binary degree-`3/4` family is not an isolated marginal construction.
Write the pair masses as `A,B,C>0`, with `A+B+C=1/2`, fix the phases

```text
(-pi/4, -3*pi/4, pi/4),
```

and set the two missing densities to `t_AB=s_BC=0`.  The four densities
`s_AB,s_AC,t_AC,t_BC` are free (subject to the density box), and the three
reflected-pair densities are

```text
v_A = (B*s_AB-C*s_AC)/A,
v_B = (A*s_AB+C*t_BC)/B,
v_C = (B*t_BC-A*s_AC)/C.
```

Whenever these values lie in the density box, all three torques vanish
identically.  The quotient stability matrices reduce exactly to

```text
E = t_AC [[C,0,-sqrt(A*C)],[0,0,0],[-sqrt(A*C),0,A]],
O = t_AC [[C,0,+sqrt(A*C)],[0,0,0],[+sqrt(A*C),0,A]].
```

Thus both matrices are rank one and positive semidefinite for `t_AC>=0`.
Each has strong eigenvalue `t_AC*(A+C)` and two zero eigenvalues.  This is an
exact six-dimensional marginal sheet with three extra nonrotation zero modes
(one even and two odd), in addition to the structural rotation zero.  Its
coordinates are
`(A,B,s_AB,s_AC,t_AC,t_BC)`.  The degree-`3/4` family is its boundary subset

```text
B=1/4,
(v_A,v_B,v_C,s_AB,s_AC,t_AC,t_BC)=(1,1,1,1,1,1,1).
```

The sheet explains why some apparently new higher-order directions remain
flat.  For example, keep the masses fixed, decrease every existing density at
rate `1` except `t_AC`, which decreases at rate `1/2`, and leave the two
missing densities zero.  This path stays on the sheet exactly: all weak
eigenvalues vanish to every order and the strong eigenvalue is `(1-t/2)/4`.

The exact verification is

```bash
python scripts/analyze_marginal_manifold.py \
  --out data/marginal_manifold.json
```

Its saved output is
[`data/marginal_manifold.json`](../data/marginal_manifold.json).

## First- and second-jet obstructions transverse to the sheet

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

The larger marginal-sheet calculation sharpens the description of the
second-flat directions.  After the first-order conditions, every critical
first jet is parameterized by

```text
dA, dB, d(s_AB), d(s_AC), d(t_AC), d(t_BC), d(theta_A),
```

with `d(theta_B)=0`, `(1/4-a)d(theta_C)=a d(theta_A)`, and the torque equations
determining `d(v_A),d(v_B),d(v_C)`.  Requiring a second-flat torque and weak
Schur-complement jet gives a rank-five linear system.  Its four exact
cokernel compatibility equations are

```text
0,
0,
4a*d(theta_A)^2/(4a-1),
0.
```

Since `0<a<1/4`, compatibility forces all three phase derivatives to vanish.
The six remaining directions are **exactly** the tangent space of the
six-dimensional marginal sheet above.  In particular, they are not limited
to varying `a` and uniformly scaling the existing densities.

This is an exact second-jet statement, not yet an all-orders exclusion.  A
promising proposed continuation is to eliminate the reflected densities by
the torque equations, use the explicit sheet as local tangent coordinates,
and apply analytic Lyapunov--Schmidt reduction to the one excess phase
direction.  The nonzero quadratic compatibility coefficient should force the
first nonzero phase coefficient of an analytic or Puiseux arc to vanish.
What remains to be checked is the valuation step when lower-order tangent
motion along the sheet is present, together with the one-sided orders of the
two missing densities.  Until that bookkeeping is complete, curve selection
does not supply a proved punctured-neighborhood theorem.

## What remains for a proof

1. Subdivide the compact five-dimensional ordered-mass/phase domain.
2. Use common-dual interval bounds to discard every box whose weak SDP value
   is below the candidate.
3. Use the existing full-model KKT certificate for the candidate neighborhood.
4. Complete the Lyapunov--Schmidt/valuation argument, or an interval
   neighborhood exclusion, transverse to the exact six-dimensional marginal
   sheet through the `3/4` family.
5. If a box instead has positive certified stability margin above the target,
   extract it as a counterexample.

The dual-box and second-order calculations are concrete progress toward this
program.  The global theorem should not be claimed until all outer boxes and
all marginal strata have been resolved.
