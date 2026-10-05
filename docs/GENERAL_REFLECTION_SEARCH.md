# General three-pair reflection search

## Outcome

The original two-density pattern is no longer assumed.  The expanded model
allows every edge density compatible with reflection symmetry to vary.  Its
nine meaningful densities are

```text
vA, vB, vC,
sAB, tAB, sAC, tAC, sBC, tBC,
```

where `v` joins a class to its reflection and `s,t` join equal or opposite
signs in different pairs.  Same-class densities are fixed to one without loss
for this optimization: raising one changes neither torque nor the
block-constant Hessian, while it increases degree and transverse stability.

A 53-start search found 17 feasible endpoints.  No endpoint exceeded the
known branch.  At stability floor `1e-10`, the three leading distinct branches
were

| Branch | Minimum degree | Character |
| --- | ---: | --- |
| 0 | `0.691553757513` | known witness branch, up to relabeling |
| 1 | `0.684369721568` | distinct active density face |
| 2 | `0.662344149453` | distinct active density face |

For the leading branch, seven of the nine densities converge to zero or one
and only two remain interior.  Thus the previously observed two-density
architecture emerges from the general local search as an active face rather
than being imposed in advance.

The stored leading branch is a reflection-preserving relabeling of the
canonical face used by the interval certificate below.  The masses, phases,
and density names therefore appear in a different order, but the two points
define the same weighted six-class system.

The search is reproducible with

```bash
python scripts/search_general_reflection_model.py \
  --near-starts 12 --random-starts 40 \
  --out data/general_reflection_multistart.json
```

The search eliminates the three torque equations algebraically, then checks
the reconstructed torques and all even, odd, and transverse stability sectors.
The independent consistency audit is

```bash
python scripts/audit_general_reflection_model.py
```

It checks reflection identities, the rotation kernel, agreement between the
full six-block spectrum and the even/odd sectors, and a finite-difference
Jacobian/Hessian identity.

## Certified local maximum in the full nine-density model

On the leading two-density face, the zero-margin endpoint satisfies a
17-dimensional KKT system: nine primal variables and eight multipliers.  The
certificate proves a unique KKT root inside a radius-`1e-25` box around

```text
mu = 0.69155376059797087278515045074165659...
```

Using outward-rounded interval arithmetic, it also verifies:

- strict Krawczyk inclusion;
- linear independence of the eight active constraint gradients;
- the required strict signs of all four active inequality multipliers;
- positive reduced Lagrangian curvature on the one-dimensional critical
  tangent space; and
- positive interval lower bounds for every inactive mass, density, spectral,
  transverse, nonsynchrony, objective, and angle constraint.

The seven density coordinates fixed at zero or one on this face require an
additional check before the result says anything about the general model.
The certificate now computes their seven KKT multipliers directly from the
full nine-density equations and proves that every multiplier has the required
strict sign.  The smallest sign margin is the lower-bound multiplier for
`ref_C`, whose entire certified interval lies below `-0.0035079204`.

These strict signs have two consequences.  First, the seven density-bound
gradients form a coordinate block; after eliminating them, the remaining
LICQ minor is exactly the already certified face minor.  Second, every
critical direction must be tangent to those seven bounds.  The full critical
cone therefore reduces to the same one-dimensional face tangent on which the
reduced Lagrangian curvature is interval-certified positive.

This proves a strict local maximum of the closed weak-stability problem in the
**full nine-density three-pair model**, not merely relative to the discovered
face.  The boundary point itself has two additional zero Hessian modes and is
not strictly stable.  The numerical positive-margin continuation approaches
it from the stable side.

Reproduce the certificate with

```bash
python scripts/certify_limiting_kkt.py \
  --out data/limiting_kkt_certificate.json
```

## A second reduction: five outer variables plus a convex SDP

There is a useful reduction that does not assume an active density face.  If
the three pair masses and three phases are fixed (two independent mass
coordinates and three angles), then all of the following are affine in the
nine densities:

- the three torque equations;
- the three degree inequalities;
- the physical even and odd Hessian blocks; and
- the three transverse stability values.

Consequently, maximizing the minimum degree over all nine densities is a
small convex semidefinite program.  The only nonlinear search is over the five
outer mass/phase variables.  This is the general reduced problem suggested by
the observed three reflected pairs.

Writing `q` for the nine-density vector and `epsilon` for the requested
stability floor, the inner problem has the form

```text
maximize    mu
subject to  T q = 0                         (equilibrium)
            d0 + D q >= mu                  (minimum degree)
            E(q) >= epsilon I_2             (physical even sector)
            O(q) >= epsilon I_3             (odd sector)
            r0 + R q >= epsilon             (transverse sectors)
            0 <= q <= 1.
```

Here `E(q)` and `O(q)` are affine symmetric matrices, and `>=` on those two
lines means positive-semidefinite ordering.  Every constraint is therefore
linear or semidefinite in `(q,mu)`.  The outer problem varies two independent
mass coordinates and three phase angles.

Install the optional solver and reproduce the fixed-outer and outer-search
experiments with

```bash
python3 -m venv .venv-sdp
.venv-sdp/bin/pip install -r requirements-sdp.txt
.venv-sdp/bin/python scripts/audit_density_sdp.py
.venv-sdp/bin/python scripts/solve_density_sdp.py --case known
.venv-sdp/bin/python scripts/solve_density_sdp.py --case marginal
.venv-sdp/bin/python scripts/solve_density_sdp.py \
  --stability-floor 1e-6 --search-iterations 80 \
  --out data/sdp_outer_search.json
```

At the certified outer point the SDP recovers the known density face and
`mu=0.69155376059...` at zero stability floor.  At a positive floor of
`1e-6`, the saved differential-evolution run evaluated 3,240 inner SDPs, of
which 1,360 were feasible, and found no point above the known continuation
branch.  The known branch was included in the initial population, so this is
a counterexample search rather than an independent global-optimality proof.
Clarabel supplies approximate floating-point SDP solutions; the saved result
is reconstructed in the independent full evaluator, but it has no rigorous
SDP duality-gap certificate.

The randomized algebra audit compares every affine sector used by the SDP
with the direct six-block evaluator on 100 sampled models.  Its largest
absolute discrepancy in the saved implementation is below `4e-16`.

## The degree-`3/4` marginal family

The general model also contains a more connected exact family.  For any
`0<a<1/4`, take pair masses

```text
(a, 1/4, 1/4-a),
```

phases

```text
(-pi/4, -3*pi/4, pi/4),
```

and densities

```text
(vA,vB,vC,sAB,tAB,sAC,tAC,sBC,tBC)
  = (1,1,1,1,0,1,1,0,1).
```

Every class has degree `3/4`, every torque is zero, and every transverse value
is `1/4`.  Nevertheless, both the even and odd characteristic polynomials are

```text
lambda^2 (lambda-1/4).
```

Besides rotation, there are therefore three additional zero quotient modes.
This is a marginal family, not a stable degree-`3/4` construction.

There is also an exact first-order obstruction.  If `tAB` and `sBC` denote the
two densities that vanish on this family, the sum of the first variations of
the one weak even mode and the two weak odd modes is

```text
-4 a d(tAB) - 4 (1/4-a) d(sBC).
```

Feasible one-sided density variations satisfy `d(tAB)>=0` and `d(sBC)>=0`.
Consequently, no first-order perturbation can make all three modes strictly
positive.  The symbolic verification is

```bash
python scripts/analyze_marginal_family.py
```

This obstruction does not rule out a branch whose eigenvalue openings begin
only at second or higher order.  Numerical probes found none, but that remains
an open question.

## Limits of the result

- Multistart local optimization is not a global search.
- The interval certificate proves local optimality in all nine density
  directions, but it does not prove global optimality.
- No general three-pair upper bound has been proved.
- The limiting KKT point is marginal; strict stability belongs to nearby
  numerical continuation points, not to the endpoint.
- Turning a weighted step solution into arbitrarily large simple unweighted
  graphs still requires the separate equitable-realization argument.
