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

## Certified boundary point on the discovered face

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

This proves a strict local maximum of the closed semidefinite problem
**relative to that two-density face**.  The boundary point itself has two
additional zero Hessian modes and is not strictly stable.  The numerical
positive-margin continuation approaches it from the stable side.

Reproduce the certificate with

```bash
python scripts/certify_limiting_kkt.py \
  --out data/limiting_kkt_certificate.json
```

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
- The interval local-optimality certificate applies to the discovered
  two-density face, not to all nine density directions.
- No general three-pair upper bound has been proved.
- The limiting KKT point is marginal; strict stability belongs to nearby
  numerical continuation points, not to the endpoint.
- Turning a weighted step solution into arbitrarily large simple unweighted
  graphs still requires the separate equitable-realization argument.
