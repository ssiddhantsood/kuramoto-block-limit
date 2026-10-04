# Optimizing three classes plus their reflections

This is the finite-dimensional experiment suggested by the repeated
six-cluster architecture. It replaces a search over a huge finite graph by an
optimization over three class masses and their reflected partners.

## Problem

The classes are

```text
A+, A-, B+, B-, C+, C-.
```

Reflection fixes the masses and phases as

```text
(a,a,b,b,c,c),             a+b+c=1/2,
(alpha,-alpha,beta,-beta,gamma,-gamma).
```

For the observed block topology, only two edge densities `x,y` remain free.
The numerical problem maximizes

```text
mu = min(dA,dB,dC)
```

subject to:

- the three independent torque equations;
- every degree being at least `mu`;
- the smallest physical even Hessian eigenvalue being at least `epsilon`;
- the smallest odd Hessian eigenvalue being at least `epsilon`;
- every transverse multiplication value being at least `epsilon`; and
- a nonsynchrony constraint on the order parameter.

The even rotation vector is projected out before imposing the even-sector
constraint. Without this projection an optimizer can incorrectly treat the
symmetry zero as a stability eigenvalue. Without the nonsynchrony constraint,
the optimizer can return the irrelevant synchronous equilibrium.

## Reproduce

```bash
python3 -m venv .venv-opt
.venv-opt/bin/pip install -r requirements-optimization.txt
.venv-opt/bin/python scripts/optimize_reflection_model.py \
  --out data/reduced_optimization_results.json
```

The optimizer uses continuation: the solution at one strict stability floor
is the initial point at the next smaller floor.

## Result

| Stability floor | Optimized `mu` | Weak even | Weak odd | Minimum transverse |
| ---: | ---: | ---: | ---: | ---: |
| `1e-5` | `0.691249589628013` | `1.0000e-5` | `1.0000e-5` | `0.2402587832` |
| `1e-7` | `0.691550676325926` | `1.0000e-7` | `1.0000e-7` | `0.2402800844` |
| `1e-9` | `0.691553729750870` | `1.0000e-9` | `1.0000e-9` | `0.2402803035` |
| `1e-11` | `0.691553760289604` | `9.9975e-12` | `9.9940e-12` | `0.2402803028` |

A two-point extrapolation to zero stability floor gives

```text
mu_* = 0.6915537605980755.
```

At the final continuation point,

```text
(a,b,c) = (0.191553760289605,
           0.236733060574374,
           0.071713179136021),

(x,y)   = (0.788448504093836,
           0.868249518995658),

(dA,dB,dC) = (0.691553760289606,
               0.711495760518028,
               0.691553760289606).
```

Thus `dA=dC` are the active degree constraints, while `dB` remains slack.
Both the physical even and odd stability constraints become active together;
the transverse constraints remain far from zero. This gives a concrete
explanation for the two weak eigenvalues in the large finite witnesses: the
search is approaching a reduced optimum at the intersection of two stability
boundaries.

## What this does and does not establish

This is a real reduction in computational complexity: one continuation run
takes less than a second on a laptop and reproduces the large-graph sequence.
It is still only a local numerical result:

- it starts from the observed two-density topology;
- SLSQP does not prove global optimality;
- the zero-floor value is extrapolated rather than interval-certified; and
- it does not rule out another three-pair topology or a model with more
  reflection pairs having a larger value.

The next mathematical targets are an interval enclosure of the limiting
Karush--Kuhn--Tucker system and a search over more general
reflection-symmetric block-density patterns.
