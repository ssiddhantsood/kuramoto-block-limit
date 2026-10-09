# Strict block equilibria give large stable simple graphs

This note supplies a direct finite-graph transfer argument for a **fixed,
strictly stable** finite-step Kuramoto solution. It applies to the exact
six-block interior point in
[`three_phase_interior_certificate.json`](../data/three_phase_interior_certificate.json)
and the stronger but more delicate point in
[`three_phase_interior_near_boundary.json`](../data/three_phase_interior_near_boundary.json).
The result is asymptotic; it does not give a useful explicit graph size.

In plain language, sampling the block model gives a finite graph whose copied
class phases are *almost* an equilibrium. Strict stability lets us correct
those phases slightly to an *exact* equilibrium; concentration of the random
edges ensures that the corrected equilibrium stays stable and retains nearly
the block model's minimum degree. The finite phases need not be exactly
constant within each class.

## The block hypotheses

Take finitely many classes with positive masses `p_s` summing to one, a
symmetric matrix `W_st` in `[0,1]`, and class phases `phi_s` that are not all
equal modulo rotation. Assume the block torque equations

```text
sum_t p_t W_st sin(phi_t-phi_s) = 0                         (1)
```

for every class. Define its normalized degree and transverse coefficient by

```text
d_s = sum_t p_t W_st,
q_s = sum_t p_t W_st cos(phi_t-phi_s).
```

The mass-orthonormal class-constant Hessian `Q` has off-diagonal entries
`Q_st=-sqrt(p_s p_t) W_st cos(phi_t-phi_s)` and diagonal entries
`Q_ss=sum_{t≠s} p_t W_st cos(phi_t-phi_s)`. Suppose `Q` is positive on the
orthogonal complement of its rotation vector `(sqrt(p_s))`, and every
`q_s>0`. These are exactly the quotient and transverse strict-stability
conditions of the step graphon.

**Theorem.** For all sufficiently large `n`, choose integer class sizes
`n_s` with `sum_s n_s=n` and `|n_s-n p_s|` bounded independently of `n`.
Independently include each edge between distinct vertices in classes `s,t`
with probability `W_st`. With probability tending to one, the resulting
simple graph has an exact nonsynchronous Kuramoto equilibrium whose phase at
each vertex differs from its class phase by
`O(sqrt(log n/n))` after fixing global rotation. Its actual `n`-vertex
Hessian is positive on the complement of rotation, so the equilibrium is
locally asymptotically stable modulo global rotation, and

```text
minimum_degree(G_n)/n = min_s d_s + O(sqrt(log n/n)).       (2)
```

In particular, every fixed target `tau<min_s d_s` is exceeded by the
normalized minimum degree of some stable finite simple graphs, for every
sufficiently large size. Equation (2) does not assert exact attainment of
`min_s d_s` at a finite `n`.

## Deterministic Newton lemma

For a finite symmetric adjacency matrix `A`, write the normalized Kuramoto
vector field and its potential Hessian at phase vector `theta` as

```text
F_i(theta) = (1/n) sum_j A_ij sin(theta_j-theta_i),
H(theta) = -DF(theta).
```

The field `F` lies in the rotation complement `V={h:sum_i h_i=0}`; `H`
annihilates the rotation vector and maps `V` into itself. At a template
phase vector `phi`, suppose

```text
H(phi)|_V >= g I,
q_i^G=(1/n) sum_j A_ij cos(phi_j-phi_i) >= q>0,
eta=||F(phi)||_infinity,
C=(1+1/g)/q.
```

If `8 C^2 eta<1` and `8 C eta<g`, there is an exact equilibrium `phi+h`
with `h∈V`, `||h||_infinity<=2 C eta`, and
`H(phi+h)|_V >= (g-8 C eta) I>0`.

To see the dimension-independent inverse bound, solve `H(phi)h=r∈V`.
The spectral gap gives `||h||_2/sqrt(n)<=||r||_infinity/g`. Writing the
`i`th equation separately, using `q_i^G>=q`, `|cos|<=1`, and
Cauchy–Schwarz gives

```text
q_i^G |h_i| <= |r_i| + (1/n) sum_j A_ij |h_j|
          <= (1+1/g) ||r||_infinity.
```

Thus `||H(phi)|_V^{-1}||_(infinity→infinity)<=C`. On a sup-norm ball of radius
`rho`, the Taylor remainder
`R(h)=F(phi+h)-F(phi)+H(phi)h` also lies in `V` and obeys
`||R(h)||_infinity<=2rho^2` and has Lipschitz constant at most `4rho`.
The map `h↦H(phi)|_V^{-1}(F(phi)+R(h))` contracts the radius
`rho=2C eta` ball under the first inequality above. Changing phases by `h`
changes the symmetric Hessian in operator norm by at most
`4||h||_infinity`, so the second inequality preserves its gap.

## Why the sampled block graph satisfies the lemma

Let `p_s^(n)=n_s/n`. The expectation of `H(phi)` splits exactly into the
class-constant quotient with masses `p_s^(n)` and, for each class `s`, a
within-class zero-sum eigenspace with eigenvalue

```text
q_s^(n)=sum_t p_t^(n) W_st cos(phi_t-phi_s).
```

Indeed, the omitted self-loop subtracts `W_ss/n` from the expected diagonal
at a vertex of class `s`. For a vector with zero sum within that class, the
sum of the off-diagonal within-class terms adds `W_ss/n` back. The two cancel,
giving exactly the displayed `q_s^(n)`. Since
`p_s^(n)=p_s+O(1/n)`, the expected Hessian has a fixed positive gap on `V`
for large `n`, and its expected signed cosine row sums tend to the positive
`q_s`. The block equilibrium equation (1) gives
`||E F(phi)||_infinity=O(1/n)`.

Hoeffding's inequality and a union bound over vertices give, with probability
`1-O(n^-7)`, deviations at most `2 sqrt(log n/n)` simultaneously for all
torque rows, cosine row sums, and normalized degrees. Each summand has range
length at most `1/n`.

For the Hessian, one random edge `e={i,j}` contributes a centered matrix

```text
X_e=(A_ij-W_st) cos(phi_j-phi_i)
    (e_i-e_j)(e_i-e_j)^T/n.
```

It has norm at most `2/n`, and the matrix variance satisfies

```text
||sum_e E X_e^2|| <= 1/(2n).
```

For the latter, `W_st(1-W_st)<=1/4` bounds each edge's variance by
`(e_i-e_j)(e_i-e_j)^T/(2n^2)` in positive-semidefinite order, and the sum of
these rank-one matrices over all unordered vertex pairs is
`n I - 1 1^T`, whose norm is `n`.

The standard self-adjoint matrix Bernstein inequality therefore yields

```text
P(||H-E H|| >= t) <= 2n exp[-n t^2/(1+4t/3)].
```

At `t=3 sqrt(log n/n)` this failure probability is `O(n^-8)`.
Consequently the sampled Hessian retains a fixed positive gap, its cosine
row sums remain positive, and `eta=O(sqrt(log n/n))`; the deterministic lemma
applies for all sufficiently large `n`. Degree concentration gives (2).

## Application and relation to the graphon paper

[`certify_three_phase_interior.py`](../scripts/certify_three_phase_interior.py)
defines an exact six-block solution from three rational-decimal phase angles
and interval-checks all the strict hypotheses. Its minimum normalized degree
lies between `0.688826561038580559095455` and
`0.688826561038580559095456`, above `11/16=0.6875`.
The [near-boundary strict point](../data/three_phase_interior_near_boundary.json)
has minimum degree in
`[0.691553760567123749076110, 0.691553760567123749076111]`.
An [independent Arb-ball calculation](../data/three_phase_arb_certificate.json)
checks the same strict point and proves the displayed lower threshold.
The theorem therefore yields arbitrarily large finite simple graphs with
stable nonsynchronous equilibria and normalized minimum degree above any
fixed threshold less than `0.691553760567123749076110`. This includes
`11/16`; the less delicate first point suffices for that threshold.
It does not assert the limiting candidate `0.69155376059797...` is attained:
that candidate has extra zero modes. Strict interior points may approach it,
but the required graph size is not uniform as their spectral gaps shrink.

Bramburger, Holzer, and Williams prove continuum-to-finite persistence under
their regularity and invertibility hypotheses
([arXiv:2402.09276, Theorems 3.1–3.2](https://arxiv.org/html/2402.09276v1#S3.SS2)).
Their [Section 6.4](https://arxiv.org/html/2402.09276v1#S6.SS4) sketches the
partition-adapted extension needed for discontinuous block graphons and
piecewise-constant equilibria. The argument above proves the transfer used
here directly for finite-step Kuramoto models, without asserting that the
published theorem applies verbatim. It also fixes the unavoidable global
rotation zero mode by working throughout on `V`.
