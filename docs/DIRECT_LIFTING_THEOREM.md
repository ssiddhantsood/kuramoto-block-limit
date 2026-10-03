# Exact lifting and stability reduction for equitable blow-ups

This note proves the finite algebraic part of Track A. It does **not** yet prove
that suitable partial biregular blocks exist for every sufficiently large set
of class sizes; that is the remaining realization theorem.

## Setup

Let a finite simple undirected graph have an equitable partition

```text
V = V_1 disjoint-union ... disjoint-union V_k,   |V_i| = n_i.
```

For every `i,j`, each vertex of `V_i` has exactly `d_ij` neighbors in `V_j`.
Undirectedness implies the edge-balance identities

```text
n_i d_ij = n_j d_ji.
```

Assign a common phase `phi_i` to every vertex in `V_i` and write

```text
c_ij = cos(phi_i-phi_j).
```

The irrelevant positive global factor in the Kuramoto vector field is omitted.

## Theorem 1 — exact equilibrium lifting

If

```text
sum_j d_ij sin(phi_j-phi_i) = 0
```

for every class `i`, then the class-constant phase assignment is an exact
equilibrium of the Kuramoto system on the full finite graph.

### Proof

Fix a vertex `v` in `V_i`. Equitability says that `v` has exactly `d_ij`
neighbors in `V_j`, and every such neighbor has phase `phi_j`. Its torque is
therefore

```text
sum_{w adjacent to v} sin(theta_w-theta_v)
  = sum_j d_ij sin(phi_j-phi_i)
  = 0.
```

This holds for every vertex. No approximation or limiting argument is used.

## Theorem 2 — mass-orthonormal quotient Hessian

Let `H` be the potential Hessian at the class-constant equilibrium:

```text
(H z)_v = sum_{w adjacent to v} cos(theta_v-theta_w)(z_v-z_w).
```

The class-constant subspace is invariant. In ordinary class coordinates its
matrix `L` is

```text
L_ii = sum_{j != i} d_ij c_ij,
L_ij = -d_ij c_ij,  i != j.
```

After the coordinate change `q_i=sqrt(n_i) z_i`, the quotient becomes the
symmetric matrix

```text
S = diag(sqrt(n_i)) L diag(1/sqrt(n_i)),
S_ij = -sqrt(d_ij d_ji)c_ij,  i != j.
```

Its rotation vector is `(sqrt(n_1),...,sqrt(n_k))`.

### Proof

For `z` constant on each class, every vertex in `V_i` has

```text
(H z)_v = sum_j d_ij c_ij(z_i-z_j),
```

which proves invariance and gives `L`. For `i != j`, edge balance gives

```text
sqrt(n_i/n_j)d_ij = sqrt(d_ij d_ji).
```

Consequently `S_ij=-sqrt(d_ij d_ji)c_ij=S_ji`, so `S` is symmetric. The
constant vertex vector is the global rotation mode of `H`; under the coordinate
change it becomes the displayed square-root size vector.

After division by the total vertex count and passage to proportions
`p_i=n_i/N` and densities `w_ij=d_ij/n_j`, this is precisely the
mass-orthonormal step-block formula

```text
S_ij/N = -sqrt(p_i p_j)w_ij cos(phi_i-phi_j).
```

## Theorem 3 — finite transverse comparison criterion

Assume every induced graph on `V_i` is a clique. For `i != j`, let `B_ij` be
the `n_i`-by-`n_j` biadjacency matrix of the block and define its centered part

```text
Btilde_ij = B_ij - (d_ij/n_j) J.
```

Let `sigma_ij` be any upper bound for the Euclidean operator norm of
`Btilde_ij`. Define the symmetric `k`-by-`k` comparison matrix `C` by

```text
C_ii = n_i + sum_{j != i} d_ij c_ij,
C_ij = -|c_ij| sigma_ij,  i != j.
```

For an absent or complete block one may take `sigma_ij=0`. If `C` is positive
definite, then the full Hessian is positive on the subspace

```text
Z = {z : sum_{v in V_i} z_v = 0 for every i}.
```

Together with positivity of the quotient away from rotation, this proves that
the full Hessian has exactly the rotation zero mode.

### Proof

Because `H` is symmetric and the class-constant subspace is invariant, its
orthogonal complement `Z` is invariant as well. Write `z_i` for the restriction
of `z` to `V_i` and `r_i=||z_i||_2`.

On the zero-sum subspace of `V_i`, the clique Laplacian is `n_i I`. A cross
block contributes `d_ij c_ij I` to the diagonal Hessian block. Its off-diagonal
quadratic contribution is

```text
-2 c_ij z_i^T B_ij z_j.
```

Since both `z_i` and `z_j` have zero sum, the all-ones part of `B_ij` vanishes,
leaving `Btilde_ij`. Cauchy--Schwarz and the operator-norm bound give

```text
-2 c_ij z_i^T Btilde_ij z_j
  >= -2 |c_ij| sigma_ij r_i r_j.
```

Summing the diagonal and cross-block contributions yields

```text
z^T H z >= r^T C r.
```

If `C` is positive definite, the right side is positive for every nonzero
`z` in `Z`. The quotient and transverse invariant subspaces are orthogonal and
span the vertex space, so adding strict quotient positivity away from rotation
leaves rotation as the unique zero mode.

## What remains

Theorems 1--3 turn the lifting problem into a realization problem:

1. choose rational class proportions and block densities satisfying the three
   torque equations with strict quotient and limiting transverse margins;
2. realize the fractional blocks by simple biregular graphs; and
3. ensure their centered norms are small enough to make `C` positive definite.

For dense blocks, a centered norm of order `sqrt(N)` is sufficient because the
diagonal entries of `C` are of order `N`. Establishing an arbitrarily-large
existence theorem with all divisibility conditions explicit is the next proof
step.

