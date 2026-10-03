# Clique blow-ups preserve strict Kuramoto stability

This gives an unconditional infinite-family result from any one finite stable
witness. It is weaker than realizing an arbitrary optimized weighted block
graph, but it requires no graphon persistence theorem.

## Construction

Let `G` be a simple graph on `N` vertices. Replace every vertex `v` by a clique
of size `m`. For every edge `{v,w}` of `G`, join the two corresponding cliques
by a complete bipartite graph. This is the lexicographic product `G[K_m]`.

Assign every clone of `v` the original phase `theta_v`.

## Theorem

Suppose `theta` is a nonsynchronous equilibrium of the homogeneous first-order
Kuramoto model on `G`, and its potential Hessian `H_G` is positive semidefinite
with kernel exactly the global rotation vector. Then, for every integer
`m >= 1`, the replicated phase assignment on `G[K_m]` is nonsynchronous and
has a positive-semidefinite Hessian whose kernel is again exactly rotation.

Consequently, every clique blow-up is locally asymptotically stable modulo
global phase rotation.

The finite-dimensional algebraic statement is machine-checked as
`cliqueBlowup_preserves_stable_equilibrium` in
[`KuramotoBlockLimit/CliqueBlowup.lean`](../KuramotoBlockLimit/CliqueBlowup.lean).
Lean verifies equilibrium replication, the actual Kuramoto Hessian matrix,
positive semidefiniteness, and the rotation-only kernel. The final standard
ODE implication to local asymptotic stability modulo rotation remains a
paper-level consequence rather than a formalized dynamical-systems theorem.

## Readable proof corresponding to the Lean theorem

### Equilibrium

Every clone of `v` has `m` neighbors over each original neighbor `w`. Edges
inside its own clique contribute zero torque because their phase difference is
zero. Therefore its torque is

```text
m sum_{w adjacent to v} sin(theta_w-theta_v) = 0.
```

Nonsynchrony is unchanged because the set of phase values is unchanged.

### Hessian splitting

Decompose the vector space on the blow-up into:

1. `U`, vectors constant on every clone clique; and
2. `Z`, vectors whose entries sum to zero inside every clone clique.

The two spaces are orthogonal, invariant, and span the full space.

On `U`, identify a vector with its original-vertex values `z_v`. Internal
clique edges vanish on a constant vector, while every original edge now has
`m` identical neighbors. Thus

```text
H_blowup restricted to U = m H_G.
```

This formula uses unnormalized Hessians. If each vector field includes its
usual `1/|V|` factor, the normalized restrictions are equal instead. The
kernel and all stability signs are unchanged.

It is positive away from the replicated rotation vector.

Now take a vector in `Z`. Write

```text
h_v = sum_{w adjacent to v} cos(theta_v-theta_w) = (H_G)_{vv}.
```

The clique Laplacian acts as `m I` on the zero-sum subspace. For a neighboring
clone clique over `w`, the off-diagonal complete-bipartite term vanishes
because the entries over `w` sum to zero, while its diagonal contribution is
`m cos(theta_v-theta_w)`. Hence the Hessian on the zero-sum directions over
`v` is scalar multiplication by

```text
m(1+h_v).
```

Positive semidefiniteness of `H_G` implies
`h_v = e_v^T H_G e_v >= 0`. Therefore every new transverse eigenvalue is at
least `m>0`. The only remaining zero is the quotient rotation mode.

## Minimum degree

If `G` has minimum degree `delta`, then `G[K_m]` has

```text
N_m     = mN,
delta_m = m(delta+1)-1.
```

Its normalized minimum degree is

```text
mu_m = [m(delta+1)-1]/[mN-1]
```

and converges to `(delta+1)/N`.

For the certified `N=80,002`, `delta=55,018` witness,

```text
16 delta_m - 11(N_m-1)
  = m[16(delta+1)-11N]-5
  = 282m-5 > 0.
```

Thus every `m >= 1` gives another finite stable nonsynchronous graph strictly
above `11/16`, and the infinite family approaches

```text
55019/80002 = 0.6877203069923252...
```

## Limitation

Uniform clique blow-up does not reach the apparent six-block boundary near
`0.691553760598`. Reaching that value requires independently tuning class
proportions and cross-block densities, which is the remaining realization
problem.
