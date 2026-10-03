# Direct proof plan

This track deliberately avoids relying on a general graphon persistence
theorem.

## Proposition A — exact equilibrium lifting

Fix class sizes `n_i` and, for each ordered pair of classes, a biregular
bipartite block with degrees `d_ij`, satisfying edge balance
`n_i d_ij = n_j d_ji`. If phases are constant on classes and

```text
sum_j d_ij sin(phi_j-phi_i) = 0
```

for every class `i`, then every vertex has zero Kuramoto torque. This is exact;
no limiting argument is needed.

## Proposition B — quotient Hessian

The class-constant subspace is invariant under the Hessian. After the change
of variables `q_i=sqrt(n_i) delta-theta_i`, its restriction is symmetric and
is the finite analogue of the mass-orthonormal block Hessian. Prove that its
only zero is rotation whenever the reduced matrix is positive away from the
rotation vector.

## Proposition C — transverse comparison theorem

Decompose the vertex space into the class-constant subspace and the direct sum
of class-zero-sum subspaces. Bound every partial biregular block by the norm of
its centered biadjacency operator. Assemble these bounds into a small
comparison matrix. Positive definiteness of that matrix implies positivity of
the full transverse Hessian.

This proposition should generalize the comparison argument already used by
the finite-witness verifiers.

## Theorem D — arbitrarily large unweighted realizations

For rational block masses and densities with strict quotient and transverse
margins, construct arbitrarily large simple unweighted biregular blow-ups.
Use either:

1. explicit biregular expanders, or
2. existence of random biregular blocks with centered norm `o(N)`.

The diagonal Hessian terms scale as `N`, while centered block fluctuations
scale sublinearly. Therefore the strict limiting transverse margin eventually
dominates the fluctuations. Biregularity preserves the equilibrium equations
exactly.

## Theorem E — within-topology optimum

Use interval Newton/Krawczyk calculations to enclose the candidate boundary
point. Then prove either:

- a constrained second-order sufficient condition for a strict local maximum;
  or
- by interval branch-and-bound, a global upper bound on the compact admissible
  region of this fixed block topology.

The boundary optimizer is expected to have two additional zero Hessian modes.
The stable construction must therefore use nearby interior points, not the
boundary point itself.

