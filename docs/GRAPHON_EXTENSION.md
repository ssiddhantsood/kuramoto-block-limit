# Finite-step graphon extension

Reference: J. J. Bramburger, M. Holzer, and J. Williams,
“Persistence of steady-states for dynamical systems on large networks,”
arXiv:2402.09276.

## Why the published theorem is not plug-and-play here

The published hypotheses use a continuous equilibrium and a row-continuity
condition on the graphon. A nonconstant finite-step graphon has jumps between
partition cells, and the class-constant Kuramoto equilibrium has the same kind
of jumps. The paper's Section 6.4 explicitly sketches replacing `C[0,1]` by a
space of functions continuous on each partition cell and says the idea should
extend to graphons with finitely many jumps, including block models.

There is a separate invertibility issue even at a strict point: Kuramoto
dynamics are invariant under global phase rotation, so the full linearization
always has one zero mode. Any persistence argument must first fix that gauge
or work on the rotation quotient. At the conjectured optimum, two *further*
eigenvalues approach zero after that quotient. The published invertibility
hypothesis therefore still cannot be used at the marginal optimizer without a
degenerate/bifurcation argument.

## Proposed theorem

The general partition-adapted persistence theorem below is still open in this
repository. A narrower [direct random-block theorem](RANDOM_BLOCK_TRANSFER_THEOREM.md)
is now proved for independent edge sampling from any **fixed strict** step
graphon: the sampled phases are corrected to an exact finite equilibrium, and
the Hessian gap and minimum degree persist. It does not cover arbitrary graph
sequences that merely converge in cut norm.

Fix a finite measurable partition of `[0,1]`. Let `W` and the equilibrium `u*`
be continuous on each partition cell (step functions are included). Assume:

1. graph approximants converge to `W` in cut norm;
2. their degree functions converge uniformly on the partition cells;
3. the Kuramoto linearization, restricted to mean-zero functions, is
   invertible and spectrally stable; and
4. the transverse multiplication coefficient is uniformly positive.

Then sufficiently large approximating graph systems possess nearby stable
equilibria modulo rotation.

## Proof tasks

1. Define the partition-adapted Banach space and its mean-zero subspace.
2. Recheck boundedness, compactness, and inverse estimates cell by cell.
3. Repeat the Newton--Kantorovich argument with constants uniform over cells.
4. Repeat spectral convergence after quotienting rotation.
5. Verify the hypotheses for a strict interior point of the six-block family.
6. Compare the conclusion with the more explicit direct blow-up theorem.
