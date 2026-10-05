# Goals and acceptance criteria

## Track A — exact block/blow-up theorem

This is the current priority.

### A1. Finite reduction

Define the reflection-paired six-block family using three class masses, two
fractional block densities, and three phase angles. Derive:

- the three independent zero-torque equations;
- the three normalized class degrees;
- the mass-orthonormal quotient Hessian;
- the even/odd reflection decomposition; and
- the within-class transverse stability quantities.

**Done when:** the formulas are stated without reference to a large adjacency
matrix and reproduce every relevant number from a packaged finite witness.

**Status:** complete at the algebraic/numerical level; interval reproduction is
part of A2.

### A2. Candidate boundary point

Solve the reduced equations near the observed ultrasmall-gap sequence and
identify the limiting point where the smallest even and odd nonrotation modes
vanish simultaneously.

**Done when:** an interval calculation encloses the candidate parameters and
the limiting value of the normalized minimum degree.

**Status:** complete locally in the full nine-density reflected model.  A
Krawczyk calculation encloses the limiting KKT point and its degree; interval
bounds on the seven density-bound multipliers, together with LICQ and reduced
curvature, extend the strict local result off the discovered face.  Global
uniqueness over other density faces is not claimed.

### A3. Direct lifting theorem

Prove that a rational, strictly stable block solution has arbitrarily large
simple unweighted equitable blow-ups with:

- an exact class-constant Kuramoto equilibrium;
- the same quotient stability signs; and
- positive transverse spectrum.

The proof should expose the precise spectral condition required of every
partial biregular block.

**Done when:** the theorem is written in ordinary mathematics, all hypotheses
are checkable, and a verifier constructs/checks finite instances.

**Status:** exact equilibrium lifting, quotient reduction, and the transverse
comparison criterion are proved. Existence of arbitrarily large block
realizations satisfying the spectral bounds remains open.

The exact equilibrium-lifting lemma is also machine-checked in Lean.

The uniform clique-blow-up special case is now completely machine-checked:
equilibrium replication, identification of the actual blown-up Kuramoto
Hessian, positive semidefiniteness, and the rotation-only kernel. This does not
settle the more general partial-biregular realization problem stated above.

### A4. Optimality within the fixed topology

Maximize normalized minimum degree subject to the torque and stability
constraints in the six-block topology.

**Minimum result:** a rigorous local maximum and a certified value.

**Stronger result:** a global upper bound within this fixed topology.

**Explicit non-goal:** this alone is not an upper bound over all graphs or all
step graphons.

**Status:** the minimum result is complete in the full nine-density
three-reflected-pair model: interval calculations certify LICQ, strict
multiplier signs (including all seven active density bounds), positive reduced
curvature, and all inactive constraints.  A complementary formulation reduces
the numerical search to five nonlinear mass/phase variables with a convex SDP
over all densities.  A stability-margin search comprising 6,060 outer
evaluations and 5,951 SDP calls found no strict counterexample, while an
exact-decimal dual plus interval arithmetic now gives
rigorous upper bounds on individual outer boxes.  A first exact-box
branch-and-bound run is reproducible, but it leaves most of the domain volume
unresolved and is not a global certificate.  The exact degree-`3/4` family is
now known to lie on a six-dimensional marginal sheet; every second-flat
critical direction is tangent to that sheet.  The global upper bound remains
open because the outer domain and all-orders approaches transverse to the
marginal sheet have not been exhaustively excluded.

## Track B — finite-step graphon persistence

### B1. Isolate the missing hypotheses

Document exactly why the published theorem does not apply verbatim: a finite
step graphon and its class-constant equilibrium have jump discontinuities, and
the limiting optimizer also loses invertibility at its stability boundary.

### B2. Partition-adapted persistence theorem

Replace the continuous function space with functions continuous on each cell
of a fixed finite partition. Prove persistence and stability after quotienting
the global phase-rotation mode.

**Done when:** the proof covers a strictly stable finite-step Kuramoto
equilibrium and graph sequences converging in cut norm with uniformly
convergent degree functions.

### B3. Application and comparison

Apply the theorem to interior points of the six-block family and compare its
conclusion with Track A's direct equitable-blow-up theorem.

**Boundary qualification:** neither persistence theorem should claim to apply
at the marginal optimizer without an additional bifurcation/degeneracy
argument.

## Order of work

1. A1 finite reduction.
2. A2 boundary candidate and interval formulation.
3. A3 direct lifting theorem.
4. A4 within-family optimality.
5. B1--B3 graphon extension.
