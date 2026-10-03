# Reflection-paired six-block model

## Variables

Order the six classes as

```text
A+, A-, B+, B-, C+, C-.
```

Each class in a reflection pair has mass `a`, `b`, or `c`, respectively, so

```text
a > 0, b > 0, c > 0,  a + b + c = 1/2.
```

Assign phases

```text
alpha, -alpha, beta, -beta, gamma, -gamma.
```

The observed limiting block-density matrix has two free densities `x,y`:

```text
      A+ A- B+ B- C+ C-
A+     1  x  y  0  1  1
A-     x  1  0  y  1  1
B+     y  0  1  1  0  1
B-     0  y  1  1  1  0
C+     1  1  0  1  1  0
C-     1  1  1  0  0  1
```

Here a diagonal `1` means a clique inside a finite blow-up class; it does not
contribute to class-constant torque or quotient motion, but it does contribute
to transverse stability.

## Degree objective

The three distinct normalized degrees are

```text
dA = a(1+x) + b y + 2c,
dB = 2b + a y + c,
dC = c + 2a + b = 1/2 + a.
```

Thus the reduced objective is

```text
maximize mu = min(dA,dB,dC).
```

The identity `dC = 1/2 + a` is an exact structural simplification: within this
topology, improving the bottleneck eventually requires increasing `a`.

## Equilibrium equations

For the homogeneous first-order Kuramoto flow, the three independent class
torques are

```text
TA = a x sin(-2 alpha)
     + b y sin(beta-alpha)
     + c[sin(gamma-alpha) + sin(-gamma-alpha)],

TB = a y sin(alpha-beta)
     + b sin(-2 beta)
     + c sin(-gamma-beta),

TC = a[sin(alpha-gamma) + sin(-alpha-gamma)]
     + b sin(-beta-gamma).
```

The partner torques are their negatives. Equilibrium is therefore equivalent
to `TA=TB=TC=0`.

## Quotient Hessian

Let `p_i`, `phi_i`, and `W_ij` denote the six masses, phases, and block
densities. In mass-orthonormal coordinates the class-constant Hessian is the
symmetric matrix

```text
H_ii = sum_{j != i} p_j W_ij cos(phi_i-phi_j),
H_ij = -sqrt(p_i p_j) W_ij cos(phi_i-phi_j),  i != j.
```

The vector `(sqrt(p_1),...,sqrt(p_6))` is the rotation zero mode. Strict
quotient stability means that all five remaining eigenvalues are positive.
Reflection symmetry decomposes these five modes into a two-dimensional even
sector (after removing rotation) and a three-dimensional odd sector.

## Transverse quantities

For a zero-mean perturbation supported inside class `i`, the step-graphon
Hessian has multiplication coefficient

```text
q_i = sum_j p_j W_ij cos(phi_i-phi_j),
```

including the diagonal block. Positivity of every `q_i` controls the essential
transverse spectrum of the ideal step graphon. A finite deterministic blow-up
also requires control of the centered adjacency operators of its partial
biregular blocks; positivity of the quotient alone is not sufficient.

## Numerical boundary conjecture

The current sequence suggests

```text
mu_* approximately 0.691553760598,
```

with the smallest nonrotation eigenvalue in both the even and odd sectors
approaching zero. The natural boundary system therefore consists of the three
torque equations, active degree constraints, and one zero-eigenvalue condition
from each reflection sector.

