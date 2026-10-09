# Exact classification of the degree-`3/4` phase slice

Consider three reflected pairs with positive pair masses `a,b,c`, where
`a+b+c=1/2`, and phases

```text
(theta_A,theta_B,theta_C) = (-pi/4,-3*pi/4,pi/4).
```

Each class has within-class density `1`.  The other densities lie in `[0,1]`:
`v_i` joins a class to its reflection; `s_ij` joins equal signs in distinct
pairs; and `t_ij` joins opposite signs.  Let `mu` be the minimum of the three
class degrees.  Weak stability means that the Kuramoto Hessian is
positive semidefinite off global rotation, including its transverse sector.

**Theorem, at these exact phases.** Every weakly stable equilibrium satisfies
`mu <= 3/4`.  Equality holds exactly when `b=1/4`, `a+c=1/4`, and

```text
(v_A,v_B,v_C,s_AB,t_AB,s_AC,t_AC,s_BC,t_BC)
  = (1,1,1,1,0,1,1,0,1).
```

Thus `0<a<1/4` parameterizes the equality family.  Every weakly stable
equilibrium at these phases has extra zero quotient modes, so none is
strictly stable.

## Proof

In the three-by-three even Hessian, the `B` diagonal is exactly

```text
E_BB = -a*t_AB - c*s_BC.
```

Positive semidefiniteness requires `E_BB >= 0`.  Since `a,c>0` and both
densities are nonnegative, `t_AB=s_BC=0`.  The three zero-torque equations
then uniquely determine the reflected-pair densities:

```text
v_A = (b*s_AB - c*s_AC)/a,
v_B = (a*s_AB + c*t_BC)/b,
v_C = (b*t_BC - a*s_AC)/c.
```

On this torque-zero sheet, the even and odd quotient Hessians are

```text
E = t_AC * [[ c, 0, -sqrt(a*c)],
            [ 0, 0,          0],
            [-sqrt(a*c), 0, a]],

O = t_AC * [[ c, 0, +sqrt(a*c)],
            [ 0, 0,          0],
            [+sqrt(a*c), 0, a]].
```

Each has rank at most one.  In particular, `E` has a zero mode besides
rotation and `O` has two zero modes.  This proves the failure of strict
quotient stability at the exact phases, regardless of the transverse sector.

The degrees simplify to

```text
d_A = a + 2*b*s_AB + c*t_AC,
d_B = b + 2*a*s_AB + 2*c*t_BC,
d_C = c + a*t_AC + 2*b*t_BC.
```

The following identities make the degree bound transparent:

```text
3*b - d_B = 2*b*(1-v_B),

3*(a+c) - d_A
  = 2*a*(1-v_A) + 2*c*(1-s_AC) + c*(1-t_AC),

3*(a+c) - d_C
  = 2*c*(1-v_C) + 2*a*(1-s_AC) + a*(1-t_AC).
```

Every term on the right is nonnegative.  Hence

```text
mu <= min(3*b, 3*(a+c))
   <= [3*b + 3*(a+c)]/2
    = 3/4.
```

For equality, both upper bounds must equal `3/4`, so `b=a+c=1/4`.
Equality in `d_B<=3b` gives `v_B=1`.  The torque formula then gives
`a*(1-s_AB)+c*(1-t_BC)=0`; positivity of `a,c` and the density box force
`s_AB=t_BC=1`.  The formula for `v_A` and its upper bound force `s_AC=1`
and `v_A=1`.  Finally `d_A>=3/4` forces `t_AC=1`, and the remaining torque
formula gives `v_C=1`.  This is the displayed binary family.  Direct
substitution gives all three degrees `3/4`, all three transverse values
`1/4`, and nonnegative rank-one even and odd sectors, so equality is attained
in the weak problem.

The identities are independently reconstructed from the six-by-six block
matrix and checked exactly by SymPy.  With `requirements-optimization.txt`
installed, run:

```bash
python scripts/certify_fixed_angle_slice.py
```

This theorem applies **only to the exact phase triple above**.  It does not
bound nearby-angle strictly stable equilibria, prove the conjectured
`0.691553760598...` global three-pair bound, or transfer a step graphon to
finite unweighted graphs.
