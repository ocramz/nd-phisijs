# N-Dimensional Rigid Body Dynamics — Technical Summary for a 4D Game Engine

Source: Marc ten Bosch, "N-Dimensional Rigid Body Dynamics", ACM Transactions on
Graphics 39(4), Article 1, July 2020.

---

## 0. How to read this document

This document is a build guide. It has three parts:

- **Part A — Mathematics.** These rules are exact. They are true for all
  dimensions `n > 1`. They have no tolerances and no time steps.
- **Part B — Numerics.** These rules are approximate. They have time steps,
  tolerances, iteration counts and stability limits.
- **Part C — 4D data.** These are the constant tables for `n = 4`.

Keep Part A and Part B in different modules. Part A must not include a
tolerance value. Part B must not include a new algebraic rule.

---

## 1. Scope and result

The paper shows one result: the equations of rigid body motion do not change
when the number of dimensions changes. Only the data types change.

To get this result, do not use the cross product and do not use axial vectors.
The cross product is correct only in 3D. Use these types in their place:

| Quantity | 3D usual type | General type |
|---|---|---|
| Position, velocity, force, impulse | vector | vector, `n` components |
| Orientation | quaternion | rotor, `2^(n-1)` components |
| Angular velocity, torque, angular momentum | axial vector | bivector, `k` components |
| Inertia | 3x3 matrix | `k` x `k` matrix |

Where:

```
k = C(n, 2) = n (n - 1) / 2      number of bivector components
r = 2^(n-1)                      number of rotor components
```

For `n = 4`: `k = 6` and `r = 8`.
For `n = 3`: `k = 3` and `r = 4`. The library thus gives the usual 3D engine
as a special case. Use this property to test the library.

---

## 2. Recommended module layout

Write the library with `n` as a compile-time constant. Do not make `n` a
run-time variable. Many table sizes come from `n`, and the compiler can
unroll the loops.

```
core/          n, k, r, the bivector basis order, index tables
algebra/       PART A. vector, bivector, rotor, multivector, products
body/          PART A. mass properties, inertia tensor, state
integrate/     PART B. time step, gyroscopic term, rotor correction
detect/        PART A + B. Minkowski difference, separating axis theorem
resolve/       PART B. impulses, friction, contact graph, shock propagation,
               constraints and joints
slice/         4D only. slice of the mesh, display, user input
```

The layers `algebra`, `body`, `detect` and `resolve` are general for all `n`.
Only `slice` is specific to 4D.

---

# PART A — MATHEMATICS

## A1. Basis and storage

A bivector has one component for each pair of different basis vectors. Put the
pairs in lexicographic order. This order is a free choice, but you must use the
same order everywhere in the library.

For `n = 4` the order is:

```
index 0: e_xy
index 1: e_xz
index 2: e_xw
index 3: e_yz
index 4: e_yw
index 5: e_zw
```

Store a bivector as a flat array of `k` numbers. Store a vector as a flat array
of `n` numbers. The memory layout is the same.

A rotor is a member of the even sub-algebra. It contains only grades 0, 2, 4,
and so on. For `n = 4` a rotor has:

- 1 scalar component (grade 0)
- 6 bivector components (grade 2)
- 1 quadvector component (grade 4)

This gives 8 numbers. Store a 4D rotor as a flat array of 8 numbers.

**Caution:** in 3D a rotor has no grade-4 part. Thus a 3D rotor keeps one plane
of rotation. In 4D a body can turn in two independent planes at the same time.
The grade-4 part is not zero in that condition. Do not copy a quaternion class
and add components. Write the even sub-algebra as its own type.

## A2. Products

The library needs these products:

| Product | Symbol | Use |
|---|---|---|
| Geometric product | `a b` | rotor composition, rotor build |
| Exterior product | `a ∧ b` | bivectors, blades, separating axes |
| Left contraction (dot) | `a · B` | point velocity, projection |
| Commutator | `A × B = (A B − B A) / 2` | gyroscopic term |
| Reverse | `R~` | rotor inverse for unit rotors |

For two vectors: `R = a b = a · b + a ∧ b`. This rotor turns a vector in the
plane of `a` and `b`, through twice the angle between `a` and `b`.

To turn a vector with a rotor, use the sandwich product:

```
x' = R x R~
```

The same sandwich product turns a bivector, and it turns any blade.

## A3. The star matrix

Define a `k` x `n` matrix `[r]*` for a vector `r`:

```
r ∧ a  = [r]* a          (a is a vector, the result is a bivector)
r · w  = [r]*^T w        (w is a bivector, the result is a vector)
```

`[r]*` is the general form of the 3D cross-product matrix `[r]x`, but it does
not take the dual. This matrix is the key tool of the paper. It connects the
vector world (positions, forces, impulses) to the bivector world (torques,
angular velocity, angular momentum). Almost all of the mechanics code uses it.

The instantaneous velocity of a point `r` on the body, in the body frame, is:

```
u = v + r · w = v + [r]*^T w
```

## A4. The rotor matrix for bivectors

Define a `k` x `k` matrix `[R]2` such that, for a bivector `B` seen as a vector
of length `k`:

```
[R]2 B = R B R~
```

Build `[R]2` one time per body per frame. Then you can transport the inertia
tensor with a matrix product, and not with an algebra product.

## A5. State and equations of motion

The state of a body is:

```
x   position          vector      n
R   orientation       rotor       2^(n-1)
v   linear velocity   vector      n
L   angular momentum  bivector    k
```

The equations of motion are:

```
dx/dt = v                 dv/dt = F / m
dR/dt = −(1/2) w R        dL/dt = tau
```

`F` is the net force (vector). `tau` is the net torque (bivector). `m` is the
mass (scalar). The product `w R` is the geometric product.

**Note:** the rotor equation has the same form as the quaternion equation of a
3D engine. You can put it in the place of the quaternion equation with no other
change to the integrator structure.

Angular momentum of a particle with momentum `p` at position `x` is the
bivector `L = x ∧ p`.

## A6. The inertia tensor

Angular momentum is a linear map from bivectors to bivectors:

```
L(w) = INT_V  r ∧ (r · w) dV
```

Do not take the dual of this map. The dual is correct only in 3D. Keep the map
as a `k` x `k` matrix:

```
I = INT_V  rho [r]* [r]*^T dV
L = I w
```

The integral is over the `n`-dimensional volume of the body. For `n = 4` this
is a 4-volume, not a 3-volume.

Call `dI = [r]* [r]*^T` the inertia density matrix. Part C gives `dI` for
`n = 2, 3, 4`.

To move the inertia tensor from the body frame to the world frame:

```
I' = [R]2 I [R]2^T
```

## A7. The Euler equation

Take the time derivative of `L = R I(R~ w R) R~`. Only the rotors change with
time. The result is the Euler equation:

```
dL/dt = I(dw/dt) − w × I(w) = tau
```

Here `×` is the commutator product of two bivectors, and not a cross product.
In the world frame:

```
I' dw/dt − w × I' w = tau
```

**Implementation note:** the sign of the commutator term depends on two
conventions: the sign in `dR/dt = −(1/2) w R`, and the order in the commutator.
Fix both conventions one time. Then test the sign with a torque-free body: the
kinetic energy and the angular momentum must stay constant.

## A8. Mass properties of a simplicial mesh

An `n`-dimensional simplicial mesh has a surface of `(n-1)`-simplices:

- 3D mesh: the surface is a set of triangles.
- 4D mesh: the surface is a set of tetrahedra.

The inertia tensor is linear. Thus you can sum the contribution of each
`n`-simplex. Build each `n`-simplex from one surface `(n-1)`-simplex and the
origin of the frame.

Do not integrate the products of inertia directly. Use the covariance method
instead. It is simpler and it extends to all `n`.

**Step 1. Canonical covariance.** The canonical `n`-simplex has these two
values:

```
C_xx = 2 / (n + 2)!
C_xy = 1 / (n + 2)!
```

Put `C_xx` on the diagonal of the `n` x `n` matrix `C`, and `C_xy` in all other
positions.

**Step 2. Transform.** Let `M` be the `n` x `n` matrix that changes the
canonical simplex vertices into the vertices of the actual simplex. Then:

```
C' = det(M) M C M^T
volume = det(M) / n!
mass   = rho * volume
```

**Step 3. Covariance to inertia.** Expand the star matrix in the coordinate
basis:

```
[r]* = SUM_i  r_i [e_i]*
```

Then:

```
I = SUM_{i,j}  C'_ij  [e_i]* [e_j]*^T
```

The matrices `[e_i]*` are constants. Build them one time at start-up.

**Step 4. Diagonal frame.** It is useful to select a local frame in which `I` is
diagonal. Do not diagonalize `I`. A `k` x `k` rotation matrix applies to
bivectors, not to vectors. You cannot use it to rotate the mesh.

Do this in its place:

1. Diagonalize the `n` x `n` covariance matrix `C'`.
2. Rotate the body vertices with that `n` x `n` rotation.
3. Calculate `I` again from the rotated mesh.

If the covariance matrix is diagonal, then the inertia tensor is also diagonal.

## A9. Collision detection — hypersphere against polytope or cylinder

Extend the usual sphere test that uses the Minkowski difference:

1. Find the point on the surface of the object that is nearest to the center of
   the hypersphere.
2. Compare the distance with the radius of the hypersphere.

The nearest point can be on an `(n-1)`-cell, or on a boundary that two or more
cells share. The dimension of the boundary depends on the number of cells that
touch it. In 4D, two 3-cells share a face, three 3-cells share an edge, and so
on. Your code must walk down through all these boundary dimensions.

To project a vector `a` onto a subspace with blade `B`, use:

```
a_parallel = (a · B) / B
```

## A10. Collision detection — two convex polytopes

Extend the separating axis theorem. An `n`-polytope contains `m`-cells for
`m < n`. Each `m`-cell `i` is spanned by `m` vectors. The exterior product of
these vectors is the `m`-blade `V_a^m(i)`.

**Rule.** For each pair of cell dimensions `(m_a, m_b)` such that
`m_a + m_b = n − 1`, form all the exterior products:

```
V_a^{m_a}(i)  ∧  V_b^{m_b}(j)     for all i, j
```

The candidate separating axes are the dual vectors of these `(n-1)`-blades.

The case `m_a = n − 1`, `m_b = 0` is the empty exterior product. It gives the
normals of the `(n-1)`-cells of polytope `a`.

Examples:

- **3D.** The pairs are (2,0), (0,2) and (1,1). This gives face normals and
  edge-against-edge axes. This is the usual 3D test.
- **4D.** The pairs are (3,0), (0,3), (2,1) and (1,2). This gives cell normals,
  and the products of an edge direction of one body with a face bivector of the
  other body. The edge-against-edge case is inside the (2,1) and (1,2) cases.
  You do not need a separate test for it.

## A11. Collision response

Apply equal and opposite impulses `j` to the two bodies:

```
v' = v ± j / m
w' = w ± I^-1 (r ∧ j) = w ± I^-1 [r]* j
```

`r` points from the center of mass of each body to the contact point.

The velocity at the contact point changes as `u' = u ± K j`, with:

```
K = delta / m + [r]*^T I^-1 [r]*
```

`delta` is the `n` x `n` identity matrix. `K` is an `n` x `n` matrix.

`K` has the same shape as the 3D case. Thus you can use your 3D formulas for
the normal impulse, for static friction and for kinetic friction with no
change.

## A12. Forces

Gravity and air friction extend to `n` dimensions with no change of form.
Gravity points along the normal of the ground, towards the ground.

## A13. Constraints and joints

The paper does not give these. This section gives the form that follows from
A3, A5 and A11.

A constraint is a list of **rows**, and each row takes away one degree of
freedom. A row is linear or angular. The error and the velocity always measure
`b` against `a`.

**A linear row.** It holds a unit direction `d` of `n` components at the world
offsets `rA` and `rB`. Its velocity and its effective mass are

```
C'   = (v_b + [rB]*^T w_b - v_a - [rA]*^T w_a) . d
K    = 1/m_a + (rA ^ d) . I'a^-1 (rA ^ d)
     + 1/m_b + (rB ^ d) . I'b^-1 (rB ^ d)
```

`K` is the value that A11 gives for a contact, along `d`. Thus a linear row is
a contact row with no limit of sign. The impulse `j = lambda d` goes to `b` at
`rB` and `-j` goes to `a` at `rA`, thus it keeps the momentum of the pair.

**An angular row.** It holds a unit bivector `x` of `k` components:

```
C'   = (w_b - w_a) . x
K    = x . I'a^-1 x + x . I'b^-1 x
```

The impulse is a torque impulse `lambda x`, and `-lambda x` goes to `a`.

**The error of the rotation.** A row needs the error of the position of the
rotation, and not only its velocity. The rotor holds that, but a logarithm of a
rotor is expensive. Take the antisymmetric part of the relative rotation matrix
instead. With `Q` the relation of the two bodies at rest,

```
T    = Rm_a Q                             the rotation that b must have
M    = Rm_b T^T                           the error rotation, in the world frame
e[p] = (M[j n + i] - M[i n + j]) / 2      for pairs[p] = [i, j]
```

`e` is the logarithm of `M` to the first order, and it is a bivector of `k`
components. It carries the same sign as `w`: a value of more than zero is a
turn from the axis `i` toward the axis `j`.

**When `e` is exact.** For a **simple** rotation, thus a rotation in one plane,
`e` holds exactly `sin(angle)` on that plane and 0 on each other plane, at any
angle. For a rotation in two planes together `e` mixes: a turn of 1.0 in the
plane `(x y)` and of 1.0 in the plane `(y z)` gives 0.354 on the plane `(x z)`,
in which nothing turned.

Thus a joint may only use `e` when the part that stays free is one plane or
nothing, thus when the joint holds `k` or `k - 1` planes. With two or more
planes free, the joint holds the angular velocity only.

**The five types, and the degrees of freedom that each one takes away.**

| type | linear rows | angular rows | n = 3 | n = 4 |
|---|---|---|---|---|
| point | `n` | 0 | 3 of 6 | 4 of 10 |
| distance | 1 | 0 | 1 of 6 | 1 of 10 |
| fixed | `n` | `k` | 6 of 6 | 10 of 10 |
| hinge | `n` | `k - 1` | 5 of 6 | 9 of 10 |
| subspace | the count that you name | the count that you name | | |

A hinge leaves one rotation **plane** free, and not one axis. In 3D that is 1
plane of 3, which is the door hinge that you know. In 4D it is 1 of 6, thus a
hinge takes away 5 angular degrees of freedom. In 2D, `k` is 1, thus a hinge
has no angular row and it is the same as a point joint.

The angular rows of a hinge are the `k - 1` unit bivectors orthogonal to the
free plane. There is no axis, and there is no cross product.

The ANGLE of a hinge does not come from the error of this section. It comes
from `atan2` in the free plane, and that is exact at any angle. See B10.

---

# PART B — NUMERICS

## B1. Time integration

The paper uses a symplectic Euler scheme, in the style of Guendelman, Bridson
and Fedkiw (2003). Other schemes are also possible. The algebra of Part A does
not change with the scheme.

Integrate the gyroscopic term of the Euler equation separately, with an
implicit Euler method (Catto, 2015). The gyroscopic term is the commutator
term `w × I' w`. It is stiff. An explicit method makes it unstable.

**Note:** the implicit method removes a small quantity of angular momentum over
long times. This is a known effect. Do not read it as a bug in the algebra.

## B2. Rotor error correction — critical

Small errors increase in the rotor `R` at each time step. In 3D you correct a
quaternion when you normalize it. **This is not sufficient in 4D and higher.**

Simple rotors, which are the geometric product of two unit vectors, are on a
sphere in the even sub-algebra. Rotors of a double rotation are not. Thus a
norm is not enough to make the rotor a correct rotation.

Do this in its place:

1. Factorize the rotor into a set of vectors whose geometric product is the
   rotor. Use the algorithm of Perwass (2009).
2. The algorithm gives normalized vectors.
3. Multiply the vectors again to get the corrected rotor.

The result is always a correct rotation.

**Caution:** the gyroscopic term makes double rotations very frequently. Thus
you cannot skip this step in a 4D engine. Put this correction in the numerics
module, and not in the algebra module.

## B3. Matrix work per step

Per body, per frame:

- Build the `k` x `k` matrix `[R]2` from the rotor.
- Calculate `I' = [R]2 I [R]2^T`.
- Calculate `I'^-1`. Or, invert `I` one time in the body frame and transform
  the inverse: `I'^-1 = [R]2 I^-1 [R]2^T`. This is less costly.

For `n = 4`, `I` is 6 x 6. This is small. A direct inverse is sufficient. Use a
Cholesky factorization, because `I` is symmetric and positive definite.

## B4. Cost

```
n     k = n(n-1)/2     rotor size 2^(n-1)     inertia matrix
2     1                2                      1 x 1
3     3                4                      3 x 3
4     6                8                      6 x 6
5     10               16                     10 x 10
6     15               32                     15 x 15
```

The rotor size increases exponentially. The inertia matrix increases as `n^4`
in element count. For a game engine at `n = 4` or `n = 5` the cost is small.
Do not plan for large `n`.

## B5. Separating axis theorem cost in 4D

A 4D hypercube has 8 cells, 24 faces, 32 edges and 16 vertices. But it is
symmetric. Thus it has only 4 unique directions (normals or edge directions),
and you can form only 6 unique bivectors.

Axis count for two oriented hyperboxes:

```
cell normals:        4 x 2                 =  8
exterior products:   4 x 6 x 2             = 48
total                                      = 56
```

Compare this with 15 axes in 3D. The 4D test costs about 4 times more. Extend
the method of Eberly (2002) for oriented boxes. Replace each scalar triple
product with an exterior product.

## B6. Contact and stacking

Use the Contact Graph and the Shock Propagation methods of Guendelman, Bridson
and Fedkiw (2003). They do not change with the dimension. They keep stacks of
bodies stable.

## B7. Tolerances to expose

Put all of these values in the numerics module:

- time step, and the number of sub-steps
- restitution and friction coefficients
- contact distance and penetration slop
- rotor factorization threshold (near-identity rotors)
- velocity thresholds for sleep and for static friction
- iteration count for the contact solver
- constraint bias, constraint slop, the largest bias of a constraint, and the
  shortest direction that a constraint accepts
- the frequency and the damping ratio of a soft joint, and the largest part of
  the rate of a substep that the frequency may use

## B8. The constraint solver

The rows of A13 go through the same sequential impulse method as the contacts.

**The bias.** A row drives its velocity toward `-(bias / dt) * error`, and not
toward zero. That gives back a part of the error in each step. Use a bias of
the constraints that is separate from the bias of the contacts: 0.2 is tuned
for the push-out of a contact, and a joint needs its own value.

**The slop.** Make the error smaller by the slop along the direction of the
whole error vector, and not component by component. A slop on each component
makes a dead zone in the shape of a box, and the joint is then not the same in
every direction.

**The limit of the bias.** The bias is an error divided by the time of a step.
A body that a teleport moves 50 units would ask for 6000 units each second.
Hold the bias at a largest value.

**The clamps.** Hold the TOTAL impulse of a row between `lower` and `upper`,
and apply only the change. An equality row has no limit. A rope is a row with
`upper` at 0, and a strut is a row with `lower` at 0.

**The warm start.** A joint keeps its impulses on its own rows, because a joint
lives between the steps. A contact needs a cache, because the narrow phase
builds a new contact in each step.

**The order in the loop.** Solve the joints before the contacts on an even
turn, and after them on an odd turn. The odd turn is then a true reversal of
the even turn, thus the change of direction still takes away the effect of the
order. A contact is an inequality, thus its clamp runs last on the even turn.

**The complement basis must be continuous.** The angular rows of a hinge are
the complement of the free plane. A basis that comes from the axes at each step
turns over when the free plane passes a tie in the order of the axes. The
impulse that the row keeps then goes along a different bivector, and that adds
energy: a jolt one time in each turn. Start from the basis of the last step,
and only build a new bivector when an old one is not usable.

**Shock propagation.** That pass makes the impulses of the contacts zero and
solves them again with false masses, thus it can pull a joint open. Run one
more pass on the joints after it. Do not put the joints in the contact graph:
that would change the stacking of every scene.

**The slop and the free directions.** Take the error into the frame of the rows
first, and only then make it smaller by the slop. A direction that the joint
leaves free carries an error that grows without limit. That error would make
the length of the error vector large, and the slop would then do nothing at
all on the directions that the joint does hold.

**What this method does not give.** The solver works at the velocity level, and
it holds the lever arms still through a step. The energy therefore falls by an
amount of the first order in `dt`. It never rises.

**Iterations and sub-steps do different work.** `iterations` brings the
VELOCITY of a row to its target. `subSteps` builds the rows again, thus it
makes the error of the POSITION smaller. Which one helps depends on the
condition of the scene. A chain of 10 links, the largest error over 10 seconds:

| the start | it = 10 | it = 100 | it = 1000 | 4 sub-steps |
|---|---|---|---|---|
| horizontal, and it whips | 0.1142 | 0.1146 | 0.1149 | 0.0137 |
| hanging, and it settles | 0.0014 | 0.000034 | 0.000000 | 0.000089 |

Thus a chain that has come to rest goes to its answer with more iterations,
and a chain that moves fast does not: there the error comes from the lever
arms that stay still through the step, and only a shorter step helps. A large
ratio of the masses also needs sub-steps.

**The drift of a lock that leaves two or more planes free.** Such a lock has no
geometric error that means anything, because `SO(n)` is not abelian: a turn of
1.0 in the plane `(x y)` and of 1.0 in the plane `(y z)` reads 0.919 in the
plane `(x z)`, and the reading changes with the order of the two turns. No
logarithm can take that away, exact or not.

Integrate the violated speed of the row in the place of it:

```
drift += ((w_b - w_a) . axis) dt
```

Add it up AFTER the solver and BEFORE the step of the position, thus `w` and
`axis` are the values that the integrator uses. The axis turns during the step,
thus the value is correct to the first order only. A warm start that is off
must not make the drift zero, and a teleport must.

## B9. Soft constraints

A raw bias of Baumgarte gives a stiffness that changes with the length of the
step, and it can put energy in. A soft constraint is a spring of a frequency
that you name, and its stiffness does not change with the step. The method is
from Erin Catto, "Soft Constraints", GDC 2011.

For each joint, one time in each substep:

```
om = 2 pi hertz
a1 = 2 zeta + dt om;  a2 = dt om a1;  a3 = 1 / (1 + a2)
biasRate = om / a1;  massScale = a2 a3;  impulseScale = a3
```

The row then solves

```
dJ = -massScale (C' - bias) / K - impulseScale (the total impulse)
```

A `hertz` of 0 must give back the rigid solver exactly: `massScale` at 1,
`impulseScale` at 0, and the old bias rate. Write the line so that a mass scale
of 1 and an impulse scale of 0 give the same numbers, and not only near ones.

The stretch of a soft joint under a weight follows the formula of a spring:
`g / (2 pi f)^2`. A mass of 1 kg on a joint of 5 Hz hangs 9.94 mm low. Use this
to test the three factors, and not only that the joint is stable.

A frequency above about `0.25 / dt` is not stable. Hold it there.

The row of a limit and the row of a motor stay RIGID. A limit that gives way is
not a limit, and a motor is a source of speed and not a spring.

## B10. The limit and the motor of a hinge

**The angle.** With `Q` the relation of the two bodies at rest, and `(u, v)` an
orthonormal basis of the free plane turned into the world by `Rm_a`:

```
M = Rm_b Q^T Rm_a^T
angle = atan2(v . (M u), u . (M u))
```

This is EXACT at any angle, and it stays exact when the two bodies also turn in
other planes. The first-order logarithm of A13 is not exact there; this is,
because the free part of a hinge is ONE plane, and a rotation in one plane
commutes with itself. The value wraps at `pi`.

Take `(u, v)` one time, when the joint starts. Build the `n` x `n` antisymmetric
matrix of the plane bivector, take the column of the largest length as `u`, and
take `v` as that matrix times `u`. A basis rebuilt at each step would jump when
the largest column changed.

The plane must be SIMPLE, thus `h ^ h = 0`. In 4 dimensions and more a bivector
can hold two planes at the same time, for example `e_xy + e_zw`. That is not a
hinge. Refuse it.

**The limit.** One row on the free plane, on one side only:

| condition | error | lowest | highest |
|---|---|---|---|
| angle below the lower limit | `angle - lower` | 0 | +inf |
| angle above the upper limit | `angle - upper` | -inf | 0 |
| between the two | the row has no mass, thus it goes away | | |

**The motor.** One row on the free plane, with the target speed as its bias and
with the clamp `maxMotorTorque * dt`. Use a TORQUE and not an impulse: an
impulse would make the motor stronger when the count of the sub-steps grew.

Put the motor BEFORE the limit. The solver walks the rows forward, thus the
limit has the last word and a motor that drives into a limit stalls.

**The two rows are always there**, and they have no mass when they are off. The
rows of a joint are named by their position, and a count that changes throws
away every impulse that the joint keeps -- at the moment that a limit engages,
which is the worst moment.

## B11. Islands and sleep

A body that is almost still for `sleepTime` seconds may sleep. Two rules make
that work:

**The solver must not reset the timer.** An impulse of a contact or of a row is
internal, and it holds the body where it is. A wake that resets the timer at
every turn of the solver means that a box on the ground can never sleep,
because its contact gives it an impulse in every turn. Give the solver a wake
that only acts on a body that IS asleep.

**The bodies sleep in islands.** A contact and a joint each join two bodies. A
whole island sleeps together, or none of it sleeps. Without that rule one box
of a stack could sleep while the box under it still moves. A joint must be an
edge of the island too, thus two bodies that a joint holds sleep at the same
time even when no contact joins them.

A static body is in no island. The ground touches everything, and it would make
one island of the whole world.

---

# PART C — 4D CONSTANT TABLES

## C1. Star matrix, 4D

Basis order: `x, y, z, w` for vectors. `xy, xz, xw, yz, yw, zw` for bivectors.

```
        |  -y   x   0   0 |
        |  -z   0   x   0 |
[r]* =  |  -w   0   0   x |     (6 rows x 4 columns)
        |   0  -z   y   0 |
        |   0  -w   0   y |
        |   0   0  -w   z |
```

## C2. Inertia density matrix, 4D

```
       | x²+y²   yz     yw    -xz    -xw     0   |
       |  yz    x²+z²   zw     xy      0    -xw  |
dI =   |  yw     zw    x²+w²    0      xy    xz  |
       | -xz     xy      0    y²+z²    zw   -yw  |
       | -xw      0     xy      zw   y²+w²   yz  |
       |   0    -xw     xz    -yw     yz   z²+w² |
```

`dI` is symmetric. `dI = [r]* [r]*^T`. Use this identity as a unit test.

## C3. Star matrix and inertia density, 2D and 3D

2D (`k = 1`):

```
[r]* = ( -y   x )                dI = ( x²+y² )
```

3D (`k = 3`, basis `xy, xz, yz`):

```
        | -y   x   0 |          | x²+y²   yz    -xz  |
[r]* =  | -z   0   x |    dI =  |  yz    x²+z²   xy  |
        |  0  -z   y |          | -xz     xy   y²+z² |
```

**Caution:** the 3D `dI` above is not the classical 3D inertia tensor. The
classical tensor is the dual form. The two give the same physics, but the
components are different. Do not mix the two forms.

## C4. Canonical simplex covariance values

```
n = 2:  C_xx = 2/24 = 1/12       C_xy = 1/24
n = 3:  C_xx = 2/120 = 1/60      C_xy = 1/120
n = 4:  C_xx = 2/720 = 1/360     C_xy = 1/720
```

## C5. The 4D Euler equation in components

For torque-free motion, the component for the plane `(i j)` is:

```
I_ij dw_ij/dt = (I_jk − I_ik) w_ik w_jk + (I_jm − I_im) w_im w_jm
```

`k` and `m` are the two other indices.

Read these properties from the equation:

- The equation has a sum of two terms. In 3D it has only one term. Thus the
  stability of a rotation in 4D does not follow only from distinct principal
  moments of inertia.
- The plane `(i j)` depends on the planes `(i k)`, `(j k)`, `(i m)` and
  `(j m)`. It does not depend on the perpendicular plane `(k m)`. This agrees
  with the property that a 4D body turns in two perpendicular planes
  independently.

---

# PART D — 4D DISPLAY AND USER INPUT

This part is specific to the 4D game. It is not part of the general library.

## D1. Slice

A 4D object has a 3-dimensional surface. The mesh of that surface is a set of
tetrahedra.

To display the object, cut the 3D surface with one 3D hyperplane. This gives a
2D surface. The slice of one tetrahedron is a triangle or a quadrilateral. All
these polygons together make the surface to display.

The user moves a slider to move the 3D slice along the fourth axis.

## D2. User input

The user can tap, drag and rotate objects only inside the current 3D slice. An
object can move out of the slice. The slider lets the user find it again.

Each position that the user sees agrees with a position in 4D. Thus the
simulation can apply a force or an impulse at that position.

To rotate an object, find the 2D plane of rotation. Then form a bivector and
apply it as a torque. For the general case of one 4D frame to another 4D frame:

1. Find a rotor that makes the first frame equal to the second frame. Apply
   successive rotors that change one pair of coordinate axes at a time.
2. Extract the bivector part of that rotor.
3. Use the bivector as the torque.

## D3. Other feedback

- On a phone, use the accelerometer to control the direction of gravity. First
  keep the gravity vector inside the visible 3D slice. Then let a button turn
  the gravity vector to point a small quantity along the fourth axis.
- On a VR controller with haptic feedback, vibrate the controller when a held
  4D object hits another object.

---

# PART E — TEST PLAN

Test the library in this order. Each test uses the layer below it.

1. **Algebra identities.** Test `r ∧ a = [r]* a` and `r · w = [r]*^T w` with
   random data. Test `dI = [r]* [r]*^T`. Test `[R]2 B = R B R~`.
2. **Rotor closure.** Multiply many random rotors. Then factorize and rebuild.
   The rotor must stay a correct rotation.
3. **3D regression.** Set `n = 3`. Compare the results with a known 3D engine.
   All the physics must agree.
4. **Mass properties.** Calculate the inertia of a 4D hypercuboid from the
   mesh. Compare it with the analytic value.
5. **Conservation.** Run a torque-free body for a long time. The angular
   momentum and the kinetic energy must stay near constant.
6. **4D Dzhanibekov effect.** Set the angular velocities in a 3D subspace, that
   is `w_ic = w_ci = 0` for one index `c`. The equations must then reduce to
   the 3D case. You must see the usual unstable rotation about the intermediate
   plane.
7. **Double rotation stability.** Start a rotation in the `(x y)` plane with
   small perturbations in the other planes. A double rotation in the `(z w)`
   plane must start. The system must go to a stable state. The moments of
   inertia have little effect on this behavior.
8. **Collision.** Test a hypersphere against a hyperbox. Then test two
   hyperboxes. Then stack three hypercubes and check that the stack is stable.
9. **The sign of the error of a rotation.** Turn a body by 0.1 in the plane
   `p`. The error of A13 must give `+sin(0.1)` on that plane and 0 on each
   other plane. Then hold the body with a weld and solve one time: the angular
   velocity on that plane must become negative. The first test alone does not
   find a sign that is wrong only in the bias.
10. **A joint keeps the momentum.** Put the two anchors of a point joint at the
   same world point and solve. The linear momentum and the angular momentum of
   the pair must not change. A lever arm that is wrong breaks the angular part.
11. **A joint does no work.** With a bias of 0 and no gravity, the energy must
   never rise, and it must fall half as much when the step is half as long.
12. **A hinge in 4D.** A hinge on the plane 0 must leave that plane free and
   hold the other five. Give a torque in every plane at the same time.

---

# PART F — RISKS AND LIMITS

- **Rotor correction.** This is the most probable source of a defect. See B2.
  A norm is not enough.
- **Dual forms.** Do not take a dual anywhere. The dual between vectors and
  bivectors is correct only in 3D. If you find a dual in the code outside of
  the separating axis theorem, it is probably an error.
- **Sign conventions.** Fix the sign of `dR/dt` and the order of the commutator
  one time, at the start. Write them in a comment. The error of the rotation of
  A13 has the same risk: hold its sign with a test on the error, and with a
  second test on the correction that comes from it.
- **Basis order.** One lexicographic order for the bivector basis must be used
  in all modules. A different order in one table gives errors that are hard to
  find.
- **Boundary cases in collision detection.** In 4D the number of boundary
  dimensions is larger than in 3D. Test each dimension.
- **Not in the paper.** Soft bodies, continuous collision detection, and
  non-euclidean spaces. The paper gives these as future work.
- **The constraints and the joints** are not in the paper either. A13 and B8
  give the form that this library uses. Their limits:
  - A lock of the rotation that leaves **two or more planes free** has no
    geometric error that means anything, thus it adds up the drift of the
    velocity in the place of one. The drift is correct to the first order
    only. See A13 and B8.
  - `effectiveMass` does not use `linearFactor` and `angularFactor`, but the
    impulse does. A joint on a body with a factor of 0 does not converge. Use a
    subspace constraint in the place of a factor.
  - The angle of a hinge wraps at `pi`, thus a limit outside `(-pi, pi)` has no
    meaning.
  - A joint sends no report to the plugin. A third report would give a third
    stride to hold in agreement.