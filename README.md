# PhysiN

N-dimensional rigid body dynamics for three.js.

PhysiN gives the interface of [Physijs](https://github.com/chandlerprall/Physijs),
but it operates in any number of dimensions `n > 1`. It has special cases for
3D and for 4D. The physics is new code. It is not a wrapper of ammo.js,
because ammo.js is correct only in 3D.

The mathematics comes from Marc ten Bosch, *N-Dimensional Rigid Body Dynamics*,
ACM Transactions on Graphics 39(4), 2020.

---

## 0. What is in this repository

| File | Content |
|---|---|
| `physin.js` | the three.js plugin. One classic script. It contains all of the library. |
| `physin_worker.js` | the same engine, for a web worker. One classic worker script. |
| `sandbox.html` | a page of bodies that fall, in 3 and in 4 dimensions. |
| `donuts.html` | a page of two donuts that make a chain link. |
| `plugin.html` | the five steps, on a page that loads the two files above. |
| `test/worker.js` | the test of the two files and of the message protocol. |
| `ND-PHYSICS.md` | the mathematics, in full. |


---

## 1. What it does

- Rigid body dynamics for `n = 2, 3, 4, 5, 6` or more.
- Collision of hyperboxes, hyperspheres, tori and half spaces.
- Stacks of bodies, with a contact graph and shock propagation.
- Friction and restitution.
- Constraints and joints: a point joint, a distance joint, a weld, a hinge
  and a lock to a subspace.
- A three.js plugin, with the same interface as Physijs.
- A web worker, or the main thread. The same engine code operates in both.
- A 4D slice, to show a 4D body on the screen.

---

## 2. Look at it first

Three pages need no build, and no server. Each page gets three.js from a CDN.
Open a page in a browser, directly from the disk.

- `sandbox.html` — bodies that fall, in 3 and in 4 dimensions. A ruler at the
  bottom of the screen moves the 3D hyperplane along the `w` axis.
- `donuts.html` — two donuts. Move one donut along `x` and along `w` to
  make a chain link. You do not cut anything. In three dimensions this is not
  possible.
- `plugin.html` — a stack of boxes, a ball that hits the stack, and an
  octahedron. This page is the five steps of section 3, and it is short.

`sandbox.html` and `donuts.html` each contain the engine, and each one sets
`PhysiN.scripts.worker = null`. Thus they do the physics in the main thread.

`plugin.html` is different: it loads `physin.js` with a `<script>` tag, from
the same folder. Thus it shows how to put the library on your own page. A
browser does not permit a worker from the disk. Thus the page uses the main
thread for a `file:` address, and the worker if a server gives the page.

## 3. The two files, and the five steps

Put these two files on your page.

| File | Use |
|---|---|
| `physin.js` | the three.js plugin. One classic script. It contains all of the library. |
| `physin_worker.js` | the same engine, for a web worker. One classic worker script. |

The interface follows physi.js:

1. Load three.js. Then load `physin.js` with a `<script>` tag.
2. Point `PhysiN.scripts.worker` at the worker file.
3. Use `PhysiN.Scene` in the place of `THREE.Scene`.
4. Use `PhysiN.BoxMesh`, `PhysiN.SphereMesh`, `PhysiN.ConvexMesh`,
   `PhysiN.TorusMesh` or `PhysiN.PlaneMesh` in the place of `THREE.Mesh`.
5. Call `scene.simulate()` at each frame.

```html
<script src="three.min.js"></script>
<script src="physin.js"></script>
<script>
  PhysiN.scripts.worker = 'physin_worker.js';   // null = the main thread

  var scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
  scene.add(new PhysiN.PlaneMesh(groundGeometry, stone, [0, 1, 0], 0));

  var box = new PhysiN.BoxMesh(new THREE.BoxGeometry(1, 1, 1), teal, 1);
  box.setPositionN([0, 5, 0]);
  scene.add(box);

  function frame() {
    requestAnimationFrame(frame);
    scene.simulate(1 / 60, 2);
    renderer.render(scene, camera);
  }
  frame();
</script>
```

`plugin.html` is this code, as a page that operates.

**Note on the global.** `physin.js` is a classic script, and not an ES module.
It has no exports. It puts `PhysiN` on the global object.

**Note on the worker.** The default value of `PhysiN.scripts.worker` is `null`.
Thus the physics operates in the main thread until you set that property. The
default worker type is `classic`, because the two files are classic scripts.
`new PhysiN.Scene({ workerType: 'module' })` selects a module worker, but you
must then supply an ES module build. This repository does not contain one.

## 4. Install and test

The library has no dependencies, and it needs no build. The two files are the
build.

```
npm test          # 20 tests of the two files and the worker protocol, then the specs
```

The example pages need no server. Open `sandbox.html`, `donuts.html` or
`plugin.html` directly from the disk.

### The files

```
physin.js          the three.js plugin, and all of the library
physin_worker.js   the same engine, for a web worker
sandbox.html       bodies that fall, in 3 and in 4 dimensions
donuts.html        two donuts that make a chain link
plugin.html        the five steps, on a page that loads the two files
test/worker.js     the test of the two files and the worker protocol
ND-PHYSICS.md      the mathematics, in full
package.json       the test command
```

---

## 5. Quick start, 3 dimensions

Load three.js and then `physin.js`, as in section 3. The two scripts give
`THREE` and `PhysiN` as globals.

```js
PhysiN.scripts.worker = 'physin_worker.js';   // null = the main thread

const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });

const ground = new PhysiN.PlaneMesh(
  new THREE.BoxGeometry(40, 0.2, 40), material, [0, 1, 0], 0);
scene.add(ground);

const box = new PhysiN.BoxMesh(new THREE.BoxGeometry(1, 1, 1), material, 1);
box.setPositionN([0, 5, 0]);
scene.add(box);

function frame() {
  requestAnimationFrame(frame);
  scene.simulate(1 / 60, 2);
  renderer.render(scene, camera);
}
```

This is the interface of Physijs. If you know Physijs, you know this.

---

## 6. Quick start, 4 dimensions

A 4D body has a 3-dimensional surface. You cannot show that surface directly.
Thus PhysiN cuts the surface with one 3D hyperplane. The cut is a 2D surface,
and three.js can show it. Move the hyperplane along the fourth axis `w` to
see the other parts of the body.

```js
const scene = new PhysiN.Scene({ dimensions: 4, gravity: [0, -9.81, 0, 0] });

const floor = new PhysiN.HyperPlaneMesh(
  [0, 1, 0, 0], 0, material, new THREE.BoxGeometry(30, 0.2, 30));
scene.add(floor);

const cube = new PhysiN.HyperBoxMesh([0.5, 0.5, 0.5, 0.5], material, 1);
cube.setPositionN([0, 3, 0, 0]);
cube.rotateInPlane(0, 3, 0.4);        // turn 0.4 rad in the (x w) plane
scene.add(cube);

scene.sliceW = 0;                     // the position of the hyperplane
```

Set `scene.sliceW` at each frame, or call `scene.refreshSlice()` after you
change it. The mesh of each 4D body then changes to the new cut. If a body is
not in the slice, PhysiN sets `visible = false`.

The page `sandbox.html` gives a ruler at the bottom of the screen. The
ruler shows the interval of `w` that each body covers. Drag the ruler to move
the hyperplane. A body that the hyperplane cuts has a colored band.

---

## 7. The web worker

```js
// Physics in a worker:
PhysiN.scripts.worker = 'physin_worker.js';

// Physics in the main thread (the default):
PhysiN.scripts.worker = null;
```

Give the address of `physin_worker.js` as the browser sees it. The plugin gives
that value to `new Worker()` with no change.

The two files contain the same engine code. Thus the results are the same in
the two conditions. Use the main thread to debug. Use the worker for a game.
`test/worker.js` compares the two conditions.

The worker sends a `Float32Array` report at each step. The report has one
record for each body. The record has this content:

```
id | position (n) | rotor (2^(n-1)) | linear velocity (n) | angular momentum (k)
```

---

## 8. Interface

### PhysiN.Scene

| Member | Description |
|---|---|
| `new PhysiN.Scene(params)` | `params`: `dimensions`, `gravity`, `sliceW`, `workerType`, `params` |
| `.add(object)` | Adds a mesh. A `THREE.Object3D` with no body is also correct. |
| `.remove(object)` | Removes the mesh and its body. |
| `.simulate(timeStep, maxSubSteps)` | Makes one step. |
| `.setGravity(g)` | `g` is an array of `n` numbers. |
| `.sliceW` | The position of the 3D hyperplane on the `w` axis. |
| `.refreshSlice()` | Calculates the slice again. |
| `.execute(cmd, params)` | Sends a command to the engine. |
| Events | `ready`, `update` |

### PhysiN.Mesh

| Member | Description |
|---|---|
| `.setPositionN(p)` | `p` is an array of `n` numbers. |
| `.getPositionN()` | Gives the position in `n` dimensions. |
| `.getRotor()` | Gives the rotor, `2^(n-1)` numbers. |
| `.rotateInPlane(i, j, angle)` | Turns the body in the plane of the axes `i` and `j`. |
| `.setLinearVelocity(v)` | `v` has `n` components. |
| `.setAngularVelocity(w)` | `w` is a **bivector**. It has `k` components. |
| `.applyCentralImpulse(j)`, `.applyImpulse(j, r)` | `j` has `n` components. |
| `.applyCentralForce(f)`, `.applyForce(f, r)` | `f` has `n` components. |
| `.applyTorque(t)` | `t` is a **bivector**. It has `k` components. |
| `.getLinearVelocity()`, `.getAngularVelocity()` | Give the last values from the engine. |
| Events | `ready`, `collision` |

`setPositionN` and `rotateInPlane` write into the mesh. Thus you can call them
before `scene.add`. Each other method in the table is a command to the engine,
and the mesh sends it through its parent. Thus you must call it **after**
`scene.add`. A command that you send before `scene.add` goes away, and the
library gives no message.

`PhysiN.Mesh` has no `setMass` method. Give the mass to the constructor. A mass
of 0 makes the body static. To change the mass after that, send the command to
the engine:

```js
scene.execute('setMass', { id: mesh._physiN.id, value: 2 });
```

### The shapes

For 3 dimensions, the mesh gets its size from the geometry:

- `new PhysiN.BoxMesh(geometry, material, mass)`
- `new PhysiN.SphereMesh(geometry, material, mass)`
- `new PhysiN.PlaneMesh(geometry, material, normal, offset)`
- `new PhysiN.TorusMesh(geometry, material, mass, majorRadius, minorRadius, plane)`
- `new PhysiN.ConvexMesh(geometry, material, mass)`

For 4 dimensions and more, you give the size directly:

- `new PhysiN.HyperBoxMesh(halfExtents, material, mass)`
- `new PhysiN.HyperSphereMesh(radius, material, mass)`
- `new PhysiN.HyperPlaneMesh(normal, offset, material, displayGeometry)`
- `new PhysiN.HyperTorusMesh(majorRadius, minorRadius, plane, material, mass, options)`
- `new PhysiN.HyperMesh(vertices, cells, material, mass)`

### The convex mesh

`PhysiN.ConvexMesh` reads the triangles of the geometry. It removes the
repeated vertices, then it calculates the volume, the center and the inertia
from those triangles.

The geometry must be **convex**, and the **origin must be inside it**. A shape
that is not convex gives wrong values, and the library cannot find that error.

A convex body touches a half space and a hypersphere. It does not touch a
hyperbox or a second convex body.

### The torus

The solid `n`-torus is the set of points at a distance of `r` or less from a
circle of radius `R`. `plane` gives the two axes of that circle.

```
( sqrt(x_i^2 + x_j^2) - R )^2  +  SUM_(other k) x_k^2  <=  r^2
```

In 3D this is the usual donut. In 4D it is a 4-dimensional solid. Its cross
section is a 3-ball, and not a disk.

The slice of a 4D torus with the hyperplane `x_w = value` is a **3D torus**.
It has the same major radius, and a smaller minor radius:

```
minor radius of the slice = sqrt(r^2 - (value - center_w)^2)
```

The plugin calculates that slice directly. It is exact, and it is much faster
than a cut of a tetrahedral mesh. It uses the mesh only if the body turns out
of the slice, or if the major plane uses the `w` axis.

The volume and the inertia are analytic for each `n`:

```
volume = 2 pi R V(n-1, r)          V(d, r) = the volume of a d-ball
C_ii = C_jj = m (R^2 + 3 r^2 / (n + 1)) / 2      i, j = the major plane
C_kk         = m r^2 / (n + 1)                    each other axis
```

`C` is the covariance matrix. `inertiaFromCovariance` changes it into the
`k` x `k` inertia tensor. For `n = 3` these values give the classical results
`m (R^2 + 3 r^2 / 4)` about the symmetry axis and `m (4 R^2 + 5 r^2) / 8`
about a diameter. The tests use this.

**Limit:** a torus touches a half space and a hypersphere. It does not touch a
hyperbox or a second torus. The library gives a warning one time if you try.

### Friction and restitution

Put the values on the three.js material:

```js
material._physiN = { friction: 0.7, restitution: 0.05 };
```

### The joints

A joint holds two bodies together, or it holds one body to the world. Give a
null second object for a joint to the world.

**Add the two objects to the scene first.** The engine builds a joint from the
ids of the two bodies. `addConstraint` throws if an object is not in the scene.

```js
scene.add(post);
scene.add(bob);
const joint = new PhysiN.PointJoint(post, bob, [0, 0, 0], [-2, 0, 0]);
scene.addConstraint(joint);
// later
scene.removeConstraint(joint);
```

**The anchors are always in the local frame of their body.** In the example the
anchor of `bob` is 2 units to its left, thus that point of `bob` holds the
center of `post`, and `bob` hangs 2 units away.

| Class | It holds | It takes away, n = 3 | It takes away, n = 4 |
|---|---|---|---|
| `PointJoint(a, b, localA, localB)` | one point on one point | 3 of 6 | 4 of 10 |
| `DistanceJoint(a, b, localA, localB, opts)` | the length between two points | 1 of 6 | 1 of 10 |
| `FixedJoint(a, b, localA, localB)` | the point and every plane | 6 of 6 | 10 of 10 |
| `HingeJoint(a, b, localA, localB, opts)` | the point, and every plane but one | 5 of 6 | 9 of 10 |
| `SubspaceJoint(a, b, opts)` | the axes and the planes that you name | | |

`DistanceJoint` takes `opts.rest` (the length; the default is the length at the
start) and `opts.mode`: `"rod"` holds the length, `"rope"` only stops it from
growing, and `"strut"` only stops it from falling.

**A hinge leaves one rotation PLANE free, and not one axis.** In 3 dimensions
that is the door hinge that you know: 1 plane of 3. In 4 dimensions it is 1
plane of 6, thus a 4D hinge holds five planes. `opts.plane` names the free
plane in the frame of `a`: a number of the lexicographic order, or a bivector
of `k` components. See section 9.

```js
// a 4D hinge that turns only in the (x y) plane
scene.addConstraint(new PhysiN.HingeJoint(a, b, [1, 0, 0, 0], [-1, 0, 0, 0], { plane: 0 }));
```

The plane must be **simple**: `e_xy + e_zw` holds two planes at the same time
and it is not a hinge, thus the engine refuses it.

**The limits and the motor of a hinge.** `opts.lowerAngle` and
`opts.upperAngle` give the two limits in radians, and `opts.motorSpeed` and
`opts.maxMotorTorque` give the motor. Change them while the world runs with
`setLimits` and `setMotor`.

```js
const hinge = new PhysiN.HingeJoint(base, arm, [0, 0, 0], [-1, 0, 0], {
  plane: 0, lowerAngle: -0.5, upperAngle: 0.5,
});
scene.addConstraint(hinge);
hinge.setMotor(2, 50);       // 2 rad/s, at 50 Nm at the most
```

The motor takes a **torque** and not an impulse, thus its strength does not
change when you change `subSteps`. A motor that drives into a limit stalls
there.

The angle of a hinge wraps at `pi`. A limit outside `(-pi, pi)` has no meaning.

`SubspaceJoint` holds a body in a subspace. Give `opts.lockAxes` or
`opts.freeAxes` for the position, and `opts.lockPlanes` or `opts.freePlanes`
for the rotation. `opts.worldFrame` at true holds the directions in the world
frame; the default turns them with `a`.

```js
// hold a 4D body on the hyperplane w = 0, and leave x, y and z free
scene.addConstraint(new PhysiN.SubspaceJoint(body, null, { lockAxes: [3], worldFrame: true }));
```

Use a subspace joint in the place of `linearFactor` and `angularFactor`. Those
two masks work on the velocity alone, and they do not correct a drift.

**A joint is soft under a load.** The error of a joint grows with the length of
a chain and with the ratio of the masses. `iterations` and `subSteps` do
different work, and which one helps depends on the scene:

| A chain of 10 links, the worst error over 10 s | it 10 | it 100 | 4 sub-steps |
|---|---|---|---|
| it starts horizontal, thus it whips | 0.114 | 0.115 | 0.014 |
| it starts hanging, thus it settles | 0.0014 | 0.00003 | 0.00009 |

`iterations` brings the **velocity** of a row to its target, thus it helps a
chain that has come to rest. `subSteps` builds the rows again, thus it makes
the error of the **position** smaller; that is the only thing that helps a
chain that moves fast. Use `subSteps`, or a shorter `fixedTimeStep`, for a
scene that moves.

A ratio of the masses of 1000 to 1 needs sub-steps. At one step of 1/120 that
chain comes apart; with `subSteps: 8` its error is 0.45.

**A soft joint.** Give `opts.hertz` and `opts.damping` to make a joint a spring
of that frequency. The stiffness then does not change with the length of the
step, and the joint cannot put energy in. The stretch under a weight is the
stretch of a spring, `g / (2 pi f)^2`: a mass of 1 kg on a joint of 5 Hz hangs
9.94 mm low. A `hertz` of 0, the default, makes the joint rigid.

```js
joint.setSoftness(5, 1);     // 5 Hz, and no overshoot
joint.setSoftness(0, 1);     // rigid again
```

**A joint that breaks.** `opts.breakForce` and `opts.breakTorque` take newtons
and newton metres. A value of 0, the default, means that the joint never
breaks. A joint that breaks sends the event `broken` and it leaves the scene.
The two bodies can touch each other after that, and they usually overlap.

```js
joint.setBreak(500, 0);
joint.addEventListener('broken', () => console.log('it came apart'));
```

The limit of a hinge counts toward the load, and the motor does not. Thus a
motor never breaks its own joint, but a motor that pushes against a limit puts
its torque through the limit and can break a joint that is not strong enough.

The two bodies of a joint do not touch each other. Give
`opts.collideConnected: true` to let them touch.

---

## 9. The conventions of `n` dimensions

**This is the part that is different from a 3D engine. Read it.**

The cross product is correct only in 3D. PhysiN does not use it, and it does
not use axial vectors. It uses these types in their place:

| Quantity | 3D type | PhysiN type | Component count |
|---|---|---|---|
| Position, velocity, force, impulse | vector | vector | `n` |
| Orientation | quaternion | rotor | `2^(n-1)`, that is `r` |
| Angular velocity, torque, angular momentum | axial vector | **bivector** | `n(n-1)/2`, that is `k` |
| Inertia | 3x3 matrix | matrix | `k` x `k` |

For `n = 4`: `k = 6` and `r = 8`.

The order of the bivector components is lexicographic. For `n = 4`:

```
index 0: (x y)   index 1: (x z)   index 2: (x w)
index 3: (y z)   index 4: (y w)   index 5: (z w)
```

Thus a torque in the `(x w)` plane is `[0, 0, 9, 0, 0, 0]`.

A rotation is not a rotation "about an axis". It is a rotation "in a plane".
In 3D each plane has one perpendicular axis, and the two ideas agree. In 4D
they do not agree. A 4D body can turn in two perpendicular planes at the same
time.

If you set `dimensions: 3`, the rotor becomes a quaternion and the bivector
becomes the usual angular velocity. Thus the library gives a usual 3D engine
as a special case. The tests use this property.

---

## 10. Why the rotor correction is necessary

**This is the most probable source of a defect in a 4D engine.**

In 3D you correct a quaternion when you normalize it. This is not sufficient
in 4D.

A simple rotor is the geometric product of two unit vectors. Simple rotors are
on a sphere in the even sub-algebra. A rotor of a double rotation is not on
that sphere. Thus the norm can be exactly 1.0 while the rotor is not a
rotation.

The test shows this. Without the correction, and with a 4D double rotation:

```
step    |R|          orthonormality error    kinetic energy
0       1.0000000    2.5e-05                 1.00
50      1.0000042    1.45e-03                1.00
400     1.0000066    grows                   3.3e+143
```

The norm gives no warning. The energy increases without limit.

PhysiN makes the frame of the rotor orthonormal, then it builds the rotor
again as a product of simple rotors. The result is always a correct rotation.
The correction is the function `rotorCorrect`, and `PhysiN.nd.rotor` gives it.
The integrator calls it when the defect is more than `rotorTolerance`.

The gyroscopic term makes double rotations very frequently. Thus you cannot
remove this step from a 4D engine.

---

## 11. The library without three.js

`physin.js` gives the general library as `PhysiN.nd`. That part has no
dependency on three.js. Use it for a server, for a test, or for a different
renderer. `physin_worker.js` uses the same code.

```js
const { World, HyperBox, HalfSpace, Body, dims } = PhysiN.nd;

const D = dims(4);                       // the tables for n = 4
const world = new World({ dimensions: 4, gravity: [0, -9.81, 0, 0] });

world.addBody(new Body(D, { shape: HalfSpace(D, [0, 1, 0, 0], 0), mass: 0 }));

const body = new Body(D, {
  shape: HyperBox(D, [0.5, 0.5, 0.5, 0.5]),
  mass: 1,
  position: [0, 3, 0, 0],
});
world.addBody(body);

for (let i = 0; i < 120; i += 1) world.step(1 / 60);
console.log(body.x, body.w, body.kineticEnergy());
```

`body.x` is the position, `n` numbers. `body.R` is the rotor, `r` numbers.
`body.v` is the linear velocity, `n` numbers. `body.w` is the angular velocity
bivector, `k` numbers. `body.L` is the angular momentum bivector.

---

## 12. Tolerances

All the tolerances are in one place: `PhysiN.nd.defaultParams`. Give new values
in `new PhysiN.Scene({ params: { ... } })`.

| Name | Default | Use |
|---|---|---|
| `fixedTimeStep` | 1/60 | the length of one step |
| `subSteps` | 1 | the number of sub-steps in one step |
| `iterations` | 10 | the iteration count of the contact solver |
| `shockIterations` | 2 | the iteration count of the shock propagation |
| `penetrationSlop` | 0.005 | the permitted penetration |
| `biasFactor` | 0.2 | the correction of the position error |
| `contactMargin` | 0.02 | the distance that makes a contact |
| `restitutionThreshold` | 0.5 | the minimum speed for a bounce |
| `maxContacts` | 0 | the limit of the contact count. 0 = no limit |
| `constraintBias` | 0.2 | the correction of the error of a joint |
| `constraintSlop` | 0.001 | the error of a joint that the solver accepts |
| `constraintMaxBias` | 10 | the largest speed that the bias of a joint asks for |
| `constraintTolerance` | 1e-6 | the shortest direction that a joint accepts |
| `constraintHertz` | 0 | the frequency of a soft joint, in Hz. 0 = rigid |
| `constraintDamping` | 1 | the damping ratio of a soft joint |
| `constraintHertzRatio` | 0.25 | the largest part of the rate of a substep that the frequency may use |
| `rotorTolerance` | 1e-9 | the limit of the rotor defect |
| `gyroscopic` | true | set it to false to remove the gyroscopic term |
| `gyroscopicIterations` | 1 | the iteration count of the implicit gyroscopic step |
| `allowSleep` | true | set it to false to keep all the bodies awake |
| `sleepLinearVelocity`, `sleepAngularVelocity`, `sleepTime` | 0.03, 0.03, 0.6 | the sleep limits |
| `useShockPropagation`, `useWarmStart` | true | set to false to compare |

The mathematics module has no tolerance. The numerics module has no new
algebraic rule. Keep this separation if you change the code.

---

## 13. Tests

This repository has one test file.

```
npm test               # the same as: node test/worker.js
```

`test/worker.js` gives 20 tests of the two files and of the message protocol.
It reads `physin.js` and `physin_worker.js` from the disk. Then it connects
them with a stub Worker. Thus it tests the protocol from end to end. It also
compares the main thread with the worker: a box comes to rest at the same
height in the two conditions. Sections 3 and 5 use four dimensions, and they
include the slice of a 4D torus and a subspace joint. Section 4 hangs a
pendulum on a point joint.

The other tests of the test plan operate on the source tree, and they are not
in this repository. The full plan has this order:

1. Algebra identities: `r ∧ a = [r]* a`, `r · w = [r]*ᵀ w`, `dI = [r]* [r]*ᵀ`,
   `[R]₂ B = R B R~`, for `n = 2` to `n = 6`.
2. Rotor closure: multiply 200 random rotors, then correct the result.
3. 3D regression: compare with the known 3D equations.
4. Mass properties: compare the mesh result with the analytic inertia.
5. Conservation: a torque-free 4D body keeps its energy.
6. The 4D Dzhanibekov effect: a 3D subspace stays closed.
7. Double rotation stability.
8. Collision: a hypersphere on a hyperbox, and a stack of three hypercubes.
9. Torus: the volume and the inertia against the classical 3D results and
   against the mesh method, the slice, and a torus in the physics.
10. Convex mesh: the volume of the cross polytope, a box as a convex mesh,
    the nearest point on a simplex, and the contact with a hypersphere.

The plugin tests use a small stub in the place of three.js. Thus they operate
in Node, with no browser. `test/worker.js` contains that stub.

### Known effect

The implicit method of the gyroscopic term removes a small quantity of angular
momentum over a long time. At `dt = 1/120` the loss is near 7%. At `dt = 1/960`
it is near 1%. This is a property of the method, and it is not a defect in the
algebra. Use a smaller time step, or set `gyroscopic: false` if your bodies
have no fast free rotation.

---

## 14. Speed

Measured with 9 boxes in a stack, on one core:

| `n` | Time for one step | Narrow phase |
|---|---|---|
| 3 | 0.86 ms | near 0 ms |
| 4 | 1.39 ms | 0.25 ms |

In 4D the separating axis theorem uses 56 axes. In 3D it uses 15. This is the
largest part of the 4D cost.

The rotor has `2^(n-1)` components. Thus the cost increases quickly with `n`.
Do not plan for a large `n`.

---

## 15. Limits

- Collision operates for hyperbox, hypersphere, torus and half space. A torus
  touches a half space and a hypersphere only. A general convex body gives
  mass properties and a display, but not a collision.
- A joint that leaves **two or more rotation planes free** has no geometric
  error that means anything, because a turn in one plane and a turn in another
  do not commute. The engine adds up the drift of the velocity in the place of
  one. That is correct to the first order only.
- The angle of a hinge wraps at `pi`. A limit outside `(-pi, pi)` has no
  meaning, and a hinge that turns past `pi` reports a jump.
- `linearFactor` and `angularFactor` do not work together with a joint on the
  same body. Use a `SubspaceJoint` in the place of the factor.
- There are no soft bodies.
- There is no continuous collision detection. A fast body can go through a
  thin body. Use more sub-steps.
- There are no vehicles.
- The examples are not tested in a browser. The tests use stubs in the place
  of three.js and the DOM. Thus they test the logic, and not the appearance.

---

## 16. License

MIT.

The interface follows Physijs by Chandler Prall (MIT). The physics code is
new. It is an implementation of the paper of Marc ten Bosch.