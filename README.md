# PhysiN

N-dimensional rigid body dynamics for three.js.

PhysiN gives the interface of [Physijs](https://github.com/chandlerprall/Physijs),
but it operates in any number of dimensions `n > 1`. It has special cases for
3D and for 4D. The physics is new code. It is not a wrapper of ammo.js,
because ammo.js is correct only in 3D.

The mathematics comes from Marc ten Bosch, *N-Dimensional Rigid Body Dynamics*,
ACM Transactions on Graphics 39(4), 2020.

---

## 1. What it does

- Rigid body dynamics for `n = 2, 3, 4, 5, 6` or more.
- Collision of hyperboxes, hyperspheres, tori and half spaces.
- Stacks of bodies, with a contact graph and shock propagation.
- Friction and restitution.
- A three.js plugin, with the same interface as Physijs.
- A web worker, or the main thread. The same engine code operates in both.
- A 4D slice, to show a 4D body on the screen.

---

## 2. Look at it first

Two pages need no build. Each page is one file, and it gets three.js from a
CDN. Open them in a browser.

- `examples/sandbox.html` — bodies that fall, in 3 and in 4 dimensions.
- `examples/plugin.html` — the five steps of the interface, with the two files.
- `examples/link.html` — two donuts. Move one donut along `x` and along `w` to
  make a chain link. You do not cut anything. In three dimensions this is not
  possible.

## 3. The two files, and the five steps

The build gives two files. Put them on your page.

| File | Use |
|---|---|
| `dist/physiN.js` | the three.js plugin. One classic script. It contains all of the library. |
| `dist/physiN_worker.js` | the same engine, for a web worker. One classic worker script. |

The interface follows physi.js:

1. Load three.js. Then load `physiN.js` with a `<script>` tag.
2. Point `PhysiN.scripts.worker` at the worker file.
3. Use `PhysiN.Scene` in the place of `THREE.Scene`.
4. Use `PhysiN.BoxMesh`, `PhysiN.SphereMesh`, `PhysiN.ConvexMesh`,
   `PhysiN.TorusMesh` or `PhysiN.PlaneMesh` in the place of `THREE.Mesh`.
5. Call `scene.simulate()` at each frame.

```html
<script src="three.min.js"></script>
<script src="physiN.js"></script>
<script>
  PhysiN.scripts.worker = 'physiN_worker.js';   // null = the main thread

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

The page `examples/plugin.html` gives these five steps in full.

**Note on the worker.** The default worker type is `classic`, because the two
files are classic scripts. Use the ES module build with
`new PhysiN.Scene({ workerType: 'module' })`.

## 4. Install

```
npm install
npm test           # 59 tests
npm run build      # makes the bundles in dist/
npm run serve      # then open http://localhost:8080/examples/4d.html
```

You must use a web server for the examples. A browser does not load an ES
module from a `file://` address.

### The files

```
src/nd/            the general library. It does not know about three.js.
  core/            dimension tables, linear algebra
  algebra/         multivector, rotor, star matrix
  body/            mass properties, shapes, the body state
  integrate/       the time step, the gyroscopic term, the rotor correction
  detect/          collision detection, the separating axis theorem
  resolve/         impulses, friction, the contact graph
  world.js         the world, the broad phase, the loop
  workerCore.js    the message protocol
src/physiN.js      the three.js plugin
src/slice/         the 4D slice, for the display
examples/          browser pages. `sandbox.html` is one file, with no build.
test/              the test plan
dist/              the bundles
```

---

## 5. Quick start, 3 dimensions

```js
import * as THREE from 'three';
import { createPhysiN } from 'physin';

const PhysiN = createPhysiN(THREE);
PhysiN.scripts.worker = new URL('physin/src/physiN_worker.js', import.meta.url).href;

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

The page `examples/4d.html` gives a ruler at the bottom of the screen. The
ruler shows the interval of `w` that each body covers. Drag the ruler to move
the hyperplane. A body that the hyperplane cuts has a colored band.

---

## 7. The web worker

```js
// Physics in a worker (the default):
PhysiN.scripts.worker = new URL('../src/physiN_worker.js', import.meta.url).href;

// Physics in the main thread:
PhysiN.scripts.worker = null;
```

The same file `src/nd/workerCore.js` operates in both conditions. Thus the
results are the same. Use the main thread to debug. Use the worker for a game.

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
| `new PhysiN.Scene(params)` | `params`: `dimensions`, `gravity`, `sliceW`, `params` |
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
| `.setMass(m)` | A mass of 0 makes the body static. |
| Events | `ready`, `collision` |

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
The correction is in `src/nd/algebra/rotor.js`, in the function `rotorCorrect`.
The integrator calls it when the defect is more than `rotorTolerance`.

The gyroscopic term makes double rotations very frequently. Thus you cannot
remove this step from a 4D engine.

---

## 11. The library without three.js

The directory `src/nd/` has no dependency. Use it for a server, for a test, or
for a different renderer.

```js
import { World, HyperBox, HalfSpace, Body, dims } from 'physin/nd';

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

All the tolerances are in one place: `defaultParams` in `src/nd/world.js`.
Give new values in `new PhysiN.Scene({ params: { ... } })`.

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
| `rotorTolerance` | 1e-9 | the limit of the rotor defect |
| `gyroscopic` | true | set it to false to remove the gyroscopic term |
| `sleepLinearVelocity`, `sleepAngularVelocity`, `sleepTime` | 0.03, 0.03, 0.6 | the sleep limits |
| `useShockPropagation`, `useWarmStart` | true | set to false to compare |

The mathematics module has no tolerance. The numerics module has no new
algebraic rule. Keep this separation if you change the code.

---

## 13. Tests

```
node test/run.js       # 46 tests of the library
node test/plugin.js    # 13 tests of the three.js plugin
node test/artifact.js  # 15 tests of the sandbox page
node test/link.js      # 14 tests of the chain link page
node test/worker.js    # 13 tests of the two files and the worker protocol
```

The tests follow the test plan in this order:

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
in Node, with no browser.

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
- There are no constraints and no joints.
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