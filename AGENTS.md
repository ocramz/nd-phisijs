# Writing style

Only write in ASD-STE100 Simplified Technical English. This applies to the code
comments, to the documents and to the messages of the commits.

# Commands

```
npm test                  # node test/worker.js (21 tests), then the specs.
node --check physin.js    # the syntax only
node tools/sync-bundles.mjs         # copies the engine into physin_worker.js
node tools/sync-bundles.mjs --check # says what differs, and changes nothing
```

**Use `tools/sync-bundles.mjs` to keep the two bundles the same.** It finds each
module by its `// src/...` marker and puts it in the place of the module of the
same name. Do not copy a range of lines: the order of the modules is not the
same in the two files, thus a range copy deletes the module that comes next.
`test/worker.js` fails when a shared module differs, thus the error cannot go
out in silence.

There is no build, and there are no dependencies. 

`test/worker.js` is one file, and it has no filter for one test. It has three
blocks that the function `section()` starts. 

The example pages need no server. Open `sandbox.html`, `donuts.html` or
`plugin.html` in a browser, directly from the disk.

# Skills

## Discover skills

npx skills list -g

## Use skills

Example : `npx skills use  dubzzz/fast-check@javascript-testing-expert`

# Architecture

## The two bundles

| File | Content |
|---|---|
| `physin.js` |  three.js plugin. It puts `PhysiN` on the global object. |
| `physin_worker.js` | the same engine, for a web worker. |

**The two files hold the same engine modules, word for word.** The worker file
leaves out the slice and the plugin only. **Make each engine change in the two
files.** A `// src/...` comment starts each module. The order of the modules is
not the same in the two files, thus find a module by its `// src/...` line and
not by its position.

There is no `src/` directory. The bundles are the source.

## The modules

```
src/nd/core/dims.js         n, k, r, the bivector basis order, the index tables
src/nd/core/linalg.js       small dense matrices
src/nd/algebra/*.js         PART A: multivector, rotor, star matrix
src/nd/body/*.js            PART A: the state, the shapes, the mass properties
src/nd/detect/*.js          the nearest point, the collision
src/nd/integrate/*.js       PART B: the time step, the gyroscopic term
src/nd/resolve/solver.js    PART B: the impulses, the friction, the shock propagation
src/nd/resolve/constraint.js  PART B: the rows, the five joints
src/nd/world.js             `World`, and `defaultParams` (all the tolerances)
src/slice/slice4.js         4D only: the cut of the surface with a 3D hyperplane
src/nd/workerCore.js        the engine behind the message protocol
src/physiN.js               the three.js plugin (in `physin.js` only)
```

`PhysiN.nd` gives the engine, and it has no dependency on three.js. `PhysiN.slice`
gives the slice.

## The rules of the design

1. **Part A and Part B stay apart.** The algebra modules (Part A) must hold no
   tolerance. The numerical modules (Part B) must hold no new algebraic rule.
   `ND-PHYSICS.md` holds the mathematics, and the code points to it by the
   number of a part, for example `A3` or `B2`.
2. **All the tolerances are in `defaultParams`**, in `src/nd/world.js`. Do not
   write a new constant in an other module.
3. **The types of `n` dimensions.** Position, velocity, force and impulse are
   vectors of `n` components. The orientation is a **rotor** of `2^(n-1)`
   components. The angular velocity, the torque and the angular momentum are
   **bivectors** of `k = n(n-1)/2` components. The inertia is a `k` x `k`
   matrix. There is no cross product, and there is no axial vector. The order
   of the bivector components is lexicographic: for `n = 4`, index 0 is `(x y)`
   and index 2 is `(x w)`.
4. **The rows of a joint are named by their position.** A count that changes
   throws away every impulse that the joint keeps. Thus the row of a motor and
   the row of a limit are ALWAYS there, and they have no mass when they are
   off. Do not build a row only when it is needed.
5. **A joint holds planes, and not axes.** A hinge leaves one rotation
   **plane** free. In 3D that is 1 plane of 3, and in 4D it is 1 of 6. The
   angular rows of a joint are bivectors, thus a hinge in 4D holds 5 planes.
   The error of a rotation comes from the antisymmetric part of the relative
   rotation matrix, and it is only exact when the free part is one plane or
   nothing. See `angularError`, and `A13`.
6. **`rotorCorrect` is necessary.** In 4D the norm of a rotor can be exactly 1
   while the rotor is not a rotation. Only a rotor that is a product of simple
   rotors is correct. Without the correction the energy increases without
   limit. See section 10 of the README, and `B2`.

## The message protocol

The plugin and the engine speak with two message types:

- A command object, `{ cmd, params }`. The commands are in `createEngine()` in
  `src/nd/workerCore.js`.
- A binary report, a `Float32Array`. The first number gives the type.

```
world report,     stride = 1 + n + r + n + k
  id | position (n) | rotor (r) | linear velocity (n) | angular velocity (k)
collision report, stride = 2 + n + n + 1
  id of a | id of b | normal (n) | point (n) | depth
```

**`workerCore.js` writes the reports and `PhysiN.Scene` reads them. The two
calculations of the stride must always agree.** If you change a report, change
the two.

A joint sends no report, and there is no third stride. `addConstraint`,
`removeConstraint` and `setConstraintParams` are command objects only. A joint
that the engine cannot build sends `constraintFailed` back, and the plugin
writes that to the console. A joint that breaks sends `constraintBroken`, and
the plugin sends the event `broken` on the joint. A command must never go away
without a message.

The same `createEngine()` code runs in the two conditions. With
`PhysiN.scripts.worker = null` the plugin calls `handle()` directly, in the main
thread. With a worker file, `postMessage` carries the messages. Thus the results
must be the same in the two conditions, and `test/worker.js` compares them.

## The plugin

The plugin does no physics. It keeps the `n`-dimensional state of an object in
`object._physiN`, it sends the commands to the engine, and it puts the answer
back onto the three.js objects. The display depends on `n`: at `n = 3` the rotor
becomes a quaternion; at `n = 4` the slice builds a new 3D surface at each frame;
at `n > 4` the plugin draws the first three axes of the position only.

`setPositionN` and `rotateInPlane` write into the mesh, thus you can call them
before `scene.add`. Each other method sends a command to the engine, thus you
must call it **after** `scene.add`. A command that you send before `scene.add`
goes away, and the library gives no message.

# Tests

`test/worker.js` reads the two bundles from the disk, then it connects them with
a `Worker` that is a stub and a `THREE` that is a stub. Thus it tests the
protocol from end to end, in Node, with no browser. It tests the logic, and not
the appearance.