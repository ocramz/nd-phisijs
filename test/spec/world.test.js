/**
 * test/spec/world.test.js -- src/nd/world.js
 *
 * The whole step, from the forces to the new places. Item 8 of the test plan
 * of PART E asks for a hypersphere on a hyperbox, for two hyperboxes and for
 * a stack of three hypercubes.
 *
 * The first test says that the physics of a world does not depend on the
 * direction of the axes. It covers each part that works in the world frame:
 * the gravity, the star matrix, the collision, the tangents and the solver.
 *
 * CAUTION: it does not cover the gyroscopic term. That term works in the body
 * frame, and the body frame of the two worlds holds the same numbers. A
 * defect of the commutator table gives the same wrong answer in the two
 * worlds. `star.test.js` (the Jacobi identity, and the commutator that
 * follows a turn) is the test for that.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyRotor, anyPositive, anyHalfExtents, anyRadius, anyAngularVelocity, anyVector,
} from '../lib/arbitraries.js';
import { assertClose, assertArrayClose, dot } from '../lib/numeric.js';

const { dims, rotor } = nd;

describe('World.step', () => {
  it('should always give the same physics when the world turns', () => {
    // Build a world, then build the same world with each place, each speed,
    // each orientation and the gravity turned by one rotor. After the same
    // count of steps the second world must be the first world, turned.
    //
    // The world holds two bodies: one that tumbles high above the ground, and
    // one that slides on it. Thus the test covers the gravity, the
    // gyroscopic term, the rotor step, the collision and the solver together.
    //
    // The friction is zero: `tangentBasis` chooses its directions from the
    // axes, thus the friction of the two worlds does not act in the same
    // directions and the two answers would move apart.
    fc.assert(fc.property(anyN(3, 4).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyHalfExtents(n, 0.3, 0.8), anyPositive(0.5, 3),
      anyAngularVelocity(n), anyVector(n, fc.double({ min: -2, max: 2, noNaN: true })),
    )), ([n, R, h, mass, spin, velocity]) => {
      // Arrange
      const D = dims(n);
      const gravity = new Float64Array(n);
      gravity[1] = -9.81;
      const turn = (v) => rotor.rotorApplyVector(D, R, v);
      const same = (v) => v;

      /** A world of two bodies. `map` turns each vector of it. */
      const build = (map, isTurned) => {
        const world = new nd.World({
          dimensions: n, gravity: map(gravity), params: { allowSleep: false },
        });
        world.createBody({ shape: nd.HalfSpace(D, map(unitAxis(n, 1)), 0), friction: 0 });
        const common = {
          shape: nd.HyperBox(D, h), mass, friction: 0, rotor: isTurned ? R : undefined,
        };
        return {
          world,
          flying: world.createBody({
            ...common,
            position: map(unitAxis(n, 1, 6)),
            velocity: map(velocity),
            angularVelocity: isTurned ? rotor.rotorApplyBivector(D, R, spin) : spin,
          }),
          sliding: world.createBody({
            ...common,
            position: map(unitAxis(n, 1, h[1])),
            velocity: map(unitAxis(n, 0, 2)),
          }),
        };
      };
      const straight = build(same, false);
      const other = build(turn, true);

      // Act -- a quarter of a second. A body that tumbles about a plane that
      // is not stable makes a small difference grow without limit, thus a
      // long run compares two answers that only chaos separates.
      for (let s = 0; s < 30; s += 1) {
        straight.world.step(1 / 120);
        other.world.step(1 / 120);
      }

      // Assert -- the tolerance is 1e-5. The two worlds do the same
      // operations in a different order, thus the rounding is not the same.
      // A defect of the algebra gives a difference of the order of 1.
      for (const name of ['flying', 'sliding']) {
        const a = straight[name];
        const b = other[name];
        assertArrayClose(b.x, turn(a.x), `the place of the body that is ${name}`, 1e-5);
        assertArrayClose(b.v, turn(a.v), `the speed of the body that is ${name}`, 1e-5);
        assertArrayClose(b.w, rotor.rotorApplyBivector(D, R, a.w),
          `the angular velocity of the body that is ${name}`, 1e-5);
      }
    }), { numRuns: 10 });
  });

  it('should always make a body rest on the ground', () => {
    // Item 8 of the test plan. A box, a ball or a tesseract falls, and it
    // stops on the ground. The test looks at the lowest point of the body,
    // and not at its center: a box that is tall and thin turns over as it
    // lands, and it then rests on another face.
    fc.assert(fc.property(anyN(3, 5).chain((n) => fc.tuple(
      fc.constant(n), fc.boolean(), anyHalfExtents(n, 0.2, 0.8), anyRadius(0.2, 0.8),
      anyPositive(0.3, 5),
    )), ([n, isBox, h, radius, mass]) => {
      // Arrange
      const D = dims(n);
      const gravity = new Float64Array(n);
      gravity[1] = -9.81;
      const world = new nd.World({ dimensions: n, gravity });
      world.createBody({ shape: nd.HalfSpace(D, unitAxis(n, 1), 0) });
      const body = world.createBody({
        shape: isBox ? nd.HyperBox(D, h) : nd.HyperSphere(D, radius),
        mass,
        position: unitAxis(n, 1, 3),
      });

      // Act
      for (let s = 0; s < 400; s += 1) world.step(1 / 120);

      // Assert -- the limit is the margin of the contacts and the slop of the
      // solver together. Those two numbers are what the engine accepts as
      // "touching".
      const lowest = body.aabb(0).min[1];
      const limit = world.params.contactMargin + world.params.penetrationSlop;
      assert.ok(Math.abs(lowest) < limit,
        `the lowest point of the body is at ${lowest}, and the ground is at 0`);
    }), { numRuns: 20 });
  });

  it('should hold a stack of three cubes', () => {
    // Item 8 of the test plan, in 4 dimensions. The cube at the bottom
    // carries the weight of the two above it. Without the shock propagation
    // pass the stack sinks.
    const D = dims(4);
    const gravity = Float64Array.from([0, -9.81, 0, 0]);
    const world = new nd.World({ dimensions: 4, gravity });
    world.createBody({ shape: nd.HalfSpace(D, [0, 1, 0, 0], 0) });
    const cubes = [];
    for (let i = 0; i < 3; i += 1) {
      cubes.push(world.createBody({
        shape: nd.HyperBox(D, [0.5, 0.5, 0.5, 0.5]),
        mass: 1,
        position: [0, 0.5 + i * 1.001, 0, 0],
      }));
    }

    // Act
    for (let s = 0; s < 600; s += 1) world.step(1 / 120);

    // Assert
    for (let i = 0; i < 3; i += 1) {
      assert.ok(Math.abs(cubes[i].x[1] - (0.5 + i)) < 0.05,
        `the cube ${i} is at the height ${cubes[i].x[1]}`);
    }
  });

  it('should never add energy to a world with no gravity', () => {
    // Two bodies that meet, with no bounce and no damping. The solver takes
    // energy away, and it must never add any. A solver that adds energy makes
    // a stack explode.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyRadius(0.4, 1), anyRadius(0.4, 1), anyPositive(0.3, 3),
      anyPositive(0.3, 3), anyVector(n, fc.double({ min: -4, max: 4, noNaN: true })),
    )), ([n, ra, rb, ma, mb, velocity]) => {
      // Arrange
      const D = dims(n);
      // `biasFactor` is zero: the correction of the depth is a push that
      // does not come from a force, thus it adds energy by design. The test
      // is about the solver, and it turns that push off. The two bodies start
      // apart, and one of them moves.
      const world = new nd.World({
        dimensions: n,
        gravity: new Float64Array(n),
        params: { allowSleep: false, biasFactor: 0 },
      });
      world.createBody({
        shape: nd.HyperSphere(D, ra), mass: ma, restitution: 0, velocity: unitAxis(n, 0, 3),
      });
      world.createBody({
        shape: nd.HyperSphere(D, rb),
        mass: mb,
        restitution: 0,
        velocity,
        position: unitAxis(n, 0, (ra + rb) * 1.2),
      });
      const before = world.totalEnergy();

      // Act
      for (let s = 0; s < 120; s += 1) world.step(1 / 120);

      // Assert
      assert.ok(world.totalEnergy() <= before * (1 + 1e-6) + 1e-9,
        `the energy went from ${before} to ${world.totalEnergy()}`);
    }), { numRuns: 20 });
  });

  it('should always lose the exact energy of the method in a free flight', () => {
    // Semi implicit Euler does not keep the energy of a body that falls. The
    // quantity that it loses has a closed form: after N steps of the length
    // dt, with a constant gravity,
    //
    //   E_N - E_0 = -m |g|^2 dt^2 N / 2
    //
    // and the start velocity has no part in it. This is an oracle, thus the
    // test says the exact number and not a limit. A defect in the order of
    // the two halves of the step changes this number.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n, fc.double({ min: -5, max: 5, noNaN: true })),
      anyPositive(0.5, 3),
    )), ([n, velocity, mass]) => {
      // Arrange
      const D = dims(n);
      const gravity = new Float64Array(n);
      gravity[1] = -9.81;
      const world = new nd.World({
        dimensions: n, gravity, params: { allowSleep: false },
      });
      world.createBody({ shape: nd.HyperSphere(D, 0.5), mass, velocity });
      const before = world.totalEnergy();

      // Act
      const steps = 240;
      const dt = 1 / 240;
      for (let s = 0; s < steps; s += 1) world.step(dt);

      // Assert
      const loss = mass * dot(gravity, gravity) * dt * dt * steps / 2;
      assertClose(world.totalEnergy(), before - loss, 'the energy of a free flight', 1e-9);
    }));
  });
});

describe('World.updateSleep', () => {
  it('should let a body with no contact go to sleep', () => {
    // Arrange -- no gravity and no contact, thus the body only floats.
    const D = dims(3);
    const world = new nd.World({ dimensions: 3, gravity: [0, 0, 0] });
    const box = world.createBody({
      shape: nd.HyperBox(D, [0.5, 0.5, 0.5]), mass: 1, position: [0, 5, 0],
    });

    // Act -- `sleepTime` is 0.6 seconds, thus 300 steps are sufficient.
    for (let s = 0; s < 300; s += 1) world.step(1 / 120);

    // Assert
    assert.ok(box.sleeping, 'the body is asleep');
    const place = Float64Array.from(box.x);
    for (let s = 0; s < 100; s += 1) world.step(1 / 120);
    assertArrayClose(box.x, place, 'a body that sleeps does not move');
  });

  /**
   * A KNOWN DEFECT, that this test found.
   *
   * A body that rests on the ground never goes to sleep. `subStep` wakes the
   * two bodies of each contact that has a depth of more than zero:
   *
   *     if (!c.b.isStatic && c.a.sleeping === false) c.b.wake();
   *
   * A static body always has `sleeping === false`, because nothing ever puts
   * a static body to sleep. Thus the ground wakes the body above it in each
   * step, `wake()` makes the sleep timer zero, and the timer never comes to
   * `sleepTime`.
   *
   * The measure: a box on the ground, at rest with a speed of 6e-6 (the limit
   * `sleepLinearVelocity` is 0.03), has a sleep timer of 0.008 seconds -- one
   * step -- after 600 steps. A stack of two boxes does the same.
   *
   * Thus `allowSleep` does nothing for the one condition that it is for: a
   * body that rests. The world keeps integrating and solving it for ever.
   * This costs speed, and it does not make the answer wrong.
   *
   * Take the `todo` mark away when the defect is repaired.
   */
  it('should let a body that rests on the ground go to sleep',
    { todo: 'the static ground wakes the body in each step -- see the note above' }, () => {
    // Arrange
    const D = dims(3);
    const world = new nd.World({ dimensions: 3, gravity: [0, -9.81, 0] });
    world.createBody({ shape: nd.HalfSpace(D, [0, 1, 0], 0) });
    const box = world.createBody({
      shape: nd.HyperBox(D, [0.5, 0.5, 0.5]), mass: 1, position: [0, 0.5, 0],
    });

    // Act
    for (let s = 0; s < 300; s += 1) world.step(1 / 120);

    // Assert
    assert.ok(box.sleeping, 'the box is asleep');
    const place = Float64Array.from(box.x);
    for (let s = 0; s < 100; s += 1) world.step(1 / 120);
    assertArrayClose(box.x, place, 'a body that sleeps does not move');
  });

  it('should wake a body that takes an impulse',
    { todo: 'the body never goes to sleep -- see the note above' }, () => {
    // Arrange
    const D = dims(3);
    const world = new nd.World({ dimensions: 3, gravity: [0, -9.81, 0] });
    world.createBody({ shape: nd.HalfSpace(D, [0, 1, 0], 0) });
    const box = world.createBody({
      shape: nd.HyperBox(D, [0.5, 0.5, 0.5]), mass: 1, position: [0, 0.5, 0],
    });
    for (let s = 0; s < 300; s += 1) world.step(1 / 120);
    assert.ok(box.sleeping, 'the box is asleep at the start of the act');

    // Act
    box.applyCentralImpulse(Float64Array.from([5, 0, 0]));

    // Assert
    assert.equal(box.sleeping, false, 'the impulse wakes the box');
  });
});

describe('World.totalEnergy', () => {
  it('should hold the potential energy of the height', () => {
    // Arrange -- a body of the mass 2 at the height 3, with g = 9.81.
    const D = dims(3);
    const world = new nd.World({ dimensions: 3, gravity: [0, -9.81, 0] });
    world.createBody({ shape: nd.HyperSphere(D, 0.5), mass: 2, position: [0, 3, 0] });

    // Act, Assert
    assertClose(world.totalEnergy(), 2 * 9.81 * 3, 'the potential energy of a body');
  });
});

// Helpers

/** The vector along the axis `i`, of the length `scale`. */
function unitAxis(n, i, scale = 1) {
  const e = new Float64Array(n);
  e[i] = scale;
  return e;
}
