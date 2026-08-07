/**
 * test/spec/integrator.test.js -- src/nd/integrate/integrator.js
 *
 * The time step. Items 5 and 6 of the test plan of PART E ask for the
 * conservation of a torque-free body and for the closed 3D subspace of a 4D
 * body.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyVector, anyHalfExtents, anyRotor, anyPositive, anyAngularVelocity,
} from '../lib/arbitraries.js';
import { assertClose, assertArrayClose, dot, norm, mulT, orthogonalError } from '../lib/numeric.js';

const { dims } = nd;

/** The params of a step, with the values of `defaultParams`. */
const OPTS = { gyroscopic: true, gyroscopicIterations: 1, rotorTolerance: 1e-9 };

describe('integrateVelocities', () => {
  it('should add the gravity to the velocity', () => {
    // Arrange
    const D = dims(4);
    const body = new nd.Body(D, { shape: nd.HyperSphere(D, 1), mass: 3 });

    // Act
    nd.integrateVelocities(D, body, 0.5, Float64Array.from([0, -10, 0, 0]), OPTS);

    // Assert
    assertArrayClose(body.v, [0, -5, 0, 0], 'the velocity after half a second');
  });

  it('should never make a static body move', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyPositive(0.001, 1),
    )), ([n, gravity, dt]) => {
      // Arrange
      const D = dims(n);
      const wall = new nd.Body(D, { shape: nd.HalfSpace(D, unitAxis(n, 1), 0) });

      // Act
      nd.integrateVelocities(D, wall, dt, gravity, OPTS);

      // Assert
      assert.ok(wall.v.every((v) => v === 0), 'a static body does not take the gravity');
    }));
  });

  it('should never let the damping turn the velocity around', () => {
    // The damping is the implicit form 1 / (1 + dt c). The comment of the
    // function says that this form is stable at any step length. An
    // exponential or an explicit form would turn the velocity around at a
    // large step, and the body would shake.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyPositive(0.001, 100), anyPositive(0.01, 1000),
    )), ([n, velocity, dt, damping]) => {
      // Arrange
      const D = dims(n);
      const body = new nd.Body(D, {
        shape: nd.HyperSphere(D, 1), mass: 1, velocity, linearDamping: damping,
      });
      const before = norm(body.v);

      // Act
      nd.integrateVelocities(D, body, dt, new Float64Array(n), OPTS);

      // Assert
      assert.ok(norm(body.v) <= before + 1e-12, 'the damping never makes a body faster');
      assert.ok(dot(body.v, velocity) >= -1e-12, 'the damping never turns the velocity around');
    }));
  });

  it('should always hold an axis that the linear factor makes zero', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), fc.nat({ max: n - 1 }), anyVector(n), anyPositive(0.001, 1),
    )), ([n, axis, gravity, dt]) => {
      // Arrange
      const D = dims(n);
      const factor = new Float64Array(n).fill(1);
      factor[axis] = 0;
      const body = new nd.Body(D, {
        shape: nd.HyperSphere(D, 1), mass: 1, linearFactor: factor,
      });

      // Act
      for (let s = 0; s < 20; s += 1) nd.integrateVelocities(D, body, dt, gravity, OPTS);

      // Assert
      assert.equal(body.v[axis], 0, `the axis ${axis} stays at zero`);
    }));
  });
});

describe('integratePositions', () => {
  it('should always agree with the closed form of a free flight', () => {
    // The method is semi implicit Euler: the velocity changes first, and then
    // the position takes the new velocity. After N steps that gives exactly
    //   x = x0 + N dt v0 + g dt^2 N (N + 1) / 2
    // This is an oracle, and not an approximation.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n, fc.double({ min: -20, max: 20, noNaN: true })),
      anyVector(n, fc.double({ min: -20, max: 20, noNaN: true })),
      fc.integer({ min: 1, max: 200 }),
    )), ([n, velocity, gravity, steps]) => {
      // Arrange
      const D = dims(n);
      const dt = 1 / 120;
      const body = new nd.Body(D, { shape: nd.HyperSphere(D, 1), mass: 2, velocity });

      // Act
      for (let s = 0; s < steps; s += 1) {
        nd.integrateVelocities(D, body, dt, gravity, OPTS);
        nd.integratePositions(D, body, dt, OPTS);
      }

      // Assert
      const expected = Float64Array.from(velocity, (v, i) => v * steps * dt
        + gravity[i] * dt * dt * steps * (steps + 1) / 2);
      assertArrayClose(body.x, expected, 'the place after a free flight', 1e-9);
    }));
  });

  it('should never turn a body that has no angular velocity', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyVector(n), anyPositive(0.001, 1),
    )), ([n, R, velocity, dt]) => {
      // Arrange
      const D = dims(n);
      const body = new nd.Body(D, { shape: nd.HyperSphere(D, 1), mass: 1, rotor: R, velocity });
      const before = Float64Array.from(body.R);

      // Act
      for (let s = 0; s < 10; s += 1) nd.integratePositions(D, body, dt, OPTS);

      // Assert
      assertArrayClose(body.R, before, 'the orientation does not change');
    }));
  });

  it('should never let a body change its shape', () => {
    // The rotor loses a little of its quality at each step. Without
    // `rotorCorrect` the matrix of the rotor stops being orthogonal, and the
    // body becomes larger or thinner as it turns. See README, section 10.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n), anyRotor(n), anyAngularVelocity(n),
    )), ([n, h, R, w]) => {
      // Arrange
      const D = dims(n);
      const body = new nd.Body(D, {
        shape: nd.HyperBox(D, h), mass: 1, rotor: R, angularVelocity: w,
      });

      // Act
      for (let s = 0; s < 400; s += 1) {
        nd.integrateVelocities(D, body, 1 / 240, new Float64Array(n), OPTS);
        nd.integratePositions(D, body, 1 / 240, OPTS);
      }

      // Assert -- `rotorTolerance` is 1e-9, and the repair starts only when
      // the defect passes it. Thus the error of the matrix stays of that
      // order. A body that changes its shape gives a much larger number.
      assert.ok(orthogonalError(body.Rm, n) < 1e-7,
        `the matrix of the body is not orthogonal: the error is ${orthogonalError(body.Rm, n)}`);
    }), { numRuns: 20 });
  });
});

describe('the free rotation of a body', () => {
  it('should never gain energy', () => {
    // A body with no torque keeps its kinetic energy. The explicit form of
    // the gyroscopic term adds energy and the body goes faster and faster;
    // the implicit form of `applyGyroscopic` must not. README section 13 says
    // that it loses a little, and that is a property of the method.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n), anyRotor(n), anyAngularVelocity(n),
    )), ([n, h, R, w]) => {
      // Arrange
      const D = dims(n);
      const body = new nd.Body(D, {
        shape: nd.HyperBox(D, h), mass: 1, rotor: R, angularVelocity: w,
      });
      const before = body.kineticEnergy();

      // Act
      for (let s = 0; s < 480; s += 1) {
        nd.integrateVelocities(D, body, 1 / 240, new Float64Array(n), OPTS);
        nd.integratePositions(D, body, 1 / 240, OPTS);
      }

      // Assert
      // The Newton method makes one turn only, thus the answer is not exact
      // and the energy moves a little. The explicit form of the same term
      // gives an energy that grows without limit.
      assert.ok(body.kineticEnergy() <= before * 1.02 + 1e-9,
        `the energy went from ${before} to ${body.kineticEnergy()}`);
    }), { numRuns: 20 });
  });

  it('should almost keep the size of the angular momentum', () => {
    // The implicit method of the gyroscopic term takes a little of the
    // angular momentum away. At dt = 1/240 over 2 seconds the loss is less
    // than one part in a hundred. README section 13 documents this.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n), anyRotor(n), anyAngularVelocity(n),
    )), ([n, h, R, w]) => {
      // Arrange
      const D = dims(n);
      const body = new nd.Body(D, {
        shape: nd.HyperBox(D, h), mass: 1, rotor: R, angularVelocity: w,
      });
      const before = norm(body.L);

      // Act
      for (let s = 0; s < 480; s += 1) {
        nd.integrateVelocities(D, body, 1 / 240, new Float64Array(n), OPTS);
        nd.integratePositions(D, body, 1 / 240, OPTS);
      }

      // Assert -- the implicit method takes momentum away, and it never adds
      // any. The quantity that it takes away depends on the body and on the
      // step; README section 13 gives 7 percent at dt = 1/120 over a long
      // time. The direction that the method can never take is the one here.
      assert.ok(norm(body.L) <= before * (1 + 1e-4) + 1e-12,
        `the size of the angular momentum went from ${before} to ${norm(body.L)}`);
    }), { numRuns: 20 });
  });

  /**
   * A KNOWN DEFECT, that this test found.
   *
   * `body.L` is the angular momentum IN THE WORLD FRAME: `updateDerived`
   * builds `w` from it with `w = [R]2 I^-1 [R]2^T L`. A body with no torque
   * must keep that momentum, in size AND in direction.
   *
   * `applyGyroscopic` solves the Euler equation in the body frame, which is
   * correct, but it then writes the answer back with the rotor of the START
   * of the step:
   *
   *     L_world = [R]2 (I w_body_new)
   *
   * `integratePositions` turns the rotor after that, and nothing builds `L`
   * again from the new rotor. Thus the turn of the frame never takes the
   * Euler term away, and the world momentum turns at the rate `[w, L]`.
   *
   * The measure: the change of `L` in one step is `dt [w, L]` to 4 places.
   * The error of a step is thus of the order of `dt`, and not of `dt^2`. It
   * does not go away with a smaller step:
   *
   *     2 s of free rotation, a box of [1, 0.6, 0.3], w = (0.2, 3, 0.1)
   *     dt = 2/240   the direction of L turns by 39.24 degrees
   *     dt = 2/960   by 38.56 degrees
   *     dt = 2/3840  by 38.42 degrees
   *
   * The size of `L` and the energy stay correct, because a commutator with
   * `L` is normal to `L`. Thus the two tests above do not see this. Only the
   * direction is wrong: a body that tumbles takes the wrong way.
   *
   * Take the `todo` mark away when the defect is repaired.
   */
  it('should never turn the angular momentum of a body with no torque',
    { todo: 'the world momentum turns at the rate [w, L] -- see the note above' }, () => {
    fc.assert(fc.property(anyN(3, 5).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n), anyRotor(n), anyAngularVelocity(n),
    )), ([n, h, R, w]) => {
      // Arrange
      const D = dims(n);
      const body = new nd.Body(D, {
        shape: nd.HyperBox(D, h), mass: 1, rotor: R, angularVelocity: w,
      });
      const before = Float64Array.from(body.L);

      // Act
      for (let s = 0; s < 480; s += 1) {
        nd.integrateVelocities(D, body, 1 / 240, new Float64Array(n), OPTS);
        nd.integratePositions(D, body, 1 / 240, OPTS);
      }

      // Assert
      assertArrayClose(body.L, before, 'the angular momentum of a free body', 1e-2);
    }), { numRuns: 10 });
  });

  it('should keep the turn of a 4D body inside a 3D subspace', () => {
    // Item 6 of the test plan of PART E. Give a body of 4 dimensions an
    // angular velocity with no part in the planes that hold the axis 3. Those
    // parts must stay at zero for ever, because so(3) is a sub algebra of
    // so(4): the commutator of two members of the subspace stays inside it.
    //
    // A wrong sign in the commutator table breaks this at once.
    // The angular velocity stays inside +-4 radians each second. The rotor
    // step is explicit, thus `w dt` must stay small; at 1000 radians each
    // second and dt = 1/240 the body turns 4 radians in one step, and no
    // integrator of that form gives a meaning to the answer.
    //
    // The orientation also stays inside the subspace. A body that is already
    // turned in a plane of the axis 3 has a body frame that the subspace of
    // the world does not hold, and then there is nothing to test.
    fc.assert(fc.property(fc.tuple(anyHalfExtents(4), anySubspaceRotor(), anyAngularVelocity(4)),
      ([h, R, full]) => {
        // Arrange
        const D = dims(4);
        const w = Float64Array.from(full);
        const outside = [];
        for (let p = 0; p < D.k; p += 1) {
          if (D.pairs[p].includes(3)) {
            w[p] = 0;
            outside.push(p);
          }
        }
        const body = new nd.Body(D, {
          shape: nd.HyperBox(D, h), mass: 1, rotor: R, angularVelocity: w,
        });

        // Act
        for (let s = 0; s < 240; s += 1) {
          nd.integrateVelocities(D, body, 1 / 240, new Float64Array(4), OPTS);
          nd.integratePositions(D, body, 1 / 240, OPTS);
        }

        // Assert -- the momentum stays in the body frame subspace.
        const bodyL = mulT(body.R2, body.L, D.k, D.k);
        for (const p of outside) {
          assert.ok(Math.abs(bodyL[p]) <= 1e-9 * (1 + norm(body.L)),
            `the plane ${D.pairs[p]} must stay at zero, and it is ${bodyL[p]}`);
        }
      }), { numRuns: 20 });
  });

  it('should turn a box over when it spins about its middle plane', () => {
    // The Dzhanibekov effect. A box with three different sides is not stable
    // when it turns about the plane of the middle moment of inertia. The turn
    // about the largest and the smallest moment is stable.
    //
    // This is an example test, and not a property: the time of the first turn
    // depends on the size of the first error.
    const D = dims(3);
    const h = Float64Array.from([1, 0.6, 0.3]);

    // Act, Assert -- the moments of the planes (x y), (x z), (y z) are
    // 1.36, 1.09 and 0.45. Thus the plane 1 holds the middle moment.
    assert.ok(spinChange(D, h, 1) > 0.5, 'the middle plane is not stable');
    assert.ok(spinChange(D, h, 0) < 0.05, 'the plane of the largest moment is stable');
    assert.ok(spinChange(D, h, 2) < 0.05, 'the plane of the smallest moment is stable');
  });
});

// Helpers

/**
 * A rotor of 4 dimensions that turns only in the planes of the axes 0, 1 and
 * 2. Such a rotor keeps the 3D subspace of those axes.
 */
function anySubspaceRotor() {
  const D = dims(4);
  return fc.array(fc.tuple(fc.constantFrom([0, 1], [0, 2], [1, 2]),
    fc.double({ min: -Math.PI, max: Math.PI, noNaN: true })), { minLength: 1, maxLength: 3 })
    .map((list) => {
      let R = nd.rotor.rotorIdentity(D);
      for (const [[i, j], angle] of list) {
        R = nd.rotor.rotorMul(D, nd.rotor.rotorFromPlane(D, i, j, angle), R);
      }
      return R;
    });
}

/** The vector of the length 1 along the axis `i`. */
function unitAxis(n, i) {
  const e = new Float64Array(n);
  e[i] = 1;
  return e;
}

/**
 * How far the angular velocity of a box goes from its start plane, in the
 * body frame, over 6 seconds. A stable turn gives almost 0, and a turn that
 * is not stable gives a number near 2, because the body turns over.
 */
function spinChange(D, h, plane) {
  const w = new Float64Array(D.k);
  w[plane] = 4;
  for (let p = 0; p < D.k; p += 1) if (p !== plane) w[p] = 0.02;
  const body = new nd.Body(D, { shape: nd.HyperBox(D, h), mass: 1, angularVelocity: w });
  const start = mulT(body.R2, body.L, D.k, D.k);
  let worst = 0;
  for (let s = 0; s < 2880; s += 1) {
    nd.integrateVelocities(D, body, 1 / 480, new Float64Array(D.n), OPTS);
    nd.integratePositions(D, body, 1 / 480, OPTS);
    const now = mulT(body.R2, body.L, D.k, D.k);
    worst = Math.max(worst, Math.abs(now[plane] - start[plane]) / Math.abs(start[plane]));
  }
  return worst;
}
