/**
 * test/spec/constraint.test.js -- src/nd/resolve/constraint.js
 *
 * The constraints and the joints. See ND-PHYSICS.md, A13 and B8.
 *
 * The laws here:
 *   - The impulse of a row is internal, thus it keeps the momentum of the
 *     pair. This is the strongest test of the Jacobian of a row.
 *   - The error of the rotation carries the same sign as the angular
 *     velocity. The sign is the one risk of the bivector form, thus two tests
 *     hold it: one on the error alone, and one on the correction.
 *   - A joint drives its error to zero, and it never adds energy.
 *   - A row that has no effective mass goes away, and it gives no NaN.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd, loadPhysiNFor } from '../lib/load.js';
import { THREE } from '../lib/three-stub.js';
import { anyN, anyRotor, anyPositive, anyVector } from '../lib/arbitraries.js';
import { assertClose, maxAbs, norm, dot } from '../lib/numeric.js';

const { dims, World, Body, HyperSphere, HyperBox, defaultParams } = nd;
const { rotorFromPlane, rotorMul, rotorMatrix, rotorDefect } = nd.rotor;

/** The 3 dimensional tables, for the tests that do not sweep `n`. */
const D3 = dims(3);

/** A ball of the radius 0.5, at `position`. */
function ball(D, mass, position) {
  return new Body(D, { shape: HyperSphere(D, 0.5), mass, position });
}

/** A world with no gravity and no sleep. A joint test wants a still world. */
function quietWorld(n, params = {}) {
  return new World({
    dimensions: n,
    gravity: new Float64Array(n),
    params: Object.assign({ allowSleep: false }, params),
  });
}

/** The total linear momentum of the bodies. */
function linearMomentum(world) {
  const { n } = world.D;
  const p = new Float64Array(n);
  for (const b of world.bodies) {
    for (let i = 0; i < n; i += 1) p[i] += b.mass * b.v[i];
  }
  return p;
}

/**
 * The total angular momentum about the origin: the spin `L` of each body plus
 * the moment `x ^ (m v)` of its linear momentum. A pair of impulses that are
 * opposite and at the same world point keeps this.
 */
function angularMomentum(world) {
  const { D } = world;
  const { n, k } = D;
  const total = new Float64Array(k);
  const mv = new Float64Array(n);
  for (const b of world.bodies) {
    for (let i = 0; i < n; i += 1) mv[i] = b.mass * b.v[i];
    const moment = nd.wedgeVec(D, b.x, mv);
    for (let p = 0; p < k; p += 1) total[p] += b.L[p] + moment[p];
  }
  return total;
}

describe('angularError', () => {
  it('should give +sin(angle) on the plane that turned, and zero on the others', () => {
    // The sign of the error of the rotation is the one real risk of the
    // bivector form, because `dR/dt = -(1/2) w R` is not the usual sign. This
    // test pins it: a turn from the axis `i` toward the axis `j` gives a
    // PLUS on that plane. See `angularError` and ND-PHYSICS.md, A13.
    for (let n = 2; n <= 5; n += 1) {
      const D = dims(n);
      for (let p = 0; p < D.k; p += 1) {
        const [i, j] = D.pairs[p];

        // Arrange: `a` at the identity, and `b` at the same place. The rest
        // relation of the joint is then the relation of no turn.
        const a = ball(D, 1, new Float64Array(n));
        const b = ball(D, 1, new Float64Array(n));
        const joint = nd.createConstraint(D, { type: 'fixed' }, a, b, defaultParams);

        // Act: turn `b` by 0.1 in the plane `p`.
        b.R.set(rotorFromPlane(D, i, j, 0.1));
        b.updateDerived();
        const e = nd.angularError(D, a, b, joint.restRel);

        // Assert
        assertClose(e[p], Math.sin(0.1),
          `n = ${n}, the plane ${p} = (${i} ${j})`, 1e-12);
        for (let q = 0; q < D.k; q += 1) {
          if (q !== p) assert.ok(Math.abs(e[q]) < 1e-12, `n = ${n}, the plane ${q} leaked`);
        }
      }
    }
  });

  it('should be exact for a simple rotation of any angle', () => {
    // A rotation in ONE plane gives exactly `sin(angle)` on that plane, at any
    // angle. That is why a hinge, whose free part is one plane, can use this
    // error although the hinge turns a long way.
    const D = dims(4);
    const a = ball(D, 1, [0, 0, 0, 0]);
    const b = ball(D, 1, [0, 0, 0, 0]);
    const joint = nd.createConstraint(D, { type: 'fixed' }, a, b, defaultParams);
    for (const angle of [0.3, 1.2, 2.5]) {
      b.R.set(rotorFromPlane(D, 0, 1, angle));
      b.updateDerived();
      const e = nd.angularError(D, a, b, joint.restRel);
      assertClose(e[0], Math.sin(angle), `the angle ${angle}`, 1e-12);
      assert.ok(maxAbs(e.slice(1)) === 1, 'no other plane holds an error');
    }
  });

  it('should mix the planes for a rotation in two planes together', () => {
    // This is the limit of the first order logarithm, and the reason that
    // `prepareRows` only takes the error when the joint holds `k` or `k - 1`
    // planes. A turn of 1.0 in (x y) and of 1.0 in (y z) shows an error in
    // (x z), in which nothing turned.
    const D = dims(4);
    const a = ball(D, 1, [0, 0, 0, 0]);
    const b = ball(D, 1, [0, 0, 0, 0]);
    const joint = nd.createConstraint(D, { type: 'fixed' }, a, b, defaultParams);
    b.R.set(rotorMul(D, rotorFromPlane(D, 1, 2, 1), rotorFromPlane(D, 0, 1, 1)));
    b.updateDerived();
    const e = nd.angularError(D, a, b, joint.restRel);
    assert.ok(Math.abs(e[1]) > 0.3,
      `the plane (x z) must show the mixing, and it holds ${e[1]}`);
  });

  it('should be zero when the joint starts, at any orientation', () => {
    // `createConstraint` takes the rest relation from the two bodies as they
    // are. Thus a weld of two bodies that already turned starts with no error
    // and it gives no jump.
    fc.assert(fc.property(anyN(2, 5), (n) => {
      const D = dims(n);
      return fc.assert(fc.property(anyRotor(n), anyRotor(n), (ra, rb) => {
        const a = ball(D, 1, new Float64Array(n));
        const b = ball(D, 1, new Float64Array(n));
        a.R.set(ra); a.updateDerived();
        b.R.set(rb); b.updateDerived();
        const joint = nd.createConstraint(D, { type: 'fixed' }, a, b, defaultParams);
        const e = nd.angularError(D, a, b, joint.restRel);
        assert.ok(maxAbs(e) === 1 && norm(e) < 1e-9,
          `n = ${n}, the error at the start is ${norm(e)}`);
      }), { numRuns: 20 });
    }), { numRuns: 4 });
  });
});

describe('solveJoint', () => {
  it('should turn the body back the other way, thus the correction has the right sign', () => {
    // The test above holds the sign of the ERROR. This one holds the sign of
    // the BIAS: a positive error must give a negative angular velocity. A sign
    // that is wrong only in the bias would pass the first test and fail here,
    // and the joint would then run away.
    for (let n = 2; n <= 5; n += 1) {
      const D = dims(n);
      for (let p = 0; p < D.k; p += 1) {
        const [i, j] = D.pairs[p];

        // Arrange: `a` is static, `b` is free and still, and `b` holds a
        // positive error of the rotation on the plane `p`.
        const a = new Body(D, { shape: HyperSphere(D, 0.5), mass: 0 });
        const b = ball(D, 1, new Float64Array(n));
        const params = Object.assign({}, defaultParams, { constraintBias: 1 });
        const joint = nd.createConstraint(D, { type: 'fixed' }, a, b, params);
        b.R.set(rotorFromPlane(D, i, j, 0.1));
        b.updateDerived();

        // Act
        nd.prepareRows(D, joint, 1 / 60, params);
        nd.solveJoint(D, joint);

        // Assert
        assert.ok(b.w[p] < -1e-6,
          `n = ${n}, the plane ${p}: the correction is ${b.w[p]} and it must be negative`);
        for (let q = 0; q < D.k; q += 1) {
          if (q !== p) {
            assert.ok(Math.abs(b.w[q]) < 1e-9, `n = ${n}, the plane ${q} must not move`);
          }
        }
      }
    }
  });

  it('should keep the momentum of the pair, because an impulse of a row is internal', () => {
    // A row gives `+j` to `b` at `rB` and `-j` to `a` at `rA`, and the two
    // points are the same point. Thus the linear momentum and the angular
    // momentum about the origin do not change. This is the strongest test of
    // the Jacobian of a row: a wrong lever arm breaks the angular part.
    fc.assert(fc.property(
      anyN(2, 5), anyPositive(0.5, 5), anyPositive(0.5, 5),
      (n, ma, mb) => {
        const D = dims(n);
        return fc.assert(fc.property(
          anyRotor(n), anyRotor(n), anyVector(n),
          (ra, rb, la) => {
            // Arrange
            const world = quietWorld(n);
            const a = world.addBody(ball(D, ma, new Float64Array(n)));
            const b = world.addBody(ball(D, mb, new Float64Array(n).fill(1)));
            a.R.set(ra); a.updateDerived();
            b.R.set(rb); b.updateDerived();
            a.setLinearVelocity(new Float64Array(n).fill(0.3));
            b.setAngularVelocity(new Float64Array(D.k).fill(0.2));
            // The two anchors must be the SAME world point. A pair of opposite
            // impulses at two different points is a couple, and a couple does
            // change the angular momentum. The law holds for a joint that is
            // closed, which is the state that the solver drives toward.
            const world0 = a.localToWorld(la, new Float64Array(n));
            const lbSame = b.worldToLocal(world0, new Float64Array(n));
            const joint = world.createConstraint(
              { type: 'point', localA: la, localB: lbSame }, a, b,
            );
            const p0 = linearMomentum(world);
            const l0 = angularMomentum(world);

            // Act
            nd.prepareRows(D, joint, 1 / 60, world.params);
            for (let t = 0; t < 5; t += 1) nd.solveJoint(D, joint);

            // Assert
            const p1 = linearMomentum(world);
            const l1 = angularMomentum(world);
            for (let i = 0; i < n; i += 1) {
              assertClose(p1[i], p0[i], `n = ${n}, the linear momentum on the axis ${i}`, 1e-9);
            }
            for (let q = 0; q < D.k; q += 1) {
              assertClose(l1[q], l0[q], `n = ${n}, the angular momentum on the plane ${q}`, 1e-9);
            }
          },
        ), { numRuns: 15 });
      },
    ), { numRuns: 4 });
  });

  it('should never add energy, and it should lose less as the step gets smaller', () => {
    // A joint holds a body, thus it does no work in the true equations. The
    // solver works at the velocity level, and it holds the lever arms still
    // through a step. That gives an error of the first order in `dt`: the
    // energy falls, and it falls half as much when `dt` is half as long.
    //
    // The first step is not part of this. A joint that starts on a state that
    // breaks it must take the two points to the same velocity, and that is an
    // impact that loses energy of its own.
    const D = dims(3);
    function drift(dt) {
      const world = quietWorld(3, { constraintBias: 0, useShockPropagation: false });
      const a = world.addBody(new Body(D, {
        shape: HyperBox(D, [0.5, 0.5, 0.5]), mass: 1, position: [-1, 0, 0],
      }));
      const b = world.addBody(new Body(D, {
        shape: HyperBox(D, [0.5, 0.5, 0.5]), mass: 3, position: [1, 0, 0],
      }));
      world.createConstraint({ type: 'point', localA: [1, 0, 0], localB: [-1, 0, 0] }, a, b);
      a.setLinearVelocity([0, 2, 0]);
      b.setAngularVelocity([0.5, 0, 0]);
      world.step(dt);
      const start = world.totalEnergy();
      for (let s = 0; s < Math.round(10 / dt); s += 1) world.step(dt);
      return { start, end: world.totalEnergy() };
    }

    const coarse = drift(1 / 240);
    const fine = drift(1 / 480);
    assert.ok(coarse.end <= coarse.start + 1e-9, 'a joint must never add energy');
    assert.ok(fine.end <= fine.start + 1e-9, 'a joint must never add energy');
    const dc = Math.abs(coarse.end - coarse.start) / coarse.start;
    const df = Math.abs(fine.end - fine.start) / fine.start;
    assert.ok(dc < 0.1, `the drift over 10 seconds is ${dc} and it must be small`);
    assert.ok(df < 0.6 * dc,
      `the drift must fall with the step: ${df} against ${dc}`);
  });
});

describe('the point joint', () => {
  it('should drive the error to zero and hold it there', () => {
    // A pendulum: the bob hangs 2 units from an anchor of the world. The
    // length must stay 2, and the error must not grow with the time.
    for (let n = 2; n <= 5; n += 1) {
      const D = dims(n);
      const gravity = new Float64Array(n);
      gravity[1] = -9.81;
      const world = new World({ dimensions: n, gravity, params: { allowSleep: false } });
      const start = new Float64Array(n);
      start[0] = 2; start[1] = 5;
      const anchorWorld = new Float64Array(n);
      anchorWorld[1] = 5;
      const localA = new Float64Array(n);
      localA[0] = -2;
      const bob = world.addBody(ball(D, 1, start));
      world.createConstraint({ type: 'point', localA, localB: anchorWorld }, bob, null);

      // The anchor is in the frame of the body, and a point joint leaves the
      // rotation free. Thus the world anchor must come from `localToWorld`,
      // and not from `x + localA`.
      const err = () => {
        const p = bob.localToWorld(localA, new Float64Array(n));
        let e = 0;
        for (let i = 0; i < n; i += 1) e += (p[i] - anchorWorld[i]) ** 2;
        return Math.sqrt(e);
      };
      // The error rises and falls with the swing, because the pull is largest
      // at the bottom. Thus take the LARGEST error of each half of the run.
      // Two single readings would only compare two points of that wave.
      const worst = (steps) => {
        let m = 0;
        for (let s = 0; s < steps; s += 1) {
          world.step(1 / 120);
          m = Math.max(m, err());
        }
        return m;
      };
      const first = worst(300);
      const second = worst(300);
      assert.ok(second < 5e-3, `n = ${n}: the error is ${second}`);
      assert.ok(second <= first * 1.5,
        `n = ${n}: the error grew with the time, ${first} to ${second}`);
    }
  });
});

describe('the distance joint', () => {
  it('should hold the length and keep the momentum', () => {
    const D = dims(3);
    const world = quietWorld(3);
    const a = world.addBody(ball(D, 1, [-1, 0, 0]));
    const b = world.addBody(ball(D, 2, [1, 0, 0]));
    a.setLinearVelocity([0, 3, 0]);
    b.setLinearVelocity([0, -1, 0]);
    world.createConstraint({ type: 'distance', rest: 2 }, a, b);
    const p0 = linearMomentum(world);

    for (let s = 0; s < 400; s += 1) world.step(1 / 240);

    const len = Math.hypot(b.x[0] - a.x[0], b.x[1] - a.x[1], b.x[2] - a.x[2]);
    assertClose(len, 2, 'the length', 1e-2);
    const p1 = linearMomentum(world);
    for (let i = 0; i < 3; i += 1) assertClose(p1[i], p0[i], `the momentum on the axis ${i}`, 1e-9);
  });

  it('should take the rest length from the two anchors when you give none', () => {
    const D = dims(3);
    const world = quietWorld(3);
    const a = world.addBody(ball(D, 1, [0, 0, 0]));
    const b = world.addBody(ball(D, 1, [3, 4, 0]));
    const joint = world.createConstraint({ type: 'distance' }, a, b);
    assertClose(joint.rest, 5, 'the rest length', 1e-12);
  });

  it('should let a rope go slack, and only stop it from growing', () => {
    const D = dims(3);
    const world = quietWorld(3);
    const a = world.addBody(new Body(D, { shape: HyperSphere(D, 0.5), mass: 0 }));
    const b = world.addBody(ball(D, 1, [3, 0, 0]));
    world.createConstraint({ type: 'distance', rest: 3, mode: 'rope' }, a, b);

    // Act: push `b` toward `a`. A rope must not stop it.
    b.setLinearVelocity([-1, 0, 0]);
    for (let s = 0; s < 120; s += 1) world.step(1 / 120);
    assert.ok(b.x[0] < 2.2, `the rope went slack, and x is ${b.x[0]}`);

    // Act: push `b` away. The rope must stop it at the rest length.
    b.setLinearVelocity([5, 0, 0]);
    for (let s = 0; s < 240; s += 1) world.step(1 / 120);
    assert.ok(b.x[0] < 3.05, `the rope must hold at 3, and x is ${b.x[0]}`);
  });
});

describe('the fixed joint', () => {
  it('should weld the body, in 3 and in 4 dimensions', () => {
    for (const n of [3, 4]) {
      const D = dims(n);
      const gravity = new Float64Array(n);
      gravity[1] = -9.81;
      const world = new World({ dimensions: n, gravity, params: { allowSleep: false } });
      const anchorPos = new Float64Array(n);
      anchorPos[1] = 3;
      const bodyPos = new Float64Array(n);
      bodyPos[0] = 1; bodyPos[1] = 3;
      const anchor = world.addBody(new Body(D, { shape: HyperSphere(D, 0.5), mass: 0, position: anchorPos }));
      const b = world.addBody(new Body(D, {
        shape: HyperBox(D, new Array(n).fill(0.5)), mass: 1, position: bodyPos,
      }));
      b.R.set(rotorFromPlane(D, 0, n - 1, 0.4));
      b.updateDerived();
      b.setAngularVelocity(new Float64Array(D.k).fill(0.3));
      const localA = new Float64Array(n);
      localA[0] = 1;
      const joint = world.createConstraint(
        { type: 'fixed', localA, localB: new Float64Array(n) }, anchor, b,
      );

      for (let s = 0; s < 400; s += 1) world.step(1 / 120);

      const e = nd.angularError(D, joint.a, joint.b, joint.restRel);
      assert.ok(norm(e) < 1e-2, `n = ${n}: the angular error is ${norm(e)}`);
      let pos = 0;
      for (let i = 0; i < n; i += 1) pos += (b.x[i] - bodyPos[i]) ** 2;
      assert.ok(Math.sqrt(pos) < 5e-3, `n = ${n}: the position error is ${Math.sqrt(pos)}`);
      assert.ok(rotorDefect(D, b.R) < 1e-6, `n = ${n}: the rotor went bad`);
    }
  });
});

describe('the hinge', () => {
  it('should leave one plane free and hold every other plane, at k = 3 and at k = 6', () => {
    // This is the test that shows what a hinge means when `k > 3`. In 4
    // dimensions the hinge holds 5 planes and it leaves 1 free, and a torque
    // in every plane at the same time only turns the free plane.
    for (const n of [3, 4]) {
      const D = dims(n);
      const world = quietWorld(n);
      const anchor = world.addBody(new Body(D, { shape: HyperSphere(D, 0.5), mass: 0 }));
      const b = world.addBody(new Body(D, {
        shape: HyperBox(D, new Array(n).fill(0.5)), mass: 1,
      }));
      world.createConstraint({ type: 'hinge', plane: 0 }, anchor, b);

      // Act: a torque in every plane at the same time.
      b.applyTorqueImpulse(new Float64Array(D.k).fill(0.4));
      for (let s = 0; s < 200; s += 1) world.step(1 / 240);

      // Assert
      assert.ok(Math.abs(b.w[0]) > 0.1, `n = ${n}: the free plane must turn, and it is ${b.w[0]}`);
      for (let p = 1; p < D.k; p += 1) {
        assert.ok(Math.abs(b.w[p]) < 1e-3, `n = ${n}: the plane ${p} must not turn, and it is ${b.w[p]}`);
      }
    }
  });

  it('should keep the complement basis near the basis of the last step', () => {
    // `refreshComplement` must not let the basis turn over. A basis that came
    // from the axes at each step would jump when the free plane passed a tie
    // in the order, and the impulse that the row keeps would then go along
    // another direction and add energy.
    const D = dims(4);
    const world = quietWorld(4);
    const a = world.addBody(new Body(D, { shape: HyperBox(D, [0.5, 0.5, 0.5, 0.5]), mass: 1 }));
    const b = world.addBody(new Body(D, { shape: HyperBox(D, [0.5, 0.5, 0.5, 0.5]), mass: 1 }));
    const joint = world.createConstraint({ type: 'hinge', plane: 0 }, a, b);
    a.setAngularVelocity([0, 0, 0, 2, 0, 0]);

    let previous = joint.angularBasis.map((u) => Float64Array.from(u));
    for (let s = 0; s < 400; s += 1) {
      world.step(1 / 240);
      assert.equal(joint.angularBasis.length, D.k - 1, 'the basis must hold k - 1 bivectors');
      for (let c = 0; c < previous.length; c += 1) {
        assert.ok(dot(joint.angularBasis[c], previous[c]) > 0.99,
          `at the step ${s} the bivector ${c} jumped`);
      }
      previous = joint.angularBasis.map((u) => Float64Array.from(u));
    }
  });
});

describe('the subspace lock', () => {
  it('should hold one axis of the position and leave the others free, in 4 dimensions', () => {
    const D = dims(4);
    const world = new World({
      dimensions: 4, gravity: [-1, -9.81, 0, 0], params: { allowSleep: false },
    });
    const b = world.addBody(new Body(D, {
      shape: HyperBox(D, [0.5, 0.5, 0.5, 0.5]), mass: 1, position: [0, 2, 0, 0],
    }));
    world.createConstraint({ type: 'subspace', lockAxes: [1], worldFrame: true }, b, null);

    for (let s = 0; s < 240; s += 1) world.step(1 / 120);

    assertClose(b.x[1], 2, 'the locked axis must not move', 1e-3);
    const t = 240 / 120;
    assertClose(b.x[0], -0.5 * t * t, 'the free axis must fall as usual', 2e-2);
  });

  it('should hold the named planes of the rotation', () => {
    const D = dims(4);
    const world = quietWorld(4);
    const b = world.addBody(new Body(D, {
      shape: HyperBox(D, [0.5, 0.7, 0.3, 0.6]), mass: 1,
    }));
    b.setAngularVelocity([0.8, 0.3, -0.5, 0.2, 0.6, -0.1]);
    world.createConstraint({ type: 'subspace', lockPlanes: [0, 3], worldFrame: true }, b, null);

    for (let s = 0; s < 240; s += 1) world.step(1 / 120);

    assert.ok(Math.abs(b.w[0]) < 1e-3, `the plane 0 must not turn, and it is ${b.w[0]}`);
    assert.ok(Math.abs(b.w[3]) < 1e-3, `the plane 3 must not turn, and it is ${b.w[3]}`);
    let free = 0;
    for (const p of [1, 2, 4, 5]) free = Math.max(free, Math.abs(b.w[p]));
    assert.ok(free > 0.1, 'the free planes must still turn');
  });

  it('should take the complement when you name the free axes', () => {
    const D = dims(4);
    const world = quietWorld(4);
    const b = world.addBody(ball(D, 1, [0, 0, 0, 0]));
    const joint = world.createConstraint(
      { type: 'subspace', freeAxes: [0], worldFrame: true }, b, null,
    );
    assert.equal(joint.linearLocal.length, 3, 'one free axis leaves three held axes');
    b.setLinearVelocity([1, 1, 1, 1]);
    for (let s = 0; s < 120; s += 1) world.step(1 / 120);
    assert.ok(b.x[0] > 0.9, 'the free axis moves');
    for (const i of [1, 2, 3]) {
      assert.ok(Math.abs(b.x[i]) < 1e-2, `the axis ${i} must be held, and it is ${b.x[i]}`);
    }
  });
});

describe('the joints and the world', () => {
  it('should give no NaN when a row has no effective mass', () => {
    // Four cases in which a row would divide by zero: two static bodies, a
    // static body joined to the world, a distance joint with the two anchors
    // at the same point, and a body that cannot turn with a weld.
    const D = dims(3);
    const world = quietWorld(3);
    const s1 = world.addBody(new Body(D, { shape: HyperSphere(D, 0.5), mass: 0 }));
    const s2 = world.addBody(new Body(D, { shape: HyperSphere(D, 0.5), mass: 0, position: [2, 0, 0] }));
    const free = world.addBody(ball(D, 1, [0, 4, 0]));
    const stiff = world.addBody(new Body(D, {
      shape: HyperSphere(D, 0.5), mass: 1, position: [0, 8, 0],
      angularFactor: new Float64Array(D.k),
    }));
    world.createConstraint({ type: 'point' }, s1, s2);
    world.createConstraint({ type: 'point' }, s1, null);
    world.createConstraint({ type: 'distance', rest: 0 }, free, free);
    world.createConstraint({ type: 'fixed' }, stiff, null);

    for (let s = 0; s < 10; s += 1) world.step(1 / 120);

    for (const b of world.bodies) {
      for (const field of ['x', 'R', 'v', 'L', 'w']) {
        for (const value of b[field]) {
          assert.ok(Number.isFinite(value), `the body ${b.id} holds a bad ${field}`);
        }
      }
    }
  });

  it('should take the joints out with the body', () => {
    const D = dims(3);
    const world = quietWorld(3);
    const a = world.addBody(ball(D, 1, [0, 0, 0]));
    const b = world.addBody(ball(D, 1, [2, 0, 0]));
    world.createConstraint({ type: 'point', localA: [1, 0, 0], localB: [-1, 0, 0] }, a, b);
    assert.equal(world.constraints.length, 1);

    world.removeBody(a);

    assert.equal(world.constraints.length, 0, 'the joint must go with the body');
    for (let s = 0; s < 10; s += 1) world.step(1 / 120);
    for (const value of b.x) assert.ok(Number.isFinite(value), 'the other body stays good');
  });

  it('should stop two bodies that a joint holds from touching', () => {
    const D = dims(3);
    function contactsWith(collideConnected) {
      const world = quietWorld(3);
      const a = world.addBody(ball(D, 1, [0, 0, 0]));
      const b = world.addBody(ball(D, 1, [0.4, 0, 0]));
      world.createConstraint({ type: 'point', collideConnected }, a, b);
      world.step(1 / 120);
      return world.contacts.length;
    }
    assert.equal(contactsWith(false), 0, 'the default must stop the contact');
    assert.ok(contactsWith(true) > 0, 'collideConnected true must let them touch');
  });

  it('should keep the joint when setMass builds a new body', () => {
    // `setMass` builds a new `Body`, because the inertia comes from the mass
    // at build time. A joint that kept the old body would pull a body that
    // nothing integrates.
    const PhysiN = loadPhysiNFor(false);
    const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
    scene.execute('addBody', {
      id: 100, shape: { type: 'sphere', n: 3, radius: 0.5 }, mass: 0, position: [0, 5, 0],
    });
    scene.execute('addBody', {
      id: 101, shape: { type: 'sphere', n: 3, radius: 0.5 }, mass: 1, position: [2, 5, 0],
    });
    scene.execute('addConstraint', {
      id: 0, type: 'point', a: 100, b: 101, localA: [0, 0, 0], localB: [-2, 0, 0],
    });
    for (let s = 0; s < 60; s += 1) scene.execute('simulate', { timeStep: 1 / 120 });
    scene.execute('setMass', { id: 101, value: 4 });
    for (let s = 0; s < 200; s += 1) scene.execute('simulate', { timeStep: 1 / 120 });

    const world = scene._engine.world;
    const bob = world.bodies.find((b) => b.id === 101);
    assert.equal(bob.mass, 4, 'the mass changed');
    assert.equal(world.constraints.length, 1, 'the joint is still there');
    assert.equal(world.constraints[0].b, bob, 'the joint holds the new body');
    const d = Math.hypot(bob.x[0], bob.x[1] - 5, bob.x[2]);
    assertClose(d, 2, 'the joint still holds the length', 2e-2);
  });

  it('should not kick a body that a teleport moves', () => {
    // `updateTransform` moves a body without a velocity. The impulses that the
    // rows keep are then wrong, thus the engine makes them zero.
    const PhysiN = loadPhysiNFor(false);
    const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
    scene.execute('addBody', {
      id: 200, shape: { type: 'sphere', n: 3, radius: 0.5 }, mass: 0, position: [0, 5, 0],
    });
    scene.execute('addBody', {
      id: 201, shape: { type: 'sphere', n: 3, radius: 0.5 }, mass: 1, position: [2, 5, 0],
    });
    scene.execute('addConstraint', {
      id: 0, type: 'point', a: 200, b: 201, localA: [0, 0, 0], localB: [-2, 0, 0],
    });
    for (let s = 0; s < 120; s += 1) scene.execute('simulate', { timeStep: 1 / 120 });
    scene.execute('updateTransform', { id: 201, position: [52, 5, 0] });
    scene.execute('simulate', { timeStep: 1 / 120 });

    const world = scene._engine.world;
    for (const b of world.bodies) {
      assert.ok(norm(b.v) < 50, `the body ${b.id} was kicked to the speed ${norm(b.v)}`);
    }
  });

  it('should tell the plugin when a joint names a body that is not there', () => {
    const PhysiN = loadPhysiNFor(false);
    const scene = new PhysiN.Scene({ dimensions: 3 });
    scene._constraints[7] = { id: 7 };
    const seen = [];
    const original = console.error;
    console.error = (m) => seen.push(m);
    try {
      scene.execute('addConstraint', { id: 7, type: 'point', a: 999, b: null });
    } finally {
      console.error = original;
    }
    assert.equal(scene._engine.world.constraints.length, 0, 'the engine took no joint');
    assert.equal(seen.length, 1, 'the plugin said what went wrong');
    assert.ok(!scene._constraints[7], 'the plugin dropped the joint');
  });

  it('should throw when a joint goes before scene.add', () => {
    const PhysiN = loadPhysiNFor(false);
    const scene = new PhysiN.Scene({ dimensions: 3 });
    const a = new PhysiN.SphereMesh(new THREE.SphereGeometry(0.5), {}, 1);
    const b = new PhysiN.SphereMesh(new THREE.SphereGeometry(0.5), {}, 1);
    scene.add(a);
    assert.throws(
      () => scene.addConstraint(new PhysiN.PointJoint(a, b, [0, 0, 0], [0, 0, 0])),
      /scene.add/,
      'a joint before scene.add must say so, and not go away in silence',
    );
  });
});

describe('the two bundles', () => {
  it('should give the same answer in the main thread and in the worker', () => {
    // `physin.js` and `physin_worker.js` must hold the same engine, word for
    // word. A change made on one side only shows here.
    function run(inWorker) {
      const PhysiN = loadPhysiNFor(inWorker);
      const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
      const post = new PhysiN.SphereMesh(new THREE.SphereGeometry(0.5), {}, 0);
      post.setPositionN([0, 6, 0]);
      scene.add(post);
      const arm = new PhysiN.SphereMesh(new THREE.SphereGeometry(0.5), {}, 1);
      arm.setPositionN([2, 6, 0]);
      scene.add(arm);
      const hand = new PhysiN.SphereMesh(new THREE.SphereGeometry(0.5), {}, 2);
      hand.setPositionN([4, 6, 0]);
      scene.add(hand);
      scene.addConstraint(new PhysiN.PointJoint(post, arm, [0, 0, 0], [-2, 0, 0]));
      scene.addConstraint(new PhysiN.HingeJoint(arm, hand, [1, 0, 0], [-1, 0, 0], { plane: 0 }));
      for (let s = 0; s < 120; s += 1) scene.simulate(1 / 120, 1);
      return [Array.from(arm.getPositionN()), Array.from(hand.getPositionN())];
    }
    const main = run(false);
    const worker = run(true);
    for (let b = 0; b < main.length; b += 1) {
      for (let i = 0; i < main[b].length; i += 1) {
        assertClose(worker[b][i], main[b][i],
          `the body ${b} on the axis ${i}`, 1e-6);
      }
    }
    assert.ok(Math.abs(main[0][1] - 6) > 0.05, 'the arm must have moved, or the test proves nothing');
  });
});

describe('the soft constraint', () => {
  it('should be the same, number for number, when hertz is 0', () => {
    // The line of `solveJoint` is
    //   -massScale * (cv - bias) / km - impulseScale * impulse
    // and it must be BIT identical to the old `-(cv - bias) / km` when the
    // mass scale is 1 and the impulse scale is 0. A multiply by 1 and a
    // subtract of 0 are exact in IEEE 754, but only in that order. This test
    // holds the order. It uses `Object.is`, and not a tolerance.
    const run = (params) => {
      const world = quietWorld(3, Object.assign({ constraintBias: 0.2 }, params));
      world.gravity[1] = -9.81;
      const anchor = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
      const b = world.addBody(ball(D3, 1, [2, 5, 0]));
      const c = world.addBody(ball(D3, 3, [4, 5, 0]));
      world.createConstraint({ type: 'point', localA: [0, 0, 0], localB: [-2, 0, 0] }, anchor, b);
      world.createConstraint({ type: 'fixed', localA: [1, 0, 0], localB: [-1, 0, 0] }, b, c);
      for (let s = 0; s < 300; s += 1) world.step(1 / 120);
      const out = [];
      for (const body of world.bodies) out.push(...body.x, ...body.v, ...body.L, ...body.R);
      return out;
    };
    const a = run({});
    const b = run({ constraintHertz: 0 });
    assert.equal(a.length, b.length);
    for (let i = 0; i < a.length; i += 1) {
      assert.ok(Object.is(a[i], b[i]), `at the index ${i}, ${a[i]} is not exactly ${b[i]}`);
    }
  });

  it('should stretch by g / (2 pi f)^2, which is the stretch of a spring', () => {
    // A soft joint is a true spring, thus its stretch under a weight follows
    // the formula of a spring. This is the test that the three factors of
    // `softness` are right, and not only that they are stable.
    for (const hertz of [3, 5, 10, 60]) {
      const world = quietWorld(3, { constraintSlop: 0, constraintHertz: hertz, constraintDamping: 1 });
      world.gravity[1] = -9.81;
      const anchor = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
      const b = world.addBody(ball(D3, 1, [0, 3, 0]));
      world.createConstraint({ type: 'distance', rest: 2 }, anchor, b);
      for (let s = 0; s < 1800; s += 1) world.step(1 / 240);
      const stretch = Math.abs(Math.hypot(b.x[0], b.x[1] - 5, b.x[2]) - 2);
      const theory = 9.81 / ((2 * Math.PI * hertz) ** 2);
      assertClose(stretch / theory, 1, `the stretch at ${hertz} Hz`, 0.05);
    }
  });

  it('should overshoot when the damping is small, and not when it is 1', () => {
    const overshoot = (zeta) => {
      const world = quietWorld(3, { constraintSlop: 0, constraintHertz: 4, constraintDamping: zeta });
      world.gravity[1] = -9.81;
      const anchor = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
      const b = world.addBody(ball(D3, 1, [0, 3, 0]));
      world.createConstraint({ type: 'distance', rest: 2 }, anchor, b);
      let peak = 0;
      for (let s = 0; s < 2400; s += 1) {
        world.step(1 / 240);
        peak = Math.max(peak, Math.abs(Math.hypot(b.x[0], b.x[1] - 5, b.x[2]) - 2));
      }
      const settled = Math.abs(Math.hypot(b.x[0], b.x[1] - 5, b.x[2]) - 2);
      return peak / settled - 1;
    };
    assert.ok(overshoot(0.2) > 0.3, 'a damping of 0.2 must overshoot by more than 30 per cent');
    assert.ok(overshoot(1) < 0.05, 'a damping of 1 must almost not overshoot');
  });

  it('should hold a frequency that is too high for the step', () => {
    // Above `constraintHertzRatio / dt` the spring is not stable. The solver
    // holds the frequency there, thus an absurd value does not blow up.
    const world = quietWorld(3, { constraintHertz: 1e5 });
    world.gravity[1] = -9.81;
    const anchor = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
    const b = world.addBody(ball(D3, 1, [0, 3, 0]));
    world.createConstraint({ type: 'distance', rest: 2 }, anchor, b);
    for (let s = 0; s < 600; s += 1) world.step(1 / 120);
    for (const v of b.x) assert.ok(Number.isFinite(v), 'the body must stay finite');
    assert.ok(Math.abs(Math.hypot(b.x[0], b.x[1] - 5, b.x[2]) - 2) < 1e-2, 'and it must still hold');
  });
});

describe('the hinge angle, the limit and the motor', () => {
  it('should give the angle exactly, at any angle, with a turn in another plane', () => {
    // The first order logarithm of `angularError` is NOT exact here. This is,
    // because the free part of a hinge is ONE plane, and a rotation in one
    // plane commutes with itself.
    const D = dims(4);
    const world = quietWorld(4);
    const a = world.addBody(ball(D, 1, new Float64Array(4)));
    const b = world.addBody(ball(D, 1, new Float64Array(4)));
    a.R.set(rotorFromPlane(D, 2, 3, 0.7));
    a.updateDerived();
    b.R.set(a.R);
    b.updateDerived();
    const joint = world.createConstraint({ type: 'hinge', plane: 0 }, a, b);
    for (const angle of [0.3, 2, 3, -1.5]) {
      b.R.set(rotorMul(D, rotorFromPlane(D, 0, 1, angle), a.R));
      b.updateDerived();
      assertClose(nd.hingeAngle(D, joint), angle, `the angle ${angle}`, 1e-12);
    }
  });

  it('should refuse a plane that holds two planes at the same time', () => {
    // `e_xy + e_zw` is not simple. It is not a hinge, and the angle of a hinge
    // would have no meaning for it. Round 1 took it in silence.
    const D = dims(4);
    const world = quietWorld(4);
    const a = world.addBody(ball(D, 1, new Float64Array(4)));
    const b = world.addBody(ball(D, 1, new Float64Array(4)));
    assert.throws(
      () => world.createConstraint({ type: 'hinge', plane: [1, 0, 0, 0, 0, 1] }, a, b),
      /not simple/,
    );
  });

  it('should hold the arm at the limit, and let it swing past without one', () => {
    const arm = (opts) => {
      const world = new World({ dimensions: 3, gravity: [0, -9.81, 0], params: { allowSleep: false } });
      const base = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
      const bar = world.addBody(new Body(D3, {
        shape: HyperBox(D3, [1, 0.1, 0.1]), mass: 1, position: [1, 5, 0],
      }));
      const joint = world.createConstraint(Object.assign(
        { type: 'hinge', plane: 0, localA: [0, 0, 0], localB: [-1, 0, 0] }, opts,
      ), base, bar);
      let lowest = 0;
      for (let s = 0; s < 720; s += 1) {
        world.step(1 / 120);
        lowest = Math.min(lowest, joint.angle);
      }
      return { joint, lowest };
    };
    const free = arm({});
    assert.ok(free.lowest < -3, `with no limit the arm must hang down, and it reached ${free.lowest}`);
    const held = arm({ lowerAngle: -0.5, upperAngle: 0.5 });
    assertClose(held.joint.angle, -0.5, 'the limit must hold the arm', 5e-3);
  });

  it('should not change the answer when the limits are far away', () => {
    // A limit that is not engaged must set `km` to 0, thus the row goes away
    // and it changes nothing at all.
    //
    // The run is short on purpose. The angle wraps at pi, thus an arm that
    // swings all the way down reaches -3.14 and a limit at -3 WOULD engage.
    // The test asserts that the angle stayed well inside the limits.
    const run = (opts) => {
      const world = new World({ dimensions: 3, gravity: [0, -9.81, 0], params: { allowSleep: false } });
      const base = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
      const bar = world.addBody(new Body(D3, {
        shape: HyperBox(D3, [1, 0.1, 0.1]), mass: 1, position: [1, 5, 0],
      }));
      const joint = world.createConstraint(Object.assign(
        { type: 'hinge', plane: 0, localA: [0, 0, 0], localB: [-1, 0, 0] }, opts,
      ), base, bar);
      let reach = 0;
      for (let s = 0; s < 60; s += 1) {
        world.step(1 / 120);
        reach = Math.max(reach, Math.abs(joint.angle));
      }
      return { place: Array.from(bar.x), reach };
    };
    const none = run({});
    const wide = run({ lowerAngle: -3, upperAngle: 3 });
    assert.ok(none.reach < 1, `the arm must stay well inside the limits, and it reached ${none.reach}`);
    for (let i = 0; i < 3; i += 1) assertClose(wide.place[i], none.place[i], `the axis ${i}`, 1e-9);
  });

  it('should keep the same row objects when a limit engages', () => {
    // The row of the motor and the row of the limit are always there. A count
    // that changed would throw away every impulse that the joint keeps, at the
    // moment that the limit engages.
    const world = new World({ dimensions: 3, gravity: [0, -9.81, 0], params: { allowSleep: false } });
    const base = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
    const bar = world.addBody(new Body(D3, {
      shape: HyperBox(D3, [1, 0.1, 0.1]), mass: 1, position: [1, 5, 0],
    }));
    const joint = world.createConstraint(
      { type: 'hinge', plane: 0, localA: [0, 0, 0], localB: [-1, 0, 0], lowerAngle: -0.5, upperAngle: 0.5 },
      base, bar,
    );
    world.step(1 / 120);
    const rows = joint.rows;
    const each = rows.slice();
    const limit = rows[rows.length - 1];
    assert.equal(limit.km, 0, 'the limit starts off');
    for (let s = 0; s < 240; s += 1) world.step(1 / 120);
    assert.ok(limit.km > 0, 'the limit engaged');
    assert.equal(joint.rows, rows, 'the array of the rows is the same object');
    for (let c = 0; c < each.length; c += 1) {
      assert.equal(joint.rows[c], each[c], `the row ${c} is the same object`);
    }
  });

  it('should reach the speed of the motor, and stall at a limit', () => {
    const build = (opts) => {
      const world = quietWorld(3);
      const base = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0 }));
      const bar = world.addBody(new Body(D3, {
        shape: HyperBox(D3, [1, 0.1, 0.1]), mass: 1, position: [1, 0, 0],
      }));
      const joint = world.createConstraint(Object.assign(
        { type: 'hinge', plane: 0, localA: [0, 0, 0], localB: [-1, 0, 0] }, opts,
      ), base, bar);
      return { world, bar, joint };
    };
    const free = build({ motorSpeed: 2, maxMotorTorque: 50 });
    for (let s = 0; s < 240; s += 1) free.world.step(1 / 120);
    assertClose(free.bar.w[0], 2, 'the motor must reach its speed', 0.02);

    const stalled = build({ motorSpeed: 2, maxMotorTorque: 50, upperAngle: 0.5 });
    for (let s = 0; s < 600; s += 1) stalled.world.step(1 / 120);
    assertClose(stalled.joint.angle, 0.5, 'the motor must stall at the limit', 5e-3);
    const rows = stalled.joint.rows;
    assert.ok(rows[rows.length - 2].impulse > 0, 'the motor pushes one way');
    assert.ok(rows[rows.length - 1].impulse < 0, 'and the limit holds the other way');
  });

  it('should give the motor the same strength at any count of sub-steps', () => {
    // The clamp of the motor is `maxMotorTorque * dt`, thus a TORQUE. An
    // impulse would make the motor stronger when `subSteps` grew.
    const stall = (subSteps) => {
      const world = new World({
        dimensions: 3, gravity: [0, -9.81, 0], params: { allowSleep: false, subSteps },
      });
      const base = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
      const bar = world.addBody(new Body(D3, {
        shape: HyperBox(D3, [1, 0.1, 0.1]), mass: 1, position: [1, 5, 0],
      }));
      const joint = world.createConstraint({
        type: 'hinge', plane: 0, localA: [0, 0, 0], localB: [-1, 0, 0],
        motorSpeed: 5, maxMotorTorque: 4,
      }, base, bar);
      for (let s = 0; s < 600; s += 1) world.step(1 / 120);
      return joint.angle;
    };
    // An impulse budget in the place of a torque budget would make the motor
    // four times as strong at four sub-steps, and the two angles would then be
    // nowhere near each other.
    assertClose(stall(4), stall(1), 'the angle where a weak motor stalls', 3e-2);
  });
});

describe('a joint that breaks', () => {
  const hang = (opts, subSteps) => {
    const world = new World({
      dimensions: 3, gravity: [0, -9.81, 0], params: { allowSleep: false, subSteps },
    });
    const anchor = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0, position: [0, 5, 0] }));
    const b = world.addBody(ball(D3, 10, [0, 3, 0]));
    const joint = world.createConstraint(Object.assign({ type: 'distance', rest: 2 }, opts), anchor, b);
    const seen = [];
    world.on('constraintBroken', (j) => seen.push(j));
    for (let s = 0; s < 240; s += 1) world.step(1 / 120);
    return { world, joint, seen };
  };

  it('should break under a load that is above the limit, one time', () => {
    const r = hang({ breakForce: 20 }, 1);
    assert.equal(r.seen.length, 1, 'the event goes out exactly one time');
    assert.equal(r.seen[0], r.joint, 'and it gives the joint');
    assert.equal(r.world.constraints.length, 0, 'the world took the joint out');
  });

  it('should not break under a load that it can carry', () => {
    const r = hang({ breakForce: 1e6 }, 1);
    assert.equal(r.seen.length, 0);
    assert.equal(r.world.constraints.length, 1);
  });

  it('should measure a force, thus the limit does not change with the sub-steps', () => {
    // `row.impulse` is the impulse of ONE substep. A limit on it would mean
    // something different at each count of sub-steps. `impulse / dt` is a
    // force, and it does not change. The load here is 10 kg * 9.81 = 98.1 N.
    for (const subSteps of [1, 8]) {
      assert.equal(hang({ breakForce: 20 }, subSteps).seen.length, 1,
        `a weak joint must break at ${subSteps} sub-steps`);
      assert.equal(hang({ breakForce: 1e6 }, subSteps).seen.length, 0,
        `a strong joint must hold at ${subSteps} sub-steps`);
    }
    const measured = hang({ breakForce: 1e6 }, 1);
    assertClose(measured.joint.linearImpulse * 120, 98.1, 'the measured force', 0.02);
  });

  it('should not let a motor break its own joint', () => {
    // The row of the motor does not count toward the load. A motor that runs
    // free carries a large impulse by design.
    const world = quietWorld(3);
    const base = world.addBody(new Body(D3, { shape: HyperSphere(D3, 0.2), mass: 0 }));
    const bar = world.addBody(new Body(D3, {
      shape: HyperBox(D3, [1, 0.1, 0.1]), mass: 1, position: [1, 0, 0],
    }));
    world.createConstraint({
      type: 'hinge', plane: 0, localA: [0, 0, 0], localB: [-1, 0, 0],
      motorSpeed: 2, maxMotorTorque: 1000, breakTorque: 1,
    }, base, bar);
    let broke = 0;
    world.on('constraintBroken', () => { broke += 1; });
    for (let s = 0; s < 600; s += 1) world.step(1 / 120);
    assert.equal(broke, 0, 'a motor that runs free must not break its joint');
  });
});

describe('the drift of a lock of the rotation', () => {
  /** A 4D body that bounces on the ground with two planes held. */
  const tumble = (bias, iterations) => {
    const D = dims(4);
    const world = new World({
      dimensions: 4, gravity: [0, -9.81, 0, 0],
      params: { allowSleep: false, constraintBias: bias, iterations },
    });
    world.addBody(new Body(D, { shape: nd.HalfSpace(D, [0, 1, 0, 0], 0), mass: 0 }));
    const b = world.addBody(new Body(D, {
      shape: HyperBox(D, [0.5, 0.7, 0.3, 0.6]), mass: 1, position: [0, 3, 0, 0],
    }));
    b.setAngularVelocity([2, 0.3, -0.5, 0.2, 0.6, -0.1]);
    const joint = world.createConstraint(
      { type: 'subspace', lockPlanes: [0, 3], worldFrame: true }, b, null,
    );
    let worst = 0;
    for (let s = 0; s < 1800; s += 1) {
      world.step(1 / 120);
      let m = 0;
      for (const row of joint.rows) if (row.kind === 'angular') m += row.drift * row.drift;
      worst = Math.max(worst, Math.sqrt(m));
    }
    return { worst, joint };
  };

  it('should use the drift when two or more planes stay free', () => {
    const r = tumble(0.2, 10);
    assert.equal(r.joint.useDrift, true, 'a lock of two planes of six must use the drift');
  });

  it('should hold the drift smaller than with no correction', () => {
    // With few iterations the solver cannot fully hold the velocity, thus a
    // drift builds up. That is the case that this feature is for.
    const off = tumble(0, 2);
    const on = tumble(0.2, 2);
    assert.ok(on.worst < off.worst / 4,
      `the drift with a bias is ${on.worst} and with no bias it is ${off.worst}`);
  });

  it('should not lose the drift when the warm start is off', () => {
    // `prepareRows` makes the impulses zero when `useWarmStart` is false. It
    // must NOT make the drift zero as well, or the correction goes away with
    // no message.
    const D = dims(4);
    const world = new World({
      dimensions: 4, gravity: [0, -9.81, 0, 0],
      params: { allowSleep: false, useWarmStart: false, iterations: 2 },
    });
    world.addBody(new Body(D, { shape: nd.HalfSpace(D, [0, 1, 0, 0], 0), mass: 0 }));
    const b = world.addBody(new Body(D, {
      shape: HyperBox(D, [0.5, 0.7, 0.3, 0.6]), mass: 1, position: [0, 3, 0, 0],
    }));
    b.setAngularVelocity([2, 0.3, -0.5, 0.2, 0.6, -0.1]);
    const joint = world.createConstraint(
      { type: 'subspace', lockPlanes: [0, 3], worldFrame: true }, b, null,
    );
    for (let s = 0; s < 600; s += 1) world.step(1 / 120);
    let m = 0;
    for (const row of joint.rows) if (row.kind === 'angular') m += Math.abs(row.drift);
    assert.ok(m > 0, 'the drift must survive a step with no warm start');
  });

  it('should keep the drift at zero when the geometric error is usable', () => {
    const world = quietWorld(4);
    const b = world.addBody(new Body(dims(4), {
      shape: HyperBox(dims(4), [0.5, 0.5, 0.5, 0.5]), mass: 1,
    }));
    b.setAngularVelocity([0.5, 0.2, 0.1, 0, 0, 0]);
    const joint = world.createConstraint({ type: 'fixed' }, b, null);
    for (let s = 0; s < 100; s += 1) world.step(1 / 120);
    assert.equal(joint.useDrift, false, 'a weld has a geometric error');
    for (const row of joint.rows) {
      if (row.kind === 'angular') assert.equal(row.drift, 0, 'and it keeps no drift');
    }
  });

  it('should lose the drift when the joint is put back on', () => {
    const D = dims(4);
    const world = quietWorld(4);
    const b = world.addBody(new Body(D, { shape: HyperBox(D, [0.5, 0.5, 0.5, 0.5]), mass: 1 }));
    const joint = world.createConstraint(
      { type: 'subspace', lockPlanes: [0, 3], worldFrame: true }, b, null,
    );
    world.step(1 / 120);
    for (const row of joint.rows) row.drift = 0.5;
    nd.resetConstraint(joint);
    for (const row of joint.rows) assert.equal(row.drift, 0, 'resetConstraint clears the drift');
  });
});
