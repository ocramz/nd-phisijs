/**
 * test/spec/body.test.js -- src/nd/body/body.js
 *
 * The state of a rigid body, and the operations that change it. The laws
 * here are the laws of mechanics: an impulse keeps the total momentum, a
 * rigid body does not change its size, and the box of the broad phase holds
 * the body.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyVector, anyUnitVector, anyBivector, anyRotor, anyPositive, anyHalfExtents,
  anyBody, anyPosition, anyAngularVelocity,
} from '../lib/arbitraries.js';
import { assertClose, assertArrayClose, dot, norm, maxAbs, assertScaled } from '../lib/numeric.js';

const { dims } = nd;

/** A body of `n` dimensions, with the tables of that `n`. */
const anyBodyOf = (min = 2, max = 5) => anyN(min, max)
  .chain((n) => anyBody(dims(n)).map((body) => ({ D: dims(n), body })));

describe('localToWorld', () => {
  it('should always go back with worldToLocal', () => {
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyPosition(D.n),
    )), ([{ D, body }, p]) => {
      // Act
      const world = body.localToWorld(p, new Float64Array(D.n));
      const back = body.worldToLocal(world, new Float64Array(D.n));

      // Assert
      assertArrayClose(back, p, 'a point that goes there and comes back', 1e-9);
    }));
  });

  it('should never change the length of a direction', () => {
    // `Rm` is a rotation matrix, thus the change of frame is an isometry.
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyVector(D.n),
    )), ([{ D, body }, v]) => {
      assertClose(norm(body.localToWorldDir(v, new Float64Array(D.n))), norm(v),
        'the length of a direction', 1e-9);
    }));
  });
});

describe('pointVelocity', () => {
  it('should give the velocity of the center at the center', () => {
    // Arrange
    const D = dims(4);
    const body = new nd.Body(D, {
      shape: nd.HyperBox(D, [1, 1, 1, 1]), mass: 2, velocity: [1, 2, 3, 4],
    });

    // Act
    const u = body.pointVelocity(new Float64Array(4));

    // Assert
    assertArrayClose(u, [1, 2, 3, 4], 'the velocity at the center of mass');
  });

  it('should never let a point move away from the center of mass', () => {
    // The part of the velocity that comes from the turn is r . w, and that is
    // always normal to r. Thus the distance from the point to the center of
    // mass does not change: the body is rigid.
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyVector(D.n),
    )), ([{ D, body }, r]) => {
      // Act
      const u = body.pointVelocity(r, new Float64Array(D.n));

      // Assert -- the tolerance is relative to `|r|^2`, and `|r|^2` goes to
      // zero for an `r` of the size 1e-160, which `anyVector` gives. The
      // relative part is then zero, and one bit of the last place of the sum
      // stops the test. Thus add an absolute floor. It is 1e-318, thus it
      // hides no error that has a physical size.
      const spin = Float64Array.from(u, (v, i) => v - body.v[i]);
      const tol = 1e-9 * norm(r) ** 2 * (1 + norm(body.w)) + Number.MIN_VALUE * 1024;
      assert.ok(Math.abs(dot(spin, r)) <= tol,
        `the turn moves the point by ${dot(spin, r)} along the offset`);
    }));
  });
});

describe('applyImpulse', () => {
  it('should never change the total momentum of a pair of bodies', () => {
    // Newton: the two bodies of a contact take the same impulse, with the
    // sign turned around, at the same point. Thus the total linear momentum
    // and the total angular momentum about any origin stay the same. The
    // whole solver is built on this. See ND-PHYSICS.md, A11.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(dims(n)), anyBody(dims(n)), anyBody(dims(n)), anyPosition(n), anyVector(n),
    )), ([D, a, b, point, j]) => {
      // Arrange
      const before = totalMomentum(D, [a, b]);
      const rA = Float64Array.from(point, (v, i) => v - a.x[i]);
      const rB = Float64Array.from(point, (v, i) => v - b.x[i]);

      // Act
      a.applyImpulse(Float64Array.from(j, (v) => -v), rA);
      b.applyImpulse(j, rB);

      // Assert
      const after = totalMomentum(D, [a, b]);
      const scale = maxAbs(j) * (1 + maxAbs(point)) * 10;
      assertScaled(after.linear, before.linear, scale, 'the total linear momentum', 1e-9);
      assertScaled(after.angular, before.angular, scale, 'the total angular momentum', 1e-9);
    }));
  });

  it('should never make a static body move', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(dims(n)), anyVector(n), anyPosition(n),
    )), ([D, j, r]) => {
      // Arrange
      const wall = new nd.Body(D, { shape: nd.HalfSpace(D, unitAxis(D.n, 1), 0) });

      // Act
      wall.applyImpulse(j, r);

      // Assert
      assert.ok(wall.v.every((v) => v === 0), 'the velocity stays zero');
      assert.ok(wall.L.every((v) => v === 0), 'the angular momentum stays zero');
    }));
  });
});

describe('applyCentralImpulse', () => {
  it('should never make a body start to turn', () => {
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyVector(D.n),
    )), ([{ body }, j]) => {
      // Arrange
      const before = Float64Array.from(body.L);

      // Act
      body.applyCentralImpulse(j);

      // Assert
      assertArrayClose(body.L, before, 'the angular momentum does not change');
    }));
  });

  it('should always change the energy by the work of the impulse', () => {
    // The work-energy law: the change of the kinetic energy is
    // j . (v_before + v_after) / 2.
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyVector(D.n, fc.double({ min: -50, max: 50, noNaN: true })),
    )), ([{ body }, j]) => {
      // Arrange
      const before = Float64Array.from(body.v);
      const energyBefore = body.kineticEnergy();

      // Act
      body.applyCentralImpulse(j);

      // Assert
      const work = dot(j, Float64Array.from(before, (v, i) => (v + body.v[i]) / 2));
      assertClose(body.kineticEnergy() - energyBefore, work, 'the work of the impulse', 1e-6);
    }));
  });
});

describe('setAngularVelocity', () => {
  it('should always give the same angular velocity back', () => {
    // It builds the momentum with L = I' w, and `updateDerived` gives w back
    // with w = I'^-1 L. The two must be the inverse of each other.
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyAngularVelocity(D.n),
    )), ([{ D, body }, w]) => {
      // Act
      body.setAngularVelocity(w);
      const back = nd.linalg.matVec(body.invInertiaWorld, body.L, D.k, D.k);

      // Assert
      assertArrayClose(back, w, 'the angular velocity from the momentum', 1e-6);
    }));
  });
});

describe('kineticEnergy', () => {
  it('should never be less than zero', () => {
    // The energy is (m v.v + w.L) / 2. The second part is w . I' w, and the
    // inertia tensor is positive definite, thus it is never negative.
    fc.assert(fc.property(anyBodyOf(), ({ body }) => {
      assert.ok(body.kineticEnergy() >= 0, `the energy is ${body.kineticEnergy()}`);
    }));
  });

  it('should never change when the body only moves to another place', () => {
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyPosition(D.n),
    )), ([{ D, body }, place]) => {
      // Arrange
      const before = body.kineticEnergy();

      // Act
      body.x.set(place);

      // Assert
      assertClose(body.kineticEnergy(), before, 'the energy of a body that moves');
    }));
  });
});

describe('angularFactor', () => {
  it('should always hold a plane at zero', () => {
    // A factor of zero in one plane must stop the body from turning in that
    // plane, whatever impulse it takes.
    fc.assert(fc.property(anyN(3, 5).chain((n) => fc.tuple(
      fc.constant(dims(n)), anyHalfExtents(n), fc.nat({ max: dims(n).k - 1 }),
      fc.array(fc.tuple(anyVector(n), anyPosition(n)), { minLength: 1 }),
    )), ([D, h, plane, hits]) => {
      // Arrange
      const factor = new Float64Array(D.k).fill(1);
      factor[plane] = 0;
      const body = new nd.Body(D, { shape: nd.HyperBox(D, h), mass: 1, angularFactor: factor });

      // Act
      for (const [j, r] of hits) body.applyImpulse(j, r);

      // Assert
      assert.equal(body.L[plane], 0, 'the momentum of that plane stays zero');
    }));
  });
});

describe('aabb', () => {
  it('should never lose a point of the body', () => {
    // The broad phase drops a pair whose boxes do not meet. Thus a support
    // point outside the box gives a lost collision. In 4 and 5 dimensions the
    // box of a body that turns is the sum of 2n terms, and only a test of
    // each direction finds a defect there.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(dims(n)),
      anyBody(dims(n)),
      anyUnitVector(n),
    )), ([D, body, dir]) => {
      // Arrange
      const box = body.aabb(0);

      // Act -- the farthest point of the body along `dir`, in the world frame.
      const local = body.shape.support(body.worldToLocalDir(dir, new Float64Array(D.n)));
      const world = body.localToWorld(local, new Float64Array(D.n));

      // Assert
      for (let i = 0; i < D.n; i += 1) {
        assert.ok(world[i] >= box.min[i] - 1e-9 && world[i] <= box.max[i] + 1e-9,
          `the axis ${i}: the point at ${world[i]} is outside ${box.min[i]}..${box.max[i]}`);
      }
    }));
  });

  it('should always grow by the margin on each side', () => {
    fc.assert(fc.property(anyBodyOf().chain(({ D, body }) => fc.tuple(
      fc.constant({ D, body }), anyPositive(0, 1),
    )), ([{ D, body }, margin]) => {
      // Act
      const tight = body.aabb(0);
      const large = body.aabb(margin);

      // Assert
      for (let i = 0; i < D.n; i += 1) {
        assertClose(large.min[i], tight.min[i] - margin, `the low side of the axis ${i}`, 1e-9);
        assertClose(large.max[i], tight.max[i] + margin, `the high side of the axis ${i}`, 1e-9);
      }
    }));
  });

  it('should have no limit for a half space', () => {
    // Arrange
    const D = dims(3);
    const wall = new nd.Body(D, { shape: nd.HalfSpace(D, [0, 1, 0], 0) });

    // Act
    const box = wall.aabb(0);

    // Assert
    assert.ok(box.min.every((v) => v === -Infinity), 'a half space fills the space');
    assert.ok(box.max.every((v) => v === Infinity), 'on each side');
  });
});

// Helpers

/** The vector of the length 1 along the axis `i`. */
function unitAxis(n, i) {
  const e = new Float64Array(n);
  e[i] = 1;
  return e;
}

/**
 * The total linear momentum and the total angular momentum about the origin
 * of a list of bodies. The angular part holds the spin `L` of each body and
 * the part `x ^ p` of its movement about the origin.
 */
function totalMomentum(D, bodies) {
  const linear = new Float64Array(D.n);
  const angular = new Float64Array(D.k);
  for (const b of bodies) {
    const p = Float64Array.from(b.v, (v) => v * b.mass);
    for (let i = 0; i < D.n; i += 1) linear[i] += p[i];
    const orbit = nd.wedgeVec(D, b.x, p);
    for (let q = 0; q < D.k; q += 1) angular[q] += b.L[q] + orbit[q];
  }
  return { linear, angular };
}
