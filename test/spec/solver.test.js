/**
 * test/spec/solver.test.js -- src/nd/resolve/solver.js
 *
 * The contact solver. The contacts come from `collide`, and not from a
 * private constructor, thus these tests do not touch an internal detail.
 *
 * The three laws here are the laws of A11: an impulse is internal, thus it
 * keeps the momentum of the pair; a contact pushes and it cannot pull; and
 * the friction stays inside the cone of Coulomb.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyRadius, anyPositive, anyUnitVector, anyVector, anyRotor, anyAngularVelocity,
} from '../lib/arbitraries.js';
import { assertClose, assertScaled, dot, norm, maxAbs } from '../lib/numeric.js';

const { dims, defaultParams } = nd;

describe('prepareContact', () => {
  it('should never give an effective mass that is less than the two inverse masses', () => {
    // kn = 1/mA + 1/mB + (rA ^ nrm) . I'^-1 (rA ^ nrm) + the same for B.
    // The inertia part is never negative, because the inverse inertia is
    // positive definite. Thus a turn only makes a body softer at a contact,
    // and never harder. See ND-PHYSICS.md, A11.
    fc.assert(fc.property(anyContact(), (spec) => {
      // Arrange
      const { D, a, b, contacts } = makeContacts(spec);
      if (contacts.length === 0) return;

      // Act
      for (const c of contacts) nd.prepareContact(D, c, 1 / 60, defaultParams);

      // Assert
      for (const c of contacts) {
        assert.ok(c.kn >= a.invMass + b.invMass - 1e-9,
          `the effective mass ${c.kn} is less than ${a.invMass + b.invMass}`);
        for (const kt of c.kt) {
          assert.ok(kt >= a.invMass + b.invMass - 1e-9, 'the same along a tangent');
        }
      }
    }));
  });

  it('should give only the two inverse masses when the offset is along the normal', () => {
    // Two balls touch along the line of their centers, thus the offset of the
    // contact is parallel to the normal and the wedge product is zero. The
    // effective mass is then the mass of two points.
    const D = dims(4);
    const a = new nd.Body(D, { shape: nd.HyperSphere(D, 1), mass: 2, position: [0, 0, 0, 0] });
    const b = new nd.Body(D, { shape: nd.HyperSphere(D, 0.5), mass: 4, position: [1.4, 0, 0, 0] });
    const contacts = [];
    nd.collide(D, a, b, contacts, defaultParams);

    // Act
    nd.prepareContact(D, contacts[0], 1 / 60, defaultParams);

    // Assert
    assertClose(contacts[0].kn, 1 / 2 + 1 / 4, 'the effective mass of two balls that meet');
  });
});

describe('solveContact', () => {
  it('should never pull the two bodies together', () => {
    // A contact can push, and it cannot pull. The solver holds the TOTAL
    // impulse at zero or more, and it applies the change only. That clamp is
    // what makes the sequential impulse method work.
    fc.assert(fc.property(anyContact(), (spec) => {
      // Arrange
      const { D, contacts } = makeContacts(spec);
      for (const c of contacts) nd.prepareContact(D, c, 1 / 60, defaultParams);

      // Act, Assert
      for (let it = 0; it < 10; it += 1) {
        for (const c of contacts) {
          nd.solveContact(D, c);
          assert.ok(c.normalImpulse >= 0, `the impulse of the normal is ${c.normalImpulse}`);
        }
      }
    }));
  });

  it('should always hold the friction inside the cone of Coulomb', () => {
    // The clamp is on the LENGTH of the whole tangent impulse, and not on
    // each tangent one by one. Thus the friction does not depend on the
    // choice of the tangents. In 4 dimensions there are 3 tangents.
    fc.assert(fc.property(anyContact(), (spec) => {
      // Arrange
      const { D, contacts } = makeContacts(spec);
      for (const c of contacts) nd.prepareContact(D, c, 1 / 60, defaultParams);

      // Act, Assert
      for (let it = 0; it < 10; it += 1) {
        for (const c of contacts) {
          nd.solveContact(D, c);
          assert.ok(norm(c.tangentImpulse) <= c.friction * c.normalImpulse + 1e-9,
            `the friction impulse ${norm(c.tangentImpulse)} passes the limit ${c.friction * c.normalImpulse}`);
        }
      }
    }));
  });

  it('should never change the total momentum of the two bodies', () => {
    // The impulse of a contact is internal: the two bodies take the same
    // impulse, with the sign turned around, at the same point. Thus the sum
    // of the momentum of the pair does not change, whatever the solver does.
    fc.assert(fc.property(anyContact(), (spec) => {
      // Arrange
      const { D, a, b, contacts } = makeContacts(spec);
      if (contacts.length === 0) return;
      for (const c of contacts) nd.prepareContact(D, c, 1 / 60, defaultParams);
      const before = totalMomentum(D, [a, b]);

      // Act
      for (let it = 0; it < 10; it += 1) for (const c of contacts) nd.solveContact(D, c);

      // Assert
      const after = totalMomentum(D, [a, b]);
      const scale = 10 * (1 + maxAbs(before.linear) + maxAbs(before.angular));
      assertScaled(after.linear, before.linear, scale, 'the total linear momentum', 1e-9);
      assertScaled(after.angular, before.angular, scale, 'the total angular momentum', 1e-9);
    }));
  });

  it('should always bring the two bodies to the target speed', () => {
    // One contact between two bodies that come together, with no friction.
    // The normal part of the solver is then one linear equation, and the
    // solver must answer it exactly.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRadius(0.5, 1.5), anyRadius(0.5, 1.5), anyUnitVector(n),
      anyPositive(0.5, 4), anyPositive(0.5, 4), anyPositive(1, 8),
    )), ([n, ra, rb, dir, ma, mb, speed]) => {
      // Arrange -- the two balls hold each other a little, and they come
      // together. A distance of exactly the sum of the radii would give a
      // depth of zero, and rounding would then take the contact away.
      const D = dims(n);
      const params = { ...defaultParams, restitutionThreshold: 1e9 };
      const a = new nd.Body(D, {
        shape: nd.HyperSphere(D, ra), mass: ma, friction: 0, restitution: 0,
        velocity: Float64Array.from(dir, (v) => v * speed),
      });
      const b = new nd.Body(D, {
        shape: nd.HyperSphere(D, rb), mass: mb, friction: 0, restitution: 0,
        position: Float64Array.from(dir, (v) => v * (ra + rb) * 0.98),
      });
      const contacts = [];
      nd.collide(D, a, b, contacts, params);
      assert.equal(contacts.length, 1, 'the two balls make one contact');
      nd.prepareContact(D, contacts[0], 1 / 60, params);

      // Act
      for (let it = 0; it < 10; it += 1) nd.solveContact(D, contacts[0]);

      // Assert -- the target holds the correction of the depth. The two
      // bodies must come apart at exactly that speed.
      const u = Float64Array.from(b.pointVelocity(contacts[0].rB),
        (v, i) => v - a.pointVelocity(contacts[0].rA)[i]);
      assertClose(dot(u, contacts[0].normal), contacts[0].target,
        'the speed along the normal after the solver', 1e-6);
    }));
  });
});

describe('warmStart', () => {
  it('should never apply an impulse to a contact that has a gap', () => {
    // A contact with a negative depth only stops the two bodies from coming
    // together. An impulse there would push two bodies apart that do not
    // touch.
    fc.assert(fc.property(anyContact(), (spec) => {
      // Arrange
      const { D, a, b, contacts } = makeContacts(spec);
      for (const c of contacts) {
        nd.prepareContact(D, c, 1 / 60, defaultParams);
        c.normalImpulse = 5;
        c.tangentImpulse.fill(3);
      }
      const before = totalMomentum(D, [a, b]);

      // Act
      for (const c of contacts) nd.warmStart(D, c);

      // Assert
      for (const c of contacts) {
        if (c.depth >= 0) return;
        assert.equal(c.normalImpulse, 0, 'the impulse of a gap starts at zero');
      }
      const after = totalMomentum(D, [a, b]);
      assertScaled(after.linear, before.linear, 10, 'the momentum does not change', 1e-9);
    }));
  });
});

describe('buildContactGraph', () => {
  it('should always give the levels of a stack, from the ground up', () => {
    // The level is the count of the contacts between a body and the nearest
    // static body. `shockPropagation` takes the contacts in that order.
    fc.assert(fc.property(fc.integer({ min: 1, max: 5 }), (count) => {
      // Arrange -- a tower of boxes on the ground.
      const D = dims(3);
      const world = new nd.World({ dimensions: 3 });
      const ground = world.createBody({ shape: nd.HalfSpace(D, [0, 1, 0], 0) });
      const tower = [];
      for (let i = 0; i < count; i += 1) {
        tower.push(world.createBody({
          shape: nd.HyperBox(D, [0.5, 0.5, 0.5]), mass: 1, position: [0, 0.5 + i, 0],
        }));
      }

      // Act
      const contacts = world.narrowPhase(world.broadPhase());
      nd.buildContactGraph(world.bodies, contacts);

      // Assert
      assert.equal(ground.level, 0, 'a static body is the ground of the levels');
      for (let i = 0; i < count; i += 1) {
        assert.equal(tower[i].level, i + 1, `the box ${i} of the tower`);
      }
    }));
  });

  it('should never give two bodies of a contact levels that differ by more than one', () => {
    fc.assert(fc.property(anyContact(), (spec) => {
      // Arrange
      const { D, a, b, contacts } = makeContacts(spec);
      if (contacts.length === 0) return;

      // Act
      nd.buildContactGraph([a, b], contacts);

      // Assert
      assert.ok(Math.abs(a.level - b.level) <= 1,
        `the levels ${a.level} and ${b.level} differ by more than one`);
    }));
  });
});

// Helpers

/**
 * The data of two bodies that touch: two balls, or a ball and a box, with
 * their masses, their speeds and the friction of the pair.
 */
function anyContact() {
  return anyN(2, 4).chain((n) => fc.record({
    n: fc.constant(n),
    kinds: fc.constantFrom(['sphere', 'sphere'], ['box', 'sphere'], ['box', 'box']),
    sizes: fc.tuple(anyRadius(0.4, 1.2), anyRadius(0.4, 1.2)),
    masses: fc.tuple(anyPositive(0.3, 4), anyPositive(0.3, 4)),
    dir: anyUnitVector(n),
    gap: fc.double({ min: 0.5, max: 1.05, noNaN: true }),
    velocities: fc.tuple(anyVector(n, fc.double({ min: -5, max: 5, noNaN: true })),
      anyVector(n, fc.double({ min: -5, max: 5, noNaN: true }))),
    spins: fc.tuple(anyAngularVelocity(n), anyAngularVelocity(n)),
    rotors: fc.tuple(anyRotor(n), anyRotor(n)),
    friction: fc.double({ min: 0, max: 1.5, noNaN: true }),
    restitution: fc.double({ min: 0, max: 1, noNaN: true }),
  }));
}

/** The two bodies of a contact specification, and their contacts. */
function makeContacts(spec) {
  const D = dims(spec.n);
  const shapeOf = (kind, i) => (kind === 'box'
    ? nd.HyperBox(D, new Float64Array(spec.n).fill(spec.sizes[i]))
    : nd.HyperSphere(D, spec.sizes[i]));
  const sa = shapeOf(spec.kinds[0], 0);
  const sb = shapeOf(spec.kinds[1], 1);
  const common = { friction: spec.friction, restitution: spec.restitution };
  const a = new nd.Body(D, {
    ...common,
    shape: sa,
    mass: spec.masses[0],
    rotor: spec.rotors[0],
    velocity: spec.velocities[0],
    angularVelocity: spec.spins[0],
  });
  const b = new nd.Body(D, {
    ...common,
    shape: sb,
    mass: spec.masses[1],
    rotor: spec.rotors[1],
    velocity: spec.velocities[1],
    angularVelocity: spec.spins[1],
    position: Float64Array.from(spec.dir,
      (v) => v * (sa.boundingRadius + sb.boundingRadius) * spec.gap),
  });
  const contacts = [];
  nd.collide(D, a, b, contacts, defaultParams);
  return { D, a, b, contacts };
}

/** The total linear and angular momentum of a list of bodies, about the origin. */
function totalMomentum(D, bodies) {
  const linear = new Float64Array(D.n);
  const angular = new Float64Array(D.k);
  for (const body of bodies) {
    const p = Float64Array.from(body.v, (v) => v * body.mass);
    for (let i = 0; i < D.n; i += 1) linear[i] += p[i];
    const orbit = nd.wedgeVec(D, body.x, p);
    for (let q = 0; q < D.k; q += 1) angular[q] += body.L[q] + orbit[q];
  }
  return { linear, angular };
}
