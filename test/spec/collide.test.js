/**
 * test/spec/collide.test.js -- src/nd/detect/collide.js
 *
 * The narrow phase. The properties here hold for EACH pair of shape types.
 * That is the reason to write them: the code turns the normal around with a
 * `flip` argument in four places, and PART F of ND-PHYSICS.md says that the
 * boundary conditions of a collision are more numerous in 4 dimensions than
 * in 3.
 *
 * An arbitrary here gives plain data (the count of dimensions, the sizes, a
 * direction), and the test builds the bodies. Thus a counterexample is short
 * and a person can read it.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyRadius, anyHalfExtents, anyUnitVector, anyRotor, anyPositive, anyPosition,
} from '../lib/arbitraries.js';
import { assertClose, assertArrayClose, dot, norm } from '../lib/numeric.js';

const { dims } = nd;

/** The params of a world with the default margin. */
const PARAMS = { contactMargin: 0.02, maxContacts: 0 };

/** The pairs of shape types that `collide` has a test for. */
const PAIRS = [
  ['halfspace', 'sphere'], ['halfspace', 'box'], ['halfspace', 'convex'], ['halfspace', 'torus'],
  ['sphere', 'sphere'], ['sphere', 'box'], ['sphere', 'convex'], ['sphere', 'torus'],
  ['box', 'box'],
];

describe('collide', () => {
  it('should give the depth of two balls that meet', () => {
    // Arrange
    const D = dims(4);
    const a = ball(D, 1, [0, 0, 0, 0]);
    const b = ball(D, 0.5, [1.2, 0, 0, 0]);

    // Act
    const out = [];
    nd.collide(D, a, b, out, PARAMS);

    // Assert
    assert.equal(out.length, 1, 'two balls give one contact');
    assertClose(out[0].depth, 1 + 0.5 - 1.2, 'the depth is the sum of the radii less the distance');
    assertArrayClose(out[0].normal, [1, 0, 0, 0], 'the normal goes from a to b');
  });

  it('should always give a normal of the length 1', () => {
    // The solver divides by the effective mass along the normal. A normal
    // that is not of the length 1 would change the mass of the contact.
    fc.assert(fc.property(anyPairNearby(), (spec) => {
      // Arrange
      const { D, a, b } = makePair(spec);

      // Act
      const out = [];
      nd.collide(D, a, b, out, PARAMS);

      // Assert
      for (const c of out) assertClose(norm(c.normal), 1, 'the length of the normal', 1e-9);
    }));
  });

  it('should always point the normal from a to b', () => {
    // Move `b` along the normal, by the depth and a little more. The two
    // bodies must then be apart. This one property tests the direction of the
    // normal, the size of the depth and each `flip` argument together.
    //
    // The pair must be convex. A torus has a hole: a ball that is larger than
    // the hole touches the tube at each place inside it, thus no move along a
    // normal can make the two apart. That is the shape, and not a defect.
    fc.assert(fc.property(anyPairNearby(), (spec) => {
      // Arrange
      if (spec.kinds.includes('torus')) return;
      const { D, a, b } = makePair(spec);
      const out = [];
      nd.collide(D, a, b, out, PARAMS);
      if (out.length === 0) return;
      const deepest = out.reduce((p, q) => (q.depth > p.depth ? q : p));
      if (deepest.depth <= 0) return;
      const step = deepest.depth + 0.05;

      // Act
      for (let i = 0; i < D.n; i += 1) b.x[i] += deepest.normal[i] * step;
      const after = [];
      nd.collide(D, a, b, after, PARAMS);

      // Assert
      for (const c of after) {
        assert.ok(c.depth <= 1e-4, `after the move of ${step} the depth is still ${c.depth}`);
      }
    }));
  });

  it('should always agree when the two bodies change place', () => {
    // collide(a, b) and collide(b, a) are the same contact, seen from the two
    // sides. The order of a pair comes from a sort in the broad phase, thus
    // the answer must not depend on it.
    //
    // The test looks at the count and at the deepest contact.
    //
    // TWO BOXES ARE DIFFERENT. `boxBoxAxis` keeps the first axis of the
    // smallest overlap. Two boxes can give the same overlap on more than one
    // axis -- for example a box of [0.3, 0.3, 0.5, 0.3] against a box of
    // [0.3, 0.5, 0.3, 0.3] -- and the first of them depends on the order of
    // the pair. The depth is the same, but the normal and the points are not.
    // Each engine must choose one axis of a tie, thus the test asks only for
    // the depth there.
    fc.assert(fc.property(anyPairNearby(), (spec) => {
      // Arrange
      const { D, a, b } = makePair(spec);

      // Act
      const forward = [];
      const backward = [];
      nd.collide(D, a, b, forward, PARAMS);
      nd.collide(D, b, a, backward, PARAMS);

      // Assert
      assert.equal(backward.length, forward.length, 'the count of the contacts');
      if (forward.length === 0) return;
      const deepest = (list) => list.reduce((p, q) => (q.depth > p.depth ? q : p));
      const one = deepest(forward);
      const other = deepest(backward);
      assertClose(other.depth, one.depth, 'the depth of the deepest contact', 1e-9);
      if (spec.kinds[0] === 'box' && spec.kinds[1] === 'box') return;
      assertArrayClose(other.normal, Float64Array.from(one.normal, (v) => -v),
        'the normal points the other way', 1e-9);
    }));
  });

  it('should never give a contact to two bodies that are far apart', () => {
    // The two bounding balls do not meet, thus the two bodies do not meet.
    fc.assert(fc.property(fc.tuple(anyPairNearby(), anyPositive(1.5, 20)), ([spec, extra]) => {
      // Arrange -- a half space has no bounding radius, thus leave it out.
      if (spec.kinds[0] === 'halfspace') return;
      const { D, a, b } = makePair(spec);
      const far = (a.shape.boundingRadius + b.shape.boundingRadius) * extra + 1;
      for (let i = 0; i < D.n; i += 1) b.x[i] = a.x[i] + spec.dir[i] * far;

      // Act
      const out = [];
      nd.collide(D, a, b, out, PARAMS);

      // Assert
      assert.equal(out.length, 0, 'two bodies that are far apart do not touch');
    }));
  });

  it('should always put the point of contact inside the two bodies', () => {
    // The solver builds the offsets rA and rB from this point. A point far
    // from the two bodies would give a torque that is not physical.
    //
    // The point must be inside the box of ONE of the two bodies. It is not
    // always inside the two: `boxBox` keeps a corner of `b` that is near the
    // face of `a`, and that corner can be outside `a` on another axis.
    fc.assert(fc.property(anyPairNearby(), (spec) => {
      // Arrange
      const { D, a, b } = makePair(spec);

      // Act
      const out = [];
      nd.collide(D, a, b, out, PARAMS);

      // Assert
      const boxA = a.aabb(PARAMS.contactMargin + 1e-6);
      const boxB = b.aabb(PARAMS.contactMargin + 1e-6);
      const holds = (box) => {
        for (let i = 0; i < D.n; i += 1) {
          if (c.point[i] < box.min[i] || c.point[i] > box.max[i]) return false;
        }
        return true;
      };
      let c = null;
      for (c of out) {
        assert.ok(holds(boxA) || holds(boxB),
          `the point ${c.point} is outside the box of a and the box of b`);
      }
    }));
  });

  it('should always give the exact depth of two balls', () => {
    // An oracle: the depth of two balls is the sum of the radii less the
    // distance between the two centers.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRadius(), anyRadius(), anyUnitVector(n),
      fc.double({ min: 0.1, max: 3, noNaN: true }),
    )), ([n, ra, rb, dir, gap]) => {
      // Arrange
      const D = dims(n);
      const distance = (ra + rb) * gap;
      const a = ball(D, ra, new Float64Array(n));
      const b = ball(D, rb, Float64Array.from(dir, (v) => v * distance));

      // Act
      const out = [];
      nd.collide(D, a, b, out, PARAMS);

      // Assert
      if (ra + rb - distance < 0) {
        assert.equal(out.length, 0, 'two balls that do not meet give no contact');
      } else {
        assertClose(out[0].depth, ra + rb - distance, 'the depth of two balls', 1e-9);
        assertArrayClose(out[0].normal, dir, 'the normal is the line of the two centers', 1e-9);
      }
    }));
  });

  it('should always give the exact depth of a ball on a half space', () => {
    // Another oracle: the depth is the radius less the distance from the
    // center to the plane.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyUnitVector(n), anyRadius(), anyPosition(n, 3),
    )), ([n, normal, radius, place]) => {
      // Arrange
      const D = dims(n);
      const ground = new nd.Body(D, { shape: nd.HalfSpace(D, normal, 0) });
      const b = ball(D, radius, place);

      // Act
      const out = [];
      nd.collide(D, ground, b, out, PARAMS);

      // Assert
      const depth = radius - dot(normal, place);
      if (depth < -PARAMS.contactMargin) {
        assert.equal(out.length, 0, 'a ball above the ground gives no contact');
      } else {
        assertClose(out[0].depth, depth, 'the depth of a ball on a half space', 1e-9);
        assertArrayClose(out[0].normal, normal, 'the normal goes from the ground to the ball', 1e-9);
      }
    }));
  });

  it('should always give a box and its own mesh the same contact', () => {
    // The same box, as an analytic box and as a convex mesh, against the same
    // ball. `boxSphere` and `convexSphere` (thus also `nearestOnSimplex`) must
    // give one answer. Item 10 of the test plan of PART E.
    //
    // The center of the ball stays outside the box: it is farther from the
    // center than the bounding radius. Inside the box the nearest point is
    // not defined, and the two functions then take different ways out.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n, 0.3, 1.5), anyRadius(0.3, 1),
      anyUnitVector(n), fc.double({ min: 0, max: 1.2, noNaN: true }), anyRotor(n),
    )), ([n, h, radius, dir, gap, R]) => {
      // Arrange
      const D = dims(n);
      const mesh = nd.hyperBoxMesh(D, h);
      const away = Float64Array.from(dir, (v) => v * (norm(h) + radius * gap));
      const box = new nd.Body(D, { shape: nd.HyperBox(D, h), mass: 1, rotor: R });
      const convex = new nd.Body(D, {
        shape: nd.ConvexMesh(D, mesh.vertices, mesh.cells), mass: 1, rotor: R,
      });

      // Act
      const fromBox = [];
      const fromMesh = [];
      nd.collide(D, box, ball(D, radius, away), fromBox, PARAMS);
      nd.collide(D, convex, ball(D, radius, away), fromMesh, PARAMS);

      // Assert
      assert.equal(fromMesh.length, fromBox.length, 'the two tests find the same count of contacts');
      for (let i = 0; i < fromBox.length; i += 1) {
        assertClose(fromMesh[i].depth, fromBox[i].depth, 'the depth', 1e-6);
        assertArrayClose(fromMesh[i].normal, fromBox[i].normal, 'the normal', 1e-6);
      }
    }));
  });
});

describe('boxBoxAxis', () => {
  it('should always give a unit axis and an overlap that is not negative', () => {
    fc.assert(fc.property(anyBoxPairSpec(), (spec) => {
      // Arrange
      const { D, a, b } = makeBoxPair(spec);

      // Act
      const res = nd.boxBoxAxis(D, a, b);

      // Assert
      if (res === null) return;
      assertClose(norm(res.axis), 1, 'the length of the axis', 1e-9);
      assert.ok(res.overlap >= 0, `the overlap is ${res.overlap}`);
    }));
  });

  it('should always give the true overlap of the two boxes along its axis', () => {
    // The separating axis theorem: the overlap along the axis is the sum of
    // the two extents less the distance between the two centers.
    fc.assert(fc.property(anyBoxPairSpec(), (spec) => {
      // Arrange
      const { D, a, b } = makeBoxPair(spec);

      // Act
      const res = nd.boxBoxAxis(D, a, b);
      if (res === null) return;

      // Assert
      const delta = Float64Array.from(b.x, (v, i) => v - a.x[i]);
      assertClose(res.overlap,
        extentAlong(D, a, res.axis) + extentAlong(D, b, res.axis) - Math.abs(dot(res.axis, delta)),
        'the overlap along the axis', 1e-6);
    }));
  });

  it('should never find an axis for two boxes that are far apart', () => {
    fc.assert(fc.property(fc.tuple(anyBoxPairSpec(), anyPositive(1.1, 10)), ([spec, extra]) => {
      // Arrange
      const { D, a, b } = makeBoxPair(spec);
      const far = (a.shape.boundingRadius + b.shape.boundingRadius) * extra + 0.1;
      for (let i = 0; i < D.n; i += 1) b.x[i] = a.x[i] + spec.dir[i] * far;

      // Act, Assert
      assert.equal(nd.boxBoxAxis(D, a, b), null, 'two boxes that are far apart are apart');
    }));
  });
});

describe('tangentBasis', () => {
  it('should always give n - 1 unit vectors that are normal to the normal', () => {
    // The friction acts in these directions, and the solver holds the length
    // of the whole tangent impulse. That is the cone of Coulomb only when the
    // directions are orthogonal and of the length 1.
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(
      fc.constant(n), anyUnitVector(n),
    )), ([n, normal]) => {
      // Arrange
      const D = dims(n);

      // Act
      const T = nd.tangentBasis(D, normal);

      // Assert
      assert.equal(T.length, n - 1, 'the count of the tangents');
      for (let i = 0; i < T.length; i += 1) {
        assertClose(norm(T[i]), 1, `the length of the tangent ${i}`, 1e-9);
        assertClose(dot(T[i], normal), 0, `the tangent ${i} is normal to the normal`, 1e-9);
        for (let j = i + 1; j < T.length; j += 1) {
          assertClose(dot(T[i], T[j]), 0, `the tangents ${i} and ${j} are orthogonal`, 1e-9);
        }
      }
    }));
  });
});

describe('World.broadPhase', () => {
  it('should never lose a pair that touches', () => {
    // The broad phase gives the pairs to the narrow phase. A pair that it
    // drops can never make a contact. Thus each pair that `collide` answers
    // must be in its list.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n),
      fc.array(fc.tuple(fc.boolean(), anyHalfExtents(n, 0.3, 1), anyRadius(0.3, 1),
        anyPosition(n, 2), anyRotor(n)), { minLength: 2, maxLength: 6 }),
    )), ([n, specs]) => {
      // Arrange
      const D = dims(n);
      const world = new nd.World({ dimensions: n });
      const bodies = specs.map(([isBox, h, r, place, R]) => world.addBody(new nd.Body(D, {
        shape: isBox ? nd.HyperBox(D, h) : nd.HyperSphere(D, r),
        mass: 1,
        position: place,
        rotor: R,
      })));

      // Act
      const found = new Set();
      for (const [a, b] of world.broadPhase()) found.add(key(a, b));

      // Assert
      for (let i = 0; i < bodies.length; i += 1) {
        for (let j = i + 1; j < bodies.length; j += 1) {
          const out = [];
          nd.collide(D, bodies[i], bodies[j], out, world.params);
          if (out.some((c) => c.depth > 0)) {
            assert.ok(found.has(key(bodies[i], bodies[j])),
              `the broad phase lost the pair ${bodies[i].id} and ${bodies[j].id}`);
          }
        }
      }
    }));
  });
});

// Helpers

/** A body with a ball shape at the place `x`. */
function ball(D, radius, x) {
  return new nd.Body(D, { shape: nd.HyperSphere(D, radius), mass: 1, position: x });
}

/** The name of a pair of bodies, in the order of the two ids. */
function key(a, b) {
  return a.id < b.id ? `${a.id}:${b.id}` : `${b.id}:${a.id}`;
}

/** The extent of a box along a direction, from its center. */
function extentAlong(D, body, axis) {
  let s = 0;
  for (let i = 0; i < D.n; i += 1) {
    let t = 0;
    for (let j = 0; j < D.n; j += 1) t += body.Rm[j * D.n + i] * axis[j];
    s += Math.abs(t) * body.shape.halfExtents[i];
  }
  return s;
}

/**
 * The data of two bodies that are near each other: the count of dimensions,
 * the two shape types, the sizes, the direction from the first to the second
 * and the part of the reach that the distance takes.
 */
function anyPairNearby() {
  return anyN(2, 4).chain((n) => fc.record({
    n: fc.constant(n),
    kinds: fc.constantFrom(...PAIRS),
    halfExtents: fc.tuple(anyHalfExtents(n, 0.3, 1.5), anyHalfExtents(n, 0.3, 1.5)),
    radii: fc.tuple(anyRadius(0.3, 1.5), anyRadius(0.3, 1.5)),
    torus: fc.tuple(anyPositive(0.8, 1.5), anyRadius(0.2, 0.4)),
    dir: anyUnitVector(n),
    gap: fc.double({ min: 0.2, max: 1.3, noNaN: true }),
    rotors: fc.tuple(anyRotor(n), anyRotor(n)),
  }));
}

/** The two bodies of a pair specification. */
function makePair(spec) {
  const D = dims(spec.n);
  const shapeOf = (kind, i) => {
    if (kind === 'halfspace') return nd.HalfSpace(D, spec.dir, 0);
    if (kind === 'box') return nd.HyperBox(D, spec.halfExtents[i]);
    if (kind === 'sphere') return nd.HyperSphere(D, spec.radii[i]);
    if (kind === 'torus') return nd.Torus(D, spec.torus[0], spec.torus[1], [0, 1]);
    const m = nd.hyperBoxMesh(D, spec.halfExtents[i]);
    return nd.ConvexMesh(D, m.vertices, m.cells);
  };
  const sa = shapeOf(spec.kinds[0], 0);
  const sb = shapeOf(spec.kinds[1], 1);
  const a = new nd.Body(D, {
    shape: sa,
    mass: spec.kinds[0] === 'halfspace' ? 0 : 1,
    rotor: spec.kinds[0] === 'halfspace' ? undefined : spec.rotors[0],
  });
  const reach = spec.kinds[0] === 'halfspace'
    ? sb.boundingRadius
    : sa.boundingRadius + sb.boundingRadius;
  const b = new nd.Body(D, {
    shape: sb,
    mass: 1,
    rotor: spec.rotors[1],
    position: Float64Array.from(spec.dir, (v) => v * reach * spec.gap),
  });
  return { D, a, b };
}

/** The data of two boxes that are near each other. */
function anyBoxPairSpec() {
  return anyN(2, 4).chain((n) => fc.record({
    n: fc.constant(n),
    halfExtents: fc.tuple(anyHalfExtents(n, 0.3, 1.5), anyHalfExtents(n, 0.3, 1.5)),
    dir: anyUnitVector(n),
    gap: fc.double({ min: 0.2, max: 1.4, noNaN: true }),
    rotors: fc.tuple(anyRotor(n), anyRotor(n)),
  }));
}

/** The two bodies of a box pair specification. */
function makeBoxPair(spec) {
  const D = dims(spec.n);
  const [ha, hb] = spec.halfExtents;
  const a = new nd.Body(D, { shape: nd.HyperBox(D, ha), mass: 1, rotor: spec.rotors[0] });
  const b = new nd.Body(D, {
    shape: nd.HyperBox(D, hb),
    mass: 1,
    rotor: spec.rotors[1],
    position: Float64Array.from(spec.dir, (v) => v * (norm(ha) + norm(hb)) * spec.gap),
  });
  return { D, a, b };
}
