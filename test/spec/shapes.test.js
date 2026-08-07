/**
 * test/spec/shapes.test.js -- src/nd/body/shapes.js
 *
 * The five shapes. They all obey one interface, thus the collision code can
 * take any of them. These tests hold each shape to the meaning of that
 * interface, and they compare the analytic values with the mesh values.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyHalfExtents, anyRadius, anyUnitVector, anyPositive, anyPlane, anyAngle,
} from '../lib/arbitraries.js';
import { assertClose, assertArrayClose, dot, norm } from '../lib/numeric.js';

const { dims } = nd;

describe('ballVolume', () => {
  it('should give the classical volumes of 2 and 3 dimensions', () => {
    assertClose(nd.ballVolume(2, 1), Math.PI, 'the area of a circle');
    assertClose(nd.ballVolume(3, 1), 4 * Math.PI / 3, 'the volume of a ball');
    assertClose(nd.ballVolume(4, 1), Math.PI ** 2 / 2, 'the volume of a 4-ball');
  });

  it('should always obey the recurrence of the ball volumes', () => {
    // V_n = V_(n-2) * 2 pi / n. The two branches of the function, for an odd
    // n and for an even n, must meet at each step of this recurrence.
    fc.assert(fc.property(anyN(4, 8), (n) => {
      assertClose(nd.ballVolume(n, 1), nd.ballVolume(n - 2, 1) * 2 * Math.PI / n,
        `the volume of a ball of ${n} dimensions`, 1e-12);
    }));
  });

  it('should always grow with the radius to the power n', () => {
    fc.assert(fc.property(fc.tuple(anyN(2, 8), anyPositive(0.1, 5)), ([n, r]) => {
      assertClose(nd.ballVolume(n, r), nd.ballVolume(n, 1) * r ** n, 'the volume of a ball', 1e-12);
    }));
  });
});

describe('HyperSphere', () => {
  it('should give the inertia of a solid ball in 3 dimensions', () => {
    // The classical value is 2 m r^2 / 5.
    const D = dims(3);
    const I = nd.HyperSphere(D, 2).inertia(3);

    // Assert -- the tensor is the same in each plane.
    for (let p = 0; p < D.k; p += 1) {
      assertClose(I[p * D.k + p], 2 * 3 * 4 / 5, `the plane ${p} of the inertia of a ball`);
    }
  });

  it('should always give a support point on its surface', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRadius(), anyUnitVector(n),
    )), ([n, r, dir]) => {
      // Act
      const p = nd.HyperSphere(dims(n), r).support(dir);

      // Assert
      assertClose(norm(p), r, 'the support point of a ball is on its surface', 1e-9);
      assertClose(dot(p, dir), r, 'and it is the farthest point along the direction', 1e-9);
    }));
  });
});

describe('HyperBox', () => {
  it('should throw when the count of half extents is not n', () => {
    assert.throws(() => nd.HyperBox(dims(4), [1, 1, 1]), /half extents/);
  });

  it('should always give the volume of the product of its sides', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyHalfExtents(n))),
      ([n, h]) => {
        assertClose(nd.HyperBox(dims(n), h).volume, h.reduce((a, b) => a * 2 * b, 1),
          'the volume of a box', 1e-12);
      }));
  });
});

describe('HalfSpace', () => {
  it('should always give a normal of the length 1', () => {
    // The solver takes the normal as a unit vector. `HalfSpace` must not
    // depend on the caller for that.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyUnitVector(n), anyPositive(0.1, 100),
    )), ([n, dir, scale]) => {
      // Arrange
      const long = Float64Array.from(dir, (v) => v * scale);

      // Act
      const shape = nd.HalfSpace(dims(n), long, 2);

      // Assert
      assertClose(norm(shape.normal), 1, 'the normal of a half space', 1e-12);
      assertArrayClose(shape.normal, dir, 'the direction does not change', 1e-9);
    }));
  });

  it('should throw when it is asked for a support point', () => {
    assert.throws(() => nd.HalfSpace(dims(3), [0, 1, 0]).support([1, 0, 0]), /support/);
  });
});

describe('ConvexMesh', () => {
  it('should always agree with HyperBox on the mesh of a box', () => {
    // Two paths to one answer: the analytic box and the general mesh. Item 10
    // of the test plan of PART E.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n), anyPositive(0.2, 5),
    )), ([n, h, mass]) => {
      // Arrange
      const D = dims(n);
      const box = nd.HyperBox(D, h);
      const mesh = nd.hyperBoxMesh(D, h);

      // Act
      const convex = nd.ConvexMesh(D, mesh.vertices, mesh.cells);

      // Assert
      assertClose(convex.volume, box.volume, 'the volume of the box as a mesh', 1e-9);
      assertClose(convex.boundingRadius, box.boundingRadius, 'the bounding radius', 1e-9);
      assertArrayClose(convex.inertia(mass), box.inertia(mass), 'the inertia of the box as a mesh', 1e-9);
    }));
  });

  it('should never change its volume when a cell turns around', () => {
    // `orientCells` gives all of the simplices one direction. Thus the
    // direction of the mesh that the caller gives does not change the result.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyHalfExtents(n), fc.array(fc.boolean(), { minLength: 1 }),
    )), ([n, h, flips]) => {
      // Arrange
      const D = dims(n);
      const mesh = nd.hyperBoxMesh(D, h);
      const mixed = Int32Array.from(mesh.cells);
      for (let c = 0; c * n < mixed.length; c += 1) {
        if (flips[c % flips.length]) {
          const t = mixed[c * n];
          mixed[c * n] = mixed[c * n + 1];
          mixed[c * n + 1] = t;
        }
      }

      // Act
      const straight = nd.ConvexMesh(D, mesh.vertices, mesh.cells);
      const turned = nd.ConvexMesh(D, mesh.vertices, mixed);

      // Assert
      assertClose(turned.volume, straight.volume, 'the volume of the mesh', 1e-9);
      assertArrayClose(turned.inertia(1), straight.inertia(1), 'the inertia of the mesh', 1e-9);
    }));
  });
});

describe('Torus', () => {
  it('should give the classical inertia of a ring in 3 dimensions', () => {
    // About the axis of symmetry:  m (R^2 + 3 r^2 / 4)
    // About a diameter:            m (R^2 / 2 + 5 r^2 / 8)
    // The major plane is (x y), thus the axis of symmetry is z, and that is
    // the bivector component 0.
    const D = dims(3);
    const R = 1.2;
    const r = 0.35;

    // Act
    const I = nd.torusInertia(D, R, r, [0, 1], 1);

    // Assert
    assertClose(I[0], R * R + 0.75 * r * r, 'about the axis of symmetry');
    assertClose(I[1 * 3 + 1], 0.5 * R * R + 0.625 * r * r, 'about a diameter');
    assertClose(I[2 * 3 + 2], 0.5 * R * R + 0.625 * r * r, 'about the other diameter');
  });

  it('should give the volume of Pappus', () => {
    // 2 pi R times the volume of the ball of the tube.
    assertClose(nd.torusVolume(3, 2, 0.5), 2 * Math.PI * 2 * Math.PI * 0.25,
      'the volume of a ring of 3 dimensions');
  });

  it('should throw when the major plane is one axis', () => {
    assert.throws(() => nd.Torus(dims(4), 1, 0.2, [2, 2]), /two different axes/);
  });

  it('should always give a mesh whose volume comes near the formula', () => {
    // The mesh is inside the torus, because its vertices are on the surface
    // and its faces are flat. Thus its volume is always a little smaller than
    // the formula, and it comes nearer as the count of the segments grows.
    fc.assert(fc.property(fc.tuple(anyPositive(0.8, 2), anyRadius(0.1, 0.4)), ([R, r]) => {
      // Arrange
      const D = dims(3);
      const shape = nd.Torus(D, R + r, r, [0, 1]);

      // Act
      const mesh = shape.mesh([64, 32]);
      const volume = nd.massProperties(D, mesh.vertices, mesh.cells, 1).volume;

      // Assert
      assert.ok(volume <= shape.volume, 'the mesh is inside the torus');
      assert.ok(volume > shape.volume * 0.97, `the mesh volume is ${volume} and the formula gives ${shape.volume}`);
    }));
  });
});

describe('support', () => {
  it('should never give a point that is behind a point of the shape', () => {
    // The definition of a support point: no point of the shape goes farther
    // along the direction. `boxBox` and the separating axis test of A10 need
    // this to be true for each direction.
    fc.assert(fc.property(anyN(2, 5).chain((n) => anyShapeAndPoint(dims(n))),
      ({ D, shape, dir, inside }) => {
        // Act
        const s = shape.support(dir);

        // Assert
        assert.ok(dot(s, dir) >= dot(inside, dir) - 1e-9 * (1 + norm(inside)),
          `the support point gives ${dot(s, dir)}, and a point of the shape gives ${dot(inside, dir)}`);
      }));
  });

  it('should never give a point outside the bounding radius', () => {
    // The broad phase builds the box of a body from `boundingRadius`. A
    // support point outside it would give a lost collision.
    fc.assert(fc.property(anyN(2, 5).chain((n) => anyShapeAndPoint(dims(n))),
      ({ shape, dir }) => {
        assert.ok(norm(shape.support(dir)) <= shape.boundingRadius * (1 + 1e-9),
          'the support point is inside the bounding radius');
      }));
  });
});

// Helpers

/**
 * A shape, a direction, and a point that is inside that shape. The point
 * comes from the definition of the shape, and not from `support`.
 */
function anyShapeAndPoint(D) {
  const { n } = D;
  const box = fc.tuple(anyHalfExtents(n), anyHalfExtents(n)).map(([h, t]) => ({
    shape: nd.HyperBox(D, h),
    // Each component of the point stays between the two half extents.
    inside: Float64Array.from(h, (v, i) => v * (2 * (t[i] / 3) - 1)),
  }));
  const ball = fc.tuple(anyRadius(), anyUnitVector(n), fc.double({ min: 0, max: 1, noNaN: true }))
    .map(([r, u, t]) => ({
      shape: nd.HyperSphere(D, r),
      inside: Float64Array.from(u, (v) => v * r * t),
    }));
  const ring = fc.tuple(anyPositive(0.8, 2), anyRadius(0.1, 0.4), anyPlane(n), anyAngle(),
    anyUnitVector(n), fc.double({ min: 0, max: 1, noNaN: true }))
    .map(([R, r, plane, angle, u, t]) => {
      // A point of the major circle, and a step of the length r or less. The
      // distance from that point to the circle is then r or less.
      const inside = Float64Array.from(u, (v) => v * r * t);
      inside[plane[0]] += R * Math.cos(angle);
      inside[plane[1]] += R * Math.sin(angle);
      return { shape: nd.Torus(D, R, r, plane), inside };
    });
  return fc.tuple(fc.oneof(box, ball, ring), anyUnitVector(n))
    .map(([part, dir]) => ({ D, dir, ...part }));
}
