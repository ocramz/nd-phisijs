/**
 * test/spec/slice4.test.js -- src/slice/slice4.js
 *
 * The cut of a 4D body with a hyperplane of 3 dimensions. A screen shows 3
 * dimensions, thus the plugin draws this cut. See ND-PHYSICS.md, D1.
 *
 * The cut of the surface of a tesseract is the surface of a 3D box. That
 * gives an oracle: the area of the triangles and the extents of the points
 * are the area and the extents of that box.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd, slice } from '../lib/load.js';
import { anyHalfExtents, anyPositive, anyAngle, anyRadius } from '../lib/arbitraries.js';
import { assertClose, soupArea, soupExtent } from '../lib/numeric.js';

const D = nd.dims(4);

/** The identity, as a rotation matrix of 4 by 4. */
const NO_TURN = Float64Array.from([1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]);

describe('sliceTetrahedra', () => {
  it('should cut a tesseract into a box of 3 dimensions', () => {
    // Arrange -- a tesseract of the half extents 1, 2, 3, 4, cut at w = 0.
    const mesh = nd.hyperBoxMesh(D, Float64Array.from([1, 2, 3, 4]));
    const world = slice.transformVertices(mesh.vertices, NO_TURN, new Float64Array(4));

    // Act
    const out = slice.sliceTetrahedra(world, mesh.cells, 0);

    // Assert -- the cut is a box of the half extents 1, 2 and 3. Its surface
    // is 8 (h0 h1 + h0 h2 + h1 h2).
    assertClose(soupArea(out.positions), 8 * (1 * 2 + 1 * 3 + 2 * 3), 'the area of the cut', 1e-9);
    for (const [axis, h] of [[0, 1], [1, 2], [2, 3]]) {
      const e = soupExtent(out.positions, axis);
      assertClose(e.min, -h, `the low side of the axis ${axis}`, 1e-9);
      assertClose(e.max, h, `the high side of the axis ${axis}`, 1e-9);
    }
  });

  it('should always give the surface of the box of the cut', () => {
    // The same oracle, for any half extents and any place of the cut inside
    // the tesseract.
    fc.assert(fc.property(anyHalfExtents(4, 0.2, 3).chain((h) => fc.tuple(
      fc.constant(h), fc.double({ min: -0.95, max: 0.95, noNaN: true }),
    )), ([h, part]) => {
      // Arrange
      const mesh = nd.hyperBoxMesh(D, h);
      const world = slice.transformVertices(mesh.vertices, NO_TURN, new Float64Array(4));

      // Act
      const out = slice.sliceTetrahedra(world, mesh.cells, h[3] * part);

      // Assert
      const area = 8 * (h[0] * h[1] + h[0] * h[2] + h[1] * h[2]);
      assertClose(soupArea(out.positions), area, 'the area of the cut', 1e-6);
      for (let axis = 0; axis < 3; axis += 1) {
        const e = soupExtent(out.positions, axis);
        assertClose(e.max - e.min, 2 * h[axis], `the size of the axis ${axis}`, 1e-6);
      }
    }));
  });

  it('should never give a triangle when the cut misses the body', () => {
    fc.assert(fc.property(anyHalfExtents(4, 0.2, 3).chain((h) => fc.tuple(
      fc.constant(h), anyPositive(1.05, 10),
    )), ([h, part]) => {
      // Arrange
      const mesh = nd.hyperBoxMesh(D, h);
      const world = slice.transformVertices(mesh.vertices, NO_TURN, new Float64Array(4));

      // Act
      const out = slice.sliceTetrahedra(world, mesh.cells, h[3] * part);

      // Assert
      assert.equal(out.positions.length, 0, 'a cut outside the body gives nothing');
    }));
  });

  it('should always move the cut with the body', () => {
    // A move of the body in x, y and z moves the triangles by the same
    // quantity. A move in w of `d`, with the cut at `value + d`, gives the
    // same triangles again.
    fc.assert(fc.property(fc.tuple(anyHalfExtents(4, 0.2, 2),
      anyPositive(-3, 3), anyPositive(-3, 3)), ([h, step, deep]) => {
      // Arrange
      const mesh = nd.hyperBoxMesh(D, h);
      const place = Float64Array.from([step, 2 * step, -step, deep]);
      const still = slice.transformVertices(mesh.vertices, NO_TURN, new Float64Array(4));
      const moved = slice.transformVertices(mesh.vertices, NO_TURN, place);

      // Act
      const first = slice.sliceTetrahedra(still, mesh.cells, 0);
      const second = slice.sliceTetrahedra(moved, mesh.cells, deep);

      // Assert
      assert.equal(second.positions.length, first.positions.length, 'the count of the points');
      for (let i = 0; i < first.positions.length; i += 3) {
        assertClose(second.positions[i], first.positions[i] + place[0], 'the axis x', 1e-9);
        assertClose(second.positions[i + 1], first.positions[i + 1] + place[1], 'the axis y', 1e-9);
        assertClose(second.positions[i + 2], first.positions[i + 2] + place[2], 'the axis z', 1e-9);
      }
    }));
  });

  it('should always turn the cut with a turn that holds w still', () => {
    // A turn in the plane (x y) does not change the axis w. Thus the cut of
    // the body that turns is the cut of the body, turned.
    fc.assert(fc.property(fc.tuple(anyHalfExtents(4, 0.2, 2), anyAngle()), ([h, angle]) => {
      // Arrange
      const mesh = nd.hyperBoxMesh(D, h);
      const R = nd.rotor.rotorFromPlane(D, 0, 1, angle);
      const M = nd.rotor.rotorMatrix(D, R);
      const still = slice.transformVertices(mesh.vertices, NO_TURN, new Float64Array(4));
      const turned = slice.transformVertices(mesh.vertices, M, new Float64Array(4));

      // Act
      const first = slice.sliceTetrahedra(still, mesh.cells, 0);
      const second = slice.sliceTetrahedra(turned, mesh.cells, 0);

      // Assert -- the area does not change with a turn, and the points of
      // the second cut are the points of the first, turned.
      assertClose(soupArea(second.positions), soupArea(first.positions), 'the area of the cut', 1e-9);
      const c = Math.cos(angle);
      const s = Math.sin(angle);
      for (let i = 0; i < first.positions.length; i += 3) {
        const x = first.positions[i];
        const y = first.positions[i + 1];
        assertClose(second.positions[i], c * x - s * y, 'the axis x of the turn', 1e-9);
        assertClose(second.positions[i + 1], s * x + c * y, 'the axis y of the turn', 1e-9);
        assertClose(second.positions[i + 2], first.positions[i + 2], 'the axis z', 1e-9);
      }
    }));
  });

  it('should always give a normal of the length 1 to each point', () => {
    fc.assert(fc.property(anyHalfExtents(4, 0.2, 3), (h) => {
      // Arrange
      const mesh = nd.hyperBoxMesh(D, h);
      const world = slice.transformVertices(mesh.vertices, NO_TURN, new Float64Array(4));

      // Act
      const out = slice.sliceTetrahedra(world, mesh.cells, 0);

      // Assert
      assert.equal(out.normals.length, out.positions.length, 'a normal for each point');
      for (let i = 0; i < out.normals.length; i += 3) {
        assertClose(Math.hypot(out.normals[i], out.normals[i + 1], out.normals[i + 2]), 1,
          'the length of a normal', 1e-9);
      }
    }));
  });
});

describe('sliceHyperSphereRadius', () => {
  it('should give the radius of the cut of a 4-ball', () => {
    // The cut of a ball of the radius 5 at the distance 3 is a ball of the
    // radius 4.
    assertClose(slice.sliceHyperSphereRadius(5, 0, 3), 4, 'the radius of the cut');
  });

  it('should always obey the theorem of Pythagoras', () => {
    fc.assert(fc.property(fc.tuple(anyRadius(0.2, 5), anyPositive(-5, 5), anyPositive(-5, 5)),
      ([radius, centerW, value]) => {
        // Act
        const cut = slice.sliceHyperSphereRadius(radius, centerW, value);

        // Assert
        const distance = value - centerW;
        if (Math.abs(distance) >= radius) {
          assert.equal(cut, -1, 'a cut outside the ball gives -1');
        } else {
          assertClose(cut * cut + distance * distance, radius * radius,
            'the radius of the cut and the distance', 1e-9);
        }
      }));
  });
});

describe('sliceTorusMinorRadius', () => {
  it('should always obey the theorem of Pythagoras', () => {
    // The tube of a 4D torus is a 2-sphere, thus its cut is a circle.
    fc.assert(fc.property(fc.tuple(anyRadius(0.1, 2), anyPositive(-3, 3), anyPositive(-3, 3)),
      ([minor, centerW, value]) => {
        // Act
        const cut = slice.sliceTorusMinorRadius(minor, centerW, value);

        // Assert
        const distance = value - centerW;
        if (Math.abs(distance) >= minor) {
          assert.equal(cut, -1, 'a cut outside the tube gives -1');
        } else {
          assertClose(cut * cut + distance * distance, minor * minor,
            'the minor radius of the cut', 1e-9);
        }
      }));
  });
});

describe('transformVertices', () => {
  it('should always give the same points as the body does', () => {
    // The slice moves the vertices with `Rm x + position`. `Body.localToWorld`
    // does the same thing for one point. The two must agree, or the cut would
    // not be at the place of the body.
    fc.assert(fc.property(fc.tuple(anyHalfExtents(4, 0.2, 2), anyAngle(), anyPositive(-3, 3)),
      ([h, angle, step]) => {
        // Arrange
        const mesh = nd.hyperBoxMesh(D, h);
        const R = nd.rotor.rotorFromPlane(D, 1, 3, angle);
        const place = Float64Array.from([step, -step, 2 * step, step / 2]);
        const body = new nd.Body(D, {
          shape: nd.HyperBox(D, h), mass: 1, rotor: R, position: place,
        });

        // Act
        const moved = slice.transformVertices(mesh.vertices, body.Rm, body.x);

        // Assert
        const count = mesh.vertices.length / 4;
        for (let m = 0; m < count; m += 1) {
          const one = body.localToWorld(mesh.vertices.subarray(m * 4, m * 4 + 4),
            new Float64Array(4));
          for (let i = 0; i < 4; i += 1) {
            assertClose(moved[m * 4 + i], one[i], `the vertex ${m}, the axis ${i}`, 1e-9);
          }
        }
      }));
  });
});
