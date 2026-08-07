/**
 * test/spec/massprops.test.js -- src/nd/body/massprops.js
 *
 * The volume, the center of mass and the inertia tensor of a mesh. Item 4 of
 * the test plan of PART E asks for the inertia of a 4D hypercuboid from the
 * mesh, against the analytic value. These tests do that for `n` of 2 to 5,
 * and they add the laws that hold for each mesh.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import { anyN, anyHalfExtents, anyRotor, anyPosition, anyPositive } from '../lib/arbitraries.js';
import {
  assertClose, assertArrayClose, assertScaled, maxAbs, isPositiveDefinite, symmetryError,
  mulMat, transpose,
} from '../lib/numeric.js';

const { dims, rotor } = nd;

/** An `n`, its half extents, and the mesh of that box. */
const anyBoxMesh = (max = 5) => anyN(2, max).chain((n) => fc.tuple(
  fc.constant(dims(n)), anyHalfExtents(n),
)).map(([D, h]) => ({ D, h, mesh: nd.hyperBoxMesh(D, h), volume: h.reduce((a, b) => a * 2 * b, 1) }));

describe('massProperties', () => {
  it('should give the volume and the center of a box of 4 dimensions', () => {
    // Arrange
    const D = dims(4);
    const h = Float64Array.from([0.5, 1, 1.5, 2]);
    const mesh = nd.hyperBoxMesh(D, h);

    // Act
    const p = nd.massProperties(D, mesh.vertices, mesh.cells, 1);

    // Assert
    assertClose(p.volume, 1 * 2 * 3 * 4, 'the volume of a tesseract is the product of the sides');
    assertArrayClose(p.center, [0, 0, 0, 0], 'the center of mass is at the origin');
  });

  it('should always agree with the analytic inertia of a box', () => {
    // Item 4 of the test plan. The mesh method walks over 2n faces and (n-1)!
    // simplices for each face; the analytic value is m (h_i^2 + h_j^2) / 3.
    fc.assert(fc.property(fc.tuple(anyBoxMesh(), anyPositive(0.2, 5)), ([box, mass]) => {
      // Arrange
      const { D, h, mesh, volume } = box;

      // Act
      const p = nd.massProperties(D, mesh.vertices, mesh.cells, mass / volume);

      // Assert
      assertClose(p.volume, volume, 'the volume of the mesh', 1e-9);
      assertClose(p.mass, mass, 'the mass of the mesh', 1e-9);
      assertArrayClose(p.inertia, nd.hyperBoxInertia(D, h, mass), 'the inertia of the mesh', 1e-9);
    }));
  });

  it('should never change the inertia when the mesh moves', () => {
    // The inertia is about the center of mass. Thus a move of the whole mesh
    // moves the center and leaves the inertia. This is the parallel axis law.
    //
    // The offset stays inside +-5. `massProperties` builds the covariance
    // about the origin, then it takes `volume * center * center` away. That
    // subtraction loses accuracy with the square of the offset. A mesh in a
    // body frame is always near its center of mass, thus the library does not
    // see a large offset.
    fc.assert(fc.property(anyBoxMesh().chain((box) => fc.tuple(
      fc.constant(box), anyPosition(box.D.n, 5),
    )), ([box, offset]) => {
      // Arrange
      const { D, mesh } = box;
      const moved = Float64Array.from(mesh.vertices);
      for (let m = 0; m < moved.length; m += D.n) {
        for (let i = 0; i < D.n; i += 1) moved[m + i] += offset[i];
      }

      // Act
      const before = nd.massProperties(D, mesh.vertices, mesh.cells, 1);
      const after = nd.massProperties(D, moved, mesh.cells, 1);

      // Assert
      assertClose(after.volume, before.volume, 'the volume does not change', 1e-9);
      assertArrayClose(after.center, offset, 'the center of mass moves with the mesh', 1e-9);
      assertScaled(after.inertia, before.inertia, maxAbs(before.inertia),
        'the inertia does not change', 1e-9);
    }));
  });

  it('should always change the inertia by [R]2 I [R]2^T when the mesh turns', () => {
    // This connects A8 (the mass properties) to A4 and A6 (the change of
    // frame that `updateDerived` makes in each step). If the two did not
    // agree, a body that turns would show the wrong inertia.
    fc.assert(fc.property(anyBoxMesh().chain((box) => fc.tuple(
      fc.constant(box), anyRotor(box.D.n),
    )), ([box, R]) => {
      // Arrange
      const { D, mesh } = box;
      const { n, k } = D;
      const M = rotor.rotorMatrix(D, R);
      const turned = new Float64Array(mesh.vertices.length);
      for (let m = 0; m < turned.length; m += n) {
        for (let i = 0; i < n; i += 1) {
          let s = 0;
          for (let j = 0; j < n; j += 1) s += M[i * n + j] * mesh.vertices[m + j];
          turned[m + i] = s;
        }
      }

      // Act
      const direct = nd.massProperties(D, turned, mesh.cells, 1).inertia;
      const body = nd.massProperties(D, mesh.vertices, mesh.cells, 1).inertia;
      const R2 = rotor.rotorBivectorMatrix(D, R);
      const changed = mulMat(mulMat(R2, body, k, k, k), transpose(R2, k, k), k, k, k);

      // Assert
      assertScaled(direct, changed, maxAbs(direct), 'the inertia of a mesh that turns', 1e-9);
    }));
  });

  it('should always scale the volume by s^n and the inertia by s^(n+2)', () => {
    // Dimensional analysis. A volume is a length to the power n; an inertia
    // is a mass times a length squared, and the mass follows the volume.
    fc.assert(fc.property(anyBoxMesh().chain((box) => fc.tuple(
      fc.constant(box), anyPositive(0.3, 3),
    )), ([box, s]) => {
      // Arrange
      const { D, mesh } = box;
      const larger = Float64Array.from(mesh.vertices, (v) => v * s);

      // Act
      const before = nd.massProperties(D, mesh.vertices, mesh.cells, 1);
      const after = nd.massProperties(D, larger, mesh.cells, 1);

      // Assert
      assertClose(after.volume, before.volume * s ** D.n, 'the volume of the larger mesh', 1e-9);
      assertArrayClose(after.inertia, Float64Array.from(before.inertia, (v) => v * s ** (D.n + 2)),
        'the inertia of the larger mesh', 1e-9);
    }));
  });

  it('should always be linear in the density', () => {
    fc.assert(fc.property(anyBoxMesh().chain((box) => fc.tuple(
      fc.constant(box), anyPositive(0.1, 8),
    )), ([box, density]) => {
      // Arrange
      const { D, mesh } = box;

      // Act
      const one = nd.massProperties(D, mesh.vertices, mesh.cells, 1);
      const more = nd.massProperties(D, mesh.vertices, mesh.cells, density);

      // Assert
      assertClose(more.mass, one.mass * density, 'the mass', 1e-9);
      assertArrayClose(more.inertia, Float64Array.from(one.inertia, (v) => v * density),
        'the inertia', 1e-9);
    }));
  });

  it('should always give an inertia tensor that a body can have', () => {
    // An inertia tensor is symmetrical, and all of its eigenvalues are more
    // than zero. Without that, `matInverseSPD` falls back to Gauss-Jordan and
    // a body can gain energy.
    fc.assert(fc.property(anyBoxMesh(), (box) => {
      // Arrange
      const { D, mesh } = box;

      // Act
      const I = nd.massProperties(D, mesh.vertices, mesh.cells, 1).inertia;

      // Assert
      assert.ok(symmetryError(I, D.k) < 1e-12, 'the inertia must be symmetrical');
      assert.ok(isPositiveDefinite(I, D.k), 'the inertia must be positive definite');
    }));
  });

  it('should throw when the mesh holds no volume', () => {
    // Arrange -- four points on one line make no volume.
    const D = dims(3);
    const flat = Float64Array.from([0, 0, 0, 1, 0, 0, 2, 0, 0]);

    // Act, Assert
    assert.throws(() => nd.massProperties(D, flat, Int32Array.from([0, 1, 2]), 1), /volume/);
  });
});

describe('inertiaFromCovariance', () => {
  it('should always give the inertia of the mass properties', () => {
    // `massProperties` gives the covariance and the inertia together. The two
    // must agree, because A8 builds the second from the first.
    fc.assert(fc.property(anyBoxMesh(), (box) => {
      // Arrange
      const { D, mesh } = box;

      // Act
      const p = nd.massProperties(D, mesh.vertices, mesh.cells, 1);

      // Assert
      assertArrayClose(nd.inertiaFromCovariance(D, p.covariance), p.inertia,
        'the inertia from the covariance', 1e-9);
    }));
  });
});
