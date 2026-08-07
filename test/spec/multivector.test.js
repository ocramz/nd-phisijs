/**
 * test/spec/multivector.test.js -- src/nd/algebra/multivector.js
 *
 * The full multivector. The step does not use one, but the tables of the
 * products come from `dims()`, and the rotor sandwich, the separating axis
 * and the mass properties all use those tables. Thus a defect here is a
 * defect everywhere.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import { anyN, anyVector } from '../lib/arbitraries.js';
import { assertClose, dot, norm } from '../lib/numeric.js';

const { dims, mv } = nd;

describe('mvGp', () => {
  it('should give the geometric product of two axes', () => {
    // Arrange
    const D = dims(3);
    const x = mv.mvFromVector(D, [1, 0, 0]);
    const y = mv.mvFromVector(D, [0, 1, 0]);

    // Act
    const P = mv.mvGp(D, x, y);

    // Assert -- e0 e1 is the plane (x y), and it has no scalar part.
    assert.equal(P[0], 0, 'two different axes give no scalar');
    assert.equal(P[(1 << 0) | (1 << 1)], 1, 'the product is the plane of the two axes');
  });

  it('should always be associative', () => {
    // (A B) C = A (B C). This is the strongest test of the sign table.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyMultivector(n), anyMultivector(n), anyMultivector(n),
    )), ([n, A, B, C]) => {
      // Arrange
      const D = dims(n);

      // Act
      const left = mv.mvGp(D, mv.mvGp(D, A, B), C);
      const right = mv.mvGp(D, A, mv.mvGp(D, B, C));

      // Assert
      assertProportional(left, right, size(A) * size(B) * size(C), '(A B) C is A (B C)');
    }));
  });

  it('should always be distributive over the sum', () => {
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyMultivector(n), anyMultivector(n), anyMultivector(n),
    )), ([n, A, B, C]) => {
      // Arrange
      const D = dims(n);
      const sum = new Float64Array(D.N);
      for (let m = 0; m < D.N; m += 1) sum[m] = B[m] + C[m];

      // Act
      const left = mv.mvGp(D, A, sum);
      const right = mv.mvGp(D, A, B);
      const other = mv.mvGp(D, A, C);

      // Assert
      for (let m = 0; m < D.N; m += 1) right[m] += other[m];
      assertProportional(left, right, size(A) * (size(B) + size(C)), 'A (B + C)');
    }));
  });

  it('should always give the square of the length for a vector', () => {
    // a a = a . a. This holds only in a euclidean space, and it is the sign
    // rule of dims() in one number.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyVector(n))),
      ([n, a]) => {
        // Arrange
        const D = dims(n);
        const A = mv.mvFromVector(D, a);

        // Act
        const P = mv.mvGp(D, A, A);

        // Assert
        assertClose(P[0], dot(a, a), 'the scalar part of a a', 1e-9);
        for (let m = 1; m < D.N; m += 1) {
          assert.ok(Math.abs(P[m]) <= 1e-9 * dot(a, a) + 1e-12,
            'a a has no part that is not a scalar');
        }
      }));
  });

  it('should always give the dot product and the wedge for two vectors', () => {
    // a b = a . b + a ^ b. See ND-PHYSICS.md, A2.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyVector(n),
    )), ([n, a, b]) => {
      // Arrange
      const D = dims(n);

      // Act
      const P = mv.mvGp(D, mv.mvFromVector(D, a), mv.mvFromVector(D, b));

      // Assert
      assertClose(P[0], dot(a, b), 'the scalar part is the dot product', 1e-9);
      const wedge = nd.wedgeVec(D, a, b);
      const part = mv.mvToBivector(D, P);
      assertProportional(part, wedge, norm(a) * norm(b), 'the grade 2 part is the wedge');
    }));
  });
});

describe('mvReverse', () => {
  it('should always turn a product around', () => {
    // (A B)~ = B~ A~. The sandwich R x R~ is only an isometry with this law.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyMultivector(n), anyMultivector(n),
    )), ([n, A, B]) => {
      // Arrange
      const D = dims(n);

      // Act
      const left = mv.mvReverse(D, mv.mvGp(D, A, B));
      const right = mv.mvGp(D, mv.mvReverse(D, B), mv.mvReverse(D, A));

      // Assert
      assertProportional(left, right, size(A) * size(B), '(A B)~ is B~ A~');
    }));
  });
});

describe('mvWedge', () => {
  it('should always give zero for a vector with itself', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyVector(n))),
      ([n, a]) => {
        // Arrange
        const D = dims(n);
        const A = mv.mvFromVector(D, a);

        // Act
        const P = mv.mvWedge(D, A, A);

        // Assert
        for (let m = 0; m < D.N; m += 1) assert.equal(P[m], 0, 'a ^ a is zero');
      }));
  });

  it('should always turn the sign around for two vectors', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyVector(n),
    )), ([n, a, b]) => {
      // Arrange
      const D = dims(n);
      const A = mv.mvFromVector(D, a);
      const B = mv.mvFromVector(D, b);

      // Act
      const ab = mv.mvWedge(D, A, B);
      const ba = mv.mvWedge(D, B, A);

      // Assert
      for (let m = 0; m < D.N; m += 1) assert.ok(ab[m] + ba[m] === 0, 'a ^ b is -(b ^ a)');
    }));
  });

  it('should always add the two grades', () => {
    // The wedge of a blade of the grade p and a blade of the grade q has the
    // grade p + q only.
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyMultivector(n), anyMultivector(n),
      fc.nat({ max: n }), fc.nat({ max: n }),
    )), ([n, A, B, p, q]) => {
      // Arrange
      const D = dims(n);
      const Ap = mv.mvGrade(D, A, p);
      const Bq = mv.mvGrade(D, B, q);

      // Act
      const P = mv.mvWedge(D, Ap, Bq);

      // Assert
      for (let m = 0; m < D.N; m += 1) {
        if (D.grade[m] !== p + q) assert.equal(P[m], 0, `the grade ${D.grade[m]} must be zero`);
      }
    }));
  });
});

describe('mvContract', () => {
  it('should always take the first grade away from the second', () => {
    fc.assert(fc.property(anyN(2, 4).chain((n) => fc.tuple(
      fc.constant(n), anyMultivector(n), anyMultivector(n),
      fc.nat({ max: n }), fc.nat({ max: n }),
    )), ([n, A, B, p, q]) => {
      // Arrange
      const D = dims(n);

      // Act
      const P = mv.mvContract(D, mv.mvGrade(D, A, p), mv.mvGrade(D, B, q));

      // Assert
      for (let m = 0; m < D.N; m += 1) {
        if (D.grade[m] !== q - p) assert.equal(P[m], 0, `the grade ${D.grade[m]} must be zero`);
      }
    }));
  });
});

describe('mvDual', () => {
  it('should always give a vector that is normal to the wedge of n - 1 vectors', () => {
    // `boxBoxAxis` builds a separating axis in this way: it wedges n - 1 axes
    // of the two boxes together, then it takes the dual. PART F of
    // ND-PHYSICS.md says that this is the one correct use of a dual.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), fc.array(anyVector(n), { minLength: n - 1, maxLength: n - 1 }),
    )), ([n, list]) => {
      // Arrange
      const D = dims(n);

      // Act
      let blade = mv.mvFromVector(D, list[0]);
      for (let i = 1; i < list.length; i += 1) {
        blade = mv.mvWedge(D, blade, mv.mvFromVector(D, list[i]));
      }
      const axis = mv.mvToVector(D, mv.mvDual(D, blade));

      // Assert
      let scale = norm(axis);
      for (const v of list) scale = Math.max(scale, norm(v));
      for (const v of list) {
        assert.ok(Math.abs(dot(axis, v)) <= 1e-9 * scale ** n,
          'the dual is normal to each vector of the blade');
      }
    }));
  });
});

// Helpers

/** A multivector of `2^n` components. */
function anyMultivector(n) {
  const D = dims(n);
  return fc.tuple(...Array.from({ length: D.N },
    () => fc.double({ min: -10, max: 10, noNaN: true }))).map((a) => Float64Array.from(a));
}

/** The largest component of a multivector. It scales the tolerance. */
function size(A) {
  let s = 1;
  for (let i = 0; i < A.length; i += 1) s = Math.max(s, Math.abs(A[i]));
  return s * A.length;
}

/**
 * Stops the test when the two arrays differ by more than `1e-9 * scale`.
 *
 * A product of multivectors adds `2^n` terms, and the terms cancel. Thus the
 * error grows with the size of the inputs, and not with the size of the
 * result. The tolerance must do the same.
 */
function assertProportional(actual, expected, scale, message) {
  for (let i = 0; i < actual.length; i += 1) {
    assert.ok(Math.abs(actual[i] - expected[i]) <= 1e-9 * scale,
      `${message}: at the index ${i}, ${actual[i]} is not ${expected[i]}`);
  }
}
