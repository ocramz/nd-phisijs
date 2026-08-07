/**
 * test/spec/linalg.test.js -- src/nd/core/linalg.js
 *
 * The small dense matrices. The step uses the fast forms (`matMulT`,
 * `matTVec`, `matInverseSPD`) and never the slow ones, thus nothing in the
 * library compares the two. These tests do that.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import { anyScalar } from '../lib/arbitraries.js';
import {
  assertClose, assertArrayClose, mulMat, mul, mulT, transpose, symmetryError,
} from '../lib/numeric.js';

const { linalg } = nd;

/** The side of a matrix. It covers `n` (2 to 6) and `k` (1 to 15). */
const anySide = () => fc.integer({ min: 1, max: 6 });

/** A matrix of `rows` by `cols`, in row major order. */
const anyMatrix = (rows, cols) => fc.tuple(
  ...Array.from({ length: rows * cols }, () => anyScalar()),
).map((a) => Float64Array.from(a));

/**
 * A matrix that is never singular. It adds the sum of the sizes of a row to
 * the diagonal, thus the matrix is diagonally dominant. This builds a valid
 * input, and it does not drop values with a filter.
 */
const anyInvertible = (m) => anyMatrix(m, m).map((A) => {
  const B = Float64Array.from(A);
  for (let i = 0; i < m; i += 1) {
    let sum = 0;
    for (let j = 0; j < m; j += 1) sum += Math.abs(B[i * m + j]);
    B[i * m + i] = sum + 1;
  }
  return B;
});

/** A matrix that is symmetrical and positive definite, as `L L^T + I`. */
const anySPD = (m) => anyMatrix(m, m).map((L) => {
  const A = new Float64Array(m * m);
  for (let i = 0; i < m; i += 1) {
    for (let j = 0; j < m; j += 1) {
      let s = i === j ? 1 : 0;
      for (let t = 0; t < m; t += 1) s += L[i * m + t] * L[j * m + t];
      A[i * m + j] = s;
    }
  }
  return A;
});

describe('matDet', () => {
  it('should give the determinant of a small matrix', () => {
    // Arrange
    const A = Float64Array.from([1, 2, 3, 4]);

    // Act, Assert
    assertClose(linalg.matDet(A, 2), -2, 'the determinant of [[1,2],[3,4]]');
  });

  it('should always give the product of the two determinants of a product', () => {
    // det(A B) = det(A) det(B). This one law finds almost every defect of a
    // determinant: a bad pivot, a lost sign, a wrong row.
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(
      fc.constant(m), anyMatrix(m, m), anyMatrix(m, m),
    )), ([m, A, B]) => {
      // Act
      const product = linalg.matDet(mulMat(A, B, m, m, m), m);
      const parts = linalg.matDet(A, m) * linalg.matDet(B, m);

      // Assert -- the tolerance is relative, thus a large determinant is safe.
      assertClose(product, parts, 'det(A B) is det(A) det(B)', 1e-6);
    }));
  });

  it('should never change the determinant when the matrix turns about its diagonal', () => {
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(fc.constant(m), anyMatrix(m, m))),
      ([m, A]) => {
        assertClose(linalg.matDet(transpose(A, m, m), m), linalg.matDet(A, m),
          'det(A^T) is det(A)', 1e-9);
      }));
  });

  it('should always turn the sign around when two rows change place', () => {
    fc.assert(fc.property(
      fc.integer({ min: 2, max: 6 }).chain((m) => fc.tuple(
        fc.constant(m), anyMatrix(m, m), fc.nat({ max: m - 1 }), fc.integer({ min: 1, max: m - 1 }),
      )),
      ([m, A, i, d]) => {
        // Arrange
        const j = (i + d) % m;
        const B = Float64Array.from(A);
        for (let c = 0; c < m; c += 1) {
          B[i * m + c] = A[j * m + c];
          B[j * m + c] = A[i * m + c];
        }

        // Act, Assert
        assertClose(linalg.matDet(B, m), -linalg.matDet(A, m),
          'an exchange of two rows turns the sign around', 1e-9);
      },
    ));
  });
});

describe('matInverse', () => {
  it('should always give a matrix that makes the identity', () => {
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(fc.constant(m), anyInvertible(m))),
      ([m, A]) => {
        // Act
        const I = mulMat(A, linalg.matInverse(A, m), m, m, m);

        // Assert
        for (let i = 0; i < m; i += 1) {
          for (let j = 0; j < m; j += 1) {
            assertClose(I[i * m + j], i === j ? 1 : 0, `A A^-1 at ${i},${j}`, 1e-9);
          }
        }
      }));
  });

  it('should throw when the matrix is singular', () => {
    assert.throws(() => linalg.matInverse(Float64Array.from([1, 2, 2, 4]), 2), /singular/);
  });
});

describe('matInverseSPD', () => {
  it('should always agree with matInverse on a positive definite matrix', () => {
    // Cholesky is the fast path of the inertia tensor, and the step uses it
    // in each frame. Gauss-Jordan is the simple form. The two must agree.
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(fc.constant(m), anySPD(m))),
      ([m, A]) => {
        assertArrayClose(linalg.matInverseSPD(A, m), linalg.matInverse(A, m),
          'the two inverses of a positive definite matrix', 1e-6);
      }));
  });

  it('should always give a symmetrical matrix', () => {
    // The inverse of a symmetrical matrix is symmetrical. An inertia tensor
    // that is not symmetrical gives a body that gains energy.
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(fc.constant(m), anySPD(m))),
      ([m, A]) => {
        assert.ok(symmetryError(linalg.matInverseSPD(A, m), m) < 1e-12,
          'the inverse must be symmetrical');
      }));
  });
});

describe('matMulT', () => {
  it('should always agree with the product with the transpose', () => {
    // `updateDerived` builds [R]2 I [R]2^T with this function, and it never
    // builds the transpose in memory.
    fc.assert(fc.property(
      fc.tuple(anySide(), anySide(), anySide()).chain(([ra, ca, rb]) => fc.tuple(
        fc.constant([ra, ca, rb]), anyMatrix(ra, ca), anyMatrix(rb, ca),
      )),
      ([[ra, ca, rb], A, B]) => {
        assertArrayClose(linalg.matMulT(A, B, ra, ca, rb),
          mulMat(A, transpose(B, rb, ca), ra, ca, rb), 'A B^T', 1e-9);
      },
    ));
  });
});

describe('matTVec', () => {
  it('should always agree with the product of the transpose and the vector', () => {
    // `worldToLocalDir` and the star matrix use this form.
    fc.assert(fc.property(
      fc.tuple(anySide(), anySide()).chain(([rows, cols]) => fc.tuple(
        fc.constant([rows, cols]), anyMatrix(rows, cols), anyMatrix(rows, 1),
      )),
      ([[rows, cols], A, x]) => {
        assertArrayClose(linalg.matTVec(A, x, rows, cols), mulT(A, x, rows, cols),
          'A^T x', 1e-9);
      },
    ));
  });
});

describe('matVec', () => {
  it('should always be associative with the product of two matrices', () => {
    fc.assert(fc.property(
      fc.tuple(anySide(), anySide(), anySide()).chain(([ra, ca, cb]) => fc.tuple(
        fc.constant([ra, ca, cb]), anyMatrix(ra, ca), anyMatrix(ca, cb), anyMatrix(cb, 1),
      )),
      ([[ra, ca, cb], A, B, x]) => {
        // Act
        const left = linalg.matVec(mulMat(A, B, ra, ca, cb), x, ra, cb);
        const right = mul(A, linalg.matVec(B, x, ca, cb), ra, ca);

        // Assert
        assertArrayClose(left, right, '(A B) x is A (B x)', 1e-6);
      },
    ));
  });
});

describe('jacobiEigen', () => {
  it('should always give values and vectors that build the matrix again', () => {
    // A V = V diag(values). The specs of the mass properties use this to see
    // that an inertia tensor is positive definite.
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(fc.constant(m), anySPD(m))),
      ([m, A]) => {
        // Act
        const { values, vectors } = linalg.jacobiEigen(A, m);

        // Assert
        for (let c = 0; c < m; c += 1) {
          const v = new Float64Array(m);
          for (let i = 0; i < m; i += 1) v[i] = vectors[i * m + c];
          const Av = mul(A, v, m, m);
          for (let i = 0; i < m; i += 1) {
            assertClose(Av[i], values[c] * v[i], `the eigenvector ${c}`, 1e-6);
          }
        }
      }));
  });

  it('should always give a sum of the values that is the trace', () => {
    fc.assert(fc.property(anySide().chain((m) => fc.tuple(fc.constant(m), anySPD(m))),
      ([m, A]) => {
        // Act
        const { values } = linalg.jacobiEigen(A, m);

        // Assert
        let sum = 0;
        let trace = 0;
        for (let i = 0; i < m; i += 1) {
          sum += values[i];
          trace += A[i * m + i];
        }
        assertClose(sum, trace, 'the sum of the eigenvalues is the trace', 1e-9);
      }));
  });
});
