/**
 * test/lib/numeric.js -- the comparisons and the small operations of the
 * tests.
 *
 * THE TOLERANCE IS RELATIVE. The engine works with numbers of a large range:
 * an inertia tensor of a body of the size 3 in 5 dimensions holds values of
 * more than 100, and a velocity can be 1e-4. An absolute tolerance would fail
 * on the first and it would accept a defect in the second. Thus each
 * comparison here is
 *
 *   |a - b| <= rel * max(1, |a|, |b|)
 *
 * The code of this file does not come from the library. A test must not use
 * the code that it tests to say what the answer is.
 */
import assert from 'node:assert/strict';

/** The tolerance of a comparison of two results of the same algebra. */
export const REL = 1e-9;

/** True when `a` and `b` are equal, within the relative tolerance `rel`. */
export function close(a, b, rel = REL) {
  const scale = Math.max(1, Math.abs(a), Math.abs(b));
  return Math.abs(a - b) <= rel * scale;
}

/** Stops the test when `a` and `b` are not equal. */
export function assertClose(a, b, message = '', rel = REL) {
  assert.ok(close(a, b, rel), `${message}: ${a} is not ${b} (tolerance ${rel})`);
}

/** Stops the test when the two arrays are not equal, element by element. */
export function assertArrayClose(a, b, message = '', rel = REL) {
  assert.equal(a.length, b.length, `${message}: the lengths are not the same`);
  for (let i = 0; i < a.length; i += 1) {
    assert.ok(close(a[i], b[i], rel),
      `${message}: at the index ${i}, ${a[i]} is not ${b[i]} (tolerance ${rel})`);
  }
}

/** The largest size of an element of an array. It is 1 at the smallest. */
export function maxAbs(A) {
  let s = 1;
  for (let i = 0; i < A.length; i += 1) s = Math.max(s, Math.abs(A[i]));
  return s;
}

/**
 * Stops the test when the two arrays differ by more than `rel * scale`.
 *
 * Use this for a matrix in which some elements are zero and others are large,
 * such as an inertia tensor. An error of 1e-14 in an element that must be
 * zero is not a defect, but a relative comparison of that one element would
 * say that it is.
 */
export function assertScaled(a, b, scale, message = '', rel = REL) {
  assert.equal(a.length, b.length, `${message}: the lengths are not the same`);
  for (let i = 0; i < a.length; i += 1) {
    assert.ok(Math.abs(a[i] - b[i]) <= rel * scale,
      `${message}: at the index ${i}, ${a[i]} is not ${b[i]} (tolerance ${rel * scale})`);
  }
}

/** The dot product of two arrays of the same length. */
export function dot(a, b) {
  let s = 0;
  for (let i = 0; i < a.length; i += 1) s += a[i] * b[i];
  return s;
}

/** The length of a vector. */
export function norm(a) {
  return Math.sqrt(dot(a, a));
}

/** `a + f b`, as a new array. */
export function addScaled(a, b, f) {
  const out = new Float64Array(a.length);
  for (let i = 0; i < a.length; i += 1) out[i] = a[i] + f * b[i];
  return out;
}

/** `f a`, as a new array. */
export function scaled(a, f) {
  const out = new Float64Array(a.length);
  for (let i = 0; i < a.length; i += 1) out[i] = a[i] * f;
  return out;
}

/** The vector of the length 1 in the direction of `a`. */
export function normalized(a) {
  const ln = norm(a);
  return scaled(a, ln > 0 ? 1 / ln : 0);
}

/** The product `A x`. `A` is `rows` by `cols`, in row major order. */
export function mul(A, x, rows, cols) {
  const y = new Float64Array(rows);
  for (let i = 0; i < rows; i += 1) {
    let s = 0;
    for (let j = 0; j < cols; j += 1) s += A[i * cols + j] * x[j];
    y[i] = s;
  }
  return y;
}

/** The product `A^T x`. `A` is `rows` by `cols`. */
export function mulT(A, x, rows, cols) {
  const y = new Float64Array(cols);
  for (let j = 0; j < cols; j += 1) {
    let s = 0;
    for (let i = 0; i < rows; i += 1) s += A[i * cols + j] * x[i];
    y[j] = s;
  }
  return y;
}

/** The product `A B`. `A` is `ra` by `ca`, and `B` is `ca` by `cb`. */
export function mulMat(A, B, ra, ca, cb) {
  const C = new Float64Array(ra * cb);
  for (let i = 0; i < ra; i += 1) {
    for (let j = 0; j < cb; j += 1) {
      let s = 0;
      for (let t = 0; t < ca; t += 1) s += A[i * ca + t] * B[t * cb + j];
      C[i * cb + j] = s;
    }
  }
  return C;
}

/** The transpose of `A`, of `rows` by `cols`. */
export function transpose(A, rows, cols) {
  const B = new Float64Array(rows * cols);
  for (let i = 0; i < rows; i += 1) for (let j = 0; j < cols; j += 1) B[j * rows + i] = A[i * cols + j];
  return B;
}

/** The largest difference between `M^T M` and the identity of `m` by `m`. */
export function orthogonalError(M, m) {
  let worst = 0;
  for (let i = 0; i < m; i += 1) {
    for (let j = 0; j < m; j += 1) {
      let s = 0;
      for (let t = 0; t < m; t += 1) s += M[t * m + i] * M[t * m + j];
      worst = Math.max(worst, Math.abs(s - (i === j ? 1 : 0)));
    }
  }
  return worst;
}

/** The largest difference between `M` and its transpose. */
export function symmetryError(M, m) {
  let worst = 0;
  for (let i = 0; i < m; i += 1) {
    for (let j = 0; j < m; j += 1) worst = Math.max(worst, Math.abs(M[i * m + j] - M[j * m + i]));
  }
  return worst;
}

/**
 * True when the symmetrical matrix `A` of `m` by `m` is positive definite.
 * The method is Cholesky: the cut `A = L L^T` is only possible when each
 * diagonal element stays more than zero.
 */
export function isPositiveDefinite(A, m) {
  const L = new Float64Array(m * m);
  for (let i = 0; i < m; i += 1) {
    for (let j = 0; j <= i; j += 1) {
      let s = A[i * m + j];
      for (let t = 0; t < j; t += 1) s -= L[i * m + t] * L[j * m + t];
      if (i === j) {
        if (!(s > 0)) return false;
        L[i * m + j] = Math.sqrt(s);
      } else {
        L[i * m + j] = s / L[j * m + j];
      }
    }
  }
  return true;
}

/**
 * The area of a set of triangles. `positions` holds 9 numbers for each
 * triangle. The area does not change with the direction of a triangle, thus
 * this is a safe measure of a surface that comes from a cut.
 */
export function soupArea(positions) {
  let area = 0;
  for (let t = 0; t + 8 < positions.length; t += 9) {
    const ux = positions[t + 3] - positions[t];
    const uy = positions[t + 4] - positions[t + 1];
    const uz = positions[t + 5] - positions[t + 2];
    const vx = positions[t + 6] - positions[t];
    const vy = positions[t + 7] - positions[t + 1];
    const vz = positions[t + 8] - positions[t + 2];
    area += 0.5 * Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
  }
  return area;
}

/** The smallest and the largest value of one axis of a set of triangles. */
export function soupExtent(positions, axis) {
  let min = Infinity;
  let max = -Infinity;
  for (let i = axis; i < positions.length; i += 3) {
    min = Math.min(min, positions[i]);
    max = Math.max(max, positions[i]);
  }
  return { min, max };
}
