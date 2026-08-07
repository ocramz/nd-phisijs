/**
 * test/spec/dims.test.js -- src/nd/core/dims.js
 *
 * The tables of `dims(n)` are the base of the whole library. A wrong sign
 * here gives a defect that no other test can explain. PART F of
 * ND-PHYSICS.md names two of these risks: the sign rule and the order of the
 * bivector basis.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import { anyN, anyVector } from '../lib/arbitraries.js';
import { assertArrayClose } from '../lib/numeric.js';

const { dims } = nd;

describe('dims', () => {
  it('should give the sizes of the three types', () => {
    // Arrange, Act
    const D = dims(4);

    // Assert
    assert.equal(D.k, 6, 'a bivector of 4 dimensions has 6 components');
    assert.equal(D.r, 8, 'a rotor of 4 dimensions has 8 components');
    assert.equal(D.N, 16, 'a multivector of 4 dimensions has 16 components');
  });

  it('should put the bivector components in lexicographic order', () => {
    // Arrange, Act
    const D = dims(4);

    // Assert -- README section 9 gives these two examples.
    assert.deepEqual(D.pairs[0], [0, 1], 'the index 0 is the plane (x y)');
    assert.deepEqual(D.pairs[2], [0, 3], 'the index 2 is the plane (x w)');
  });

  it('should throw when n is not an integer of 2 or more', () => {
    assert.throws(() => dims(1), /2 or more/);
    assert.throws(() => dims(2.5), /2 or more/);
  });

  it('should always give the same tables for the same n', () => {
    assert.equal(dims(4), dims(4), 'the tables are in a cache');
  });

  it('should always square an axis to +1', () => {
    // The space is euclidean. This one sign decides the sign of each product.
    fc.assert(fc.property(anyN(2, 6), (n) => {
      // Arrange
      const D = dims(n);

      // Act, Assert
      for (let i = 0; i < n; i += 1) {
        const idx = D.vecBlade[i] * D.N + D.vecBlade[i];
        assert.equal(D.gpBlade[idx], 0, 'the product of an axis with itself is a scalar');
        assert.equal(D.gpSign[idx], 1, `e${i} e${i} must be +1`);
      }
    }));
  });

  it('should always turn the sign around when two different axes change place', () => {
    fc.assert(fc.property(anyN(2, 6), (n) => {
      // Arrange
      const D = dims(n);

      // Act, Assert
      for (let i = 0; i < n; i += 1) {
        for (let j = 0; j < n; j += 1) {
          if (i === j) continue;
          const ab = D.vecBlade[i] * D.N + D.vecBlade[j];
          const ba = D.vecBlade[j] * D.N + D.vecBlade[i];
          assert.equal(D.gpBlade[ab], D.gpBlade[ba], 'the two products give the same blade');
          assert.equal(D.gpSign[ab], -D.gpSign[ba], `e${i} e${j} = -e${j} e${i}`);
        }
      }
    }));
  });

  it('should always give a commutator table that is antisymmetric', () => {
    // The commutator is (A B - B A) / 2. Thus [A, B] = -[B, A] for each pair
    // of components. The gyroscopic term of A7 depends on it.
    fc.assert(fc.property(anyN(2, 6), (n) => {
      // Arrange
      const { k, comm } = dims(n);

      // Act, Assert
      for (let p = 0; p < k; p += 1) {
        for (let q = 0; q < k; q += 1) {
          for (let s = 0; s < k; s += 1) {
            // The sum, and not `equal`: the table holds -0, and -0 is 0 here.
            assert.ok(comm[(p * k + q) * k + s] + comm[(q * k + p) * k + s] === 0,
              `comm[${p}][${q}][${s}] must be the opposite of comm[${q}][${p}][${s}]`);
          }
        }
      }
    }));
  });

  it('should always build the star matrix from the sum of r[i] eStar[i]', () => {
    // The comment of dims.js states this, and `starProducts` in massprops.js
    // depends on it: it builds the inertia tensor from `eStar` alone.
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(fc.constant(n), anyVector(n))),
      ([n, r]) => {
        // Arrange
        const D = dims(n);

        // Act
        const S = nd.starMatrix(D, r);
        const sum = new Float64Array(D.k * n);
        for (let i = 0; i < n; i += 1) {
          for (let t = 0; t < D.k * n; t += 1) sum[t] += r[i] * D.eStar[i][t];
        }

        // Assert
        assertArrayClose(S, sum, 'the star matrix is the sum of the axis matrices');
      }));
  });

  it('should always give a bivector index that finds its own pair again', () => {
    fc.assert(fc.property(anyN(2, 6), (n) => {
      // Arrange
      const D = dims(n);

      // Act, Assert
      for (let p = 0; p < D.k; p += 1) {
        const [i, j] = D.pairs[p];
        assert.ok(i < j, 'the first axis of a pair is the smaller one');
        assert.equal(D.biOfBlade[D.biBlade[p]], p, 'the blade of a pair finds the pair again');
        assert.equal(D.biBlade[p], (1 << i) | (1 << j), 'the blade holds the two axes');
      }
    }));
  });
});
