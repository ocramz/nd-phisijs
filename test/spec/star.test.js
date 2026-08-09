/**
 * test/spec/star.test.js -- src/nd/algebra/star.js
 *
 * The star matrix connects the vectors (position, force, impulse) to the
 * bivectors (torque, angular velocity, momentum). Almost all of the mechanics
 * uses it. In 3 dimensions it is the cross product matrix; at `n` of 4 and
 * more there is no cross product, and only these identities say that the code
 * is correct.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import { anyN, anyVector, anyBivector, anyRotor } from '../lib/arbitraries.js';
import { assertClose, assertArrayClose, dot, mul, mulT, norm } from '../lib/numeric.js';

const { dims, rotor } = nd;

describe('starMatrix', () => {
  it('should give the cross product matrix in 3 dimensions', () => {
    // Arrange -- the bivector basis of 3 dimensions is (x y), (x z), (y z).
    const D = dims(3);
    const r = Float64Array.from([1, 2, 3]);
    const a = Float64Array.from([4, 5, 6]);

    // Act
    const B = mul(nd.starMatrix(D, r), a, D.k, D.n);

    // Assert -- the components of r ^ a, in the order of the basis.
    assertArrayClose(B, [1 * 5 - 2 * 4, 1 * 6 - 3 * 4, 2 * 6 - 3 * 5], 'r ^ a of two vectors');
  });

  it('should always give the wedge product of the two vectors', () => {
    // [r]* a = r ^ a. Item 1 of the test plan of PART E.
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyVector(n),
    )), ([n, r, a]) => {
      // Arrange
      const D = dims(n);

      // Act
      const matrix = mul(nd.starMatrix(D, r), a, D.k, n);
      const wedge = nd.wedgeVec(D, r, a);

      // Assert
      assertArrayClose(matrix, wedge, '[r]* a is r ^ a', 1e-9);
    }));
  });

  it('should always give the contraction when it is transposed', () => {
    // [r]*^T W = r . W. The velocity of a point uses this form.
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyBivector(n),
    )), ([n, r, W]) => {
      // Arrange
      const D = dims(n);

      // Act
      const matrix = mulT(nd.starMatrix(D, r), W, D.k, n);
      const contraction = nd.contractVecBi(D, r, W);

      // Assert
      assertArrayClose(matrix, contraction, '[r]*^T W is r . W', 1e-9);
    }));
  });
});

describe('wedgeVec', () => {
  it('should always give zero for a vector with itself', () => {
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(fc.constant(n), anyVector(n))),
      ([n, a]) => {
        const B = nd.wedgeVec(dims(n), a, a);
        for (let p = 0; p < B.length; p += 1) assert.equal(B[p], 0, 'a ^ a is zero');
      }));
  });

  it('should always turn the sign around when the two vectors change place', () => {
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyVector(n),
    )), ([n, a, b]) => {
      // Arrange
      const D = dims(n);

      // Act
      const ab = nd.wedgeVec(D, a, b);
      const ba = nd.wedgeVec(D, b, a);

      // Assert
      for (let p = 0; p < D.k; p += 1) assert.ok(ab[p] + ba[p] === 0, 'a ^ b is -(b ^ a)');
    }));
  });

  it('should always be linear in the first vector', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyVector(n), anyVector(n),
    )), ([n, a, b, c]) => {
      // Arrange
      const D = dims(n);
      const sum = new Float64Array(n);
      for (let i = 0; i < n; i += 1) sum[i] = a[i] + b[i];

      // Act
      const together = nd.wedgeVec(D, sum, c);
      const apart = nd.wedgeVec(D, a, c);
      const other = nd.wedgeVec(D, b, c);

      // Assert
      for (let p = 0; p < D.k; p += 1) apart[p] += other[p];
      assertArrayClose(together, apart, '(a + b) ^ c', 1e-9);
    }));
  });
});

describe('contractVecBi', () => {
  it('should always give a vector that is normal to the offset', () => {
    // The velocity of a point of a body that turns is v + r . w. The part
    // r . w is always normal to r, thus the point does not move away from the
    // center of mass: the body is rigid.
    fc.assert(fc.property(anyN(2, 6).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n), anyBivector(n),
    )), ([n, r, W]) => {
      // Act
      const u = nd.contractVecBi(dims(n), r, W);

      // Assert
      assert.ok(Math.abs(dot(u, r)) <= 1e-9 * norm(r) ** 2 * (1 + norm(W)),
        `(r . W) . r must be zero, and it is ${dot(u, r)}`);
    }));
  });
});

describe('commutator', () => {
  it('should give the cross product of 3 dimensions, with the order turned around', () => {
    // The basis of 3 dimensions is (x y), (x z), (y z). `contractVecBi` shows
    // that the bivector (a, b, c) moves a point as the axial vector
    // (c, -b, a) does.
    //
    // With that map, [A, B] is the axial vector b x a, and not a x b. This is
    // the sign that makes the Euler equation of ND-PHYSICS.md, A7,
    // `I dw/dt - w x I w = tau`, the same as the classical
    // `I dw/dt + w x I w = tau`. PART F says to fix this sign one time and to
    // write it down. Here it is.
    const D = dims(3);
    const A = Float64Array.from([0, 0, 1]); // the axial vector x
    const B = Float64Array.from([1, 0, 0]); // the axial vector z

    // Act
    const C = nd.commutator(D, A, B);

    // Assert -- z x x is +y, and the bivector of +y is (0, -1, 0).
    assertArrayClose(C, [0, -1, 0], 'the commutator of two planes');
  });

  it('should always turn the sign around when the two bivectors change place', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyBivector(n), anyBivector(n),
    )), ([n, A, B]) => {
      // Arrange
      const D = dims(n);

      // Act
      const ab = nd.commutator(D, A, B);
      const ba = nd.commutator(D, B, A);

      // Assert
      assertArrayClose(ab, Float64Array.from(ba, (v) => -v), '[A, B] is -[B, A]', 1e-9);
    }));
  });

  it('should always obey the Jacobi identity', () => {
    // [A,[B,C]] + [B,[C,A]] + [C,[A,B]] = 0. The bivectors of n dimensions
    // are the Lie algebra so(n), and the gyroscopic term of A7 works only
    // when the commutator table obeys this law.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyBivector(n), anyBivector(n), anyBivector(n),
    )), ([n, A, B, C]) => {
      // Arrange
      const D = dims(n);

      // Act
      const one = nd.commutator(D, A, nd.commutator(D, B, C));
      const two = nd.commutator(D, B, nd.commutator(D, C, A));
      const three = nd.commutator(D, C, nd.commutator(D, A, B));

      // Assert -- the terms are products of three inputs, thus the tolerance
      // grows with the cube of the size.
      const scale = (norm(A) * norm(B) * norm(C)) || 1;
      for (let p = 0; p < D.k; p += 1) {
        assert.ok(Math.abs(one[p] + two[p] + three[p]) <= 1e-9 * scale,
          `the Jacobi identity at the component ${p}`);
      }
    }));
  });

  it('should always follow a turn of the two bivectors', () => {
    // [R B R~, R C R~] = R [B, C] R~. A rotor is an automorphism of the Lie
    // algebra. Without this the gyroscopic term would depend on the frame.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyBivector(n), anyBivector(n),
    )), ([n, R, B, C]) => {
      // Arrange
      const D = dims(n);

      // Act
      const turnedFirst = nd.commutator(D,
        rotor.rotorApplyBivector(D, R, B), rotor.rotorApplyBivector(D, R, C));
      const turnedLast = rotor.rotorApplyBivector(D, R, nd.commutator(D, B, C));

      // Assert
      const scale = (norm(B) * norm(C)) || 1;
      for (let p = 0; p < D.k; p += 1) {
        assert.ok(Math.abs(turnedFirst[p] - turnedLast[p]) <= 1e-9 * scale,
          `the commutator follows the turn, at the component ${p}`);
      }
    }));
  });
});

describe('commutatorMatrix', () => {
  it('should always agree with the commutator', () => {
    // `gyroscopicSpin` needs the matrix form for its Newton method.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyBivector(n), anyBivector(n),
    )), ([n, X, B]) => {
      // Arrange
      const D = dims(n);

      // Act
      const matrix = mul(nd.commutatorMatrix(D, X), B, D.k, D.k);
      const direct = nd.commutator(D, X, B);

      // Assert
      assertArrayClose(matrix, direct, 'the matrix of the commutator', 1e-9);
    }));
  });

  it('should always give zero on the diagonal of the energy', () => {
    // B . [X, B] = 0 for each B. The commutator with X is antisymmetric as a
    // matrix, thus the gyroscopic term does no work. This is the reason that
    // a body with no torque keeps its energy.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyBivector(n), anyBivector(n),
    )), ([n, X, B]) => {
      // Act
      const C = nd.commutator(dims(n), X, B);

      // Assert
      assertClose(dot(B, C) / ((norm(X) * norm(B) ** 2) || 1), 0, 'B . [X, B]', 1e-9);
    }));
  });
});
