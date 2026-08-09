/**
 * test/spec/rotor.test.js -- src/nd/algebra/rotor.js
 *
 * The orientation of a body. This is the most important spec of the suite.
 * PART F of ND-PHYSICS.md says that the rotor correction is the most probable
 * source of a defect, and README section 10 says why: at `n` of 4 and more, a
 * length of 1 does not make a rotation.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { nd } from '../lib/load.js';
import {
  anyN, anyVector, anyBivector, anyRotor, anyEvenUnit, anyAngle, anyPlane, anyUnitVector,
} from '../lib/arbitraries.js';
import {
  assertClose, assertArrayClose, dot, norm, mul, mulMat, orthogonalError,
} from '../lib/numeric.js';

const { dims, rotor, linalg } = nd;

describe('rotorFromPlane', () => {
  it('should turn the first axis on to the second axis', () => {
    // Arrange -- this example fixes the sign of the convention. PART F says
    // to fix the sign one time and to write it down.
    const D = dims(4);
    const angle = 0.3;

    // Act
    const R = rotor.rotorFromPlane(D, 0, 1, angle);
    const x = rotor.rotorApplyVector(D, R, Float64Array.from([1, 0, 0, 0]));

    // Assert
    assertArrayClose(x, [Math.cos(angle), Math.sin(angle), 0, 0],
      'the turn goes from the axis 0 to the axis 1');
  });

  it('should throw when the two axes are the same', () => {
    assert.throws(() => rotor.rotorFromPlane(dims(3), 1, 1, 0.5), /different/);
  });

  it('should always add the angles of two turns in the same plane', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyPlane(n), anyAngle(), anyAngle(), anyVector(n),
    )), ([n, [i, j], a, b, x]) => {
      // Arrange
      const D = dims(n);

      // Act
      const two = rotor.rotorMul(D, rotor.rotorFromPlane(D, i, j, a), rotor.rotorFromPlane(D, i, j, b));
      const one = rotor.rotorFromPlane(D, i, j, a + b);

      // Assert -- R and -R are the same rotation, thus compare the results.
      assertArrayClose(rotor.rotorApplyVector(D, two, x), rotor.rotorApplyVector(D, one, x),
        'two turns in one plane make one turn', 1e-9);
    }));
  });
});

describe('rotorMul', () => {
  it('should always be associative', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyRotor(n), anyRotor(n),
    )), ([n, A, B, C]) => {
      // Arrange
      const D = dims(n);

      // Act
      const left = rotor.rotorMul(D, rotor.rotorMul(D, A, B), C);
      const right = rotor.rotorMul(D, A, rotor.rotorMul(D, B, C));

      // Assert
      assertArrayClose(left, right, '(A B) C is A (B C)', 1e-12);
    }));
  });

  it('should always turn a product around when it takes the reverse', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyRotor(n),
    )), ([n, A, B]) => {
      // Arrange
      const D = dims(n);

      // Act
      const left = rotor.rotorReverse(D, rotor.rotorMul(D, A, B));
      const right = rotor.rotorMul(D, rotor.rotorReverse(D, B), rotor.rotorReverse(D, A));

      // Assert
      assertArrayClose(left, right, '(A B)~ is B~ A~', 1e-12);
    }));
  });
});

describe('rotorApplyVector', () => {
  it('should never change the length of a vector', () => {
    // A rotation is an isometry. This is the definition, and it is the one
    // property that a body needs: a body does not change its size when it
    // turns.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyVector(n),
    )), ([n, R, x]) => {
      // Act
      const y = rotor.rotorApplyVector(dims(n), R, x);

      // Assert
      assertClose(norm(y), norm(x), 'the length after the turn', 1e-9);
    }));
  });

  it('should never change the angle between two vectors', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyVector(n), anyVector(n),
    )), ([n, R, a, b]) => {
      // Arrange
      const D = dims(n);

      // Act
      const ra = rotor.rotorApplyVector(D, R, a);
      const rb = rotor.rotorApplyVector(D, R, b);

      // Assert
      assertClose(dot(ra, rb), dot(a, b), 'the dot product after the turn', 1e-9);
    }));
  });

  it('should always agree with the rotation matrix', () => {
    // The step uses the matrix, and the algebra uses the sandwich product.
    // Nothing in the library compares the two.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyVector(n),
    )), ([n, R, x]) => {
      // Arrange
      const D = dims(n);

      // Act
      const sandwich = rotor.rotorApplyVector(D, R, x);
      const matrix = mul(rotor.rotorMatrix(D, R), x, n, n);

      // Assert
      assertArrayClose(sandwich, matrix, 'the sandwich and the matrix', 1e-9);
    }));
  });

  it('should always go back with the inverse', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyVector(n),
    )), ([n, R, x]) => {
      // Arrange
      const D = dims(n);

      // Act
      const there = rotor.rotorApplyVector(D, R, x);
      const back = rotor.rotorApplyVectorInverse(D, R, there);

      // Assert
      assertArrayClose(back, x, 'the turn and the turn back', 1e-9);
    }));
  });
});

describe('rotorMatrix', () => {
  it('should always give an orthogonal matrix with the determinant +1', () => {
    // A rotation, and not a reflection. A determinant of -1 would turn a body
    // inside out.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyRotor(n))),
      ([n, R]) => {
        // Act
        const M = rotor.rotorMatrix(dims(n), R);

        // Assert
        assert.ok(orthogonalError(M, n) < 1e-9, 'the matrix must be orthogonal');
        assertClose(linalg.matDet(M, n), 1, 'the determinant of a rotation', 1e-9);
      }));
  });

  it('should always make a product of rotors a product of matrices', () => {
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyRotor(n),
    )), ([n, A, B]) => {
      // Arrange
      const D = dims(n);

      // Act
      const left = rotor.rotorMatrix(D, rotor.rotorMul(D, A, B));
      const right = mulMat(rotor.rotorMatrix(D, A), rotor.rotorMatrix(D, B), n, n, n);

      // Assert
      assertArrayClose(left, right, 'the matrix of a product', 1e-9);
    }));
  });
});

describe('rotorBivectorMatrix', () => {
  it('should always agree with the sandwich product of a bivector', () => {
    // [R]2 B = R B R~. Item 1 of the test plan of PART E. `updateDerived`
    // changes the frame of the inertia tensor with this matrix.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyBivector(n),
    )), ([n, R, B]) => {
      // Arrange
      const D = dims(n);

      // Act
      const matrix = mul(rotor.rotorBivectorMatrix(D, R), B, D.k, D.k);
      const sandwich = rotor.rotorApplyBivector(D, R, B);

      // Assert
      assertArrayClose(matrix, sandwich, '[R]2 B is R B R~', 1e-9);
    }));
  });

  it('should always give an orthogonal matrix', () => {
    // The change of frame of the inertia tensor is [R]2 I [R]2^T. That keeps
    // the eigenvalues of the tensor only when [R]2 is orthogonal.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyRotor(n))),
      ([n, R]) => {
        const D = dims(n);
        assert.ok(orthogonalError(rotor.rotorBivectorMatrix(D, R), D.k) < 1e-9,
          '[R]2 must be orthogonal');
      }));
  });
});

describe('rotorApplyBivector', () => {
  it('should always make the wedge of the two turned vectors', () => {
    // R (a ^ b) R~ = (R a R~) ^ (R b R~). This is the reason that a bivector
    // is the correct type for the angular velocity, and it is a strong test
    // of the order of the bivector basis (PART F).
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyVector(n), anyVector(n),
    )), ([n, R, a, b]) => {
      // Arrange
      const D = dims(n);

      // Act
      const turned = rotor.rotorApplyBivector(D, R, nd.wedgeVec(D, a, b));
      const wedge = nd.wedgeVec(D,
        rotor.rotorApplyVector(D, R, a), rotor.rotorApplyVector(D, R, b));

      // Assert
      assertArrayClose(turned, wedge, 'the wedge follows the turn', 1e-8);
    }));
  });
});

describe('rotorExp', () => {
  it('should always give a rotor that is correct', () => {
    // The exponential of a bivector is always a member of the spin group.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyBivector(n))),
      ([n, B]) => {
        // Arrange
        const D = dims(n);

        // Act
        const R = rotor.rotorExp(D, B);

        // Assert
        assert.ok(rotor.rotorDefect(D, R) < 1e-9, `the defect is ${rotor.rotorDefect(D, R)}`);
        assert.ok(orthogonalError(rotor.rotorMatrix(D, R), n) < 1e-9, 'the matrix is orthogonal');
      }));
  });

  it('should always cut a double rotation into its two planes', () => {
    // In 4 dimensions a body can turn in two planes at the same time. The two
    // planes are orthogonal, thus the two bivectors commute, and
    // exp(B1 + B2) = exp(B1) exp(B2). A rotor of 4 components could not hold
    // this rotation.
    fc.assert(fc.property(fc.tuple(anyAngle(), anyAngle()), ([a, b]) => {
      // Arrange
      const D = dims(4);
      const B1 = new Float64Array(D.k);
      const B2 = new Float64Array(D.k);
      B1[D.biOfBlade[(1 << 0) | (1 << 1)]] = a;
      B2[D.biOfBlade[(1 << 2) | (1 << 3)]] = b;
      const sum = new Float64Array(D.k);
      for (let p = 0; p < D.k; p += 1) sum[p] = B1[p] + B2[p];

      // Act
      const together = rotor.rotorExp(D, sum);
      const apart = rotor.rotorMul(D, rotor.rotorExp(D, B1), rotor.rotorExp(D, B2));

      // Assert
      assertArrayClose(together, apart, 'exp(B1 + B2) is exp(B1) exp(B2)', 1e-9);
    }));
  });
});

describe('rotorFromBivectorAngle', () => {
  it('should agree with rotorFromPlane on one plane', () => {
    // `angle` is the turn of the BODY in the two functions. The exponent holds
    // the half angle. A full angle in the exponent made `Math.PI` give the
    // identity, and that made the half turn branch of `rotorBetweenVectors` do
    // nothing. See the note above `rotorCorrect` below.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyPlane(n), anyAngle(),
    )), ([n, [i, j], angle]) => {
      // Arrange
      const D = dims(n);
      const B = new Float64Array(D.k);
      B[D.biOfBlade[(1 << i) | (1 << j)]] = i < j ? 1 : -1;

      // Act
      const R = rotor.rotorFromBivectorAngle(D, B, angle);

      // Assert
      assertArrayClose(R, rotor.rotorFromPlane(D, i, j, angle),
        'the two ways to name a plane give the same rotor', 1e-12);
    }));
  });

  it('should give a half turn at pi', () => {
    // A half turn is the condition that `rotorCorrect` must hold. Take the
    // vector back and front: the two must be the same vector with the sign
    // changed, and not the vector itself.
    const D = dims(4);
    const B = new Float64Array(D.k);
    B[D.biOfBlade[(1 << 0) | (1 << 1)]] = 1;

    // Act
    const R = rotor.rotorFromBivectorAngle(D, B, Math.PI);
    const x = rotor.rotorApplyVector(D, R, Float64Array.from([0.6, 0.8, 0, 0]));

    // Assert
    assertArrayClose(x, [-0.6, -0.8, 0, 0], 'pi turns the vector against itself', 1e-12);
  });
});

describe('rotorBetweenVectors', () => {
  it('should always take the first vector on to the second', () => {
    // The direct form divides by `sqrt(2 (1 + dot))`, thus it loses all of its
    // accuracy as the two vectors come near opposite. `rotorCorrect` is the
    // one caller, and it meets that condition each time a body passes through
    // a half turn.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyUnitVector(n), anyUnitVector(n),
    )), ([n, a, b]) => {
      // Arrange
      const D = dims(n);

      // Act
      const R = rotor.rotorBetweenVectors(D, a, b);

      // Assert
      assertClose(rotor.rotorNorm(D, R), 1, 'the length of the rotor', 1e-12);
      assertArrayClose(rotor.rotorApplyVector(D, R, a), b, 'a comes on to b', 1e-9);
    }));
  });

  it('should always take a vector on to a vector that is almost opposite', () => {
    // The hard condition, one axis at a time. The angle comes to pi from an
    // error of 1e-3 down to 0, and it takes the exact half turn as well.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyUnitVector(n), anyPlane(n),
      fc.constantFrom(0, 1e-16, 1e-12, 1e-9, 1e-6, 1e-3),
    )), ([n, a, [i, j], gap]) => {
      // Arrange -- turn `a` through `pi - gap` in the plane of the axes i, j.
      const D = dims(n);
      const b = rotor.rotorApplyVector(D, rotor.rotorFromPlane(D, i, j, Math.PI - gap), a);

      // Act
      const R = rotor.rotorBetweenVectors(D, a, b);

      // Assert
      assertClose(rotor.rotorNorm(D, R), 1, 'the length of the rotor', 1e-12);
      assertArrayClose(rotor.rotorApplyVector(D, R, a), b, 'a comes on to b', 1e-9);
    }));
  });
});

describe('rotorDefect', () => {
  it('should give zero for the rotor that turns nothing', () => {
    assert.equal(rotor.rotorDefect(dims(4), rotor.rotorIdentity(dims(4))), 0);
  });

  it('should find the bad rotor that has the length 1', () => {
    // README section 10. (1 + e0 e1 e2 e3) / sqrt(2) has the length 1, but the
    // sandwich product with it is not a rotation. A 3D engine tests the
    // length of its quaternion; here that test would accept this element.
    const D = dims(4);
    const R = new Float64Array(D.r);
    R[D.scalarSlot] = Math.SQRT1_2;
    R[D.slotOfBlade[0b1111]] = Math.SQRT1_2;

    // Act, Assert
    assertClose(rotor.rotorNorm(D, R), 1, 'the length is exactly 1');
    assert.ok(rotor.rotorDefect(D, R) > 0.5, 'the defect finds it');
    assert.ok(orthogonalError(rotor.rotorMatrix(D, R), 4) > 0.5, 'the matrix is not a rotation');
  });

  it('should never say that a bad element of the length 1 is good', () => {
    // For each member of the even sub algebra of the length 1 at n of 4 and
    // 5: when the defect is more than 0.01, the matrix of the sandwich is not
    // orthogonal. Thus the length alone can never be the test.
    fc.assert(fc.property(fc.integer({ min: 4, max: 5 }).chain((n) => fc.tuple(
      fc.constant(n), anyEvenUnit(n),
    )), ([n, R]) => {
      // Arrange
      const D = dims(n);

      // Act
      const defect = rotor.rotorDefect(D, R);

      // Assert
      if (defect > 0.01) {
        assert.ok(orthogonalError(rotor.rotorMatrix(D, R), n) > 1e-9,
          `the defect is ${defect}, thus the matrix must not be orthogonal`);
      }
    }));
  });
});

/**
 * THE HALF TURN. These tests found three defects, and they hold the repair.
 *
 *   1. `rotorFromBivectorAngle` put the full angle in the exponent, thus the
 *      body turned through two times the angle. `Math.PI` gave the identity,
 *      and the "almost opposite" branch of `rotorBetweenVectors` did nothing.
 *   2. `rotorCorrect` did not move its work frame `cur` when the length of
 *      `q` was zero, but it had already put the turn into the product. The
 *      next column then took the same turn a second time.
 *   3. The divisor `sqrt(2 (1 + dot))` of `rotorBetweenVectors` loses all of
 *      its accuracy as `dot` comes near -1. About one of 2000 random rotors
 *      came back with an orthogonal error of more than 1e-9, which is the
 *      tolerance that `integratePositions` works with.
 *
 * A half turn is thus the condition to test, and not a rare one: a body that
 * turns passes through it two times in each revolution.
 */
describe('rotorCorrect', () => {
  it('should keep a half turn', () => {
    // Arrange
    const D = dims(3);
    const R = rotor.rotorFromPlane(D, 0, 1, Math.PI);
    const x = Float64Array.from([0.3, 0.5, 0.8]);

    // Act
    const C = rotor.rotorCorrect(D, R);

    // Assert
    assertArrayClose(rotor.rotorApplyVector(D, C, x), rotor.rotorApplyVector(D, R, x),
      'a half turn must stay a half turn');
  });

  it('should always give a rotor whose matrix is orthogonal', () => {
    // From any element of the even sub algebra, good or bad.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyEvenUnit(n))),
      ([n, R]) => {
        // Arrange
        const D = dims(n);

        // Act
        const C = rotor.rotorCorrect(D, R);

        // Assert
        assert.ok(orthogonalError(rotor.rotorMatrix(D, C), n) < 1e-9, 'the matrix is orthogonal');
        assert.ok(rotor.rotorDefect(D, C) < 1e-9, `the defect after the repair is ${rotor.rotorDefect(D, C)}`);
      }));
  });

  it('should almost do nothing to a rotor that is already correct', () => {
    // The repair runs after a step. It must not move a body that is good.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyRotor(n), anyUnitVector(n),
    )), ([n, R, x]) => {
      // Arrange
      const D = dims(n);

      // Act
      const C = rotor.rotorCorrect(D, R);

      // Assert -- compare the turn, because R and -R are the same rotation.
      assertArrayClose(rotor.rotorApplyVector(D, C, x), rotor.rotorApplyVector(D, R, x),
        'the repair does not turn the body', 1e-9);
    }));
  });

  it('should not move an axis that is in its place when a turn is near pi', () => {
    // The hard condition, and the reason for the orthogonalization of `q` in
    // step 3. `anyAngle` is uniform, thus it almost never gives an angle at
    // 1e-11 from pi, and a uniform test does not find this. Here the angles
    // come to pi on a ladder.
    //
    // With no orthogonalization the repair moved the body by 1e-4 at n = 4.
    fc.assert(fc.property(fc.integer({ min: 2, max: 5 }).chain((n) => fc.tuple(
      fc.constant(n),
      fc.array(fc.tuple(
        anyPlane(n),
        fc.constantFrom(0, 1e-16, 1e-13, 1e-11, 1e-9, 1e-6, 1e-3, 0.05, 0.2),
        fc.boolean(),
      ), { minLength: 1, maxLength: 4 }),
      anyUnitVector(n),
    )), ([n, list, x]) => {
      // Arrange
      const D = dims(n);
      let R = rotor.rotorIdentity(D);
      for (const [[i, j], gap, back] of list) {
        const angle = back ? -(Math.PI - gap) : Math.PI - gap;
        R = rotor.rotorMul(D, rotor.rotorFromPlane(D, i, j, angle), R);
      }

      // Act
      const C = rotor.rotorCorrect(D, R);

      // Assert
      assert.ok(orthogonalError(rotor.rotorMatrix(D, C), n) < 1e-9, 'the matrix is orthogonal');
      assertArrayClose(rotor.rotorApplyVector(D, C, x), rotor.rotorApplyVector(D, R, x),
        'the repair does not turn the body', 1e-9);
    }), { numRuns: 500 });
  });

  it('should never change the sign of a rotor', () => {
    // R and -R are the same rotation, but the plugin sends R to the display.
    // A change of the sign would make the body jump.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(fc.constant(n), anyRotor(n))),
      ([n, R]) => {
        // Arrange
        const D = dims(n);

        // Act
        const C = rotor.rotorCorrect(D, R);

        // Assert
        let d = 0;
        for (let i = 0; i < D.r; i += 1) d += C[i] * R[i];
        assert.ok(d >= 0, `the two rotors point the same way, and not ${d}`);
      }));
  });

  it('should keep a rotor correct after 200 products', () => {
    // Item 2 of the test plan of PART E: the rotor closure. Each product adds
    // a little error. The repair must take all of it away.
    fc.assert(fc.property(fc.integer({ min: 3, max: 5 }).chain((n) => fc.tuple(
      fc.constant(n),
      fc.array(fc.tuple(anyPlane(n), anyAngle()), { minLength: 200, maxLength: 200 }),
    )), ([n, list]) => {
      // Arrange
      const D = dims(n);

      // Act
      let R = rotor.rotorIdentity(D);
      for (const [[i, j], angle] of list) {
        R = rotor.rotorMul(D, rotor.rotorFromPlane(D, i, j, angle), R);
      }
      const C = rotor.rotorCorrect(D, R);

      // Assert -- the matrix must be orthogonal to 1e-9. The defect of the
      // rebuilt rotor is larger, because the rebuild of `rotorCorrect` is the
      // step that loses the accuracy, and not the 200 products.
      assert.ok(orthogonalError(rotor.rotorMatrix(D, C), n) < 1e-9, 'the matrix is orthogonal');
      assert.ok(rotor.rotorDefect(D, C) < 1e-7, `the defect after 200 products is ${rotor.rotorDefect(D, C)}`);
    }), { numRuns: 20 });
  });
});
