/**
 * test/lib/arbitraries.js -- the fast-check arbitraries of the types of the
 * library.
 *
 * A generator here builds a VALID value of the type. A rotor is not a list of
 * random numbers: it is a product of simple rotors, thus it is always a true
 * rotation. `anyEven` gives the random numbers, and only the tests of
 * `rotorDefect` and `rotorCorrect` use it.
 *
 * THE LIMITS OF THE RANGES. `anyScalar` stays between -1e3 and 1e3, and a
 * half extent stays between 0.1 and 3. This is a requirement of the code and
 * not a convenience: `defaultParams` holds ABSOLUTE tolerances
 * (`penetrationSlop` is 5e-3 and `contactMargin` is 0.02), thus the library
 * is only defined on a physical range of length. A body of the size 1e300 has
 * no meaning for the solver.
 *
 * WHEN A TEST FAILS, fast-check writes the seed and the path. Run that one
 * case again with:
 *   fc.assert(property, { seed: 1234, path: '0:0' })
 */
import fc from 'fast-check';
import { nd } from './load.js';

const { dims, rotor } = nd;

/**
 * The count of dimensions. The tables of `dims(n)` have the size `2^n`, thus
 * a large `n` is slow. ND-PHYSICS.md asks for 2 to 6; the specs use 2 to 5,
 * and `dims.test.js` adds 6 by itself.
 */
export function anyN(min = 2, max = 5) {
  return fc.integer({ min, max });
}

/** One number of a physical size. */
export function anyScalar() {
  return fc.double({ min: -1e3, max: 1e3, noNaN: true });
}

/** One number that is not near zero, for a divisor or a mass. */
export function anyPositive(min = 0.1, max = 10) {
  return fc.double({ min, max, noNaN: true });
}

/** An angle in radians. */
export function anyAngle() {
  return fc.double({ min: -Math.PI, max: Math.PI, noNaN: true });
}

/** An array of `count` numbers from `item`, as a `Float64Array`. */
function fixedArray(count, item) {
  return fc.tuple(...Array.from({ length: count }, () => item)).map((a) => Float64Array.from(a));
}

/** A vector of `n` components. */
export function anyVector(n, item = anyScalar()) {
  return fixedArray(n, item);
}

/**
 * A vector of the length 1. It divides a vector by its length. When the
 * vector is almost zero it gives an axis, thus the generator drops no value.
 */
export function anyUnitVector(n) {
  return fc.tuple(anyVector(n), fc.nat({ max: n - 1 })).map(([v, axis]) => {
    let ln = 0;
    for (let i = 0; i < n; i += 1) ln += v[i] * v[i];
    ln = Math.sqrt(ln);
    if (ln < 1e-8) {
      const e = new Float64Array(n);
      e[axis] = 1;
      return e;
    }
    const u = new Float64Array(n);
    for (let i = 0; i < n; i += 1) u[i] = v[i] / ln;
    return u;
  });
}

/** A bivector of `k = n (n-1) / 2` components. */
export function anyBivector(n, item = anyScalar()) {
  return fixedArray(dims(n).k, item);
}

/** A bivector of a small size, for a rotation of a body in a step. */
export function anyAngularVelocity(n) {
  return anyBivector(n, fc.double({ min: -4, max: 4, noNaN: true }));
}

/** Two different axes, in order. */
export function anyPlane(n) {
  return fc.tuple(fc.nat({ max: n - 1 }), fc.integer({ min: 1, max: n - 1 }))
    .map(([i, d]) => {
      const j = (i + d) % n;
      return i < j ? [i, j] : [j, i];
    });
}

/**
 * A rotor. It is the product of 1 to 4 simple rotors. Each rotation of `n`
 * dimensions is a product of `floor(n / 2)` simple rotations, thus 4 covers
 * each `n` to 8, and the product is always a true rotation.
 */
export function anyRotor(n) {
  const D = dims(n);
  const simple = fc.tuple(anyPlane(n), anyAngle());
  return fc.array(simple, { minLength: 1, maxLength: 4 }).map((list) => {
    let R = rotor.rotorIdentity(D);
    for (const [[i, j], angle] of list) {
      R = rotor.rotorMul(D, rotor.rotorFromPlane(D, i, j, angle), R);
    }
    return R;
  });
}

/**
 * A member of the even sub algebra of the length 1. This is NOT a rotor: at
 * `n` of 4 and more, the length 1 does not make a rotation. See README,
 * section 10.
 */
export function anyEvenUnit(n) {
  const D = dims(n);
  return fc.tuple(...Array.from({ length: D.r }, () => anyScalar())).map((a) => {
    const R = Float64Array.from(a);
    let ln = 0;
    for (let i = 0; i < D.r; i += 1) ln += R[i] * R[i];
    ln = Math.sqrt(ln);
    if (ln < 1e-8) return rotor.rotorIdentity(D);
    for (let i = 0; i < D.r; i += 1) R[i] /= ln;
    return R;
  });
}

/** The half extents of a box, all more than zero. */
export function anyHalfExtents(n, min = 0.1, max = 3) {
  return anyVector(n, fc.double({ min, max, noNaN: true }));
}

/** The radius of a ball or of the tube of a torus. */
export function anyRadius(min = 0.2, max = 2) {
  return fc.double({ min, max, noNaN: true });
}

/** A position in a box of the side `2 * reach` about the origin. */
export function anyPosition(n, reach = 5) {
  return anyVector(n, fc.double({ min: -reach, max: reach, noNaN: true }));
}

/**
 * A shape with a mesh: a box, or the same box as a convex mesh. The two must
 * give the same mass properties, thus a test can compare them.
 */
export function anyMeshShape(D) {
  return anyHalfExtents(D.n).chain((h) => fc.constantFrom(
    nd.HyperBox(D, h),
    (() => {
      const m = nd.hyperBoxMesh(D, h);
      return nd.ConvexMesh(D, m.vertices, m.cells);
    })(),
  ));
}

/**
 * A shape that can move: a box, a ball or a torus. A half space is always
 * static, thus it is not here. `withTorus` is false at `n` of 5 and more,
 * because a torus has no mesh there.
 */
export function anyDynamicShape(D, withTorus = D.n === 3 || D.n === 4) {
  const parts = [
    anyHalfExtents(D.n).map((h) => nd.HyperBox(D, h)),
    anyRadius().map((r) => nd.HyperSphere(D, r)),
  ];
  if (withTorus) {
    parts.push(fc.tuple(anyRadius(0.8, 2), anyRadius(0.2, 0.5), anyPlane(D.n))
      .map(([R, r, plane]) => nd.Torus(D, R, r, plane)));
  }
  return fc.oneof(...parts);
}

/** A body that can move, with a shape, a place, a turn and a speed. */
export function anyBody(D, opts = {}) {
  return fc.record({
    shape: opts.shape || anyDynamicShape(D),
    mass: anyPositive(0.2, 5),
    position: anyPosition(D.n, opts.reach === undefined ? 3 : opts.reach),
    rotor: anyRotor(D.n),
    velocity: anyVector(D.n, fc.double({ min: -5, max: 5, noNaN: true })),
    angularVelocity: anyAngularVelocity(D.n),
  }).map((o) => new nd.Body(D, o));
}

/** The `D` of each `n` of `min` to `max`, one at a time. */
export function anyDims(min = 2, max = 5) {
  return anyN(min, max).map((n) => dims(n));
}
