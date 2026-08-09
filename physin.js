/**
 * PhysiN 0.1.0 -- N-dimensional rigid body dynamics for three.js.
 * https://github.com/  (MIT)
 *
 * Single file containing algebra, shapes, mass properties, the
 * integrator, collision detection, the solver, the 4D slice and the
 * three.js plugin.
 *
 * HOW TO USE IT, as with physi.js:
 *   1. Load three.js, then load this file with a <script> tag.
 *   2. Point PhysiN.scripts.worker at the worker file:
 *        PhysiN.scripts.worker = 'physiN_worker.js';
 *      Leave it as null to run the physics in the main thread.
 *   3. Use PhysiN.Scene in the place of THREE.Scene.
 *   4. Use PhysiN.BoxMesh, PhysiN.SphereMesh, PhysiN.ConvexMesh,
 *      PhysiN.TorusMesh or PhysiN.PlaneMesh in the place of THREE.Mesh.
 *      For four dimensions and more, use PhysiN.HyperBoxMesh,
 *      PhysiN.HyperSphereMesh, PhysiN.HyperTorusMesh, PhysiN.HyperMesh or
 *      PhysiN.HyperPlaneMesh.
 *   5. Call scene.simulate(timeStep, maxSubSteps) at each frame.
 *
 * CAUTION for four dimensions and more: the angular velocity, the torque and
 * the angular momentum are BIVECTORS with n(n-1)/2 components. They are not
 * vectors. The orientation is a ROTOR with 2^(n-1) components. It is not a
 * quaternion. See the README.
 *
 * HOW TO READ THIS FILE. The `// src/...` lines show the modules. Three type
 * blocks give the names that all of the other code uses:
 *   `Dims`     above `dims()`, in src/nd/core/dims.js
 *   `Shape`    above `HyperBox()`, in src/nd/body/shapes.js
 *   `Contact`  above `makeContact()`, in src/nd/detect/collide.js
 * ND-PHYSICS.md holds the mathematics, and the blocks below point to it by
 * the number of its parts, for example "A3" or "B2".
 *
 * CAUTION for a change: physiN_worker.js holds the same engine modules as
 * this file, word for word. It leaves out only the code that a worker does
 * not use. Make the same change in the two files.
 */
(() => {
  // The two helpers of the bundler. They build the namespace objects, such as
  // `PhysiN.nd` and `PhysiN.slice`, from a table of getters.
  var __defProp = Object.defineProperty;
  /** Puts a getter on `target` for each name in `all`. */
  var __export = (target, all) => {
    for (var name in all)
      __defProp(target, name, { get: all[name], enumerable: true });
  };

  // src/nd/index.js
  var nd_exports = {};
  __export(nd_exports, {
    Body: () => Body,
    ConvexMesh: () => ConvexMesh,
    HalfSpace: () => HalfSpace,
    HyperBox: () => HyperBox,
    HyperSphere: () => HyperSphere,
    Torus: () => Torus,
    World: () => World,
    accumulateDrift: () => accumulateDrift,
    angularError: () => angularError,
    angularMass: () => angularMass,
    ballVolume: () => ballVolume,
    boxBoxAxis: () => boxBoxAxis,
    boxVertices: () => boxVertices,
    buildContactGraph: () => buildContactGraph,
    collide: () => collide,
    commutator: () => commutator,
    commutatorMatrix: () => commutatorMatrix,
    contractVecBi: () => contractVecBi,
    createConstraint: () => createConstraint,
    defaultParams: () => defaultParams,
    dims: () => dims,
    hingeAngle: () => hingeAngle,
    hyperBoxInertia: () => hyperBoxInertia,
    hyperBoxMesh: () => hyperBoxMesh,
    hyperSphereInertia: () => hyperSphereInertia,
    inertiaFromCovariance: () => inertiaFromCovariance,
    integratePositions: () => integratePositions,
    integrateVelocities: () => integrateVelocities,
    linalg: () => linalg_exports,
    massProperties: () => massProperties,
    mv: () => multivector_exports,
    orthoComplement: () => orthoComplement,
    prepareContact: () => prepareContact,
    prepareRows: () => prepareRows,
    refreshComplement: () => refreshComplement,
    resetConstraint: () => resetConstraint,
    rotor: () => rotor_exports,
    shockPropagation: () => shockPropagation,
    solveContact: () => solveContact,
    solveJoint: () => solveJoint,
    starMatrix: () => starMatrix,
    tangentBasis: () => tangentBasis,
    torusCovariance: () => torusCovariance,
    torusInertia: () => torusInertia,
    torusMesh: () => torusMesh,
    torusVolume: () => torusVolume,
    warmStart: () => warmStart,
    warmStartJoint: () => warmStartJoint,
    wedgeVec: () => wedgeVec
  });

  // src/nd/core/dims.js
  /**
   * The constant tables of the geometric algebra of `n` dimensions.
   *
   * A multivector has `N = 2^n` components. The index of a component is a bit
   * mask of the axes of its blade. Bit `i` is set when the axis `i` is in the
   * blade. Thus index 0 is the scalar, index `1 << i` is the axis `i`, and
   * index `(1 << i) | (1 << j)` is the plane of the axes `i` and `j`.
   *
   * Three lengths are in use through the library. Each array is flat.
   *   `n`  a vector: the position, the velocity, the force
   *   `k`  a bivector: the angular velocity, the torque, the momentum
   *   `r`  a rotor: the orientation
   * A matrix is flat and row major. See ND-PHYSICS.md, A1.
   *
   * @typedef {object} Dims
   * @property {number} n the count of dimensions. It is 2 or more.
   * @property {number} k the count of bivector components, `n (n-1) / 2`
   * @property {number} r the count of rotor components, `2^(n-1)`
   * @property {number} N the count of multivector components, `2^n`
   * @property {Int8Array} grade the grade of each blade, of length `N`
   * @property {Int8Array} reverseSign the sign that the reverse gives to a blade
   * @property {Int32Array} gpBlade the blade of a geometric product, `N` by `N`
   * @property {Int8Array} gpSign the sign of a geometric product, `N` by `N`
   * @property {Int32Array} vecBlade the blade of each axis, of length `n`
   * @property {number[][]} pairs the two axes of each bivector component
   * @property {number[]} biBlade the blade of each bivector component
   * @property {Int32Array} biOfBlade the bivector component of a blade, or -1
   * @property {number[]} evenBlade the blade of each rotor component
   * @property {Int32Array} slotOfBlade the rotor component of a blade, or -1
   * @property {Int32Array} evenGpSlot the component of a rotor product, `r` by `r`
   * @property {Int8Array} evenGpSign the sign of a rotor product, `r` by `r`
   * @property {number} scalarSlot the rotor component that holds the scalar
   * @property {Int32Array} biSlot the rotor component of a bivector component
   * @property {Float64Array} comm the commutator table, `k` by `k` by `k`
   * @property {Float64Array[]} eStar the star matrix of each axis, each `k` by `n`
   * @property {number} cXX a diagonal term of the canonical simplex covariance
   * @property {number} cXY an off diagonal term of the same covariance
   * @property {number} simplexVolumeDiv the divisor `n!` of a simplex volume
   * @property {number} simplexMomentDiv the divisor `(n+1)!` of a simplex moment
   */
  /** The tables of each `n` that `dims()` built before. */
  var cache = /* @__PURE__ */ new Map();
  /** The count of the bits that are set in `x`. This is the grade of a blade. */
  function popcount(x) {
    let c = 0;
    while (x !== 0) {
      x &= x - 1;
      c += 1;
    }
    return c;
  }
  /** The factorial `m!`. It is 1 when `m` is 0 or less. */
  function factorial(m) {
    let f = 1;
    for (let i = 2; i <= m; i += 1) f *= i;
    return f;
  }
  /**
   * The tables of the algebra of `n` dimensions. Build them one time, then
   * give the same object `D` to each function of the library.
   *
   * The sign of a product of two blades comes from the count of the swaps that
   * put the axes in order. The commutator table `comm` and the star matrices
   * `eStar` come from those signs. See ND-PHYSICS.md, A1, A2 and C4.
   *
   * The result is in a cache. Two calls with the same `n` give the same
   * object. Do not change the tables.
   *
   * @param {number} n the count of dimensions. It must be an integer, 2 or more.
   * @returns {Dims} the tables
   * @throws {Error} when `n` is not an integer of 2 or more
   */
  function dims(n) {
    if (!Number.isInteger(n) || n < 2) throw new Error("dims: n must be an integer of 2 or more");
    if (cache.has(n)) return cache.get(n);
    const N = 1 << n;
    // The grade of a blade is the count of its axes. The reverse turns the
    // order of the axes around. That gives the sign (-1)^(g (g-1) / 2).
    const grade = new Int8Array(N);
    const reverseSign = new Int8Array(N);
    for (let m = 0; m < N; m += 1) {
      const g = popcount(m);
      grade[m] = g;
      reverseSign[m] = g * (g - 1) / 2 & 1 ? -1 : 1;
    }
    // The geometric product of two blades. The blade of the result is the
    // exclusive or of the two masks. The sign comes from the count of the
    // swaps that put the axes back in order.
    const gpBlade = new Int32Array(N * N);
    const gpSign = new Int8Array(N * N);
    for (let a = 0; a < N; a += 1) {
      for (let b = 0; b < N; b += 1) {
        let m = a >> 1;
        let swaps = 0;
        while (m !== 0) {
          swaps += popcount(m & b);
          m >>= 1;
        }
        gpBlade[a * N + b] = a ^ b;
        gpSign[a * N + b] = swaps & 1 ? -1 : 1;
      }
    }
    const vecBlade = new Int32Array(n);
    for (let i = 0; i < n; i += 1) vecBlade[i] = 1 << i;
    // The bivector components, in lexicographic order of the axis pair. This
    // order is a free choice, but the library keeps it everywhere.
    const pairs = [];
    const biBlade = [];
    const biOfBlade = new Int32Array(N).fill(-1);
    for (let i = 0; i < n; i += 1) {
      for (let j = i + 1; j < n; j += 1) {
        biOfBlade[1 << i | 1 << j] = pairs.length;
        pairs.push([i, j]);
        biBlade.push(1 << i | 1 << j);
      }
    }
    const k = pairs.length;
    // A rotor holds only the blades of even grade: 0, 2, 4, and so on. Give
    // each one a slot. In 4 dimensions this gives 8 slots, not 4.
    const evenBlade = [];
    const slotOfBlade = new Int32Array(N).fill(-1);
    for (let m = 0; m < N; m += 1) {
      if ((grade[m] & 1) === 0) {
        slotOfBlade[m] = evenBlade.length;
        evenBlade.push(m);
      }
    }
    const r = evenBlade.length;
    // The geometric product of two rotors stays in the even sub algebra. Cut
    // the full `N` by `N` table down to `r` by `r`, and `rotorMul` is fast.
    const evenGpSlot = new Int32Array(r * r);
    const evenGpSign = new Int8Array(r * r);
    for (let a = 0; a < r; a += 1) {
      for (let b = 0; b < r; b += 1) {
        const idx = evenBlade[a] * N + evenBlade[b];
        evenGpSlot[a * r + b] = slotOfBlade[gpBlade[idx]];
        evenGpSign[a * r + b] = gpSign[idx];
      }
    }
    const scalarSlot = slotOfBlade[0];
    const biSlot = new Int32Array(k);
    for (let p = 0; p < k; p += 1) biSlot[p] = slotOfBlade[biBlade[p]];
    // The commutator `A x B = (A B - B A) / 2` of two bivectors. The result is
    // again a bivector. `comm[(p * k + q) * k + s]` is the part of the product
    // of the components `p` and `q` that goes to the component `s`. The
    // gyroscopic term uses this table. See ND-PHYSICS.md, A2 and A7.
    const comm = new Float64Array(k * k * k);
    for (let p = 0; p < k; p += 1) {
      for (let q = 0; q < k; q += 1) {
        const ab = biBlade[p] * N + biBlade[q];
        const ba = biBlade[q] * N + biBlade[p];
        if (gpBlade[ab] === gpBlade[ba]) {
          const s = biOfBlade[gpBlade[ab]];
          if (s >= 0) {
            comm[(p * k + q) * k + s] += (gpSign[ab] - gpSign[ba]) / 2;
          }
        }
      }
    }
    // The star matrix of each axis vector. `starMatrix` of any vector `r` is
    // the sum of `r[i] * eStar[i]`. `starProducts` uses these to build the
    // inertia tensor. See ND-PHYSICS.md, A3.
    const eStar = [];
    for (let i = 0; i < n; i += 1) {
      const S = new Float64Array(k * n);
      for (let p = 0; p < k; p += 1) {
        const [a, b] = pairs[p];
        if (a === i) S[p * n + b] += 1;
        if (b === i) S[p * n + a] -= 1;
      }
      eStar.push(S);
    }
    const D = {
      n,
      k,
      r,
      N,
      grade,
      reverseSign,
      gpBlade,
      gpSign,
      vecBlade,
      pairs,
      biBlade,
      biOfBlade,
      evenBlade,
      slotOfBlade,
      evenGpSlot,
      evenGpSign,
      scalarSlot,
      biSlot,
      comm,
      eStar,
      // Canonical simplex covariance values of Part A, step A8.
      cXX: 2 / factorial(n + 2),
      cXY: 1 / factorial(n + 2),
      simplexVolumeDiv: factorial(n),
      simplexMomentDiv: factorial(n + 1)
    };
    cache.set(n, D);
    return D;
  }

  // src/nd/core/linalg.js
  //
  // Small dense matrix work. The sides are `n` (2 to 6) or `k` (1 to 15), thus
  // a plain loop is faster than a library.
  //
  // The convention of the whole library: a matrix is a flat `Float64Array` in
  // row major order. The element of the row `i` and the column `j` of a matrix
  // of `cols` columns is at `A[i * cols + j]`. An array does not hold its own
  // shape. The caller gives the row count and the column count.
  //
  // Each function takes an optional `out` buffer. It writes into that buffer
  // and gives it back. Give a buffer to keep the step free of garbage.
  var linalg_exports = {};
  __export(linalg_exports, {
    jacobiEigen: () => jacobiEigen,
    matDet: () => matDet,
    matIdentity: () => matIdentity,
    matInverse: () => matInverse,
    matInverseSPD: () => matInverseSPD,
    matMul: () => matMul,
    matMulT: () => matMulT,
    matTVec: () => matTVec,
    matTranspose: () => matTranspose,
    matVec: () => matVec,
    matZero: () => matZero
  });
  /** A new matrix of `rows` by `cols`. All of its elements are zero. */
  function matZero(rows, cols) {
    return new Float64Array(rows * cols);
  }
  /** The identity matrix of `m` by `m`. It writes into `out` when you give it. */
  function matIdentity(m, out) {
    const A = out || matZero(m, m);
    A.fill(0);
    for (let i = 0; i < m; i += 1) A[i * m + i] = 1;
    return A;
  }
  /**
   * The product `A B`. `A` is `ra` by `ca`, and `B` is `ca` by `cb`.
   * The loop skips a zero element of `A`, because the inertia tensor and the
   * rotor matrices hold many zeros.
   * @returns {Float64Array} the product, `ra` by `cb`
   */
  function matMul(A, B, ra, ca, cb, out) {
    const C = out || matZero(ra, cb);
    C.fill(0);
    for (let i = 0; i < ra; i += 1) {
      for (let t = 0; t < ca; t += 1) {
        const a = A[i * ca + t];
        if (a === 0) continue;
        for (let j = 0; j < cb; j += 1) C[i * cb + j] += a * B[t * cb + j];
      }
    }
    return C;
  }
  /**
   * The product `A B^T`. `A` is `ra` by `ca`, and `B` is `rb` by `ca`.
   * With `matMul` this gives the change of frame `[R]2 I [R]2^T` of the
   * inertia tensor. See ND-PHYSICS.md, A6.
   * @returns {Float64Array} the product, `ra` by `rb`
   */
  function matMulT(A, B, ra, ca, rb, out) {
    const C = out || matZero(ra, rb);
    C.fill(0);
    for (let i = 0; i < ra; i += 1) {
      for (let j = 0; j < rb; j += 1) {
        let s = 0;
        for (let t = 0; t < ca; t += 1) s += A[i * ca + t] * B[j * ca + t];
        C[i * rb + j] = s;
      }
    }
    return C;
  }
  /** The product `A x`. `A` is `rows` by `cols`. It gives a vector of `rows`. */
  function matVec(A, x, rows, cols, out) {
    const y = out || new Float64Array(rows);
    for (let i = 0; i < rows; i += 1) {
      let s = 0;
      for (let j = 0; j < cols; j += 1) s += A[i * cols + j] * x[j];
      y[i] = s;
    }
    return y;
  }
  /**
   * The product `A^T x`. `A` is `rows` by `cols`, and `x` has `rows`. It gives
   * a vector of `cols`. This does not transpose `A` in memory.
   */
  function matTVec(A, x, rows, cols, out) {
    const y = out || new Float64Array(cols);
    y.fill(0);
    for (let i = 0; i < rows; i += 1) {
      const xi = x[i];
      if (xi === 0) continue;
      for (let j = 0; j < cols; j += 1) y[j] += A[i * cols + j] * xi;
    }
    return y;
  }
  /** The transpose of `A`. `A` is `rows` by `cols`, the result `cols` by `rows`. */
  function matTranspose(A, rows, cols, out) {
    const B = out || matZero(cols, rows);
    for (let i = 0; i < rows; i += 1) for (let j = 0; j < cols; j += 1) B[j * rows + i] = A[i * cols + j];
    return B;
  }
  /**
   * The determinant of the matrix `A` of `m` by `m`. It uses Gauss removal
   * with a partial pivot. It does not change `A`.
   *
   * `massProperties` and `orientCells` use the sign of the determinant to find
   * the direction of a simplex. See ND-PHYSICS.md, A8.
   *
   * @returns {number} the determinant. It is 0 when the matrix is singular.
   */
  function matDet(A, m) {
    const a = Float64Array.from(A);
    let det2 = 1;
    for (let c = 0; c < m; c += 1) {
      let piv = c;
      let best = Math.abs(a[c * m + c]);
      for (let i = c + 1; i < m; i += 1) {
        const v = Math.abs(a[i * m + c]);
        if (v > best) {
          best = v;
          piv = i;
        }
      }
      if (best === 0) return 0;
      if (piv !== c) {
        for (let j = 0; j < m; j += 1) {
          const t = a[c * m + j];
          a[c * m + j] = a[piv * m + j];
          a[piv * m + j] = t;
        }
        det2 = -det2;
      }
      const d = a[c * m + c];
      det2 *= d;
      for (let i = c + 1; i < m; i += 1) {
        const f = a[i * m + c] / d;
        if (f === 0) continue;
        for (let j = c; j < m; j += 1) a[i * m + j] -= f * a[c * m + j];
      }
    }
    return det2;
  }
  /**
   * The inverse of the matrix `A` of `m` by `m`. It uses Gauss-Jordan removal
   * with a partial pivot. It does not change `A`.
   *
   * Use `matInverseSPD` for an inertia tensor. It is symmetrical and positive
   * definite, and Cholesky is faster and more stable.
   *
   * @throws {Error} when the matrix is singular
   */
  function matInverse(A, m, out) {
    const a = Float64Array.from(A);
    const inv = matIdentity(m, out);
    for (let c = 0; c < m; c += 1) {
      let piv = c;
      let best = Math.abs(a[c * m + c]);
      for (let i = c + 1; i < m; i += 1) {
        const v = Math.abs(a[i * m + c]);
        if (v > best) {
          best = v;
          piv = i;
        }
      }
      if (best < 1e-300) throw new Error("matInverse: the matrix is singular");
      if (piv !== c) {
        for (let j = 0; j < m; j += 1) {
          let t = a[c * m + j];
          a[c * m + j] = a[piv * m + j];
          a[piv * m + j] = t;
          t = inv[c * m + j];
          inv[c * m + j] = inv[piv * m + j];
          inv[piv * m + j] = t;
        }
      }
      const d = a[c * m + c];
      for (let j = 0; j < m; j += 1) {
        a[c * m + j] /= d;
        inv[c * m + j] /= d;
      }
      for (let i = 0; i < m; i += 1) {
        if (i === c) continue;
        const f = a[i * m + c];
        if (f === 0) continue;
        for (let j = 0; j < m; j += 1) {
          a[i * m + j] -= f * a[c * m + j];
          inv[i * m + j] -= f * inv[c * m + j];
        }
      }
    }
    return inv;
  }
  /**
   * The inverse of a symmetrical positive definite matrix of `m` by `m`. This
   * is the usual condition of an inertia tensor.
   *
   * The three steps: cut `A` into `L L^T` with Cholesky, invert the lower
   * triangle `L`, then give `A^-1 = L^-T L^-1`. The result is symmetrical.
   *
   * A body that is flat in one axis can give a matrix that is not positive
   * definite. In that condition the function falls back to `matInverse`.
   */
  function matInverseSPD(A, m, out) {
    const L = matZero(m, m);
    for (let i = 0; i < m; i += 1) {
      for (let j = 0; j <= i; j += 1) {
        let s = A[i * m + j];
        for (let t = 0; t < j; t += 1) s -= L[i * m + t] * L[j * m + t];
        if (i === j) {
          if (s <= 0) return matInverse(A, m, out);
          L[i * m + j] = Math.sqrt(s);
        } else {
          L[i * m + j] = s / L[j * m + j];
        }
      }
    }
    const Li = matZero(m, m);
    for (let i = 0; i < m; i += 1) {
      Li[i * m + i] = 1 / L[i * m + i];
      for (let j = 0; j < i; j += 1) {
        let s = 0;
        for (let t = j; t < i; t += 1) s += L[i * m + t] * Li[t * m + j];
        Li[i * m + j] = -s * Li[i * m + i];
      }
    }
    const inv = out || matZero(m, m);
    inv.fill(0);
    for (let i = 0; i < m; i += 1) {
      for (let j = 0; j <= i; j += 1) {
        let s = 0;
        for (let t = i; t < m; t += 1) s += Li[t * m + i] * Li[t * m + j];
        inv[i * m + j] = s;
        inv[j * m + i] = s;
      }
    }
    return inv;
  }
  /**
   * The eigenvalues and the eigenvectors of a symmetrical matrix of `m` by
   * `m`. It uses the cyclic Jacobi method. Each turn of the loop makes one
   * pair of off diagonal elements zero. It stops when the sum of the squares
   * of the off diagonal elements is very small, or after `sweeps` turns.
   *
   * This is a tool for a test and for analysis. The step does not use it.
   *
   * @returns {{values: Float64Array, vectors: Float64Array}} the eigenvalues,
   *   of length `m`, and the eigenvectors as the columns of an `m` by `m`
   *   matrix. The values are not in order.
   */
  function jacobiEigen(A, m, sweeps = 60) {
    const a = Float64Array.from(A);
    const V = matIdentity(m);
    for (let s = 0; s < sweeps; s += 1) {
      let off = 0;
      for (let i = 0; i < m; i += 1) for (let j = i + 1; j < m; j += 1) off += a[i * m + j] * a[i * m + j];
      if (off < 1e-30) break;
      for (let p = 0; p < m; p += 1) {
        for (let q = p + 1; q < m; q += 1) {
          const apq = a[p * m + q];
          if (Math.abs(apq) < 1e-300) continue;
          const theta = (a[q * m + q] - a[p * m + p]) / (2 * apq);
          const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          const c = 1 / Math.sqrt(t * t + 1);
          const sn = t * c;
          for (let i = 0; i < m; i += 1) {
            const aip = a[i * m + p];
            const aiq = a[i * m + q];
            a[i * m + p] = c * aip - sn * aiq;
            a[i * m + q] = sn * aip + c * aiq;
          }
          for (let i = 0; i < m; i += 1) {
            const api = a[p * m + i];
            const aqi = a[q * m + i];
            a[p * m + i] = c * api - sn * aqi;
            a[q * m + i] = sn * api + c * aqi;
          }
          for (let i = 0; i < m; i += 1) {
            const vip = V[i * m + p];
            const viq = V[i * m + q];
            V[i * m + p] = c * vip - sn * viq;
            V[i * m + q] = sn * vip + c * viq;
          }
        }
      }
    }
    const values = new Float64Array(m);
    for (let i = 0; i < m; i += 1) values[i] = a[i * m + i];
    return { values, vectors: V };
  }

  // src/nd/algebra/multivector.js
  //
  // The full multivector of `N = 2^n` components. The index of a component is
  // the bit mask of its blade, as `dims()` gives it.
  //
  // A multivector is large. The step does not use one. The three small types
  // (vector, bivector, rotor) do the work of each frame. Use a multivector
  // only to build a table, or for an operation that mixes the grades, such as
  // the wedge product of the separating axis test.
  //
  // The `mvFrom*` and `mvTo*` functions move between the small types and the
  // full multivector. See ND-PHYSICS.md, A1 and A2.
  var multivector_exports = {};
  __export(multivector_exports, {
    mvContract: () => mvContract,
    mvDual: () => mvDual,
    mvFromBivector: () => mvFromBivector,
    mvFromRotor: () => mvFromRotor,
    mvFromVector: () => mvFromVector,
    mvGp: () => mvGp,
    mvGrade: () => mvGrade,
    mvNorm: () => mvNorm,
    mvReverse: () => mvReverse,
    mvToBivector: () => mvToBivector,
    mvToRotor: () => mvToRotor,
    mvToVector: () => mvToVector,
    mvWedge: () => mvWedge,
    mvZero: () => mvZero
  });
  /** A new multivector of `N` components. All of them are zero. */
  function mvZero(D) {
    return new Float64Array(D.N);
  }
  /** Puts the vector `v`, of length `n`, into the grade 1 blades. */
  function mvFromVector(D, v, out) {
    const M = out || mvZero(D);
    M.fill(0);
    for (let i = 0; i < D.n; i += 1) M[D.vecBlade[i]] = v[i];
    return M;
  }
  /** Puts the bivector `B`, of length `k`, into the grade 2 blades. */
  function mvFromBivector(D, B, out) {
    const M = out || mvZero(D);
    M.fill(0);
    for (let p = 0; p < D.k; p += 1) M[D.biBlade[p]] = B[p];
    return M;
  }
  /** Puts the rotor `R`, of length `r`, into the blades of even grade. */
  function mvFromRotor(D, R, out) {
    const M = out || mvZero(D);
    M.fill(0);
    for (let s = 0; s < D.r; s += 1) M[D.evenBlade[s]] = R[s];
    return M;
  }
  /** Takes the grade 1 part of `M` out, as a vector of length `n`. */
  function mvToVector(D, M, out) {
    const v = out || new Float64Array(D.n);
    for (let i = 0; i < D.n; i += 1) v[i] = M[D.vecBlade[i]];
    return v;
  }
  /** Takes the grade 2 part of `M` out, as a bivector of length `k`. */
  function mvToBivector(D, M, out) {
    const B = out || new Float64Array(D.k);
    for (let p = 0; p < D.k; p += 1) B[p] = M[D.biBlade[p]];
    return B;
  }
  /** Takes the even part of `M` out, as a rotor of length `r`. */
  function mvToRotor(D, M, out) {
    const R = out || new Float64Array(D.r);
    for (let s = 0; s < D.r; s += 1) R[s] = M[D.evenBlade[s]];
    return R;
  }
  /**
   * The geometric product `A B` of two multivectors. This is the basic
   * product of the algebra. For two vectors it gives `a b = a . b + a ^ b`,
   * thus the scalar part and the bivector part together.
   *
   * The rotor sandwich `R x R~` uses two of these products.
   *
   * @param {Dims} D the tables from `dims(n)`
   * @param {Float64Array} A a multivector of `N`
   * @param {Float64Array} B a multivector of `N`
   * @param {Float64Array} [out] a buffer of `N`. It must not be `A` or `B`.
   * @returns {Float64Array} the product
   * @throws {Error} when `out` is one of the inputs
   */
  function mvGp(D, A, B, out) {
    const C = out || mvZero(D);
    if (C === A || C === B) throw new Error("mvGp: the output must not be an input");
    C.fill(0);
    const N = D.N;
    for (let a = 0; a < N; a += 1) {
      const va = A[a];
      if (va === 0) continue;
      for (let b = 0; b < N; b += 1) {
        const vb = B[b];
        if (vb === 0) continue;
        const idx = a * N + b;
        C[D.gpBlade[idx]] += D.gpSign[idx] * va * vb;
      }
    }
    return C;
  }
  /**
   * The exterior product `A ^ B`. It is the geometric product without the
   * terms that hold a common axis, thus the test `(a & b) !== 0`.
   *
   * The wedge of two vectors gives the plane through them. `boxBoxAxis` wedges
   * `n - 1` axes together, then takes the dual, to build a separating axis.
   * See ND-PHYSICS.md, A10.
   *
   * @throws {Error} when `out` is one of the inputs
   */
  function mvWedge(D, A, B, out) {
    const C = out || mvZero(D);
    if (C === A || C === B) throw new Error("mvWedge: the output must not be an input");
    C.fill(0);
    const N = D.N;
    for (let a = 0; a < N; a += 1) {
      const va = A[a];
      if (va === 0) continue;
      for (let b = 0; b < N; b += 1) {
        const vb = B[b];
        if (vb === 0) continue;
        if ((a & b) !== 0) continue;
        const idx = a * N + b;
        C[a ^ b] += D.gpSign[idx] * va * vb;
      }
    }
    return C;
  }
  /**
   * The left contraction `A . B`. It keeps only the terms in which the blade
   * of `A` is fully inside the blade of `B`, thus the test `(a & b) === a`.
   * The grade of the result is the grade of `B` less the grade of `A`.
   *
   * The velocity of a point uses the contraction `r . w` of a vector and a
   * bivector. `contractVecBi` does that one operation, and it is faster.
   *
   * @throws {Error} when `out` is one of the inputs
   */
  function mvContract(D, A, B, out) {
    const C = out || mvZero(D);
    if (C === A || C === B) throw new Error("mvContract: the output must not be an input");
    C.fill(0);
    const N = D.N;
    for (let a = 0; a < N; a += 1) {
      const va = A[a];
      if (va === 0) continue;
      for (let b = 0; b < N; b += 1) {
        const vb = B[b];
        if (vb === 0) continue;
        if ((a & b) !== a) continue;
        const idx = a * N + b;
        C[a ^ b] += D.gpSign[idx] * va * vb;
      }
    }
    return C;
  }
  /**
   * The reverse `A~`. It turns the order of the axes of each blade around.
   * That gives the sign `(-1)^(g (g-1) / 2)` to a blade of the grade `g`.
   * For a unit rotor the reverse is the inverse. See ND-PHYSICS.md, A2.
   */
  function mvReverse(D, A, out) {
    const C = out || mvZero(D);
    for (let m = 0; m < D.N; m += 1) C[m] = D.reverseSign[m] * A[m];
    return C;
  }
  /** Keeps the part of `A` of the grade `g`. All other components become zero. */
  function mvGrade(D, A, g, out) {
    const C = out || mvZero(D);
    C.fill(0);
    for (let m = 0; m < D.N; m += 1) if (D.grade[m] === g) C[m] = A[m];
    return C;
  }
  /**
   * The dual `A I^-1`. `I` is the pseudoscalar, thus the blade of all `n`
   * axes, at the index `N - 1`. The dual changes a blade of the grade `g`
   * into a blade of the grade `n - g`.
   *
   * `boxBoxAxis` takes the dual of a blade of the grade `n - 1` to get the
   * vector normal to it, and that vector is a separating axis.
   */
  function mvDual(D, A, out) {
    const C = out || mvZero(D);
    C.fill(0);
    const I = D.N - 1;
    const sInv = D.reverseSign[I];
    const N = D.N;
    for (let a = 0; a < N; a += 1) {
      const va = A[a];
      if (va === 0) continue;
      const idx = a * N + I;
      C[D.gpBlade[idx]] += D.gpSign[idx] * sInv * va;
    }
    return C;
  }
  /** The square root of the sum of the squares of all `N` components. */
  function mvNorm(D, A) {
    let s = 0;
    for (let m = 0; m < D.N; m += 1) s += A[m] * A[m];
    return Math.sqrt(s);
  }

  // src/nd/algebra/rotor.js
  //
  // The orientation of a body is a rotor. A rotor is a member of the even sub
  // algebra: it holds only the blades of the grade 0, 2, 4, and so on. It is a
  // flat `Float64Array` of `r = 2^(n-1)` components, in the order of
  // `D.evenBlade`. `D.scalarSlot` is the component of the scalar.
  //
  // CAUTION: a rotor is not a quaternion. In 3 dimensions it has 4 components,
  // and it acts as a quaternion does. In 4 dimensions it has 8 components: 1
  // scalar, 6 bivector and 1 quadvector. The quadvector part is not zero when
  // the body turns in two planes at the same time. Do not take a quaternion
  // class and add components to it. See ND-PHYSICS.md, A1.
  //
  // A rotor turns a vector or a blade with the sandwich product `x' = R x R~`.
  // A unit rotor has `R R~ = 1`, thus the reverse is the inverse. Rounding
  // makes a rotor lose its unit length. `rotorDefect` measures that loss, and
  // `rotorCorrect` repairs it. See ND-PHYSICS.md, B2.
  var rotor_exports = {};
  __export(rotor_exports, {
    rotorAddScaled: () => rotorAddScaled,
    rotorApplyBivector: () => rotorApplyBivector,
    rotorApplyVector: () => rotorApplyVector,
    rotorApplyVectorInverse: () => rotorApplyVectorInverse,
    rotorBetweenVectors: () => rotorBetweenVectors,
    rotorBivectorMatrix: () => rotorBivectorMatrix,
    rotorCorrect: () => rotorCorrect,
    rotorDefect: () => rotorDefect,
    rotorExp: () => rotorExp,
    rotorFromBivectorAngle: () => rotorFromBivectorAngle,
    rotorFromPlane: () => rotorFromPlane,
    rotorIdentity: () => rotorIdentity,
    rotorMatrix: () => rotorMatrix,
    rotorMul: () => rotorMul,
    rotorNorm: () => rotorNorm,
    rotorReverse: () => rotorReverse,
    rotorScale: () => rotorScale
  });
  /**
   * A length below this one counts as zero. This is not a tolerance of the
   * world, thus it is not in `defaultParams`: it says where the plane of two
   * vectors stops to have a meaning, and it comes from the accuracy of the
   * double and not from a choice of the user. See design rule 2.
   */
  var ROTOR_ZERO = 1e-12;
  /** The rotor that turns nothing. Its scalar component is 1. */
  function rotorIdentity(D, out) {
    const R = out || new Float64Array(D.r);
    R.fill(0);
    R[D.scalarSlot] = 1;
    return R;
  }
  /**
   * The geometric product `A B` of two rotors. The result is the rotor that
   * puts `B` first and then `A`.
   *
   * This uses the small `r` by `r` table `D.evenGpSlot`, and not the full `N`
   * by `N` table. Thus it does not build a multivector.
   *
   * @param {Float64Array} [out] a buffer of `r`. It must not be `A` or `B`.
   * @throws {Error} when `out` is one of the inputs
   */
  function rotorMul(D, A, B, out) {
    const r = D.r;
    const C = out || new Float64Array(r);
    if (C === A || C === B) throw new Error("rotorMul: the output must not be an input");
    C.fill(0);
    for (let a = 0; a < r; a += 1) {
      const va = A[a];
      if (va === 0) continue;
      for (let b = 0; b < r; b += 1) {
        const vb = B[b];
        if (vb === 0) continue;
        const idx = a * r + b;
        C[D.evenGpSlot[idx]] += D.evenGpSign[idx] * va * vb;
      }
    }
    return C;
  }
  /** The reverse `R~`. For a unit rotor this is the inverse, thus the turn back. */
  function rotorReverse(D, R, out) {
    const C = out || new Float64Array(D.r);
    for (let s = 0; s < D.r; s += 1) C[s] = D.reverseSign[D.evenBlade[s]] * R[s];
    return C;
  }
  /** The length of `R`. A rotor that turns a body has the length 1. */
  function rotorNorm(D, R) {
    let s = 0;
    for (let i = 0; i < D.r; i += 1) s += R[i] * R[i];
    return Math.sqrt(s);
  }
  /** Multiplies each component of `R` by `f`. It gives a new rotor. */
  function rotorScale(D, R, f, out) {
    const C = out || new Float64Array(D.r);
    for (let i = 0; i < D.r; i += 1) C[i] = R[i] * f;
    return C;
  }
  /** Adds `f` times `S` to `R`, in place. It gives `R`. */
  function rotorAddScaled(D, R, S, f) {
    for (let i = 0; i < D.r; i += 1) R[i] += S[i] * f;
    return R;
  }
  /**
   * The rotor of a turn of `angle` radians in the plane of the axes `i` and
   * `j`. The turn goes from the axis `i` to the axis `j`.
   *
   * The half angle is correct: the sandwich product applies the rotor two
   * times, thus the body turns through the full angle.
   *
   * @param {number} i the first axis, 0 to `n - 1`
   * @param {number} j the second axis. It must not be `i`.
   * @param {number} angle the angle in radians
   * @throws {Error} when `i` and `j` are the same axis
   */
  function rotorFromPlane(D, i, j, angle, out) {
    const R = rotorIdentity(D, out);
    const p = D.biOfBlade[1 << i | 1 << j];
    if (p < 0) throw new Error("rotorFromPlane: i and j must be different");
    const sign = i < j ? 1 : -1;
    R[D.scalarSlot] = Math.cos(angle / 2);
    R[D.biSlot[p]] = -sign * Math.sin(angle / 2);
    return R;
  }
  /**
   * The exponential `exp(B)` of a bivector. This is the general way to make a
   * rotor from a plane and an angle. In 4 dimensions `B` can hold two planes
   * at the same time, and one series solves both.
   *
   * The method is scale and square:
   *   1. Halve `B` until its length is 0.25 or less.
   *   2. Sum 16 terms of the Taylor series of the exponential.
   *   3. Square the result one time for each halving.
   * A small input keeps the series short and the error low.
   *
   * @param {Float64Array} B a bivector of length `k`
   * @returns {Float64Array} the rotor, of length `r`
   */
  function rotorExp(D, B, out) {
    let mag = 0;
    for (let p = 0; p < D.k; p += 1) mag += B[p] * B[p];
    mag = Math.sqrt(mag);
    let steps = 0;
    let scale = 1;
    while (mag * scale > 0.25) {
      scale /= 2;
      steps += 1;
    }
    const X = new Float64Array(D.r);
    for (let p = 0; p < D.k; p += 1) X[D.biSlot[p]] = B[p] * scale;
    let term = rotorIdentity(D);
    let acc = rotorIdentity(D);
    let tmp = new Float64Array(D.r);
    for (let m = 1; m <= 16; m += 1) {
      rotorMul(D, term, X, tmp);
      for (let i = 0; i < D.r; i += 1) tmp[i] /= m;
      const t = term;
      term = tmp;
      tmp = t;
      for (let i = 0; i < D.r; i += 1) acc[i] += term[i];
    }
    for (let s = 0; s < steps; s += 1) {
      const sq = rotorMul(D, acc, acc);
      acc = sq;
    }
    const R = out || new Float64Array(D.r);
    R.set(acc);
    return R;
  }
  /**
   * The rotor of a turn of `angle` radians in the plane of the bivector `B`.
   * It makes `B` a unit bivector first, thus the length of `B` has no effect.
   * It gives the identity when `B` is zero.
   *
   * `angle` is the turn of the BODY, as in `rotorFromPlane`. The exponent
   * holds the half angle, because the sandwich product applies the rotor two
   * times, and it holds the minus of `dR/dt = -(1/2) w R`. Thus
   * `rotorFromBivectorAngle(D, e_i ^ e_j, angle)` and
   * `rotorFromPlane(D, i, j, angle)` give the same rotor.
   */
  function rotorFromBivectorAngle(D, B, angle, out) {
    let mag = 0;
    for (let p = 0; p < D.k; p += 1) mag += B[p] * B[p];
    mag = Math.sqrt(mag);
    if (mag < 1e-300) return rotorIdentity(D, out);
    const S = new Float64Array(D.k);
    for (let p = 0; p < D.k; p += 1) S[p] = -B[p] / mag * (angle / 2);
    return rotorExp(D, S, out);
  }
  /** The scratch multivectors of each `D`, for the sandwich products. */
  var scratch = /* @__PURE__ */ new WeakMap();
  /** Four scratch multivectors of the algebra `D`. They keep the step free of garbage. */
  function pad(D) {
    let s = scratch.get(D);
    if (!s) {
      s = { a: mvZero(D), b: mvZero(D), c: mvZero(D), d: mvZero(D) };
      scratch.set(D, s);
    }
    return s;
  }
  /**
   * Turns the vector `x` with the sandwich product `R x R~`. It gives the
   * vector in the world frame. See ND-PHYSICS.md, A2.
   *
   * For many vectors, build the matrix with `rotorMatrix` one time, then use
   * `matVec`. That is much faster than one sandwich for each vector.
   *
   * @param {Float64Array} x a vector of length `n`
   * @returns {Float64Array} the vector after the turn, of length `n`
   */
  function rotorApplyVector(D, R, x, out) {
    const s = pad(D);
    mvFromRotor(D, R, s.a);
    mvFromVector(D, x, s.b);
    mvGp(D, s.a, s.b, s.c);
    for (let m = 0; m < D.N; m += 1) s.b[m] = D.reverseSign[m] * s.a[m];
    mvGp(D, s.c, s.b, s.d);
    return mvToVector(D, s.d, out);
  }
  /**
   * Turns the bivector `B` with the same sandwich product `R B R~`. The
   * angular velocity, the torque and the angular momentum are bivectors, thus
   * they change frame with this function and not with `rotorApplyVector`.
   * `rotorBivectorMatrix` builds the matrix form.
   */
  function rotorApplyBivector(D, R, B, out) {
    const s = pad(D);
    mvFromRotor(D, R, s.a);
    mvFromBivector(D, B, s.b);
    mvGp(D, s.a, s.b, s.c);
    for (let m = 0; m < D.N; m += 1) s.b[m] = D.reverseSign[m] * s.a[m];
    mvGp(D, s.c, s.b, s.d);
    return mvToBivector(D, s.d, out);
  }
  /** Turns the vector `x` back, with `R~ x R`. It goes from world to body. */
  function rotorApplyVectorInverse(D, R, x, out) {
    const Rr = rotorReverse(D, R);
    return rotorApplyVector(D, Rr, x, out);
  }
  /**
   * The rotation matrix of the rotor `R`, of `n` by `n`. Its column `i` is the
   * axis `i` after the turn. `Body.Rm` holds this matrix, and the body builds
   * it one time in each step.
   *
   * `M x` changes a vector from the body frame to the world frame. `M^T x`
   * changes it back, because the matrix is orthogonal.
   *
   * @returns {Float64Array} the matrix, `n` by `n`, row major
   */
  function rotorMatrix(D, R, out) {
    const n = D.n;
    const M = out || new Float64Array(n * n);
    const e = new Float64Array(n);
    const f = new Float64Array(n);
    for (let i = 0; i < n; i += 1) {
      e.fill(0);
      e[i] = 1;
      rotorApplyVector(D, R, e, f);
      for (let j = 0; j < n; j += 1) M[j * n + i] = f[j];
    }
    return M;
  }
  /**
   * The matrix `[R]2` of the rotor `R`, of `k` by `k`. It does to a bivector
   * what `rotorMatrix` does to a vector: `[R]2 B = R B R~`.
   *
   * `Body.R2` holds this matrix. With it, the change of frame of the inertia
   * tensor is a matrix product `[R]2 I [R]2^T`, and not an algebra product.
   * See ND-PHYSICS.md, A4.
   *
   * @returns {Float64Array} the matrix, `k` by `k`, row major
   */
  function rotorBivectorMatrix(D, R, out) {
    const k = D.k;
    const M = out || new Float64Array(k * k);
    const e = new Float64Array(k);
    const f = new Float64Array(k);
    for (let p = 0; p < k; p += 1) {
      e.fill(0);
      e[p] = 1;
      rotorApplyBivector(D, R, e, f);
      for (let q = 0; q < k; q += 1) M[q * k + p] = f[q];
    }
    return M;
  }
  /**
   * The shortest rotor that turns the unit vector `a` onto the unit vector
   * `b`. The turn stays in the plane of `a` and `b`.
   *
   * The direct form is `1 + b a`, made a unit rotor. It is the half angle
   * form. It is accurate while `a . b` is not near -1, and it loses ALL of
   * its accuracy at -1, because `1 + a . b` and `b ^ a` are then two small
   * differences of large numbers. Thus there are three conditions:
   *
   *   - `a . b` is 0 or more: the direct form.
   *   - `a . b` is less than 0: two turns. `P` takes `a` to `-a` in the plane
   *     of `a` and `b`, and the direct form takes `-a` to `b`. The dot of the
   *     second turn is `-a . b`, thus more than 0, thus the form is accurate.
   *     The two turns are in the same plane, thus the answer is still the
   *     shortest turn.
   *   - `b` is `-a`: the plane has no meaning. It makes a vector normal to
   *     `a`, then it turns through pi radians.
   *
   * The last step makes the length of the rotor 1 with the true euclidean
   * norm, and not with the divisor `sqrt(2 (1 + a . b))` of the algebra. The
   * rotor is a scalar and ONE 2-blade, thus a euclidean length of 1 makes it
   * an exact rotor, whatever the rounding did.
   *
   * `rotorCorrect` does NOT use this function. It builds its own plane,
   * because it must also keep that plane normal to the axes that are already
   * in their place. See the note there.
   *
   * @param {Float64Array} a a unit vector of length `n`
   * @param {Float64Array} b a unit vector of length `n`
   */
  function rotorBetweenVectors(D, a, b, out) {
    const n = D.n;
    let dot = 0;
    for (let i = 0; i < n; i += 1) dot += a[i] * b[i];
    if (dot < 0) {
      const t = new Float64Array(n);
      let ln = 0;
      for (let i = 0; i < n; i += 1) {
        t[i] = b[i] - dot * a[i];
        ln += t[i] * t[i];
      }
      ln = Math.sqrt(ln);
      if (ln < ROTOR_ZERO) {
        let best = 0;
        for (let i = 1; i < n; i += 1) if (Math.abs(a[i]) < Math.abs(a[best])) best = i;
        t.fill(0);
        t[best] = 1;
        let d2 = 0;
        for (let i = 0; i < n; i += 1) d2 += t[i] * a[i];
        for (let i = 0; i < n; i += 1) t[i] -= d2 * a[i];
        ln = 0;
        for (let i = 0; i < n; i += 1) ln += t[i] * t[i];
        ln = Math.sqrt(ln);
      }
      for (let i = 0; i < n; i += 1) t[i] /= ln;
      const B = new Float64Array(D.k);
      for (let p = 0; p < D.k; p += 1) {
        const [i, j] = D.pairs[p];
        B[p] = a[i] * t[j] - a[j] * t[i];
      }
      const P = rotorFromBivectorAngle(D, B, Math.PI);
      const back = new Float64Array(n);
      for (let i = 0; i < n; i += 1) back[i] = -a[i];
      return rotorMul(D, rotorBetweenVectors(D, back, b), P, out);
    }
    const R = out || new Float64Array(D.r);
    R.fill(0);
    R[D.scalarSlot] = 1 + dot;
    for (let p = 0; p < D.k; p += 1) {
      const [i, j] = D.pairs[p];
      R[D.biSlot[p]] = b[i] * a[j] - b[j] * a[i];
    }
    let len = 0;
    for (let i = 0; i < D.r; i += 1) len += R[i] * R[i];
    len = Math.sqrt(len);
    for (let i = 0; i < D.r; i += 1) R[i] /= len;
    return R;
  }
  /**
   * Repairs a rotor that rounding made bad. This function is critical in 4
   * dimensions and more. See ND-PHYSICS.md, B2, and README, section 10.
   *
   * A 3D engine keeps a quaternion good with a divide by its length. That is
   * not sufficient here. A rotor of 8 components can hold an error that the
   * length does not show: the sandwich product then stops being a rotation,
   * and the body becomes larger or thinner as it turns.
   *
   * The method has four steps:
   *   1. Build the matrix `F` of the rotor with `rotorMatrix`.
   *   2. Make the columns of `F` orthogonal and of the length 1, with the
   *      Gram-Schmidt method. Now `F` is a true rotation matrix.
   *   3. Build a new rotor from `F`. Take each axis `c` in turn, and turn the
   *      column `c` of a work frame `cur` on to the column `c` of `F`. The
   *      plane holds `p`, the column of `cur`, and `q`, the part of the column
   *      of `F` that is normal to `p`. The angle is `atan2(|q|, p . tgt)`.
   *      Multiply the rotors, and turn `cur` with the SAME rotor, thus `cur`
   *      is always the matrix of the product.
   *   4. Give the new rotor the same sign as the old one. `R` and `-R` are the
   *      same rotation, and a change of the sign would make the body jump.
   *
   * A HALF TURN IS THE HARD CONDITION. A turn must not move an axis that is
   * already in its place. `q` is a small difference of large numbers when the
   * angle comes near pi: at an angle of `pi - 1e-11` the part of `q` along a
   * placed column can grow to 1e-4, and a turn of almost pi radians then moves
   * that column by two times as much. Thus step 3 makes `q` normal to the
   * columns `0` to `c - 1` of `F` again, before it makes the length of `q` 1.
   *
   * When `q` is zero, the column of `F` is the column of `cur`, or its
   * opposite. The plane then has no meaning. If the two agree, there is
   * nothing to do. If they are opposite, take the plane of `p` and the column
   * `c + 1` of `cur`: that column is normal to `p` and to every placed column,
   * by construction. The column `c + 1` is always there, because the loop
   * stops at `n - 2`.
   *
   * It gives the identity when the matrix is degenerate.
   *
   * @param {Float64Array} R the rotor to repair. It does not change.
   * @returns {Float64Array} the repaired rotor, of length `r`
   */
  function rotorCorrect(D, R, out) {
    const n = D.n;
    const F = rotorMatrix(D, R);
    for (let c = 0; c < n; c += 1) {
      for (let p2 = 0; p2 < c; p2 += 1) {
        let d = 0;
        for (let i = 0; i < n; i += 1) d += F[i * n + c] * F[i * n + p2];
        for (let i = 0; i < n; i += 1) F[i * n + c] -= d * F[i * n + p2];
      }
      let ln = 0;
      for (let i = 0; i < n; i += 1) ln += F[i * n + c] * F[i * n + c];
      ln = Math.sqrt(ln);
      if (ln < ROTOR_ZERO) {
        return rotorIdentity(D, out);
      }
      for (let i = 0; i < n; i += 1) F[i * n + c] /= ln;
    }
    let acc = rotorIdentity(D);
    const cur = new Float64Array(n * n);
    for (let i = 0; i < n; i += 1) cur[i * n + i] = 1;
    const p = new Float64Array(n);
    const q = new Float64Array(n);
    const tgt = new Float64Array(n);
    const col = new Float64Array(n);
    const turned = new Float64Array(n);
    for (let c = 0; c < n - 1; c += 1) {
      for (let i = 0; i < n; i += 1) {
        p[i] = cur[i * n + c];
        tgt[i] = F[i * n + c];
      }
      let dot2 = 0;
      for (let i = 0; i < n; i += 1) dot2 += p[i] * tgt[i];
      for (let i = 0; i < n; i += 1) q[i] = tgt[i] - dot2 * p[i];
      for (let e = 0; e < c; e += 1) {
        let d = 0;
        for (let i = 0; i < n; i += 1) d += q[i] * F[i * n + e];
        for (let i = 0; i < n; i += 1) q[i] -= d * F[i * n + e];
      }
      let ln = 0;
      for (let i = 0; i < n; i += 1) ln += q[i] * q[i];
      ln = Math.sqrt(ln);
      let angle;
      if (ln < ROTOR_ZERO) {
        if (dot2 > 0) continue;
        for (let i = 0; i < n; i += 1) q[i] = cur[i * n + c + 1];
        angle = Math.PI;
      } else {
        for (let i = 0; i < n; i += 1) q[i] /= ln;
        angle = Math.atan2(ln, dot2);
      }
      const B = new Float64Array(D.k);
      for (let t = 0; t < D.k; t += 1) {
        const [i, j] = D.pairs[t];
        B[t] = p[i] * q[j] - p[j] * q[i];
      }
      const S = rotorFromBivectorAngle(D, B, angle);
      acc = rotorMul(D, S, acc);
      for (let cc = 0; cc < n; cc += 1) {
        for (let i = 0; i < n; i += 1) col[i] = cur[i * n + cc];
        rotorApplyVector(D, S, col, turned);
        for (let i = 0; i < n; i += 1) cur[i * n + cc] = turned[i];
      }
    }
    let dot = 0;
    for (let i = 0; i < D.r; i += 1) dot += acc[i] * R[i];
    const C = out || new Float64Array(D.r);
    const s = dot < 0 ? -1 : 1;
    for (let i = 0; i < D.r; i += 1) C[i] = acc[i] * s;
    return C;
  }
  /**
   * How far the rotor `R` is from a good rotor. A good rotor has `R R~ = 1`.
   * The function multiplies `R` by its reverse, then it adds the error of the
   * scalar part to the size of all other parts.
   *
   * `integratePositions` calls this after each step. It calls `rotorCorrect`
   * only when the result is more than `params.rotorTolerance`. Thus the
   * repair, which costs much, does not run in each step.
   *
   * @returns {number} the error. It is 0 for a good rotor.
   */
  function rotorDefect(D, R) {
    const Rr = rotorReverse(D, R);
    const P = rotorMul(D, R, Rr);
    let e = Math.abs(P[D.scalarSlot] - 1);
    for (let i = 0; i < D.r; i += 1) if (i !== D.scalarSlot) e += Math.abs(P[i]);
    return e;
  }

  // src/nd/algebra/star.js
  //
  // The star matrix is the key tool of the library. It connects the vector
  // world (position, force, impulse) to the bivector world (torque, angular
  // velocity, angular momentum). Almost all of the mechanics uses it.
  //
  // In 3 dimensions the star matrix is the cross product matrix `[r]x`. In 4
  // dimensions and more there is no cross product, but the star matrix stays
  // correct. See ND-PHYSICS.md, A3.
  /**
   * The star matrix `[r]*` of the vector `r`. It maps a vector to a bivector,
   * and its transpose maps a bivector back to a vector:
   *   [r]* a = r ^ a          (a is a vector, the result is a bivector)
   *   [r]*^T w = r . w        (w is a bivector, the result is a vector)
   *
   * The row `p` holds the axis pair `D.pairs[p]`. Two examples of its use:
   * the velocity of a point is `v + [r]*^T w`, and an impulse `j` at the
   * offset `r` changes the angular momentum by `[r]* j`.
   *
   * @param {Dims} D the tables from `dims(n)`
   * @param {Float64Array} r the offset, of length `n`
   * @param {Float64Array} [out] a buffer of `k * n`
   * @returns {Float64Array} the matrix, `k` by `n`, row major
   */
  function starMatrix(D, r, out) {
    const n = D.n;
    const S = out || new Float64Array(D.k * n);
    S.fill(0);
    for (let p = 0; p < D.k; p += 1) {
      const [i, j] = D.pairs[p];
      S[p * n + j] = r[i];
      S[p * n + i] = -r[j];
    }
    return S;
  }
  /**
   * The bivector `a ^ b` of two vectors. This is the plane through them, and
   * its size is the area that they hold.
   *
   * It is the same as `[a]* b`, but it does not build the matrix. The solver
   * calls this for each contact, thus the speed is of value.
   *
   * @returns {Float64Array} the bivector, of length `k`
   */
  function wedgeVec(D, a, b, out) {
    const B = out || new Float64Array(D.k);
    for (let p = 0; p < D.k; p += 1) {
      const [i, j] = D.pairs[p];
      B[p] = a[i] * b[j] - a[j] * b[i];
    }
    return B;
  }
  /**
   * The vector `v . W` of a vector and a bivector. It is the same as
   * `[v]*^T W`, but it does not build the matrix.
   *
   * @param {Float64Array} v a vector of length `n`
   * @param {Float64Array} W a bivector of length `k`
   * @returns {Float64Array} the vector, of length `n`
   */
  function contractVecBi(D, v, W, out) {
    const n = D.n;
    const u = out || new Float64Array(n);
    u.fill(0);
    for (let p = 0; p < D.k; p += 1) {
      const [i, j] = D.pairs[p];
      const w = W[p];
      u[j] += v[i] * w;
      u[i] -= v[j] * w;
    }
    return u;
  }
  /**
   * The commutator `A x B = (A B - B A) / 2` of two bivectors. The result is
   * again a bivector. It uses the table `D.comm`.
   *
   * CAUTION: this is not a cross product. In 3 dimensions the two are the
   * same, but in 4 dimensions and more there is no cross product.
   *
   * The Euler equation `I dw/dt - w x I w = tau` uses this product. It is the
   * gyroscopic term. See ND-PHYSICS.md, A7.
   *
   * @returns {Float64Array} the bivector, of length `k`
   */
  function commutator(D, A, B, out) {
    const k = D.k;
    const C = out || new Float64Array(k);
    C.fill(0);
    const T = D.comm;
    for (let p = 0; p < k; p += 1) {
      const a = A[p];
      if (a === 0) continue;
      for (let q = 0; q < k; q += 1) {
        const b = B[q];
        if (b === 0) continue;
        const base = (p * k + q) * k;
        const ab = a * b;
        for (let s = 0; s < k; s += 1) {
          const t = T[base + s];
          if (t !== 0) C[s] += t * ab;
        }
      }
    }
    return C;
  }
  /**
   * The matrix of `k` by `k` such that `M B = X x B` for each bivector `B`.
   * Thus it is the commutator with `X`, in matrix form.
   *
   * `gyroscopicSpin` needs this form. Its Newton method must differentiate
   * the commutator term, and a matrix makes that possible.
   *
   * @returns {Float64Array} the matrix, `k` by `k`, row major
   */
  function commutatorMatrix(D, X, out) {
    const k = D.k;
    const M = out || new Float64Array(k * k);
    M.fill(0);
    const T = D.comm;
    for (let p = 0; p < k; p += 1) {
      const x = X[p];
      if (x === 0) continue;
      for (let q = 0; q < k; q += 1) {
        const base = (p * k + q) * k;
        for (let s = 0; s < k; s += 1) {
          const t = T[base + s];
          if (t !== 0) M[s * k + q] += t * x;
        }
      }
    }
    return M;
  }

  // src/nd/body/body.js
  /** The next automatic body id. The worker gives its own ids. */
  var nextId = 1;
  /**
   * A rigid body of `n` dimensions.
   *
   * The state, as ND-PHYSICS.md, A5, gives it:
   *   `x`  the position of the center of mass   vector    `n`
   *   `R`  the orientation                      rotor     `r`
   *   `v`  the linear velocity                  vector    `n`
   *   `L`  the angular momentum                 bivector  `k`
   *
   * The library integrates the momentum `L`, and not the angular velocity `w`.
   * `w` comes from `L` with `w = I'^-1 L` at each change. This keeps the
   * momentum correct when there is no torque.
   *
   * These fields come from `R`, and `updateDerived()` builds them:
   *   `Rm`               the rotation matrix, `n` by `n`
   *   `R2`               the bivector matrix `[R]2`, `k` by `k`
   *   `invInertiaWorld`  the inverse inertia in the world frame, `k` by `k`
   *   `w`                the angular velocity, bivector `k`
   *
   * The other fields:
   *   `inertia`         the inertia tensor in the body frame, `k` by `k`
   *   `invInertia`      its inverse
   *   `force`           the force of this step, vector `n`
   *   `torque`          the torque of this step, bivector `k`
   *   `linearFactor`    a multiplier of each component of the impulse, `n`
   *   `angularFactor`   a multiplier of each component of the torque, `k`
   *   `isStatic`        true for a body that no force can move
   *   `sleeping`        true for a body that the world does not integrate
   *   `level`           the distance to a static body. `shockPropagation`
   *                     sets it.
   *
   * A body with the mass 0, or with the shape `halfspace`, is static.
   */
  var Body = class {
    /**
     * @param {Dims} D the tables from `dims(n)`
     * @param {object} [opts] the options
     * @param {object} opts.shape a shape, for example `HyperBox(D, [1, 1, 1])`
     * @param {number} [opts.mass] the mass. 0 makes the body static. Default 1.
     * @param {number} [opts.id] the id. It gives an automatic id if you do not.
     * @param {string} [opts.name] a name, for your own use
     * @param {ArrayLike<number>} [opts.position] the position, `n`
     * @param {ArrayLike<number>} [opts.rotor] the orientation, `r`
     * @param {ArrayLike<number>} [opts.velocity] the linear velocity, `n`
     * @param {ArrayLike<number>} [opts.angularVelocity] the angular velocity, `k`
     * @param {number} [opts.friction] the friction. Default 0.5.
     * @param {number} [opts.restitution] the bounce, 0 to 1. Default 0.1.
     * @param {number} [opts.linearDamping] the linear damping. Default 0.
     * @param {number} [opts.angularDamping] the angular damping. Default 0.
     * @param {ArrayLike<number>} [opts.linearFactor] a multiplier of each axis, `n`
     * @param {ArrayLike<number>} [opts.angularFactor] a multiplier of each plane, `k`
     * @param {boolean} [opts.allowSleep] false keeps the body awake. Default true.
     */
    constructor(D, opts = {}) {
      const { n, k, r } = D;
      this.D = D;
      this.id = opts.id !== void 0 ? opts.id : nextId++;
      this.shape = opts.shape;
      this.name = opts.name || "";
      this.mass = opts.mass !== void 0 ? opts.mass : 1;
      this.isStatic = this.mass <= 0 || this.shape.type === "halfspace";
      if (this.isStatic) this.mass = 0;
      this.invMass = this.isStatic ? 0 : 1 / this.mass;
      this.x = new Float64Array(n);
      if (opts.position) this.x.set(opts.position);
      this.R = rotorIdentity(D);
      if (opts.rotor) this.R.set(opts.rotor);
      this.v = new Float64Array(n);
      if (opts.velocity) this.v.set(opts.velocity);
      this.L = new Float64Array(k);
      this.w = new Float64Array(k);
      this.inertia = this.isStatic ? matZero(k, k) : this.shape.inertia(this.mass);
      this.invInertia = this.isStatic ? matZero(k, k) : matInverseSPD(this.inertia, k);
      this.Rm = matZero(n, n);
      this.R2 = matZero(k, k);
      this.invInertiaWorld = matZero(k, k);
      this.force = new Float64Array(n);
      this.torque = new Float64Array(k);
      this.friction = opts.friction !== void 0 ? opts.friction : 0.5;
      this.restitution = opts.restitution !== void 0 ? opts.restitution : 0.1;
      this.linearDamping = opts.linearDamping || 0;
      this.angularDamping = opts.angularDamping || 0;
      this.linearFactor = new Float64Array(n).fill(1);
      if (opts.linearFactor) this.linearFactor.set(opts.linearFactor);
      this.angularFactor = new Float64Array(k).fill(1);
      if (opts.angularFactor) this.angularFactor.set(opts.angularFactor);
      this.sleeping = false;
      this.sleepTimer = 0;
      this.allowSleep = opts.allowSleep !== false;
      this.level = 0;
      this._island = 0;
      this._tmpN = new Float64Array(n);
      this._tmpK = new Float64Array(k);
      this._star = new Float64Array(k * n);
      this.updateDerived();
      if (opts.angularVelocity) this.setAngularVelocity(opts.angularVelocity);
    }
    /**
     * Builds the fields that come from the orientation `R`: the rotation
     * matrix `Rm`, the bivector matrix `R2`, the inverse inertia in the world
     * frame `I'^-1 = [R]2 I^-1 [R]2^T`, and the angular velocity `w = I'^-1 L`.
     *
     * Call this after each change of `R`. The integrator calls it one time in
     * each step. See ND-PHYSICS.md, A4 and A6.
     */
    updateDerived() {
      const { D } = this;
      const { n, k } = D;
      rotorMatrix(D, this.R, this.Rm);
      rotorBivectorMatrix(D, this.R, this.R2);
      if (!this.isStatic) {
        const T = matMul(this.R2, this.invInertia, k, k, k);
        matMulT(T, this.R2, k, k, k, this.invInertiaWorld);
        matVec(this.invInertiaWorld, this.L, k, k, this.w);
      } else {
        this.w.fill(0);
      }
    }
    /**
     * The inertia tensor in the world frame, `[R]2 I [R]2^T`, of `k` by `k`.
     * It makes a new matrix at each call. The step does not need it, because
     * `invInertiaWorld` holds the inverse. Use it for a test or for the energy.
     */
    inertiaWorld() {
      const { k } = this.D;
      const T = matMul(this.R2, this.inertia, k, k, k);
      return matMulT(T, this.R2, k, k, k);
    }
    /**
     * Sets the angular velocity, and builds the angular momentum `L = I' w`
     * from it. Give a bivector of `k` components, and not a vector.
     * @param {ArrayLike<number>} w the angular velocity, of length `k`
     */
    setAngularVelocity(w) {
      const { k } = this.D;
      this.w.set(w);
      if (this.isStatic) {
        this.L.fill(0);
        return;
      }
      const I = this.inertiaWorld();
      matVec(I, this.w, k, k, this.L);
    }
    /** Sets the linear velocity, of length `n`, and wakes the body. */
    setLinearVelocity(v) {
      this.v.set(v);
      this.wake();
    }
    /** Change a vector from the body frame to the world frame. */
    localToWorldDir(a, out) {
      return matVec(this.Rm, a, this.D.n, this.D.n, out);
    }
    /** Change a vector from the world frame to the body frame. */
    worldToLocalDir(a, out) {
      return matTVec(this.Rm, a, this.D.n, this.D.n, out);
    }
    /** Changes a point from the body frame to the world frame. It adds `x`. */
    localToWorld(a, out) {
      const p = this.localToWorldDir(a, out);
      for (let i = 0; i < this.D.n; i += 1) p[i] += this.x[i];
      return p;
    }
    /**
     * Changes a point from the world frame to the body frame. It takes `x`
     * away first. The result goes into the scratch buffer of the body when
     * you do not give `out`. Thus copy it before the next call.
     */
    worldToLocal(a, out) {
      const { n } = this.D;
      const t = this._tmpN;
      for (let i = 0; i < n; i += 1) t[i] = a[i] - this.x[i];
      return this.worldToLocalDir(t, out);
    }
    /**
     * Velocity of the point at the world offset `r` from the center of mass.
     *   u = v + r . w = v + [r]*^T w
     */
    pointVelocity(r, out) {
      const { n, k } = this.D;
      const u = out || new Float64Array(n);
      const S = starMatrix(this.D, r, this._star);
      matTVec(S, this.w, k, n, u);
      for (let i = 0; i < n; i += 1) u[i] += this.v[i];
      return u;
    }
    /**
     * Applies an impulse `j` at the world offset `r` from the center of mass.
     * The linear part is `v += j / m`. The angular part is `L += [r]* j`, and
     * then `w` comes again from `L`. See ND-PHYSICS.md, A11.
     *
     * A static body does not change. The function wakes the body.
     *
     * @param {Float64Array} j the impulse, a vector of length `n`
     * @param {Float64Array} r the offset, a vector of length `n`
     * @param {boolean} [fromSolver] true when the solver calls this. The body
     *   then keeps its sleep timer. See `touch`. A caller from outside must
     *   not give this, thus an impulse of a program always wakes the body.
     */
    applyImpulse(j, r, fromSolver) {
      if (this.isStatic) return;
      const { n, k } = this.D;
      for (let i = 0; i < n; i += 1) this.v[i] += this.invMass * j[i] * this.linearFactor[i];
      const S = starMatrix(this.D, r, this._star);
      const dL = matVec(S, j, k, n, this._tmpK);
      for (let p = 0; p < k; p += 1) this.L[p] += dL[p] * this.angularFactor[p];
      matVec(this.invInertiaWorld, this.L, k, k, this.w);
      if (fromSolver) this.touch();
      else this.wake();
    }
    /** Applies an impulse at the center of mass. The body does not start to turn. */
    applyCentralImpulse(j) {
      if (this.isStatic) return;
      for (let i = 0; i < this.D.n; i += 1) this.v[i] += this.invMass * j[i] * this.linearFactor[i];
      this.wake();
    }
    /**
     * Adds `dL` to the angular momentum. Give a bivector of length `k`.
     * @param {boolean} [fromSolver] see `applyImpulse`
     */
    applyTorqueImpulse(dL, fromSolver) {
      if (this.isStatic) return;
      const { k } = this.D;
      for (let p = 0; p < k; p += 1) this.L[p] += dL[p] * this.angularFactor[p];
      matVec(this.invInertiaWorld, this.L, k, k, this.w);
      if (fromSolver) this.touch();
      else this.wake();
    }
    /**
     * Adds a force at the center of mass. The force holds until the end of the
     * step, and then `World.step` makes it zero again.
     */
    applyCentralForce(f) {
      for (let i = 0; i < this.D.n; i += 1) this.force[i] += f[i];
      this.wake();
    }
    /**
     * Adds a force `f` at the world offset `r`. It also adds the torque
     * `[r]* f`. Both hold until the end of the step.
     */
    applyForce(f, r) {
      const { n, k } = this.D;
      for (let i = 0; i < n; i += 1) this.force[i] += f[i];
      const S = starMatrix(this.D, r, this._star);
      const t = matVec(S, f, k, n, this._tmpK);
      for (let p = 0; p < k; p += 1) this.torque[p] += t[p];
      this.wake();
    }
    /** Adds a torque. Give a bivector of length `k`, and not a vector. */
    applyTorque(t) {
      for (let p = 0; p < this.D.k; p += 1) this.torque[p] += t[p];
      this.wake();
    }
    /** Makes the force and the torque zero. `World.step` calls this at the end. */
    clearForces() {
      this.force.fill(0);
      this.torque.fill(0);
    }
    /** Wakes the body and makes its sleep timer zero. */
    wake() {
      if (this.sleeping) {
        this.sleeping = false;
      }
      this.sleepTimer = 0;
    }
    /**
     * Wakes a body that sleeps, and does NOT touch the timer of a body that is
     * already awake.
     *
     * The solver uses this, and not `wake()`. An impulse of a contact or of a
     * row is internal: it holds the body where it is, and it must not stop the
     * body from going to sleep. `wake()` makes `sleepTimer` zero at every call,
     * thus a box that rests on the ground could never sleep, because the
     * contact gives it an impulse in every turn of the solver.
     */
    touch() {
      if (this.sleeping) {
        this.sleeping = false;
        this.sleepTimer = 0;
      }
    }
    /**
     * The kinetic energy `(m v.v + w.L) / 2`. The angular part uses `w` and
     * `L` together, thus it does not need the inertia tensor.
     *
     * A body with no torque must keep this value. That is the test of the
     * gyroscopic term. See ND-PHYSICS.md, A7.
     */
    kineticEnergy() {
      const { n, k } = this.D;
      if (this.isStatic) return 0;
      let e = 0;
      for (let i = 0; i < n; i += 1) e += this.v[i] * this.v[i];
      e *= 0.5 * this.mass;
      let a = 0;
      for (let p = 0; p < k; p += 1) a += this.w[p] * this.L[p];
      return e + 0.5 * a;
    }
    /**
     * The axis aligned box of the body in the world frame. The broad phase
     * uses it. There are three conditions:
     *   - A half space fills all of the space, thus the box has no limit.
     *   - A box turns with the body. The extent on the axis `i` is the sum of
     *     `|Rm[i][j]| * halfExtents[j]` over all `j`.
     *   - All other shapes use the bounding radius. The box is then larger
     *     than the shape, but it is correct.
     *
     * @param {number} [margin] a length to add on each side. Default 0.
     * @returns {{min: Float64Array, max: Float64Array}} the two corners, each `n`
     */
    aabb(margin = 0) {
      const { n } = this.D;
      const min = new Float64Array(n);
      const max = new Float64Array(n);
      if (this.shape.type === "halfspace") {
        min.fill(-Infinity);
        max.fill(Infinity);
        return { min, max };
      }
      if (this.shape.type === "box") {
        const h = this.shape.halfExtents;
        for (let i = 0; i < n; i += 1) {
          let e = 0;
          for (let j = 0; j < n; j += 1) e += Math.abs(this.Rm[i * n + j]) * h[j];
          min[i] = this.x[i] - e - margin;
          max[i] = this.x[i] + e + margin;
        }
        return { min, max };
      }
      const rad = this.shape.boundingRadius + margin;
      for (let i = 0; i < n; i += 1) {
        min[i] = this.x[i] - rad;
        max[i] = this.x[i] + rad;
      }
      return { min, max };
    }
  };

  // src/nd/body/massprops.js
  //
  // The mass, the center of mass and the inertia tensor of a simplicial mesh.
  //
  // The method is the covariance method of ND-PHYSICS.md, A8. Do not integrate
  // the products of inertia one by one. Integrate the covariance matrix `C` of
  // `n` by `n` instead, then change it into the inertia tensor of `k` by `k`.
  // The covariance method is shorter, and it works for each `n`.
  /** The `starProducts` tables of each `D`. */
  var productCache = /* @__PURE__ */ new WeakMap();
  /**
   * The `n` by `n` grid of matrices `P[i][j] = eStar[i] eStar[j]^T`, each of
   * `k` by `k`. `inertiaFromCovariance` adds these together.
   *
   * The result is in a cache, because it costs much and it does not change.
   */
  function starProducts(D) {
    let P = productCache.get(D);
    if (P) return P;
    const { n, k, eStar } = D;
    P = [];
    for (let i = 0; i < n; i += 1) {
      P.push([]);
      for (let j = 0; j < n; j += 1) {
        const M = matZero(k, k);
        for (let p = 0; p < k; p += 1) {
          for (let q = 0; q < k; q += 1) {
            let s = 0;
            for (let t = 0; t < n; t += 1) s += eStar[i][p * n + t] * eStar[j][q * n + t];
            M[p * k + q] = s;
          }
        }
        P[i].push(M);
      }
    }
    productCache.set(D, P);
    return P;
  }
  /**
   * The inertia tensor of `k` by `k`, from the covariance matrix `C` of `n` by
   * `n`. The covariance holds the integral of `rho x_i x_j` over the body.
   *
   * The inertia tensor is `I = INT rho [r]* [r]*^T dV`. The star matrix is
   * linear in `r`. Thus the integral becomes a sum of `C[i][j]` times the
   * constant matrix `eStar[i] eStar[j]^T`. See ND-PHYSICS.md, A6 and A8.
   *
   * @param {Float64Array} C the covariance, `n` by `n`, about the center of mass
   * @returns {Float64Array} the inertia tensor, `k` by `k`
   */
  function inertiaFromCovariance(D, C, out) {
    const { n, k } = D;
    const P = starProducts(D);
    const I = out || matZero(k, k);
    I.fill(0);
    for (let i = 0; i < n; i += 1) {
      for (let j = 0; j < n; j += 1) {
        const c = C[i * n + j];
        if (c === 0) continue;
        const M = P[i][j];
        for (let t = 0; t < k * k; t += 1) I[t] += c * M[t];
      }
    }
    return I;
  }
  /**
   * The volume, the mass, the center of mass and the inertia tensor of a
   * closed simplicial mesh.
   *
   * The surface of the mesh is a set of `(n-1)`-simplices: triangles in 3
   * dimensions, tetrahedra in 4. Each surface simplex, together with the
   * origin, makes an `n`-simplex. The integrals are linear, thus the function
   * adds the part of each `n`-simplex.
   *
   * For one simplex, with the matrix `M` of its `n` vertices as the columns:
   *   volume    = det(M) / n!
   *   moment    = det(M) * (sum of the columns) / (n+1)!
   *   covariance = det(M) * M C0 M^T,  with C0 the canonical value of A8
   * The determinant is signed. Thus a part outside the body is negative, and
   * the sum over a closed surface gives the true body. The direction of the
   * surface simplices must be the same everywhere. `orientCells` does that.
   *
   * The function turns the sign around when the total volume is negative. At
   * the end it moves the covariance to the center of mass.
   *
   * See ND-PHYSICS.md, A8.
   *
   * @param {Dims} D the tables from `dims(n)`
   * @param {Float64Array} vertices the vertices, `n` numbers for each
   * @param {ArrayLike<number>} cells the surface simplices, `n` indices for each
   * @param {number} [density] the density. Default 1.
   * @returns {{volume: number, mass: number, center: Float64Array,
   *   covariance: Float64Array, inertia: Float64Array}} the properties. The
   *   center has `n`, the covariance is `n` by `n`, the inertia `k` by `k`.
   *   Both matrices are about the center of mass.
   * @throws {Error} when the volume is zero
   */
  function massProperties(D, vertices, cells, density = 1) {
    const { n, cXX, cXY, simplexVolumeDiv, simplexMomentDiv } = D;
    const cellCount = cells.length / n;
    const M = matZero(n, n);
    const C0 = matZero(n, n);
    for (let i = 0; i < n; i += 1) {
      for (let j = 0; j < n; j += 1) C0[i * n + j] = i === j ? cXX : cXY;
    }
    let volume = 0;
    const moment = new Float64Array(n);
    const cov = matZero(n, n);
    const MC = matZero(n, n);
    for (let c = 0; c < cellCount; c += 1) {
      for (let col = 0; col < n; col += 1) {
        const vi = cells[c * n + col] * n;
        for (let row = 0; row < n; row += 1) M[row * n + col] = vertices[vi + row];
      }
      const det2 = matDet(M, n);
      if (det2 === 0) continue;
      volume += det2 / simplexVolumeDiv;
      for (let i = 0; i < n; i += 1) {
        let s = 0;
        for (let j = 0; j < n; j += 1) s += M[i * n + j];
        moment[i] += det2 * s / simplexMomentDiv;
      }
      MC.fill(0);
      for (let i = 0; i < n; i += 1) {
        for (let t = 0; t < n; t += 1) {
          const a = M[i * n + t];
          if (a === 0) continue;
          for (let j = 0; j < n; j += 1) MC[i * n + j] += a * C0[t * n + j];
        }
      }
      for (let i = 0; i < n; i += 1) {
        for (let j = 0; j < n; j += 1) {
          let s = 0;
          for (let t = 0; t < n; t += 1) s += MC[i * n + t] * M[j * n + t];
          cov[i * n + j] += det2 * s;
        }
      }
    }
    if (volume < 0) {
      volume = -volume;
      for (let i = 0; i < n; i += 1) moment[i] = -moment[i];
      for (let i = 0; i < n * n; i += 1) cov[i] = -cov[i];
    }
    if (volume <= 0) throw new Error("massProperties: the volume is zero");
    const center = new Float64Array(n);
    for (let i = 0; i < n; i += 1) center[i] = moment[i] / volume;
    const covC = matZero(n, n);
    for (let i = 0; i < n; i += 1) {
      for (let j = 0; j < n; j += 1) {
        covC[i * n + j] = density * (cov[i * n + j] - volume * center[i] * center[j]);
      }
    }
    return {
      volume,
      mass: density * volume,
      center,
      covariance: covC,
      inertia: inertiaFromCovariance(D, covC)
    };
  }

  // src/nd/body/shapes.js
  //
  // The five shapes. Each one is a plain object, and not a class. All of them
  // obey the same interface, thus the collision code can use any of them.
  /**
   * The interface of a shape. A shape does not know its body. All of its
   * values are in the body frame, and the center of mass is at the origin.
   *
   * @typedef {object} Shape
   * @property {string} type one of `box`, `sphere`, `halfspace`, `convex`,
   *   `torus`. The dispatch of `collide()` uses this value.
   * @property {number} n the count of dimensions
   * @property {number} boundingRadius the radius of a ball, at the origin,
   *   that holds all of the shape. It is `Infinity` for a half space.
   * @property {number} volume the `n`-volume
   * @property {function(Float64Array, Float64Array=): Float64Array} support
   *   the point of the shape that is farthest along the direction `dir`. A
   *   half space throws, because it has no such point.
   * @property {function(): ?{vertices: Float64Array, cells: Int32Array}} mesh
   *   the surface as simplices, or null when the shape has no mesh. The 4D
   *   slice needs a mesh.
   * @property {function(number): Float64Array} inertia the inertia tensor of
   *   the shape with the given mass, `k` by `k`, about the center of mass.
   */
  /** All of the orders of the members of `list`. `hyperBoxMesh` uses it. */
  function permutations(list) {
    if (list.length <= 1) return [list.slice()];
    const out = [];
    for (let i = 0; i < list.length; i += 1) {
      const rest = list.slice(0, i).concat(list.slice(i + 1));
      for (const p of permutations(rest)) out.push([list[i]].concat(p));
    }
    return out;
  }
  /**
   * The surface of a box as simplices. In 4 dimensions this gives tetrahedra,
   * and the slice needs them.
   *
   * The method is the Kuhn cut. Start at one corner of a face, then move along
   * the axes of that face, one axis at a time. Each order of the axes gives
   * one simplex. Thus each of the `2n` faces gives `(n-1)!` simplices.
   *
   * At the end the function turns each simplex that has a negative
   * determinant around, thus all of them have the same direction. The mass
   * properties need that.
   *
   * @param {Float64Array} h the half extents, of length `n`
   * @returns {{vertices: Float64Array, cells: Int32Array}} the `2^n` corners,
   *   and the simplices as `n` indices for each
   */
  function hyperBoxMesh(D, h) {
    const n = D.n;
    const corners = 1 << n;
    const vertices = new Float64Array(corners * n);
    for (let m = 0; m < corners; m += 1) {
      for (let i = 0; i < n; i += 1) vertices[m * n + i] = m >> i & 1 ? h[i] : -h[i];
    }
    const cells = [];
    for (let d = 0; d < n; d += 1) {
      const others = [];
      for (let i = 0; i < n; i += 1) if (i !== d) others.push(i);
      const perms = permutations(others);
      for (const s of [0, 1]) {
        for (const perm of perms) {
          let mask = s ? 1 << d : 0;
          const idx = [mask];
          for (const ax of perm) {
            mask |= 1 << ax;
            idx.push(mask);
          }
          cells.push(idx);
        }
      }
    }
    const M = matZero(n, n);
    const flat = new Int32Array(cells.length * n);
    for (let c = 0; c < cells.length; c += 1) {
      const idx = cells[c];
      for (let col = 0; col < n; col += 1) {
        for (let row = 0; row < n; row += 1) M[row * n + col] = vertices[idx[col] * n + row];
      }
      if (matDet(M, n) < 0) {
        const t = idx[0];
        idx[0] = idx[1];
        idx[1] = t;
      }
      for (let col = 0; col < n; col += 1) flat[c * n + col] = idx[col];
    }
    return { vertices, cells: flat };
  }
  /**
   * The inertia tensor of a box, `k` by `k`. It is diagonal: the axes of the
   * box are its principal axes. The element of the plane of the axes `i` and
   * `j` is `m (h_i^2 + h_j^2) / 3`.
   */
  function hyperBoxInertia(D, h, mass) {
    const { n, k } = D;
    const I = matZero(k, k);
    for (let p = 0; p < k; p += 1) {
      const [i, j] = D.pairs[p];
      I[p * k + p] = mass / 3 * (h[i] * h[i] + h[j] * h[j]);
    }
    return I;
  }
  /**
   * The inertia tensor of a solid ball, `k` by `k`. Each diagonal element is
   * `2 m R^2 / (n + 2)`. A ball is the same in each direction, thus the tensor
   * is a multiple of the identity. In 3 dimensions this gives `2 m R^2 / 5`.
   */
  function hyperSphereInertia(D, radius, mass) {
    const { n, k } = D;
    const I = matZero(k, k);
    const v = 2 * mass * radius * radius / (n + 2);
    for (let p = 0; p < k; p += 1) I[p * k + p] = v;
    return I;
  }
  /**
   * The volume of a ball of `n` dimensions. The formula has two conditions,
   * because the gamma function of a half integer is different:
   *   n even:  pi^(n/2) R^n / (n/2)!
   *   n odd:   2^((n+1)/2) pi^((n-1)/2) R^n / (1 * 3 * 5 * ... * n)
   * In 3 dimensions this gives `4 pi R^3 / 3`, and in 4 `pi^2 R^4 / 2`.
   */
  function ballVolume(n, radius) {
    let v;
    if (n % 2 === 0) {
      const m = n / 2;
      let f = 1;
      for (let i = 2; i <= m; i += 1) f *= i;
      v = Math.PI ** m / f;
    } else {
      const m = (n - 1) / 2;
      let num = 2 ** (m + 1);
      let den = 1;
      for (let i = 1; i <= n; i += 2) den *= i;
      v = num * Math.PI ** m / den;
    }
    return v * radius ** n;
  }
  /**
   * A box of `n` dimensions, at the origin, in line with the axes of the body.
   * In 4 dimensions this is a tesseract.
   * @param {ArrayLike<number>} halfExtents the half length on each axis, `n`
   * @returns {Shape} the shape, of the type `box`
   * @throws {Error} when the count of half extents is not `n`
   */
  function HyperBox(D, halfExtents) {
    const h = Float64Array.from(halfExtents);
    if (h.length !== D.n) throw new Error("HyperBox: wrong number of half extents");
    let bound = 0;
    for (let i = 0; i < D.n; i += 1) bound += h[i] * h[i];
    return {
      type: "box",
      n: D.n,
      halfExtents: h,
      boundingRadius: Math.sqrt(bound),
      volume: h.reduce((a, b) => a * 2 * b, 1),
      support(dir, out) {
        const p = out || new Float64Array(D.n);
        for (let i = 0; i < D.n; i += 1) p[i] = dir[i] >= 0 ? h[i] : -h[i];
        return p;
      },
      mesh() {
        return hyperBoxMesh(D, h);
      },
      inertia(mass) {
        return hyperBoxInertia(D, h, mass);
      }
    };
  }
  /**
   * A solid ball of `n` dimensions, at the origin. Its `mesh()` gives null,
   * because the slice of a ball is again a ball. `sliceHyperSphereRadius`
   * gives the radius of that slice directly.
   * @returns {Shape} the shape, of the type `sphere`
   */
  function HyperSphere(D, radius) {
    return {
      type: "sphere",
      n: D.n,
      radius,
      boundingRadius: radius,
      volume: ballVolume(D.n, radius),
      support(dir, out) {
        const p = out || new Float64Array(D.n);
        let ln = 0;
        for (let i = 0; i < D.n; i += 1) ln += dir[i] * dir[i];
        ln = Math.sqrt(ln) || 1;
        for (let i = 0; i < D.n; i += 1) p[i] = dir[i] / ln * radius;
        return p;
      },
      mesh() {
        return null;
      },
      inertia(mass) {
        return hyperSphereInertia(D, radius, mass);
      }
    };
  }
  /**
   * The half space `normal . x <= offset`. This is the ground, or a wall.
   *
   * A body with this shape is always static, and it has no mass and no
   * inertia. The normal points away from the solid part. `support()` throws,
   * because a half space has no farthest point.
   *
   * The function makes the normal a unit vector.
   *
   * @param {ArrayLike<number>} normal the normal, of length `n`
   * @param {number} [offset] the distance from the origin. Default 0.
   * @returns {Shape} the shape, of the type `halfspace`
   */
  function HalfSpace(D, normal, offset = 0) {
    const nv = Float64Array.from(normal);
    let ln = 0;
    for (let i = 0; i < D.n; i += 1) ln += nv[i] * nv[i];
    ln = Math.sqrt(ln);
    for (let i = 0; i < D.n; i += 1) nv[i] /= ln;
    return {
      type: "halfspace",
      n: D.n,
      normal: nv,
      offset,
      boundingRadius: Infinity,
      volume: Infinity,
      support() {
        throw new Error("HalfSpace: a half space has no support point");
      },
      mesh() {
        return null;
      },
      inertia() {
        return matZero(D.k, D.k);
      }
    };
  }
  /**
   * A convex body from a closed surface mesh. The surface is a set of
   * `(n-1)`-simplices: triangles in 3 dimensions, tetrahedra in 4.
   *
   * The origin of the mesh must be inside the body. `orientCells` gives all of
   * the simplices the same direction about the origin, thus the direction of
   * the input does not matter. `massProperties` then gives the volume and the
   * inertia tensor.
   *
   * `support()` looks at each vertex in turn. Thus a mesh of many vertices is
   * slow. A collision with a torus or with another convex mesh has no test.
   * See `collide()`.
   *
   * @param {ArrayLike<number>} vertices the vertices, `n` numbers for each
   * @param {ArrayLike<number>} cells the surface simplices, `n` indices for each
   * @returns {Shape} the shape, of the type `convex`
   */
  function ConvexMesh(D, vertices, cells) {
    const V = Float64Array.from(vertices);
    const grouped = [];
    for (let c = 0; c < cells.length; c += D.n) {
      const idx = [];
      for (let i = 0; i < D.n; i += 1) idx.push(cells[c + i]);
      grouped.push(idx);
    }
    const ORIGIN = new Float64Array(D.n);
    const C = orientCells(V, grouped, D.n, () => ORIGIN);
    const count = V.length / D.n;
    let bound = 0;
    for (let m = 0; m < count; m += 1) {
      let s = 0;
      for (let i = 0; i < D.n; i += 1) s += V[m * D.n + i] * V[m * D.n + i];
      if (s > bound) bound = s;
    }
    const props = massProperties(D, V, C, 1);
    return {
      type: "convex",
      n: D.n,
      vertices: V,
      cells: C,
      boundingRadius: Math.sqrt(bound),
      volume: props.volume,
      support(dir, out) {
        const p = out || new Float64Array(D.n);
        let best = -Infinity;
        let bi = 0;
        for (let m = 0; m < count; m += 1) {
          let s = 0;
          for (let i = 0; i < D.n; i += 1) s += V[m * D.n + i] * dir[i];
          if (s > best) {
            best = s;
            bi = m;
          }
        }
        for (let i = 0; i < D.n; i += 1) p[i] = V[bi * D.n + i];
        return p;
      },
      mesh() {
        return { vertices: V, cells: C };
      },
      inertia(mass) {
        const p = massProperties(D, V, C, mass / props.volume);
        return p.inertia;
      }
    };
  }
  /**
   * The volume of a torus of `n` dimensions. Pappus gives it: the length of
   * the major circle, `2 pi R`, times the volume of the `(n-1)`-ball of the
   * minor radius `r`.
   */
  function torusVolume(n, R, r) {
    return 2 * Math.PI * R * ballVolume(n - 1, r);
  }
  /**
   * The covariance matrix of a torus, `n` by `n`. It is diagonal. The two axes
   * of the major plane hold `m (R^2 + 3 r^2 / (n+1)) / 2`, and each other axis
   * holds `m r^2 / (n+1)`.
   * @param {number[]} plane the two axes of the major circle
   */
  function torusCovariance(D, R, r, plane, mass) {
    const { n } = D;
    const C = matZero(n, n);
    const major = mass * (R * R + 3 * r * r / (n + 1)) / 2;
    const minor = mass * r * r / (n + 1);
    for (let i = 0; i < n; i += 1) {
      C[i * n + i] = i === plane[0] || i === plane[1] ? major : minor;
    }
    return C;
  }
  /** The inertia tensor of a torus, `k` by `k`, through its covariance. */
  function torusInertia(D, R, r, plane, mass) {
    return inertiaFromCovariance(D, torusCovariance(D, R, r, plane, mass));
  }
  /**
   * The surface of a torus as simplices. Only 3 and 4 dimensions have a mesh.
   * @param {number[]} [segments] the counts of the divisions. The default is
   *   `[32, 16]` in 3 dimensions and `[16, 8, 12]` in 4.
   * @throws {Error} when `n` is not 3 and not 4
   */
  function torusMesh(D, R, r, plane, segments) {
    const { n } = D;
    if (n === 3) return torusMesh3(D, R, r, plane, segments || [32, 16]);
    if (n === 4) return torusMesh4(D, R, r, plane, segments || [16, 8, 12]);
    throw new Error("torusMesh: only n = 3 and n = 4 have a mesh");
  }
  /** The axes that are not in `plane`. A torus in 4 dimensions has two. */
  function otherAxes(n, plane) {
    const out = [];
    for (let i = 0; i < n; i += 1) if (i !== plane[0] && i !== plane[1]) out.push(i);
    return out;
  }
  /**
   * Gives all of the surface simplices the same direction, then makes them
   * flat. `massProperties` needs one direction, because it uses the sign of
   * the determinant.
   *
   * For each simplex the function takes a point inside the body, near the
   * simplex. It then makes the matrix of the vertices about that point. A
   * negative determinant shows the wrong direction, and the function changes
   * the first two vertices for each other.
   *
   * @param {function(Float64Array): Float64Array} insidePoint gives a point
   *   inside the body, from the center of the simplex. A convex mesh uses the
   *   origin. A torus uses `insideTorus`, because the origin of a torus is
   *   not inside it.
   * @returns {Int32Array} the simplices, flat, `n` indices for each
   */
  function orientCells(vertices, cells, n, insidePoint) {
    const M = matZero(n, n);
    const g = new Float64Array(n);
    for (let c = 0; c < cells.length; c += 1) {
      const idx = cells[c];
      g.fill(0);
      for (let col = 0; col < n; col += 1) {
        for (let row = 0; row < n; row += 1) g[row] += vertices[idx[col] * n + row] / n;
      }
      const o = insidePoint(g);
      for (let col = 0; col < n; col += 1) {
        for (let row = 0; row < n; row += 1) {
          M[row * n + col] = vertices[idx[col] * n + row] - o[row];
        }
      }
      if (matDet(M, n) < 0) {
        const t = idx[0];
        idx[0] = idx[1];
        idx[1] = t;
      }
    }
    const flat = new Int32Array(cells.length * n);
    for (let c = 0; c < cells.length; c += 1) {
      for (let col = 0; col < n; col += 1) flat[c * n + col] = cells[c][col];
    }
    return flat;
  }
  /**
   * The surface of a 3D torus as triangles. It is a grid of `nt` by `nb`: `nt`
   * steps along the major circle, and `nb` around the tube. Each square of the
   * grid gives two triangles.
   */
  function torusMesh3(D, R, r, plane, seg) {
    const n = 3;
    const [nt, nb] = seg;
    const [i, j] = plane;
    const m = otherAxes(n, plane)[0];
    const vertices = new Float64Array(nt * nb * n);
    const at = (t, b) => t % nt * nb + b % nb;
    for (let t = 0; t < nt; t += 1) {
      const th = 2 * Math.PI * t / nt;
      for (let b = 0; b < nb; b += 1) {
        const be = 2 * Math.PI * b / nb;
        const rad = R + r * Math.cos(be);
        const base = at(t, b) * n;
        vertices[base + i] = rad * Math.cos(th);
        vertices[base + j] = rad * Math.sin(th);
        vertices[base + m] = r * Math.sin(be);
      }
    }
    const cells = [];
    for (let t = 0; t < nt; t += 1) {
      for (let b = 0; b < nb; b += 1) {
        const v00 = at(t, b);
        const v10 = at(t + 1, b);
        const v01 = at(t, b + 1);
        const v11 = at(t + 1, b + 1);
        cells.push([v00, v10, v11]);
        cells.push([v00, v11, v01]);
      }
    }
    return { vertices, cells: orientCells(vertices, cells, n, insideTorus(n, R, plane)) };
  }
  /**
   * Gives a function that finds a point inside a torus, for `orientCells`.
   * The origin of a torus is in the hole, thus it is not inside the body. The
   * function moves to the nearest point of the major circle, then it goes half
   * of the way back to the given point.
   */
  function insideTorus(n, R, plane) {
    const o = new Float64Array(n);
    return (g) => {
      o.fill(0);
      const mag = Math.hypot(g[plane[0]], g[plane[1]]);
      let ci = 1;
      let cj = 0;
      if (mag > 1e-12) {
        ci = g[plane[0]] / mag;
        cj = g[plane[1]] / mag;
      }
      o[plane[0]] = R * ci;
      o[plane[1]] = R * cj;
      for (let a = 0; a < n; a += 1) o[a] = o[a] + (g[a] - o[a]) * 0.5;
      return o;
    };
  }
  /** The six orders of three axes. They cut a cube into six tetrahedra. */
  var KUHN3 = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0]
  ];
  /**
   * The surface of a 4D torus as tetrahedra. The surface is a grid of three
   * angles: `nt` steps along the major circle, and a 2-sphere of `na` by `nb`
   * around the tube. Each box of the grid gives six tetrahedra, with the Kuhn
   * cut of `KUHN3`.
   *
   * The 4D slice needs these tetrahedra. See `sliceTetrahedra`.
   */
  function torusMesh4(D, R, r, plane, seg) {
    const n = 4;
    const [nt, na, nb] = seg;
    const [i, j] = plane;
    const others = otherAxes(n, plane);
    const p = others[0];
    const q = others[1];
    const rows = na + 1;
    const vertices = new Float64Array(nt * rows * nb * n);
    const at = (t, a, b) => (t % nt * rows + a) * nb + b % nb;
    for (let t = 0; t < nt; t += 1) {
      const th = 2 * Math.PI * t / nt;
      for (let a = 0; a < rows; a += 1) {
        const al = Math.PI * a / na;
        for (let b = 0; b < nb; b += 1) {
          const be = 2 * Math.PI * b / nb;
          const da = Math.sin(al) * Math.cos(be);
          const db = Math.sin(al) * Math.sin(be);
          const dc = Math.cos(al);
          const rad = R + r * da;
          const base = at(t, a, b) * n;
          vertices[base + i] = rad * Math.cos(th);
          vertices[base + j] = rad * Math.sin(th);
          vertices[base + p] = r * db;
          vertices[base + q] = r * dc;
        }
      }
    }
    const cells = [];
    for (let t = 0; t < nt; t += 1) {
      for (let a = 0; a < na; a += 1) {
        for (let b = 0; b < nb; b += 1) {
          const corner = (dt, da, db) => at(t + dt, a + da, b + db);
          for (const perm of KUHN3) {
            const step = [0, 0, 0];
            const idx = [corner(0, 0, 0)];
            for (const ax of perm) {
              step[ax] = 1;
              idx.push(corner(step[0], step[1], step[2]));
            }
            cells.push(idx);
          }
        }
      }
    }
    return { vertices, cells: orientCells(vertices, cells, n, insideTorus(n, R, plane)) };
  }
  /**
   * A torus of `n` dimensions. All of the points at the distance
   * `minorRadius` from a circle of the radius `majorRadius`. That circle is
   * the major circle, and it lies in the plane of two axes.
   *
   * In 4 dimensions the tube is a 2-sphere, and not a circle.
   *
   * The shape has two functions of its own, `corePoint` and `coreDistance`.
   * The collision with a hypersphere uses them.
   *
   * @param {number} majorRadius the radius of the major circle
   * @param {number} minorRadius the radius of the tube
   * @param {number[]} [plane] the two axes of the major plane. Default `[0, 1]`.
   * @param {number[]} [meshSegments] the divisions of the mesh
   * @returns {Shape} the shape, of the type `torus`
   * @throws {Error} when the two axes of the plane are not different and valid
   */
  function Torus(D, majorRadius, minorRadius, plane = [0, 1], meshSegments = null) {
    const { n } = D;
    const R = majorRadius;
    const r = minorRadius;
    const pl = [plane[0], plane[1]];
    if (pl[0] === pl[1] || pl[0] < 0 || pl[1] < 0 || pl[0] >= n || pl[1] >= n) {
      throw new Error("Torus: the major plane needs two different axes");
    }
    const others = otherAxes(n, pl);
    return {
      type: "torus",
      n,
      majorRadius: R,
      minorRadius: r,
      plane: pl,
      otherAxes: others,
      boundingRadius: R + r,
      volume: torusVolume(n, R, r),
      /**
       * The support point. The direction of the largest extent is `R` times the
       * unit vector of the part of `dir` in the major plane, plus `r` times the
       * unit vector of `dir`.
       */
      support(dir, out) {
        const s = out || new Float64Array(n);
        let ln = 0;
        for (let a = 0; a < n; a += 1) ln += dir[a] * dir[a];
        ln = Math.sqrt(ln) || 1;
        const mi = dir[pl[0]] / ln;
        const mj = dir[pl[1]] / ln;
        let mag = Math.hypot(mi, mj);
        let ci = 1;
        let cj = 0;
        if (mag > 1e-12) {
          ci = mi / mag;
          cj = mj / mag;
        }
        for (let a = 0; a < n; a += 1) s[a] = dir[a] / ln * r;
        s[pl[0]] += R * ci;
        s[pl[1]] += R * cj;
        return s;
      },
      /**
       * The point on the major circle that is nearest to `p`, and the distance
       * from `p` to that point. `p` must be in the frame of the body.
       */
      corePoint(p, out) {
        const c = out || new Float64Array(n);
        c.fill(0);
        const mag = Math.hypot(p[pl[0]], p[pl[1]]);
        if (mag > 1e-12) {
          c[pl[0]] = p[pl[0]] / mag * R;
          c[pl[1]] = p[pl[1]] / mag * R;
        } else {
          c[pl[0]] = R;
        }
        return c;
      },
      /** The distance from a point in the body frame to the major circle. */
      coreDistance(p) {
        const mag = Math.hypot(p[pl[0]], p[pl[1]]);
        let s = (mag - R) * (mag - R);
        for (const a of others) s += p[a] * p[a];
        return Math.sqrt(s);
      },
      mesh(segments) {
        return torusMesh(D, R, r, pl, segments || meshSegments);
      },
      inertia(mass) {
        return torusInertia(D, R, r, pl, mass);
      }
    };
  }

  // src/nd/detect/nearest.js
  //
  // The nearest point of a mesh to a given point. The collision of a
  // hypersphere with a convex mesh uses it. See ND-PHYSICS.md, A9.
  /**
   * The point of a simplex that is nearest to `q`.
   *
   * The method: put the point on the plane of the simplex, with the Gram
   * system of the edge vectors. When all of the barycentric weights are
   * between 0 and 1, that point is inside the simplex, and it is the result.
   *
   * If it is not inside, the nearest point is on a face. The function then
   * drops each vertex in turn, and it calls itself on the smaller simplex.
   * The nearest of those results is the answer. The recursion stops at one
   * vertex.
   *
   * @param {ArrayLike<number>[]} pts the vertices, each of length `n`
   * @param {Float64Array} q the point, of length `n`
   * @param {number} n the count of dimensions
   * @returns {Float64Array} the nearest point, of length `n`
   */
  function nearestOnSimplex(pts, q, n) {
    const m = pts.length - 1;
    if (m === 0) return pts[0];
    const E = [];
    for (let i = 1; i <= m; i += 1) {
      const e = new Float64Array(n);
      for (let a = 0; a < n; a += 1) e[a] = pts[i][a] - pts[0][a];
      E.push(e);
    }
    const G = new Float64Array(m * m);
    const b = new Float64Array(m);
    for (let i = 0; i < m; i += 1) {
      for (let j = 0; j < m; j += 1) {
        let s2 = 0;
        for (let a = 0; a < n; a += 1) s2 += E[i][a] * E[j][a];
        G[i * m + j] = s2;
      }
      let s = 0;
      for (let a = 0; a < n; a += 1) s += (q[a] - pts[0][a]) * E[i][a];
      b[i] = s;
    }
    const Gi = matInverse(G, m);
    if (Gi) {
      const l = new Float64Array(m);
      let sum = 0;
      let ok = true;
      for (let i = 0; i < m; i += 1) {
        let s = 0;
        for (let j = 0; j < m; j += 1) s += Gi[i * m + j] * b[j];
        l[i] = s;
        sum += s;
        if (s < -1e-12) ok = false;
      }
      if (ok && sum <= 1 + 1e-12) {
        const p = new Float64Array(n);
        for (let a = 0; a < n; a += 1) {
          let s = pts[0][a];
          for (let i = 0; i < m; i += 1) s += l[i] * E[i][a];
          p[a] = s;
        }
        return p;
      }
    }
    let best = null;
    let bestD = Infinity;
    for (let drop = 0; drop <= m; drop += 1) {
      const sub = pts.filter((_, i) => i !== drop);
      const p = nearestOnSimplex(sub, q, n);
      let d = 0;
      for (let a = 0; a < n; a += 1) d += (p[a] - q[a]) * (p[a] - q[a]);
      if (d < bestD) {
        bestD = d;
        best = p;
      }
    }
    return best;
  }
  /**
   * The unit normal of a surface simplex of `n` vertices. Each component is a
   * signed minor of the matrix of the `n - 1` edge vectors. This is the
   * general form of the 3D cross product of two edges.
   *
   * The function makes the normal point away from the origin. Thus the origin
   * of the mesh must be inside the body.
   *
   * @returns {?Float64Array} the normal, of length `n`, or null when the
   *   simplex is flat
   */
  function facetNormal(pts, n) {
    const E = [];
    for (let i = 1; i < n; i += 1) {
      const e = new Float64Array(n);
      for (let a = 0; a < n; a += 1) e[a] = pts[i][a] - pts[0][a];
      E.push(e);
    }
    const N = new Float64Array(n);
    const M = new Float64Array((n - 1) * (n - 1));
    for (let j = 0; j < n; j += 1) {
      let c = 0;
      for (let a = 0; a < n; a += 1) {
        if (a === j) continue;
        for (let i = 0; i < n - 1; i += 1) M[i * (n - 1) + c] = E[i][a];
        c += 1;
      }
      N[j] = (j % 2 ? -1 : 1) * det(M, n - 1);
    }
    let ln = 0;
    for (let a = 0; a < n; a += 1) ln += N[a] * N[a];
    ln = Math.sqrt(ln);
    if (ln < 1e-14) return null;
    let dot = 0;
    for (let a = 0; a < n; a += 1) {
      N[a] /= ln;
      dot += N[a] * pts[0][a];
    }
    if (dot < 0) for (let a = 0; a < n; a += 1) N[a] = -N[a];
    return N;
  }
  /**
   * The determinant of a small matrix of `n` by `n`, by the minors of the
   * first row. It has a direct formula for `n` of 1, 2 and 3.
   *
   * `matDet` is faster for a large matrix, but `facetNormal` calls this one
   * many times on a very small matrix.
   */
  function det(M, n) {
    if (n === 1) return M[0];
    if (n === 2) return M[0] * M[3] - M[1] * M[2];
    if (n === 3) {
      return M[0] * (M[4] * M[8] - M[5] * M[7]) - M[1] * (M[3] * M[8] - M[5] * M[6]) + M[2] * (M[3] * M[7] - M[4] * M[6]);
    }
    let s = 0;
    const S = new Float64Array((n - 1) * (n - 1));
    for (let j = 0; j < n; j += 1) {
      let c = 0;
      for (let a = 0; a < n; a += 1) {
        if (a === j) continue;
        for (let i = 1; i < n; i += 1) S[(i - 1) * (n - 1) + c] = M[i * n + a];
        c += 1;
      }
      s += (j % 2 ? -1 : 1) * M[j] * det(S, n - 1);
    }
    return s;
  }
  /**
   * The point of the surface of a mesh that is nearest to `q`, and the
   * distance to it. It looks at each surface simplex in turn.
   *
   * It also says if `q` is inside the body: `q` is inside when it is on the
   * inner side of the plane of each simplex.
   *
   * @param {Float64Array} vertices the vertices, `n` numbers for each
   * @param {ArrayLike<number>} cells the simplices, `n` indices for each
   * @param {Float64Array} q the point, of length `n`, in the body frame
   * @param {?Float64Array[]} normals the normals from `meshNormals`
   * @returns {{point: Float64Array, distance: number, inside: boolean}}
   */
  function nearestOnMesh(n, vertices, cells, q, normals) {
    const count = cells.length / n;
    let best = null;
    let bestD = Infinity;
    let inside = true;
    const pts = new Array(n);
    for (let c = 0; c < count; c += 1) {
      for (let i = 0; i < n; i += 1) {
        const v = cells[c * n + i] * n;
        pts[i] = vertices.subarray(v, v + n);
      }
      const p = nearestOnSimplex(pts.slice(), q, n);
      let d = 0;
      for (let a = 0; a < n; a += 1) d += (p[a] - q[a]) * (p[a] - q[a]);
      if (d < bestD) {
        bestD = d;
        best = p;
      }
      const N = normals[c];
      if (N) {
        let s = 0;
        for (let a = 0; a < n; a += 1) s += N[a] * (q[a] - pts[0][a]);
        if (s > 0) inside = false;
      }
    }
    return { point: best, distance: Math.sqrt(bestD), inside };
  }
  /**
   * The normal of each surface simplex. A member is null for a flat simplex.
   * `convexSphere` keeps the result in a cache, because it does not change.
   */
  function meshNormals(n, vertices, cells) {
    const count = cells.length / n;
    const out = new Array(count);
    const pts = new Array(n);
    for (let c = 0; c < count; c += 1) {
      for (let i = 0; i < n; i += 1) {
        const v = cells[c * n + i] * n;
        pts[i] = vertices.subarray(v, v + n);
      }
      out[c] = facetNormal(pts, n);
    }
    return out;
  }

  // src/nd/detect/collide.js
  //
  // The narrow phase. Two bodies come in, and a list of contacts goes out.
  // `collide()` at the end of the module sends each pair of shape types to its
  // own test. See ND-PHYSICS.md, A9 and A10.
  /** The separating axis tables of each `D`. */
  var tableCache = /* @__PURE__ */ new WeakMap();
  /** All of the groups of `m` members of `list`, in order. */
  function combinations(list, m) {
    if (m === 0) return [[]];
    if (m > list.length) return [];
    const out = [];
    for (let i = 0; i <= list.length - m; i += 1) {
      for (const rest of combinations(list.slice(i + 1), m - 1)) out.push([list[i]].concat(rest));
    }
    return out;
  }
  /**
   * The list of the axis groups of the separating axis test of two boxes.
   *
   * In 3 dimensions the edge axes are the cross products of one edge of `a`
   * and one edge of `b`. In `n` dimensions there is no cross product. The
   * general form: take `ma` axes of `a` and `mb` axes of `b`, with
   * `ma + mb = n - 1`, wedge them together, then take the dual. That gives a
   * vector normal to all of them.
   *
   * The count of the groups grows fast with `n`. In 4 dimensions there are
   * two families, `(1, 2)` and `(2, 1)`, and 48 groups. See ND-PHYSICS.md, B5.
   *
   * The result is in a cache, because it does not change.
   *
   * @returns {number[][][]} each member is a pair: the axes of `a`, the axes of `b`
   */
  function satTable(D) {
    let T = tableCache.get(D);
    if (T) return T;
    const { n } = D;
    const axes = [];
    for (let i = 0; i < n; i += 1) axes.push(i);
    T = [];
    for (let ma = 1; ma <= n - 2; ma += 1) {
      const mb = n - 1 - ma;
      if (mb < 1) continue;
      for (const A of combinations(axes, ma)) {
        for (const B of combinations(axes, mb)) T.push([A, B]);
      }
    }
    tableCache.set(D, T);
    return T;
  }
  /**
   * The `n - 1` unit vectors that are normal to `normal`, and normal to each
   * other. The friction of a contact acts in these directions. In 3
   * dimensions there are 2 of them, and in 4 dimensions there are 3.
   *
   * The function takes the axis with the smallest part along the normal
   * first. That axis is the farthest from the normal, thus it gives the most
   * stable result. It then removes the parts along the normal and along each
   * tangent that it already has.
   *
   * @param {Float64Array} normal a unit vector of length `n`
   * @returns {Float64Array[]} the tangents, each of length `n`
   */
  function tangentBasis(D, normal) {
    const { n } = D;
    const order = [];
    for (let i = 0; i < n; i += 1) order.push(i);
    order.sort((a, b) => Math.abs(normal[a]) - Math.abs(normal[b]));
    const T = [];
    for (const idx of order) {
      if (T.length === n - 1) break;
      const t = new Float64Array(n);
      t[idx] = 1;
      let d = 0;
      for (let i = 0; i < n; i += 1) d += t[i] * normal[i];
      for (let i = 0; i < n; i += 1) t[i] -= d * normal[i];
      for (const u of T) {
        let e = 0;
        for (let i = 0; i < n; i += 1) e += t[i] * u[i];
        for (let i = 0; i < n; i += 1) t[i] -= e * u[i];
      }
      let ln = 0;
      for (let i = 0; i < n; i += 1) ln += t[i] * t[i];
      ln = Math.sqrt(ln);
      if (ln < 1e-6) continue;
      for (let i = 0; i < n; i += 1) t[i] /= ln;
      T.push(t);
    }
    return T;
  }
  /**
   * One point of contact between two bodies.
   *
   * The narrow phase fills `a`, `b`, `normal`, `point` and `depth`. The solver
   * fills all of the other fields in `prepareContact`.
   *
   * @typedef {object} Contact
   * @property {Body} a the first body
   * @property {Body} b the second body
   * @property {Float64Array} normal the unit normal, of length `n`. It points
   *   from `a` to `b`.
   * @property {Float64Array} point the point of contact in the world frame, `n`
   * @property {number} depth how deep the two bodies are one in the other. A
   *   negative value shows a gap, and the solver then only stops the two
   *   bodies from coming together.
   * @property {Float64Array} rA the offset from the center of mass of `a`, `n`
   * @property {Float64Array} rB the offset from the center of mass of `b`, `n`
   * @property {number} normalImpulse the impulse along the normal so far. It
   *   is never negative, because a contact can push but it cannot pull.
   * @property {?Float64Array} tangentImpulse the impulse along each tangent
   * @property {?Float64Array[]} tangents the tangents from `tangentBasis`
   * @property {number} kn the effective mass along the normal
   * @property {?number[]} kt the effective mass along each tangent
   * @property {number} target the velocity along the normal that the solver
   *   works for. It holds the correction of the depth and the bounce.
   */
  /** A new contact. Only the fields of the narrow phase have a value. */
  function makeContact(D, a, b, normal, point, depth) {
    return {
      a,
      b,
      normal: Float64Array.from(normal),
      point: Float64Array.from(point),
      depth,
      rA: new Float64Array(D.n),
      rB: new Float64Array(D.n),
      normalImpulse: 0,
      tangentImpulse: null,
      tangents: null,
      kn: 0,
      kt: null,
      target: 0
    };
  }
  /**
   * Two hyperspheres. They touch when the distance between the two centers is
   * less than the sum of the two radii. The normal is the line between the two
   * centers. When the two centers are at the same point, it uses the axis 0.
   */
  function sphereSphere(D, a, b, out) {
    const { n } = D;
    const d = new Float64Array(n);
    let len = 0;
    for (let i = 0; i < n; i += 1) {
      d[i] = b.x[i] - a.x[i];
      len += d[i] * d[i];
    }
    len = Math.sqrt(len);
    const sum = a.shape.radius + b.shape.radius;
    const depth = sum - len;
    if (depth < -0) return;
    if (len < 1e-12) {
      d[0] = 1;
      len = 1;
    }
    for (let i = 0; i < n; i += 1) d[i] /= len;
    const p = new Float64Array(n);
    for (let i = 0; i < n; i += 1) p[i] = a.x[i] + d[i] * (a.shape.radius - depth / 2);
    out.push(makeContact(D, a, b, d, p, depth));
  }
  /**
   * A box and a hypersphere. It puts the center of the sphere into the frame
   * of the box, then it holds each component between the half extents. That
   * gives the nearest point of the box.
   *
   * When the center is inside the box, there is no direction from the nearest
   * point. The function then finds the nearest face, and it uses the normal of
   * that face.
   *
   * @param {boolean} flip true when the sphere is the body `a` of the pair.
   *   The function then turns the normal around, thus the normal always points
   *   from `a` to `b`.
   */
  function boxSphere(D, bx, s, out, flip) {
    const { n } = D;
    const local = bx.worldToLocal(s.x);
    const h = bx.shape.halfExtents;
    const clamped = new Float64Array(n);
    let inside = true;
    for (let i = 0; i < n; i += 1) {
      clamped[i] = Math.max(-h[i], Math.min(h[i], local[i]));
      if (clamped[i] !== local[i]) inside = false;
    }
    const normalLocal = new Float64Array(n);
    let dist2 = 0;
    for (let i = 0; i < n; i += 1) {
      normalLocal[i] = local[i] - clamped[i];
      dist2 += normalLocal[i] * normalLocal[i];
    }
    let depth;
    if (inside) {
      let best = Infinity;
      let axis = 0;
      let sign = 1;
      for (let i = 0; i < n; i += 1) {
        const d = h[i] - Math.abs(local[i]);
        if (d < best) {
          best = d;
          axis = i;
          sign = local[i] >= 0 ? 1 : -1;
        }
      }
      normalLocal.fill(0);
      normalLocal[axis] = sign;
      clamped[axis] = sign * h[axis];
      depth = s.shape.radius + best;
    } else {
      const dist = Math.sqrt(dist2);
      depth = s.shape.radius - dist;
      if (depth < 0) return;
      for (let i = 0; i < n; i += 1) normalLocal[i] /= dist;
    }
    const normal = bx.localToWorldDir(normalLocal);
    const point = bx.localToWorld(clamped);
    if (flip) {
      for (let i = 0; i < n; i += 1) normal[i] = -normal[i];
      out.push(makeContact(D, s, bx, normal, point, depth));
    } else {
      out.push(makeContact(D, bx, s, normal, point, depth));
    }
  }
  /**
   * The points of a torus that can touch a half space. A torus lies flat on
   * the ground along a full circle, and not at one point. One contact would
   * let it fall over.
   *
   * The function gives the support point, and then 12 points around the major
   * circle. Each one is at the distance `r` from the circle, in the direction
   * `dir`. `halfSpaceOther` drops the points that are too far away.
   *
   * @param {Float64Array} dir the direction into the ground, in the world frame
   * @returns {Float64Array[]} the points, in the world frame
   */
  function torusGroundPoints(D, body, dir) {
    const { n } = D;
    const sh = body.shape;
    const local = body.worldToLocalDir(dir);
    const pts = [];
    pts.push(sh.support(local));
    const [i, j] = sh.plane;
    const R = sh.majorRadius;
    const r = sh.minorRadius;
    let ln = 0;
    for (let a = 0; a < n; a += 1) ln += local[a] * local[a];
    ln = Math.sqrt(ln) || 1;
    const SAMPLES = 12;
    for (let t = 0; t < SAMPLES; t += 1) {
      const th = 2 * Math.PI * t / SAMPLES;
      const p = new Float64Array(n);
      for (let a = 0; a < n; a += 1) p[a] = local[a] / ln * r;
      p[i] += R * Math.cos(th);
      p[j] += R * Math.sin(th);
      pts.push(p);
    }
    return pts.map((p) => body.localToWorld(p));
  }
  /**
   * A torus and a hypersphere. The surface of a torus is the set of the points
   * at the distance `r` from the major circle. Thus the test is simple: find
   * the nearest point of the major circle with `corePoint`, then compare the
   * distance to `r` plus the radius of the sphere.
   */
  function torusSphere(D, to, s, out, flip) {
    const { n } = D;
    const sh = to.shape;
    const local = to.worldToLocal(s.x);
    const core = sh.corePoint(local);
    const v = new Float64Array(n);
    let dv = 0;
    for (let a = 0; a < n; a += 1) {
      v[a] = local[a] - core[a];
      dv += v[a] * v[a];
    }
    dv = Math.sqrt(dv);
    if (dv < 1e-12) {
      v.fill(0);
      v[sh.plane[0]] = 1;
      dv = 1;
    } else {
      for (let a = 0; a < n; a += 1) v[a] /= dv;
    }
    const depth = s.shape.radius + sh.minorRadius - dv;
    if (depth < 0) return;
    const surface = new Float64Array(n);
    for (let a = 0; a < n; a += 1) surface[a] = core[a] + v[a] * sh.minorRadius;
    const normal = to.localToWorldDir(v);
    const point = to.localToWorld(surface);
    if (flip) {
      for (let a = 0; a < n; a += 1) normal[a] = -normal[a];
      out.push(makeContact(D, s, to, normal, point, depth));
    } else {
      out.push(makeContact(D, to, s, normal, point, depth));
    }
  }
  /** The facet normals of each convex shape. `meshNormals` costs much. */
  var normalCache = /* @__PURE__ */ new WeakMap();
  /**
   * A convex mesh and a hypersphere. It puts the center of the sphere into the
   * frame of the mesh, then `nearestOnMesh` gives the nearest point of the
   * surface.
   *
   * The direction of the normal depends on the side. When the center is
   * outside, the normal points from the surface to the center. When it is
   * inside, the normal points the other way, and the depth is then the radius
   * plus the distance.
   */
  function convexSphere(D, cv, s, out, flip) {
    const { n } = D;
    const sh = cv.shape;
    let normals = normalCache.get(sh);
    if (!normals) {
      normals = meshNormals(n, sh.vertices, sh.cells);
      normalCache.set(sh, normals);
    }
    const local = cv.worldToLocal(s.x);
    const near = nearestOnMesh(n, sh.vertices, sh.cells, local, normals);
    const v = new Float64Array(n);
    let depth;
    if (near.inside) {
      for (let a = 0; a < n; a += 1) v[a] = near.point[a] - local[a];
      depth = s.shape.radius + near.distance;
    } else {
      for (let a = 0; a < n; a += 1) v[a] = local[a] - near.point[a];
      depth = s.shape.radius - near.distance;
      if (depth < 0) return;
    }
    let ln = 0;
    for (let a = 0; a < n; a += 1) ln += v[a] * v[a];
    ln = Math.sqrt(ln);
    if (ln < 1e-12) {
      v.fill(0);
      v[0] = 1;
    } else {
      for (let a = 0; a < n; a += 1) v[a] /= ln;
    }
    const normal = cv.localToWorldDir(v);
    const point = cv.localToWorld(near.point);
    if (flip) {
      for (let a = 0; a < n; a += 1) normal[a] = -normal[a];
      out.push(makeContact(D, s, cv, normal, point, depth));
    } else {
      out.push(makeContact(D, cv, s, normal, point, depth));
    }
  }
  /** True after the first warning. One warning for each page is sufficient. */
  var warned = false;
  /** Writes one warning about a pair of shapes that has no test. */
  function warnNoPair(ta, tb) {
    if (warned) return;
    warned = true;
    if (typeof console !== "undefined" && console.warn) {
      console.warn(`PhysiN: no collision test for the pair (${ta}, ${tb}). A torus and a convex mesh touch a half space and a hypersphere only.`);
    }
  }
  /**
   * A half space and any other shape. The normal of the contact is always the
   * normal of the half space, thus the test is only a distance.
   *
   * There are three conditions:
   *   - a torus: it uses the points of `torusGroundPoints`.
   *   - a hypersphere: one point, at the radius below the center.
   *   - a box or a convex mesh: each vertex.
   * It drops a point that is farther than `margin` above the plane.
   */
  function halfSpaceOther(D, hs, other, out, flip, margin) {
    const { n } = D;
    const nrm = hs.shape.normal;
    const off = hs.shape.offset;
    const push = (point, depth) => {
      if (depth < -margin) return;
      if (flip) {
        const inv = new Float64Array(n);
        for (let i = 0; i < n; i += 1) inv[i] = -nrm[i];
        out.push(makeContact(D, other, hs, inv, point, depth));
      } else {
        out.push(makeContact(D, hs, other, nrm, point, depth));
      }
    };
    if (other.shape.type === "torus") {
      const dir = new Float64Array(n);
      for (let i = 0; i < n; i += 1) dir[i] = -nrm[i];
      for (const p of torusGroundPoints(D, other, dir)) {
        let d = -off;
        for (let i = 0; i < n; i += 1) d += nrm[i] * p[i];
        push(p, -d);
      }
      return;
    }
    if (other.shape.type === "sphere") {
      let d = -off;
      for (let i = 0; i < n; i += 1) d += nrm[i] * other.x[i];
      const depth = other.shape.radius - d;
      const p = new Float64Array(n);
      for (let i = 0; i < n; i += 1) p[i] = other.x[i] - nrm[i] * other.shape.radius;
      push(p, depth);
      return;
    }
    const verts = boxVertices(D, other);
    for (const p of verts) {
      let d = -off;
      for (let i = 0; i < n; i += 1) d += nrm[i] * p[i];
      push(p, -d);
    }
  }
  /** The body frame vertices of each shape. A box in 4 dimensions has 16. */
  var vertexCache = /* @__PURE__ */ new WeakMap();
  /**
   * The vertices of a box or of a convex mesh, in the world frame. The body
   * frame vertices are in a cache, because they do not change. A box of `n`
   * dimensions has `2^n` of them.
   */
  function boxVertices(D, body) {
    const { n } = D;
    let local = vertexCache.get(body.shape);
    if (!local) {
      local = [];
      if (body.shape.type === "box") {
        const h = body.shape.halfExtents;
        for (let m = 0; m < 1 << n; m += 1) {
          const v = new Float64Array(n);
          for (let i = 0; i < n; i += 1) v[i] = m >> i & 1 ? h[i] : -h[i];
          local.push(v);
        }
      } else if (body.shape.vertices) {
        const V = body.shape.vertices;
        for (let m = 0; m < V.length / n; m += 1) local.push(V.subarray(m * n, m * n + n));
      }
      vertexCache.set(body.shape, local);
    }
    return local.map((v) => body.localToWorld(v));
  }
  /**
   * How far a body goes from its center along `axis`. For a sphere this is the
   * radius. For a box it is the sum of `|axis . column_j| * h_j` over the
   * columns of the rotation matrix.
   */
  function boxExtent(D, body, axis) {
    const { n } = D;
    if (body.shape.type === "sphere") return body.shape.radius;
    const h = body.shape.halfExtents;
    let e = 0;
    for (let j = 0; j < n; j += 1) {
      let d = 0;
      for (let i = 0; i < n; i += 1) d += body.Rm[i * n + j] * axis[i];
      e += Math.abs(d) * h[j];
    }
    return e;
  }
  /**
   * The separating axis test of two boxes. See ND-PHYSICS.md, A10 and B5.
   *
   * Two convex bodies are apart when there is one axis on which their two
   * shadows do not meet. For two boxes it is sufficient to test a small set of
   * axes: the `n` axes of `a`, the `n` axes of `b`, and the axes of
   * `satTable(D)`. The last group comes from the wedge product of `n - 1`
   * axes, and then the dual.
   *
   * The function stops at the first axis with no overlap: the two bodies are
   * then apart. If all of the axes overlap, it gives the axis of the smallest
   * overlap. That axis is the normal of the contact.
   *
   * `owner` says where that best axis came from: `a`, `b`, or `mixed` for an
   * axis of the table. `boxBox` uses it to choose the face to work with.
   *
   * @returns {?{axis: Float64Array, overlap: number, owner: string}} the
   *   result, or null when the two boxes are apart. The axis points from `a`
   *   to `b`.
   */
  function boxBoxAxis(D, a, b) {
    const { n } = D;
    const T = satTable(D);
    const delta = new Float64Array(n);
    for (let i = 0; i < n; i += 1) delta[i] = b.x[i] - a.x[i];
    let bestOverlap = Infinity;
    let bestOwner = "mixed";
    const bestAxis = new Float64Array(n);
    const axis = new Float64Array(n);
    let owner = "mixed";
    const test = () => {
      let ln = 0;
      for (let i = 0; i < n; i += 1) ln += axis[i] * axis[i];
      if (ln < 1e-12) return true;
      ln = Math.sqrt(ln);
      for (let i = 0; i < n; i += 1) axis[i] /= ln;
      let d = 0;
      for (let i = 0; i < n; i += 1) d += axis[i] * delta[i];
      const overlap = boxExtent(D, a, axis) + boxExtent(D, b, axis) - Math.abs(d);
      if (overlap < 0) return false;
      if (overlap < bestOverlap) {
        bestOverlap = overlap;
        bestOwner = owner;
        const s = d < 0 ? -1 : 1;
        for (let i = 0; i < n; i += 1) bestAxis[i] = axis[i] * s;
      }
      return true;
    };
    for (const [body, tag] of [[a, "a"], [b, "b"]]) {
      owner = tag;
      for (let j = 0; j < n; j += 1) {
        for (let i = 0; i < n; i += 1) axis[i] = body.Rm[i * n + j];
        if (!test()) return null;
      }
    }
    owner = "mixed";
    const w1 = mvZero(D);
    const w2 = mvZero(D);
    const w3 = mvZero(D);
    const col = new Float64Array(n);
    for (const [A, B] of T) {
      let acc = null;
      let dst = w2;
      for (const [body, list] of [[a, A], [b, B]]) {
        for (const j of list) {
          for (let i = 0; i < n; i += 1) col[i] = body.Rm[i * n + j];
          if (acc === null) {
            acc = mvFromVector(D, col, w1);
          } else {
            mvFromVector(D, col, w3);
            mvWedge(D, acc, w3, dst);
            const prev = acc;
            acc = dst;
            dst = prev;
          }
        }
      }
      if (!acc) continue;
      const du = mvDual(D, acc, w3);
      mvToVector(D, du, axis);
      if (!test()) return null;
    }
    if (!Number.isFinite(bestOverlap)) return null;
    return { axis: bestAxis, overlap: bestOverlap, owner: bestOwner };
  }
  /**
   * Two boxes. `boxBoxAxis` gives the normal, and then this function finds the
   * points.
   *
   * The method: take the face of `a` that faces `b`, and keep each vertex of
   * `b` that is near that face. Do the same with the parts changed around.
   * A box that rests flat on another box then gives `2^(n-1)` points, and the
   * stack is stable.
   *
   * When no vertex is near, the contact is edge on edge. The function then
   * uses the middle of the two support points.
   *
   * At the end it puts the points in order of the depth, and it keeps the
   * `maxContacts` deepest.
   */
  function boxBox(D, a, b, out, margin, maxContacts) {
    const res = boxBoxAxis(D, a, b);
    if (!res) return;
    const { n } = D;
    const nrm = res.axis;
    const dot = (p) => {
      let d = 0;
      for (let i = 0; i < n; i += 1) d += p[i] * nrm[i];
      return d;
    };
    const ca = dot(a.x);
    const cb = dot(b.x);
    const found = [];
    if (res.owner === "a" || res.owner === "mixed") {
      const faceA = ca + boxExtent(D, a, nrm);
      for (const p of boxVertices(D, b)) {
        const depth = faceA - dot(p);
        if (depth > -margin && depth <= res.overlap + margin) found.push([p, depth]);
      }
    }
    if (res.owner === "b" || res.owner === "mixed") {
      const faceB = cb - boxExtent(D, b, nrm);
      for (const p of boxVertices(D, a)) {
        const depth = dot(p) - faceB;
        if (depth > -margin && depth <= res.overlap + margin) found.push([p, depth]);
      }
    }
    if (found.length === 0) {
      const inv = new Float64Array(n);
      for (let i = 0; i < n; i += 1) inv[i] = -nrm[i];
      const pa = a.localToWorld(a.shape.support(a.worldToLocalDir(nrm)));
      const pb = b.localToWorld(b.shape.support(b.worldToLocalDir(inv)));
      const mid = new Float64Array(n);
      for (let i = 0; i < n; i += 1) mid[i] = (pa[i] + pb[i]) / 2;
      found.push([mid, res.overlap]);
    }
    found.sort((p, q) => q[1] - p[1]);
    for (let i = 0; i < Math.min(found.length, maxContacts); i += 1) {
      out.push(makeContact(D, a, b, nrm, found[i][0], found[i][1]));
    }
  }
  /**
   * The narrow phase of one pair of bodies. It adds a contact to `out` for
   * each point where the two bodies touch. It adds nothing when they do not.
   *
   * The table of the pairs. A "-" shows a pair with no test:
   *
   *              halfspace  sphere  box     convex  torus
   *   halfspace  none       yes     yes     yes     yes
   *   sphere                yes     yes     yes     yes
   *   box                           yes     -       -
   *   convex                                -       -
   *   torus                                         -
   *
   * A torus and a convex mesh touch a half space and a hypersphere only. A
   * pair with no test writes one warning, and then it does nothing. Two half
   * spaces do nothing, and they write no warning.
   *
   * @param {Body} a the first body
   * @param {Body} b the second body
   * @param {Contact[]} out the list to add to
   * @param {object} [params] `contactMargin` and `maxContacts` of the world
   */
  function collide(D, a, b, out, params = {}) {
    const margin = params.contactMargin !== void 0 ? params.contactMargin : 0.02;
    const maxContacts = params.maxContacts || 1 << D.n - 1;
    const ta = a.shape.type;
    const tb = b.shape.type;
    if (ta === "halfspace" && tb === "halfspace") return;
    if (ta === "halfspace") return halfSpaceOther(D, a, b, out, false, margin);
    if (tb === "halfspace") return halfSpaceOther(D, b, a, out, true, margin);
    if (ta === "sphere" && tb === "sphere") return sphereSphere(D, a, b, out);
    if (ta === "convex" && tb === "sphere") return convexSphere(D, a, b, out, false);
    if (tb === "convex" && ta === "sphere") return convexSphere(D, b, a, out, true);
    if (ta === "convex" || tb === "convex") return warnNoPair(ta, tb);
    if (ta === "torus" && tb === "sphere") return torusSphere(D, a, b, out, false);
    if (tb === "torus" && ta === "sphere") return torusSphere(D, b, a, out, true);
    if (ta === "torus" || tb === "torus") return warnNoPair(ta, tb);
    if (ta === "sphere" && tb !== "sphere") return boxSphere(D, b, a, out, true);
    if (tb === "sphere" && ta !== "sphere") return boxSphere(D, a, b, out, false);
    return boxBox(D, a, b, out, margin, maxContacts);
  }

  // src/nd/integrate/integrator.js
  //
  // The time integration, in two halves. `integrateVelocities` runs first, and
  // the solver runs between the two. `integratePositions` runs last. See
  // ND-PHYSICS.md, B1.
  /**
   * The first half of the step: the forces change the velocities.
   *
   *   v += dt (F / m + g)
   *   L += dt tau
   *   w  = I'^-1 L
   *
   * The damping is not an exponential. It is the implicit form
   * `1 / (1 + dt c)`. That form is stable at any step length.
   *
   * `linearFactor` and `angularFactor` multiply the change, thus you can hold
   * a body on one axis or in one plane.
   *
   * The gyroscopic term is NOT here. It is in `integratePositions`, because
   * the state holds `L` and not `w`. See `gyroscopicSpin`.
   * A static body or a sleeping body does not change.
   *
   * @param {number} dt the length of the step in seconds
   * @param {Float64Array} gravity the acceleration of gravity, of length `n`
   * @param {object} opts the params of the world
   */
  function integrateVelocities(D, body, dt, gravity, opts) {
    if (body.isStatic || body.sleeping) return;
    const { n, k } = D;
    for (let i = 0; i < n; i += 1) {
      body.v[i] += dt * (body.invMass * body.force[i] + gravity[i]) * body.linearFactor[i];
    }
    if (body.linearDamping > 0) {
      const f = 1 / (1 + dt * body.linearDamping);
      for (let i = 0; i < n; i += 1) body.v[i] *= f;
    }
    for (let p = 0; p < k; p += 1) body.L[p] += dt * body.torque[p] * body.angularFactor[p];
    if (body.angularDamping > 0) {
      const f = 1 / (1 + dt * body.angularDamping);
      for (let p = 0; p < k; p += 1) body.L[p] *= f;
    }
    matVec(body.invInertiaWorld, body.L, k, k, body.w);
  }
  /**
   * The gyroscopic term of the Euler equation. It gives the angular velocity
   * at the END of the step, in the world frame. It changes NO field of the
   * body: the caller uses the answer for the rotor step only.
   *
   * WHY IT CHANGES NOTHING. The state holds `L` in the world frame, and
   * `updateDerived` builds `w` again from `L` with the new rotor. Thus the
   * turn of the frame is already in the state, and a free body keeps its
   * momentum exactly. A write of `L` here would count that turn a second
   * time, and the momentum would turn at the rate `[w, L]`. What is left for
   * this function is the accuracy of the rotor step: `w` at the end of the
   * step, and not at the start.
   *
   * The Euler equation in the body frame is `I dw/dt - w x I w = tau`. Here
   * `x` is the commutator of two bivectors, and not a cross product. The
   * explicit form of this term adds energy and it goes bad. Thus the function
   * solves the implicit form with the Newton method:
   *
   *   f(w2)  = I (w2 - w1) - dt (w2 x I w2) = 0
   *   J      = I + dt [I w2]comm - dt [w2]comm I
   *   w2    -= J^-1 f(w2)
   *
   * `opts.gyroscopicIterations` gives the count of the turns of the loop. One
   * turn is sufficient at a normal step length. The loop stops early when the
   * change is very small, or when `J` is singular.
   *
   * The work is in the body frame, where the inertia tensor is constant. The
   * function changes the frame with `R2` at the start and at the end.
   *
   * See ND-PHYSICS.md, A7, B1 and C5.
   *
   * @returns {Float64Array} the angular velocity in the world frame, of
   *   length `k`
   */
  function gyroscopicSpin(D, body, dt, opts, out) {
    const { k } = D;
    const iters = opts.gyroscopicIterations || 1;
    const w1 = matTVec(body.R2, body.w, k, k);
    const w2 = Float64Array.from(w1);
    const I = body.inertia;
    for (let it = 0; it < iters; it += 1) {
      const Iw = matVec(I, w2, k, k);
      const cross = commutator(D, w2, Iw);
      const f = new Float64Array(k);
      for (let p = 0; p < k; p += 1) {
        let s = 0;
        for (let q = 0; q < k; q += 1) s += I[p * k + q] * (w2[q] - w1[q]);
        f[p] = s - dt * cross[p];
      }
      const CmIw = commutatorMatrix(D, Iw);
      const Cmw = commutatorMatrix(D, w2);
      const CmwI = matMul(Cmw, I, k, k, k);
      const J = matZero(k, k);
      for (let t = 0; t < k * k; t += 1) J[t] = I[t] + dt * CmIw[t] - dt * CmwI[t];
      let Jinv;
      try {
        Jinv = matInverse(J, k);
      } catch (e) {
        break;
      }
      const d = matVec(Jinv, f, k, k);
      let change = 0;
      for (let p = 0; p < k; p += 1) {
        w2[p] -= d[p];
        change += d[p] * d[p];
      }
      if (change < 1e-24) break;
    }
    return matVec(body.R2, w2, k, k, out || new Float64Array(k));
  }
  /**
   * The second half of the step: the velocities change the position and the
   * orientation.
   *
   *   x += dt v
   *   R += -0.5 dt w R      (w R is the geometric product)
   *
   * The rotor equation has the same form as the quaternion equation of a 3D
   * engine. See ND-PHYSICS.md, A5.
   *
   * The `w` of the rotor step comes from `gyroscopicSpin`, if
   * `opts.gyroscopic` is not false. That is the angular velocity at the END
   * of the step. The function runs here, and not in `integrateVelocities`,
   * for two reasons: each impulse builds `body.w` again from `body.L`, thus a
   * value that the first half wrote would go away as soon as the solver
   * touches the body; and the answer must not go into `body.L`. See
   * `gyroscopicSpin`.
   *
   * That step takes the rotor a little away from a true rotation. The
   * function measures the error with `rotorDefect`, and it calls
   * `rotorCorrect` only when the error is more than `opts.rotorTolerance`.
   * This repair is necessary in 4 dimensions: without it the body becomes
   * larger or thinner as it turns. See ND-PHYSICS.md, B2.
   *
   * At the end it calls `body.updateDerived()`, because `R` changed.
   */
  function integratePositions(D, body, dt, opts) {
    if (body.isStatic || body.sleeping) return;
    const { n, r } = D;
    for (let i = 0; i < n; i += 1) body.x[i] += dt * body.v[i];
    const w = opts.gyroscopic !== false ? gyroscopicSpin(D, body, dt, opts) : body.w;
    const W = new Float64Array(r);
    for (let p = 0; p < D.k; p += 1) W[D.biSlot[p]] = w[p];
    const dR = rotorMul(D, W, body.R);
    for (let i = 0; i < r; i += 1) body.R[i] += -0.5 * dt * dR[i];
    const tol = opts.rotorTolerance !== void 0 ? opts.rotorTolerance : 1e-9;
    if (rotorDefect(D, body.R) > tol) {
      const C = rotorCorrect(D, body.R);
      body.R.set(C);
    }
    body.updateDerived();
  }

  // src/nd/resolve/solver.js
  //
  // The contact solver. It uses sequential impulses: it takes each contact in
  // turn, it applies an impulse that makes that one contact correct, and it
  // repeats for `params.iterations` turns. The result comes near to the true
  // answer of all of the contacts together.
  //
  // The order of the work in each step:
  //   1. `prepareContact`  builds the effective mass and the target velocity.
  //   2. `warmStart`       applies the impulse of the last step again.
  //   3. `solveContact`    the main loop, many turns.
  //   4. `shockPropagation` an extra pass for a tall stack.
  //
  // See ND-PHYSICS.md, A11 and B6.
  /** The scratch vectors of each `D`. They keep the solver free of garbage. */
  var scratchCache = /* @__PURE__ */ new WeakMap();
  /** The scratch vectors of the algebra `D`. Do not hold the result. */
  function scratch2(D) {
    let s = scratchCache.get(D);
    if (!s) {
      s = {
        ua: new Float64Array(D.n),
        ub: new Float64Array(D.n),
        u: new Float64Array(D.n),
        j: new Float64Array(D.n),
        neg: new Float64Array(D.n),
        bi: new Float64Array(D.k),
        iv: new Float64Array(D.k)
      };
      scratchCache.set(D, s);
    }
    return s;
  }
  /**
   * The inverse of the mass that a body shows at the offset `r`, along the
   * direction `dir`:
   *
   *   1 / m + (r ^ dir) . I'^-1 (r ^ dir)
   *
   * A large value means that the body moves easily there. The impulse of a
   * contact is the change of velocity divided by the sum of the two values.
   * See ND-PHYSICS.md, A11.
   *
   * @param {boolean} treatStatic true makes the body act as a static body,
   *   thus the result is 0. `shockPropagation` uses this.
   */
  function effectiveMass(D, body, r, dir, treatStatic) {
    if (treatStatic || body.isStatic) return 0;
    const { k } = D;
    const s = scratch2(D);
    const bi = wedgeVec(D, r, dir, s.bi);
    const iv = matVec(body.invInertiaWorld, bi, k, k, s.iv);
    let m = body.invMass;
    for (let p = 0; p < k; p += 1) m += iv[p] * bi[p];
    return m;
  }
  /**
   * The velocity of the point of `b` less the velocity of the point of `a`, at
   * the contact. It writes into a scratch vector when you do not give `out`.
   */
  function relativeVelocity(D, c, out) {
    const { n } = D;
    const s = scratch2(D);
    const ua = c.a.pointVelocity(c.rA, s.ua);
    const ub = c.b.pointVelocity(c.rB, s.ub);
    const u = out || s.u;
    for (let i = 0; i < n; i += 1) u[i] = ub[i] - ua[i];
    return u;
  }
  /**
   * Builds the fields of a contact that the solver needs. Call this one time
   * in each step, before `solveContact`.
   *
   * It builds the offsets `rA` and `rB`, the effective mass `kn` along the
   * normal, the tangents and their effective mass `kt`, and the friction and
   * the restitution of the pair. The friction of a pair is the geometric mean
   * of the two values, and the restitution is the smaller of the two.
   *
   * The target velocity `target` holds two parts:
   *   - The correction of the depth: `beta (depth - slop) / dt`. This pushes
   *     the two bodies apart slowly. `slop` is a small depth that the solver
   *     accepts, and it stops the contact from shaking.
   *   - The bounce: `-restitution * vn`, but only when the two bodies come
   *     together faster than `restitutionThreshold`. Without that limit a body
   *     that rests would never stop bouncing.
   * A contact with a gap, thus a negative depth, only stops the two bodies
   * from coming together in this step.
   *
   * @param {Contact} c the contact
   * @param {number} dt the length of the step
   * @param {object} params the params of the world
   * @param {boolean} [aStatic] true makes `a` act as a static body
   * @param {boolean} [bStatic] true makes `b` act as a static body
   */
  function prepareContact(D, c, dt, params, aStatic = false, bStatic = false) {
    const { n } = D;
    for (let i = 0; i < n; i += 1) {
      c.rA[i] = c.point[i] - c.a.x[i];
      c.rB[i] = c.point[i] - c.b.x[i];
    }
    c.aStatic = aStatic;
    c.bStatic = bStatic;
    c.kn = effectiveMass(D, c.a, c.rA, c.normal, aStatic) + effectiveMass(D, c.b, c.rB, c.normal, bStatic);
    c.tangents = tangentBasis(D, c.normal);
    c.kt = c.tangents.map((t) => effectiveMass(D, c.a, c.rA, t, aStatic) + effectiveMass(D, c.b, c.rB, t, bStatic));
    if (!c.tangentImpulse || c.tangentImpulse.length !== c.tangents.length) {
      c.tangentImpulse = new Float64Array(c.tangents.length);
    }
    const u = relativeVelocity(D, c);
    let vn = 0;
    for (let i = 0; i < n; i += 1) vn += u[i] * c.normal[i];
    c.vnInitial = vn;
    const slop = params.penetrationSlop;
    const beta = params.biasFactor;
    const rest = Math.min(c.a.restitution, c.b.restitution);
    c.friction = Math.sqrt(c.a.friction * c.b.friction);
    if (c.depth < 0) {
      c.target = c.depth / dt;
    } else {
      let t = beta * Math.max(0, c.depth - slop) / dt;
      if (vn < -params.restitutionThreshold) t += -rest * vn;
      c.target = t;
    }
  }
  /**
   * Applies the impulse `j` to `b`, and the opposite impulse to `a`. This
   * keeps the total momentum the same. It does not touch a body that the
   * contact treats as static.
   */
  function applyPair(D, c, j) {
    const { n } = D;
    if (!c.aStatic) {
      const neg = scratch2(D).neg;
      for (let i = 0; i < n; i += 1) neg[i] = -j[i];
      c.a.applyImpulse(neg, c.rA, true);
    }
    if (!c.bStatic) c.b.applyImpulse(j, c.rB, true);
  }
  /**
   * Applies the impulse of the last step again, before the main loop. A stack
   * of boxes then holds still with far fewer turns of the loop, because the
   * solver starts near the answer.
   *
   * `World.applyWarmStartCache` finds the impulse of the last step, from the
   * contact at almost the same point.
   *
   * A contact with a gap starts at zero.
   */
  function warmStart(D, c) {
    const { n } = D;
    const buf = scratch2(D).j;
    if (c.depth < 0) {
      c.normalImpulse = 0;
      c.tangentImpulse.fill(0);
      return;
    }
    if (c.normalImpulse === 0 && !c.tangentImpulse.some((x) => x !== 0)) return;
    const j = buf;
    for (let i = 0; i < n; i += 1) j[i] = c.normal[i] * c.normalImpulse;
    for (let t = 0; t < c.tangents.length; t += 1) {
      const tv = c.tangents[t];
      const m = c.tangentImpulse[t];
      for (let i = 0; i < n; i += 1) j[i] += tv[i] * m;
    }
    applyPair(D, c, j);
  }
  /** The friction work buffer of each contact. */
  var wantedCache = /* @__PURE__ */ new WeakMap();
  /**
   * One turn of the solver on one contact. Call it many times.
   *
   * The normal part. The impulse is `-(vn - target) / kn`. The total impulse
   * of the contact can never be negative, because a contact can push but it
   * cannot pull. Thus the function holds the total at 0 or more, and it
   * applies only the change. That clamp on the total, and not on the change,
   * is what makes the sequential impulse method work.
   *
   * The friction part. It finds the impulse that would stop the movement along
   * each tangent. It then holds the length of that vector at
   * `friction * normalImpulse`. This is the friction cone of Coulomb. The
   * clamp is on the length of the whole vector, and not on each tangent one by
   * one. Thus the friction does not depend on the choice of the tangents.
   *
   * See ND-PHYSICS.md, A11.
   */
  function solveContact(D, c) {
    const { n } = D;
    const j = scratch2(D).j;
    j.fill(0);
    if (c.kn > 0) {
      const u = relativeVelocity(D, c);
      let vn = 0;
      for (let i = 0; i < n; i += 1) vn += u[i] * c.normal[i];
      let dJ = -(vn - c.target) / c.kn;
      const old = c.normalImpulse;
      c.normalImpulse = Math.max(0, old + dJ);
      dJ = c.normalImpulse - old;
      if (dJ !== 0) {
        for (let i = 0; i < n; i += 1) j[i] = c.normal[i] * dJ;
        applyPair(D, c, j);
      }
    }
    const limit = c.friction * c.normalImpulse;
    const u2 = relativeVelocity(D, c);
    let wanted = wantedCache.get(c);
    if (!wanted || wanted.length !== c.tangents.length) {
      wanted = new Float64Array(c.tangents.length);
      wantedCache.set(c, wanted);
    }
    for (let t = 0; t < c.tangents.length; t += 1) {
      if (c.kt[t] <= 0) {
        wanted[t] = c.tangentImpulse[t];
        continue;
      }
      let vt = 0;
      for (let i = 0; i < n; i += 1) vt += u2[i] * c.tangents[t][i];
      wanted[t] = c.tangentImpulse[t] - vt / c.kt[t];
    }
    let mag = 0;
    for (let t = 0; t < wanted.length; t += 1) mag += wanted[t] * wanted[t];
    mag = Math.sqrt(mag);
    if (mag > limit && mag > 0) {
      const f = limit / mag;
      for (let t = 0; t < wanted.length; t += 1) wanted[t] *= f;
    }
    j.fill(0);
    let any = false;
    for (let t = 0; t < wanted.length; t += 1) {
      const d = wanted[t] - c.tangentImpulse[t];
      c.tangentImpulse[t] = wanted[t];
      if (d !== 0) {
        any = true;
        for (let i = 0; i < n; i += 1) j[i] += c.tangents[t][i] * d;
      }
    }
    if (any) applyPair(D, c, j);
  }
  /**
   * Gives each body its `level`: the count of the contacts between it and the
   * nearest static body. A static body has the level 0, a box on the ground
   * has 1, the box on that box has 2, and so on.
   *
   * The method is a breadth first walk from the static bodies. A body that no
   * contact connects to a static body gets the level 0.
   *
   * `shockPropagation` uses the levels.
   */
  function buildContactGraph(bodies, contacts) {
    const neighbours = /* @__PURE__ */ new Map();
    for (const b of bodies) {
      b.level = b.isStatic ? 0 : Infinity;
      neighbours.set(b, []);
    }
    for (const c of contacts) {
      neighbours.get(c.a).push(c.b);
      neighbours.get(c.b).push(c.a);
    }
    const queue = bodies.filter((b) => b.isStatic);
    let head = 0;
    while (head < queue.length) {
      const b = queue[head];
      head += 1;
      for (const o of neighbours.get(b) || []) {
        if (o.level > b.level + 1) {
          o.level = b.level + 1;
          queue.push(o);
        }
      }
    }
    for (const b of bodies) if (!Number.isFinite(b.level)) b.level = 0;
  }
  /**
   * An extra pass that stops a tall stack from sinking. See ND-PHYSICS.md, B6.
   *
   * A stack of boxes is difficult for the sequential impulse method: the box
   * at the bottom carries all of the weight, and the main loop does not have
   * sufficient turns to hold it.
   *
   * The method: take the contacts in the order of the level, from the ground
   * up. At each contact, treat the body of the lower level as a static body.
   * The weight then goes to the ground in one pass, and it does not go back up.
   *
   * The pass does not keep its impulses for the warm start, because it works
   * with false masses.
   */
  function shockPropagation(D, contacts, dt, params) {
    const sorted = contacts.slice().sort((p, q) => Math.min(p.a.level, p.b.level) - Math.min(q.a.level, q.b.level));
    for (const c of sorted) {
      const aLower = c.a.level <= c.b.level;
      prepareContact(D, c, dt, params, aLower, !aLower);
      c.normalImpulse = 0;
      c.tangentImpulse.fill(0);
      for (let i = 0; i < params.shockIterations; i += 1) solveContact(D, c);
    }
  }

  // src/nd/resolve/constraint.js
  //
  // The constraints and the joints. A constraint holds a list of rows, and each
  // row takes away one degree of freedom.
  //
  // A row is a linear row or an angular row:
  //   - A linear row holds a unit direction `dir` of `n` components, and the
  //     two world offsets `rA` and `rB`. It works on the velocity of a point.
  //   - An angular row holds a unit bivector `axis` of `k` components. It works
  //     on the angular velocity.
  //
  // The five types:
  //   point     `n` linear rows. It holds one point of `a` on one point of `b`.
  //   distance  1 linear row. It holds the length between the two points.
  //   fixed     `n` linear rows and `k` angular rows. It welds the two bodies.
  //   hinge     `n` linear rows and `k - 1` angular rows. One plane stays free.
  //   subspace  the rows that you name. It holds a body in a subspace.
  //
  // A hinge leaves one rotation PLANE free, and not one axis. In 3 dimensions
  // that is 1 plane free of 3. In 4 dimensions it is 1 free of 6, thus a hinge
  // takes away 5 angular degrees of freedom. In 2 dimensions `k` is 1, thus
  // there is no angular row and a hinge is the same as a point joint.
  //
  // The sign: the error and the velocity always measure `b` against `a`. A
  // positive impulse pushes `b` along `dir`, or it turns `b` along `axis`.
  //
  // The method is the method of the contact solver: take each row in turn,
  // apply the impulse that makes that one row correct, and repeat. The solver
  // holds the TOTAL impulse of the row between `lower` and `upper`, and it
  // applies only the change. The total impulse lives between the steps, thus
  // the warm start costs nothing.
  //
  // See ND-PHYSICS.md, A13 and B8.
  /** The scratch vectors of the joints of each `D`. Do not hold the result. */
  var jointScratchCache = /* @__PURE__ */ new WeakMap();
  /**
   * The scratch of the joints. It is not the pool of the contact solver,
   * because `effectiveMass` and `pointVelocity` write into that pool.
   */
  function jointScratch(D) {
    let s = jointScratchCache.get(D);
    if (!s) {
      s = {
        pa: new Float64Array(D.n),
        pb: new Float64Array(D.n),
        d: new Float64Array(D.n),
        ja: new Float64Array(D.n),
        jb: new Float64Array(D.n),
        va: new Float64Array(D.n),
        vb: new Float64Array(D.n),
        err: new Float64Array(D.k),
        dl: new Float64Array(D.k),
        hi: new Float64Array(D.k),
        proj: new Float64Array(Math.max(D.n, D.k)),
        iv: new Float64Array(D.k),
        m1: new Float64Array(D.n * D.n),
        m2: new Float64Array(D.n * D.n),
        uw: new Float64Array(D.n),
        vw: new Float64Array(D.n),
        t1: new Float64Array(D.n),
        t2: new Float64Array(D.n),
        mu: new Float64Array(D.n),
        soft: { biasRate: 0, massScale: 1, impulseScale: 0 }
      };
      jointScratchCache.set(D, s);
    }
    return s;
  }
  /** The static body of the world of each `D`. */
  var worldBodyCache = /* @__PURE__ */ new WeakMap();
  /**
   * The body that a joint to the world uses in the place of `b`. It is static,
   * it is at the origin, and its rotor is the identity. `effectiveMass` gives 0
   * for it, and `applyImpulse` does nothing to it.
   *
   * The body is not in `world.bodies`, thus the broad phase, the contact graph
   * and the sleep never see it. It takes every test for a null body out of the
   * solver.
   */
  function worldBody(D) {
    let b = worldBodyCache.get(D);
    if (!b) {
      b = new Body(D, {
        id: -1,
        name: "world",
        mass: 0,
        shape: { type: "point", n: D.n, boundingRadius: 0 }
      });
      worldBodyCache.set(D, b);
    }
    return b;
  }
  /**
   * The inverse of the inertia that a body shows along the bivector `axis`:
   *
   *   axis . (I'^-1 axis)
   *
   * A static body gives 0. This is the angular form of `effectiveMass`.
   * See ND-PHYSICS.md, A13.
   */
  function angularMass(D, body, axis) {
    if (body.isStatic) return 0;
    const { k } = D;
    const iv = matVec(body.invInertiaWorld, axis, k, k, jointScratch(D).iv);
    let m = 0;
    for (let p = 0; p < k; p += 1) m += iv[p] * axis[p];
    return m;
  }
  /**
   * The rest relation of the two bodies: `Q = Rm_a^T Rm_b`, an `n` x `n`
   * matrix. Take it one time, when the joint starts. The loop is here because
   * `matTranspose` is not in the worker bundle.
   */
  function restRelation(D, a, b) {
    const { n } = D;
    const Q = new Float64Array(n * n);
    for (let i = 0; i < n; i += 1) {
      for (let j = 0; j < n; j += 1) {
        let sum = 0;
        for (let m = 0; m < n; m += 1) sum += a.Rm[m * n + i] * b.Rm[m * n + j];
        Q[i * n + j] = sum;
      }
    }
    return Q;
  }
  /**
   * The error of the rotation of `b` against `a`, as a bivector of `k`
   * components. `Q` is the rest relation from `restRelation`.
   *
   *   T = Rm_a Q      the rotation that `b` must have
   *   M = Rm_b T^T    the error rotation, in the world frame
   *   e[p] = (M[j n + i] - M[i n + j]) / 2   for D.pairs[p] = [i, j]
   *
   * This is the logarithm of `M` to the first order.
   *
   * WHEN IT IS EXACT. For a SIMPLE rotation, thus a rotation in one plane, the
   * result is exactly `sin(angle)` on that plane and 0 on each other plane, at
   * any angle. For a rotation in two planes together the result mixes: a
   * rotation of 1.0 in the plane `(x y)` and of 1.0 in the plane `(y z)` gives
   * 0.354 on the plane `(x z)`, in which nothing turned.
   *
   * Thus a joint must only use this error when the part that stays free is one
   * plane or nothing. `prepareRows` tests that: it uses the error when the
   * joint holds `k` or `k - 1` planes, and it holds the velocity only when the
   * joint leaves two or more planes free. See ND-PHYSICS.md, A13.
   *
   * THE SIGN. The error carries the same sign as `body.w`: a component that is
   * more than zero is a turn from the axis `i` toward the axis `j`. Two
   * measurements give this, and `test/spec/constraint.test.js` holds them:
   *   - `rotorFromPlane(D, i, j, 0.1)` on `b` gives `e[p] = +sin(0.1)`.
   *   - `w[p] = 1` for 0.1 seconds gives `e[p] = +0.09975`.
   * The minus sign of `dR/dt = -(1/2) w R` and the minus sign of the reverse
   * take each other away. See ND-PHYSICS.md, A13.
   *
   * The frame is the world frame, because `invInertiaWorld` and `body.w` are in
   * the world frame.
   */
  function angularError(D, a, b, Q, out) {
    const { n, k } = D;
    const s = jointScratch(D);
    const T = matMul(a.Rm, Q, n, n, n, s.m1);
    const M = matMulT(b.Rm, T, n, n, n, s.m2);
    const e = out || s.err;
    for (let p = 0; p < k; p += 1) {
      const i = D.pairs[p][0];
      const j = D.pairs[p][1];
      e[p] = (M[j * n + i] - M[i * n + j]) / 2;
    }
    return e;
  }
  /** Takes the part of `t` that goes along the unit vector `u` away from `t`. */
  function projectOff(t, dim, u) {
    let d = 0;
    for (let i = 0; i < dim; i += 1) d += t[i] * u[i];
    for (let i = 0; i < dim; i += 1) t[i] -= d * u[i];
  }
  /**
   * Makes `t` orthogonal to each vector of `spans` and of `rest`, then makes it
   * a unit vector. It gives false when `t` becomes shorter than `tol`.
   */
  function orthonormalize(t, dim, spans, rest, tol) {
    for (let c = 0; c < spans.length; c += 1) projectOff(t, dim, spans[c]);
    for (let c = 0; c < rest.length; c += 1) projectOff(t, dim, rest[c]);
    let ln = 0;
    for (let i = 0; i < dim; i += 1) ln += t[i] * t[i];
    ln = Math.sqrt(ln);
    if (ln < tol) return false;
    for (let i = 0; i < dim; i += 1) t[i] /= ln;
    return true;
  }
  /** The largest absolute value that the vectors of `spans` hold on the axis `i`. */
  function axisWeight(spans, i) {
    let m = 0;
    for (let c = 0; c < spans.length; c += 1) {
      const a = Math.abs(spans[c][i]);
      if (a > m) m = a;
    }
    return m;
  }
  /**
   * `want` unit vectors of the length `dim` that are orthogonal to each vector
   * of `spans`. `want` is `dim - spans.length` when you do not give it.
   *
   * The method is the method of `tangentBasis`: try the axes, and start with
   * the axis that `spans` uses the least. `tangentBasis` only works at the
   * length `n`, and this function also works at the length `k`.
   */
  function orthoComplement(dim, spans, tol, want) {
    const need = want !== void 0 ? want : dim - spans.length;
    const T = [];
    const order = [];
    for (let i = 0; i < dim; i += 1) order.push(i);
    order.sort((p, q) => axisWeight(spans, p) - axisWeight(spans, q));
    for (let c = 0; c < order.length; c += 1) {
      if (T.length >= need) break;
      const t = new Float64Array(dim);
      t[order[c]] = 1;
      if (orthonormalize(t, dim, spans, T, tol)) T.push(t);
    }
    return T;
  }
  /**
   * The same as `orthoComplement`, but the result stays near the result of the
   * last step. Give the vectors of the last step in `prev`.
   *
   * THIS IS NECESSARY, and it is not a nicety. A basis that comes from the axes
   * at each step TURNS OVER when `spans` moves past a tie in the order of the
   * axes. The impulse that the row keeps then goes along a different direction,
   * and that puts energy into the world: a jolt one time in each turn of the
   * hinge. Thus the function starts from the vectors of the last step, and it
   * only builds a new vector when an old vector is not usable.
   *
   * The function writes into the vectors of `prev`, thus it makes no garbage
   * after the first step.
   */
  function refreshComplement(dim, spans, prev, tol) {
    const need = dim - spans.length;
    const T = [];
    if (prev) {
      for (let c = 0; c < prev.length; c += 1) {
        if (T.length >= need) break;
        const t = prev[c];
        if (orthonormalize(t, dim, spans, T, tol)) T.push(t);
      }
    }
    if (T.length < need) {
      const more = orthoComplement(dim, spans.concat(T), tol, need - T.length);
      for (let c = 0; c < more.length; c += 1) T.push(more[c]);
    }
    return T;
  }
  /**
   * Makes the error smaller by `slop`, along the direction of the error. It
   * gives a zero vector when the error is not larger than `slop`.
   *
   * The function works on the length of the whole vector, and not on each
   * component. A slop on each component makes a dead zone in the shape of a
   * box, and the joint is then not the same in every direction.
   *
   * GIVE ONLY THE PART OF THE ERROR THAT THE JOINT HOLDS. A direction that the
   * joint leaves free carries an error that grows without limit, and that
   * error would make the length large and the slop would then do nothing. See
   * `projectError`.
   */
  function shrinkError(e, dim, slop) {
    let m = 0;
    for (let i = 0; i < dim; i += 1) m += e[i] * e[i];
    m = Math.sqrt(m);
    if (m <= slop) {
      for (let i = 0; i < dim; i += 1) e[i] = 0;
      return e;
    }
    const f = (m - slop) / m;
    for (let i = 0; i < dim; i += 1) e[i] *= f;
    return e;
  }
  /**
   * The three factors of the softness of a joint, for one substep.
   *
   * A joint with `hertz` at more than 0 is a spring of that frequency and of
   * that damping ratio. The stiffness then does NOT change with the length of
   * the step, which the raw bias of Baumgarte cannot give. The method comes
   * from Erin Catto, "Soft Constraints", GDC 2011.
   *
   *   om = 2 pi hertz
   *   a1 = 2 zeta + dt om;  a2 = dt om a1;  a3 = 1 / (1 + a2)
   *   biasRate = om / a1;  massScale = a2 a3;  impulseScale = a3
   *
   * `hertz` at 0 makes a RIGID joint: the function gives back the old bias
   * rate of `constraintBias / dt`, a mass scale of 1 and an impulse scale of
   * 0. The solver is then the same, number for number.
   *
   * A frequency above `constraintHertzRatio / dt` is not stable, thus the
   * function holds it there. See ND-PHYSICS.md, B9.
   */
  function softness(joint, dt, params, out) {
    const hertz = joint.hertz >= 0 ? joint.hertz : params.constraintHertz;
    if (!(hertz > 0)) {
      out.biasRate = params.constraintBias / dt;
      out.massScale = 1;
      out.impulseScale = 0;
      return out;
    }
    const zeta = joint.damping >= 0 ? joint.damping : params.constraintDamping;
    const limit = params.constraintHertzRatio / dt;
    const hz = hertz > limit ? limit : hertz;
    const om = 2 * Math.PI * hz;
    const a1 = 2 * zeta + dt * om;
    const a2 = dt * om * a1;
    const a3 = 1 / (1 + a2);
    out.biasRate = om / a1;
    out.massScale = a2 * a3;
    out.impulseScale = a3;
    return out;
  }
  /**
   * The error `e` in the frame of the rows: `proj[c] = e . basis[c]`. This
   * takes away every direction that the joint leaves free, thus the slop then
   * works on the part that the joint holds and on nothing else.
   */
  function projectError(D, e, basis, dim, out) {
    for (let c = 0; c < basis.length; c += 1) {
      let sum = 0;
      for (let i = 0; i < dim; i += 1) sum += e[i] * basis[c][i];
      out[c] = sum;
    }
    return out;
  }
  /**
   * Holds the bias between `-limit` and `+limit`. A teleport gives an error of
   * many units, and the bias is the error divided by the time of the step. The
   * limit stops that from throwing the body away.
   */
  function clampBias(bias, limit) {
    if (bias > limit) return limit;
    if (bias < -limit) return -limit;
    return bias;
  }
  /**
   * An orthonormal basis `(u, v)` of the plane of a SIMPLE unit bivector `h`.
   *
   * It builds the `n` x `n` antisymmetric matrix `H` of `h`, takes the column
   * of the largest length as `u`, and takes `v = H u`. For a unit simple
   * bivector `|H u| = |u|` and `v` is orthogonal to `u`, thus the two vectors
   * span the plane.
   *
   * Take the basis ONE TIME, when the joint starts. A basis that came from the
   * columns at each step would jump when the largest column changed, and the
   * angle of the hinge would jump with it.
   *
   * @throws {Error} when the bivector is too short, or when it is not simple
   */
  function planeBasis(D, h, tol) {
    const { n, k } = D;
    const H = new Float64Array(n * n);
    for (let p = 0; p < k; p += 1) {
      const i = D.pairs[p][0];
      const j = D.pairs[p][1];
      H[j * n + i] += h[p];
      H[i * n + j] -= h[p];
    }
    let best = -1;
    let bestLen = 0;
    for (let c = 0; c < n; c += 1) {
      let ln = 0;
      for (let i = 0; i < n; i += 1) ln += H[i * n + c] * H[i * n + c];
      if (ln > bestLen) {
        bestLen = ln;
        best = c;
      }
    }
    if (best < 0 || Math.sqrt(bestLen) < tol) {
      throw new Error("createConstraint: the plane of the joint is too short");
    }
    const u = new Float64Array(n);
    for (let i = 0; i < n; i += 1) u[i] = H[i * n + best];
    let ln = Math.sqrt(bestLen);
    for (let i = 0; i < n; i += 1) u[i] /= ln;
    const v = matVec(H, u, n, n, new Float64Array(n));
    ln = 0;
    for (let i = 0; i < n; i += 1) ln += v[i] * v[i];
    ln = Math.sqrt(ln);
    if (ln < tol) throw new Error("createConstraint: the plane of the joint is too short");
    for (let i = 0; i < n; i += 1) v[i] /= ln;
    return { u, v };
  }
  /**
   * True when the bivector `h` is simple, thus when it is the plane of one
   * rotation. A bivector is simple when `h ^ h = 0`.
   *
   * In 4 dimensions and more a bivector can hold two planes at the same time,
   * for example `e_xy + e_zw`. That is not a hinge, and the angle of a hinge
   * has no meaning for it.
   */
  function isSimpleBivector(D, h) {
    if (D.k < 6) return true;
    const A = mvFromBivector(D, h, mvZero(D));
    const W = mvWedge(D, A, A, mvZero(D));
    let m = 0;
    for (let i = 0; i < W.length; i += 1) m = Math.max(m, Math.abs(W[i]));
    return m < 1e-9;
  }
  /**
   * The angle of a hinge, in radians, in the interval `(-pi, pi]`.
   *
   *   M = Rm_b Q^T Rm_a^T                 the error rotation, in the world
   *   angle = atan2(v . (M u), u . (M u))
   *
   * `(u, v)` is the basis of the free plane, turned into the world by `Rm_a`.
   *
   * THIS IS EXACT at any angle, and it stays exact when the two bodies also
   * turn in other planes. A test holds it at 0.3, 2.0, 3.0 and -1.5 radians
   * with a turn of 0.7 in another plane, to 1e-12. The first order logarithm
   * of `angularError` is NOT exact in that condition; this is, because the
   * free part of a hinge is one plane and a rotation in one plane commutes
   * with itself.
   *
   * The value wraps at `pi`. A limit outside `(-pi, pi)` has no meaning.
   * See ND-PHYSICS.md, B10.
   */
  function hingeAngle(D, joint) {
    const { n } = D;
    const s = jointScratch(D);
    const a = joint.a;
    const b = joint.b;
    const u = matVec(a.Rm, joint.planeU, n, n, s.uw);
    const v = matVec(a.Rm, joint.planeV, n, n, s.vw);
    const t1 = matTVec(a.Rm, u, n, n, s.t1);
    const t2 = matTVec(joint.restRel, t1, n, n, s.t2);
    const mu = matVec(b.Rm, t2, n, n, s.mu);
    let cu = 0;
    let cv = 0;
    for (let i = 0; i < n; i += 1) {
      cu += u[i] * mu[i];
      cv += v[i] * mu[i];
    }
    return Math.atan2(cv, cu);
  }
  /** A new row. `kind` is "linear" or "angular". */
  function makeRow(D, kind) {
    const linear = kind === "linear";
    return {
      kind,
      dir: linear ? new Float64Array(D.n) : null,
      rA: linear ? new Float64Array(D.n) : null,
      rB: linear ? new Float64Array(D.n) : null,
      axis: linear ? null : new Float64Array(D.k),
      km: 0,
      bias: 0,
      lower: -Infinity,
      upper: Infinity,
      impulse: 0,
      massScale: 1,
      impulseScale: 0,
      drift: 0
    };
  }
  /**
   * Makes the impulse of each row zero, and keeps the drift.
   *
   * The warm start uses this. A joint that runs with `useWarmStart` at false
   * must still keep the drift that it has measured, or the correction of the
   * drift goes away without a message.
   */
  function resetImpulses(joint) {
    for (let c = 0; c < joint.rows.length; c += 1) joint.rows[c].impulse = 0;
  }
  /**
   * Makes the impulse AND the drift of each row zero. Use this when the joint
   * loses its history: a teleport of a body, or a new body from `setMass`.
   */
  function resetConstraint(joint) {
    for (let c = 0; c < joint.rows.length; c += 1) {
      joint.rows[c].impulse = 0;
      joint.rows[c].drift = 0;
    }
  }
  /** A list of unit vectors of the length `dim`, from a list of axis numbers. */
  function axisVectors(dim, list) {
    const out = [];
    for (let c = 0; c < list.length; c += 1) {
      const v = new Float64Array(dim);
      if (typeof list[c] === "number") {
        if (list[c] < 0 || list[c] >= dim) {
          throw new Error("createConstraint: the axis " + list[c] + " is not in the range");
        }
        v[list[c]] = 1;
      } else {
        v.set(list[c]);
        let ln = 0;
        for (let i = 0; i < dim; i += 1) ln += v[i] * v[i];
        ln = Math.sqrt(ln);
        if (ln === 0) throw new Error("createConstraint: a direction has the length 0");
        for (let i = 0; i < dim; i += 1) v[i] /= ln;
      }
      out.push(v);
    }
    return out;
  }
  /**
   * Builds a joint. `a` is a body, and `b` is a body or null. A null `b` makes
   * a joint to the world.
   *
   * The anchors `localA` and `localB` are ALWAYS in the local frame of their
   * body. With a null `b`, `localB` takes the world position of the anchor of
   * `a` when you do not give it, thus the joint holds the body where it is.
   *
   * @param {object} def the definition
   * @param {string} def.type point, distance, fixed, hinge or subspace
   * @param {number} [def.id] the id of the joint
   * @param {ArrayLike<number>} [def.localA] the anchor on `a`, `n`
   * @param {ArrayLike<number>} [def.localB] the anchor on `b`, `n`
   * @param {number} [def.rest] the rest length             (distance)
   * @param {string} [def.mode] rod, rope or strut          (distance)
   * @param {number|ArrayLike<number>} [def.plane] the free plane, in the frame
   *   of `a`: a number of `D.pairs`, or a bivector of `k`   (hinge)
   * @param {Array} [def.lockAxes] the directions that the joint holds  (subspace)
   * @param {Array} [def.freeAxes] the directions that stay free        (subspace)
   * @param {Array} [def.lockPlanes] the planes that the joint holds    (subspace)
   * @param {Array} [def.freePlanes] the planes that stay free          (subspace)
   * @param {boolean} [def.worldFrame] true holds the directions in the world
   *   frame. The default holds them in the frame of `a`.       (subspace)
   * @param {boolean} [def.collideConnected] true lets the two bodies touch.
   *   The default is false.
   * @param {boolean} [def.enabled] false leaves the joint out. Default true.
   * @param {object} [params] the params of the world. It gives
   *   `constraintTolerance` to the complement basis.
   */
  function createConstraint(D, def, a, b, params) {
    const { n, k } = D;
    const tol = (params || defaultParams).constraintTolerance;
    const B = b || worldBody(D);
    const type = def.type;
    const joint = {
      id: def.id !== void 0 ? def.id : 0,
      type,
      a,
      b: B,
      toWorld: !b,
      localA: new Float64Array(n),
      localB: new Float64Array(n),
      rest: def.rest !== void 0 ? def.rest : 0,
      mode: def.mode || "rod",
      length: 0,
      restRel: null,
      planeLocal: null,
      linearMode: "none",
      angularMode: "none",
      linearLocal: null,
      angularLocal: null,
      linearBasis: [],
      angularBasis: [],
      useDrift: false,
      planeU: null,
      planeV: null,
      freeAxis: null,
      angle: 0,
      hasLimit: false,
      lowerAngle: def.lowerAngle !== void 0 ? def.lowerAngle : -Infinity,
      upperAngle: def.upperAngle !== void 0 ? def.upperAngle : Infinity,
      motorSpeed: def.motorSpeed !== void 0 ? def.motorSpeed : 0,
      maxMotorTorque: def.maxMotorTorque !== void 0 ? def.maxMotorTorque : 0,
      breakForce: def.breakForce !== void 0 ? def.breakForce : 0,
      breakTorque: def.breakTorque !== void 0 ? def.breakTorque : 0,
      linearImpulse: 0,
      angularImpulse: 0,
      hertz: def.hertz !== void 0 ? def.hertz : -1,
      damping: def.damping !== void 0 ? def.damping : -1,
      worldFrame: def.worldFrame === true,
      collideConnected: def.collideConnected === true,
      enabled: def.enabled !== false,
      rows: []
    };
    if (def.localA) joint.localA.set(def.localA);
    if (def.localB) joint.localB.set(def.localB);
    else if (!b) a.localToWorld(joint.localA, joint.localB);
    if (type === "point") {
      joint.linearMode = "axes";
    } else if (type === "distance") {
      joint.linearMode = "distance";
      if (def.rest === void 0) {
        const pA = a.localToWorld(joint.localA, jointScratch(D).pa);
        const pB = B.localToWorld(joint.localB, jointScratch(D).pb);
        let len = 0;
        for (let i = 0; i < n; i += 1) len += (pB[i] - pA[i]) * (pB[i] - pA[i]);
        joint.rest = Math.sqrt(len);
      }
    } else if (type === "fixed") {
      joint.linearMode = "axes";
      joint.angularMode = "axes";
      joint.restRel = restRelation(D, a, B);
    } else if (type === "hinge") {
      joint.linearMode = "axes";
      joint.angularMode = k > 1 ? "hinge" : "none";
      joint.restRel = restRelation(D, a, B);
      joint.planeLocal = new Float64Array(k);
      const plane = def.plane !== void 0 ? def.plane : 0;
      if (typeof plane === "number") {
        if (plane < 0 || plane >= k) throw new Error("createConstraint: the plane " + plane + " is not in the range");
        joint.planeLocal[plane] = 1;
      } else {
        joint.planeLocal.set(plane);
        let ln = 0;
        for (let p = 0; p < k; p += 1) ln += joint.planeLocal[p] * joint.planeLocal[p];
        ln = Math.sqrt(ln);
        if (ln < tol) throw new Error("createConstraint: the plane of the hinge has the length 0");
        for (let p = 0; p < k; p += 1) joint.planeLocal[p] /= ln;
      }
      if (!isSimpleBivector(D, joint.planeLocal)) {
        throw new Error("createConstraint: the plane of the hinge is not simple");
      }
      const basis = planeBasis(D, joint.planeLocal, tol);
      joint.planeU = basis.u;
      joint.planeV = basis.v;
      joint.hasLimit = joint.lowerAngle > -Infinity || joint.upperAngle < Infinity;
    } else if (type === "subspace") {
      if (def.lockAxes) {
        joint.linearMode = "list";
        joint.linearLocal = axisVectors(n, def.lockAxes);
      } else if (def.freeAxes) {
        joint.linearMode = "list";
        joint.linearLocal = orthoComplement(n, axisVectors(n, def.freeAxes), tol);
      }
      if (def.lockPlanes) {
        joint.angularMode = "list";
        joint.angularLocal = axisVectors(k, def.lockPlanes);
      } else if (def.freePlanes) {
        joint.angularMode = "list";
        joint.angularLocal = orthoComplement(k, axisVectors(k, def.freePlanes), tol);
      }
      if (joint.angularMode !== "none") joint.restRel = restRelation(D, a, B);
    } else {
      throw new Error("createConstraint: the type " + type + " is not known");
    }
    return joint;
  }
  /** Builds the world directions of the linear rows of this step. */
  function updateLinearBasis(D, joint, d, tol) {
    const { n } = D;
    const mode = joint.linearMode;
    if (mode === "none") {
      joint.linearBasis.length = 0;
      return;
    }
    if (mode === "distance") {
      let len = 0;
      for (let i = 0; i < n; i += 1) len += d[i] * d[i];
      len = Math.sqrt(len);
      joint.length = len;
      if (len < tol) {
        joint.linearBasis.length = 0;
        return;
      }
      if (joint.linearBasis.length !== 1) joint.linearBasis = [new Float64Array(n)];
      const u = joint.linearBasis[0];
      for (let i = 0; i < n; i += 1) u[i] = d[i] / len;
      return;
    }
    if (mode === "axes") {
      if (joint.linearBasis.length !== n) {
        joint.linearBasis = [];
        for (let i = 0; i < n; i += 1) {
          const u = new Float64Array(n);
          u[i] = 1;
          joint.linearBasis.push(u);
        }
      }
      return;
    }
    const list = joint.linearLocal;
    if (joint.linearBasis.length !== list.length) {
      joint.linearBasis = [];
      for (let c = 0; c < list.length; c += 1) joint.linearBasis.push(new Float64Array(n));
    }
    for (let c = 0; c < list.length; c += 1) {
      if (joint.worldFrame) joint.linearBasis[c].set(list[c]);
      else matVec(joint.a.Rm, list[c], n, n, joint.linearBasis[c]);
    }
  }
  /** Builds the world bivectors of the angular rows of this step. */
  function updateAngularBasis(D, joint, tol) {
    const { k } = D;
    const mode = joint.angularMode;
    if (mode === "none") {
      joint.angularBasis.length = 0;
      return;
    }
    if (mode === "axes") {
      if (joint.angularBasis.length !== k) {
        joint.angularBasis = [];
        for (let p = 0; p < k; p += 1) {
          const u = new Float64Array(k);
          u[p] = 1;
          joint.angularBasis.push(u);
        }
      }
      return;
    }
    if (mode === "hinge") {
      const h = joint.freeAxis;
      if (h === null) {
        joint.angularBasis.length = 0;
        return;
      }
      joint.angularBasis = refreshComplement(k, [h], joint.angularBasis, tol);
      return;
    }
    const list = joint.angularLocal;
    if (joint.angularBasis.length !== list.length) {
      joint.angularBasis = [];
      for (let c = 0; c < list.length; c += 1) joint.angularBasis.push(new Float64Array(k));
    }
    for (let c = 0; c < list.length; c += 1) {
      if (joint.worldFrame) joint.angularBasis[c].set(list[c]);
      else matVec(joint.a.R2, list[c], k, k, joint.angularBasis[c]);
    }
  }
  /**
   * Builds the rows of a joint. Call this one time in each step, before
   * `solveJoint`.
   *
   * The bias of a row is `-(constraintBias / dt) * error`. It gives back a part
   * of the error in each step. The solver drives the velocity of the row toward
   * the bias, thus the error goes to zero. See ND-PHYSICS.md, B8.
   *
   * The function keeps the rows when their count does not change, because the
   * impulse of a row must live between the steps for the warm start.
   *
   * An angular row only takes a bias when the joint holds `k` or `k - 1`
   * planes. With two or more planes free the first-order logarithm of the
   * relative rotation mixes the planes, and a bias from it would push the body
   * in a plane in which nothing turned. The row then holds the velocity only,
   * and it does not correct a drift. See `angularError`.
   */
  function prepareRows(D, joint, dt, params) {
    const { n, k } = D;
    const s = jointScratch(D);
    const a = joint.a;
    const b = joint.b;
    const tol = params.constraintTolerance;
    const slop = params.constraintSlop;
    const maxBias = params.constraintMaxBias;
    const soft = softness(joint, dt, params, s.soft);
    const gain = soft.biasRate;
    const pA = a.localToWorld(joint.localA, s.pa);
    const pB = b.localToWorld(joint.localB, s.pb);
    const d = s.d;
    for (let i = 0; i < n; i += 1) d[i] = pB[i] - pA[i];
    updateLinearBasis(D, joint, d, tol);
    updateFreeAxis(D, joint, tol);
    updateAngularBasis(D, joint, tol);
    const nl = joint.linearBasis.length;
    const na = joint.angularBasis.length;
    // The row of the motor and the row of the limit are ALWAYS there, and they
    // have `km` at 0 when they are off. A row that came and went would change
    // the count, and a change of the count throws away every impulse that the
    // joint keeps -- at the moment that a limit engages, which is the worst
    // moment. The test "the row identity holds" holds this.
    const extra = joint.planeLocal !== null ? 2 : 0;
    if (joint.rows.length !== nl + na + extra) {
      joint.rows = [];
      for (let c = 0; c < nl; c += 1) joint.rows.push(makeRow(D, "linear"));
      for (let c = 0; c < na + extra; c += 1) joint.rows.push(makeRow(D, "angular"));
    }
    if (!params.useWarmStart) resetImpulses(joint);
    const rawLength = joint.length;
    const proj = s.proj;
    if (joint.linearMode !== "distance") {
      shrinkError(projectError(D, d, joint.linearBasis, n, proj), nl, slop);
    }
    for (let c = 0; c < nl; c += 1) {
      const row = joint.rows[c];
      row.dir.set(joint.linearBasis[c]);
      for (let i = 0; i < n; i += 1) {
        row.rA[i] = pA[i] - a.x[i];
        row.rB[i] = pB[i] - b.x[i];
      }
      row.km = effectiveMass(D, a, row.rA, row.dir, false) + effectiveMass(D, b, row.rB, row.dir, false);
      row.lower = -Infinity;
      row.upper = Infinity;
      row.massScale = soft.massScale;
      row.impulseScale = soft.impulseScale;
      let err = 0;
      if (joint.linearMode === "distance") {
        const raw = rawLength - joint.rest;
        if (joint.mode === "rope") {
          row.lower = -Infinity;
          row.upper = 0;
          if (raw <= 0) {
            row.km = 0;
            row.impulse = 0;
          }
        } else if (joint.mode === "strut") {
          row.lower = 0;
          row.upper = Infinity;
          if (raw >= 0) {
            row.km = 0;
            row.impulse = 0;
          }
        }
        err = raw > slop ? raw - slop : raw < -slop ? raw + slop : 0;
      } else {
        err = proj[c];
      }
      row.bias = clampBias(-gain * err, maxBias);
    }
    const useError = na > 0 && na >= k - 1;
    joint.useDrift = na > 0 && !useError;
    if (na > 0) {
      if (useError) {
        shrinkError(projectError(D, angularError(D, a, b, joint.restRel, s.err),
          joint.angularBasis, k, proj), na, slop);
      } else {
        for (let c = 0; c < na; c += 1) proj[c] = joint.rows[nl + c].drift;
        shrinkError(proj, na, slop);
      }
      for (let c = 0; c < na; c += 1) {
        const row = joint.rows[nl + c];
        row.axis.set(joint.angularBasis[c]);
        row.km = angularMass(D, a, row.axis) + angularMass(D, b, row.axis);
        row.lower = -Infinity;
        row.upper = Infinity;
        row.massScale = soft.massScale;
        row.impulseScale = soft.impulseScale;
        if (useError) row.drift = 0;
        row.bias = clampBias(-gain * proj[c], maxBias);
      }
    }
    if (extra === 0) return;
    prepareHingeRows(D, joint, nl + na, dt, params);
  }
  /**
   * The row of the motor and the row of the limit of a hinge. The two rows sit
   * on the same free plane `h`, thus they work against each other.
   *
   * THE TWO ROWS STAY RIGID. They do not take the softness of the joint: a
   * limit that gives way is not a limit, and a motor is a source of speed and
   * not a spring.
   *
   * The motor comes FIRST and the limit comes SECOND. `solveJoint` walks the
   * rows forward, thus the limit has the last word in each turn. A motor that
   * drives into a limit then stalls, and the joint holds at the limit.
   *
   * The motor takes a TORQUE and not an impulse. The clamp is
   * `maxMotorTorque * dt`, thus the strength of the motor does not change when
   * you change `subSteps`. See ND-PHYSICS.md, B10.
   */
  function prepareHingeRows(D, joint, at, dt, params) {
    const { k } = D;
    const slop = params.constraintSlop;
    const maxBias = params.constraintMaxBias;
    const rigidRate = params.constraintBias / dt;
    const h = joint.freeAxis;
    const km = h === null ? 0 : angularMass(D, joint.a, h) + angularMass(D, joint.b, h);
    joint.angle = h === null ? 0 : hingeAngle(D, joint);
    const motor = joint.rows[at];
    const limit = joint.rows[at + 1];
    motor.massScale = 1;
    motor.impulseScale = 0;
    motor.drift = 0;
    limit.massScale = 1;
    limit.impulseScale = 0;
    limit.drift = 0;
    if (h !== null) {
      motor.axis.set(h);
      limit.axis.set(h);
    }
    if (joint.maxMotorTorque > 0 && km > 0) {
      const cap = joint.maxMotorTorque * dt;
      motor.km = km;
      motor.bias = joint.motorSpeed;
      motor.lower = -cap;
      motor.upper = cap;
    } else {
      motor.km = 0;
      motor.impulse = 0;
      motor.bias = 0;
      motor.lower = 0;
      motor.upper = 0;
    }
    let err = 0;
    let lower = 0;
    let upper = 0;
    if (joint.hasLimit && km > 0) {
      if (joint.angle < joint.lowerAngle) {
        err = joint.angle - joint.lowerAngle;
        lower = 0;
        upper = Infinity;
      } else if (joint.angle > joint.upperAngle) {
        err = joint.angle - joint.upperAngle;
        lower = -Infinity;
        upper = 0;
      }
    }
    if (err !== 0) {
      err = err > slop ? err - slop : err < -slop ? err + slop : 0;
      limit.km = km;
      limit.lower = lower;
      limit.upper = upper;
      limit.bias = clampBias(-rigidRate * err, maxBias);
    } else {
      limit.km = 0;
      limit.impulse = 0;
      limit.bias = 0;
      limit.lower = 0;
      limit.upper = 0;
    }
  }
  /**
   * The free plane of a hinge, in the world frame. A hinge of 2 dimensions has
   * no angular row, thus `updateAngularBasis` never runs for it; the motor and
   * the limit still need the plane.
   */
  function updateFreeAxis(D, joint, tol) {
    const { k } = D;
    if (joint.planeLocal === null) return;
    if (joint.freeAxis === null) joint.freeAxis = new Float64Array(k);
    const h = joint.freeAxis;
    matVec(joint.a.R2, joint.planeLocal, k, k, h);
    let ln = 0;
    for (let p = 0; p < k; p += 1) ln += h[p] * h[p];
    ln = Math.sqrt(ln);
    if (ln < tol) {
      joint.freeAxis = null;
      return;
    }
    for (let p = 0; p < k; p += 1) h[p] /= ln;
  }
  /**
   * Applies the impulse `dJ` of one row. A linear row gives `+dJ dir` to `b` at
   * `rB` and `-dJ dir` to `a` at `rA`, thus the total momentum does not change.
   *
   * The impulse goes through `applyImpulse` and `applyTorqueImpulse`, because
   * `updateDerived` builds `w` again from `L` at the end of each step. A write
   * into `w` would go away.
   */
  function applyRow(D, joint, row, dJ) {
    const { n, k } = D;
    const s = jointScratch(D);
    if (row.kind === "linear") {
      const jb = s.jb;
      const ja = s.ja;
      for (let i = 0; i < n; i += 1) {
        jb[i] = row.dir[i] * dJ;
        ja[i] = -jb[i];
      }
      joint.b.applyImpulse(jb, row.rB, true);
      joint.a.applyImpulse(ja, row.rA, true);
    } else {
      const dl = s.dl;
      for (let p = 0; p < k; p += 1) dl[p] = row.axis[p] * dJ;
      joint.b.applyTorqueImpulse(dl, true);
      for (let p = 0; p < k; p += 1) dl[p] = -dl[p];
      joint.a.applyTorqueImpulse(dl, true);
    }
  }
  /**
   * Applies the impulse of the last step again, before the main loop. A joint
   * then holds with far fewer turns of the loop.
   *
   * A joint keeps its impulses on its own rows, thus it needs no cache. A
   * contact needs `World.applyWarmStartCache`, because the narrow phase builds
   * a new contact in each step.
   */
  function warmStartJoint(D, joint) {
    for (let c = 0; c < joint.rows.length; c += 1) {
      const row = joint.rows[c];
      if (row.km > 0 && row.impulse !== 0) applyRow(D, joint, row, row.impulse);
    }
  }
  /**
   * One turn of the solver on one joint. Call it many times.
   *
   * For each row it finds the velocity of the row, and it applies the impulse
   * that brings that velocity to `bias`. It holds the TOTAL impulse between
   * `lower` and `upper`, and it applies only the change. That clamp on the
   * total, and not on the change, is what makes the method work.
   *
   * A row with `km` at 0 goes away. That is the guard against a division by
   * zero: two static bodies, a joint from a static body to the world, an
   * angular row on a body with no inertia, or a distance joint with the two
   * anchors at the same point. One NaN would poison the whole world in one
   * step.
   */
  /**
   * Adds the violated angular speed of each row to the drift of that row. Call
   * this one time in each substep, AFTER the solver and BEFORE
   * `integratePositions`, thus `w` and `axis` are the values that the solver
   * made.
   *
   * The rotor step of `integratePositions` uses the gyroscopic value of `w`,
   * which is not exactly this one. The difference is of the order of `dt`
   * times the commutator. The drift stays on `body.w`, because that is the
   * value that agrees with the momentum `body.L`.
   *
   * WHY A DRIFT AND NOT A GEOMETRIC ERROR. A lock that leaves two or more
   * planes free has no geometric error that means anything. `SO(n)` is not
   * abelian, thus "the rotation in the plane p" is not defined for a general
   * rotation: a turn of 1.0 in the plane `(x y)` and of 1.0 in the plane
   * `(y z)` reads 0.919 in the plane `(x z)`, and the reading changes with the
   * order of the two turns. No logarithm can take that away. The integral of
   * the violated speed IS defined, and it is what the constraint means.
   *
   * The axis turns during the step, thus the value is correct to the first
   * order only. See ND-PHYSICS.md, B8.
   */
  function accumulateDrift(D, joint, dt) {
    if (!joint.useDrift) return;
    const { k } = D;
    const a = joint.a;
    const b = joint.b;
    for (let c = 0; c < joint.rows.length; c += 1) {
      const row = joint.rows[c];
      if (row.kind !== "angular" || row.km <= 0) continue;
      if (joint.planeLocal !== null && c >= joint.rows.length - 2) continue;
      let cv = 0;
      for (let p = 0; p < k; p += 1) cv += (b.w[p] - a.w[p]) * row.axis[p];
      row.drift += cv * dt;
    }
  }
  /**
   * Measures the load that a joint carries, and says if the joint must break.
   *
   * It gives the length of the impulse of the linear rows and the length of the
   * impulse of the angular rows, and it puts the two values on the joint.
   *
   * IT COMPARES A FORCE AND A TORQUE, and not an impulse. `row.impulse` is the
   * impulse of ONE substep, thus a limit on it would mean something different
   * when you change `subSteps`. `impulse / dt` is a force, and it does not
   * change.
   *
   * The two values stay apart because an impulse and an angular impulse do not
   * have the same units, and one length of the two together has no meaning. A
   * lock of the rotation has no linear row, and it needs `breakTorque`.
   *
   * THE ROW OF THE MOTOR DOES NOT COUNT. A motor that pushes against its own
   * limit carries a large impulse by design, and it would break its own joint.
   * The row of the limit does count, because it is a true load.
   *
   * A value of 0 means "never break", thus no `Infinity` goes in a message.
   */
  function jointStress(D, joint, dt) {
    const rows = joint.rows;
    const motorAt = joint.planeLocal !== null ? rows.length - 2 : -1;
    let lin = 0;
    let ang = 0;
    for (let c = 0; c < rows.length; c += 1) {
      if (c === motorAt) continue;
      const j = rows[c].impulse;
      if (rows[c].kind === "linear") lin += j * j;
      else ang += j * j;
    }
    joint.linearImpulse = Math.sqrt(lin);
    joint.angularImpulse = Math.sqrt(ang);
    return joint.breakForce > 0 && joint.linearImpulse > joint.breakForce * dt || joint.breakTorque > 0 && joint.angularImpulse > joint.breakTorque * dt;
  }
  function solveJoint(D, joint) {
    const { n, k } = D;
    const s = jointScratch(D);
    const a = joint.a;
    const b = joint.b;
    for (let c = 0; c < joint.rows.length; c += 1) {
      const row = joint.rows[c];
      if (row.km <= 0) continue;
      let cv = 0;
      if (row.kind === "linear") {
        const ua = a.pointVelocity(row.rA, s.va);
        const ub = b.pointVelocity(row.rB, s.vb);
        for (let i = 0; i < n; i += 1) cv += (ub[i] - ua[i]) * row.dir[i];
      } else {
        for (let p = 0; p < k; p += 1) cv += (b.w[p] - a.w[p]) * row.axis[p];
      }
      let dJ = -row.massScale * (cv - row.bias) / row.km - row.impulseScale * row.impulse;
      const old = row.impulse;
      let total = old + dJ;
      if (total < row.lower) total = row.lower;
      if (total > row.upper) total = row.upper;
      row.impulse = total;
      dJ = total - old;
      if (dJ !== 0) applyRow(D, joint, row, dJ);
    }
  }
  // src/nd/world.js
  /**
   * The tolerances of the world. Give your own values in
   * `new World({ params: { ... } })`, or in `scene.setSolverParams({ ... })`.
   * See ND-PHYSICS.md, B7, and README, section 12.
   */
  var defaultParams = {
    // time
    /** The length of a step in seconds, when the caller does not give one. */
    fixedTimeStep: 1 / 60,
    /** The count of the parts of one step. More parts give more accuracy. */
    subSteps: 1,
    // solver
    /** The turns of the main solver loop. More turns make a stack more stable. */
    iterations: 10,
    /** The turns of the solver in the shock propagation pass. */
    shockIterations: 2,
    /** True runs the shock propagation pass. Necessary for a tall stack. */
    useShockPropagation: true,
    /** True applies the impulse of the last step again. It makes a stack rest. */
    useWarmStart: true,
    // contact
    /** A depth that the solver accepts. It stops a contact from shaking. */
    penetrationSlop: 5e-3,
    /**
     * How much of the depth the solver corrects in one step, 0 to 1. A large
     * value pushes the bodies apart fast, but it can add energy.
     */
    biasFactor: 0.2,
    /** The extra length of the boxes of the broad phase, and of the contacts. */
    contactMargin: 0.02,
    /** Below this speed a contact does not bounce. It stops small bounces. */
    restitutionThreshold: 0.5,
    /** The most contacts of one pair. `World` changes 0 into `2^(n-1)`. */
    maxContacts: 0,
    // zero means 2^(n-1)
    // constraint
    /**
     * How much of the error of a joint the solver corrects in one step, 0 to 1.
     * A large value closes a joint fast, but it can add energy. This is not
     * `biasFactor`: that value is for the push-out of a contact.
     */
    constraintBias: 0.2,
    /** An error of a joint that the solver accepts. It stops a joint from shaking. */
    constraintSlop: 1e-3,
    /**
     * The largest speed that the bias of a joint asks for. Without this limit a
     * body that a teleport moves a long way gets a very large velocity, because
     * the bias is the error divided by the time of the step.
     */
    constraintMaxBias: 10,
    /**
     * Below this length a direction of a joint is not defined, and the solver
     * leaves that row out. The value is also the shortest vector that the
     * complement basis of a hinge accepts.
     */
    constraintTolerance: 1e-6,
    /**
     * The frequency of the softness of a joint, in Hz. 0 makes a rigid joint,
     * and the joint then uses `constraintBias`. A small value makes a soft
     * joint that can stretch, and the stiffness does not change with the
     * length of the step. See ND-PHYSICS.md, B9.
     */
    constraintHertz: 0,
    /** The damping ratio of the softness. 1 gives no overshoot. */
    constraintDamping: 1,
    /**
     * The largest part of the rate of a substep that the softness may use. A
     * frequency above this is not stable, and the solver holds it here.
     */
    constraintHertzRatio: 0.25,
    // rotor
    /**
     * The error of the rotor that starts a repair. See `rotorCorrect`. A large
     * value is faster, but the body then changes its shape as it turns.
     */
    rotorTolerance: 1e-9,
    // gyroscopic term
    /** True applies the gyroscopic term. Necessary for a body that tumbles. */
    gyroscopic: true,
    /** The turns of the Newton loop of the gyroscopic term. */
    gyroscopicIterations: 1,
    // sleep
    /** True lets a body that is almost still go to sleep. */
    allowSleep: true,
    /** Below this linear speed a body can go to sleep. */
    sleepLinearVelocity: 0.03,
    /** Below this angular speed a body can go to sleep. */
    sleepAngularVelocity: 0.03,
    /** The time in seconds that a body must be still before it sleeps. */
    sleepTime: 0.6
  };
  /**
   * The world of the physics. It holds the bodies, and it moves them.
   *
   * This class does not need three.js. Use it directly for a test, or for a
   * program with no display. See README, section 11.
   *
   *   const world = new World({ dimensions: 4 });
   *   world.createBody({ shape: HyperSphere(world.D, 1), position: [0, 5, 0, 0] });
   *   world.step(1 / 60);
   *
   * The event `collision` gives `(a, b, contact)` for each contact with a
   * depth of more than zero.
   */
  var World = class {
    /**
     * @param {object} [opts] the options
     * @param {number} [opts.dimensions] the count of dimensions. Default 3.
     * @param {ArrayLike<number>} [opts.gravity] the gravity, of length `n`.
     *   The default is -9.81 on the axis 1.
     * @param {object} [opts.params] your own values. See `defaultParams`.
     */
    constructor(opts = {}) {
      const n = opts.dimensions || 3;
      this.D = dims(n);
      this.n = n;
      this.params = Object.assign({}, defaultParams, opts.params || {});
      if (!this.params.maxContacts) this.params.maxContacts = 1 << n - 1;
      this.gravity = new Float64Array(n);
      if (opts.gravity) this.gravity.set(opts.gravity);
      else this.gravity[1] = -9.81;
      this.bodies = [];
      this.contacts = [];
      this.constraints = [];
      this.manifolds = /* @__PURE__ */ new Map();
      this.time = 0;
      this.listeners = { collision: [], constraintBroken: [] };
      this._noCollide = /* @__PURE__ */ new Set();
      this._islandRoot = new Int32Array(0);
      this._islandTimer = new Float64Array(0);
      this._islandOk = new Uint8Array(0);
    }
    /** Adds a body that you built. It gives that body. */
    addBody(body) {
      this.bodies.push(body);
      return body;
    }
    /**
     * Takes a body out of the world. It does nothing when the body is not in it.
     *
     * It also takes out each constraint that holds the body. A constraint that
     * stays would pull a body that nothing integrates, and the momentum would
     * come from nothing.
     */
    removeBody(body) {
      const i = this.bodies.indexOf(body);
      if (i >= 0) this.bodies.splice(i, 1);
      let dropped = false;
      for (let c = this.constraints.length - 1; c >= 0; c -= 1) {
        const j = this.constraints[c];
        if (j.a === body || j.b === body) {
          this.constraints.splice(c, 1);
          dropped = true;
        }
      }
      if (dropped) this._rebuildNoCollide();
    }
    /** Makes a body with the options of `Body`, and adds it. It gives the body. */
    createBody(opts) {
      return this.addBody(new Body(this.D, opts));
    }
    /**
     * Adds a constraint that `createConstraint` built. It gives that constraint.
     * See ND-PHYSICS.md, A13.
     */
    addConstraint(joint) {
      this.constraints.push(joint);
      if (!joint.collideConnected) this._rebuildNoCollide();
      joint.a.wake();
      joint.b.wake();
      return joint;
    }
    /** Makes a constraint with the options of `createConstraint`, and adds it. */
    createConstraint(def, a, b) {
      return this.addConstraint(createConstraint(this.D, def, a, b, this.params));
    }
    /** Takes a constraint out of the world. */
    removeConstraint(joint) {
      const i = this.constraints.indexOf(joint);
      if (i < 0) return;
      this.constraints.splice(i, 1);
      this._rebuildNoCollide();
    }
    /**
     * Puts `next` in the place of `old`, in the bodies and in the constraints.
     *
     * `setMass` builds a new `Body`, because the inertia comes from the mass at
     * build time. A constraint that keeps the old body would then pull a body
     * that is not in the world.
     */
    replaceBody(old, next) {
      const i = this.bodies.indexOf(old);
      if (i >= 0) this.bodies[i] = next;
      else this.bodies.push(next);
      for (const j of this.constraints) {
        if (j.a === old) j.a = next;
        if (j.b === old) j.b = next;
        if (j.a === next || j.b === next) resetConstraint(j);
      }
      this._rebuildNoCollide();
    }
    /**
     * The constraints that the solver must work on in this step. It leaves out
     * a constraint that is off, and a constraint whose bodies are all static or
     * asleep. It wakes the partner of a body that is awake.
     */
    activeConstraints() {
      const out = [];
      for (const j of this.constraints) {
        if (!j.enabled) continue;
        const liveA = !j.a.isStatic && !j.a.sleeping;
        const liveB = !j.b.isStatic && !j.b.sleeping;
        if (!liveA && !liveB) continue;
        if (j.a.sleeping) j.a.wake();
        if (j.b.sleeping) j.b.wake();
        out.push(j);
      }
      return out;
    }
    /** Builds the set of the pairs of ids that must not touch. */
    _rebuildNoCollide() {
      this._noCollide.clear();
      for (const j of this.constraints) {
        if (j.collideConnected) continue;
        const lo = Math.min(j.a.id, j.b.id);
        const hi = Math.max(j.a.id, j.b.id);
        this._noCollide.add(`${lo}:${hi}`);
      }
    }
    /** Sets the gravity, of length `n`, and wakes all of the bodies. */
    setGravity(g) {
      this.gravity.set(g);
      for (const b of this.bodies) b.wake();
    }
    /**
     * Adds a listener. The world sends `collision`, `stepStart` and
     * `constraintBroken`. `constraintBroken` gives the joint, and the world has
     * already taken that joint out when the event goes out.
     */
    on(name, fn) {
      (this.listeners[name] = this.listeners[name] || []).push(fn);
    }
    /** Sends an event to each listener of that name. */
    emit(name, ...args) {
      for (const fn of this.listeners[name] || []) fn(...args);
    }
    /**
     * The broad phase: it finds the pairs of bodies that can touch.
     *
     * The method is sort and sweep on the axis 0. It puts the boxes in the
     * order of their lowest point on that axis. For each box it then looks
     * only at the boxes that start before the end of that box, and it stops
     * the inner loop at the first box that starts after it.
     *
     * It drops a pair of two static bodies, and a pair in which no body is
     * awake. It also drops a pair that a constraint holds together, because
     * `collideConnected` is false by default.
     *
     * @returns {Body[][]} the pairs whose boxes overlap
     */
    broadPhase() {
      const { n } = this;
      const margin = this.params.contactMargin;
      const noCollide = this._noCollide.size > 0 ? this._noCollide : null;
      const items = [];
      for (const b of this.bodies) {
        const box = b.aabb(margin);
        items.push({ b, min: box.min, max: box.max });
      }
      items.sort((p, q) => p.min[0] - q.min[0]);
      const pairs = [];
      for (let i = 0; i < items.length; i += 1) {
        for (let j = i + 1; j < items.length; j += 1) {
          if (items[j].min[0] > items[i].max[0]) break;
          const A = items[i];
          const B = items[j];
          if (A.b.isStatic && B.b.isStatic) continue;
          if (A.b.sleeping && B.b.sleeping) continue;
          if (A.b.sleeping && B.b.isStatic || B.b.sleeping && A.b.isStatic) continue;
          if (noCollide) {
            const lo = Math.min(A.b.id, B.b.id);
            const hi = Math.max(A.b.id, B.b.id);
            if (noCollide.has(`${lo}:${hi}`)) continue;
          }
          let hit = true;
          for (let t = 1; t < n; t += 1) {
            if (A.min[t] > B.max[t] || B.min[t] > A.max[t]) {
              hit = false;
              break;
            }
          }
          if (hit) pairs.push([A.b, B.b]);
        }
      }
      return pairs;
    }
    /** Calls `collide` on each pair. It gives the list of the contacts. */
    narrowPhase(pairs) {
      const out = [];
      for (const [a, b] of pairs) collide(this.D, a, b, out, this.params);
      return out;
    }
    /**
     * Copies the impulses of the last step onto the new contacts. The narrow
     * phase makes new contact objects at each step, thus the solver would
     * start at zero without this.
     *
     * The key of a manifold is the pair of the ids of the two bodies. Inside a
     * manifold the function looks for the old contact that is nearest to the
     * new point, within 0.01 of length. It then keeps the new manifolds for
     * the next step.
     */
    applyWarmStartCache(contacts) {
      const { n } = this;
      const next = /* @__PURE__ */ new Map();
      for (const c of contacts) {
        const key = `${c.a.id}:${c.b.id}`;
        const old = this.manifolds.get(key);
        if (old) {
          let best = null;
          let bestD = 0.01 * 0.01;
          for (const o of old) {
            let d = 0;
            for (let i = 0; i < n; i += 1) {
              const t = o.point[i] - c.point[i];
              d += t * t;
            }
            if (d < bestD) {
              bestD = d;
              best = o;
            }
          }
          if (best) {
            c.normalImpulse = best.normalImpulse;
            c.oldTangent = best.tangentImpulse;
          }
        }
        if (!next.has(key)) next.set(key, []);
        next.get(key).push(c);
      }
      this.manifolds = next;
    }
    /**
     * Moves the world ahead by `dt` seconds. It cuts `dt` into
     * `params.subSteps` parts, and it makes the forces zero at the end.
     * @param {number} dt the time in seconds
     */
    step(dt) {
      const { D, params } = this;
      const sub = Math.max(1, params.subSteps | 0);
      const h = dt / sub;
      for (let s = 0; s < sub; s += 1) this.subStep(h);
      for (const b of this.bodies) b.clearForces();
    }
    /**
     * One part of a step. The order of the work does not change:
     *
     *   1. The forces change the velocities.       `integrateVelocities`
     *   2. Find the pairs that can touch.          `broadPhase`
     *   3. Find the contacts.                      `narrowPhase`
     *   4. Wake the bodies, and send `collision`.
     *   5. Copy the impulses of the last step.     `applyWarmStartCache`
     *   6. Build the rows of the joints.           `prepareRows`
     *   7. Build the contacts for the solver.      `prepareContact`
     *   8. Apply the old impulses again.           `warmStartJoint`, `warmStart`
     *   9. The main solver loop.            `solveJoint`, `solveContact`
     *  10. The extra pass for a stack.             `shockPropagation`
     *  11. One more pass on the joints.            `solveJoint`
     *  12. Add up the drift of the joints.          `accumulateDrift`
     *  13. Break the joints that carry too much.    `jointStress`
     *  14. The velocities change the positions.    `integratePositions`
     *  15. Sleep.
     *
     * The main loop changes its direction at each turn. That takes away the
     * effect of the order of the contacts, and a stack then rests level.
     *
     * The joints go before the contacts on an even turn, and after them on an
     * odd turn. The odd turn is then a true reversal of the even turn. A block
     * of joints that always stays at the same end would leave the same bias
     * that the change of direction takes away. And a contact is an inequality,
     * thus its clamp runs last on the even turn: a body that a joint pushes
     * into the ground goes back out there.
     *
     * Step 11 is necessary because `shockPropagation` makes the impulses of the
     * contacts zero and solves them again with false masses. That pass can pull
     * a joint open, after the main loop has closed it. The extra pass keeps the
     * impulses of the joints, thus it costs almost nothing.
     *
     * With no joint the work is the same as it was before the joints existed.
     */
    subStep(dt) {
      const { D, params } = this;
      for (const b of this.bodies) integrateVelocities(D, b, dt, this.gravity, params);
      const pairs = this.broadPhase();
      const contacts = this.narrowPhase(pairs);
      this.contacts = contacts;
      for (const c of contacts) {
        if (c.depth > 0) {
          if (c.a.sleeping && !c.b.isStatic && !c.b.sleeping) c.a.wake();
          if (c.b.sleeping && !c.a.isStatic && !c.a.sleeping) c.b.wake();
          this.emit("collision", c.a, c.b, c);
        }
      }
      const joints = this.constraints.length > 0 ? this.activeConstraints() : this.constraints;
      if (params.useWarmStart) this.applyWarmStartCache(contacts);
      for (let i = 0; i < joints.length; i += 1) prepareRows(D, joints[i], dt, params);
      for (const c of contacts) prepareContact(D, c, dt, params);
      if (params.useWarmStart) {
        for (let i = 0; i < joints.length; i += 1) warmStartJoint(D, joints[i]);
        for (const c of contacts) {
          if (c.oldTangent && c.oldTangent.length === c.tangentImpulse.length) {
            c.tangentImpulse.set(c.oldTangent);
          }
          warmStart(D, c);
        }
      }
      for (let it = 0; it < params.iterations; it += 1) {
        if ((it & 1) === 0) {
          for (let i = 0; i < joints.length; i += 1) solveJoint(D, joints[i]);
          for (let i = 0; i < contacts.length; i += 1) solveContact(D, contacts[i]);
        } else {
          for (let i = contacts.length - 1; i >= 0; i -= 1) solveContact(D, contacts[i]);
          for (let i = joints.length - 1; i >= 0; i -= 1) solveJoint(D, joints[i]);
        }
      }
      if (params.useShockPropagation && contacts.length > 0) {
        buildContactGraph(this.bodies, contacts);
        shockPropagation(D, contacts, dt, params);
        for (let i = 0; i < joints.length; i += 1) solveJoint(D, joints[i]);
      }
      for (let i = 0; i < joints.length; i += 1) accumulateDrift(D, joints[i], dt);
      let broken = null;
      for (let i = 0; i < joints.length; i += 1) {
        if (jointStress(D, joints[i], dt)) {
          if (broken === null) broken = [];
          broken.push(joints[i]);
        }
      }
      if (broken !== null) {
        for (let i = 0; i < broken.length; i += 1) {
          const j = broken[i];
          this.removeConstraint(j);
          j.a.wake();
          j.b.wake();
          this.emit("constraintBroken", j);
        }
      }
      for (const b of this.bodies) integratePositions(D, b, dt, params);
      if (params.allowSleep) this.updateSleep(dt, contacts);
      this.time += dt;
    }
    /**
     * Puts a body to sleep after it is almost still for `params.sleepTime`
     * seconds. A sleeping body does not move, and the broad phase drops it.
     * A contact with an awake body, or any force, wakes it again.
     *
     * THE BODIES SLEEP IN ISLANDS. A contact and a joint each join two bodies,
     * and a whole island sleeps together or none of it sleeps. Without that
     * rule one box of a stack could sleep while the box under it still moves,
     * and the stack would come apart. A joint is an edge of the island too,
     * thus two bodies that a joint holds always sleep at the same time, even
     * when no contact joins them.
     *
     * A static body is in no island. The ground touches everything, and it
     * would make one island of the whole world.
     */
    updateSleep(dt, contacts) {
      const { n, k } = this.D;
      const bodies = this.bodies;
      const count = bodies.length;
      if (this._islandRoot.length < count) {
        this._islandRoot = new Int32Array(count);
        this._islandTimer = new Float64Array(count);
        this._islandOk = new Uint8Array(count);
      }
      const root = this._islandRoot;
      const timer = this._islandTimer;
      const ok = this._islandOk;
      for (let i = 0; i < count; i += 1) {
        bodies[i]._island = i;
        root[i] = i;
      }
      const find = (x) => {
        let r = x;
        while (root[r] !== r) {
          root[r] = root[root[r]];
          r = root[r];
        }
        return r;
      };
      const union = (x, y) => {
        const rx = find(x);
        const ry = find(y);
        if (rx !== ry) root[rx] = ry;
      };
      // A contact and a joint each join two bodies into one island. A static
      // body joins nothing: the ground would then make one island of the whole
      // world, and nothing would ever sleep.
      if (contacts) {
        for (let c = 0; c < contacts.length; c += 1) {
          const p = contacts[c];
          if (!p.a.isStatic && !p.b.isStatic) union(p.a._island, p.b._island);
        }
      }
      for (let c = 0; c < this.constraints.length; c += 1) {
        const j = this.constraints[c];
        if (!j.enabled) continue;
        if (!j.a.isStatic && !j.b.isStatic) union(j.a._island, j.b._island);
      }
      // The timer of each body, then the smallest timer of each island.
      const lv = this.params.sleepLinearVelocity ** 2;
      const av = this.params.sleepAngularVelocity ** 2;
      for (let i = 0; i < count; i += 1) {
        timer[i] = Infinity;
        ok[i] = 1;
      }
      for (let i = 0; i < count; i += 1) {
        const b = bodies[i];
        if (b.isStatic) continue;
        let sv = 0;
        for (let x = 0; x < n; x += 1) sv += b.v[x] * b.v[x];
        let sw = 0;
        for (let p = 0; p < k; p += 1) sw += b.w[p] * b.w[p];
        b.sleepTimer = sv < lv && sw < av ? b.sleepTimer + dt : 0;
        const r = find(i);
        if (b.sleepTimer < timer[r]) timer[r] = b.sleepTimer;
        if (!b.allowSleep) ok[r] = 0;
      }
      // A whole island sleeps together, or none of it sleeps.
      const need = this.params.sleepTime;
      for (let i = 0; i < count; i += 1) {
        const b = bodies[i];
        if (b.isStatic) continue;
        const r = find(i);
        const rest = ok[r] === 1 && timer[r] > need;
        if (rest && !b.sleeping) {
          b.sleeping = true;
          b.v.fill(0);
          b.L.fill(0);
          b.w.fill(0);
        } else if (!rest && b.sleeping) {
          b.sleeping = false;
        }
      }
    }
    /**
     * The kinetic energy plus the potential energy of all of the bodies that
     * are not static. This is a tool for a test: with no contact and no
     * damping, the value must stay the same. See ND-PHYSICS.md, E.
     */
    totalEnergy() {
      let e = 0;
      for (const b of this.bodies) {
        if (b.isStatic) continue;
        e += b.kineticEnergy();
        for (let i = 0; i < this.n; i += 1) e -= b.mass * this.gravity[i] * b.x[i];
      }
      return e;
    }
  };

  // src/slice/slice4.js
  //
  // The display of 4 dimensions. A screen shows 3 dimensions, thus the plugin
  // cuts the world at one value of the fourth axis `w`, and it draws the cut.
  // This is the same as a 3D program that cuts a body with a plane and draws
  // the section. See ND-PHYSICS.md, D1.
  //
  // The cut of a 4D body is a normal 3D surface. Thus three.js draws it with
  // no change. `PhysiN.Scene.sliceW` holds the value of the cut, and
  // `refreshSlice()` draws it again after a change.
  //
  // A body that does not meet the cut is not visible.
  var slice4_exports = {};
  __export(slice4_exports, {
    sliceHyperSphereRadius: () => sliceHyperSphereRadius,
    sliceTetrahedra: () => sliceTetrahedra,
    sliceTorusMinorRadius: () => sliceTorusMinorRadius,
    torusSliceTriangles: () => torusSliceTriangles,
    transformVertices: () => transformVertices
  });
  /**
   * Changes the vertices of a 4D mesh from the body frame to the world frame,
   * with `Rm x + position`. It uses the buffer `out` again when it is of the
   * correct length, thus a frame makes no garbage.
   *
   * @param {Float64Array} vertices the vertices, 4 numbers for each
   * @param {Float64Array} Rm the rotation matrix, 4 by 4
   * @param {Float64Array} x the position, of length 4
   */
  function transformVertices(vertices, Rm, x, out) {
    const count = vertices.length / 4;
    const V = out && out.length === vertices.length ? out : new Float64Array(vertices.length);
    for (let m = 0; m < count; m += 1) {
      const b = m * 4;
      for (let i = 0; i < 4; i += 1) {
        let s = x[i];
        for (let j = 0; j < 4; j += 1) s += Rm[i * 4 + j] * vertices[b + j];
        V[b + i] = s;
      }
    }
    return V;
  }
  /** Four 3D points that `sliceTetrahedra` uses again at each cut. */
  var POINT_POOL = [
    new Float64Array(3),
    new Float64Array(3),
    new Float64Array(3),
    new Float64Array(3)
  ];
  /**
   * Cuts a 4D mesh at `w = value`, and gives the 3D triangles of the cut. This
   * is the marching tetrahedra method, one dimension higher.
   *
   * For each tetrahedron the function looks at the sign of `w - value` at its
   * four vertices. A tetrahedron with all four on one side gives nothing. In
   * all other conditions the cut goes through the edges that change sign, and
   * there are only two forms:
   *   - 1 vertex on one side, 3 on the other: 3 cut points, thus 1 triangle.
   *   - 2 vertices on each side: 4 cut points, thus 2 triangles.
   * The point on an edge comes from a straight interpolation on `w`.
   *
   * @param {Float64Array} worldVertices the vertices in the world frame, 4 for each
   * @param {ArrayLike<number>} cells the tetrahedra, 4 indices for each
   * @param {number} value the value of `w` of the cut
   * @param {{positions: number[], normals: number[]}} [sink] the arrays to add to
   * @returns {{positions: number[], normals: number[]}} the triangles, 9
   *   numbers for each, and a normal for each vertex
   */
  function sliceTetrahedra(worldVertices, cells, value, sink) {
    const cellCount = cells.length / 4;
    const out = sink || { positions: [], normals: [] };
    const pos = out.positions;
    const nrm = out.normals;
    const d = new Float64Array(4);
    const idx = new Int32Array(4);
    const above = [];
    const below = [];
    const pts = [];
    for (let c = 0; c < cellCount; c += 1) {
      above.length = 0;
      below.length = 0;
      for (let t = 0; t < 4; t += 1) {
        idx[t] = cells[c * 4 + t];
        d[t] = worldVertices[idx[t] * 4 + 3] - value;
        if (d[t] >= 0) above.push(t);
        else below.push(t);
      }
      if (above.length === 0 || below.length === 0) continue;
      pts.length = 0;
      const cut = (t0, t1) => {
        const a = idx[t0] * 4;
        const b = idx[t1] * 4;
        const f = d[t0] / (d[t0] - d[t1]);
        const p = POINT_POOL[pts.length];
        for (let i = 0; i < 3; i += 1) {
          p[i] = worldVertices[a + i] + f * (worldVertices[b + i] - worldVertices[a + i]);
        }
        pts.push(p);
      };
      if (above.length === 1 || below.length === 1) {
        const lone = above.length === 1 ? above[0] : below[0];
        const other = above.length === 1 ? below : above;
        cut(lone, other[0]);
        cut(lone, other[1]);
        cut(lone, other[2]);
        emitTriangle(pos, nrm, pts[0], pts[1], pts[2]);
      } else {
        const [a0, a1] = above;
        const [b0, b1] = below;
        cut(a0, b0);
        cut(a0, b1);
        cut(a1, b1);
        cut(a1, b0);
        emitTriangle(pos, nrm, pts[0], pts[1], pts[2]);
        emitTriangle(pos, nrm, pts[0], pts[2], pts[3]);
      }
    }
    return out;
  }
  /**
   * Adds one triangle, with its normal at each of the three vertices. The
   * normal is the cross product of two edges. It drops a triangle that is
   * almost a line, because such a triangle has no normal.
   */
  function emitTriangle(pos, nrm, p0, p1, p2) {
    const ux = p1[0] - p0[0];
    const uy = p1[1] - p0[1];
    const uz = p1[2] - p0[2];
    const vx = p2[0] - p0[0];
    const vy = p2[1] - p0[1];
    const vz = p2[2] - p0[2];
    let nx = uy * vz - uz * vy;
    let ny = uz * vx - ux * vz;
    let nz = ux * vy - uy * vx;
    const ln = Math.hypot(nx, ny, nz);
    if (ln < 1e-12) return;
    nx /= ln;
    ny /= ln;
    nz /= ln;
    pos.push(p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], p2[0], p2[1], p2[2]);
    nrm.push(nx, ny, nz, nx, ny, nz, nx, ny, nz);
  }
  /**
   * The radius of the cut of a hypersphere. The cut of a 4-ball is a 3-ball,
   * thus a mesh is not necessary. The radius is `sqrt(R^2 - d^2)`, with `d`
   * the distance from the center to the cut.
   * @returns {number} the radius, or -1 when the ball does not meet the cut
   */
  function sliceHyperSphereRadius(radius, centerW, value) {
    const d = value - centerW;
    const s = radius * radius - d * d;
    return s > 0 ? Math.sqrt(s) : -1;
  }
  /**
   * The minor radius of the cut of a 4D torus. The tube of a 4D torus is a
   * 2-sphere, thus its cut is a circle of the radius `sqrt(r^2 - d^2)`. The
   * cut of the whole body is then a normal 3D torus.
   *
   * This is correct only when the torus does not turn out of the cut.
   * `torusFastPath` tests that.
   *
   * @returns {number} the radius, or -1 when the torus does not meet the cut
   */
  function sliceTorusMinorRadius(minorRadius, centerW, value) {
    const d = value - centerW;
    const s = minorRadius * minorRadius - d * d;
    return s > 0 ? Math.sqrt(s) : -1;
  }
  /**
   * The triangles of the cut of a 4D torus, without a mesh. This is the fast
   * path: it makes a smooth 3D torus of the minor radius `minorRadius`
   * directly, and it needs no tetrahedra.
   *
   * It writes into `Float32Array` buffers that it holds in `out`, and it uses
   * them again at each frame. Thus a frame makes no garbage.
   *
   * @param {object} out the buffers of the last call, or `{}`
   * @param {Shape} shape the torus
   * @param {Float64Array} M the rotation matrix of the body, 4 by 4
   * @param {Float64Array} x the position of the body, of length 4
   * @param {number} deltaW the distance from the body to the cut, on `w`
   * @param {number} minorRadius the radius from `sliceTorusMinorRadius`
   * @param {number[]} [segments] the divisions. Default `[40, 20]`.
   */
  function torusSliceTriangles(out, shape, M, x, deltaW, minorRadius, segments) {
    const [i, j] = shape.plane;
    const p = shape.otherAxes[0] === 3 ? shape.otherAxes[1] : shape.otherAxes[0];
    const R = shape.majorRadius;
    const rho = minorRadius;
    const nt = segments && segments[0] || 40;
    const nb = segments && segments[1] || 20;
    const tris = nt * nb * 2;
    const need = tris * 9;
    if (!out.buf || out.buf.length < need) {
      out.buf = new Float32Array(need);
      out.nbuf = new Float32Array(need);
    }
    const P = out.buf;
    const N = out.nbuf;
    const L = new Float64Array(4);
    const NL = new Float64Array(4);
    const wp = new Float64Array(3);
    const wn = new Float64Array(3);
    const point = (t, b, dst, dstN, at) => {
      const th = 2 * Math.PI * t / nt;
      const be = 2 * Math.PI * b / nb;
      const cb = Math.cos(be);
      const sb = Math.sin(be);
      const ct = Math.cos(th);
      const st = Math.sin(th);
      const rad = R + rho * cb;
      L[0] = 0;
      L[1] = 0;
      L[2] = 0;
      L[3] = deltaW;
      L[i] = rad * ct;
      L[j] = rad * st;
      L[p] = rho * sb;
      NL[0] = 0;
      NL[1] = 0;
      NL[2] = 0;
      NL[3] = 0;
      NL[i] = cb * ct;
      NL[j] = cb * st;
      NL[p] = sb;
      for (let a = 0; a < 3; a += 1) {
        let s = x[a];
        let sn = 0;
        for (let c = 0; c < 4; c += 1) {
          s += M[a * 4 + c] * L[c];
          sn += M[a * 4 + c] * NL[c];
        }
        dst[at + a] = s;
        dstN[at + a] = sn;
      }
    };
    let o = 0;
    for (let t = 0; t < nt; t += 1) {
      for (let b = 0; b < nb; b += 1) {
        point(t, b, P, N, o);
        point(t + 1, b, P, N, o + 3);
        point(t + 1, b + 1, P, N, o + 6);
        point(t, b, P, N, o + 9);
        point(t + 1, b + 1, P, N, o + 12);
        point(t, b + 1, P, N, o + 15);
        o += 18;
      }
    }
    out.positions = P.subarray(0, need);
    out.normals = N.subarray(0, need);
    return out;
  }

  // src/nd/workerCore.js
  //
  // The engine behind the message protocol. The plugin sends a command object
  // in, and the engine sends a report out.
  //
  // The same code runs in the two conditions. With a worker file, this module
  // runs in the worker, and `postMessage` carries the messages. With
  // `PhysiN.scripts.worker = null`, it runs in the main thread, and the
  // plugin calls `handle()` directly.
  //
  // THE BINARY LAYOUT. A report is a `Float32Array`. The plugin repeats these
  // two strides in `PhysiN.Scene`, and the two must always agree.
  //
  //   the world report, stride = 1 + n + r + n + k
  //     [0]    the type, WORLDREPORT
  //     [1]    the count of the bodies
  //     then, for each body:  id, position (n), rotor (r), velocity (n),
  //                           angular velocity (k)
  //
  //   the collision report, contactStride = 2 + n + n + 1
  //     [0]    the type, COLLISIONREPORT
  //     [1]    the count of the contacts
  //     then, for each contact:  id of a, id of b, normal (n), point (n), depth
  //
  // A JOINT SENDS NO REPORT. A third report would give a third stride to keep
  // in agreement, and the plugin does not need the impulse of a joint. The
  // state of a body that a joint holds comes back in the world report.
  //
  // A message that is not a `Float32Array` is a command object,
  // `{ cmd, params }`.
  /** The first number of a binary report. It says which report it is. */
  var MESSAGE_TYPES = { WORLDREPORT: 0, COLLISIONREPORT: 1 };
  /**
   * Makes the engine. It holds a `World`, and it obeys the commands.
   *
   * The commands: `init`, `addBody`, `removeBody`, `addConstraint`,
   * `removeConstraint`, `setConstraintParams`, `updateTransform`,
   * `setGravity`, `setFixedTimeStep`, `setParams`, `setLinearVelocity`,
   * `setAngularVelocity`, `applyCentralImpulse`, `applyImpulse`,
   * `applyCentralForce`, `applyForce`, `applyTorque`, `setMass`, `simulate`.
   *
   * `simulate` sends the two reports back. `init` and `addBody` send
   * `worldReady` and `objectReady`. A joint that the engine cannot build sends
   * `constraintFailed`, and a joint that breaks sends `constraintBroken`. An
   * unknown command sends `unknown`.
   *
   * @param {function(*): void} post sends a message back to the plugin
   * @returns {{world: World, handle: function(object): void}} the engine
   */
  function createEngine(post) {
    let world = null;
    let D = null;
    let stride = 0;
    let contactStride = 0;
    let worldReport = null;
    let collisionReport = null;
    const bodies = /* @__PURE__ */ new Map();
    const constraints = /* @__PURE__ */ new Map();
    const collisions = [];
    let fixedTimeStep = 1 / 60;
    /**
     * Builds a shape from the plain object that the plugin sent. The plugin
     * cannot send a shape object, because a message holds no functions.
     * @throws {Error} when the type of the shape is not known
     */
    function makeShape(def) {
      switch (def.type) {
        case "box":
          return HyperBox(D, def.halfExtents);
        case "sphere":
          return HyperSphere(D, def.radius);
        case "halfspace":
          return HalfSpace(D, def.normal, def.offset || 0);
        case "torus":
          return Torus(D, def.majorRadius, def.minorRadius, def.plane || [0, 1]);
        case "convex":
          return ConvexMesh(D, def.vertices, def.cells);
        default:
          throw new Error(`workerCore: the shape "${def.type}" is not known`);
      }
    }
    /**
     * Sends the state of each body, in the world report layout above. It uses
     * the same buffer again at each step, and it makes a larger one only when
     * the count of the bodies grows.
     */
    function reportWorld() {
      const { n, k, r } = D;
      const list = world.bodies;
      const need = 2 + list.length * stride;
      if (!worldReport || worldReport.length < need) worldReport = new Float32Array(need);
      worldReport[0] = MESSAGE_TYPES.WORLDREPORT;
      worldReport[1] = list.length;
      let o = 2;
      for (const b of list) {
        worldReport[o] = b.id;
        o += 1;
        for (let i = 0; i < n; i += 1) worldReport[o + i] = b.x[i];
        o += n;
        for (let i = 0; i < r; i += 1) worldReport[o + i] = b.R[i];
        o += r;
        for (let i = 0; i < n; i += 1) worldReport[o + i] = b.v[i];
        o += n;
        for (let i = 0; i < k; i += 1) worldReport[o + i] = b.w[i];
        o += k;
      }
      post(worldReport.subarray(0, need));
    }
    /**
     * Sends the contacts of this step, in the collision report layout above.
     * It sends nothing when there is no contact. It makes the list empty
     * again at the end.
     */
    function reportCollisions() {
      const { n } = D;
      if (collisions.length === 0) {
        collisions.length = 0;
        return;
      }
      const need = 2 + collisions.length * contactStride;
      if (!collisionReport || collisionReport.length < need) collisionReport = new Float32Array(need);
      collisionReport[0] = MESSAGE_TYPES.COLLISIONREPORT;
      collisionReport[1] = collisions.length;
      let o = 2;
      for (const c of collisions) {
        collisionReport[o] = c.a.id;
        collisionReport[o + 1] = c.b.id;
        o += 2;
        for (let i = 0; i < n; i += 1) collisionReport[o + i] = c.normal[i];
        o += n;
        for (let i = 0; i < n; i += 1) collisionReport[o + i] = c.point[i];
        o += n;
        collisionReport[o] = c.depth;
        o += 1;
      }
      collisions.length = 0;
      post(collisionReport.subarray(0, need));
    }
    /** The commands. `handle()` looks the name up here. */
    const commands = {
      /**
       * Makes the world, and works out the two strides of the reports. It
       * keeps only the first contact of each pair in one step, thus the
       * plugin does not send the same collision event many times.
       */
      init(params) {
        world = new World({
          dimensions: params.dimensions || 3,
          gravity: params.gravity,
          params: params.params
        });
        D = world.D;
        stride = 1 + D.n + D.r + D.n + D.k;
        contactStride = 2 + D.n + D.n + 1;
        fixedTimeStep = params.params && params.params.fixedTimeStep || 1 / 60;
        const seen = /* @__PURE__ */ new Set();
        world.on("collision", (a, b, c) => {
          const key = `${a.id}:${b.id}`;
          if (seen.has(key)) return;
          seen.add(key);
          collisions.push(c);
        });
        world.on("stepStart", () => seen.clear());
        world.on("constraintBroken", (j) => {
          constraints.delete(j.id);
          post({ cmd: "constraintBroken", params: { id: j.id } });
        });
        post({ cmd: "worldReady" });
      },
      addBody(def) {
        const body = new Body(D, {
          id: def.id,
          shape: makeShape(def.shape),
          mass: def.mass,
          position: def.position,
          rotor: def.rotor,
          velocity: def.velocity,
          angularVelocity: def.angularVelocity,
          friction: def.friction,
          restitution: def.restitution,
          linearDamping: def.linearDamping,
          angularDamping: def.angularDamping,
          linearFactor: def.linearFactor,
          angularFactor: def.angularFactor,
          allowSleep: def.allowSleep
        });
        bodies.set(def.id, body);
        world.addBody(body);
        post({ cmd: "objectReady", params: def.id });
      },
      removeBody(params) {
        const b = bodies.get(params.id);
        if (b) {
          world.removeBody(b);
          bodies.delete(params.id);
          for (const [cid, j] of constraints) {
            if (j.a === b || j.b === b) constraints.delete(cid);
          }
        }
      },
      /**
       * Adds a joint. `def.a` and `def.b` are the ids of the bodies, and a `b`
       * that is null makes a joint to the world. See `createConstraint` for the
       * other fields of `def`.
       *
       * The command sends `constraintFailed` when a body is not in the world.
       * A joint that goes before the body would be lost with no message.
       */
      addConstraint(def) {
        const a = bodies.get(def.a);
        const b = def.b === void 0 || def.b === null ? null : bodies.get(def.b);
        if (!a || def.b !== void 0 && def.b !== null && !b) {
          post({ cmd: "constraintFailed", params: { id: def.id, reason: "the body is not in the world" } });
          return;
        }
        let joint;
        try {
          joint = createConstraint(D, def, a, b, world.params);
        } catch (e) {
          post({ cmd: "constraintFailed", params: { id: def.id, reason: e.message } });
          return;
        }
        constraints.set(def.id, joint);
        world.addConstraint(joint);
      },
      removeConstraint(params) {
        const j = constraints.get(params.id);
        if (j) {
          world.removeConstraint(j);
          constraints.delete(params.id);
        }
      },
      /** Changes `enabled`, `collideConnected`, `rest` or `mode` of a joint. */
      setConstraintParams(params) {
        const j = constraints.get(params.id);
        if (!j) return;
        if (params.enabled !== void 0) {
          if (params.enabled && !j.enabled) resetConstraint(j);
          j.enabled = params.enabled;
        }
        if (params.rest !== void 0) j.rest = params.rest;
        if (params.mode !== void 0) j.mode = params.mode;
        if (params.lowerAngle !== void 0) j.lowerAngle = params.lowerAngle;
        if (params.upperAngle !== void 0) j.upperAngle = params.upperAngle;
        if (params.lowerAngle !== void 0 || params.upperAngle !== void 0) {
          j.hasLimit = j.lowerAngle > -Infinity || j.upperAngle < Infinity;
        }
        if (params.motorSpeed !== void 0) j.motorSpeed = params.motorSpeed;
        if (params.maxMotorTorque !== void 0) j.maxMotorTorque = params.maxMotorTorque;
        if (params.hertz !== void 0) j.hertz = params.hertz;
        if (params.damping !== void 0) j.damping = params.damping;
        if (params.breakForce !== void 0) j.breakForce = params.breakForce;
        if (params.breakTorque !== void 0) j.breakTorque = params.breakTorque;
        if (params.collideConnected !== void 0) {
          j.collideConnected = params.collideConnected;
          world._rebuildNoCollide();
        }
        j.a.wake();
        j.b.wake();
      },
      updateTransform(params) {
        const b = bodies.get(params.id);
        if (!b) return;
        if (params.position) b.x.set(params.position);
        if (params.rotor) b.R.set(params.rotor);
        b.updateDerived();
        b.wake();
        for (const j of constraints.values()) {
          if (j.a === b || j.b === b) resetConstraint(j);
        }
      },
      setGravity(g) {
        world.setGravity(g);
      },
      setFixedTimeStep(v) {
        fixedTimeStep = v;
      },
      setParams(p) {
        Object.assign(world.params, p);
      },
      setLinearVelocity(params) {
        const b = bodies.get(params.id);
        if (b) b.setLinearVelocity(params.value);
      },
      setAngularVelocity(params) {
        const b = bodies.get(params.id);
        if (b) {
          b.setAngularVelocity(params.value);
          b.wake();
        }
      },
      applyCentralImpulse(params) {
        const b = bodies.get(params.id);
        if (b) b.applyCentralImpulse(params.value);
      },
      applyImpulse(params) {
        const b = bodies.get(params.id);
        if (b) b.applyImpulse(params.value, params.offset);
      },
      applyCentralForce(params) {
        const b = bodies.get(params.id);
        if (b) b.applyCentralForce(params.value);
      },
      applyForce(params) {
        const b = bodies.get(params.id);
        if (b) b.applyForce(params.value, params.offset);
      },
      applyTorque(params) {
        const b = bodies.get(params.id);
        if (b) b.applyTorque(params.value);
      },
      /**
       * Changes the mass of a body. The inertia tensor comes from the mass at
       * build time, thus the command makes a new body with the same state,
       * and it puts that body in the place of the old one.
       */
      setMass(params) {
        const b = bodies.get(params.id);
        if (!b) return;
        const opts = { shape: b.shape, mass: params.value };
        const nb = new Body(D, Object.assign({ id: b.id }, opts));
        nb.x.set(b.x);
        nb.R.set(b.R);
        nb.v.set(b.v);
        nb.updateDerived();
        world.replaceBody(b, nb);
        bodies.set(b.id, nb);
      },
      /**
       * Moves the world ahead, then sends the two reports. `maxSubSteps` cuts
       * the time into that many parts.
       */
      simulate(params) {
        const step = params && params.timeStep ? params.timeStep : fixedTimeStep;
        const maxSub = params && params.maxSubSteps || 1;
        const h = step / maxSub;
        for (let i = 0; i < maxSub; i += 1) {
          world.emit("stepStart");
          world.step(h);
        }
        reportWorld();
        reportCollisions();
      }
    };
    return {
      get world() {
        return world;
      },
      handle(message) {
        if (!message || !message.cmd) return;
        const fn = commands[message.cmd];
        if (!fn) {
          post({ cmd: "unknown", params: message.cmd });
          return;
        }
        fn(message.params);
      }
    };
  }

  // src/physiN.js
  //
  // The three.js plugin. It gives the interface of physi.js: use
  // `PhysiN.Scene` in the place of `THREE.Scene`, and a `PhysiN` mesh class in
  // the place of `THREE.Mesh`. See README, section 8.
  //
  // The plugin does no physics. It holds the `n`-dimensional state of each
  // object in `object._physiN`, it sends the commands to the engine, and it
  // puts the answer back onto the three.js objects.
  //
  // The display depends on `n`:
  //   n = 3  the position goes to `object.position`, and the rotor becomes a
  //          quaternion.
  //   n = 4  the slice builds a new 3D surface at each frame.
  //   n > 4  it draws the first three axes of the position only.
  /**
   * Makes the `PhysiN` namespace for a given three.js.
   * @param {object} THREE the three.js namespace
   * @returns {object} the namespace
   * @throws {Error} when `THREE` is missing
   */
  function createPhysiN(THREE) {
    if (!THREE) throw new Error("PhysiN: give the THREE namespace");
    const PhysiN = {};
    PhysiN.scripts = { worker: null };
    PhysiN.version = "0.1.0";
    PhysiN.nd = nd_exports;
    PhysiN.slice = slice4_exports;
    let nextId2 = 0;
    let nextJointId = 0;
    /**
     * The event system of physi.js: `addEventListener`,
     * `removeEventListener` and `dispatchEvent`. `Eventable.make(Klass)` puts
     * the three functions on a class.
     *
     * A scene sends `ready` and `update`. An object sends `ready` and
     * `collision`.
     */
    function Eventable() {
      this._eventListeners = {};
    }
    Eventable.prototype.addEventListener = function(name, fn) {
      (this._eventListeners[name] = this._eventListeners[name] || []).push(fn);
    };
    Eventable.prototype.removeEventListener = function(name, fn) {
      const list = this._eventListeners[name];
      if (!list) return false;
      const i = list.indexOf(fn);
      if (i < 0) return false;
      list.splice(i, 1);
      return true;
    };
    Eventable.prototype.dispatchEvent = function(name, ...args) {
      const list = this._eventListeners[name];
      if (list) for (const fn of list.slice()) fn.apply(this, args);
    };
    Eventable.make = function(Klass) {
      Klass.prototype.addEventListener = Eventable.prototype.addEventListener;
      Klass.prototype.removeEventListener = Eventable.prototype.removeEventListener;
      Klass.prototype.dispatchEvent = Eventable.prototype.dispatchEvent;
    };
    PhysiN.Eventable = Eventable;
    /**
     * A copy of a three.js material, with the friction and the bounce on it.
     * A mesh that you build with this material takes the two values.
     * @param {number} [friction] the friction. Default 0.5.
     * @param {number} [restitution] the bounce, 0 to 1. Default 0.1.
     */
    PhysiN.createMaterial = function(material, friction, restitution) {
      const m = material.clone();
      m._physiN = {
        friction: friction === void 0 ? 0.5 : friction,
        restitution: restitution === void 0 ? 0.1 : restitution
      };
      return m;
    };
    function defaultGravity(n) {
      const g = new Array(n).fill(0);
      g[1] = -9.81;
      return g;
    }
    const _m4 = new THREE.Matrix4();
    /**
     * Builds the bounding sphere of a geometry that the slice made. three.js
     * needs it to know if the object is in view. The slice changes the
     * geometry at each frame, thus this runs at each frame too.
     */
    function setBoundingSphere(geometry, positions, count) {
      let cx = 0;
      let cy = 0;
      let cz = 0;
      const verts = count / 3;
      for (let i = 0; i < count; i += 3) {
        cx += positions[i];
        cy += positions[i + 1];
        cz += positions[i + 2];
      }
      cx /= verts;
      cy /= verts;
      cz /= verts;
      let r2 = 0;
      for (let i = 0; i < count; i += 3) {
        const dx = positions[i] - cx;
        const dy = positions[i + 1] - cy;
        const dz = positions[i + 2] - cz;
        const d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > r2) r2 = d2;
      }
      if (!geometry.boundingSphere) geometry.boundingSphere = new THREE.Sphere();
      geometry.boundingSphere.center.set(cx, cy, cz);
      geometry.boundingSphere.radius = Math.sqrt(r2);
    }
    /**
     * A three.js scene that also holds a world of physics. Use it in the place
     * of `THREE.Scene`.
     *
     * It starts a worker when `PhysiN.scripts.worker` names a file. If that is
     * null, the engine runs in the main thread. The interface is the same in
     * the two conditions.
     *
     * Call `simulate(timeStep, maxSubSteps)` at each frame. The answer comes
     * back later, and the scene then sends `update`.
     */
    PhysiN.Scene = class PhysiNScene extends THREE.Scene {
      /**
       * @param {object} [params] the options
       * @param {number} [params.dimensions] the count of dimensions. Default 3.
       * @param {number[]} [params.gravity] the gravity, of length `n`
       * @param {object} [params.params] the tolerances. See `defaultParams`.
       * @param {number} [params.sliceW] the value of `w` of the 4D cut. Default 0.
       * @param {string} [params.workerType] `classic` or `module`. Default `classic`.
       */
      constructor(params = {}) {
        super();
        Eventable.call(this);
        const self = this;
        this._n = params.dimensions || 3;
        this._D = dims(this._n);
        this._objects = {};
        this._constraints = {};
        this._isSimulating = false;
        this._stride = 1 + this._D.n + this._D.r + this._D.n + this._D.k;
        this._contactStride = 2 + this._D.n + this._D.n + 1;
        this.sliceW = params.sliceW || 0;
        const onMessage = (data) => self._onMessage(data);
        if (PhysiN.scripts.worker) {
          this._worker = new Worker(PhysiN.scripts.worker, { type: params.workerType || "classic" });
          this._worker.onmessage = (event) => onMessage(event.data);
          this._post = (msg, transfer) => self._worker.postMessage(msg, transfer || []);
        } else {
          this._engine = createEngine((data) => onMessage(data));
          this._post = (msg) => self._engine.handle(msg);
        }
        this.execute("init", {
          dimensions: this._n,
          gravity: params.gravity || defaultGravity(this._n),
          params: params.params || {}
        });
      }
    };
    Eventable.make(PhysiN.Scene);
    /** Sends a command to the engine, in the worker or in the main thread. */
    PhysiN.Scene.prototype.execute = function(cmd, params) {
      this._post({ cmd, params });
    };
    /**
     * Takes a message from the engine. A `Float32Array` is a binary report,
     * and its first number says which one. All other messages are command
     * objects.
     *
     * The `ready` events go out later, and not at once. Thus a listener that
     * the caller adds after the constructor still gets the event.
     */
    PhysiN.Scene.prototype._onMessage = function(data) {
      if (data instanceof Float32Array || data && data.buffer && data.BYTES_PER_ELEMENT === 4) {
        if (data[0] === MESSAGE_TYPES.WORLDREPORT) this._updateScene(data);
        else if (data[0] === MESSAGE_TYPES.COLLISIONREPORT) this._updateCollisions(data);
        return;
      }
      const later = (fn) => typeof queueMicrotask === "function" ? queueMicrotask(fn) : setTimeout(fn, 0);
      if (data && data.cmd === "worldReady") {
        later(() => this.dispatchEvent("ready"));
        return;
      }
      if (data && data.cmd === "objectReady") {
        later(() => {
          const o = this._objects[data.params];
          if (o) o.dispatchEvent("ready");
        });
        return;
      }
      if (data && data.cmd === "constraintBroken") {
        const j = this._constraints[data.params.id];
        if (j) {
          delete this._constraints[data.params.id];
          j.scene = null;
          later(() => j.dispatchEvent("broken"));
        }
        return;
      }
      if (data && data.cmd === "constraintFailed") {
        const j = this._constraints[data.params.id];
        if (j) delete this._constraints[data.params.id];
        console.error("PhysiN: the engine did not take the joint " + data.params.id + ": " + data.params.reason);
      }
    };
    /**
     * Reads a world report, and puts the state of each body onto its object.
     * The layout is the one of `workerCore`, and `this._stride` must always
     * agree with the stride there.
     */
    PhysiN.Scene.prototype._updateScene = function(data) {
      const { n, k, r } = this._D;
      const count = data[1];
      for (let m = 0; m < count; m += 1) {
        let o = 2 + m * this._stride;
        const id = data[o];
        o += 1;
        const object = this._objects[id];
        if (!object) {
          continue;
        }
        const st = object._physiN;
        for (let i = 0; i < n; i += 1) st.position[i] = data[o + i];
        o += n;
        for (let i = 0; i < r; i += 1) st.rotor[i] = data[o + i];
        o += r;
        for (let i = 0; i < n; i += 1) st.velocity[i] = data[o + i];
        o += n;
        for (let i = 0; i < k; i += 1) st.angularVelocity[i] = data[o + i];
        this._placeObject(object);
      }
      this._isSimulating = false;
      this.dispatchEvent("update");
    };
    /**
     * Reads a collision report, and sends a `collision` event to the two
     * objects. The normal that the second object gets points the other way,
     * thus each object sees the normal that comes to it.
     */
    PhysiN.Scene.prototype._updateCollisions = function(data) {
      const { n } = this._D;
      const count = data[1];
      for (let m = 0; m < count; m += 1) {
        const o = 2 + m * this._contactStride;
        const a = this._objects[data[o]];
        const b = this._objects[data[o + 1]];
        if (!a || !b) continue;
        const normal = Array.from(data.subarray(o + 2, o + 2 + n));
        const point = Array.from(data.subarray(o + 2 + n, o + 2 + 2 * n));
        const depth = data[o + 2 + 2 * n];
        a.dispatchEvent("collision", b, normal, point, depth);
        b.dispatchEvent("collision", a, normal.map((x) => -x), point, depth);
      }
    };
    /**
     * Puts the state of one object onto its three.js object. In 3 dimensions
     * it sets the position and the quaternion. In 4 dimensions it cuts the
     * body. In 5 dimensions and more it sets the first three axes of the
     * position only.
     */
    PhysiN.Scene.prototype._placeObject = function(object) {
      const { n } = this._D;
      const st = object._physiN;
      if (n === 3) {
        object.position.set(st.position[0], st.position[1], st.position[2]);
        const M = rotorMatrix(this._D, st.rotor, st.matrix);
        _m4.set(
          M[0],
          M[1],
          M[2],
          0,
          M[3],
          M[4],
          M[5],
          0,
          M[6],
          M[7],
          M[8],
          0,
          0,
          0,
          0,
          1
        );
        object.quaternion.setFromRotationMatrix(_m4);
        return;
      }
      if (n === 4) {
        this._sliceObject(object);
        return;
      }
      object.position.set(st.position[0], st.position[1], st.position[2]);
    };
    /**
     * Builds the 3D cut of a 4D body, at `this.sliceW`. There are four
     * conditions, from the fastest to the slowest:
     *   - a hypersphere: the cut is a ball. It only scales the geometry.
     *   - a half space: it draws the display geometry as it is.
     *   - a torus on the fast path: `torusSliceTriangles` makes a smooth 3D
     *     torus, and no mesh is necessary.
     *   - all other shapes: it changes the frame of the mesh, then
     *     `sliceTetrahedra` cuts it.
     * A body that does not meet the cut becomes not visible.
     */
    PhysiN.Scene.prototype._sliceObject = function(object) {
      const st = object._physiN;
      const D = this._D;
      if (st.shape.type === "sphere") {
        const rad = sliceHyperSphereRadius(st.shape.radius, st.position[3], this.sliceW);
        if (rad < 0) {
          object.visible = false;
          return;
        }
        object.visible = true;
        object.position.set(st.position[0], st.position[1], st.position[2]);
        object.scale.setScalar(rad);
        return;
      }
      if (st.shape.type === "halfspace") {
        object.visible = true;
        return;
      }
      if (st.shape.type === "torus" && (torusFastPath(D, st) || !st.mesh)) {
        const rho = sliceTorusMinorRadius(st.shape.minorRadius, st.position[3], this.sliceW);
        if (rho <= 0) {
          object.visible = false;
          return;
        }
        const M2 = rotorMatrix(D, st.rotor, st.matrix);
        st.sink = torusSliceTriangles(
          st.sink || {},
          st.shape,
          M2,
          st.position,
          this.sliceW - st.position[3],
          rho,
          st.sliceSegments
        );
        object.visible = true;
        object.position.set(0, 0, 0);
        object.quaternion.set(0, 0, 0, 1);
        writeSliceBuffers(object, st.sink);
        return;
      }
      if (!st.mesh) {
        object.visible = false;
        return;
      }
      const M = rotorMatrix(D, st.rotor, st.matrix);
      st.world = transformVertices(st.mesh.vertices, M, st.position, st.world);
      const sink = { positions: [], normals: [] };
      sliceTetrahedra(st.world, st.mesh.cells, this.sliceW, sink);
      if (sink.positions.length === 0) {
        object.visible = false;
        return;
      }
      object.visible = true;
      object.position.set(0, 0, 0);
      object.quaternion.set(0, 0, 0, 1);
      writeSliceBuffers(object, sink);
    };
    /** True after the first warning about a torus with no mesh. */
    let torusWarned = false;
    /**
     * True when the smooth cut of a torus is exact. Two conditions must hold:
     * the major plane must not hold the axis `w`, and the body must not turn
     * out of the cut, thus the corner `M[15]` of the rotation matrix must be 1.
     *
     * When the two do not hold, the exact cut needs a mesh. The function
     * writes one warning if the torus has no mesh.
     */
    function torusFastPath(D, st) {
      const pl = st.shape.plane;
      if (pl[0] > 2 || pl[1] > 2) return false;
      const M = rotorMatrix(D, st.rotor, st.matrix);
      if (Math.abs(M[15] - 1) < 1e-6) return true;
      if (!st.mesh && !torusWarned) {
        torusWarned = true;
        if (typeof console !== "undefined" && console.warn) {
          console.warn("PhysiN: this torus turns out of the slice, and it has no mesh. The display is not exact. Build it with { mesh: true }.");
        }
      }
      return false;
    }
    /**
     * Writes the triangles of a cut into the geometry of an object. It makes
     * a larger buffer only when the count of the triangles grows, and it
     * rounds the size up to a multiple of 192. Thus a frame makes almost no
     * garbage. The draw range then hides the part of the buffer that the cut
     * does not use.
     */
    function writeSliceBuffers(object, sink) {
      const need = sink.positions.length;
      const g = object.geometry;
      let attr = g.getAttribute("position");
      if (!attr || attr.array.length < need) {
        const cap = Math.ceil(Math.max(need, 192) / 192) * 192;
        g.setAttribute("position", new THREE.BufferAttribute(new Float32Array(cap), 3));
        g.setAttribute("normal", new THREE.BufferAttribute(new Float32Array(cap), 3));
        attr = g.getAttribute("position");
      }
      attr.array.set(sink.positions);
      g.getAttribute("normal").array.set(sink.normals);
      attr.needsUpdate = true;
      g.getAttribute("normal").needsUpdate = true;
      g.setDrawRange(0, need / 3);
      setBoundingSphere(g, sink.positions, need);
    }
    /**
     * Draws each object again at the value `sliceW` of the cut. Call it after
     * you change `sliceW`, thus the display agrees with the new value before
     * the next step.
     */
    PhysiN.Scene.prototype.refreshSlice = function() {
      for (const id in this._objects) {
        if (Object.prototype.hasOwnProperty.call(this._objects, id)) {
          this._placeObject(this._objects[id]);
        }
      }
    };
    /**
     * Adds an object to the scene. When the object has physics, it also sends
     * `addBody` to the engine. A normal three.js object goes in as usual.
     *
     * A shape that gives no count of dimensions, as a ball does, takes the
     * count of the scene here.
     *
     * @throws {Error} when the count of the mesh is not the count of the scene
     */
    PhysiN.Scene.prototype.add = function(object) {
      THREE.Scene.prototype.add.call(this, object);
      if (!object._physiN) return;
      const st = resolveDims(object._physiN, this._n);
      if (st.D.n !== this._n) {
        throw new Error(
          `PhysiN.Scene.add: the mesh has ${st.D.n} dimensions and the scene has ${this._n}`
        );
      }
      this._objects[st.id] = object;
      this.execute("addBody", {
        id: st.id,
        shape: st.shape,
        mass: st.mass,
        position: Array.from(st.position),
        rotor: Array.from(st.rotor),
        friction: st.friction,
        restitution: st.restitution,
        linearDamping: st.linearDamping,
        angularDamping: st.angularDamping,
        allowSleep: st.allowSleep
      });
    };
    /**
     * Takes an object out of the scene, and its body out of the engine. It also
     * takes out each joint that holds the object, because a joint that stays
     * would pull a body that is not there.
     */
    PhysiN.Scene.prototype.remove = function(object) {
      THREE.Scene.prototype.remove.call(this, object);
      if (!object._physiN) return;
      const id = object._physiN.id;
      delete this._objects[id];
      for (const key of Object.keys(this._constraints)) {
        const j = this._constraints[key];
        if (j.a === object || j.b === object) {
          delete this._constraints[key];
          j.scene = null;
        }
      }
      this.execute("removeBody", { id });
    };
    /**
     * Adds a joint. Give a `PhysiN.PointJoint`, `PhysiN.DistanceJoint`,
     * `PhysiN.FixedJoint`, `PhysiN.HingeJoint` or `PhysiN.SubspaceJoint`.
     *
     * CALL THIS AFTER `scene.add` OF THE TWO OBJECTS. The engine builds a joint
     * from the ids of the bodies, and it does not hold a joint that comes
     * first.
     *
     * @throws {Error} when an object is not in this scene
     */
    PhysiN.Scene.prototype.addConstraint = function(joint) {
      const a = joint.a;
      const b = joint.b;
      if (!a || !a._physiN || this._objects[a._physiN.id] !== a) {
        throw new Error("PhysiN.Scene.addConstraint: call scene.add of the first object first");
      }
      if (b && (!b._physiN || this._objects[b._physiN.id] !== b)) {
        throw new Error("PhysiN.Scene.addConstraint: call scene.add of the second object first");
      }
      joint.scene = this;
      this._constraints[joint.id] = joint;
      this.execute("addConstraint", joint.toMessage());
      return joint;
    };
    /** Takes a joint out of the scene and out of the engine. */
    PhysiN.Scene.prototype.removeConstraint = function(joint) {
      if (this._constraints[joint.id] !== joint) return;
      delete this._constraints[joint.id];
      joint.scene = null;
      this.execute("removeConstraint", { id: joint.id });
    };
    PhysiN.Scene.prototype.setGravity = function(g) {
      this.execute("setGravity", g);
    };
    PhysiN.Scene.prototype.setFixedTimeStep = function(v) {
      this.execute("setFixedTimeStep", v);
    };
    /** Changes the tolerances while the world runs. See `defaultParams`. */
    PhysiN.Scene.prototype.setSolverParams = function(p) {
      this.execute("setParams", p);
    };
    /**
     * Moves the physics ahead. Call it one time at each frame.
     *
     * It first sends the transform of each object that the program moved by
     * hand. It then sends `simulate`. The answer comes back later, thus the
     * new positions are not ready when this function gives its result.
     *
     * @param {number} [timeStep] the time in seconds. The default is the fixed
     *   time step of the world.
     * @param {number} [maxSubSteps] the count of the parts of the step
     * @returns {boolean} false when a step is still running. The function then
     *   does nothing.
     */
    PhysiN.Scene.prototype.simulate = function(timeStep, maxSubSteps) {
      if (this._isSimulating) return false;
      this._isSimulating = true;
      for (const id in this._objects) {
        if (!Object.prototype.hasOwnProperty.call(this._objects, id)) continue;
        const object = this._objects[id];
        const st = object._physiN;
        if (st.dirty) {
          this.execute("updateTransform", {
            id: st.id,
            position: Array.from(st.position),
            rotor: Array.from(st.rotor)
          });
          st.dirty = false;
        }
      }
      this.execute("simulate", { timeStep, maxSubSteps });
      return true;
    };
    /**
     * Builds `object._physiN`, the `n`-dimensional state of a mesh. three.js
     * holds 3 dimensions only, thus the plugin keeps the true state here.
     *
     * The fields:
     *   `id`               the id that the engine uses
     *   `shape`            the plain shape object, for the message
     *   `mass`             the mass. 0 makes the body static.
     *   `position`         the position, `n`
     *   `rotor`            the orientation, `r`
     *   `velocity`         the last velocity from the engine, `n`
     *   `angularVelocity`  the last angular velocity, `k`
     *   `matrix`           a work buffer for the rotation matrix, `n` by `n`
     *   `mesh`             the tetrahedra for the 4D cut, or null
     *   `world`            a work buffer for the vertices in the world frame
     *   `dirty`            true when the program moved the object by hand.
     *                      `simulate()` then sends the new transform.
     *   `D`                the tables of `n`, or null while the count waits
     *   `turns`            the turns that came before the count, or null
     *
     * A SHAPE THAT GIVES NO COUNT. A hyperbox gives `n` with its half extents,
     * and a half space with its normal. A ball gives nothing: a ball of 3
     * dimensions and a ball of 5 have the same radius. Thus `shape.n` can be
     * absent, and the state then waits: `D` is null, and there are no arrays.
     * `setPositionN` gives the count with the length of its argument, and
     * `Scene.add` gives the count of the scene. `resolveDims` builds the rest.
     *
     * @returns {object} the state
     */
    function initState(object, shape, mass, material) {
      const mat = material && material._physiN || {};
      object._physiN = {
        id: nextId2++,
        shape,
        mass: mass === void 0 ? 1 : mass,
        position: null,
        rotor: null,
        velocity: null,
        angularVelocity: null,
        matrix: null,
        friction: mat.friction,
        restitution: mat.restitution,
        linearDamping: 0,
        angularDamping: 0,
        allowSleep: true,
        dirty: false,
        mesh: null,
        world: null,
        turns: null,
        D: null
      };
      if (typeof shape.n === "number") resolveDims(object._physiN, shape.n);
      return object._physiN;
    }
    /**
     * Gives the state its count of dimensions, and builds the arrays of that
     * count. It plays back the turns that `rotateInPlane` held, and it writes
     * the count into the shape, thus the `addBody` message always names it.
     *
     * It does nothing when the state has the count already.
     *
     * @param {object} st the state, `object._physiN`
     * @param {number} n the count of dimensions
     */
    function resolveDims(st, n) {
      if (st.D) return st;
      const D = dims(n);
      st.D = D;
      st.shape.n = n;
      st.position = new Float64Array(n);
      st.rotor = rotorIdentity(D);
      st.velocity = new Float64Array(n);
      st.angularVelocity = new Float64Array(D.k);
      st.matrix = new Float64Array(n * n);
      if (st.turns) {
        for (const [i, j, angle] of st.turns) {
          st.rotor.set(rotorMul(D, rotorFromPlane(D, i, j, angle), st.rotor));
        }
        st.turns = null;
      }
      return st;
    }
    /**
     * The state of a mesh, with its count of dimensions.
     * @throws {Error} when the count is not known yet
     */
    function readyState(object, what) {
      const st = object._physiN;
      if (!st.D) {
        throw new Error(
          `PhysiN.${what}: the count of dimensions of this mesh is not known. ` +
          "Call setPositionN, or add the mesh to a scene, first."
        );
      }
      return st;
    }
    /**
     * The base of all of the meshes with physics. Use one of the classes below
     * it, or give your own shape object here.
     *
     * The functions that apply a force or an impulse go to the parent scene.
     * Thus the mesh must be in a scene before you call them.
     */
    PhysiN.Mesh = class PhysiNMesh extends THREE.Mesh {
      /**
       * @param {object} geometry the three.js geometry for the display
       * @param {object} material the material. `PhysiN.createMaterial` puts
       *   the friction and the bounce on it.
       * @param {number} [mass] the mass. 0 makes the body static. Default 1.
       * @param {object} shape the shape, as a plain object
       */
      constructor(geometry, material, mass, shape) {
        super(geometry, material);
        Eventable.call(this);
        initState(this, shape, mass, material);
        this.castShadow = true;
        this.receiveShadow = true;
      }
    };
    Eventable.make(PhysiN.Mesh);
    /**
     * Sets the position in `n` dimensions. `simulate()` sends the new value to
     * the engine at the next frame.
     *
     * When the shape gives no count of dimensions, as a ball does, the length
     * of `p` gives it.
     *
     * @param {ArrayLike<number>} p the position, of length `n`
     */
    PhysiN.Mesh.prototype.setPositionN = function(p) {
      const st = resolveDims(this._physiN, p.length);
      st.position.set(p);
      st.dirty = true;
    };
    /** The position in `n` dimensions. Do not change the array that it gives. */
    PhysiN.Mesh.prototype.getPositionN = function() {
      return readyState(this, "getPositionN").position;
    };
    /** The orientation as a rotor, of length `r`. It is not a quaternion. */
    PhysiN.Mesh.prototype.getRotor = function() {
      return readyState(this, "getRotor").rotor;
    };
    /**
     * Turns the body through `angle` radians in the plane of the axes `i` and
     * `j`. This is the general form of a turn about an axis: in 4 dimensions a
     * turn has a plane, and not an axis.
     *
     * When the shape gives no count of dimensions, the state holds the turn
     * until the count comes, and `resolveDims` plays it back then.
     */
    PhysiN.Mesh.prototype.rotateInPlane = function(i, j, angle) {
      const st = this._physiN;
      if (!st.D) {
        if (!st.turns) st.turns = [];
        st.turns.push([i, j, angle]);
        st.dirty = true;
        return;
      }
      const S = rotorFromPlane(st.D, i, j, angle);
      st.rotor.set(rotorMul(st.D, S, st.rotor));
      st.dirty = true;
    };
    /**
     * Makes a function that sends one command to the engine, for a mesh. The
     * functions that apply a force, an impulse and a velocity all use it.
     * The function does nothing when the mesh has no parent scene.
     */
    function meshCommand(name) {
      return function(value, offset) {
        if (!this.parent || !this.parent.execute) return;
        this.parent.execute(name, {
          id: this._physiN.id,
          value: Array.from(value),
          offset: offset ? Array.from(offset) : void 0
        });
      };
    }
    PhysiN.Mesh.prototype.applyCentralImpulse = meshCommand("applyCentralImpulse");
    PhysiN.Mesh.prototype.applyImpulse = meshCommand("applyImpulse");
    PhysiN.Mesh.prototype.applyCentralForce = meshCommand("applyCentralForce");
    PhysiN.Mesh.prototype.applyForce = meshCommand("applyForce");
    PhysiN.Mesh.prototype.applyTorque = meshCommand("applyTorque");
    PhysiN.Mesh.prototype.setLinearVelocity = meshCommand("setLinearVelocity");
    PhysiN.Mesh.prototype.setAngularVelocity = meshCommand("setAngularVelocity");
    /** The linear velocity of the last report, of length `n`. */
    PhysiN.Mesh.prototype.getLinearVelocity = function() {
      return readyState(this, "getLinearVelocity").velocity;
    };
    /** The angular velocity of the last report. It is a bivector of `k`. */
    PhysiN.Mesh.prototype.getAngularVelocity = function() {
      return readyState(this, "getAngularVelocity").angularVelocity;
    };
    /**
     * A 3D box. The half extents come from the bounding box of the geometry.
     * Use it in the place of `Physijs.BoxMesh`.
     */
    PhysiN.BoxMesh = class PhysiNBoxMesh extends PhysiN.Mesh {
      constructor(geometry, material, mass) {
        geometry.computeBoundingBox();
        const b = geometry.boundingBox;
        super(geometry, material, mass, {
          type: "box",
          n: 3,
          halfExtents: [(b.max.x - b.min.x) / 2, (b.max.y - b.min.y) / 2, (b.max.z - b.min.z) / 2]
        });
      }
    };
    /** A 3D ball. The radius comes from the bounding sphere of the geometry. */
    PhysiN.SphereMesh = class PhysiNSphereMesh extends PhysiN.Mesh {
      constructor(geometry, material, mass) {
        geometry.computeBoundingSphere();
        super(
          geometry,
          material,
          mass,
          { type: "sphere", n: 3, radius: geometry.boundingSphere.radius }
        );
      }
    };
    /**
     * A 3D half space, for the ground or a wall. It is always static. The
     * default normal points up, on the axis 1.
     */
    PhysiN.PlaneMesh = class PhysiNPlaneMesh extends PhysiN.Mesh {
      constructor(geometry, material, normal, offset) {
        super(
          geometry,
          material,
          0,
          { type: "halfspace", n: 3, normal: normal || [0, 1, 0], offset: offset || 0 }
        );
      }
    };
    const SphereGeometryClass = THREE.SphereGeometry || THREE.SphereBufferGeometry;
    /**
     * An empty geometry, with a position attribute and a normal attribute.
     * Each shape that the 4D cut draws starts with one, and
     * `writeSliceBuffers` then fills it at each frame.
     */
    function emptyGeometry() {
      const g = new THREE.BufferGeometry();
      g.setAttribute("position", new THREE.Float32BufferAttribute([], 3));
      g.setAttribute("normal", new THREE.Float32BufferAttribute([], 3));
      return g;
    }
    /**
     * A box of 4 dimensions and more. In 4 dimensions this is a tesseract, and
     * the class builds the tetrahedra that the cut needs.
     *
     * The geometry starts empty, thus give the half extents and not a
     * geometry. `frustumCulled` is false, because the cut changes the
     * geometry at each frame.
     */
    PhysiN.HyperBoxMesh = class PhysiNHyperBoxMesh extends PhysiN.Mesh {
      constructor(halfExtents, material, mass) {
        const n = halfExtents.length;
        super(
          emptyGeometry(),
          material,
          mass,
          { type: "box", n, halfExtents: Array.from(halfExtents) }
        );
        if (n === 4) this._physiN.mesh = hyperBoxMesh(dims(4), Float64Array.from(halfExtents));
        this.frustumCulled = false;
      }
    };
    /**
     * A ball of 4 dimensions and more. The cut of a 4-ball is a 3-ball, thus
     * the class uses a normal sphere geometry of the radius 1, and the cut
     * only scales it. No mesh is necessary.
     *
     * A ball gives no count of dimensions: a ball of 3 dimensions and a ball
     * of 5 have the same radius. Thus the shape holds no `n`, and the count
     * comes from the first `setPositionN` or from `scene.add`. A ball works
     * in a scene of any count.
     */
    PhysiN.HyperSphereMesh = class PhysiNHyperSphereMesh extends PhysiN.Mesh {
      constructor(radius, material, mass, segments = 24) {
        super(
          new SphereGeometryClass(1, segments, Math.max(2, segments / 2)),
          material,
          mass,
          { type: "sphere", radius }
        );
        this.frustumCulled = false;
      }
    };
    /**
     * Changes a three.js geometry into the vertices and the triangles that
     * `ConvexMesh` needs. It puts the vertices that are at almost the same
     * point together, at 1e-5 of length, because a geometry holds a separate
     * vertex for each corner of each face. It drops a triangle that has two
     * equal corners.
     * @throws {Error} when the geometry has no position attribute
     */
    function geometryToMesh(geometry) {
      const attr = geometry.getAttribute && geometry.getAttribute("position");
      if (!attr) throw new Error("PhysiN.ConvexMesh: the geometry has no position attribute");
      const arr = attr.array;
      const index = geometry.index ? geometry.index.array : null;
      const count = index ? index.length : arr.length / 3;
      const map = /* @__PURE__ */ new Map();
      const vertices = [];
      const cells = [];
      const idOf = (v) => {
        const x = arr[v * 3];
        const y = arr[v * 3 + 1];
        const z = arr[v * 3 + 2];
        const key = `${Math.round(x * 1e5)},${Math.round(y * 1e5)},${Math.round(z * 1e5)}`;
        let id = map.get(key);
        if (id === void 0) {
          id = vertices.length / 3;
          vertices.push(x, y, z);
          map.set(key, id);
        }
        return id;
      };
      for (let t = 0; t + 2 < count; t += 3) {
        const a = idOf(index ? index[t] : t);
        const b = idOf(index ? index[t + 1] : t + 1);
        const c = idOf(index ? index[t + 2] : t + 2);
        if (a === b || b === c || a === c) continue;
        cells.push(a, b, c);
      }
      return { vertices, cells };
    }
    /**
     * A 3D convex body from a three.js geometry. The geometry must be convex,
     * and its origin must be inside it. This shape touches a half space and a
     * hypersphere only. See `collide()`.
     */
    PhysiN.ConvexMesh = class PhysiNConvexMesh extends PhysiN.Mesh {
      constructor(geometry, material, mass) {
        const { vertices, cells } = geometryToMesh(geometry);
        super(geometry, material, mass, {
          type: "convex",
          n: 3,
          vertices,
          cells
        });
      }
    };
    /**
     * A 3D torus. Give the same radii as the geometry, thus the physics and
     * the display agree. The default major plane is the plane of the axes 0
     * and 1.
     */
    PhysiN.TorusMesh = class PhysiNTorusMesh extends PhysiN.Mesh {
      constructor(geometry, material, mass, majorRadius, minorRadius, plane = [0, 1]) {
        super(geometry, material, mass, {
          type: "torus",
          n: 3,
          majorRadius,
          minorRadius,
          plane: [plane[0], plane[1]],
          otherAxes: [0, 1, 2].filter((a) => a !== plane[0] && a !== plane[1])
        });
      }
    };
    /**
     * A torus of 4 dimensions and more.
     *
     * The cut has two paths. The fast path makes a smooth 3D torus, and it is
     * exact while the body does not turn out of the cut. The mesh path cuts
     * the tetrahedra, and it is always exact but it costs more. The class
     * builds the mesh by default, thus the display is always correct.
     *
     * @param {object} [options] the options
     * @param {number} [options.dimensions] the count of dimensions. Default 4.
     * @param {boolean} [options.mesh] false leaves the mesh out. The display
     *   is then not exact while the torus turns out of the cut.
     * @param {number[]} [options.meshSegments] the divisions of the mesh.
     *   Default `[12, 6, 8]`. A large value costs much.
     * @param {number[]} [options.sliceSegments] the divisions of the fast
     *   path. Default `[40, 20]`.
     */
    PhysiN.HyperTorusMesh = class PhysiNHyperTorusMesh extends PhysiN.Mesh {
      constructor(majorRadius, minorRadius, plane, material, mass, options = {}) {
        const n = options.dimensions || 4;
        const pl = [plane[0], plane[1]];
        const others = [];
        for (let a = 0; a < n; a += 1) if (a !== pl[0] && a !== pl[1]) others.push(a);
        super(emptyGeometry(), material, mass, {
          type: "torus",
          n,
          majorRadius,
          minorRadius,
          plane: pl,
          otherAxes: others
        });
        const wantMesh = options.mesh !== void 0 ? options.mesh : true;
        if (wantMesh && n === 4) {
          this._physiN.mesh = torusMesh(
            dims(4),
            majorRadius,
            minorRadius,
            pl,
            options.meshSegments || [12, 6, 8]
          );
        }
        this._physiN.sliceSegments = options.sliceSegments || [40, 20];
        this.frustumCulled = false;
      }
    };
    /**
     * A half space of 4 dimensions and more, for the ground. It is always
     * static. Give a display geometry, because a half space has no mesh and
     * the cut cannot build one.
     */
    PhysiN.HyperPlaneMesh = class PhysiNHyperPlaneMesh extends PhysiN.Mesh {
      constructor(normal, offset, material, displayGeometry) {
        super(
          displayGeometry || emptyGeometry(),
          material,
          0,
          { type: "halfspace", n: normal.length, normal: Array.from(normal), offset }
        );
        this._physiN.mesh = null;
      }
    };
    /**
     * A convex body of 4 dimensions, from your own mesh of tetrahedra. The
     * same mesh does the physics and the cut.
     * @param {ArrayLike<number>} vertices the vertices, 4 numbers for each
     * @param {ArrayLike<number>} cells the tetrahedra, 4 indices for each
     */
    PhysiN.HyperMesh = class PhysiNHyperMesh extends PhysiN.Mesh {
      constructor(vertices, cells, material, mass) {
        super(emptyGeometry(), material, mass, {
          type: "convex",
          n: 4,
          vertices: Array.from(vertices),
          cells: Array.from(cells)
        });
        this._physiN.mesh = {
          vertices: Float64Array.from(vertices),
          cells: Int32Array.from(cells)
        };
        this.frustumCulled = false;
      }
    };
    /**
     * The joint that the five classes below share. It holds the two objects and
     * the fields of the message, and `Scene.addConstraint` sends it.
     *
     * THE ANCHORS ARE ALWAYS IN THE LOCAL FRAME OF THEIR OBJECT. A null second
     * object makes a joint to the world, and the anchor of the world then takes
     * the world position of the anchor of the first object.
     *
     * The two objects must be in the scene before you add the joint. See
     * `Scene.addConstraint`.
     */
    PhysiN.Constraint = class PhysiNConstraint {
      /**
       * @param {string} type point, distance, fixed, hinge or subspace
       * @param {object} a the first object
       * @param {object} [b] the second object, or null for the world
       * @param {object} [opts] the other fields. See `createConstraint`.
       */
      constructor(type, a, b, opts = {}) {
        Eventable.call(this);
        this.id = nextJointId++;
        this.type = type;
        this.a = a;
        this.b = b || null;
        this.scene = null;
        this.opts = opts;
      }
      /** The command object that `Scene.addConstraint` sends to the engine. */
      toMessage() {
        const m = Object.assign({}, this.opts);
        m.id = this.id;
        m.type = this.type;
        m.a = this.a._physiN.id;
        m.b = this.b ? this.b._physiN.id : null;
        if (m.localA) m.localA = Array.from(m.localA);
        if (m.localB) m.localB = Array.from(m.localB);
        return m;
      }
      /** False leaves the joint out of the solver. It does not take it away. */
      setEnabled(value) {
        this.opts.enabled = value;
        if (this.scene) this.scene.execute("setConstraintParams", { id: this.id, enabled: value });
      }
      /** Changes the rest length of a distance joint. */
      setRestLength(value) {
        this.opts.rest = value;
        if (this.scene) this.scene.execute("setConstraintParams", { id: this.id, rest: value });
      }
      /**
       * The two angles of a hinge, in radians. The angle wraps at pi, thus a
       * limit outside `(-pi, pi)` has no meaning. Give `-Infinity` and
       * `Infinity` to take the limits away.
       */
      setLimits(lowerAngle, upperAngle) {
        this.opts.lowerAngle = lowerAngle;
        this.opts.upperAngle = upperAngle;
        if (this.scene) this.scene.execute("setConstraintParams", { id: this.id, lowerAngle, upperAngle });
      }
      /**
       * The motor of a hinge. `speed` is the angular speed of `b` against `a`
       * in the free plane, and `maxTorque` is the largest torque that the
       * motor gives. A torque of 0 turns the motor off.
       *
       * The value is a TORQUE and not an impulse, thus the strength of the
       * motor does not change when you change `subSteps`.
       */
      setMotor(speed, maxTorque) {
        this.opts.motorSpeed = speed;
        this.opts.maxMotorTorque = maxTorque;
        if (this.scene) {
          this.scene.execute("setConstraintParams", { id: this.id, motorSpeed: speed, maxMotorTorque: maxTorque });
        }
      }
      /**
       * Makes the joint a spring of `hertz` and of that damping ratio. A
       * `hertz` of 0 makes the joint rigid again. The stiffness does not
       * change with the length of the step. See README, section 12.
       */
      setSoftness(hertz, damping) {
        this.opts.hertz = hertz;
        this.opts.damping = damping;
        if (this.scene) this.scene.execute("setConstraintParams", { id: this.id, hertz, damping });
      }
      /**
       * Breaks the joint when the load goes above `force` newtons or above
       * `torque` newton metres. A value of 0 means that the joint never breaks.
       *
       * A joint that breaks sends the event `broken` and it goes out of the
       * scene. The two bodies can then touch each other, and they usually
       * overlap at that moment.
       */
      setBreak(force, torque) {
        this.opts.breakForce = force;
        this.opts.breakTorque = torque;
        if (this.scene) {
          this.scene.execute("setConstraintParams", { id: this.id, breakForce: force, breakTorque: torque });
        }
      }
    };
    Eventable.make(PhysiN.Constraint);
    /**
     * Holds one point of `a` on one point of `b`. This is the joint that
     * attaches two solids. It takes away `n` degrees of freedom, and it leaves
     * every rotation free.
     *
     *   scene.add(post); scene.add(bob);
     *   scene.addConstraint(new PhysiN.PointJoint(post, bob, [0, 0, 0], [-2, 0, 0]));
     */
    PhysiN.PointJoint = class PhysiNPointJoint extends PhysiN.Constraint {
      constructor(a, b, localA, localB, opts = {}) {
        super("point", a, b, Object.assign({}, opts, { localA, localB }));
      }
    };
    /**
     * Holds the length between one point of `a` and one point of `b`. It takes
     * away 1 degree of freedom.
     *
     * `opts.rest` is the length. The default is the length that the two points
     * have when the joint starts. `opts.mode` is "rod" (the default, it holds
     * the length), "rope" (it only stops the length from growing) or "strut"
     * (it only stops the length from falling).
     */
    PhysiN.DistanceJoint = class PhysiNDistanceJoint extends PhysiN.Constraint {
      constructor(a, b, localA, localB, opts = {}) {
        super("distance", a, b, Object.assign({}, opts, { localA, localB }));
      }
    };
    /**
     * Welds `a` to `b`. It takes away `n + k` degrees of freedom, thus 6 in 3
     * dimensions and 10 in 4 dimensions. The two bodies keep the relation that
     * they have when the joint starts.
     */
    PhysiN.FixedJoint = class PhysiNFixedJoint extends PhysiN.Constraint {
      constructor(a, b, localA, localB, opts = {}) {
        super("fixed", a, b, Object.assign({}, opts, { localA, localB }));
      }
    };
    /**
     * Holds one point, and leaves ONE ROTATION PLANE free. It takes away
     * `n + k - 1` degrees of freedom: 5 of 6 in 3 dimensions, and 9 of 10 in 4
     * dimensions.
     *
     * `opts.plane` names the free plane in the frame of `a`: a number of the
     * lexicographic order of the planes, or a bivector of `k` components. In 4
     * dimensions the plane 0 is `(x y)` and the plane 2 is `(x w)`. See README,
     * section 9. The plane must be SIMPLE: `e_xy + e_zw` holds two planes and
     * it is not a hinge, thus the engine refuses it.
     *
     * `opts.lowerAngle` and `opts.upperAngle` give the limits, in radians.
     * `opts.motorSpeed` and `opts.maxMotorTorque` give the motor. Use
     * `setLimits` and `setMotor` to change them while the world runs.
     */
    PhysiN.HingeJoint = class PhysiNHingeJoint extends PhysiN.Constraint {
      constructor(a, b, localA, localB, opts = {}) {
        super("hinge", a, b, Object.assign({}, opts, { localA, localB }));
      }
    };
    /**
     * Holds a body in a subspace of the position, of the rotation, or of the
     * two. With a null `b` the subspace is a subspace of the world.
     *
     * Give `opts.lockAxes` or `opts.freeAxes` for the position, and
     * `opts.lockPlanes` or `opts.freePlanes` for the rotation. A member of a
     * list is a number of an axis or of a plane, or a vector.
     *
     *   // hold a 4D body on the hyperplane w = 0, and leave x, y and z free
     *   scene.addConstraint(new PhysiN.SubspaceJoint(body, null, { lockAxes: [3] }));
     *
     * `opts.worldFrame` at true holds the directions in the world frame. The
     * default turns them with `a`.
     *
     * A lock of the rotation that leaves TWO OR MORE PLANES FREE holds the
     * angular velocity only, and it does not correct a drift. See
     * `angularError`.
     */
    PhysiN.SubspaceJoint = class PhysiNSubspaceJoint extends PhysiN.Constraint {
      constructor(a, b, opts = {}) {
        super("subspace", a, b, opts);
      }
    };
    return PhysiN;
  }

  // src/global.js
  //
  // The entry of the classic script. It puts `PhysiN` on the global object.
  // Load three.js first.
  if (typeof globalThis.THREE === "undefined") {
    console.error("PhysiN: load three.js before physiN.js");
  } else {
    globalThis.PhysiN = createPhysiN(globalThis.THREE);
  }
})();