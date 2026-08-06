/**
 * PhysiN 0.1.0 -- the physics engine, for a web worker.  (MIT)
 *
 * One file. It is a classic worker script. Point the plugin at it:
 *   PhysiN.scripts.worker = 'physiN_worker.js';
 *
 * It contains the same engine code as physiN.js. The plugin runs that code in
 * the main thread when PhysiN.scripts.worker is null. Thus the results are
 * the same in the two conditions.
 */
(() => {
  // src/nd/core/dims.js
  var cache = /* @__PURE__ */ new Map();
  function popcount(x) {
    let c = 0;
    while (x !== 0) {
      x &= x - 1;
      c += 1;
    }
    return c;
  }
  function factorial(m) {
    let f = 1;
    for (let i = 2; i <= m; i += 1) f *= i;
    return f;
  }
  function dims(n) {
    if (!Number.isInteger(n) || n < 2) throw new Error("dims: n must be an integer of 2 or more");
    if (cache.has(n)) return cache.get(n);
    const N = 1 << n;
    const grade = new Int8Array(N);
    const reverseSign = new Int8Array(N);
    for (let m = 0; m < N; m += 1) {
      const g = popcount(m);
      grade[m] = g;
      reverseSign[m] = g * (g - 1) / 2 & 1 ? -1 : 1;
    }
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
    const evenBlade = [];
    const slotOfBlade = new Int32Array(N).fill(-1);
    for (let m = 0; m < N; m += 1) {
      if ((grade[m] & 1) === 0) {
        slotOfBlade[m] = evenBlade.length;
        evenBlade.push(m);
      }
    }
    const r = evenBlade.length;
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

  // src/nd/algebra/multivector.js
  function mvZero(D) {
    return new Float64Array(D.N);
  }
  function mvFromVector(D, v, out) {
    const M = out || mvZero(D);
    M.fill(0);
    for (let i = 0; i < D.n; i += 1) M[D.vecBlade[i]] = v[i];
    return M;
  }
  function mvFromBivector(D, B, out) {
    const M = out || mvZero(D);
    M.fill(0);
    for (let p = 0; p < D.k; p += 1) M[D.biBlade[p]] = B[p];
    return M;
  }
  function mvFromRotor(D, R, out) {
    const M = out || mvZero(D);
    M.fill(0);
    for (let s = 0; s < D.r; s += 1) M[D.evenBlade[s]] = R[s];
    return M;
  }
  function mvToVector(D, M, out) {
    const v = out || new Float64Array(D.n);
    for (let i = 0; i < D.n; i += 1) v[i] = M[D.vecBlade[i]];
    return v;
  }
  function mvToBivector(D, M, out) {
    const B = out || new Float64Array(D.k);
    for (let p = 0; p < D.k; p += 1) B[p] = M[D.biBlade[p]];
    return B;
  }
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

  // src/nd/algebra/rotor.js
  function rotorIdentity(D, out) {
    const R = out || new Float64Array(D.r);
    R.fill(0);
    R[D.scalarSlot] = 1;
    return R;
  }
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
  function rotorReverse(D, R, out) {
    const C = out || new Float64Array(D.r);
    for (let s = 0; s < D.r; s += 1) C[s] = D.reverseSign[D.evenBlade[s]] * R[s];
    return C;
  }
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
  function rotorFromBivectorAngle(D, B, angle, out) {
    let mag = 0;
    for (let p = 0; p < D.k; p += 1) mag += B[p] * B[p];
    mag = Math.sqrt(mag);
    if (mag < 1e-300) return rotorIdentity(D, out);
    const S = new Float64Array(D.k);
    for (let p = 0; p < D.k; p += 1) S[p] = B[p] / mag * angle;
    return rotorExp(D, S, out);
  }
  var scratch = /* @__PURE__ */ new WeakMap();
  function pad(D) {
    let s = scratch.get(D);
    if (!s) {
      s = { a: mvZero(D), b: mvZero(D), c: mvZero(D), d: mvZero(D) };
      scratch.set(D, s);
    }
    return s;
  }
  function rotorApplyVector(D, R, x, out) {
    const s = pad(D);
    mvFromRotor(D, R, s.a);
    mvFromVector(D, x, s.b);
    mvGp(D, s.a, s.b, s.c);
    for (let m = 0; m < D.N; m += 1) s.b[m] = D.reverseSign[m] * s.a[m];
    mvGp(D, s.c, s.b, s.d);
    return mvToVector(D, s.d, out);
  }
  function rotorApplyBivector(D, R, B, out) {
    const s = pad(D);
    mvFromRotor(D, R, s.a);
    mvFromBivector(D, B, s.b);
    mvGp(D, s.a, s.b, s.c);
    for (let m = 0; m < D.N; m += 1) s.b[m] = D.reverseSign[m] * s.a[m];
    mvGp(D, s.c, s.b, s.d);
    return mvToBivector(D, s.d, out);
  }
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
  function rotorBetweenVectors(D, a, b, out) {
    const n = D.n;
    let dot = 0;
    for (let i = 0; i < n; i += 1) dot += a[i] * b[i];
    if (dot > 1 - 1e-12) return rotorIdentity(D, out);
    if (dot < -1 + 1e-12) {
      let best = 0;
      for (let i = 1; i < n; i += 1) if (Math.abs(a[i]) < Math.abs(a[best])) best = i;
      const t = new Float64Array(n);
      t[best] = 1;
      let d2 = 0;
      for (let i = 0; i < n; i += 1) d2 += t[i] * a[i];
      for (let i = 0; i < n; i += 1) t[i] -= d2 * a[i];
      let ln = 0;
      for (let i = 0; i < n; i += 1) ln += t[i] * t[i];
      ln = Math.sqrt(ln);
      for (let i = 0; i < n; i += 1) t[i] /= ln;
      const B = new Float64Array(D.k);
      for (let p = 0; p < D.k; p += 1) {
        const [i, j] = D.pairs[p];
        B[p] = a[i] * t[j] - a[j] * t[i];
      }
      return rotorFromBivectorAngle(D, B, Math.PI, out);
    }
    const R = out || new Float64Array(D.r);
    R.fill(0);
    R[D.scalarSlot] = 1 + dot;
    for (let p = 0; p < D.k; p += 1) {
      const [i, j] = D.pairs[p];
      R[D.biSlot[p]] = b[i] * a[j] - b[j] * a[i];
    }
    const f = 1 / Math.sqrt(2 * (1 + dot));
    for (let i = 0; i < D.r; i += 1) R[i] *= f;
    return R;
  }
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
      if (ln < 1e-12) {
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
    for (let c = 0; c < n - 1; c += 1) {
      for (let i = 0; i < n; i += 1) {
        p[i] = cur[i * n + c];
        tgt[i] = F[i * n + c];
      }
      let dot2 = 0;
      for (let i = 0; i < n; i += 1) dot2 += p[i] * tgt[i];
      if (dot2 > 1 - 1e-15) continue;
      const S = rotorBetweenVectors(D, p, tgt);
      acc = rotorMul(D, S, acc);
      let ln = 0;
      for (let i = 0; i < n; i += 1) {
        q[i] = tgt[i] - dot2 * p[i];
        ln += q[i] * q[i];
      }
      ln = Math.sqrt(ln);
      if (ln < 1e-15) continue;
      for (let i = 0; i < n; i += 1) q[i] /= ln;
      const cs = dot2;
      const sn = ln;
      for (let col = c; col < n; col += 1) {
        let a = 0;
        let b = 0;
        for (let i = 0; i < n; i += 1) {
          a += cur[i * n + col] * p[i];
          b += cur[i * n + col] * q[i];
        }
        const na = a * cs - b * sn;
        const nb = a * sn + b * cs;
        for (let i = 0; i < n; i += 1) {
          cur[i * n + col] += (na - a) * p[i] + (nb - b) * q[i];
        }
      }
    }
    let dot = 0;
    for (let i = 0; i < D.r; i += 1) dot += acc[i] * R[i];
    const C = out || new Float64Array(D.r);
    const s = dot < 0 ? -1 : 1;
    for (let i = 0; i < D.r; i += 1) C[i] = acc[i] * s;
    return C;
  }
  function rotorDefect(D, R) {
    const Rr = rotorReverse(D, R);
    const P = rotorMul(D, R, Rr);
    let e = Math.abs(P[D.scalarSlot] - 1);
    for (let i = 0; i < D.r; i += 1) if (i !== D.scalarSlot) e += Math.abs(P[i]);
    return e;
  }

  // src/nd/core/linalg.js
  function matZero(rows, cols) {
    return new Float64Array(rows * cols);
  }
  function matIdentity(m, out) {
    const A = out || matZero(m, m);
    A.fill(0);
    for (let i = 0; i < m; i += 1) A[i * m + i] = 1;
    return A;
  }
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
  function matVec(A, x, rows, cols, out) {
    const y = out || new Float64Array(rows);
    for (let i = 0; i < rows; i += 1) {
      let s = 0;
      for (let j = 0; j < cols; j += 1) s += A[i * cols + j] * x[j];
      y[i] = s;
    }
    return y;
  }
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

  // src/nd/algebra/star.js
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
  function wedgeVec(D, a, b, out) {
    const B = out || new Float64Array(D.k);
    for (let p = 0; p < D.k; p += 1) {
      const [i, j] = D.pairs[p];
      B[p] = a[i] * b[j] - a[j] * b[i];
    }
    return B;
  }
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
  var nextId = 1;
  var Body = class {
    /**
     * @param {object} D dimension tables
     * @param {object} opts shape, mass, position, rotor, material
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
      this._tmpN = new Float64Array(n);
      this._tmpK = new Float64Array(k);
      this._star = new Float64Array(k * n);
      this.updateDerived();
      if (opts.angularVelocity) this.setAngularVelocity(opts.angularVelocity);
    }
    /** Build `[R]2`, the world inverse inertia and the angular velocity. */
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
    /** World inertia tensor. Use it only for a test or for the energy. */
    inertiaWorld() {
      const { k } = this.D;
      const T = matMul(this.R2, this.inertia, k, k, k);
      return matMulT(T, this.R2, k, k, k);
    }
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
    localToWorld(a, out) {
      const p = this.localToWorldDir(a, out);
      for (let i = 0; i < this.D.n; i += 1) p[i] += this.x[i];
      return p;
    }
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
    /** Apply an impulse `j` at the world offset `r`. See item A11. */
    applyImpulse(j, r) {
      if (this.isStatic) return;
      const { n, k } = this.D;
      for (let i = 0; i < n; i += 1) this.v[i] += this.invMass * j[i] * this.linearFactor[i];
      const S = starMatrix(this.D, r, this._star);
      const dL = matVec(S, j, k, n, this._tmpK);
      for (let p = 0; p < k; p += 1) this.L[p] += dL[p] * this.angularFactor[p];
      matVec(this.invInertiaWorld, this.L, k, k, this.w);
      this.wake();
    }
    applyCentralImpulse(j) {
      if (this.isStatic) return;
      for (let i = 0; i < this.D.n; i += 1) this.v[i] += this.invMass * j[i] * this.linearFactor[i];
      this.wake();
    }
    applyTorqueImpulse(dL) {
      if (this.isStatic) return;
      const { k } = this.D;
      for (let p = 0; p < k; p += 1) this.L[p] += dL[p] * this.angularFactor[p];
      matVec(this.invInertiaWorld, this.L, k, k, this.w);
      this.wake();
    }
    applyCentralForce(f) {
      for (let i = 0; i < this.D.n; i += 1) this.force[i] += f[i];
      this.wake();
    }
    applyForce(f, r) {
      const { n, k } = this.D;
      for (let i = 0; i < n; i += 1) this.force[i] += f[i];
      const S = starMatrix(this.D, r, this._star);
      const t = matVec(S, f, k, n, this._tmpK);
      for (let p = 0; p < k; p += 1) this.torque[p] += t[p];
      this.wake();
    }
    applyTorque(t) {
      for (let p = 0; p < this.D.k; p += 1) this.torque[p] += t[p];
      this.wake();
    }
    clearForces() {
      this.force.fill(0);
      this.torque.fill(0);
    }
    wake() {
      if (this.sleeping) {
        this.sleeping = false;
      }
      this.sleepTimer = 0;
    }
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
    /** Axis aligned box in the world frame. It gives `min` and `max`. */
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

  // src/nd/detect/nearest.js
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
  var tableCache = /* @__PURE__ */ new WeakMap();
  function combinations(list, m) {
    if (m === 0) return [[]];
    if (m > list.length) return [];
    const out = [];
    for (let i = 0; i <= list.length - m; i += 1) {
      for (const rest of combinations(list.slice(i + 1), m - 1)) out.push([list[i]].concat(rest));
    }
    return out;
  }
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
  var normalCache = /* @__PURE__ */ new WeakMap();
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
  var warned = false;
  function warnNoPair(ta, tb) {
    if (warned) return;
    warned = true;
    if (typeof console !== "undefined" && console.warn) {
      console.warn(`PhysiN: no collision test for the pair (${ta}, ${tb}). A torus and a convex mesh touch a half space and a hypersphere only.`);
    }
  }
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
  var vertexCache = /* @__PURE__ */ new WeakMap();
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

  // src/nd/resolve/solver.js
  var scratchCache = /* @__PURE__ */ new WeakMap();
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
  function relativeVelocity(D, c, out) {
    const { n } = D;
    const s = scratch2(D);
    const ua = c.a.pointVelocity(c.rA, s.ua);
    const ub = c.b.pointVelocity(c.rB, s.ub);
    const u = out || s.u;
    for (let i = 0; i < n; i += 1) u[i] = ub[i] - ua[i];
    return u;
  }
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
  function applyPair(D, c, j) {
    const { n } = D;
    if (!c.aStatic) {
      const neg = scratch2(D).neg;
      for (let i = 0; i < n; i += 1) neg[i] = -j[i];
      c.a.applyImpulse(neg, c.rA);
    }
    if (!c.bStatic) c.b.applyImpulse(j, c.rB);
  }
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
  var wantedCache = /* @__PURE__ */ new WeakMap();
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

  // src/nd/integrate/integrator.js
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
    if (opts.gyroscopic !== false) applyGyroscopic(D, body, dt, opts);
  }
  function applyGyroscopic(D, body, dt, opts) {
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
    matVec(body.R2, w2, k, k, body.w);
    const Ib = matVec(I, w2, k, k);
    matVec(body.R2, Ib, k, k, body.L);
  }
  function integratePositions(D, body, dt, opts) {
    if (body.isStatic || body.sleeping) return;
    const { n, r } = D;
    for (let i = 0; i < n; i += 1) body.x[i] += dt * body.v[i];
    const W = new Float64Array(r);
    for (let p = 0; p < D.k; p += 1) W[D.biSlot[p]] = body.w[p];
    const dR = rotorMul(D, W, body.R);
    for (let i = 0; i < r; i += 1) body.R[i] += -0.5 * dt * dR[i];
    const tol = opts.rotorTolerance !== void 0 ? opts.rotorTolerance : 1e-9;
    if (rotorDefect(D, body.R) > tol) {
      const C = rotorCorrect(D, body.R);
      body.R.set(C);
    }
    body.updateDerived();
  }

  // src/nd/world.js
  var defaultParams = {
    // time
    fixedTimeStep: 1 / 60,
    subSteps: 1,
    // solver
    iterations: 10,
    shockIterations: 2,
    useShockPropagation: true,
    useWarmStart: true,
    // contact
    penetrationSlop: 5e-3,
    biasFactor: 0.2,
    contactMargin: 0.02,
    restitutionThreshold: 0.5,
    maxContacts: 0,
    // zero means 2^(n-1)
    // rotor
    rotorTolerance: 1e-9,
    // gyroscopic term
    gyroscopic: true,
    gyroscopicIterations: 1,
    // sleep
    allowSleep: true,
    sleepLinearVelocity: 0.03,
    sleepAngularVelocity: 0.03,
    sleepTime: 0.6
  };
  var World = class {
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
      this.manifolds = /* @__PURE__ */ new Map();
      this.time = 0;
      this.listeners = { collision: [] };
    }
    addBody(body) {
      this.bodies.push(body);
      return body;
    }
    removeBody(body) {
      const i = this.bodies.indexOf(body);
      if (i >= 0) this.bodies.splice(i, 1);
    }
    createBody(opts) {
      return this.addBody(new Body(this.D, opts));
    }
    setGravity(g) {
      this.gravity.set(g);
      for (const b of this.bodies) b.wake();
    }
    on(name, fn) {
      (this.listeners[name] = this.listeners[name] || []).push(fn);
    }
    emit(name, ...args) {
      for (const fn of this.listeners[name] || []) fn(...args);
    }
    /** Sort and sweep on axis 0. It gives the pairs whose boxes overlap. */
    broadPhase() {
      const { n } = this;
      const margin = this.params.contactMargin;
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
    narrowPhase(pairs) {
      const out = [];
      for (const [a, b] of pairs) collide(this.D, a, b, out, this.params);
      return out;
    }
    /** Copy the impulses of the last step onto the new contacts. */
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
    step(dt) {
      const { D, params } = this;
      const sub = Math.max(1, params.subSteps | 0);
      const h = dt / sub;
      for (let s = 0; s < sub; s += 1) this.subStep(h);
      for (const b of this.bodies) b.clearForces();
    }
    subStep(dt) {
      const { D, params } = this;
      for (const b of this.bodies) integrateVelocities(D, b, dt, this.gravity, params);
      const pairs = this.broadPhase();
      const contacts = this.narrowPhase(pairs);
      this.contacts = contacts;
      for (const c of contacts) {
        if (c.depth > 0) {
          if (!c.a.isStatic && c.b.sleeping === false) c.a.wake();
          if (!c.b.isStatic && c.a.sleeping === false) c.b.wake();
          this.emit("collision", c.a, c.b, c);
        }
      }
      if (params.useWarmStart) this.applyWarmStartCache(contacts);
      for (const c of contacts) prepareContact(D, c, dt, params);
      if (params.useWarmStart) {
        for (const c of contacts) {
          if (c.oldTangent && c.oldTangent.length === c.tangentImpulse.length) {
            c.tangentImpulse.set(c.oldTangent);
          }
          warmStart(D, c);
        }
      }
      for (let it = 0; it < params.iterations; it += 1) {
        if ((it & 1) === 0) {
          for (let i = 0; i < contacts.length; i += 1) solveContact(D, contacts[i]);
        } else {
          for (let i = contacts.length - 1; i >= 0; i -= 1) solveContact(D, contacts[i]);
        }
      }
      if (params.useShockPropagation && contacts.length > 0) {
        buildContactGraph(this.bodies, contacts);
        shockPropagation(D, contacts, dt, params);
      }
      for (const b of this.bodies) integratePositions(D, b, dt, params);
      if (params.allowSleep) this.updateSleep(dt);
      this.time += dt;
    }
    updateSleep(dt) {
      const { n, k } = this.D;
      const lv = this.params.sleepLinearVelocity ** 2;
      const av = this.params.sleepAngularVelocity ** 2;
      for (const b of this.bodies) {
        if (b.isStatic || !b.allowSleep) continue;
        let s = 0;
        for (let i = 0; i < n; i += 1) s += b.v[i] * b.v[i];
        let a = 0;
        for (let p = 0; p < k; p += 1) a += b.w[p] * b.w[p];
        if (s < lv && a < av) {
          b.sleepTimer += dt;
          if (b.sleepTimer > this.params.sleepTime) {
            b.sleeping = true;
            b.v.fill(0);
            b.L.fill(0);
            b.w.fill(0);
          }
        } else {
          b.sleepTimer = 0;
          b.sleeping = false;
        }
      }
    }
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

  // src/nd/body/massprops.js
  var productCache = /* @__PURE__ */ new WeakMap();
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
  function permutations(list) {
    if (list.length <= 1) return [list.slice()];
    const out = [];
    for (let i = 0; i < list.length; i += 1) {
      const rest = list.slice(0, i).concat(list.slice(i + 1));
      for (const p of permutations(rest)) out.push([list[i]].concat(p));
    }
    return out;
  }
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
  function hyperBoxInertia(D, h, mass) {
    const { n, k } = D;
    const I = matZero(k, k);
    for (let p = 0; p < k; p += 1) {
      const [i, j] = D.pairs[p];
      I[p * k + p] = mass / 3 * (h[i] * h[i] + h[j] * h[j]);
    }
    return I;
  }
  function hyperSphereInertia(D, radius, mass) {
    const { n, k } = D;
    const I = matZero(k, k);
    const v = 2 * mass * radius * radius / (n + 2);
    for (let p = 0; p < k; p += 1) I[p * k + p] = v;
    return I;
  }
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
  function torusVolume(n, R, r) {
    return 2 * Math.PI * R * ballVolume(n - 1, r);
  }
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
  function torusInertia(D, R, r, plane, mass) {
    return inertiaFromCovariance(D, torusCovariance(D, R, r, plane, mass));
  }
  function torusMesh(D, R, r, plane, segments) {
    const { n } = D;
    if (n === 3) return torusMesh3(D, R, r, plane, segments || [32, 16]);
    if (n === 4) return torusMesh4(D, R, r, plane, segments || [16, 8, 12]);
    throw new Error("torusMesh: only n = 3 and n = 4 have a mesh");
  }
  function otherAxes(n, plane) {
    const out = [];
    for (let i = 0; i < n; i += 1) if (i !== plane[0] && i !== plane[1]) out.push(i);
    return out;
  }
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
  var KUHN3 = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0]
  ];
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

  // src/nd/workerCore.js
  var MESSAGE_TYPES = { WORLDREPORT: 0, COLLISIONREPORT: 1 };
  function createEngine(post) {
    let world = null;
    let D = null;
    let stride = 0;
    let contactStride = 0;
    let worldReport = null;
    let collisionReport = null;
    const bodies = /* @__PURE__ */ new Map();
    const collisions = [];
    let fixedTimeStep = 1 / 60;
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
    const commands = {
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
        }
      },
      updateTransform(params) {
        const b = bodies.get(params.id);
        if (!b) return;
        if (params.position) b.x.set(params.position);
        if (params.rotor) b.R.set(params.rotor);
        b.updateDerived();
        b.wake();
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
      setMass(params) {
        const b = bodies.get(params.id);
        if (!b) return;
        const opts = { shape: b.shape, mass: params.value };
        const nb = new Body(D, Object.assign({ id: b.id }, opts));
        nb.x.set(b.x);
        nb.R.set(b.R);
        nb.v.set(b.v);
        nb.updateDerived();
        world.removeBody(b);
        world.addBody(nb);
        bodies.set(b.id, nb);
      },
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

  // src/physiN_worker.js
  var engine = createEngine((data) => {
    if (data && data.buffer && data.BYTES_PER_ELEMENT === 4) {
      const copy = new Float32Array(data);
      self.postMessage(copy, [copy.buffer]);
    } else {
      self.postMessage(data);
    }
  });
  self.onmessage = (event) => engine.handle(event.data);
})();