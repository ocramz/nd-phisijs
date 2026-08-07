/**
 * test/lib/three-stub.js -- a small part of three.js, sufficient for the
 * plugin.
 *
 * `physin.js` builds its classes from `THREE.Mesh`, `THREE.Object3D` and the
 * other names below. It builds them one time, when the file loads. Thus a
 * test must give these names, but the tests do not draw anything. The stub
 * holds the data only.
 *
 * `test/worker.js` and `test/lib/load.js` use this file.
 */

/** A point of 3 dimensions. */
export class V3 {
  constructor(x = 0, y = 0, z = 0) { this.x = x; this.y = y; this.z = z; }
  set(x, y, z) { this.x = x; this.y = y; this.z = z; return this; }
  setScalar(s) { return this.set(s, s, s); }
}

/** A surface, as a table of attributes. */
export class Geom {
  constructor() { this.attributes = {}; this.drawRange = { start: 0, count: Infinity }; }
  setAttribute(k, a) { this.attributes[k] = a; return this; }
  getAttribute(k) { return this.attributes[k]; }
  setDrawRange(s, c) { this.drawRange = { start: s, count: c }; }
  computeBoundingBox() {
    this.boundingBox = { min: new V3(-0.5, -0.5, -0.5), max: new V3(0.5, 0.5, 0.5) };
  }
  computeBoundingSphere() { this.boundingSphere = { center: new V3(), radius: 1 }; }
}

/** A node of the tree of the scene. */
export class Obj3D {
  constructor() {
    this.position = new V3();
    this.quaternion = { x: 0, y: 0, z: 0, w: 1, set() { return this; }, setFromRotationMatrix() { return this; } };
    this.scale = new V3(1, 1, 1); this.children = []; this.visible = true;
    this.parent = null; this.userData = {};
  }
  add(o) { this.children.push(o); o.parent = this; return this; }
  remove(o) { const i = this.children.indexOf(o); if (i >= 0) this.children.splice(i, 1); }
}

/** The namespace that `createPhysiN` takes. */
export const THREE = {
  Vector3: V3,
  Quaternion: class { constructor() { this.w = 1; } set() { return this; } setFromRotationMatrix() { return this; } },
  Matrix4: class { constructor() { this.elements = new Float64Array(16); } set() { return this; } },
  Sphere: class { constructor() { this.center = new V3(); this.radius = 0; } },
  BufferAttribute: class { constructor(a, i) { this.array = a; this.itemSize = i; } },
  Float32BufferAttribute: class { constructor(a, i) { this.array = Float32Array.from(a); this.itemSize = i; } },
  BufferGeometry: Geom,
  SphereGeometry: class extends Geom {
    constructor(radius = 1) { super(); this.parameters = { radius }; }
    computeBoundingSphere() {
      this.boundingSphere = { center: new V3(), radius: this.parameters.radius };
    }
    computeBoundingBox() {
      const r = this.parameters.radius;
      this.boundingBox = { min: new V3(-r, -r, -r), max: new V3(r, r, r) };
    }
  },
  BoxGeometry: class extends Geom {},
  Object3D: Obj3D,
  Scene: class extends Obj3D {},
  Mesh: class extends Obj3D { constructor(g, m) { super(); this.geometry = g; this.material = m; } },
};
