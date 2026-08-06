/**
 * test/worker.js -- tests the two files that go to a page:
 *   dist/physiN.js         the three.js plugin, a classic script
 *   dist/physiN_worker.js  the engine, a classic worker script
 *
 * The test loads the two bundles, and it connects them with a Worker that is
 * a stub. Thus it tests the message protocol from end to end.
 */
import fs from 'fs';

let pass = 0; let fail = 0;
function check(name, ok, detail = '') {
  if (ok) { pass += 1; console.log(`  ok    ${name}${detail ? `  (${detail})` : ''}`); } else {
    fail += 1; console.log(`  FAIL  ${name}  ${detail}`);
  }
}

// ------------------------------------------------------------ THREE stub
class V3 {
  constructor(x = 0, y = 0, z = 0) { this.x = x; this.y = y; this.z = z; }
  set(x, y, z) { this.x = x; this.y = y; this.z = z; return this; }
  setScalar(s) { return this.set(s, s, s); }
}
class Geom {
  constructor() { this.attributes = {}; this.drawRange = { start: 0, count: Infinity }; }
  setAttribute(k, a) { this.attributes[k] = a; return this; }
  getAttribute(k) { return this.attributes[k]; }
  setDrawRange(s, c) { this.drawRange = { start: s, count: c }; }
  computeBoundingBox() {
    this.boundingBox = { min: new V3(-0.5, -0.5, -0.5), max: new V3(0.5, 0.5, 0.5) };
  }
  computeBoundingSphere() { this.boundingSphere = { center: new V3(), radius: 1 }; }
}
class Obj3D {
  constructor() {
    this.position = new V3();
    this.quaternion = { x: 0, y: 0, z: 0, w: 1, set() { return this; }, setFromRotationMatrix() { return this; } };
    this.scale = new V3(1, 1, 1); this.children = []; this.visible = true;
    this.parent = null; this.userData = {};
  }
  add(o) { this.children.push(o); o.parent = this; return this; }
  remove(o) { const i = this.children.indexOf(o); if (i >= 0) this.children.splice(i, 1); }
}
const THREE = {
  Vector3: V3,
  Quaternion: class { constructor() { this.w = 1; } set() { return this; } setFromRotationMatrix() { return this; } },
  Matrix4: class { constructor() { this.elements = new Float64Array(16); } set() { return this; } },
  Sphere: class { constructor() { this.center = new V3(); this.radius = 0; } },
  BufferAttribute: class { constructor(a, i) { this.array = a; this.itemSize = i; } },
  Float32BufferAttribute: class { constructor(a, i) { this.array = Float32Array.from(a); this.itemSize = i; } },
  BufferGeometry: Geom,
  SphereGeometry: class extends Geom {},
  BoxGeometry: class extends Geom {},
  Object3D: Obj3D,
  Scene: class extends Obj3D {},
  Mesh: class extends Obj3D { constructor(g, m) { super(); this.geometry = g; this.material = m; } },
};

// ------------------------------------------------- the worker, as a stub
const workerCode = fs.readFileSync(new URL('../dist/physiN_worker.js', import.meta.url), 'utf8');
const mainCode = fs.readFileSync(new URL('../dist/physiN.js', import.meta.url), 'utf8');

let messagesToWorker = 0;
let messagesFromWorker = 0;

/** One instance of the worker script, with its own `self`. */
function makeWorkerScript(onPost) {
  const self = {
    onmessage: null,
    postMessage(data) { messagesFromWorker += 1; onPost(data); },
  };
  // eslint-disable-next-line no-new-func
  new Function('self', workerCode)(self);
  return self;
}

class FakeWorker {
  constructor(url, opts) {
    this.url = url;
    this.type = opts && opts.type;
    this.onmessage = null;
    this._self = makeWorkerScript((data) => {
      // A real worker sends the message later. Send it now, so that the test
      // stays simple; the order of the messages is the same.
      if (this.onmessage) this.onmessage({ data });
    });
  }
  postMessage(data) {
    messagesToWorker += 1;
    if (this._self.onmessage) this._self.onmessage({ data });
  }
}

// ---------------------------------------------------------- load the two
global.THREE = THREE;
global.Worker = FakeWorker;
// eslint-disable-next-line no-new-func
new Function('THREE', 'Worker', 'globalThis', `${mainCode}\nreturn globalThis.PhysiN;`);
const load = new Function('THREE', 'Worker', `
  var globalThis = { THREE: THREE };
  ${mainCode}
  return globalThis.PhysiN;
`);
const PhysiN = load(THREE, FakeWorker);

console.log('\nThe two consolidated files');
check('physiN.js gives a global PhysiN', !!PhysiN && !!PhysiN.Scene);
check('it has the shapes of the physi.js interface',
  !!(PhysiN.BoxMesh && PhysiN.SphereMesh && PhysiN.ConvexMesh && PhysiN.PlaneMesh));
check('it has the shapes of four dimensions and more',
  !!(PhysiN.HyperBoxMesh && PhysiN.HyperSphereMesh && PhysiN.HyperTorusMesh
    && PhysiN.HyperMesh && PhysiN.HyperPlaneMesh));
check('the worker script is null at the start', PhysiN.scripts.worker === null);

// ---- 1. the physics in the main thread
{
  PhysiN.scripts.worker = null;
  const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
  scene.add(new PhysiN.PlaneMesh(new Geom(), {}, [0, 1, 0], 0));
  const box = new PhysiN.BoxMesh(new Geom(), {}, 1);
  box.setPositionN([0, 4, 0]);
  scene.add(box);
  for (let s = 0; s < 500; s += 1) scene.simulate(1 / 120, 1);
  check('the main thread makes the box rest on the ground',
    Math.abs(box.position.y - 0.5) < 0.02, `y = ${box.position.y.toFixed(4)}`);
}

// ---- 2. the same with the worker file
{
  PhysiN.scripts.worker = 'physiN_worker.js';
  const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
  check('the plugin starts a classic worker', scene._worker.type === 'classic',
    String(scene._worker.type));
  check('the plugin uses the file of PhysiN.scripts.worker',
    scene._worker.url === 'physiN_worker.js');
  scene.add(new PhysiN.PlaneMesh(new Geom(), {}, [0, 1, 0], 0));
  const box = new PhysiN.BoxMesh(new Geom(), {}, 1);
  box.setPositionN([0, 4, 0]);
  scene.add(box);
  const ball = new PhysiN.SphereMesh(new THREE.SphereGeometry(0.4), {}, 1);
  ball.setPositionN([2, 4, 0]);
  scene.add(ball);
  for (let s = 0; s < 500; s += 1) scene.simulate(1 / 120, 1);
  check('the worker makes the box rest on the ground',
    Math.abs(box.position.y - 0.5) < 0.02, `y = ${box.position.y.toFixed(4)}`);
  check('the worker sends the state of each body',
    Math.abs(ball.position.y - 1) < 0.02 || Math.abs(ball.position.y - 0.4) < 0.02,
    `y = ${ball.position.y.toFixed(4)}`);
  check('the two files exchange messages',
    messagesToWorker > 500 && messagesFromWorker > 500,
    `${messagesToWorker} to the worker, ${messagesFromWorker} from it`);
}

// ---- 3. four dimensions through the worker
{
  PhysiN.scripts.worker = 'physiN_worker.js';
  const scene = new PhysiN.Scene({ dimensions: 4, gravity: [0, -9.81, 0, 0] });
  scene.add(new PhysiN.HyperPlaneMesh([0, 1, 0, 0], 0, {}, new Geom()));
  const cube = new PhysiN.HyperBoxMesh([0.5, 0.5, 0.5, 0.5], {}, 1);
  cube.setPositionN([0, 3, 0, 0]);
  scene.add(cube);
  const ring = new PhysiN.HyperTorusMesh(1.2, 0.3, [0, 2], {}, 1, { mesh: false });
  ring.setPositionN([3, 3, 0, 0]);
  scene.add(ring);
  scene.sliceW = 0;
  for (let s = 0; s < 600; s += 1) scene.simulate(1 / 120, 1);
  check('the worker makes a tesseract rest on the 4D ground',
    Math.abs(cube.getPositionN()[1] - 0.5) < 0.02,
    `y = ${cube.getPositionN()[1].toFixed(4)}`);
  check('the worker makes a 4D torus rest on the 4D ground',
    Math.abs(ring.getPositionN()[1] - 0.3) < 0.02,
    `y = ${ring.getPositionN()[1].toFixed(4)}`);
  check('the slice of the torus gives triangles',
    ring.visible && ring.geometry.drawRange.count > 0,
    `${ring.geometry.drawRange.count} vertices`);
}

console.log(`\n${pass} tests pass, ${fail} tests fail`);
if (fail > 0) process.exitCode = 1;