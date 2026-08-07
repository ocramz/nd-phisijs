/**
 * test/worker.js -- tests the two files that go to a page:
 *   physin.js         the three.js plugin, a classic script
 *   physin_worker.js  the engine, a classic worker script
 *
 * The test loads the two bundles, and it connects them with a Worker that is
 * a stub. Thus it tests the message protocol from end to end.
 */
import fs from 'fs';
import { THREE, Geom } from './lib/three-stub.js';

let pass = 0; let fail = 0;
function check(name, ok, detail = '') {
  if (ok) { pass += 1; console.log(`  ok    ${name}${detail ? `  (${detail})` : ''}`); } else {
    fail += 1; console.log(`  FAIL  ${name}  ${detail}`);
  }
}

/** Keeps one bad section from stopping the file. */
function section(name, fn) {
  try {
    fn();
  } catch (e) {
    fail += 1; console.log(`  FAIL  ${name}  ${e && e.message}`);
  }
}

// The THREE stub is in test/lib/three-stub.js, and the specs use it too.

// ------------------------------------------------- the worker, as a stub
const workerCode = fs.readFileSync(new URL('../physin_worker.js', import.meta.url), 'utf8');
const mainCode = fs.readFileSync(new URL('../physin.js', import.meta.url), 'utf8');

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
// eslint-disable-next-line no-new-func
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
section('the physics in the main thread', () => {
  PhysiN.scripts.worker = null;
  const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
  scene.add(new PhysiN.PlaneMesh(new Geom(), {}, [0, 1, 0], 0));
  const box = new PhysiN.BoxMesh(new Geom(), {}, 1);
  box.setPositionN([0, 4, 0]);
  scene.add(box);
  for (let s = 0; s < 500; s += 1) scene.simulate(1 / 120, 1);
  check('the main thread makes the box rest on the ground',
    Math.abs(box.position.y - 0.5) < 0.02, `y = ${box.position.y.toFixed(4)}`);
});

// ---- 2. the same with the worker file
section('the physics in the worker', () => {
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
    Math.abs(ball.position.y - 0.4) < 0.02, `y = ${ball.position.y.toFixed(4)}`);
  check('the two files exchange messages',
    messagesToWorker > 500 && messagesFromWorker > 500,
    `${messagesToWorker} to the worker, ${messagesFromWorker} from it`);
});

// ---- 3. four dimensions through the worker
section('four dimensions through the worker', () => {
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
});

console.log(`\n${pass} tests pass, ${fail} tests fail`);
if (fail > 0) process.exitCode = 1;