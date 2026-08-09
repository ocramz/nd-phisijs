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

// ---- 4. the two bundles hold the same engine
section('the two bundles hold the same engine', () => {
  // `physin.js` and `physin_worker.js` must hold the same engine modules, word
  // for word. A change made on one side only is the classic failure of this
  // repository, and it is silent. This test makes it an error that a person
  // sees at once.
  //
  // Four algebra modules are NOT in this list. The worker bundle leaves out the
  // functions that the worker never calls, thus `linalg`, `multivector`,
  // `rotor` and `star` are shorter there. Every function that stays is the same
  // word for word.
  const shared = [
    'src/nd/core/dims.js',
    'src/nd/body/body.js',
    'src/nd/body/massprops.js',
    'src/nd/body/shapes.js',
    'src/nd/detect/nearest.js',
    'src/nd/detect/collide.js',
    'src/nd/integrate/integrator.js',
    'src/nd/resolve/solver.js',
    'src/nd/resolve/constraint.js',
    'src/nd/world.js',
    'src/nd/workerCore.js',
  ];

  /** The text of each module of a bundle, by its `// src/...` marker. */
  const modulesOf = (code) => {
    const lines = code.split('\n');
    const marks = [];
    for (let i = 0; i < lines.length; i += 1) {
      if (/^\s*\/\/ src\/\S+\.js$/.test(lines[i])) marks.push([i, lines[i].trim().slice(3)]);
    }
    const out = {};
    for (let m = 0; m < marks.length; m += 1) {
      const end = m + 1 < marks.length ? marks[m + 1][0] : lines.length;
      out[marks[m][1]] = lines.slice(marks[m][0], end).join('\n');
    }
    return out;
  };

  const a = modulesOf(mainCode);
  const b = modulesOf(workerCode);
  const missing = shared.filter((m) => !a[m] || !b[m]);
  const differ = shared.filter((m) => a[m] && b[m] && a[m] !== b[m]);
  check('every shared engine module is in the two bundles',
    missing.length === 0, missing.join(', '));
  check('every shared engine module is the same in the two bundles',
    differ.length === 0, differ.length ? `these differ: ${differ.join(', ')}` : `${shared.length} modules`);
});

// ---- 5. a joint through the worker
section('a joint through the worker', () => {
  PhysiN.scripts.worker = 'physiN_worker.js';
  const scene = new PhysiN.Scene({ dimensions: 3, gravity: [0, -9.81, 0] });
  // The geometry stub always gives a box of the size 1, thus the post and the
  // bob are unit cubes. The post has the mass 0, thus it is static.
  const post = new PhysiN.BoxMesh(new Geom(), {}, 0);
  post.setPositionN([0, 5, 0]);
  scene.add(post);
  const bob = new PhysiN.BoxMesh(new Geom(), {}, 1);
  bob.setPositionN([2, 5, 0]);
  scene.add(bob);
  // The anchor of the bob is 2 units to its left, thus it holds the center of
  // the post. The bob then hangs 2 units below the post, and it swings.
  scene.addConstraint(new PhysiN.PointJoint(post, bob, [0, 0, 0], [-2, 0, 0]));
  // The bob swings, thus its height rises and falls. Take the LOWEST point of
  // the whole run, and the WORST length. One reading at the end would only
  // give one point of that wave, and the bob can be at the top there.
  let lowest = Infinity;
  let worstLength = 0;
  for (let s = 0; s < 600; s += 1) {
    scene.simulate(1 / 120, 1);
    const q = bob.getPositionN();
    lowest = Math.min(lowest, q[1]);
    worstLength = Math.max(worstLength, Math.abs(Math.hypot(q[0], q[1] - 5, q[2]) - 2));
  }
  const p = bob.getPositionN();
  check('the point joint holds the length through the worker',
    worstLength < 0.02, `the largest error of the length = ${worstLength.toFixed(4)}`);
  check('the pendulum swings down', lowest < 3.5, `the lowest y = ${lowest.toFixed(4)}`);
  check('the joint does not blow up',
    Array.from(p).every((v) => Number.isFinite(v) && Math.abs(v) < 10),
    `position = ${Array.from(p).map((v) => v.toFixed(3)).join(', ')}`);
});

// ---- 6. a subspace lock in four dimensions
section('a subspace lock in four dimensions', () => {
  PhysiN.scripts.worker = 'physiN_worker.js';
  // There is no ground here. Only the joint holds the cube against gravity.
  const scene = new PhysiN.Scene({ dimensions: 4, gravity: [0, -9.81, 0, 0] });
  const cube = new PhysiN.HyperBoxMesh([0.5, 0.5, 0.5, 0.5], {}, 1);
  cube.setPositionN([0, 2, 0, 0]);
  scene.add(cube);
  scene.addConstraint(new PhysiN.SubspaceJoint(cube, null, {
    lockAxes: [1],
    worldFrame: true,
  }));
  for (let s = 0; s < 300; s += 1) scene.simulate(1 / 120, 1);
  check('the lock holds the body against gravity in four dimensions',
    Math.abs(cube.getPositionN()[1] - 2) < 0.01,
    `y = ${cube.getPositionN()[1].toFixed(4)}`);
  cube.setLinearVelocity([1, 0, 0, 0]);
  for (let s = 0; s < 120; s += 1) scene.simulate(1 / 120, 1);
  check('the axes that the lock leaves free still move',
    cube.getPositionN()[0] > 0.9, `x = ${cube.getPositionN()[0].toFixed(4)}`);
});

console.log(`\n${pass} tests pass, ${fail} tests fail`);
if (fail > 0) process.exitCode = 1;