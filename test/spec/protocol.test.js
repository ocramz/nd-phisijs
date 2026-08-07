/**
 * test/spec/protocol.test.js -- src/nd/workerCore.js and the plugin
 *
 * The message protocol. `workerCore.js` writes a binary report, and
 * `PhysiN.Scene` reads it. The two calculate the stride of that report
 * separately, and CLAUDE.md says that the two must always agree.
 *
 * The two bundles hold the same engine, word for word. Thus the same commands
 * must give the same answer in the main thread and in the worker. That is the
 * test that finds a change that went into one file only.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { loadPhysiNFor, FakeWorker } from '../lib/load.js';
import { anyN, anyPositive, anyVector, anyRadius } from '../lib/arbitraries.js';
import { assertArrayClose } from '../lib/numeric.js';

describe('the world report', () => {
  it('should always give each body its own state back', () => {
    // The stride of the report is 1 + n + r + n + k. `workerCore.js` and
    // `PhysiN.Scene` work it out separately. A stride that does not agree
    // mixes the bodies: the state of one body goes to another. Thus put each
    // body at a place of its own, run a step with no gravity, and see that
    // each body has its own place back.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), fc.integer({ min: 1, max: 8 }),
    )), ([n, count]) => {
      // Arrange
      const PhysiN = loadPhysiNFor(false);
      const scene = new PhysiN.Scene({ dimensions: n, gravity: new Array(n).fill(0) });
      const places = [];
      const boxes = [];
      for (let i = 0; i < count; i += 1) {
        // 10 of distance on the axis 0, thus no two bodies ever touch.
        const place = new Array(n).fill(0).map((_, axis) => (axis === 0 ? i * 10 : axis + i));
        const box = new PhysiN.HyperBoxMesh(new Array(n).fill(0.5), {}, 1);
        box.setPositionN(place);
        scene.add(box);
        places.push(place);
        boxes.push(box);
      }

      // Act
      scene.simulate(1 / 60, 1);

      // Assert
      for (let i = 0; i < count; i += 1) {
        assertArrayClose(boxes[i].getPositionN(), places[i],
          `the body ${i} of ${count}, in ${n} dimensions`, 1e-6);
      }
    }));
  });

  it('should always give the state of a body that moves', () => {
    // The same test, with a body that moves: the report must carry the place
    // and the velocity of each step.
    fc.assert(fc.property(anyN(2, 5).chain((n) => fc.tuple(
      fc.constant(n), anyVector(n, fc.double({ min: -3, max: 3, noNaN: true })),
      fc.integer({ min: 1, max: 30 }),
    )), ([n, velocity, steps]) => {
      // Arrange
      const PhysiN = loadPhysiNFor(false);
      const scene = new PhysiN.Scene({ dimensions: n, gravity: new Array(n).fill(0) });
      const ball = new PhysiN.HyperBoxMesh(new Array(n).fill(0.5), {}, 1);
      scene.add(ball);
      ball.setLinearVelocity(Array.from(velocity));

      // Act
      for (let s = 0; s < steps; s += 1) scene.simulate(1 / 60, 1);

      // Assert -- the report is a Float32Array, thus the tolerance is the
      // accuracy of a float of 32 bits.
      const expected = Array.from(velocity, (v) => v * steps / 60);
      assertArrayClose(ball.getPositionN(), expected, 'the place of a body that moves', 1e-5);
      assertArrayClose(ball.getLinearVelocity(), velocity, 'the velocity of the body', 1e-6);
    }));
  });
});

describe('the two bundles', () => {
  it('should always give the same answer in the main thread and in the worker', () => {
    // physin.js and physin_worker.js hold the same engine, word for word.
    // The plugin calls it directly when `PhysiN.scripts.worker` is null, and
    // it sends messages when that names a file. The two conditions must give
    // the same numbers, and not almost the same numbers: it is the same code
    // and the same order of operations.
    fc.assert(fc.property(anyScript(), (script) => {
      // Act
      const one = runScript(loadPhysiNFor(false), script);
      const other = runScript(loadPhysiNFor(true), script);

      // Assert
      assert.equal(other.length, one.length, 'the count of the bodies');
      for (let i = 0; i < one.length; i += 1) {
        assert.deepEqual(other[i].position, one[i].position, `the place of the body ${i}`);
        assert.deepEqual(other[i].rotor, one[i].rotor, `the orientation of the body ${i}`);
        assert.deepEqual(other[i].velocity, one[i].velocity, `the velocity of the body ${i}`);
        assert.deepEqual(other[i].angularVelocity, one[i].angularVelocity,
          `the angular velocity of the body ${i}`);
      }
    }), { numRuns: 20 });
  });

  it('should start a classic worker at the file of PhysiN.scripts.worker', () => {
    // Arrange
    const PhysiN = loadPhysiNFor(true);

    // Act
    const scene = new PhysiN.Scene({ dimensions: 4 });

    // Assert
    assert.equal(scene._worker.url, 'physin_worker.js', 'the file of the worker');
    assert.equal(scene._worker.type, 'classic', 'the type of the worker');
  });
});

/**
 * A KNOWN DEFECT, that this test found.
 *
 * README section 8 gives `PhysiN.HyperSphereMesh` under the head "For 4
 * dimensions and more". The class holds `n: 4` in its shape, and it takes no
 * count of dimensions:
 *
 *     { type: "sphere", n: 4, radius }
 *
 * Thus the plugin builds the state of the body with 4 dimensions, whatever
 * the scene has. At 2, 3 and 5 dimensions the first `setPositionN` or
 * `scene.add` stops with `RangeError: offset is out of bounds`, from inside
 * `Body` or `setPositionN`. The message does not name the true cause.
 *
 * `HyperBoxMesh` takes the count from the length of its half extents, and
 * `HyperPlaneMesh` from the length of its normal. `HyperTorusMesh` takes
 * `options.dimensions`. Only the ball has no way.
 *
 * Take the `todo` mark away when the defect is repaired.
 */
describe('PhysiN.HyperSphereMesh', () => {
  it('should work in each count of dimensions',
    { todo: 'the class holds n = 4 -- see the note above' }, () => {
    fc.assert(fc.property(anyN(2, 5), (n) => {
      // Arrange
      const PhysiN = loadPhysiNFor(false);
      const scene = new PhysiN.Scene({ dimensions: n, gravity: new Array(n).fill(0) });

      // Act
      const ball = new PhysiN.HyperSphereMesh(0.5, {}, 1);
      ball.setPositionN(new Array(n).fill(1));
      scene.add(ball);
      scene.simulate(1 / 60, 1);

      // Assert
      assertArrayClose(ball.getPositionN(), new Array(n).fill(1),
        `the place of a ball of ${n} dimensions`, 1e-6);
    }));
  });
});

describe('the engine of the worker', () => {
  it('should answer a command that it does not know, and go on', () => {
    // Arrange
    const answers = [];
    const worker = new FakeWorker('physin_worker.js', { type: 'classic' });
    worker.onmessage = (event) => answers.push(event.data);

    // Act
    worker.postMessage({ cmd: 'fly', params: {} });
    worker.postMessage({ cmd: 'init', params: { dimensions: 3, gravity: [0, -9.81, 0] } });

    // Assert
    assert.deepEqual(answers[0], { cmd: 'unknown', params: 'fly' }, 'it names the command');
    assert.deepEqual(answers[1], { cmd: 'worldReady' }, 'the engine still works');
  });
});

// Helpers

/**
 * A short list of commands for a scene: the count of dimensions, the bodies
 * and the impulses. The same list runs in the two conditions.
 */
function anyScript() {
  return anyN(3, 4).chain((n) => fc.record({
    n: fc.constant(n),
    bodies: fc.array(fc.record({
      isBox: fc.boolean(),
      radius: anyRadius(0.3, 0.8),
      height: anyPositive(1, 6),
      impulse: anyVector(n, fc.double({ min: -4, max: 4, noNaN: true })),
      // A torque is a bivector of k = n (n - 1) / 2 components.
      torque: anyVector((n * (n - 1)) / 2, fc.double({ min: -2, max: 2, noNaN: true })),
    }), { minLength: 1, maxLength: 4 }),
    steps: fc.integer({ min: 1, max: 40 }),
  }));
}

/**
 * Builds the scene of a script, runs it, and gives the state of each body.
 * The state comes from the reports, thus it holds what the plugin shows.
 */
function runScript(PhysiN, script) {
  const { n } = script;
  const gravity = new Array(n).fill(0);
  gravity[1] = -9.81;
  const scene = new PhysiN.Scene({ dimensions: n, gravity });
  const ground = new PhysiN.HyperBoxMesh(new Array(n).fill(0.5).map((v, i) => (i === 1 ? 0.5 : 8)),
    {}, 0);
  ground.setPositionN(new Array(n).fill(0).map((v, i) => (i === 1 ? -0.5 : 0)));
  scene.add(ground);

  const bodies = [];
  script.bodies.forEach((spec, i) => {
    // `HyperSphereMesh` only works at n = 4. See the test of that below.
    const body = spec.isBox || n !== 4
      ? new PhysiN.HyperBoxMesh(new Array(n).fill(0.4), {}, 1)
      : new PhysiN.HyperSphereMesh(spec.radius, {}, 1);
    body.setPositionN(new Array(n).fill(0).map((v, axis) => {
      if (axis === 0) return i * 3;
      if (axis === 1) return spec.height;
      return 0;
    }));
    scene.add(body);
    body.applyCentralImpulse(Array.from(spec.impulse));
    body.applyTorque(Array.from(spec.torque));
    bodies.push(body);
  });

  for (let s = 0; s < script.steps; s += 1) scene.simulate(1 / 60, 1);

  return bodies.map((body) => ({
    position: Array.from(body.getPositionN()),
    rotor: Array.from(body.getRotor()),
    velocity: Array.from(body.getLinearVelocity()),
    angularVelocity: Array.from(body.getAngularVelocity()),
  }));
}
