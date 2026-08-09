/**
 * test/spec/stress.test.js -- how well a joint holds under a load.
 *
 * These tests do not check that the mathematics is correct;
 * `constraint.test.js` does that. They check HOW WELL the solver holds, and
 * they catch a change that makes a joint softer.
 *
 * THE NUMBERS ARE IN A BAND OF A FACTOR OF TWO. An exact number would break on
 * a rounding change and say nothing. The band catches a real regression. The
 * measurement of each band is written next to it, with its date.
 *
 * The tests that matter most are the ORDER tests at the end. They hold the
 * relations that must be true of any build, and they do not depend on a
 * measurement at all.
 */
import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import { nd } from '../lib/load.js';

const { dims, World, Body, HyperSphere, defaultParams } = nd;

const D = dims(3);

/**
 * A chain of `links` point joints, from a static anchor at (0, 10, 0). Each
 * link is 1 long, and the chain STARTS HORIZONTAL, thus it whips as it falls.
 * That is the hard case: the lever arms of a row stay still through the step,
 * and a body that moves fast leaves the row behind.
 *
 * @returns {number} the largest gap between the two anchors of any joint, over
 *   the whole run. This is the true error of the constraint.
 */
function chainError(links, endMass, params) {
  const world = new World({
    dimensions: 3,
    gravity: [0, -9.81, 0],
    params: Object.assign({}, defaultParams, { allowSleep: false }, params),
  });
  let prev = world.addBody(new Body(D, {
    shape: HyperSphere(D, 0.2), mass: 0, position: [0, 10, 0],
  }));
  const joints = [];
  for (let i = 1; i <= links; i += 1) {
    const b = world.addBody(new Body(D, {
      shape: HyperSphere(D, 0.2), mass: i === links ? endMass : 1, position: [i, 10, 0],
    }));
    joints.push(world.createConstraint(
      { type: 'point', localA: [0.5, 0, 0], localB: [-0.5, 0, 0] }, prev, b,
    ));
    prev = b;
  }
  let worst = 0;
  const pA = new Float64Array(3);
  const pB = new Float64Array(3);
  for (let s = 0; s < 1200; s += 1) {
    world.step(1 / 120);
    for (const j of joints) {
      j.a.localToWorld(j.localA, pA);
      j.b.localToWorld(j.localB, pB);
      worst = Math.max(worst, Math.hypot(pA[0] - pB[0], pA[1] - pB[1], pA[2] - pB[2]));
    }
  }
  return worst;
}

/** Stops the test when `value` is not inside a band of a factor of two. */
function assertBand(value, expected, message) {
  assert.ok(value > expected / 2 && value < expected * 2,
    `${message}: ${value.toFixed(5)} is not near ${expected} (a factor of two)`);
}

describe('a chain of point joints', () => {
  it('should hold to about 0.11 with the values of the defaults', () => {
    // Measured 2026-08-08 on the build of round 2: 0.11418.
    assertBand(chainError(10, 1, {}), 0.11, 'a chain of 10 links');
  });

  it('should not get better with more iterations, in the case that moves fast', () => {
    // THIS IS THE IMPORTANT ONE. `iterations` brings the VELOCITY of a row to
    // its target. It does nothing for the error of the POSITION of a chain that
    // whips, because that error comes from the lever arms that stay still
    // through the step. Only a shorter step helps there.
    //
    // A SETTLED chain is the other case, and it DOES get better with more
    // iterations: the same chain, started hanging straight down, goes from
    // 0.0014 at 10 iterations to 0.00003 at 100. Do not read this test as
    // "iterations never help".
    const ten = chainError(10, 1, {});
    const hundred = chainError(10, 1, { iterations: 100 });
    assert.ok(Math.abs(hundred - ten) / ten < 0.1,
      `10 iterations gave ${ten.toFixed(5)} and 100 gave ${hundred.toFixed(5)}, `
      + 'thus the iterations changed the answer. The error of the position is '
      + 'no longer limited by the bias, and this test must be looked at again.');
  });

  it('should get much better with sub-steps', () => {
    // Measured 2026-08-08: 0.01373 with four sub-steps, against 0.11418.
    const one = chainError(10, 1, {});
    const four = chainError(10, 1, { subSteps: 4 });
    assertBand(four, 0.0137, 'a chain of 10 links with four sub-steps');
    assert.ok(four < one / 4, `four sub-steps gave ${four.toFixed(5)}, and one gave ${one.toFixed(5)}`);
  });
});

describe('a large ratio of the masses', () => {
  it('should need sub-steps, and should hold with them', () => {
    // A body of 1000 at the end of a chain of bodies of 1. At one step of 1/120
    // the chain comes apart; measured 2026-08-08 as 180. With eight sub-steps
    // it holds to 0.44647.
    const one = chainError(5, 1000, {});
    const eight = chainError(5, 1000, { subSteps: 8 });
    assert.ok(one > 10,
      `a mass ratio of 1000 to 1 must come apart at one sub-step, and the error was ${one}`);
    assertBand(eight, 0.45, 'a mass ratio of 1000 to 1 with eight sub-steps');
    assert.ok(Number.isFinite(eight), 'the error must stay finite with sub-steps');
  });
});

describe('a soft joint', () => {
  it('should stretch more than a rigid joint', () => {
    // This holds for any build: a spring gives way, and a rigid joint does not.
    const rigid = chainError(10, 1, {});
    const soft = chainError(10, 1, { constraintHertz: 5 });
    assert.ok(soft > rigid,
      `a soft joint gave ${soft.toFixed(5)} and a rigid joint gave ${rigid.toFixed(5)}`);
  });
});
