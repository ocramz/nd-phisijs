/**
 * test/lib/load.js -- loads the two bundles into Node.
 *
 * The bundles are classic scripts, and not modules. `physin.js` puts `PhysiN`
 * on the global object, but only when `globalThis.THREE` is there.
 * `physin_worker.js` waits for `self.onmessage`. Thus a test cannot import
 * them. It builds a function from the text of the file, and it gives that
 * function its own `globalThis` or its own `self`.
 *
 * `test/worker.js` does the same thing. This file holds the method one time.
 */
import fs from 'fs';
import { THREE } from './three-stub.js';

const mainCode = fs.readFileSync(new URL('../../physin.js', import.meta.url), 'utf8');
const workerCode = fs.readFileSync(new URL('../../physin_worker.js', import.meta.url), 'utf8');

/**
 * A new `PhysiN`, with its own state. Each call builds the bundle again, thus
 * one test cannot change `PhysiN.scripts.worker` for another test.
 *
 * @param {Function} [Worker] the class of the worker. The plugin only calls
 *   it when `PhysiN.scripts.worker` names a file.
 * @returns {object} `PhysiN`. `PhysiN.nd` is the engine, and `PhysiN.slice`
 *   is the cut of 4 dimensions.
 */
export function loadPhysiN(Worker) {
  // eslint-disable-next-line no-new-func
  const build = new Function('THREE', 'Worker', `
    var globalThis = { THREE: THREE };
    ${mainCode}
    return globalThis.PhysiN;
  `);
  return build(THREE, Worker);
}

/** One `PhysiN` for all of the specs of one file. */
const shared = loadPhysiN();

/**
 * The engine, `PhysiN.nd`. Almost every spec starts with this.
 *
 * CAUTION: use this one object, and do not build a second engine. The tables
 * of `dims(n)` are in a cache, and some functions keep their scratch memory
 * in a `WeakMap` on those tables. A body from one engine and a function from
 * another engine do not work together.
 */
export const nd = shared.nd;

/** The cut of 4 dimensions, `PhysiN.slice`. */
export const slice = shared.slice;

/**
 * One instance of the worker script, with its own `self`.
 * @param {function(*): void} onPost takes each message that the worker sends
 */
function makeWorkerScript(onPost) {
  const self = { onmessage: null, postMessage(data) { onPost(data); } };
  // eslint-disable-next-line no-new-func
  new Function('self', workerCode)(self);
  return self;
}

/**
 * A `Worker` that is a stub. It runs `physin_worker.js` in this thread, and
 * it sends each message at once. A true worker sends the message later, but
 * the order of the messages is the same.
 *
 * `counts` holds the count of the messages of the two directions, for a test
 * that must see that the two files speak.
 */
export class FakeWorker {
  constructor(url, opts) {
    this.url = url;
    this.type = opts && opts.type;
    this.onmessage = null;
    this.counts = { toWorker: 0, fromWorker: 0 };
    this._self = makeWorkerScript((data) => {
      this.counts.fromWorker += 1;
      if (this.onmessage) this.onmessage({ data });
    });
  }
  postMessage(data) {
    this.counts.toWorker += 1;
    if (this._self.onmessage) this._self.onmessage({ data });
  }
}

/**
 * A `PhysiN` that does the physics in the given condition.
 * @param {boolean} inWorker true puts the engine in the worker file, and
 *   false keeps it in the main thread
 */
export function loadPhysiNFor(inWorker) {
  const PhysiN = loadPhysiN(FakeWorker);
  PhysiN.scripts.worker = inWorker ? 'physin_worker.js' : null;
  return PhysiN;
}
