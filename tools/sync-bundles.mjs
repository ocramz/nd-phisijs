/**
 * tools/sync-bundles.mjs -- copies the shared engine modules of `physin.js`
 * into `physin_worker.js`.
 *
 * The two bundles hold the same engine modules, word for word, but NOT in the
 * same order. Thus a copy must find each module by its `// src/...` marker and
 * put it in the place of the module of the same name. A copy that works on a
 * range of lines deletes the module that comes next, because the module that
 * comes next is not the same in the two files.
 *
 * The four algebra modules are not in the list. The worker bundle leaves out
 * the functions that the worker never calls, thus those modules are shorter
 * there. `test/worker.js` holds the same list, and it fails when a module of
 * the list is not the same in the two files.
 *
 *   node tools/sync-bundles.mjs          copies, and says what changed
 *   node tools/sync-bundles.mjs --check  says what differs, and changes nothing
 */
import fs from 'fs';

const SHARED = [
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

/** The first and the last line of each module of a bundle. */
function moduleSpans(code) {
  const lines = code.split('\n');
  const marks = [];
  for (let i = 0; i < lines.length; i += 1) {
    if (/^\s*\/\/ src\/\S+\.js$/.test(lines[i])) marks.push([i, lines[i].trim().slice(3)]);
  }
  const out = new Map();
  for (let m = 0; m < marks.length; m += 1) {
    const end = m + 1 < marks.length ? marks[m + 1][0] : lines.length;
    out.set(marks[m][1], { start: marks[m][0], end, text: lines.slice(marks[m][0], end).join('\n') });
  }
  return { lines, out };
}

const check = process.argv.includes('--check');
const mainPath = new URL('../physin.js', import.meta.url);
const workerPath = new URL('../physin_worker.js', import.meta.url);
const main = moduleSpans(fs.readFileSync(mainPath, 'utf8'));
const worker = moduleSpans(fs.readFileSync(workerPath, 'utf8'));

const missing = SHARED.filter((m) => !main.out.has(m) || !worker.out.has(m));
if (missing.length > 0) {
  console.error(`a module of the list is not in a bundle: ${missing.join(', ')}`);
  process.exit(1);
}

const differ = SHARED.filter((m) => main.out.get(m).text !== worker.out.get(m).text);
if (differ.length === 0) {
  console.log('the two bundles already hold the same engine.');
  process.exit(0);
}
if (check) {
  console.error(`these modules differ: ${differ.join(', ')}`);
  process.exit(1);
}

// Work from the last module to the first, thus an earlier span stays correct.
const spans = differ
  .map((m) => ({ name: m, span: worker.out.get(m), text: main.out.get(m).text }))
  .sort((p, q) => q.span.start - p.span.start);
let lines = worker.lines;
for (const s of spans) {
  lines = lines.slice(0, s.span.start).concat(s.text.split('\n'), lines.slice(s.span.end));
  console.log(`copied ${s.name}`);
}
fs.writeFileSync(workerPath, lines.join('\n'));
console.log(`${differ.length} module(s) copied into physin_worker.js`);
