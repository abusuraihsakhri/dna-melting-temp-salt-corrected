/* Reproducible static browser-model smoke tests with Node built-ins only. */
"use strict";
const fs = require("node:fs");
const vm = require("node:vm");
const assert = require("node:assert/strict");
const html = fs.readFileSync("web/index.html", "utf8");
const script = html.match(/<script>([\s\S]*?)<\/script>/);
assert.ok(script, "Inline browser calculator script was not found");
const context = {
  document: { getElementById: () => ({ addEventListener: () => {} }) },
  URL: {}, Blob, setTimeout,
};
vm.runInNewContext(script[1] + "\nglobalThis.api = {calculate, parseInput};", context);
const { calculate, parseInput } = context.api;
const buffer = { na: 50, mg: 1.5, ct: 250 };
for (const [seq, expected] of [
  ["ATCGATCGATCGATCGATCG", 59.72],
  ["GCGCGCGCGCGCGCGCGCGC", 86.30],
  ["ATATATATATATATATATAT", 31.46],
  ["AGTCAGTCAGTCAGTCAGTC", 59.07],
]) {
  assert.ok(Math.abs(calculate(seq, buffer).combined - expected) < 0.03, seq);
}
assert.equal(Math.round(calculate("ATCGATCGATCGATCGATCG", { ...buffer, mg: 0 }).combined * 100) / 100, 53.09);
assert.throws(() => calculate("ATCGN", buffer), /A\/C\/G\/T/);
assert.equal(parseInput("id,sequence\nP1,ATCG\nP2,GCAT")[1].seq, "GCAT");
assert.equal(parseInput(">P1\nATCG\n>P2\nGCAT").length, 2);
console.log("Browser scientific calculations, input validation, CSV and FASTA parsing passed.");
