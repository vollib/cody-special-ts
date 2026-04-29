import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import { erfCody, erfcCody, erfcxCody, inverseNormCdf, normCdf, normPdf, erfInv } from "../src/index.js";

function close(actual: number, expected: number, tolerance = 1e-14): void {
  assert.ok(Math.abs(actual - expected) <= tolerance || Math.abs(actual - expected) <= tolerance * Math.abs(expected), `${actual} != ${expected}`);
}

type NormalDistributionVectors = {
  tolerancePolicy: {
    absolute: number;
    relative: number;
  };
  cases: Array<{
    input: {
      z: number;
    };
    output: {
      normCdf: number;
    };
  }>;
};

const normalDistributionVectors = JSON.parse(
  readFileSync("node_modules/vollib-test-vectors/vectors/lets-be-rational/binary64/normal-distribution.json", "utf8")
) as NormalDistributionVectors;

test("Cody erf functions match reference values", () => {
  assert.equal(erfCody(0), 0);
  close(erfCody(1), 0.8427007929497149, 1e-15);
  close(erfCody(-1), -erfCody(1), 1e-15);
  assert.equal(erfcCody(0), 1);
  assert.ok(erfcCody(10) < 1e-40);
  close(erfcxCody(1), 0.427583576155807, 1e-14);
});

test("normal PDF/CDF and inverse CDF round trip", () => {
  close(normPdf(0), 0.3989422804014327, 1e-15);
  close(normPdf(1), normPdf(-1), 1e-15);
  close(normCdf(0), 0.5, 1e-15);
  close(normCdf(10), 1, 1e-15);
  assert.ok(normCdf(-10) < 1e-20);
  close(inverseNormCdf(0.5), 0, 1e-15);
  for (const x of [-2, -1, 0, 1, 2]) close(inverseNormCdf(normCdf(x)), x, 2e-14);
  for (const p of [0.1, 0.25, 0.5, 0.75, 0.9]) close(normCdf(inverseNormCdf(p)), p, 2e-14);
});

test("normal CDF matches shared C++ binary64 vectors", () => {
  const tolerance = Math.max(normalDistributionVectors.tolerancePolicy.absolute, normalDistributionVectors.tolerancePolicy.relative);
  for (const row of normalDistributionVectors.cases) {
    close(normCdf(row.input.z), row.output.normCdf, tolerance);
  }
});

test("PJ-2024 inverse normal tail values are stable", () => {
  close(inverseNormCdf(1e-20), -9.262340089798409, 1e-14);
  close(inverseNormCdf(1 - 1e-12), 7.0344869100478356, 1e-12);
  close(erfInv(0), 0, 1e-15);
  close(erfInv(0.5), 0.4769362762044699, 1e-14);
});
