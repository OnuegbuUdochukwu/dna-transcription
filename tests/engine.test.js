import test from "node:test";
import assert from "node:assert/strict";
import { processSequence } from "../src/engine.js";

test("Template DNA TACGGC transcribes to AUGCCG", () => {
    const result = processSequence({ dna: "TACGGC", strandType: "template" });
    assert.equal(result.error, null);
    assert.equal(result.mrna, "AUGCCG");
});

test("Translation stops at internal stop codon", () => {
    const result = processSequence({
        dna: "ATGAAATAATTT",
        strandType: "coding",
    });
    assert.equal(result.error, null);
    assert.equal(result.mrna, "AUGAAAUAAUUU");
    assert.equal(result.proteinChain, "Met-Lys");
    assert.equal(result.encounteredStop, true);
});

test("Batch-like large sequence input remains valid", () => {
    const dna = "ATGC".repeat(200);
    const result = processSequence({ dna, strandType: "coding" });
    assert.equal(result.error, null);
    assert.ok(result.mrna.length > 0);
});
