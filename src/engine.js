import {
    AMINO_TO_ONE_LETTER,
    CODON_TO_AMINO,
    DNA_CODING_TO_MRNA,
    DNA_TEMPLATE_TO_MRNA,
    PROTEIN_MOTIFS,
    RESIDUE_MASS_DA,
    STOP_CODONS,
} from "./geneticCode.js";

const DNA_VALID_PATTERN = /^[ATCG\s]+$/i;

export function normalizeDNA(input) {
    return (input || "").toUpperCase().replace(/\s+/g, "");
}

export function validateDNA(input) {
    const trimmed = (input || "").trim();
    if (!trimmed) {
        return { valid: false, error: "DNA sequence is required." };
    }
    if (!DNA_VALID_PATTERN.test(trimmed)) {
        return { valid: false, error: "DNA must only contain A, T, C, and G." };
    }
    return { valid: true, error: null };
}

export function transcribeDNA(dna, strandType = "template") {
    const normalized = normalizeDNA(dna);
    const map =
        strandType === "coding" ? DNA_CODING_TO_MRNA : DNA_TEMPLATE_TO_MRNA;
    let mrna = "";

    for (const base of normalized) {
        mrna += map[base] ?? "N";
    }

    return mrna;
}

export function splitIntoCodons(mrna) {
    const codons = [];
    for (let i = 0; i + 2 < mrna.length; i += 3) {
        codons.push(mrna.slice(i, i + 3));
    }
    return codons;
}

export function translateMRNA(mrna) {
    const codons = splitIntoCodons(mrna);
    const codonDetails = [];
    const aminoAcids = [];

    for (const codon of codons) {
        const amino = CODON_TO_AMINO[codon] ?? "?";
        const type =
            codon === "AUG"
                ? "start"
                : STOP_CODONS.has(codon)
                  ? "stop"
                  : amino === "?"
                    ? "invalid"
                    : "normal";

        codonDetails.push({ codon, amino, type });

        if (type === "stop") {
            break;
        }
        if (amino !== "?") {
            aminoAcids.push(amino);
        }
    }

    return {
        codonDetails,
        aminoAcids,
        proteinChain: aminoAcids.join("-"),
        encounteredStop: codonDetails.some((x) => x.type === "stop"),
    };
}

export function gcContentPercent(dna) {
    const normalized = normalizeDNA(dna);
    if (!normalized.length) return 0;
    const gc = (normalized.match(/[GC]/g) || []).length;
    return (gc / normalized.length) * 100;
}

export function estimateMolecularWeight(aminoAcids) {
    if (!aminoAcids.length) return 0;

    const residueTotal = aminoAcids.reduce(
        (sum, aa) => sum + (RESIDUE_MASS_DA[aa] || 0),
        0,
    );
    const waterLoss = (aminoAcids.length - 1) * 18.015;
    return Math.max(0, residueTotal - waterLoss);
}

export function toOneLetterProtein(aminoAcids) {
    return aminoAcids.map((aa) => AMINO_TO_ONE_LETTER[aa] || "X").join("");
}

export function findProteinMotifs(oneLetterProtein) {
    const hits = [];
    for (const motif of PROTEIN_MOTIFS) {
        if (motif.regex.test(oneLetterProtein)) {
            hits.push(motif);
        }
    }
    return hits;
}

export function processSequence({ dna, strandType }) {
    const validation = validateDNA(dna);
    if (!validation.valid) {
        return { error: validation.error };
    }

    const normalizedDNA = normalizeDNA(dna);
    const mrna = transcribeDNA(normalizedDNA, strandType);
    const translation = translateMRNA(mrna);
    const oneLetter = toOneLetterProtein(translation.aminoAcids);
    const motifs = findProteinMotifs(oneLetter);

    return {
        error: null,
        dna: normalizedDNA,
        mrna,
        ...translation,
        oneLetterProtein: oneLetter,
        gcContent: gcContentPercent(normalizedDNA),
        molecularWeight: estimateMolecularWeight(translation.aminoAcids),
        motifs,
    };
}
