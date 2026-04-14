import re

from src.genetic_code import (
    AMINO_TO_ONE_LETTER,
    CODON_TO_AMINO,
    DNA_CODING_TO_MRNA,
    DNA_TEMPLATE_TO_MRNA,
    PROTEIN_MOTIFS,
    RESIDUE_MASS_DA,
    STOP_CODONS,
)

_DNA_VALID_PATTERN = re.compile(r"^[ATCG\s]+$", re.IGNORECASE)


def normalize_dna(raw_input):
    return re.sub(r"\s+", "", (raw_input or "")).upper()


def validate_dna(raw_input):
    trimmed = (raw_input or "").strip()
    if not trimmed:
        return {"valid": False, "error": "DNA sequence is required."}
    if not _DNA_VALID_PATTERN.match(trimmed):
        return {"valid": False, "error": "DNA must only contain A, T, C, and G."}
    return {"valid": True, "error": None}


def transcribe_dna(dna, strand_type="template"):
    normalized = normalize_dna(dna)
    mapping = DNA_CODING_TO_MRNA if strand_type == "coding" else DNA_TEMPLATE_TO_MRNA
    return "".join(mapping.get(base, "N") for base in normalized)


def split_into_codons(mrna):
    return [mrna[i : i + 3] for i in range(0, len(mrna) - 2, 3)]


def translate_mrna(mrna):
    codons = split_into_codons(mrna)
    codon_details = []
    amino_acids = []

    for codon in codons:
        amino = CODON_TO_AMINO.get(codon, "?")
        if codon == "AUG":
            codon_type = "start"
        elif codon in STOP_CODONS:
            codon_type = "stop"
        elif amino == "?":
            codon_type = "invalid"
        else:
            codon_type = "normal"

        codon_details.append({"codon": codon, "amino": amino, "type": codon_type})

        if codon_type == "stop":
            break
        if amino != "?":
            amino_acids.append(amino)

    return {
        "codonDetails": codon_details,
        "aminoAcids": amino_acids,
        "proteinChain": "-".join(amino_acids),
        "encounteredStop": any(d["type"] == "stop" for d in codon_details),
    }


def gc_content_percent(dna):
    normalized = normalize_dna(dna)
    if not normalized:
        return 0
    gc = sum(1 for ch in normalized if ch in "GC")
    return (gc / len(normalized)) * 100


def estimate_molecular_weight(amino_acids):
    if not amino_acids:
        return 0
    residue_total = sum(RESIDUE_MASS_DA.get(aa, 0) for aa in amino_acids)
    water_loss = (len(amino_acids) - 1) * 18.015
    return max(0, residue_total - water_loss)


def to_one_letter_protein(amino_acids):
    return "".join(AMINO_TO_ONE_LETTER.get(aa, "X") for aa in amino_acids)


def find_protein_motifs(one_letter_protein):
    hits = []
    for motif in PROTEIN_MOTIFS:
        if re.search(motif["regex"], one_letter_protein):
            hits.append({"name": motif["name"], "description": motif["description"]})
    return hits


def process_sequence(dna, strand_type="template"):
    validation = validate_dna(dna)
    if not validation["valid"]:
        return {"error": validation["error"]}

    normalized_dna = normalize_dna(dna)
    mrna = transcribe_dna(normalized_dna, strand_type)
    translation = translate_mrna(mrna)
    one_letter = to_one_letter_protein(translation["aminoAcids"])
    motifs = find_protein_motifs(one_letter)

    return {
        "error": None,
        "dna": normalized_dna,
        "mrna": mrna,
        **translation,
        "oneLetterProtein": one_letter,
        "gcContent": gc_content_percent(normalized_dna),
        "molecularWeight": estimate_molecular_weight(translation["aminoAcids"]),
        "motifs": motifs,
    }
