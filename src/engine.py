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

_DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}


def normalize_dna(raw_input):
    return re.sub(r"\s+", "", (raw_input or "")).upper()


def validate_dna(raw_input):
    trimmed = (raw_input or "").strip()
    if not trimmed:
        return {"valid": False, "error": "DNA sequence is required."}
    if not _DNA_VALID_PATTERN.match(trimmed):
        return {"valid": False, "error": "DNA must only contain A, T, C, and G."}
    return {"valid": True, "error": None}


def clean_sequence(raw):
    """
    Clean a raw sequence string:
    - Skip FASTA header lines (starting with '>')
    - Remove whitespace, digits, and all non-letter characters
    - Convert to UPPERCASE
    """
    lines = (raw or "").strip().splitlines()
    cleaned = []
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            continue
        line = re.sub(r"[^A-Za-z]", "", line)
        cleaned.append(line.upper())
    return "".join(cleaned)


def validate_sequence(seq):
    """
    Validate sequence length, empty status, and ensure it doesn't
    contain both T and U or invalid characters.

    Returns: {"valid": True/False, "error": None/message}
    """
    if not seq:
        return {"valid": False, "error": "Sequence is empty after cleaning."}
    if len(seq) < 3:
        return {"valid": False, "error": "Sequence too short (minimum 3 bases)."}

    has_t = "T" in seq
    has_u = "U" in seq
    if has_t and has_u:
        return {"valid": False, "error": "Contains both T and U — not valid DNA or RNA."}

    valid = set("ATGCN") if (has_t or not has_u) else set("AUGCN")
    invalid = set(seq) - valid
    if invalid:
        return {
            "valid": False,
            "error": f'Invalid characters: {", ".join(sorted(invalid))}',
        }
    return {"valid": True, "error": None}


def detect_sequence_type(seq):
    """
    Rule: T present (no U) -> DNA; U present (no T) -> RNA.
    Returns: 'DNA' or 'RNA'
    """
    has_t = "T" in seq
    has_u = "U" in seq
    if has_t and has_u:
        raise ValueError("Ambiguous: sequence contains both T and U.")
    elif has_u:
        return "RNA"
    else:
        return "DNA"


def reverse_transcribe(mrna):
    """
    Reverse transcription: convert mRNA back to the DNA coding strand.
    Rule: replace every U with T.
    """
    return mrna.replace("U", "T")


def get_reverse_complement(seq):
    """
    Generate the reverse complement of a DNA sequence.
    Step 1: Complement each base. Step 2: Reverse the sequence.
    """
    complement = "".join(_DNA_COMPLEMENT.get(b, b) for b in seq)
    return complement[::-1]


def calculate_base_percentages(seq):
    """
    Returns a dictionary mapping each base letter to its percentage (rounded to 2 dp).
    """
    if not seq:
        raise ValueError("Cannot calculate base percentages of empty sequence.")
    total = len(seq)
    return {base: round(seq.count(base) / total * 100, 2) for base in sorted(set(seq))}


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
    dna = clean_sequence(dna)
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
        "reverseTranscription": reverse_transcribe(mrna),
        "reverseComplement": get_reverse_complement(normalized_dna),
        "basePercentages": calculate_base_percentages(normalized_dna),
    }
