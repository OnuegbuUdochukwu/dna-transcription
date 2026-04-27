import os
import sys
import re

from src.engine import (
    clean_sequence,
    detect_sequence_type,
    validate_sequence,
    transcribe_dna,
    reverse_transcribe,
    get_reverse_complement,
    translate_mrna,
    gc_content_percent,
    calculate_base_percentages,
    to_one_letter_protein,
    estimate_molecular_weight,
    find_protein_motifs,
)


def get_sequence_input():
    """
    Accept a nucleotide sequence from the user via:
    1. Direct typing or pasting
    2. File upload (.txt or .fasta)
    Returns: raw sequence string
    Raises: ValueError, FileNotFoundError
    """
    print("How would you like to provide your sequence?")
    print("  1. Type or paste directly")
    print("  2. Upload a file (.txt or .fasta)")
    choice = input("Enter choice (1 or 2): ").strip()

    if choice == "1":
        sequence = input("Enter your sequence: ").strip()
        if not sequence:
            raise ValueError("No sequence entered.")
        return sequence
    elif choice == "2":
        path = input("Enter file path: ").strip()
        if not os.path.exists(path):
            raise FileNotFoundError(f"File not found: {path}")
        with open(path, "r") as f:
            content = f.read()
        if not content.strip():
            raise ValueError("File is empty.")
        return content
    else:
        raise ValueError("Invalid choice. Enter 1 or 2.")


def get_strand_type():
    """
    Ask user to specify if the DNA provided is Template or Non-Template.
    """
    print("Which DNA strand have you provided?")
    print("  1. Non-Template Strand (Coding/Sense strand)")
    print("  2. Template Strand (Antisense strand)")
    choice = input("Enter choice (1 or 2): ").strip()
    if choice == "1":
        return "coding"
    elif choice == "2":
        return "template"
    else:
        raise ValueError("Invalid choice. Enter 1 or 2.")


def display_results(seq, seq_type, strand_type, mrna, translation, gc, base_pcts, mw, motifs):
    """
    Print a formatted summary including detection, transcription,
    and translation details.
    """
    LINE = "=" * 60
    one_letter = to_one_letter_protein(translation["aminoAcids"])

    print(f"\n{LINE}")
    print("  DNA TO PROTEIN — RESULTS")
    print(LINE)

    print(f"\n[DETECTION]")
    print(f"  Type:   {seq_type}")
    print(f"  Length: {len(seq)} bases")
    if seq_type == "DNA":
        print(f"  Strand: {strand_type}")

    print(f"\n[TRANSCRIPTION]")
    print(f"  mRNA:               {mrna}")
    if seq_type == "DNA":
        print(f"  Reverse transcription (mRNA→DNA): {reverse_transcribe(mrna)}")
        print(f"  Reverse complement (DNA):          {get_reverse_complement(seq)}")

    print(f"\n[TRANSLATION]")
    print(f"  Protein chain:      {translation['proteinChain'] or '(none)'}")
    print(f"  Protein (1-letter): {one_letter or '(none)'}")
    print(f"  Amino acid count:   {len(translation['aminoAcids'])}")
    print(f"  Stop codon found:   {'Yes' if translation['encounteredStop'] else 'No'}")
    print(f"  Est. mol. weight:   {mw:.2f} Da")

    print(f"\n[STATISTICS]")
    print(f"  GC content: {gc:.2f}%")
    print(f"  Base percentages:")
    for base, pct in sorted(base_pcts.items()):
        print(f"    {base}: {pct:.2f}%")

    if motifs:
        print(f"\n[MOTIFS]")
        for m in motifs:
            print(f"  - {m['name']}: {m['description']}")

    print(f"\n{LINE}\n")


def main():
    """
    Orchestrates the full CLI pipeline with exhaustive error handling.
    """
    try:
        raw = get_sequence_input()
        seq = clean_sequence(raw)

        seq_validation = validate_sequence(seq)
        if not seq_validation["valid"]:
            raise ValueError(seq_validation["error"])

        seq_type = detect_sequence_type(seq)

        if seq_type == "DNA":
            strand_type = get_strand_type()
            mrna = transcribe_dna(seq, strand_type)
        else:
            strand_type = "rna_input"
            mrna = seq

        translation = translate_mrna(mrna)
        gc = gc_content_percent(seq)
        base_pcts = calculate_base_percentages(seq)
        one_letter = to_one_letter_protein(translation["aminoAcids"])
        mw = estimate_molecular_weight(translation["aminoAcids"])
        motifs = find_protein_motifs(one_letter)

        display_results(seq, seq_type, strand_type, mrna, translation, gc, base_pcts, mw, motifs)

    except (ValueError, FileNotFoundError) as e:
        print(f"ERROR: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
