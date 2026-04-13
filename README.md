# DNA Transcription & Translation Suite

A lightweight browser-based bioinformatics tool that:

- Transcribes DNA to mRNA from either **template** or **coding (non-template)** strands.
- Translates mRNA codons into amino-acid chains using the standard genetic code.
- Stops translation at `UAA`, `UAG`, or `UGA`.
- Highlights start (`AUG`) and stop codons in the UI.
- Computes quick metrics like **GC-content** and estimated **protein molecular weight**.
- Supports **CSV batch processing** with selectable DNA column.
- Adds a direct outbound link to NCBI BLAST for downstream protein identification.

## Biological Rules Implemented

### DNA → mRNA

- Coding strand: replace `T` with `U`.
- Template strand complement:
    - `A → U`
    - `T → A`
    - `C → G`
    - `G → C`

### Translation

- mRNA is split into codons (triplets) in reading frame from position 1.
- Codons map to amino acids using the standard codon table.
- Translation terminates at the first stop codon.

## Run Locally

This project is plain HTML + JavaScript and has no build step.

1. Open `index.html` in a browser, or
2. Serve with a local server (recommended for module loading).

Example using Python:

```bash
python3 -m http.server 8080
```

Then visit `http://localhost:8080`.

## Suggested Next Enhancements

- Add configurable reading-frame selection.
- Integrate ORF detection and longest ORF extraction.
- Expand motif database with curated protein signatures.
- Add export options (CSV/FASTA) for processed output.
- Add automated tests for known benchmark sequences.
