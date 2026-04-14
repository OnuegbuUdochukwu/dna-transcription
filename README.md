# DNA Transcription & Translation Suite

A lightweight browser-based bioinformatics tool with a **Python (Flask) backend** that:

- Transcribes DNA to mRNA from either **template** or **coding (non-template)** strands.
- Translates mRNA codons into amino-acid chains using the standard genetic code.
- Stops translation at `UAA`, `UAG`, or `UGA`.
- Highlights start (`AUG`) and stop codons in the UI.
- Computes quick metrics like **GC-content** and estimated **protein molecular weight**.
- Supports **CSV batch processing** with selectable DNA column.
- Adds a direct outbound link to NCBI BLAST for downstream protein identification.

## Architecture

The core transcription/translation logic is implemented in **Python**:

- `src/genetic_code.py` – codon tables, mappings, and constants.
- `src/engine.py` – DNA validation, transcription, translation, and statistics.
- `src/server.py` – Flask API server that exposes `/api/process` and `/api/process-batch` endpoints, and serves the static frontend.

The frontend (`index.html` + `src/app.js`) is a lightweight HTML/CSS/JS interface that calls the Python API via `fetch`.

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

### Prerequisites

- Python 3.10+
- pip

### Setup

```bash
pip install -r requirements.txt
```

### Start the Server

```bash
python3 src/server.py
```

Then visit `http://localhost:8080`.

### Run Tests

```bash
python3 -m pytest tests/ -v
```

## Suggested Next Enhancements

- Add configurable reading-frame selection.
- Integrate ORF detection and longest ORF extraction.
- Expand motif database with curated protein signatures.
- Add export options (CSV/FASTA) for processed output.
- Add automated tests for known benchmark sequences.
