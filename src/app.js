import { processSequence } from "./engine.js";

const dnaInput = document.getElementById("dnaInput");
const strandType = document.getElementById("strandType");
const processBtn = document.getElementById("processBtn");
const inputError = document.getElementById("inputError");

const resultsSection = document.getElementById("resultsSection");
const dnaOut = document.getElementById("dnaOut");
const mrnaOut = document.getElementById("mrnaOut");
const codonOut = document.getElementById("codonOut");
const proteinOut = document.getElementById("proteinOut");
const gcOut = document.getElementById("gcOut");
const aaCountOut = document.getElementById("aaCountOut");
const mwOut = document.getElementById("mwOut");
const blastLink = document.getElementById("blastLink");
const motifOut = document.getElementById("motifOut");

const csvInput = document.getElementById("csvInput");
const csvColumn = document.getElementById("csvColumn");
const processCsvBtn = document.getElementById("processCsvBtn");
const csvInfo = document.getElementById("csvInfo");
const csvResultsBody = document.getElementById("csvResultsBody");

let csvRows = [];

function setError(message) {
    if (!message) {
        inputError.classList.add("hidden");
        inputError.textContent = "";
        return;
    }

    inputError.textContent = message;
    inputError.classList.remove("hidden");
}

function codonBadge(detail) {
    const classByType = {
        start: "bg-green-100 text-green-900 border-green-300",
        stop: "bg-red-100 text-red-900 border-red-300",
        normal: "bg-slate-100 text-slate-900 border-slate-300",
        invalid: "bg-amber-100 text-amber-900 border-amber-300",
    };

    const span = document.createElement("span");
    span.className = `px-2 py-0.5 rounded border text-xs ${classByType[detail.type] || classByType.normal}`;
    span.textContent = `${detail.codon} (${detail.amino})`;
    return span;
}

function renderSingleResult(result) {
    resultsSection.classList.remove("hidden");

    dnaOut.textContent = result.dna;
    mrnaOut.textContent = result.mrna;

    codonOut.innerHTML = "";
    result.codonDetails.forEach((detail) => {
        codonOut.appendChild(codonBadge(detail));
    });

    proteinOut.textContent =
        result.proteinChain || "No translated amino acids detected.";
    gcOut.textContent = `${result.gcContent.toFixed(2)}%`;
    aaCountOut.textContent = `${result.aminoAcids.length}`;
    mwOut.textContent = `${result.molecularWeight.toFixed(2)} Da`;

    blastLink.href = result.oneLetterProtein
        ? `https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp&QUERY=${encodeURIComponent(result.oneLetterProtein)}`
        : "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins&PROGRAM=blastp";

    if (!result.oneLetterProtein) {
        motifOut.textContent = "No motif scan: protein chain is empty.";
    } else if (!result.motifs.length) {
        motifOut.textContent = "Motif scan: no quick-match motifs found.";
    } else {
        motifOut.textContent = `Motif scan: ${result.motifs.map((m) => m.name).join(", ")}.`;
    }
}

processBtn.addEventListener("click", () => {
    const result = processSequence({
        dna: dnaInput.value,
        strandType: strandType.value,
    });
    if (result.error) {
        setError(result.error);
        return;
    }

    setError(null);
    renderSingleResult(result);
});

csvInput.addEventListener("change", () => {
    const [file] = csvInput.files || [];
    if (!file) return;

    Papa.parse(file, {
        header: true,
        skipEmptyLines: true,
        complete: (results) => {
            csvRows = results.data || [];
            const headers = results.meta?.fields || [];

            csvColumn.innerHTML = "";
            headers.forEach((header) => {
                const option = document.createElement("option");
                option.value = header;
                option.textContent = header;
                csvColumn.appendChild(option);
            });

            csvInfo.textContent = `Loaded ${csvRows.length} rows.`;
        },
        error: (err) => {
            csvRows = [];
            csvInfo.textContent = `CSV parse error: ${err.message}`;
        },
    });
});

processCsvBtn.addEventListener("click", () => {
    if (!csvRows.length) {
        csvInfo.textContent = "Please upload a CSV file first.";
        return;
    }

    const column = csvColumn.value;
    if (!column) {
        csvInfo.textContent = "No DNA column selected.";
        return;
    }

    const t0 = performance.now();

    const rows = csvRows.map((row, i) => {
        const dna = row[column] || "";
        const result = processSequence({ dna, strandType: strandType.value });

        return {
            index: i + 1,
            dna: String(dna).trim(),
            mrna: result.error ? "" : result.mrna,
            protein: result.error ? "" : result.proteinChain,
            gc: result.error ? "" : `${result.gcContent.toFixed(2)}%`,
            status: result.error || "OK",
        };
    });

    csvResultsBody.innerHTML = "";
    rows.forEach((row) => {
        const tr = document.createElement("tr");
        tr.className = "border-t border-slate-200";
        tr.innerHTML = `
      <td class="p-2 align-top">${row.index}</td>
      <td class="p-2 align-top font-mono">${row.dna}</td>
      <td class="p-2 align-top font-mono">${row.mrna}</td>
      <td class="p-2 align-top">${row.protein}</td>
      <td class="p-2 align-top">${row.gc}</td>
      <td class="p-2 align-top ${row.status === "OK" ? "text-emerald-700" : "text-red-700"}">${row.status}</td>
    `;
        csvResultsBody.appendChild(tr);
    });

    const elapsed = performance.now() - t0;
    csvInfo.textContent = `Processed ${rows.length} rows in ${elapsed.toFixed(1)} ms.`;
});
