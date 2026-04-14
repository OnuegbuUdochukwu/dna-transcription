from flask import Flask, jsonify, request, send_from_directory
import os

from src.engine import process_sequence

app = Flask(__name__, static_folder=os.path.join(os.path.dirname(__file__), ".."))


@app.route("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.route("/<path:path>")
def static_files(path):
    return send_from_directory(app.static_folder, path)


@app.route("/api/process", methods=["POST"])
def api_process():
    data = request.get_json(force=True)
    dna = data.get("dna", "")
    strand_type = data.get("strandType", "template")
    result = process_sequence(dna, strand_type)
    return jsonify(result)


@app.route("/api/process-batch", methods=["POST"])
def api_process_batch():
    data = request.get_json(force=True)
    strand_type = data.get("strandType", "template")
    sequences = data.get("sequences", [])

    results = []
    for seq in sequences:
        result = process_sequence(seq, strand_type)
        results.append(result)

    return jsonify(results)


if __name__ == "__main__":
    app.run(host="0.0.0.0", port=8080)
