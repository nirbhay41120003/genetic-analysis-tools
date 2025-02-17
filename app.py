from flask import Flask, request, render_template
import os
import subprocess
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from sequence_alignment.aligner import align_sequences

app = Flask(__name__)
UPLOAD_FOLDER = "uploads"
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/align", methods=["POST"])
def align():
    file = request.files["file"]
    filepath = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(filepath)

    sequences = list(SeqIO.parse(filepath, "fasta"))
    if len(sequences) < 2:
        return "Upload a FASTA file with at least two sequences."

    alignment = align_sequences(sequences[0].seq, sequences[1].seq, method="global")
    return f"<pre>{alignment}</pre>"

@app.route("/variants", methods=["POST"])
def call_variants():
    file = request.files["bam_file"]
    filepath = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(filepath)

    result = subprocess.run(["bash", "variant_calling/pipeline.sh"], capture_output=True, text=True)
    return f"<pre>{result.stdout}</pre>"

@app.route("/gene_expression", methods=["POST"])
def gene_expression():
    df = pd.read_csv("gene_expression/gene_counts_tpm.csv")
    fig = px.imshow(df.iloc[:, 1:], labels=dict(x="Sample", y="Gene", color="Expression Level"))
    return fig.to_html()

if __name__ == "__main__":
    app.run(debug=True)
