from flask import Flask, request, render_template
from aligner import align_sequences
from Bio import SeqIO
import io

app = Flask(__name__)

@app.route("/", methods=["GET", "POST"])
def upload_file():
    if request.method == "POST":
        file = request.files["file"]
        sequences = list(SeqIO.parse(io.StringIO(file.read().decode("utf-8")), "fasta"))

        if len(sequences) < 2:
            return "Upload a FASTA file with at least two sequences."

        alignment = align_sequences(sequences[0].seq, sequences[1].seq, method="global")
        return f"<pre>{alignment}</pre>"

    return render_template("upload.html")

if __name__ == "__main__":
    app.run(debug=True)
