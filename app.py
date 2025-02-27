from flask import Flask, request, render_template, jsonify, abort
import os
import subprocess
import pandas as pd
import plotly.express as px
from Bio import SeqIO
from sequence_alignment.aligner import align_sequences
from werkzeug.utils import secure_filename
import logging
from pathlib import Path

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

app = Flask(__name__)

# Configuration
class Config:
    UPLOAD_FOLDER = Path("uploads")
    ALLOWED_EXTENSIONS = {
        'fasta': ['.fasta', '.fa'],
        'bam': ['.bam'],
        'csv': ['.csv']
    }
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB max file size

app.config.from_object(Config)

# Create upload directory if it doesn't exist
Config.UPLOAD_FOLDER.mkdir(exist_ok=True)

def allowed_file(filename, file_type):
    return '.' in filename and \
           any(filename.lower().endswith(ext) for ext in Config.ALLOWED_EXTENSIONS[file_type])

@app.errorhandler(Exception)
def handle_error(error):
    logger.error(f"An error occurred: {str(error)}")
    return jsonify({"error": str(error)}), getattr(error, 'code', 500)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/align", methods=["POST"])
def align():
    try:
        if "file" not in request.files:
            abort(400, description="No file provided")
        
        file = request.files["file"]
        if not file or not allowed_file(file.filename, 'fasta'):
            abort(400, description="Invalid file format. Please upload a FASTA file.")

        filename = secure_filename(file.filename)
        filepath = Config.UPLOAD_FOLDER / filename
        file.save(filepath)

        try:
            sequences = list(SeqIO.parse(filepath, "fasta"))
            if len(sequences) < 2:
                raise ValueError("FASTA file must contain at least two sequences")

            alignment = align_sequences(sequences[0].seq, sequences[1].seq, method="global")
            return jsonify({
                "status": "success",
                "alignment": str(alignment)
            })
        finally:
            # Clean up uploaded file
            filepath.unlink(missing_ok=True)

    except Exception as e:
        logger.error(f"Alignment error: {str(e)}")
        return jsonify({"error": str(e)}), 400

@app.route("/variants", methods=["POST"])
def call_variants():
    try:
        if "bam_file" not in request.files:
            abort(400, description="No BAM file provided")
        
        file = request.files["bam_file"]
        if not file or not allowed_file(file.filename, 'bam'):
            abort(400, description="Invalid file format. Please upload a BAM file.")

        filename = secure_filename(file.filename)
        filepath = Config.UPLOAD_FOLDER / filename
        file.save(filepath)

        try:
            result = subprocess.run(
                ["bash", "variant_calling/pipeline.sh", str(filepath)],
                capture_output=True,
                text=True,
                check=True
            )
            return jsonify({
                "status": "success",
                "variants": result.stdout
            })
        finally:
            # Clean up uploaded file
            filepath.unlink(missing_ok=True)

    except subprocess.CalledProcessError as e:
        logger.error(f"Variant calling error: {e.stderr}")
        return jsonify({"error": "Variant calling failed", "details": e.stderr}), 400
    except Exception as e:
        logger.error(f"Variant calling error: {str(e)}")
        return jsonify({"error": str(e)}), 400

@app.route("/gene_expression", methods=["POST"])
def gene_expression():
    try:
        if "file" not in request.files:
            abort(400, description="No file provided")
        
        file = request.files["file"]
        if not file or not allowed_file(file.filename, 'csv'):
            abort(400, description="Invalid file format. Please upload a CSV file.")

        filename = secure_filename(file.filename)
        filepath = Config.UPLOAD_FOLDER / filename
        file.save(filepath)

        try:
            df = pd.read_csv(filepath)
            if df.empty:
                raise ValueError("The CSV file is empty")

            fig = px.imshow(
                df.iloc[:, 1:],
                labels=dict(x="Sample", y="Gene", color="Expression Level"),
                title="Gene Expression Heatmap"
            )
            return jsonify({
                "status": "success",
                "plot": fig.to_json()
            })
        finally:
            # Clean up uploaded file
            filepath.unlink(missing_ok=True)

    except Exception as e:
        logger.error(f"Gene expression analysis error: {str(e)}")
        return jsonify({"error": str(e)}), 400

if __name__ == "__main__":
    app.run(debug=False)  # Set debug=False for production
