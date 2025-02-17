# Genetic Analysis Tools

## Overview
This project provides a suite of genetic analysis tools for:
1. **Sequence Alignment** (using Biopython)
2. **Variant Calling** (using SAMtools & BCFtools)
3. **Gene Expression Analysis** (using DESeq2/edgeR)

The tools include a **web-based interface** (Flask) for easier interaction.

---

## Features
### ✅ Sequence Alignment Tool
- Implements **Needleman-Wunsch** & **Smith-Waterman** algorithms.
- Uses **Biopython** for pairwise alignment.
- Web interface for uploading **FASTA** sequences.
- Visualizes alignments using **Matplotlib**.

### ✅ Variant Calling Tool
- Accepts **BAM/FASTQ** input files.
- Uses **SAMtools** for preprocessing.
- Calls variants with **BCFtools**.
- Saves output in **VCF format**.

### ✅ Gene Expression Analysis Tool
- Processes **RNA-Seq FASTQ** data.
- Normalizes using **TPM/RPKM/FPKM**.
- Performs **differential expression analysis** with **DESeq2**.
- Visualizes data with **volcano plots & heatmaps**.

---

## Installation
### 🔹 Prerequisites
- **Python 3.8+**
- **Flask, Biopython, Matplotlib**
- **SAMtools, BCFtools, GATK (for Variant Calling)**
- **R with DESeq2/edgeR (for Gene Expression Analysis)**

### 🔹 Install Python Dependencies
```bash
pip install flask biopython matplotlib plotly pandas
```

### 🔹 Install System Dependencies
```bash
# Ubuntu/Debian
sudo apt install samtools bcftools

# macOS (Homebrew)
brew install samtools bcftools
```

---

## Usage
### **1️⃣ Sequence Alignment**
Run the Flask web app:
```bash
python app.py
```
Then, open `http://127.0.0.1:5000` in your browser to upload **FASTA files**.

---

### **2️⃣ Variant Calling**
Place your **BAM file** in `variant_calling/`, then run:
```bash
python variant_calling.py sample.bam
```
This will generate `output.vcf` with detected variants.

---

### **3️⃣ Gene Expression Analysis**
Run the analysis script with your RNA-Seq data:
```bash
python gene_expression.py input_counts.csv
```
This will generate **volcano plots** & **heatmaps** of gene expression changes.

---

## File Structure
```
📂 genetic_analysis/
├── 📂 sequence_alignment/
│   ├── app.py  # Flask Web App
│   ├── templates/
│   │   ├── index.html
│   ├── static/
│   │   ├── style.css
├── 📂 variant_calling/
│   ├── variant_calling.py
│   ├── pipeline.sh  # Shell script for SAMtools/BCFtools
├── 📂 gene_expression/
│   ├── gene_expression.py
│   ├── sample_counts.csv  # Sample input
├── README.md
├── requirements.txt
```

---

## Example BAM File Download
You can download a test BAM file using:
```bash
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/wgEncodeUwRepliSeqBg02esG1AlnRep1.bam -O sample.bam
```
Then, run the variant calling tool.

---

## Contributing
Feel free to submit **pull requests** or open an **issue** for improvements.

---

## License
This project is **open-source** under the MIT License.

