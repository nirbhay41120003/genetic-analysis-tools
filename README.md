1️⃣ Sequence Alignment Tool
✔ Implements Needleman-Wunsch (Global) and Smith-Waterman (Local) using Biopython.
✔ Upload FASTA files for pairwise alignment.
✔ Displays results in text format (can be extended with Matplotlib/D3.js).

How It Works
Upload a FASTA file containing two sequences.
The backend extracts sequences and aligns them using Biopython.PairwiseAligner.
Displays the aligned sequences on the webpage.
🔹 Try it:

Run python app.py
Open http://127.0.0.1:5000/
Upload a FASTA file with two sequences.
2️⃣ Variant Calling Tool
✔ Uses SAMtools for preprocessing BAM files.
✔ Calls variants using a basic pipeline (Can be extended with GATK/FreeBayes).
✔ Displays detected variants.

How It Works
Upload a BAM file.
The pipeline script (variant_calling/pipeline.sh) runs:
samtools mpileup to generate an intermediate variant file.
bcftools call -mv to call variants and generate a VCF file.
The results are displayed in plain text (Can be visualized using Plotly).
🔹 Try it:

Upload a BAM file via the web interface.
The backend processes it and outputs variant calls.
3️⃣ Gene Expression Analysis Tool
✔ Accepts RNA-Seq data (CSV format for now, can be expanded to FASTQ).
✔ Normalizes gene expression using TPM/RPKM/FPKM.
✔ Generates volcano plots & heatmaps using Seaborn.

How It Works
Upload an RNA-Seq gene expression dataset (CSV format).
The script processes the data and normalizes it using TPM/RPKM.
A heatmap is displayed using Seaborn.
🔹 Try it:

Run python app.py
Click "Generate Heatmap" in the web UI.
🛠 Improvements & Enhancements
✅ To Add BLAST+ Support (for Advanced Alignment):
Modify aligner.py to use blastn (NCBI BLAST+) for more complex alignments.

✅ To Integrate GATK/FreeBayes (for Variant Calling):
Update pipeline.sh to use:

bash
Copy
Edit
gatk HaplotypeCaller -R reference.fasta -I sample.bam -O output.vcf
✅ To Add Differential Expression Analysis:
Use DESeq2 in R for better gene expression comparison.