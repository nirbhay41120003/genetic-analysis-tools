#!/bin/bash
echo "Running Variant Calling Pipeline..."
samtools mpileup -f reference.fasta sample.bam | bcftools call -mv -Ov > variants.vcf
echo "Variant calling complete. Output saved in variants.vcf"


chmod +x variant_calling/pipeline.sh
