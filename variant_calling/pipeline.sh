#!/bin/bash

# Index the reference genome
bwa index reference.fasta

# Align reads to reference genome
bwa mem reference.fasta input.fastq > aligned.sam

# Convert SAM to BAM, sort, and index
samtools view -S -b aligned.sam > aligned.bam
samtools sort aligned.bam -o sorted.bam
samtools index sorted.bam

# Call variants using FreeBayes
freebayes -f reference.fasta sorted.bam > variants.vcf
