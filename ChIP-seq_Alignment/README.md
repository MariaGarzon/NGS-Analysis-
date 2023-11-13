# ChIP-seq Analysis Project

![Alt text](https://github.com/MariaGarzon/NGS-Analysis-/blob/main/ChIP-seq_Alignment/workflow.png)

## Overview
In this project, I analyzed ChIP-seq data, aiming to characterize the locations of the Androgen Receptor (AR) transcription factor binding in prostate cancer tumors. 

## Background 
The primary objective of a transcription factor ChIP-seq experiment is to determine the genomic regions that are bound by a 
specific transcription factor. This information can be used to identify potential target genes regulated by the transcription 
factor and to elucidate its role in gene expression. By analyzing the DNA sequences associated with the binding sites, 
this technique can also provide insight into the DNA sequence motifs recognized by the transcription factor.

## Methods

### Data Retrieval

I retrieved the ChIP-seq data from the Sequence Read Archive (SRA) using the `sra-tools` software. The data consisted of 
single-end (1 x 65 bp) FASTQ files sequenced on an Illumina Hi-Seq sequencer. I obtained ChIP samples for Androgen Receptor (AR) 
from two patients' tumors (P1_AR_DSG and P2_AR_DSG) and an "input" control sample (P_Input_DSG).

### Read Processing

To process the raw reads, I utilized the `fastp` tool. This step involved adapter removal and quality filtering to ensure that 
only high-quality reads were used in subsequent analyses. I set an appropriate minimum length for processed reads, considering 
the 1 x 65 bp sequencing length of the data.

### Alignment to Reference Genome

For aligning the processed reads to the human reference genome (GRCh38), I employed the BWA (Burrows-Wheeler Aligner) tool. 
Using the `bwa mem` command, I aligned the reads, specifying parameters like multi-threading for efficient processing.

### BAM File Preparation

After alignment, I prepared BAM files for downstream analysis. This involved three key steps for each SAM file produced 
in the previous alignment step:

1. **Conversion to BAM Format**: I used `samtools view` to convert SAM files to BAM format while retaining header information.

2. **Coordinate Sorting**: I utilized Picard tools to perform coordinate sorting of reads in the BAM files.

3. **BAM Indexing**: To enable efficient random access to BAM files, I created BAM index files (.bai) for each coordinate-sorted BAM file.

### Filtering Low-Quality Reads

In ChIP-seq analysis, it's common to remove reads with low mapping quality. I employed `samtools view` once again, this time using 
the `-q` option to filter out reads with mapping quality less than 20. This step ensured that only high-quality reads were retained 
for subsequent analyses.

### Peak Calling

The final step involved identifying regions of significant enrichment, known as peaks, using the MACS2 (Model-based Analysis of 
ChIP-Seq) tool. I ran MACS2 separately for each Androgen Receptor ChIP sample, providing the corresponding control sample as input. 
This step helped identify binding sites of the AR transcription factor in prostate cancer tumors.

