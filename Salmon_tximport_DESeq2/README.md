# RNA-seq Analysis Project

## Overview
In this project, I have performed differential gene expression analysis using a workflow based on Salmon, tximport, and DESeq2. 
This analysis aims to identify genes that are differentially expressed between two groups of date palm varieties based on their 
sucrose content. 

I have used a "pseudo-alignment" based method, which eliminates the need for read alignment and the creation of SAM/BAM files. 
Instead, I employed Salmon to estimate transcript abundances as Transcripts Per Million (TPMs), which were then converted to 
gene-level counts for differential gene expression (DGE) analysis using DESeq2.

## Analysis Workflow
### Task 1: Running Salmon
In this task, I used Salmon to estimate transcript abundances for the RNA-seq data.
Salmon requires a reference transcriptome fasta file for pseudo-alignment.
Index files were created for the reference transcriptome before running Salmon.
### Task 2: Using tximport
tximport is used to import the Salmon output, including TPMs, into DESeq2.
A mapping file (tx2gene) was used to convert TPMs to gene-level counts.
A DESeq2DataSet was created for downstream analysis.
### Task 3: Running DESeq2
DESeq2 was used to perform differential gene expression analysis.
Genes with fewer than 10 reads were removed from the analysis.
Log fold-change shrinkage was applied.
The results were ordered based on adjusted p-values (FDR).

## Results and Interpretation
### Task 3: DESeq2 Results
The DESeq2 analysis provided a list of differentially expressed genes between high and low sucrose content date palm varieties.
The results include log fold-changes, p-values, and adjusted p-values.
Further interpretation and visualization of the results were performed.

### Questions and Discussion
The project includes questions for discussion and interpretation of the results, such as the distribution of p-values, differentially expressed genes, and dispersion-mean plots.
These discussions help in gaining insights into the biological significance of the findings.
### Conclusion
This RNA-seq analysis project demonstrates the use of a Salmon + tximport + DESeq2 workflow for differential gene expression analysis. 

