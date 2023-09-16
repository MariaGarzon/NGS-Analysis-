# RNA-seq Differential Gene Expression Analysis Assignment

## Overview
This assignment focuses on conducting differential gene expression analysis using the DESeq2 package in R. 
I worked with RNA-seq data from date palm fruit and investigated differential gene expression between date palm varieties with 
high sucrose content and low sucrose content.

## Dataset
Source: RNA-seq data from date palm fruit.

## Objective: 
Determine if a group of linked invertase enzymes, identified by Genome Wide Association Study (GWAS), exhibits differential 
gene expression between date palm varieties with high sucrose content (n=4) and those with low sucrose content (n=4).

## Tasks
### Task 1: Create a DESeqDataSet object
Construct a sample table with sample information.
Read htseq-count files using DESeqDataSetFromHTSeqCount.
Explore the DESeqDataSet object.
### Task 2: Pre-filter low count genes, normalize counts, and run DESeq
Remove genes with 10 or fewer reads.
Run DESeq to estimate size factors, dispersions, and conduct a Wald Test for DGE.
Summarize the DESeqDataSet object.
### Task 3: Analyze DGE Data
Perform hierarchical clustering and PCA on normalized counts.
Interpret clustering and PCA results.
Investigate the presence of batch effects.
Extract results of differential expression analysis.
Report log2 fold-changes for candidate genes.
Report normalized counts for candidate genes.
