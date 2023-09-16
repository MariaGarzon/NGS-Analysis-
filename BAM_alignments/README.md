# Coverage Depth and Copy Number Variant Detection

## Assignment Overview
Use sequencing coverage and the insert size distribution from BAM alignments from paired-end sequencing to predict structural variants.

## Background 
### Deletions: 
In paired-end sequencing, the insert size represents the distance between the two sequenced ends of a DNA fragment. 
When you align these paired-end reads to a reference genome, the expected insert size should match the library preparation protocol.
Deletions result in a decrease in the insert size of the aligned reads. If a portion of the genome is deleted, paired-end reads that 
should span the deleted region will have a smaller insert size than expected. By analyzing the insert size distribution, you can 
identify outliers with significantly smaller insert sizes than the median or mean. These outliers may indicate the presence of 
deletions.
### Insertions:
Conversely, insertions can lead to an increase in the insert size of aligned reads. When an extra piece of DNA is inserted into the 
genome, the paired-end reads that span the insertion site will have a larger insert size than expected.

## Assignment Tasks
Here's a brief summary of the key tasks covered in this assignment:

### Task 1: Summarize Coverage and Insert Size Distribution
Calculated the genome-wide coverage depth for a specific sample (CR2342) using Samtools. 
Delved into the insert size distribution analysis, which is vital for identifying structural variants like deletions in the 
sample genome.

### Task 2: Coverage Depth in Genomic Regions and CNV Discovery
Calculated coverage depth within specific genomic intervals (500 bp intervals on chromosome_1) for two different samples 
(CR407 and CR2342). The ultimate goal was to identify copy number variants (CNVs) based on these coverage depth values.

## Tools and Software
Throughout the assignment, we utilize fundamental bioinformatics tools such as Samtools for processing BAM files and R for data visualization. 


