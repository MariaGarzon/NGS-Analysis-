## Genomic Variant Calling and Filtering

### Description:
This project involves processing genomic data from the 1000 Genomes Project, specifically focusing on single nucleotide polymorphism 
(SNP) calling and filtering. The goal is to extract high-quality SNPs from a multi-sample Variant Call Format (VCF) file using 
the Genome Analysis Toolkit (GATK) pipeline. 

### The project is divided into three main tasks:

#### Task 1: Call SNPs and Genotypes with GenotypeGVCFs
The GATK GenotypeGVCFs tool is employed to call SNPs and genotypes from a combined multi-sample .gvcf file generated 
from previous alignment and processing steps. 

#### Task 2: Subset SNPs from VCF
Subsetting the VCF file to extract SNPs only using GATK SelectVariants 

#### Task 3: Hard Filtering SNPs
Hard filtering is applied to the SNPs obtained from Task 2. Although not as comprehensive as Variant Quality Score Recalibration (VQSR), 
hard filtering helps eliminate outlier SNPs by setting thresholds on various criteria such as read depth and strand bias. 
