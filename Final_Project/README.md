
# Introduction

Differential Gene Expression (DGE) analysis through RNA-sequencing (RNA-seq) is a fundamental technique in molecular biology, 

particularly in experimental biomedicine. In this project, we delve into DGE analysis to compare gene expression patterns between 

control breast cancer cell lines and those subjected to RNA interference (RNAi) targeting the NRDE2 gene. 

The objective is to identify genes whose expression is significantly altered when NRDE2 is silenced, providing valuable insights 

into its functional role.


Our dataset comprises RNA-seq libraries from three biological replicates of control cell lines and three replicates of cell lines 

where NRDE2 has been silenced using RNAi. These libraries were generated as single-end (SE) reads on an Illumina NextSeq platform, 

offering an opportunity to investigate the regulatory impact of NRDE2.


This report outlines a systematic analysis workflow using tools such as Salmon, tximport, and DESeq2. We begin by preprocessing raw 

Fastq files, addressing issues like adapter contamination and platform-specific artifacts. Subsequently, we conduct quality control 

assessments using FastQC and MultiQC to ensure data reliability. Then transcript quantification using Salmon, enabling the 

computation of gene-level counts essential for DEG analysis with DESeq2.


We present read statistics, mapping rates, and insights into differentially expressed genes. Accompanied by figures and tables, this 

report serves as a practical guide for those interested in exploring DGE analysis with RNA-seq data.


Our objective is to uncover the biological implications of NRDE2 silencing in breast cancer cell lines, contributing to a deeper 

understanding of gene regulation in this context.


# Methods 


## Step 1: Data Processing and Quality Control:

The following steps were performed for data processing and quality control:
  
### 1. Trimming of FastQ Files:

The raw FastQ files were subjected to trimming using the fastp tool. This tool automatically removed adapters from the single-end 

reads and eliminates polyG sequences introduced during sequencing on the Illumina NextSeq platform. 

The bash script used was designed to process multiple samples in parallel using array indexing. 

Each sample's FastQ file was trimmed and subjected to quality control analysis individually. 

Array Indexing and Variable Assignment:

The script uses SLURM_ARRAY_TASK_ID environmental variable to identify the current array index.

It extracts the corresponding sample name and FastQ file name from the table using the "head", "tail", and "cut" commands.

The sample name is assigned to the variable "sample", and the FastQ file name is assigned to the variable "fq".

#### Directory Setup:

The script creates a directory for the current sample using the "mkdir" command and navigates to it using "cd".

#### Output Filename:

The script defines the output filename for the trimmed FastQ file by removing the file extension and adding ".fP.fastq.gz" (forward 

paired) suffix.

The output filename is assigned to the variable "fq_fastp".

#### Trimming with fastp:

The script runs the fastp tool on the sample FastQ file using the "fastp" command.

It specifies input and output files, sets trimming parameters (length_required, n_base_limit), and generates an HTML report and a 

JSON file for further analysis.

After trimming, the script runs FastQC on the processed RNA-seq read using the "fastqc" command.

The trimmed FastQ file is specified as input.

### 2. Quality Control Assessment with FastQC:

The processed RNA-seq reads from each sample were individually subjected to QC analysis using the FastQC tool providing information 

about various quality metrics, including sequence length distribution, per-base quality scores, GC content and the presence of any 

remaining sequencing artifacts or biases. 

### 3. MultiQC Report Generation:

To summarize the quality control results across all samples, a MultiQC report was generated. 

The MultiQC tool combined the FastQC reports from each sample and provided a consolidated overview of the quality control metrics. 

This facilitated the identification of common quality patterns or issues across the control and NRDE2 knockdown samples.

Overall, the GC content remained consistent before and after filtering.

Per Sequence Quality Scores representing the distribution of quality scores for all sequences all show a the majority of sequences 

having a high-quality scores, overall above a threshold of 20 or 30.

No samples found with any adapter contamination > 0.1%.

6 samples had less than 1% of reads made up of overrepresented sequences.

All 6 reads failed the Per Base Sequence Content metric in the MultiQC report, suggesting that there may be issues with the sequence 

content across different positions in the reads. A failure in the Per Base Sequence Content metric could indicate several potential 

problems such as biases or skewness in the nucleotide composition if the distribution of nucleotides is highly skewed or imbalanced 

at certain positions, it could suggest biases in the library preparation or sequencing process. 

For example, it could be caused by biased PCR amplification or nucleotide misincorporation during sequencing. Also, if the sequencing 

process introduces systematic errors or biases, it can result in irregular patterns in the nucleotide distribution. 

### 4. Downloaded the reference genome: 

I obtained the "Homo_sapiens.GRCh38.cdna.all.fa.gz" reference genome file containing the transcript sequences of the Homo sapiens 

(human) genome.

Unzipped the transcripts fasta and ran Picard tools NormalizeFasta. Specifically the NormalizeFasta tool, to process the transcript 

sequences. The NormalizeFasta tool removes everything after the transcript ID, resulting in a simplified version of the transcript 

sequences.

This step ensures consistency in the transcript ID format and reduces any potential complications during downstream analysis.

### 5. Generated a Salmon index: 

Using the simplified transcript sequences, I created a Salmon index. 

The Salmon index is a pre-built data structure that enables efficient mapping and quantification of RNA-seq reads to the reference 

transcriptome. This step involved running the Salmon command with the following parameters: 

  -t normalized_sequence.fasta specifies the input file containing the normalized transcript sequences 
  
  -i transcripts_index defines the output directory or prefix for the Salmon index 
  
  -k 31 specifies the length of k-mers used for indexing.


### 6. Running Salmon in mapping-based mode using a command appropriate for single-end data.

a. The script defines several variables to store relevant information: 
  
  The table variable holds the path to a 2-column (tab-delimited) table file that contains sample names and corresponding fastq file    names. 
  
  The line variable extracts a specific line from the table file based on the current array index. 
  
  The sample variable stores the sample name extracted from the line, and the fq1 variable stores the corresponding fastq file name.
  
  The fqdir variable holds the path to the directory where the fastq files are located. 
  
  The salmon_index_dir variable stores the directory path of the Salmon index that was generated in a previous step. These paths are 
  
  necessary for providing the input files and output directory for the salmon quantification process.

b. A new directory is created with the name of the current sample. This ensures that the output files from Salmon quantification are 

organized and stored separately for each sample.

#### Changing to the sample directory: 

The script changes the current working directory to the newly created sample directory using the cd command. This ensures that the 

subsequent Salmon command is executed within the correct directory.

c. The Salmon quant command is executed with the specified parameters: 
  
  -i ${salmon_index_dir} specifies the path to the Salmon index directory
  
  -l A indicates that the library type is "A" (unstranded)
  
  -r $fqdir/$fq1 specifies the input fastq file for the sample
  
  --validateMappings enables mapping validation
  
  --gcBias enables GC bias correction
  
  --threads ${SLURM_CPUS_PER_TASK} specifies the number of threads to be used for parallel processing
  
  -o $sample.transcripts_quant defines the output directory or prefix for the Salmon quantification results


#### Salmon inferred the library to be a stranded single-end protocol where the reads from come the reverse strand. 

## Step 2: Differential Expression Analysis

### 7. After converting the Salmon TPMs to gene-level counts using tximport, I conducted differential gene expression (DGE) analysis with DESeq2.

a. Defined the sample names (sample_names) and corresponding conditions (sample_condition) in vectors.

b. Specified the folder path (folder_path) where the quant.sf files are located. Then, created the file paths for each sample by concatenating the folder path with the sample names and the file extension.

c. Read the transcript-to-gene mapping information from the tx2gene.csv file using the read.table function.

d. Used the tximport function to import the transcript quantification data. I provided the file paths, specified the data type as "salmon," and passed the tx2gene mapping information.

e. Created a data frame (samples) to store the sample names and conditions. The row names of the data frame were set to the sample names.

f. The DESeqDataSetFromTximport function was used to create a DESeqDataSet object (ddsTxi). I provided the imported transcript quantification data (txi), the sample metadata (colData), and the design formula (design) specifying the relationship between the conditions.

g. Filtered out genes with low counts by keeping only those genes with a minimum sum of counts across all samples.

h. Performed differential expression analysis using the DESeq function on the filtered DESeqDataSet object (ddsTxi).

i. The results function was applied to the ddsTxi object to extract the differential expression results. 

j. The results were sorted in ascending order based on the p-value, and the top 10 most significant genes were viewed. 

k. To obtain more stable estimates of differential expression, I performed log-fold change shrinkage using the lfcShrink function on 
the ddsTxi object. The coef parameter specified the contrast coefficient, and the shrinkage method was set to 'apeglm'.

l. The shrunken results were ordered based on the adjusted p-value (res.shrunk[order(res.shrunk$pvalue),]).

### The key statistical approaches used are as follows:

#### Differential Expression Analysis: 

DESeq2 performs differential expression analysis by fitting a generalized linear model to the count data. 

The design formula ~ condition specifies the comparison of interest, which is the contrast between the "Treated" and "Control" conditions.

#### Statistical Test: 

DESeq2 employs a Wald test to assess the significance of the differences in gene expression between the two conditions. 

The log-fold change estimates and associated p-values are calculated using this test.

#### Multiple-Test Correction: 

To account for the issue of multiple hypothesis testing, the padj column in the DESeq2 results represents adjusted p-values. 

These p-values are corrected for multiple testing using the Benjamini-Hochberg procedure (FDR correction). By applying the log-fold change shrinkage using the lfcShrink function with the type='apeglm' argument, the fold-change estimates are further stabilized and improved, enhancing the interpretation of the results. The shrunken log fold-change estimates are provided in theres.shrunkOrdered variable.

