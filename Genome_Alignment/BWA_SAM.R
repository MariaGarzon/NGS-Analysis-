Purpose: Aligned short reads to the human reference genome and used samtools to process SAM/BAM short read alignments. 
This is part of multi-week “re-sequencing” workflow where the aim is to call single nucleotide polymorphisms (SNPs) from Illumina 
short read sequencing of human genomes.

# 1. Create a FASTA index and a bwa index for the normalized reference genome
Many NGS software require that input files be indexed prior to including them in work flows. 
Index files permit rapid lookups of genome coordinates in large files and dramatically speed up computation times.
The re-sequencing workflow we will run to align reads from samples in the 1,000 genomes project followed by snp-calling requires 
the following index files:
  
  a. a FASTA index
b. a set of index files required by the BWA-MEM aligner
c. a dictionary file required by the Genome Analysis Toolkit (GATK)

### Create a directory and copy to it the hg38 reference genome FASTA as follows:
cd $SCRATCH
mkdir hg38
cd hg38
cp /scratch/work/courses/BI7653/hw3.2023/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa .

### Create a slurm job script that will create the FASTA index files in the same directory as the reference genome FASTA. 

#!/bin/bash 
#
#SBATCH –nodes=1 
#SBATCH –tasks-per-node=1 
#SBATCH –cpus-per-task=1 
#SBATCH –time=5:00:00 
#SBATCH –mem=32GB 
#SBATCH –job-name=slurm_template 
#SBATCH –mail-type=FAIL 
#SBATCH –mail-user=mpg8596@nyu.edu

module purge 
# Load the most recent versions of the samtools and bwa modules
module load samtools/intel/1.14 module load bwa/intel/0.7.17

# Navigate to the directory with the reference fasta file
cd /scratch/mpg8596/hg38

# Run samtools faidx
samtools faidx Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa

# Run bwa index
bwa index -a bwtsw Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa

slurm_template.sh (END)

# ls -al in  hg38 directory: 
[mpg8596@ga001 hg38]$ ls -al total 8355588 drwxrwsr-x. 2 mpg8596 mpg8596 4096 Feb 20 10:47 . drwx–S—. 6 mpg8596 
mpg8596 4096 Feb 20 10:14 .. -rw-rwx—. 1 mpg8596 mpg8596 3130750435 Feb 20 08:36 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa -rw-rw-r–. 1 
mpg8596 mpg8596 18172 Feb 20 10:32 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.amb -rw-rw-r–. 1 
mpg8596 mpg8596 7418 Feb 20 10:32 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.ann -rw-rw-r–. 1 
mpg8596 mpg8596 3099750792 Feb 20 10:31 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.bwt -rw-rw-r–. 1 
mpg8596 mpg8596 6793 Feb 20 09:55 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.fai -rw-rw-r–. 1 
mpg8596 mpg8596 774937681 Feb 20 10:32 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.pac -rw-rw-r–. 1 
mpg8596 mpg8596 1549875408 Feb 20 10:47 Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa.sa -rw-rw-r–. 1 
mpg8596 mpg8596 881 Feb 20 08:52 slurm-30375296.out -rw-rw-r–. 1 mpg8596 mpg8596 905 Feb 20 09:12 slurm-30375365.out -rw-rw-r–. 1 
mpg8596 mpg8596 447 Feb 20 09:27 slurm-30375429.out -rw-rw-r–. 1 mpg8596 mpg8596 447 Feb 20 09:31 slurm-30375437.out -rw-rw-r–. 1 
mpg8596 mpg8596 447 Feb 20 09:46 slurm-30375489.out -rw-rw-r–. 1 mpg8596 mpg8596 6107 Feb 20 10:47 slurm-30375512.out -rw-rwx—. 1
mpg8596 mpg8596 637 Feb 20 09:55 slurm_template.sh

### 2. Align paired-end fastq data from 30 human samples to the human reference genomes using BWA-MEM
Purpose: Samples represent 4 populations in the 1000 genomes project.

The script executes a job array that will find the paired-end fastqs to be aligned using a tab-delimited table with columns 
sample name, read 1 fastq file name, and read 2 fastq filename.

#!/bin/bash
#
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=bwamem_array
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=<mpg8596@nyu.edu>
#SBATCH --array=1-30

module purge

# This variable holds the path to the normalized hg18 reference genome FASTA 
ref=/scratch/mpg8596/hg38/Homo_sapiens.GRCh38.dna_sm.primary_assembly.normalized.fa 

# This variable contains the path to a 3-column tab-delimited table (hw3_fastqs.processed.txt) that specifies 
# sample names and the corresponding paired-end FASTQ files for each sample. 
# This table is used to determine which samples to process in the job array.
table=/scratch/work/courses/BI7653/hw3.2023/fastqs.processed/hw3_fastqs.processed.txt
fqdir=/scratch/work/courses/BI7653/hw3.2023/fastqs.processed

# The following code defines sample, fq1 and fq2 variables for current array index
# note: SLURM_ARRAY_TASK_ID environmental variable will contain a single value corresponding to the current array index
# For each array index (sample), the script extracts information (sample name, FASTQ file names) from the sample table using 
# the head, tail, and cut commands based on the current SLURM_ARRAY_TASK_ID.
line="$(head -n $SLURM_ARRAY_TASK_ID $table | tail -n 1)"
sample="$(printf "%s" "${line}" | cut -f1)"
fq1="$(printf "%s" "${line}" | cut -f2)"
fq2="$(printf "%s" "${line}" | cut -f3)"

# Print to standard out the array index and the sample name
echo Processing array index: $SLURM_ARRAY_TASK_ID sample: $sample

# Make a directory for the sample and cd to it
mkdir $sample
cd $sample

# Load the bwa mem module
module load bwa/intel/0.7.17

bwa mem \
-M \
-t $SLURM_CPUS_PER_TASK \
-R "@RG\tID:${sample}.id\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}.lb" \ #Adds a read group header line to the output SAM file, providing sample-specific information.
"${ref}" \ #Specifies the path to the reference genome FASTA file.
$fqdir/$fq1 \ # paths to the input paired-end FASTQ files.
$fqdir/$fq2 > $sample.sam

echo _ESTATUS_ [ bwa mem for $sample ]: $?
  
echo _END_ [ hw3_bwamem.slurm for $sample ]: $(date)

### Report

There are 29 .sam files in total.

The exit statuses of all 30 subjobs are 1, which indicates successful completion of the BWA MEM alignment step for each sample.

# SAM/BAM format and Samtools
 
The @HD tag line has the following format:
  
@HD:The character is a tab delimiter that separates the fields of the header line. 
The VN field indicates the version number of the SAM/BAM format specification. 
SO field indicates the sorting order of the alignments. 
The sorting order is important for downstream analysis tools that require the alignments to be in a specific order, such as for variant calling or de novo assembly.

The @HD tag line is: @HD VN:1.0 SO:coordinate

This indicates that the file is in SAM/BAM format version 1.0 and that the alignments are sorted by genomic coordinate.

# Number of unmapped reads in the BAM: 7247
                                                    
# Number of mapped reads are  in the BAM?
[mpg8596@cs095 ~]$ samtools view -c -F 4 /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam
2924253
# The percentage mapping rate (total mapped reads / total reads in the alignment) for this sample?
total mapped reads: 2924253 total reads: 2931500 % mapping rate = 99.75
                                                    
# A hypothetical SAM file has alignment records with the bitwise flag values that include 4, 147, 113, 99 on the decimal scale. 
# What are the binary and hexadecimal representations of each of the these values?
Decimal Binary Hexadecimal 4 0000 0100 0x0004 147 1001 0011 0x0093 113 0111 0001 0x0071 99 0110 0011 0x0063

# Picard has an online tool for determining the meanings of bitwise flag values such as those in Q3.3:
# Using this tool, what are the characteristics of the reads in each these four flags in Q3.3 (4,147,113,99)?
                                                    
Flag value 4: This indicates that the read is unmapped, meaning it does not align to the reference genome.
                                                    
Flag value 147: This indicates that the read is both paired-end and mapped in a proper pair, meaning that it has a mate that is also mapped and in the expected orientation and distance. The read is on the reverse strand, meaning that it aligns to the reverse complement of the reference genome. The mate is on the forward strand, and it is also mapped.
                                                    
Flag value 113: This indicates that the read is both paired-end and mapped in a proper pair. The read is on the forward strand, and the mate is on the reverse strand.
                                                    
Flag value 99: This indicates that the read is both paired-end and mapped in a proper pair. The read is on the forward strand, and the mate is on the reverse strand. The read is the first in the pair.

SAM specification allows for 3 types of alignment records: primary, secondary and supplementary alignments. 
Depending on the alignment software and command line used, secondary or supplementary alignments may also exist in a SAM/BAM file.
BWA-MEM will typically add both primary and secondary alignments. Note that by using the -M option above we instructed bwa to set all supplementary alignments to secondary (0x100).
                                                    
For Q3.5, use samtools view with appropriate -c, -f, -F options to count the following in the BAM file at /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam. For each answer, provide the number of reads and the command line you used.
                                                    
How many alignments are primary?
[mpg8596@ga001 ~]$ samtools view -c -F 256 /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 2931500
                                                    
How many alignments are secondary?
[mpg8596@ga001 ~]$ samtools view -c -f 0x100 /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 0
                                                    
How many alignments are supplementary?
[mpg8596@ga001 ~]$ samtools view -c -f 2048 /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 0
                                                    
What is the number of reads excluding unmapped reads, supplementary reads, secondary reads and PCR duplicates?
                                                      
[mpg8596@ga001 ~]$ samtools view -F 0x904 -f 0x2 -F 0x400 -c /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 2849394

Use samtools view to subset the BAM from chromosome 20 position 1 to 2000000 (i.e., 2 Mb), while also retaining the header. 
Note that to perform this type of operation, the BAM must be coordinate-sorted, (which it is).
samtools view -c <subsetted bam>
[mpg8596@ga001 ~]$ samtools view -b -h /scratch/work/courses/BI7653/hw3.2023/HG00096.chrom20.ILLUMINA.bwa.GBR.low_coverage.20120522.bam 
20:1-2000000 > subset.bam 
[mpg8596@ga001 ~]$ samtools view -c subset.bam 95338
                                                    
