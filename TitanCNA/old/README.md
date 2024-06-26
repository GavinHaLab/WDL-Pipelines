# TITAN


Overview of running TITAN

Description
This workflow will run the TITAN a set of tumour-normal pairs, starting from the BAM files and generating TitanCNA outputs. It will also perform model selection at the end of the workflow to choose the optimal ploidy and clonal cluster solutions.

Contact
Gavin Ha
Fred Hutchinson Cancer Research Center
contact: gavinha@gmail.com or gha@fredhutch.org
Date: May 11, 2019
Website: GavinHaLab.org

Requirements
Software packages or libraries
R-3.5
TitanCNA (v1.15.0+)
TitanCNA imports: GenomicRanges, dplyr, data.table, doMC
ichorCNA (v0.1.0)
HMMcopy
optparse
stringr
SNPchip
Python 3.4
snakemake-3.12.0
PySAM-0.11.2.1
PyYAML-3.12
samtools v1.3.1
bcftools v1.1
HMMcopy Suite.
-In particular, readCounter is used.
Scripts/executables
readCounter (C++ executable; HMMcopy Suite)
ichorCNA.R (ichorCNA tool for normalizing and correcting read coverage)
countPysam.py (generates input allele counts)
titanCNA.R (main R script to run TitanCNA)
selectSolution.R (R script to select optimal solution for each sample)
Tumour-Normal sample list
The list of tumour-normal paired samples should be defined in a YAML file. See config/samples.yaml for an example. Both fields samples and pairings must to be provided. pairings key must match the tumour sample while the value must match the normal sample.

samples:
  tumor_sample_1:  /path/to/bam/tumor.bam
  normal_sample_1:  /path/to/bam/normal.bam


pairings:
tumor_sample_1:  normal_sample_1
snakefiles
ichorCNA.snakefile
getAlleleCounts.snakefile
TitanCNA.snakefile
config.yaml
See below for details about config/config.yaml

## OPTIONAL TODO Description of files for running on cromwell

a. qsub
There are 2 separate files in use for qsub, which are provided as a template: config/cluster_qsub.sh - This file contains other qsub parameters. Note that these settings are used for the Broad's UGER cluster so users will need to modify this for their own clusters.
config/cluster_qsub.yaml - This file contains the memory, runtime, and number of cores for certain tasks.

To invoke the snakemake pipeline for qsub:

b. slurm
There is only one file in use for slurm: config/cluster_slurm.yaml - This file contains the memory, runtime, and number of cores for certain tasks. To invoke the snakemake pipeline for qsub:


1. Path to tools
samTools is used by (getAlleleCounts.snakefile)[getAlleleCounts.snakefile].

samTools:  /path/to/samtools
2. Path to scripts
These are provided in this (TitanCNA)[https://github.com/gavinha/TitanCNA] repo either under code/ or R_scripts.

readCounterScript:  /path/to/readCounter
ichorCNA_rscript:  /path/to/ichorCNA.R
pyCountScript:  code/countPysam.py
TitanCNA_rscript: ../R_scripts/titanCNA.R
TitanCNA_combineTitanIchorCNA:  code/combineTITAN-ichor.R
TitanCNA_selectSolutionRscript: ../R_scripts/selectSolution.R
See the ichorCNA repo for the ichorCNA R script.

3. Path to R package files
Specify the directory in which TitanCNA and ichorCNA are installed.
Set these if the R files in these libraries have been modified or updated but not yet installed or updated in R.

TitanCNA_libdir:  ../../R/
ichorCNA_libdir:  /path/to/ichorCNA/ ## optional
4. Reference files and settings
Global reference files used by many of the snakefiles and scripts.

snpVCF you can download the HapMap file (used for filtering heterozygous SNPs) here: https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0
genomeStyle specifies the chromosome naming convention to used for output files. Input files can be any convention as long as it is the same genome build. Only use UCSC (e.g. chr1) or NCBI (e.g. 1).
sex set to male or female, otherwise None if both females and males are in sample set.
cytobandFile Only used for hg38! Coordinates for plotting the idiogram. The file is found in the (ichorCNA repo)[https://github.com/broadinstitute/ichorCNA/blob/master/inst/extdata/cytoBand_hg38.txt].
genomeBuild: hg38 # Use "None" if hg19
genomeStyle:  UCSC
refFasta: /path/to/Homo_sapiens_assembly38.fasta
snpVCF:  /path/to/hapmap_3.3.hg38.vcf.gz 
cytobandFile:  data/cytoBand_hg38.txt # only need if hg38
centromere:  /path/to/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
sex:  male   # use "None" if both females and males are in sample set
5. readCounter settings
Settings for the computing read coverage using readCounter.

chrs Used by readCounter command line tool. Please comment out the chromosome naming convention that you do not wish to use.
binSize Various bin sizes can be used in the analysis.
## readCounter params ##
# use this for NCBI chr naming
chrs: 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y 
# use this for UCSC chr naming
#chrs: chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY 
binSize:  10000
6. ichorCNA.snakefile settings
Settings for the analysis of read coverage using ichorCNA.

ichorCNA_chrs specifies the chromosomes to analyze in the R script; users do not need to be concerned about chromosome naming convention here as the code will handle it based on the genomeStyle set in the reference settings above.
The GC and Map wig files must have match the binSize above. The various sizes supported will depend on which GC and Map wig files are generated. We provide a few to select from in the ichorCNA repo. Users can also create their own by following instructions from the HMMcopy repo.
## ichorCNA params ##
ichorCNA_gcWig: /path/to/gc_hg19_10kb.wig
ichorCNA_mapWig:  /path/to/map_hg19_10kb.wig
ichorCNA_chrs:  c(1:22, \"X\")
ichorCNA_normal:  c(0.5)  
ichorCNA_ploidy:  c(2,3)  
ichorCNA_estimateNormal:  TRUE
ichorCNA_estimatePloidy:  TRUE
ichorCNA_estimateClonality: TRUE
ichorCNA_scStates:  c(1,3)
ichorCNA_maxCN:  8
ichorCNA_includeHOMD: FALSE
ichorCNA_txnE:  0.9999
ichorCNA_txnStrength:  10000
ichorCNA_plotFileType:  png
ichorCNA_plotYlim:  c(-2,4) 
7. getAlleleCounts.snakefile settings: Tumor allelic counts
Minimum thresholds used when extracting read counts from the tumor BAM file at heterozygous SNP sites. In addition, the vcf_quality of the heterozygous SNP site from the matched normal sample is also used in the filtering. Note: Users must modify getAlleleCounts.snakefile to set the appropriate chromosome naming convention for samtools and the countPysam.py script. Comment out line 6 if using UCSC or comment out line 8 when using NCBI.

# USERS MUST MODIFY getAlleleCounts.snakefile to use the correct CHRS naming
map_quality:  10
base_quality: 10
vcf_quality:  100
8. TitanCNA.snakefile settings
Most settings can be left as default.

TitanCNA_maxNumClonalClusters specifies the maximum number of clonal clusters to consider. For example, if set to 5, then 5 solutions are generated, each one considering a different number of cluster(s).
TitanCNA_maxPloidy specifies the maximum ploidy to initialize. This be set to either 2 (only considers diploid solutions), 3 (considers diploid and triploid, and usually accounts for tetraploid), or 4 (for diploid, triploid, tetraploid or higher ploidies). Usually, 3 is suitable for most samples unless you know that your samples are tetraploid or even higher. For example, if set to 3, then solutions for diploid and triploid will be generated. code/selectSolution.R will try to select the optimal solution; however, users should inspect to make sure results are accurate.
TitanCNA_numCores specifies the number of cores to use on a single machine. Ff using a cluster, then must match the settings for number of cpus in config/cluster_qsub.yaml or config/cluster_slurm.yaml.
TitanCNA_maxNumClonalClusters: 2
TitanCNA_chrs:  c(1:22, \"X\")
TitanCNA_normalInit: 0.5
TitanCNA_maxPloidy: 3
TitanCNA_estimateNormal:  map
TitanCNA_estimatePloidy:  TRUE
TitanCNA_estimateClonality: TRUE
TitanCNA_alleleModel: binomial
TitanCNA_alphaK:  10000
TitanCNA_alphaR:  10000
TitanCNA_txnExpLen: 1e15
TitanCNA_plotYlim:  c(-2,4)
TitanCNA_solutionThreshold: 0.05
TitanCNA_numCores: 1 
