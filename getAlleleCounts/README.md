# getAlleleCounts.wdl
This was originally made for TitanCNA, but can be run standalone. I've also provided an input json for if you only want to run getAlleleCounts.wdl on any samples.

## Inputs:
For each sample:

"`sampleName`"  
"`sex`"- use "None" if both male and female are in sample      
"`tumorBam`"- Path to tumor bam file.  
"`tumorBai`"- Path to tumor bam.bai file- if there is none, leave it as "".  
"`normalBam`"- Path to normal bam file  
"`normalBai`"- Path to normal bam.bai file. Can be null  
"`normalPanel`"- Median corrected depth from panel of normals. Default: null.  
"`genomeBuild`"- hg19 or hg38 only, capitalization matters  
"`genomeStyle`"- "NCBI" when hg19 or "UCSC" when hg38 

Minimum thresholds used when extracting read counts from the tumor BAM file at heterozygous SNP sites. In addition, the `vcf_quality` of the heterozygous SNP site from the matched normal sample is also used in the filtering.  
"`getAlleleCounts.mapQ`"- Integer, i.e. 10  
"`getAlleleCounts.baseQ`"- Integer, i.e. 10  
"`getAlleleCounts.vcfQ`"- Integer, i.e. 100  

Files:  
"`getAlleleCounts.refFasta`"- "path/to/Homo_sapiens_assembly38.fasta"  
"`getAlleleCounts.snpVCF`"- "path/to/hapmap_3.3.hg38.vcf.gz"  
"`getAlleleCounts.countScript`"- "https://raw.githubusercontent.com/gavinha/TitanCNA/master/scripts/snakemake/code/countPysam.py"  

Tool commands:  
"`getAlleleCounts.samtools`"- "samtools"  
"`getAlleleCounts.bcftools`"- "bcftools"

## Outputs:
hetSites = "{tumorName}.chr{chr}.vcf" for each chr in each tumor  
alleleCounts = "{tumorName}.tumCounts.chr{chr}.txt" for each chr in each tumor  
concatenatedCounts = "{tumorName}.tumCounts.txt" for each tumor  
