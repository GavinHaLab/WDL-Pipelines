# TitanCNA
`TitanCNA.wdl` is a workflow that uses ichorCNA.wdl and getAlleleCounts.wdl before running TitanCNA.

Source code [here](https://github.com/gavinha/TitanCNA/tree/master), where complete info on TitanCNA is held.

## Inputs- wip

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

"`TitanCNA.ichorCNA_exons`"- Path to bed file containing exon regions. Default: "NULL"  
"`TitanCNA.ichorCNA_binSize`"- "10kb" This must match binSizeNumeric  
"`TitanCNA.ichorCNA_binSizeNumeric`"- 10000, but must match the other binSize   
"`TitanCNA.ichorCNA_qual`"- Integer, i.e. 20  
"`TitanCNA.ichorCNA_normal`"- Initial normal contamination; can be more than one value if additional normal initializations are desired. Default: "0.5"    
"`TitanCNA.ichorCNA_ploidy`"- Initial tumour ploidy; can be more than one value if additional ploidy initializations are desired. Default: "2"  
"`TitanCNA.ichorCNA_estimateNormal`"- Estimate normal. These need to be quoted strings in all caps instead of Boolean b/c of R. Default: "TRUE"   
"`TitanCNA.ichorCNA_estimatePloidy`"- Estimate tumour ploidy. Default: "TRUE"     
"`TitanCNA.ichorCNA_estimateClonality`"- Estimate clonality. Default: "TRUE"      
"`TitanCNA.ichorCNA_scStates`"- Subclonal states to consider. Default: "NULL"   
"`TitanCNA.ichorCNA_maxCN`"- Total clonal CN states. Default: "7"     
"`TitanCNA.ichorCNA_scPenalty`"- penalize subclonal events - n-fold multiplier; n=1 for no penalty     
"`TitanCNA.ichorCNA_includeHOMD`"- If FALSE, then exclude HOMD state. Useful when using large bins (e.g. 1Mb). Default: "FALSE" 
"`TitanCNA.ichorCNA_plotFileType`"- File format for output plots. pdf or png. Default: "pdf"   
"`TitanCNA.ichorCNA_plotYlim`"- "ylim to use for chromosome plots. Default: "c(-2,2)"  
"`TitanCNA.ichorCNA_likModel`"- Default: "t". if multisample, use "gauss"    
"`TitanCNA.ichorCNA_minMapScore`"- control segmentation - higher (e.g. 0.9999999) leads to higher specificity and fewer segments     
"`TitanCNA.ichorCNA_maxFracGenomeSubclone`"- Exclude solutions with subclonal genome fraction greater than this value. Default: "0.5"   
"`TitanCNA.ichorCNA_maxFracCNASubclone`"- Exclude solutions with fraction of subclonal events greater than this value. Default: "0.7"  
"`TitanCNA.ichorCNA_normal2IgnoreSC`"- Ignore subclonal analysis when initial normal setting >= this value    
"`TitanCNA.ichorCNA_txnE`"- Self-transition probability. Increase to decrease number of segments. Lower (e.g. 0.99) leads to higher sensitivity and more segments. Default: "0.9999999". Higher (e.g. 10000000) leads to higher specificity and fewer segments and lower (e.g. 100) leads to higher sensitivity and more segments     
"`TitanCNA.ichorCNA_txnStrength`"- Transition pseudo-counts. Exponent should be the same as the number of decimal places of txnE. Default: "1e+07"    
"`TitanCNA.ichorCNA_fracReadsInChrYForMale`"- Threshold for fraction of reads in chrY to assign as male. Default: "0.001"

Minimum thresholds used when extracting read counts from the tumor BAM file at heterozygous SNP sites. In addition, the `vcf_quality` of the heterozygous SNP site from the matched normal sample is also used in the filtering.  
"`TitanCNA.getAlleleCounts_mapQuality`"- Integer, i.e. 10  
"`TitanCNA.getAlleleCounts_baseQuality`"- Integer, i.e. 10  
"`TitanCNA.getAlleleCounts_vcfQuality`"- Integer, i.e. 100  

"`TitanCNA.maxClusters`"- specifies the maximum number of clonal clusters to consider. For example, if set to 5, then 5 solutions are generated, each one considering a different number of cluster(s).  
"`TitanCNA.maxPloidy`"- specifies the maximum ploidy to initialize. This be set to either 2 (only considers diploid solutions), 3 (considers diploid and triploid, and usually accounts for tetraploid), or 4 (for diploid, triploid, tetraploid or higher ploidies). Usually, 3 is suitable for most samples unless you know that your samples are tetraploid or even higher. For example, if set to 3, then solutions for diploid and triploid will be generated. code/selectSolution.R will try to select the optimal solution; however, users should inspect to make sure results are accurate.  
"`TitanCNA.numCores`"- specifies the number of cores to use on a single machine.  

"`TitanCNA.alphaK`"- Integer, i.e. 10000  
"`TitanCNA.alphaR`"- Integer, i.e. 10000  
"`TitanCNA.txnExpLen`"- Formatted like: 1e15  
"`TitanCNA.mergeIchorHOMD`"- "FALSE", consider setting to "TRUE" when working with pure tumor   
"`TitanCNA.normalInit`"- Float, i.e. 0.5  
"`TitanCNA.estimateNormal`"- "map",  
"`TitanCNA.estimatePloidy`"- "TRUE",  
"`TitanCNA.estimateClonality`"- "FALSE",  
"`TitanCNA.alleleModel`"- "binomial",  
"`TitanCNA.plotYlim`"- "ylim to use for chromosome plots. Default: "c(-2,2)"   
"`TitanCNA.threshold`"- 0.05,  
"`TitanCNA.outputDirectory`"- "path/to/outdir",

"`TitanCNA.refFasta`"- "path/to/Homo_sapiens_assembly38.fasta",     
"`TitanCNA.snpVCF`"- "path/to/hapmap_3.3.b37.vcf.gz",   
"`TitanCNA.cytobandFile`"- If hg38, then: "https://raw.githubusercontent.com/GavinHaLab/ichorCNA/master/inst/extdata/cytoBand_hg38.txt". If hg19, then "None"   
"`TitanCNA.centromere`"- "https://raw.githubusercontent.com/GavinHaLab/ichorCNA/master/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt"