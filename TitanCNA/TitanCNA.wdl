version 1.0

import "https://raw.githubusercontent.com/argage/wdl-repo/main/getAlleleCounts/getAlleleCounts.wdl" as getAlleleCounts
import "https://raw.githubusercontent.com/argage/wdl-repo/main/ichorCNA/ichorCNA.wdl" as ichorCNA

struct sampleData {
  ## All sample information including any paired normal bams, or a normal panel, and genome build information for individual bams
  String sampleName 
  String sex 
  File tumorBam 
  File tumorBai 
  File? normalBam 
  File? normalBai
  String? normalPanel  # this is listed as a string for a path when it's a file found in the ichor container, but if it's an external file, it needs to be changed to a type File?
  String genomeBuild # hg38 or hg19, only, capitalization matters
  String genomeStyle  #"NCBI" # or NCBI, only
}

workflow TitanCNA {
    input {
        Array[sampleData] batchSamples

        # ichor inputs
        ## Batch level params
        File? ichorCNA_exons  
        Int ichorCNA_binSizeNumeric ## 10000, but must match the other binSize below
        String ichorCNA_binSize  #  "10kb" This must match binSizeNumeric!!!
        Int ichorCNA_qual 
        String ichorCNA_normal 
        String ichorCNA_ploidy 
        String ichorCNA_estimateNormal  # these need to be quoted strings in all caps instead of Boolean b/c of R
        String ichorCNA_estimatePloidy 
        String ichorCNA_estimateClonality 
        String ichorCNA_scStates  # states to use for subclonal CN
        Int ichorCNA_maxCN 
        Float ichorCNA_scPenalty  # penalize subclonal events - n-fold multiplier; n=1 for no penalty,  
        String ichorCNA_includeHOMD 
        String ichorCNA_plotFileType # "pdf" # "png"
        String ichorCNA_plotYlim 
        String ichorCNA_likModel  # "t" # if multisample, use "gauss"
        Float ichorCNA_minMapScore  # control segmentation - higher (e.g. 0.9999999) leads to higher specificity and fewer segments
        Float ichorCNA_maxFracGenomeSubclone 
        Float ichorCNA_maxFracCNASubclone 
        Float ichorCNA_normal2IgnoreSC # Ignore subclonal analysis when initial normal setting >= this value
        Float ichorCNA_txnE # lower (e.g. 0.99) leads to higher sensitivity and more segments
        Float ichorCNA_txnStrength  # control segmentation - higher (e.g. 10000000) leads to higher specificity and fewer segments
        # lower (e.g. 100) leads to higher sensitivity and more segments
        Float ichorCNA_fracReadsInChrYForMale 
        
        # allele inputs
        Int getAlleleCounts_mapQuality
        Int getAlleleCounts_baseQuality
        Int getAlleleCounts_vcfQuality

        # titan inputs
        Int maxClusters
        Int maxPloidy
        Int numCores
        Float normalInit
        String estimatePloidy
        String estimateClonality
        String estimateNormal
        Float alphaK
        Float? alphaR
        Float txnExpLen
        String plotYlim
        String mergeIchorHOMD
        Float threshold
        String outputDirectory
        String? alleleModel

        #path inputs
        String refFasta #https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references
        String snpVCF
        String cytobandFile
        String centromere
    }
    # this is only changed if docker is changed
    # in wdl instead of inputs as it should be an uncommon change
    
    String titan_docker = "argage/titancna:v1.23.1"
    
    String samTools = "samtools"
    String bcfTools = "bcftools"

    String getAlleleCounts_pyCountScript = "/root/scripts/snakemake/code/countPysam.py"#"https://raw.githubusercontent.com/gavinha/TitanCNA/master/scripts/snakemake/code/countPysam.py"
    String rscript = "/root/scripts/R_scripts/titanCNA.R"
    String combineTitanIchorCNA = "/root/scripts/snakemake/code/combineTITAN-ichor.R"
    String selectSolutionRscript = "/root/scripts/R_scripts/selectSolution.R"
    String libdir = "/root/R/"

    

    Array[String] ucscChrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                           "chr22", "chrX", "chrY"]

    Array[String] ncbiChrs = ["1", "2", "3", "4", "5", "6", "7", 
                          "8", "9", "10", "11", "12", "13", "14", 
                          "15", "16", "17", "18", "19", "20", "21", 
                          "22", "X", "Y"]

    # call ichorcna
    # scatter included in ichor
    call ichorCNA.ichorCNA {
        input:
            batchSamples = batchSamples,
            exons = ichorCNA_exons, 
            binSize = ichorCNA_binSize, 
            binSizeNumeric = ichorCNA_binSizeNumeric, 
            qual = ichorCNA_qual, 
            normal = ichorCNA_normal, 
            ploidy = ichorCNA_ploidy,
            estimateNormal = ichorCNA_estimateNormal,
            estimatePloidy = ichorCNA_estimatePloidy,
            estimateClonality = ichorCNA_estimateClonality,
            scStates = ichorCNA_scStates,
            maxCN = ichorCNA_maxCN,
            scPenalty = ichorCNA_scPenalty,
            includeHOMD = ichorCNA_includeHOMD,
            plotFileType = ichorCNA_plotFileType,
            plotYlim = ichorCNA_plotYlim,
            likModel = ichorCNA_likModel,
            minMapScore = ichorCNA_minMapScore,
            maxFracGenomeSubclone = ichorCNA_maxFracGenomeSubclone,
            maxFracCNASubclone = ichorCNA_maxFracCNASubclone,
            normal2IgnoreSC = ichorCNA_normal2IgnoreSC,
            txnE = ichorCNA_txnE,
            txnStrength = ichorCNA_txnStrength,
            fracReadsInChrYForMale = ichorCNA_fracReadsInChrYForMale
    }
    # call getAlleleCounts subworkflow
    # scatter included in allele
    call getAlleleCounts.getAlleleCounts {
        input:
            tumors = batchSamples,
            refFasta = refFasta,
            snpVCF = snpVCF,
            samtools = samTools,
            bcftools = bcfTools,
            countScript = getAlleleCounts_pyCountScript,
            baseQ = getAlleleCounts_baseQuality,
            mapQ = getAlleleCounts_mapQuality,
            vcfQ = getAlleleCounts_vcfQuality
    }

    # Call runTitanCNA for each tumor
    scatter (tumor in batchSamples) {
        Array[String] chrs = if tumor.genomeStyle == "NCBI" then ncbiChrs else ucscChrs

        # Assigning the tumor specific file

        ###### Workflow error:
        # [1] "Failed to process workflow definition 'TitanCNA' (reason 1 of 2): Failed to process 'call runTitanCNA' (reason 1 of 2): 
        # Failed to process expression 'ichorCNA.corrDepth[tumor]' (reason 1 of 1): Type evaluation failed for IndexAccess(IdentifierMemberAccess(ichorCNA,corrDepth,Vector()),IdentifierLookup(tumor)).
        # Cannot dereference a Array[File] value using a WomCompositeType 
        # {\n normalBam -> File?\ngenomeStyle -> String\ngenomeBuild -> String\nsampleName -> String\nnormalPanel -> String?\ntumorBai -> File\nsex -> String\ntumorBam -> File\nnormalBai -> File? \n} key"

        File alleleCountsFile = getAlleleCounts.concatenatedCounts[tumor]
        File corrDepthFile = ichorCNA.corrDepth[tumor]
        File ichorSegFile = ichorCNA.segTxt[tumor]
        File ichorBinFile = ichorCNA.seg[tumor]
        File ichorParamFile = ichorCNA.params[tumor]
        
        call runTitanCNA {
            input:
                alleleCounts = getAlleleCounts.concatenatedCounts[tumor],
                corrDepth = ichorCNA.corrDepth[tumor],
                tumorName = tumor.sampleName, 
                clustNum = maxClusters, 
                ploidy = maxPloidy,
                titanRscript = rscript,
                libdir = libdir,
                numCores = numCores,
                normal = normalInit,
                chrs = chrs,
                sex = tumor.sex,
                genomeStyle = tumor.genomeStyle,
                genomeBuild = tumor.genomeBuild,
                cytobandFile = cytobandFile,
                estimatePloidy = estimatePloidy,
                estimateClonality = estimateClonality,
                estimateNormal = estimateNormal,
                centromere = centromere,
                alphaK = alphaK,
                txnExpLen = txnExpLen,
                plotYlim = plotYlim,
                titan_docker = titan_docker
        }

        # Call combineTitanAndIchorCNA
        call combineTitanAndIchorCNA {
            input:
                titanSeg = runTitanCNA.segTxt,
                titanBin = runTitanCNA.titanOutput,
                titanParam = runTitanCNA.param,
                ichorSeg = ichorCNA.segTxt[tumor], 
                ichorBin = ichorCNA.seg[tumor], 
                ichorParam = ichorCNA.params[tumor],  
                combineScript = combineTitanIchorCNA,
                libdir = libdir,
                mergeIchorHOMD = mergeIchorHOMD,
                sex = tumor.sex,
                centromere = centromere,
                log = "combineTitanAndIchorCNA.log",
                titan_docker = titan_docker
        }
    }
    # Call selectSolution
    call selectSolution {
        input:
            resultFiles = runTitanCNA.titanOutput,
            solutionRscript = selectSolutionRscript,
            threshold = threshold,
            log = "selectSolution.log",
            titan_docker = titan_docker
    }

    # Call copyOptSolution
    call copyOptSolution {
        input:
            optimalClusterSolutionFile = selectSolution.optimalClusterSolution,
            outputDirectory = outputDirectory,
            log = "copyOptSolution.log"
    }

    # Define workflow outputs here
    output {
        Array[File] combinedSegFile = combineTitanAndIchorCNA.segFile
        Array[File] combinedBinFile = combineTitanAndIchorCNA.binFile
        File optimalClusterSolution = selectSolution.optimalClusterSolution
        String copiedSolution = copyOptSolution.copiedSolution
    }
}

task runTitanCNA {
    input {
        File alleleCounts
        File corrDepth
        String tumorName
        Int clustNum
        Int ploidy
        String titanRscript
        String libdir
        Int numCores
        Float normal
        Array[String] chrs
        String sex
        String genomeStyle
        String genomeBuild
        String cytobandFile
        Boolean estimatePloidy
        Boolean estimateClonality
        Boolean estimateNormal
        String centromere
        Float alphaK
        Float txnExpLen
        String plotYlim
        String titan_docker
    }
    command <<<
        Rscript ~{titanRscript} --hetFile ~{alleleCounts} --cnFile ~{corrDepth} \
        --outFile results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.titan.txt \
        --outSeg results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.segs.txt \
        --outParam results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.params.txt \
        --outIGV results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.seg \
        --outPlotDir results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}/ \
        --libdir ~{libdir} --id ~{tumorName} --numClusters ~{clustNum} --numCores ~{numCores} \
        --normal_0 ~{normal} --ploidy_0 ~{ploidy} --genomeStyle ~{genomeStyle} \
        --genomeBuild ~{genomeBuild} --cytobandFile ~{cytobandFile} --chrs "~{chrs}" \
        --sex ~{sex} --estimatePloidy ~{estimatePloidy} --estimateClonality ~{estimateClonality} \
        --estimateNormal ~{estimateNormal} --centromere ~{centromere} --alphaK ~{alphaK} \
        --txnExpLen ~{txnExpLen} --plotYlim ~{plotYlim}
    >>>
    runtime {
        docker: titan_docker
        }
    output {
        File titanOutput = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.titan.txt"
        File segTxt = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.segs.txt"
        File param = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.params.txt"
        File seg = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.seg"
    }
}

task combineTitanAndIchorCNA {
    input {
        File titanSeg
        File titanBin
        File titanParam
        File ichorSeg
        File ichorBin
        File ichorParam
        String combineScript
        String libdir
        Boolean mergeIchorHOMD
        String sex
        String centromere
        String log
        String titan_docker
    }
    command <<<
        Rscript ~{combineScript} --libdir ~{libdir} --titanSeg ~{titanSeg} --titanBin ~{titanBin} \
        --titanParam ~{titanParam} --ichorSeg ~{ichorSeg} --ichorBin ~{ichorBin} --ichorParam ~{ichorParam} \
        --mergeIchorHOMD ~{mergeIchorHOMD} --sex ~{sex} --outSegFile results/titan/combinedSegFile.txt \
        --outBinFile results/titan/combinedBinFile.txt --centromere ~{centromere} > ~{log} 2> ~{log}
    >>>
    runtime {
        docker: titan_docker
    }
    output {
        File segFile = "results/titan/combinedSegFile.txt"
        File binFile = "results/titan/combinedBinFile.txt"
    }
}

task selectSolution {
    input {
        Array[File] resultFiles
        String solutionRscript
        Float threshold
        String log
        String titan_docker
    }
    command <<<
        # Assuming the R script handles multiple result files
        Rscript ~{solutionRscript} --resultFiles ~{sep=' ' resultFiles} --threshold ~{threshold} --outFile results/titan/hmm/optimalClusterSolution.txt > ~{log} 2> ~{log}
    >>>
    runtime {
        docker: titan_docker
    }
    output {
        File optimalClusterSolution = "results/titan/hmm/optimalClusterSolution.txt"
    }
}

task copyOptSolution {
    input {
        File optimalClusterSolutionFile
        String outputDirectory
        String log
    }
    command <<<
        set -e
        curDir=$(pwd)
        while IFS= read -r line; do
            if [ "$line" != "path" ]; then
                echo "Copying $curDir/$line to ~{outputDirectory}"
                cp -r "$curDir/$line"* "~{outputDirectory}"
            fi
        done < <(cut -f11 ~{optimalClusterSolutionFile} | grep -v "path")
    >>>
    output {
        String copiedSolution = "~{outputDirectory}"
    }
}