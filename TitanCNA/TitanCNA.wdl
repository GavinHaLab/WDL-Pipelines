version 1.0

import getAlleleCounts.wdl as getAlleleCounts
import "https://raw.githubusercontent.com/GavinHaLab/WDL_Pipelines/59321d241165993d126f277466911cbd4ebe2be7/ichorCNA/ichorCNA.wdl?token=GHSAT0AAAAAACLT5MBAMQGXUWLN7MVCDFUOZMDOKGA" as ichorCNA

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
        # ## allele counts - samtools, pysam ##
        Int getAlleleCounts_mapQuality
        Int getAlleleCounts_baseQuality
        Int getAlleleCounts_vcfQuality

        # titan inputs
        Int maxClusters
        Int maxPloidy
        String libdir
        Int numCores
        Float normalInit
        String sex
        String genomeStyle
        String genomeBuild
        Boolean estimatePloidy
        Boolean estimateClonality
        Boolean estimateNormal
        Float alphaK
        Float txnExpLen
        Float plotYlim
        String combineScript
        Boolean mergeIchorHOMD
        String solutionRscript
        Float threshold
        String outputDirectory

        #path inputs
        String refFasta
        String snpVCF
        String cytobandFile
        String centromere
    }

    String titanDocker = "argage/titancna:vsomething.something.hopefully"
    
    String samTools = "/path/to/samtools"
    String bcfTools = "/path/to/bcftools"

    String getAlleleCounts_pyCountScript = "code/countPysam.py"
    String rscript = "../R_scripts/titanCNA.R"
    String combineTitanIchorCNA = "code/combineTITAN-ichor.R"
    String selectSolutionRscript = "../R_scripts/selectSolution.R"
    String libdir = "../../R/"

    

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
    call ichorCNA {
        input:
            batchSamples = Samples, #there might be issues as importing this as an array previously and ichor making it an array may make an Array[Array[]]
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
    call getAlleleCounts {
        input:
            tumors = batchSamples
            refFasta = refFasta
            snpDB = snpDB
            samtools = samTools
            bcftools = bcfTools
            countScript = getAlleleCounts_pyCountScript
            baseQ = getAlleleCounts_baseQuality
            mapQ = getAlleleCounts_mapQuality
            vcfQ = getAlleleCounts_vcfQuality
    }

    # Call runTitanCNA for each tumor
    Array[File] runTitanCNA_outputs = [] # defining outputs array before first call of scatter

    Array[File] concatenatedCounts = getAlleleCounts.concatenatedCounts
    Array[File] corrDepth = ichorCNA.corrDepth
    scatter (tumor in Samples) {
        Array[String] chrs = if tumor.genomeStyle == "NCBI" then ncbiChrs else ucscChrs

        # this is to get the specific file from concatenatedCounts from getAlleleCounts and corrDepth from ichor for the tumor in the scatter
        # taking all files that match the name (just one)
        # I'm hoping all the below code works in wdl?
        Array[File] alleleCountsFiles = select_all(concatenatedCounts, file -> file.name.contains(tumor.sampleName + ".tumCounts.txt")) 
        Array[File] corrDepthFiles = select_all(corrDepth, file -> file.name.contains(tumor.sampleName + ".correctedDepth.txt"))

        # turning array of one file to file
        File alleleCountsFile = alleleCountsFiles[0] 
        File corrDepthFile = corrDepthFiles[0]
        call runTitanCNA {
            input:
                alleleCounts = alleleCountsFile
                corrDepth = corrDepthFile,
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
                plotYlim = plotYlim
        }
        # append tumor's output- I could just move this to the runTitanCNA task if there's issues. Just would require more inputting and outputting
        Array[File] runTitanCNA_outputs = runTitanCNA_outputs + [runTitanCNA.titanOutput]

        Array[File] ichorSeg = ichorCNA.segTxt
        Array[File] ichorBin = ichorCNA.seg
        Array[File] ichorParam = ichorCNA.params

        #unpacking array to get tumor files from ichor
        Array[File] ichorSegFiles = select_all(ichorSeg, file -> file.name.contains(tumor.sampleName + ".seg.txt")) 
        Array[File] ichorBinFiles = select_all(ichorBin, file -> file.name.contains(tumor.sampleName + ".cna.seg")) 
        Array[File] ichorParamFiles = select_all(ichorParam, file -> file.name.contains(tumor.sampleName + ".params.txt")) 

        File ichorSegFile = ichorSegFiles[0]
        File ichorBinFile = ichorBinFiles[0]
        File ichorParamFile = ichorParamFiles[0]
        
        # Call combineTitanAndIchorCNA
        call combineTitanAndIchorCNA {
            input:
                titanSeg = runTitanCNA.segTxt,
                titanBin = runTitanCNA.titanOutput,
                titanParam = runTitanCNA.param,
                ichorSeg = ichorSegFile, 
                ichorBin = ichorBinFile, 
                ichorParam = ichorParamFile,  
                combineScript = combineTitanIchorCNA,
                libdir = libdir,
                mergeIchorHOMD = mergeIchorHOMD,
                sex = sex,
                centromere = centromere,
                log = "combineTitanAndIchorCNA.log"
        }
    }
    # Call selectSolution
    call selectSolution {
        input:
            resultFiles = runTitanCNA_outputs,
            solutionRscript = selectSolutionRscript,
            threshold = threshold,
            log = "selectSolution.log"
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
        File combinedSegFile = combineTitanAndIchorCNA.segFile
        File combinedBinFile = combineTitanAndIchorCNA.binFile
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
        String chrs
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
        Float plotYlim
    }
    command <<<
        Rscript ~{titanRscript} --hetFile ~{alleleCounts} --cnFile ~{corrDepth} \
        --outFile results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.titan.txt \
        --outSeg results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorName}_cluster~{clustNum}.segs.txt \
        --outParam results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumorname}_cluster~{clustNum}.params.txt \
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
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3
        }
    output {
        File titanOutput = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumor}_cluster~{clustNum}.titan.txt"
        File segTxt = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumor}_cluster~{clustNum}.segs.txt"
        File param = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumor}_cluster~{clustNum}.params.txt"
        File seg = "results/titan/hmm/titanCNA_ploidy~{ploidy}/~{tumor}_cluster~{clustNum}.seg"
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