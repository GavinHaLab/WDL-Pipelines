version 1.0

## workflow description
workflow ichorCNA {
  input {
    ## Future ref struct
    String genomeBuild = "hg19" # or hg19, only, capitalization matters
    String genomeStyle = "NCBI" # or NCBI, only
    File? exons  
    Pair[Int, String] binSize = (10000, "10kb") # These must match, to be improved later

    ## future sample struct
    # File? normalBam 
    # File? normalBai
    File tumorBam = "s3://fh-pi-ha-g-eco/Collaborator_Data/Cascadia/BCCRC_pilot/cfDNA_ULPWGS/RP-2598/WGS/TNBC032_P/v1/TNBC032_P.bam"
    File tumorBai = "s3://fh-pi-ha-g-eco/Collaborator_Data/Cascadia/BCCRC_pilot/cfDNA_ULPWGS/RP-2598/WGS/TNBC032_P/v1/TNBC032_P.bam.bai"
    
    String? normalPanel = "/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds"
    
    String sampleName = "TNBC032_P"
    String sex = "None"

    ## future params struct?
    Int qual = 20
    String normal = "c(0.5)"
    String ploidy = "c(2,3,4)"
    String estimateNormal = "TRUE" # these need to be quoted strings in all caps instead of Boolean b/c of R
    String estimatePloidy = "TRUE"
    String estimateClonality = "TRUE"
    String scStates = "c(1,3)" # states to use for subclonal CN
    Int maxCN = 8
    Float scPenalty = 1 # penalize subclonal events - n-fold multiplier; n=1 for no penalty,  
    String includeHOMD = "TRUE"
    String plotFileType = "pdf" # "png"
    String plotYlim = "c(-2,4)"
    String likModel = "t" # if multisample, use "gauss"
    Float minMapScore = 0.75 # control segmentation - higher (e.g. 0.9999999) leads to higher specificity and fewer segments
    Float maxFracGenomeSubclone = 0.5
    Float maxFracCNASubclone = 0.6
    Float normal2IgnoreSC = 0.90 # Ignore subclonal analysis when initial normal setting >= this value
    Float txnE = 0.9999 # lower (e.g. 0.99) leads to higher sensitivity and more segments
    Float txnStrength = 10000 # control segmentation - higher (e.g. 10000000) leads to higher specificity and fewer segments
    # lower (e.g. 100) leads to higher sensitivity and more segments
    Float fracReadsChrYMale = 0.001

  }
    ## Workflow and docker level params
    String ichorDocker = "vortexing/ichorcna:v0.4.9"
    
    Array[String] ucscChrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                           "chr22", "chrX"]

    Array[String] ncbiChrs = ["1", "2", "3", "4", "5", "6", "7", 
                          "8", "9", "10", "11", "12", "13", "14", 
                          "15", "16", "17", "18", "19", "20", "21", 
                          "22", "X"]

    Array[String] chrs = if genomeStyle == "NCBI" then ncbiChrs else ucscChrs

    String? centromere = if genomeBuild == "hg38" then "/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt" else "/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt"

  call read_counter as read_counter_tumor {
    input:
      bamFile = tumorBam,
      baiFile = tumorBai,
      sampleName = sampleName + "_tumor",
      binSize = binSize.left, 
      qual = qual,
      chrs = chrs,
      taskDocker = ichorDocker
  }

  ## Amy will revisit this, making it conditional based on whether it's run in paired mode or single mode
  # call read_counter as read_counter_normal {
  #   input:
  #     bamFile = normalBam,
  #     baiFile = normalBai,
  #     sampleName = sampleName + "_normal",
  #     binSize = binSize.left, 
  #     qual = qual,
  #     chrs = chrs,
  #     taskDocker = ichorDocker
  # }
    ## then call ichorCNA using that information
  call run_ichorCNA {
    input:
      tumorWig = read_counter_tumor.readDepth,
      #normalWig = read_counter_normal.readDepth,
      normalPanel = normalPanel,
      sampleId = sampleName,
      binSizeName = binSize.right,
      ichorChrs = chrs,
      sex = sex,
      normal = normal,
      ploidy = ploidy,
      genomeStyle = genomeStyle,
      genomeBuild = genomeBuild,
      estimateNormal = estimateNormal,
      estimatePloidy = estimatePloidy,
      estimateClonality = estimateClonality,
      scStates = scStates,
      maxCN = maxCN,
      centromere = centromere,
      exons = exons,
      includeHOMD = includeHOMD,
      txnE = txnE,
      txnStrength = txnStrength,
      plotFileType = plotFileType,
      plotYlim = plotYlim,
      likModel = likModel,
      minMapScore = minMapScore,
      maxFracGenomeSubclone = maxFracGenomeSubclone,
      maxFracCNASubclone = maxFracCNASubclone,
      normal2IgnoreSC = normal2IgnoreSC,
      scPenalty = scPenalty,
      fracReadsChrYMale = fracReadsChrYMale,
      taskDocker = ichorDocker
  }

  output {
    File tumorWig = read_counter_tumor.readDepth
    File corrDepth = run_ichorCNA.corrDepth
    File cna = run_ichorCNA.cna
    File segTxt = run_ichorCNA.segTxt
    File seg = run_ichorCNA.seg
    File rdata = run_ichorCNA.rdata
  }
}


## Task definitions
## Files input should be in form file.bam and file.bam.bai and in the same directory in their original locations so they are localized together by Cromwell too. 
task read_counter {
  input {
    File bamFile
    File baiFile
    String sampleName
    Int binSize
    Int qual
    Array[String] chrs
    String taskDocker
  }
  command {
    set -eo pipefail

    readCounter ~{bamFile} \
      -c ~{sep="," chrs} \
      -w ~{binSize} \
      -q ~{qual} > ~{sampleName}.bin~{binSize}.wig
  }
  runtime {
    memory: "4G"
    cpu: 2
    docker: taskDocker
  }
  output {
    File readDepth = "~{sampleName}.bin~{binSize}.wig"
  }
}

task run_ichorCNA {
  input {
    File tumorWig
    File? normalWig
    String? normalPanel
    String sampleId
    String binSizeName
    Array[String] ichorChrs
    String sex
    String normal
    String ploidy
    String genomeStyle
    String genomeBuild
    String estimateNormal
    String estimatePloidy
    String estimateClonality
    String scStates
    Int maxCN
    String? centromere
    File? exons
    String includeHOMD
    Float txnE
    Float txnStrength
    String plotFileType
    String plotYlim
    String likModel
    Float minMapScore
    Float maxFracGenomeSubclone
    Float maxFracCNASubclone
    Float normal2IgnoreSC
    Float scPenalty
    Float fracReadsChrYMale
    String taskDocker
  }
    Int taskCPU = 5
    ## Notes:  So, anything that R wants to be a number needs to be unquoted.  We cannot use the Boolean type here to my knowledge
    ## for the logical params b/c R gets confused when they aren't all caps TRUE/FALSE that are also unquoted.  
  command {
    set -o pipefail
    read -r -d '' VAR <<'EOF'
    library(ichorCNA); run_ichorCNA(id = '~{sampleId}',
    tumor_wig = '~{tumorWig}',
    repTimeWig = '/ichorCNA/inst/extdata/RepTiming_~{genomeBuild}_~{binSizeName}.wig',
    sex = '~{sex}',
    gcWig = '/ichorCNA/inst/extdata/gc_~{genomeBuild}_~{binSizeName}.wig',
    mapWig = '/ichorCNA/inst/extdata/map_~{genomeBuild}_~{binSizeName}.wig',
    ploidy = '~{ploidy}',
    normal = '~{normal}',
    maxCN = ~{maxCN},
    minMapScore = ~{minMapScore},
    chrs = 'c("~{sep='\", \"' ichorChrs}")',
    includeHOMD = ~{includeHOMD},
    genomeStyle = '~{genomeStyle}',
    genomeBuild = '~{genomeBuild}',
    estimateNormal = ~{estimateNormal},
    estimatePloidy = ~{estimatePloidy},
    estimateScPrevalence = ~{estimateClonality},
    scStates = '~{scStates}',
    likModel = '~{likModel}',
    maxFracGenomeSubclone = ~{maxFracGenomeSubclone},
    maxFracCNASubclone = ~{maxFracCNASubclone},
    normal2IgnoreSC = ~{normal2IgnoreSC},
    scPenalty = ~{scPenalty},
    txnE = ~{txnE},
    txnStrength = ~{txnStrength},
    fracReadsInChrYForMale = ~{fracReadsChrYMale},
    plotFileType = '~{plotFileType}',
    plotYLim = '~{plotYlim}',
    outDir = '.',
    cores  = ~{taskCPU} ~{", centromere = '" + centromere + "'"} ~{", exons.bed = '" + exons + "'"} ~{", normal_wig = '" + normalWig + "'"} ~{", normal_panel = '" + normalPanel + "'"})
    EOF

    echo "$VAR" | Rscript -
  }
  runtime {
    memory: "10G"
    cpu: taskCPU
    docker: taskDocker
  }
  output {
    File corrDepth = "~{sampleId}.correctedDepth.txt"
    File cna = "~{sampleId}.cna.seg"
    File params = "~{sampleId}.params.txt"
    File segTxt = "~{sampleId}.seg.txt"
    File seg = "~{sampleId}.seg"
    File rdata = "~{sampleId}.RData"
  }
}
