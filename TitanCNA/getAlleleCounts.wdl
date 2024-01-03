version 1.0

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

workflow getAlleleCounts {
    input {
        Array[sampleData] tumors
        String refFasta
        String snpDB
        String samtools
        String bcftools
        String countScript
        Int baseQ
        Int mapQ
        Int vcfQ
    }
    Array[String] ucscChrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                           "chr22", "chrX", "chrY"]

    Array[String] ncbiChrs = ["1", "2", "3", "4", "5", "6", "7", 
                          "8", "9", "10", "11", "12", "13", "14", 
                          "15", "16", "17", "18", "19", "20", "21", 
                          "22", "X", "Y"]

    ## Process each tumor and chromosome
    scatter (tumor in tumors) {
        Array[String] chrs = if tumor.genomeStyle == "NCBI" then ncbiChrs else ucscChrs
        scatter (chr in chrs) {
            call getHETsites {
                input:
                    sample = tumor.normalBam,
                    tumorName = tumor.sampleName,
                    chr = chr,
                    refFasta = refFasta,
                    snpDB = snpDB,
                    samtools = samtools,
                    bcftools = bcftools
            }

            call getAlleleCountsByChr {
                input:
                    hetSites = getHETsites.hetSites,
                    tumor = tumor.tumorBam,
                    tumorName = tumor.sampleName,
                    chr = chr,
                    countScript = countScript,
                    baseQ = baseQ,
                    mapQ = mapQ,
                    vcfQ = vcfQ
            }
        }

        ## Concatenate Allele Count Files
        call catAlleleCountFiles {
            input:
                alleleCountFiles = getAlleleCountsByChr.alleleCounts, #array here?
                tumorName = tumor.sampleName
        }
    }
    output {
        Array[Array[File]] hetSites = getHETsites.hetSites
        Array[Array[File]] alleleCounts = getAlleleCountsByChr.alleleCounts
        Array[File] concatenatedCounts = catAlleleCountFiles.concatenatedCounts
    }
}

task getHETsites {
    input {
        #lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
        # ^ appears to just use the normal sample paired with the tumor
        String tumorName
        File? sample #this needs to be the path to the normal bam
        String refFasta
        String snpDB
        String samtools
        String bcftools
        String chr
    }
    command {
        ~{samtools} mpileup -uv -I -f ~{refFasta} -r ~{chr} -l ~{snpDB} ~{sample} | ~{bcftools} call -v -c - | grep -e '0/1' -e '#' > hetSites.vcf
    }
    output {
        #"results/titan/hetPosns/{tumor}/{tumor}.chr{chr}.vcf"
        File hetSites = "~{tumorName}.chr~{chr}.vcf" #not sure if I need to add results/titan/hetPosns/{tumor}/ or if cromwell already puts into its own folder
    }
}

task getAlleleCountsByChr {
    input {
        #hetSites="results/titan/hetPosns/{tumor}/{tumor}.chr{chr}.vcf",
		#tumBam=lambda wildcards: config["samples"][wildcards.tumor]
    
        File hetSites
        File tumor
        String tumorName
        String chr
        String countScript
        Int baseQ
        Int mapQ
        Int vcfQ
    }
    command {
        python ~{countScript} ~{chr} ~{hetSites} ~{tumor} ~{baseQ} ~{mapQ} ~{vcfQ} > alleleCounts.txt
    }
    output {
        #"results/titan/tumCounts/{tumor}/{tumor}.tumCounts.chr{chr}.txt"
        File alleleCounts = "~{tumorName}.tumCounts.chr{chr}.txt"
    }
}

task catAlleleCountFiles {
    input {
        #expand("results/titan/tumCounts/{{tumor}}/{{tumor}}.tumCounts.chr{chr}.txt", chr=CHRS)
        Array[File] alleleCountFiles
        String tumorName
    }
    command {
        cat ~{sep=' ' alleleCountFiles} | grep -v Chr > concatenatedCounts.txt
    }
    output {
        #"results/titan/tumCounts/{tumor}.tumCounts.txt"
        File concatenatedCounts = "~{tumorName}.tumCounts.txt"
    }
}


