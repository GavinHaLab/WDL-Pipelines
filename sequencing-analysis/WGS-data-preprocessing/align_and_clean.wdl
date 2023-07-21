version 1.0
# struct for the input files for a given sample, if downloadCache task is used
struct preprocessSample {
  String sample_name
  String dataset_id
  File unmappedCramBam
  File unmappedCraiBai
}

# struct for all the reference data needed for the run
struct wgsReferenceData {
  String ref_name
  File ref_fasta
  File ref_fasta_index
  File ref_dict
  # This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
  # listing the reference contigs that are "alternative". Leave blank in JSON for legacy
  # references such as b37 and hg19.
  File? ref_alt
  File ref_amb
  File ref_ann
  File ref_bwt
  File ref_pac
  File ref_sa
  File dbSNP_vcf
  File dbSNP_vcf_index
  Array[File] known_indels_sites_VCFs
  Array[File] known_indels_sites_indices
  File gnomad
  File gnomad_index
  File wgs_intervals
}

## WGS data preprocessing workflow for downstream variant calling.
## Input requirements:
## - Paired end sequencing data in unmapped BAM (uBAM) or unmapped CRAM (uCRAM) format that comply with the following requirements:
## - - files must pass validation by ValidateSamFile (a Picard tool)
##
## Output Files:
## - An analysis-ready recalibrated bam and it's index
## 
workflow WGS_preprocess_for_variants {
  input {
    Array[preprocessSample] batchInputs
    wgsReferenceData referenceDataSet
  }
    # Docker containers this workflow has been designed for
    String GATKdocker = "broadinstitute/gatk:4.2.2.0"
    String GATKBWAdocker = "fredhutch/gatk-bwa:4.2.20-0.7.17"

    Int scatter_count = 30
    Int bwaThreads = 24

  call SplitIntervals {
    input:
      intervals = referenceDataSet.wgs_intervals,
      ref_fasta = referenceDataSet.ref_fasta,
      ref_fasta_index = referenceDataSet.ref_fasta_index,
      ref_dict = referenceDataSet.ref_dict,
      scatter_count = scatter_count,
      taskDocker = GATKdocker
  }

  scatter (job in batchInputs) {
    # Incorporate a sample name, dataset id and the reference genome name into the filenames
    String base_file_name = job.sample_name + "_" + job.dataset_id + "." + referenceDataSet.ref_name

    # Due to large file sizes, using this task can substantially increase speed
    # call createDownloadCache {
    #   input:
    #     fileToCache = job.unmappedCram
    # }

  call MarkIlluminaAdapters {
    input:
      input_bam = job.unmappedCramBam,
      input_bai = job.unmappedCraiBai,
      base_file_name = base_file_name,
      taskDocker = GATKdocker
  }
    ## task took ~19 hours for a ~200GB bam file with 24 threads. 
    call bamtoBWAtoMergeAlignment {
    input:
      input_bam = MarkIlluminaAdapters.output_bam,
      base_file_name = base_file_name,
      ref_fasta = referenceDataSet.ref_fasta,
      ref_fasta_index = referenceDataSet.ref_fasta_index,
      ref_dict = referenceDataSet.ref_dict,
      ref_alt = referenceDataSet.ref_alt,
      ref_amb = referenceDataSet.ref_amb,
      ref_ann = referenceDataSet.ref_ann,
      ref_bwt = referenceDataSet.ref_bwt,
      ref_pac = referenceDataSet.ref_pac,
      ref_sa = referenceDataSet.ref_sa,
      bwaThreads = bwaThreads,
      taskDocker = GATKBWAdocker
  }

    # Aggregate aligned+merged flowcell BAM files and mark duplicates
  call MarkDuplicatesSpark {
    input:
      input_bam = bamtoBWAtoMergeAlignment.output_bam,
      output_bam_basename = base_file_name + ".aligned.duplicates_marked",
      metrics_filename = base_file_name + ".duplicate_metrics",
      taskDocker = GATKdocker
  }


  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter (subinterval in SplitIntervals.interval_files) {

  # Generate the recalibration model by interval and apply it
  call ApplyBaseRecalibrator {
    input:
      input_bam = MarkDuplicatesSpark.output_bam,
      input_bam_index = MarkDuplicatesSpark.output_bai,
      base_file_name = base_file_name,
      dbSNP_vcf = referenceDataSet.dbSNP_vcf,
      dbSNP_vcf_index = referenceDataSet.dbSNP_vcf_index,
      known_indels_sites_VCFs = referenceDataSet.known_indels_sites_VCFs,
      known_indels_sites_indices = referenceDataSet.known_indels_sites_indices,
      ref_dict = referenceDataSet.ref_dict,
      ref_fasta = referenceDataSet.ref_fasta,
      ref_fasta_index = referenceDataSet.ref_fasta_index,
      intervals = subinterval,
      taskDocker = GATKdocker
    }
  } # End of interval scatter
  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GatherBamFiles {
    input:
      input_bams = ApplyBaseRecalibrator.recalibrated_bam,
      base_file_name = base_file_name,
      taskDocker = GATKdocker
  }

  ##  Consider which type(s) of qc metrics you are interested in and assemble those tasks here. 
    # call CollectAlignmentSummaryMetrics {
    #   input: 
    #     input_bam = GatherBamFiles.output_bam,
    #     base_file_name = base_file_name,
    #     ref_fasta = referenceDataSet.ref_fasta,
    #     ref_fasta_index = referenceDataSet.ref_fasta_index,
    #     docker = GATKdocker
    # }
    call CollectWgsMetrics {
      input: 
        input_bam = GatherBamFiles.output_bam,
        input_bam_index = GatherBamFiles.output_bai,
        base_file_name = base_file_name,
        ref_fasta = referenceDataSet.ref_fasta,
        ref_fasta_index = referenceDataSet.ref_fasta_index,
        taskDocker = GATKdocker
    }
    
  } # End of job scatter

  # Outputs that will be retained when execution is complete
  output {
    Array[File] analysisReadyBam = GatherBamFiles.output_bam
    Array[File] analysisReadyIndex = GatherBamFiles.output_bai
    Array[File] wgsMetrics = CollectWgsMetrics.metrics
    #Array[File] alignmentMetrics = CollectAlignmentSummaryMetrics.out
  }
}# End workflow



#### TASK DEFINITIONS

# Generate Base Quality Score Recalibration (BQSR) model and apply it
task ApplyBaseRecalibrator {
  input {
    File input_bam
    File intervals 
    File input_bam_index
    String base_file_name
    File dbSNP_vcf
    File dbSNP_vcf_index
    Array[File] known_indels_sites_VCFs
    Array[File] known_indels_sites_indices
    File ref_dict
    File ref_fasta
    File ref_fasta_index
    String taskDocker
  }
  command {
  set -eo pipefail

  samtools index ~{input_bam}

  gatk --java-options "-Xms8g" \
      BaseRecalibrator \
      -R ~{ref_fasta} \
      -I ~{input_bam} \
      -O ~{base_file_name}.recal_data.csv \
      --known-sites ~{dbSNP_vcf} \
      --known-sites ~{sep=" --known-sites " known_indels_sites_VCFs} \
      --intervals ~{intervals} \
      --interval-padding 100 

  gatk --java-options "-Xms8g" \
      ApplyBQSR \
      -bqsr ~{base_file_name}.recal_data.csv \
      -I ~{input_bam} \
      -O ~{base_file_name}.recal.bam \
      -R ~{ref_fasta} \
      --intervals ~{intervals} \
      --interval-padding 100 

  #finds the current sort order of this bam file
  samtools view -H ~{base_file_name}.recal.bam | grep @SQ | sed 's/@SQ\tSN:\|LN://g' > ~{base_file_name}.sortOrder.txt

  }
  output {
    File recalibrated_bam = "~{base_file_name}.recal.bam"
    File recalibrated_bai = "~{base_file_name}.recal.bai"
    File sortOrder = "~{base_file_name}.sortOrder.txt"
  }
  runtime {
    memory: "36 GB"
    cpu: 2
    docker: taskDocker
  }
}

# Note these tasks will break if the read lengths in the bam are greater than 250.
task CollectWgsMetrics {
  input {
    File input_bam
    File input_bam_index
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    String taskDocker
  }

  command {
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms1g" \
      CollectWgsMetrics \
        --INPUT ~{input_bam} \
        --VALIDATION_STRINGENCY SILENT \
        --REFERENCE_SEQUENCE ~{ref_fasta} \
        --INCLUDE_BQ_HISTOGRAM true \
        --COUNT_UNPAIRED true \
        --OUTPUT ~{base_file_name}_wgs_metrics.txt \
        --USE_FAST_ALGORITHM true 
  }
  runtime {
    dockerSL: taskDocker
    cpu: 2
  }
  output {
    File metrics = "~{base_file_name}_wgs_metrics.txt"
  }
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
task GatherBamFiles {
  input {
    Array[File] input_bams
    String base_file_name
    String taskDocker
  }

  command {
    set -eo pipefail

    gatk --java-options "-Dsamjdk.compression_level=5 -Xms1g" \
      GatherBamFiles \
      --INPUT ~{sep=' --INPUT ' input_bams} \
      --OUTPUT ~{base_file_name}.bam \
      --CREATE_INDEX true \
      --VERBOSITY WARNING 
    
   mv ~{base_file_name}.bai ~{base_file_name}.bam.bai
  }
  runtime {
    dockerSL: taskDocker
    cpu: 2
  }
  output {
    File output_bam = "~{base_file_name}.bam"
    File output_bai = "~{base_file_name}.bam.bai"
  }
}

# this task is a hybrid task suggested by the Broad for computational efficiency as the pipes result in far less work than separating the tools
task bamtoBWAtoMergeAlignment {
  input {
    File input_bam
    String base_file_name
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa
    String taskDocker
    Int bwaThreads
  }
  command {
    set -eo pipefail

    mkdir tmp4pipe

    gatk --java-options "-Xmx16G" \
      SamToFastq \
        --INPUT ~{input_bam} \
        --FASTQ /dev/stdout \
        --CLIPPING_ATTRIBUTE XT \
        --CLIPPING_ACTION 2 \
        --INTERLEAVE true \
        --INCLUDE_NON_PF_READS true \
        --VERBOSITY WARNING | \
      bwa mem -M -t ~{bwaThreads-2} -p ~{ref_fasta} /dev/stdin | \
      gatk --java-options "-Xmx16G" \
        MergeBamAlignment \
          --ALIGNED_BAM /dev/stdin \
          --UNMAPPED_BAM ~{input_bam} \
          --OUTPUT ~{base_file_name}.aligned.merged.bam \
          --REFERENCE_SEQUENCE ~{ref_fasta} \
          --ADD_MATE_CIGAR true \
          --CLIP_ADAPTERS false \
          --CLIP_OVERLAPPING_READS true \
          --INCLUDE_SECONDARY_ALIGNMENTS true \
          --MAX_INSERTIONS_OR_DELETIONS -1 \
          --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
          --ATTRIBUTES_TO_RETAIN XS \
          --VERBOSITY WARNING \
          --SORT_ORDER queryname

  }
  output {
    File output_bam = "~{base_file_name}.aligned.merged.bam"
  }
  runtime {
    cpu: bwaThreads
    dockerSL: taskDocker
  }
}


task MarkIlluminaAdapters {
  input {
    File input_bam
    File input_bai
    String base_file_name
    String taskDocker
  }
  command {
    set -eo pipefail

    gatk --java-options "-Xmx8G" \
      MarkIlluminaAdapters \
        --INPUT ~{input_bam} \
        --OUTPUT ~{base_file_name}.markilluminaadapters.bam \
        --METRICS ~{base_file_name}.markilluminaadapters_metrics.txt 
  }
  output {
    File output_bam = "~{base_file_name}.markilluminaadapters.bam"
    File output_metrics = "~{base_file_name}.markilluminaadapters_metrics.txt"
  }
  runtime {
    cpu: 2
    dockerSL: taskDocker
  }
}




#### This is much faster than Cromwell localizing files for large, many-part files (one part is 5GB, so files need to be >30-40GB for this to become an issue).  
#### 
task createDownloadCache {
  input {
    String fileToCache
  }
    String outputFileName = basename(fileToCache)
  command {
    set -eo pipefail
    aws s3 cp --quiet ~{fileToCache} .
  }
  runtime {
    cpu: 12
    modules: "awscli/2.1.37-GCCcore-10.2.0"
  }
  output {
    File file = "~{outputFileName}"
  }
}



task MarkDuplicatesSpark {
  input {
    File input_bam
    String output_bam_basename
    String metrics_filename
    String taskDocker
  }
  # Later use: --VERBOSITY WARNING
  # Scales linearly up to 16 cores
  # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly.
  # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
  # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
  command {
    set -eo pipefail
    gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=4 -Dsamjdk.compression_level=5 -Xms32g" \
      MarkDuplicatesSpark \
      --input ~{input_bam} \
      --output ~{output_bam_basename}.bam \
      --metrics-file ~{metrics_filename} \
      --optical-duplicate-pixel-distance 2500 \
      --verbosity WARNING
  }
  runtime {
    dockerSL: taskDocker
    cpu: 16
    walltime: '12-0'
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bai = "~{output_bam_basename}.bam.bai"
    File duplicate_metrics = "~{metrics_filename}"
  }
}

task SplitIntervals {
    input {
      File intervals
      File ref_fasta
      File ref_fasta_index
      File ref_dict
      Int scatter_count
      String taskDocker
    }
    command {
        set -eo pipefail
        mkdir interval-files1
        gatk --java-options "-Xms2g" \
          SplitIntervals \
            -R ~{ref_fasta} \
            -L ~{intervals} \
            -scatter ~{scatter_count} \
            -O interval-files1 \
            --subdivision-mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
            --verbosity WARNING

        cp interval-files1/*.interval_list .
        ls *.interval_list > globGetAround.txt
    }
    runtime {
        dockerSL: taskDocker
        cpu: 2
    }
    output {
        Array[File] interval_files = read_lines("globGetAround.txt")
    }
}