version 1.0
# This workflow converts an aligned bam into an unmapped cram for storage.
# Creates unmapped crams with SAMPLE_ALIAS=~{sampleName} but incorporates
# a dataset_id via the TGR into the filename for future use. 

struct sampleInformation {
  String sample_name
  String dataset_id
  String alignedBam
}

workflow AlignedToUnMappedCram {
  input {
    Array[sampleInformation] samplesToUnmap
  }

  ## Docker containers this workflow has been validated with
  String GATKDocker = "broadinstitute/gatk:4.2.2.0"
  String samtoolsDocker = "vortexing/samtools:1.10"

scatter (job in samplesToUnmap) {

    String base_file_name = job.sample_name + "_" + job.dataset_id
    call createDownloadCache {
      input:
        fileToCache = job.alignedBam
    }
    call RevertBam {
      input: 
        alignedBam = createDownloadCache.file,
        base_file_name = base_file_name,
        sampleName = job.sample_name, 
        taskDocker = GATKDocker
    }
    call AddOrReplaceReadGroups {
      input: 
        bam = RevertBam.unmappedbam,
        libraryName = base_file_name,
        platformName = "ILLUMINA",
        platformUnit = job.dataset_id,
        sampleName = job.sample_name,
        taskDocker = GATKDocker
    }
    call ValidateBam {
      input: 
        unmappedBam = AddOrReplaceReadGroups.repairedBam,
        base_file_name = base_file_name,
        taskDocker = GATKDocker
    }
    call CramABam {
      input:
        validationResult = ValidateBam.validation,
        bamtocram = AddOrReplaceReadGroups.repairedBam,
        base_file_name = base_file_name,
        taskDocker = samtoolsDocker,
        threads = 6
    }


 } # End Scatter
  # Outputs that will be retained when execution is complete
  output {
    Array[File] unmappedCram = CramABam.cram
    Array[File] unmappedCramIndex = CramABam.cramIndex
    Array[File] validation = ValidateBam.validation
    }
# End workflow
}

#### TASK DEFINITIONS


#### This is much faster for large, many-part files (one part is 5GB, so files need to be >30-40GB for this to become an issue).  
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
task CramABam {
  input { 
    File validationResult
    File bamtocram
    String base_file_name
    String taskDocker
    Int threads
  }
  command {
    set -eo pipefail

    samtools view -@~{threads-1} ~{bamtocram} -o ~{base_file_name}.cram
    samtools index -@~{threads-1} ~{base_file_name}.cram
    }
  runtime {
    dockerSL: taskDocker
    cpu: threads
  }
  output {
    File cram = "~{base_file_name}.cram"
    File cramIndex = "~{base_file_name}.cram.crai"
  }
}
# Reverts an aligned bam for a sample into an unmapped bam and validates it
task RevertBam {
  input {
    String base_file_name
    String sampleName
    File alignedBam
    String taskDocker
  }
  command {
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms4g" \
      RevertSam \
        --INPUT ~{alignedBam} \
        --OUTPUT ~{base_file_name}.unmapped.bam \
        --SAMPLE_ALIAS ~{sampleName} 
  }
  runtime {
    cpu: 2
    dockerSL: taskDocker
  }
  output {
    File unmappedbam = "~{base_file_name}.unmapped.bam"
  }
}


# Updatess read groups and tags if needed
task AddOrReplaceReadGroups {
  input {
    File bam
    String libraryName
    String platformName
    String platformUnit
    String sampleName
    String taskDocker
  }
    String filename = basename(bam, ".bam")
  command {
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms1g" \
      AddOrReplaceReadGroups \
        --INPUT ~{bam} \
        --OUTPUT ~{filename}.rgrepaired.bam \
        --RGLB ~{libraryName} \
        --RGPL ~{platformName} \
        --RGPU ~{platformUnit} \
        --RGSM ~{sampleName}
  }
  runtime {
    cpu: 2
    dockerSL: taskDocker
  }
  output {
    File repairedBam = "~{filename}.rgrepaired.bam"
  }
}


# Validates cram/bam files for formatting issues. 
task ValidateBam {
  input {
    File unmappedBam
    String base_file_name
    String taskDocker
  }
  command {
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2g" \
      ValidateSamFile \
        --INPUT ~{unmappedBam} \
        --MODE SUMMARY \
        --IGNORE_WARNINGS true > ~{base_file_name}.validation.txt
  }
  runtime {
    cpu: 2
    dockerSL: taskDocker
  }
  output {
    File validation = "~{base_file_name}.validation.txt"
  }
}


