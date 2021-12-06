version 1.0
# Note: As of 11/15/21, this still needs to be re-tested with current edits

# This workflow converts an aligned bam into an unmapped cram for storage.
# Creates unmapped crams with SAMPLE_ALIAS=~{sampleName} but incorporates
# a dataset_id via the TGR into the filename for future use. 

struct sampleInformation {
  String sampleName
  String dataset_id
  File alignedBam
}
workflow AlignedToUnMappedCram {
  input {
    Array[sampleInformation] samplesToUnmap
  }

  ## Docker containers this workflow has been validated with
  String GATKDocker = "broadinstitute/gatk:4.2.2.0"
  String samtoolsDocker = "vortexing/samtools:1.10"

scatter (job in samplesToUnmap) {

    String base_file_name = job.sampleName + "_" + job.dataset_id

    call RevertBam {
      input: 
        alignedBam = job.alignedBam,
        base_file_name = base_file_name,
        sampleName = job.sampleName, 
        taskDocker = GATKDocker
    }
    call CramABam {
      input:
        bamtocram = RevertBam.unmappedbam,
        base_file_name = base_file_name,
        taskDocker = samtoolsDocker,
        threads = 6
  }

    call ValidateCramorBam {
      input: 
        unmappedCramorBam = CramABam.cram,
        base_file_name = base_file_name,
        taskDocker = GATKDocker
  }
 } # End Scatter
  # Outputs that will be retained when execution is complete
  output {
    Array[File] unmappedCram = CramABam.cram
    Array[File] unmappedCramIndex = CramABam.cramIndex
    Array[File] validation = ValidateCramorBam.validation
    }
# End workflow
}

#### TASK DEFINITIONS
task CramABam {
  input {
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
    docker: taskDocker
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
        --INPUT=~{alignedBam} \
        --OUTPUT=~{base_file_name}.unmapped.bam \
        --SAMPLE_ALIAS=~{sampleName} 
  }
  runtime {
    cpu: 2
    memory: "32 GB"
    docker: taskDocker
  }
  output {
    File unmappedbam = "~{base_file_name}.unmapped.bam"
  }
}

# Validates cram/bam files for formatting issues. 
task ValidateCramorBam {
  input {
    File unmappedCramorBam
    String base_file_name
    String taskDocker
  }
  command {
    set -eo pipefail
    gatk --java-options "-Dsamjdk.compression_level=5 -Xms2g" \
      ValidateSamFile \
        --INPUT=~{unmappedCramorBam} \
        --MODE=SUMMARY \
        --IGNORE_WARNINGS=true > ~{base_file_name}.validation.txt
  }
  runtime {
    cpu: 2
    memory: "4 GB"
    docker: taskDocker
  }
  output {
    File validation = "~{base_file_name}.validation.txt"
  }
}


