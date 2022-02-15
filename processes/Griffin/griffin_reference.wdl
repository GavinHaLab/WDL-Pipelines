version 1.0

## Griffin workflow reference file bundle preparation (in this case for hg38)
workflow griffinReferencePrep {
  ## Inputs are here instead of a separate json as this is rarely run and easier to contain in one file since it's such a small "workflow".
  File ref_fasta = "/fh/scratch/delete90/paguirigan_a/apaguiri/cromwell-executions/ReferenceDataSets/Homo_sapiens_assembly38.fasta"
  File chrom_sizes = "/fh/scratch/delete90/paguirigan_a/apaguiri/cromwell-executions/ReferenceDataSets/hg38.standard.chrom.sizes"
  File mappable_regions = "/fh/scratch/delete90/paguirigan_a/apaguiri/cromwell-executions/ReferenceDataSets/k100_minus_exclusion_lists.mappable_regions.hg38.bed"
  Int read_length = 500

  ## Docker container validated for Griffin
  String griffinDocker = "vortexing/griffin:v0.9"
 
  scatter (frag in range(read_length)) {
    call calc_GC_frequency {
      input:
        mappable_regions = mappable_regions,
        ref_fasta = ref_fasta,
        chrom_sizes = chrom_sizes,
        fragment_length = frag,
        read_length = read_length,
        taskDocker = griffinDocker
    }
  }
  ## Could just concatenate them all into one file in the future, otherwise create a bundle for use as is.
    String regions_name = basename(mappable_regions)
  call createTar {
    input:
      files = calc_GC_frequency.out,
      tarName = regions_name + ".maxfrag_size." + read_length + ".GC_frequency.tar.gz"
  }
  output {
    File tar = createTar.tar
  }
}


## Task definitions
task calc_GC_frequency {
  input {
    File mappable_regions
    File ref_fasta
    File chrom_sizes
    Int read_length
    Int fragment_length
    String taskDocker
  }
    String outfilename = basename(mappable_regions, ".bed")
  command {
    set -eo pipefail

    python3 /Griffin/scripts/griffin_calc_GC_frequency.py \
      --mappable_regions_path ~{mappable_regions} \
      --ref_seq ~{ref_fasta} \
      --chrom_sizes ~{chrom_sizes} \
      --out_dir . \
      --read_length ~{read_length} \
      --fragment_length ~{fragment_length + 1} 
  }
  runtime {
    docker: taskDocker
  }
  output {
    File out = "~{outfilename}.~{fragment_length + 1}bp.GC_frequency.txt"
  }
}

task createTar {
  input {
    Array[File] files
    String tarName
  }
  command {
    set -eo
    tar -czf ~{tarName} ~{sep=" " files}
  }
  runtime {
    docker: "ubuntu:bionic"
  }
  output {
    File tar = "~{tarName}"
  }
}