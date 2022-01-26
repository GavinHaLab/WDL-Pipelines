version 1.0
## Griffin workflow reference file preparation for hg38
workflow griffinReferencePrep {
  input {
    File ref_fasta
    File mappable_regions
    File chrom_sizes
    Int read_length
  }
  Pair[Int, Int] fragment_size_range = (1, 500)
  #Array[Int] fragment_size = ["1"]
  String griffinDocker = "vortexing/griffin:v0.1"
  call createRange {
    input:
      fragment_size_range = fragment_size_range
  }

  scatter (frag in createRange.seq) {
    #scatter (frag in fragment_size) {
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
    String regions_name = basename(mappable_regions)
  call createTar {
    input:
      files = calc_GC_frequency.out,
      tarName = regions_name + ".GC_frequency.tar.gz"
  }
  output {
    File tar = createTar.tar
  }
}


task createRange {
  input {
    Pair[Int, Int] fragment_size_range
  }
  command {
    set -eo
    seq ~{fragment_size_range.left} ~{fragment_size_range.right} > sequence.txt
  }
  output {
    Array[Int] seq = read_lines("sequence.txt")
  }
}
task calc_GC_frequency {
  input {
    File mappable_regions
    File ref_fasta
    File chrom_sizes
    Int read_length
    Int fragment_length
    String taskDocker
  }
    String outfilename = basename(mappable_regions)
  command {
    set -eo pipefail

    python3 /Griffin/scripts/griffin_calc_GC_frequency.py \
      --mappable_regions_path ~{mappable_regions} \
      --ref_seq ~{ref_fasta} \
      --chrom_sizes ~{chrom_sizes} \
      --out_dir . \
      --read_length ~{read_length} \
      --fragment_length ~{fragment_length} 
  }
  runtime {
    dockerSL: taskDocker
  }
  output {
    File out = "~{outfilename}.~{fragment_length}bp.GC_frequency.txt"
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
  output {
    File tar = "~{tarName}"
  }
}