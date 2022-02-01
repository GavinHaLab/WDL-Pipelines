version 1.0
# Griffin WDL workflow 0.0 (an initial wiring)
## Author:  Amy Paguirigan (GitHub @vortexing)

## Sample-level inputs
## An array of the 'griffinInputs' struct will create a parallel analysis of many samples
## using the same set of parameters.
struct griffinInput {
  File bam_file
  File bam_index
  String sample_name
}

## The reference data bundle here can be defined in a struct for the entire batch.  
## Housing this input in a separate input json will allow for easy switching between 
## hg19 and hg38 as needed.  
struct referenceData {
    File chrom_sizes
    Array[File] excluded_regions
    File mappability_bw
    File ref_fasta
    File mappable_regions
    File genome_GC_freq_tar
    File sites_yaml
    File sites_tar
}

## Workflow definition
workflow runGriffin {
  input {
    Array[griffinInput] griffinBatch
    referenceData griffinReferences
  }
    ## Defined variables for this workflow, not likely to be changed per-run by the user, 
    ## and validated for use by the workflow writer

    String griffinDocker = "vortexing/griffin:v0.8"
    
    Array[String] chroms = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
                           "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", 
                           "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", 
                           "chr22"]

    ## Note: these params may be something users specify in their input json with the griffinBatch inputs
    ## OR, you can leave them here b/c these are the "best" way for users to run this particular workflow.
    Int map_quality = 20
    Pair[Int, Int] GC_bias_size_range = (15, 500)
    Pair[Int, Int] profiling_size_range = (100, 200)
    
    String chrom_column = "Chrom"
    String position_column = "position"
    String strand_column = "Strand"
    String number_of_sites = "none"
    String sort_by = "none"
    String ascending = "none"

    Pair[Int, Int] norm_window = (-5000, 5000)
    Pair[Int, Int] save_window = (-1000, 1000)
    Pair[Int, Int] center_window = (-30, 30)
    Pair[Int, Int] fft_window = (-960, 960)
    Int fft_index = 10
    Int smoothing_length = 165
    Int step_size = 15
    Boolean CNA_normalization = false
    Boolean individual = false
    Boolean smoothing = true
    Boolean exclude_zero_mappability = true
    Boolean exclude_outliers = true

## For each of the samples in the array of griffinBatch, perform the workflow tasks
scatter (sample in griffinBatch) {

    call compute_GC_counts {
    input:
      bam_file = sample.bam_file,
      bam_index = sample.bam_index,
      sample_name = sample.sample_name,
      mappable_regions = griffinReferences.mappable_regions,
      ref_fasta = griffinReferences.ref_fasta,
      chrom_sizes = griffinReferences.chrom_sizes,
      map_quality = map_quality,
      size_range = GC_bias_size_range,
      taskDocker = griffinDocker
  }
  
  call compute_GC_bias {
    input:
      sample_name = sample.sample_name,
      GC_counts = compute_GC_counts.out,
      mappable_name = "repeat_masker.mapable.k50.Umap.hg38", # b/c I can't get the reference workflow to run
      genome_GC_freq_tar  = griffinReferences.genome_GC_freq_tar,
      size_range = GC_bias_size_range,
      taskDocker = griffinDocker
  }

  call compute_mappability_bias {
    input:
      bam_file = sample.bam_file,
      bam_index = sample.bam_index,
      sample_name = sample.sample_name,
      mappability_bw = griffinReferences.mappability_bw,
      excluded_regions = griffinReferences.excluded_regions,
      chrom_sizes = griffinReferences.chrom_sizes,
      chroms = chroms,
      map_quality = map_quality,
      taskDocker = griffinDocker
  }

  call calculate_coverage {
    input:
      bam_file = sample.bam_file,
      bam_index = sample.bam_index,
      GC_bias = compute_GC_bias.bias,
      mappability_bias = compute_mappability_bias.bias,
      sample_name = sample.sample_name,
      ref_fasta = griffinReferences.ref_fasta,
      mappability_bw = griffinReferences.mappability_bw,
      chroms = chroms,
      chrom_sizes = griffinReferences.chrom_sizes,
      chrom_col = chrom_column,
      position_col = position_column,
      strand_col = strand_column,
      sites_yaml = griffinReferences.sites_yaml,
      sites_tar = griffinReferences.sites_tar,
      norm_window = norm_window,
      size_range = profiling_size_range,
      map_quality = map_quality,
      sort_by = sort_by,
      ascending = ascending,
      number_of_sites = number_of_sites,
      taskDocker = griffinDocker
  }

  call merge_sites {
    input:
      sample_name = sample.sample_name,
      mappability_bw = griffinReferences.mappability_bw,
      uncorrected_bw = calculate_coverage.uncorrected_bw,
      GC_corrected_bw = calculate_coverage.GC_corrected_bw,
      GC_map_corrected_bw = calculate_coverage.GC_map_corrected_bw,
      chroms = chroms,
      excluded_regions = griffinReferences.excluded_regions,
      chrom_sizes = griffinReferences.chrom_sizes,
      chrom_col = chrom_column,
      position_col = position_column,
      strand_col = strand_column,
      sites_yaml = griffinReferences.sites_yaml,
      sites_tar = griffinReferences.sites_tar,
      norm_window = norm_window,
      save_window = save_window,
      center_window = center_window,
      step_size = step_size,
      fft_window = fft_window,
      fft_index = fft_index,
      smoothing_length = smoothing_length,
      individual = individual,
      smoothing = smoothing,
      exclude_zero_mappability = exclude_zero_mappability,
      exclude_outliers = exclude_outliers,
      CNA_normalization = CNA_normalization,
      sort_by = sort_by,
      ascending = ascending,
      number_of_sites = number_of_sites,
      taskDocker = griffinDocker
  }

} ## End of Sample Scatter
  output {}
} ## End of workflow

## Individual task descriptions
task compute_mappability_bias {
  input {
    File bam_file
    File bam_index
    String sample_name
    File mappability_bw
    Array[File] excluded_regions
    File chrom_sizes
    Array[String] chroms
    Int map_quality
    String taskDocker
  }
    Int taskcpu = 1
  command {
    set -eo pipefail
    touch ~{bam_index}
    mkdir tmp
    python3 /Griffin/scripts/griffin_mappability_correction.py \
      --bam_file ~{bam_file} \
      --bam_file_name ~{sample_name} \
      --output ~{sample_name}.mappability_bias.txt \
      --output_plot ~{sample_name}.mappability_bias.pdf \
      --mappability ~{mappability_bw} \
      --exclude_paths ~{sep=' ' excluded_regions} \
      --chrom_sizes ~{chrom_sizes} \
      --chroms ~{sep=" " chroms} \
      --map_quality ~{map_quality} \
      --CPU ~{taskcpu} \
      --tmp_dir tmp/
  }
  runtime {
    cpu: taskcpu
    docker: taskDocker
  }
  output {
    File bias = "~{sample_name}.mappability_bias.txt"
    File plot = "~{sample_name}.mappability_bias.pdf"
  }
}



task compute_GC_counts {
  input {
    File bam_file
    File bam_index
    String sample_name
    File mappable_regions
    File ref_fasta
    File chrom_sizes
    Int map_quality
    Pair[Int, Int] size_range
    String taskDocker
  }
  Int taskcpu = 1
  command {
    set -eo pipefail
    touch ~{bam_index}
    python3 /Griffin/scripts/griffin_GC_counts.py \
      --bam_file ~{bam_file} \
      --bam_file_name ~{sample_name} \
      --mappable_regions ~{mappable_regions} \
      --ref_seq ~{ref_fasta} \
      --chrom_sizes ~{chrom_sizes} \
      --out_dir . \
      --map_q ~{map_quality} \
      --size_range ~{size_range.left} ~{size_range.right} \
      --CPU ~{taskcpu}
  }
  runtime {
    cpu: taskcpu
    docker: taskDocker
  }
  output {
    File out = "GC_counts/~{sample_name}.GC_counts.txt"
  }
}


task compute_GC_bias {
  input {
    String sample_name
    File GC_counts
    String mappable_name 
    File genome_GC_freq_tar 
    Pair[Int, Int] size_range
    String taskDocker
    }
    
  command {
    set -eo pipefail
    mkdir GC_counts
    cp ~{GC_counts} ./GC_counts
    mkdir genome_dir
    tar -xf ~{genome_GC_freq_tar} -C ./genome_dir
    python3 /Griffin/scripts/griffin_GC_bias.py \
      --bam_file_name ~{sample_name} \
      --mappable_name ~{mappable_name} \
      --genome_GC_frequency genome_dir/genome_GC_frequency \
      --out_dir . \
      --size_range ~{size_range.left} ~{size_range.right}
  }
  runtime {
    docker: taskDocker
  }
  output {
    File bias = "GC_bias/~{sample_name}.GC_bias.txt"
    File plot = "GC_plots/~{sample_name}.GC_bias.pdf"
    File summary = "GC_plots/~{sample_name}.GC_bias.summary.pdf"
    File key_lengths = "GC_plots/~{sample_name}.GC_bias.key_lengths.pdf"
    File countsout = "GC_counts/~{sample_name}.GC_counts.txt"
  }

}

## calculate_coverage is subject to the hacky approach of tar'ing a sites list file set and yaml'ing the listing of them

task calculate_coverage {
  input {
    File bam_file
    File bam_index
    File GC_bias
    File mappability_bias
    String sample_name
    File ref_fasta
    File mappability_bw
    Array[String] chroms
    File chrom_sizes
    String chrom_col
    String position_col
    String strand_col
    File sites_yaml
    File sites_tar
    Pair[Int, Int] norm_window
    Pair[Int, Int] size_range
    Int map_quality
    String sort_by
    String ascending
    String number_of_sites
    String taskDocker
  }
  Int taskcpu = 1
  command {
    set -eo pipefail
    touch ~{bam_index}
    mkdir sites
    tar -xzvf ~{sites_tar} -C sites

    python3 /Griffin/scripts/griffin_coverage.py \
      --sample_name ~{sample_name} \
      --bam ~{bam_file} \
      --GC_bias ~{GC_bias} \
      --mappability_bias ~{mappability_bias} \
      --tmp_dir . \
      --reference_genome ~{ref_fasta} \
      --mappability_bw ~{mappability_bw} \
      --chrom_sizes_path ~{chrom_sizes} \
      --griffin_scripts /Griffin/scripts \
      --sites_yaml ~{sites_yaml} \
      --chrom_column ~{chrom_col} \
      --position_column ~{position_col} \
      --strand_column ~{strand_col} \
      --chroms ~{sep=" " chroms} \
      --norm_window ~{norm_window.left} ~{norm_window.right} \
      --size_range ~{size_range.left} ~{size_range.right} \
      --map_quality ~{map_quality} \
      --number_of_sites ~{number_of_sites} \
      --sort_by ~{sort_by} \
      --ascending ~{ascending} \
      --CPU ~{taskcpu} 
  }
  runtime {
    cpu: taskcpu
    docker: taskDocker
  }
  output {
    File uncorrected_bw = "~{sample_name}/tmp_bigWig/~{sample_name}.uncorrected.bw"
    File GC_corrected_bw = "~{sample_name}/tmp_bigWig/~{sample_name}.GC_corrected.bw"
    File GC_map_corrected_bw = "~{sample_name}/tmp_bigWig/~{sample_name}.GC_map_corrected.bw"
  }
}

## merge_sites is subject to the hacky approach of tar'ing a sites list file set and yaml'ing the listing of them

task merge_sites {
  input {
    String sample_name
    File mappability_bw
    File uncorrected_bw
    File GC_corrected_bw
    File GC_map_corrected_bw
    Array[String] chroms
    Array[File] excluded_regions
    File chrom_sizes
    String chrom_col
    String position_col
    String strand_col
    File sites_yaml
    File sites_tar
    Pair[Int, Int] norm_window
    Pair[Int, Int] save_window
    Pair[Int, Int] center_window
    Int step_size
    Pair[Int, Int] fft_window
    Int fft_index
    Int smoothing_length
    Boolean individual
    Boolean smoothing
    Boolean exclude_zero_mappability
    Boolean exclude_outliers
    Boolean CNA_normalization
    String sort_by
    String ascending
    String number_of_sites
    String taskDocker

  }
  Int taskcpu = 1
  command {
    set -eo pipefail
    mkdir tmp_dir
    mkdir results
    mkdir sites
    tar -xzvf ~{sites_tar} -C sites
    python3 /Griffin/scripts/griffin_merge_sites.py \
      --sample_name ~{sample_name} \
      --uncorrected_bw_path ~{uncorrected_bw} \
      --GC_corrected_bw_path ~{GC_corrected_bw} \
      --GC_map_corrected_bw_path ~{GC_map_corrected_bw} \
      --tmp_dir tmp_dir \
      --results_dir results \
      --mappability_bw ~{mappability_bw} \
      --chrom_sizes_path ~{chrom_sizes} \
      --griffin_scripts /Griffin/scripts \
      --chrom_column ~{chrom_col} \
      --position_column ~{position_col} \
      --strand_column ~{strand_col} \
      --chroms ~{sep=" " chroms} \
      --sites_yaml ~{sites_yaml} \
      --norm_window ~{norm_window.left} ~{norm_window.right} \
      --save_window ~{save_window.left} ~{save_window.right} \
      --center_window ~{center_window.left} ~{center_window.right} \
      --fft_window ~{fft_window.left} ~{fft_window.right} \
      --fft_index ~{fft_index} \
      --smoothing_length ~{smoothing_length} \
      --exclude_paths ~{sep=' ' excluded_regions} \
      --step ~{step_size} \
      --CNA_normalization ~{CNA_normalization} \
      --individual ~{individual} \
      --smoothing ~{smoothing} \
      --exclude_outliers ~{exclude_outliers} \
      --exclude_zero_mappability ~{exclude_zero_mappability} \
      --number_of_sites ~{number_of_sites} \
      --sort_by ~{sort_by} \
      --ascending ~{ascending} \
      --CPU ~{taskcpu}
  }
  runtime {
    cpu: taskcpu
    docker: taskDocker
  }
  output {}
}

## This task is too dependent on file structures and snakemake yamls so needs adjustment of the code itself
task generate_plots {
  input {
    Pair[Int, Int] save_window
    String step_size
    String individual
    String taskDocker
  }
  command {
    set -eo pipefail
    mkdir out_dir
    python3 /Griffin/scripts/griffin_plot.py \
      --in_dir in_dir \ ## HMMMMMMM
      --samples_yaml samples_yaml \ ##### YIKES
      --save_window ~{save_window.left} ~{save_window.right} \
      --step ~{step_size} \
      --individual ~{individual} \
      --out_dir out_dir 
  }
  runtime {
    docker: taskDocker
  }
  output {}
}