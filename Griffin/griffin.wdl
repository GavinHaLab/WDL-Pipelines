version 1.0
struct sampleData {
    String sample_name
    File bam_file
}

workflow griffin {
    input {
        File mappability_bw
        File encode_exclude
        File centromeres
        File gaps
        File patches
        File alternative_haplotypes
        File chrom_sizes
        Int map_quality  # minimum mapping quality to keep a read
        String griffin_docker
        File mappable_regions
        File reference_genome
        File genome_GC_frequency
        Array[Int] GC_bias_size_range
        
        File sites_yaml
        String chrom_column # column containing the chromosome in your sites file
        String position_column # column containing the site position in your sites file
        String strand_column # column for indicating site direction in your sites file. If this column doesn't exist, the script will assume non-directional sites.
        Array[String] chroms
        Array[Int] norm_window # window around each site for normalizing to 1 (-5000 5000 bp for typical TFBS WGS analysis)
        Array[Int] size_range # range of fragment lengths to be used for analysis. 100 to 200 captures nucleosome sized fragments. 15 to 500 also okay.
        Int number_of_sites = 0 # how many sites to analyze. use 0 to analyze all sites
        String sort_by = 'none' # column to use for sorting sites use 'none' if analyzing all sites
        String ascending = 'none' # whether to sort sites in ascending order, use 'none' to analyze all sites
        Array[Int] save_window # window around each site to save to outputs
        Array[Int] center_window # range of positions used to calculate the central coverage feature
        Array[Int] fft_window
        Int fft_index
        Int smoothing_length # approximately the fragment length
        # bed files contraining regions to exclude from analysis
        Int step
        Boolean CNA_normalization = false
        Boolean individual = false #indivudal=true does not work yet
        Boolean smoothing = true
        Boolean exclude_zero_mappability = true
        Boolean exclude_outliers = true
        Boolean call_nucleosome_profiling
        Array[sampleData] samples
    }
    scatter (sample in samples) {
        call GC_counts {
            input:
                griffin_docker = griffin_docker,
                bam_file = sample.bam_file,
                sample_name = sample.sample_name,
                mappable_regions = mappable_regions,
                reference_genome = reference_genome,
                chrom_sizes = chrom_sizes,
                map_q = map_quality,
                GC_bias_size_range = GC_bias_size_range
        }

        call GC_bias {
            input:
                griffin_docker = griffin_docker,
                sample_name = sample.sample_name,
                GC_counts_file = GC_counts.GC_counts_file,
                mappable_name = basename(mappable_regions, ".bed"),
                genome_GC_frequency = genome_GC_frequency,
                GC_bias_size_range = GC_bias_size_range
        }

        if (call_nucleosome_profiling) {
            call calc_cov {
                input:
                    griffin_docker = griffin_docker,
                    bam_file = sample.bam_file,
                    GC_bias_file = GC_bias.GC_bias_file,
                    sample_name = sample.sample_name,
                    reference_genome = reference_genome,
                    chrom_sizes = chrom_sizes,
                    sites_yaml = sites_yaml,
                    chrom_column = chrom_column,
                    position_column = position_column,
                    strand_column = strand_column,
                    chroms = chroms,
                    norm_window = norm_window,
                    size_range = size_range,
                    map_quality = map_quality,
                    number_of_sites = number_of_sites,
                    sort_by = sort_by,
                    ascending = ascending
            }

            call merge_sites {
                input:
                    griffin_docker = griffin_docker,
                    sample_name = sample.sample_name,
                    uncorrected_bw = calc_cov.uncorrected_bw,
                    GC_corrected_bw = calc_cov.GC_corrected_bw,
                    mappability_bw = mappability_bw,
                    chrom_sizes = chrom_sizes,
                    sites_yaml = sites_yaml,
                    chrom_column = chrom_column,
                    position_column = position_column,
                    strand_column = strand_column,
                    chroms = chroms,
                    norm_window = norm_window,
                    save_window = save_window,
                    center_window = center_window,
                    fft_window = fft_window,
                    fft_index = fft_index,
                    smoothing_length = smoothing_length,
                    encode_exclude = encode_exclude,
                    centromeres = centromeres,
                    gaps = gaps,
                    patches = patches,
                    alternative_haplotypes = alternative_haplotypes,
                    step = step,
                    CNA_normalization = CNA_normalization,
                    individual = individual,
                    smoothing = smoothing,
                    exclude_zero_mappability = exclude_zero_mappability,
                    exclude_outliers = exclude_outliers,
                    number_of_sites = number_of_sites,
                    sort_by = sort_by,
                    ascending = ascending
            }

            call generate_plots {
                input:
                    griffin_docker = griffin_docker,
                    sample_name = sample.sample_name,
                    uncorrected_cov = merge_sites.uncorrected_cov,
                    GC_corrected_cov = merge_sites.GC_corrected_cov,
                    GC_bias_file = GC_bias.GC_bias_file,
                    bam_file = sample.bam_file,
                    save_window = save_window,
                    step = step,
                    individual = individual
            }
    }
    }

    output {
        Array[File] GC_counts_file = GC_counts.GC_counts_file #results/GC_counts/~{sample_name}.GC_counts.txt
        Array[File] GC_bias_file = GC_bias.GC_bias_file #results/GC_bias/~{sample_name}.GC_bias.txt
        Array[File] GC_plots = GC_bias.GC_plots #results/GC_plots/~{sample_name}.GC_bias.summary.pdf
        Array[File] plot = GC_bias.plot
        Array[File] key_lengths = GC_bias.key_lengths

        Array[File?] uncorrected_bw = calc_cov.uncorrected_bw #results/calc_cov/temp/~{sample_name}/tmp_bigWig/~{sample_name}.uncorrected.bw
        Array[File?] GC_corrected_bw = calc_cov.GC_corrected_bw #results/calc_cov/temp/~{sample_name}/tmp_bigWig/~{sample_name}.GC_corrected.bw
        Array[File?] uncorrected_cov = merge_sites.uncorrected_cov #results/merge_sites/~{sample_name}/~{sample_name}.uncorrected.coverage.tsv
        Array[File?] GC_corrected_cov = merge_sites.GC_corrected_cov #results/merge_sites/~{sample_name}/~{sample_name}.GC_corrected.coverage.tsv
        Array[File?] output_plots = generate_plots.output_plots #glob("results/plots/*.pdf")
        Array[File?] samples_yaml = generate_plots.samples_yaml #griffin_nucleosome_profiling_files/config/samples.GC.yaml
    }
    meta {
        
    }
}


task GC_counts {
    input {
        String griffin_docker
        File bam_file
        String sample_name
        File mappable_regions
        File reference_genome
        File chrom_sizes
        Int map_q
        Array[Int] GC_bias_size_range
        Int cpu_num = 8
        Int mem_size = 8
        Int disk_size = 100

    }

    command <<<
        set -e

        mkdir results
        samtools index ~{bam_file}

        python3 /Griffin/scripts/griffin_GC_counts.py \
        --bam_file ~{bam_file} \
        --bam_file_name ~{sample_name} \
        --mappable_regions_path ~{mappable_regions} \
        --ref_seq ~{reference_genome} \
        --chrom_sizes ~{chrom_sizes} \
        --out_dir results \
        --map_q ~{map_q} \
        --size_range ~{sep=" " GC_bias_size_range} \
        --CPU ~{cpu_num} \

    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3

}
    output {
        File GC_counts_file = "results/GC_counts/~{sample_name}.GC_counts.txt"
    }
}

task GC_bias {
    input {
        String griffin_docker
        String sample_name
        String mappable_name
        File GC_counts_file
        File genome_GC_frequency
        Array[Int] GC_bias_size_range
    }
    command <<<
        set -e

        # Uncompress the genome_GC_frequency.tar.gz into a direcotry named genome_GC_frequency
        tar -xvf ~{genome_GC_frequency} --directory ./

        # Copy the GC_count_file that were generated from GC_Count task
        mkdir -p results/GC_counts/
        cp ~{GC_counts_file} results/GC_counts/

        # Run griffin_GC_bias
        #conda run --no-capture-output \
        #-n griffin_env \
        python3 /Griffin/scripts/griffin_GC_bias.py \
        --bam_file_name ~{sample_name} \
        --mappable_name ~{mappable_name} \
        --genome_GC_frequency ./genome_GC_frequency \
        --out_dir results \
        --size_range ~{sep=" " GC_bias_size_range}

    >>>
    runtime {
        docker: griffin_docker
        preemptible: 2
        bootDiskSizeGb: 12
        cpu: 2
        memory: "4 GB"
        disks: "local-disk " + 25 + " HDD"
        preemptible: 2
        maxRetries: 3
}
    output {
        File GC_bias_file = "results/GC_bias/~{sample_name}.GC_bias.txt"
        File GC_plots = "results/GC_plots/~{sample_name}.GC_bias.summary.pdf"
        File plot = "results/GC_plots/~{sample_name}.GC_bias.pdf"
        File key_lengths = "results/GC_plots/~{sample_name}.GC_bias.key_lengths.pdf"
    }
}

task calc_cov {
    input {
        String griffin_docker
        File bam_file
        File GC_bias_file
        String sample_name
        File reference_genome
        File chrom_sizes
        File sites_yaml
        String chrom_column
        String position_column
        String strand_column
        Array[String] chroms
        Array[Int] norm_window
        Array[Int] size_range
        Int map_quality
        Int number_of_sites
        String sort_by
        String ascending
        # If cpu, mem, and disk size were not specified, use 8 cores, 10GB, and 100 GB as default
        Int cpu_num = 8
        Int mem_size = 32
        Int disk_size = 100
        }
    command <<<
        set -e
        # Index the input bam file and ref genome
        samtools index ~{bam_file}
        samtools faidx ~{reference_genome}

        # Make temporary directory
        mkdir -p results/calc_cov/temp/

        # Run griffin_coverage_script to calculate coverage
        
        python3 /Griffin/scripts/griffin_coverage.py \
        --sample_name ~{sample_name} \
        --bam ~{bam_file} \
        --GC_bias ~{GC_bias_file} \
        --mappability_bias none \
        --mappability_correction False \
        --tmp_dir results/calc_cov/temp \
        --reference_genome ~{reference_genome} \
        --mappability_bw none \
        --chrom_sizes_path ~{chrom_sizes} \
        --sites_yaml ~{sites_yaml} \
        --griffin_scripts /Griffin/scripts/ \
        --chrom_column ~{chrom_column} \
        --position_column ~{position_column} \
        --strand_column ~{strand_column} \
        --chroms ~{sep=" " chroms} \
        --norm_window ~{sep=" " norm_window} \
        --size_range ~{sep=" " size_range} \
        --map_quality ~{map_quality} \
        --number_of_sites ~{number_of_sites} \
        --sort_by ~{sort_by} \
        --ascending ~{ascending} \
        --CPU ~{cpu_num}
    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3
        }
    output {
        File uncorrected_bw = "results/calc_cov/temp/~{sample_name}/tmp_bigWig/~{sample_name}.uncorrected.bw"
        File GC_corrected_bw = "results/calc_cov/temp/~{sample_name}/tmp_bigWig/~{sample_name}.GC_corrected.bw"
        }
}

task merge_sites {
    input {
        String griffin_docker
        String sample_name
        File uncorrected_bw
        File GC_corrected_bw
        File mappability_bw
        File chrom_sizes
        File sites_yaml
        String chrom_column
        String position_column
        String strand_column
        Array[String] chroms
        Array[Int] norm_window
        Array[Int] save_window
        Array[Int] center_window
        Array[Int] fft_window
        Int fft_index
        Int smoothing_length
        File encode_exclude
        File centromeres
        File gaps
        File patches
        File alternative_haplotypes
        Int step
        Boolean CNA_normalization
        Boolean individual
        Boolean smoothing
        Boolean exclude_zero_mappability
        Boolean exclude_outliers
        Int number_of_sites
        String sort_by
        String ascending
        # If cpu, mem, and disk size were not specified, use 8 cores, 10GB, and 100 GB as default
        Int cpu_num = 8
        Int mem_size = 32
        Int disk_size = 100
        }
    command <<<
        set -e
        # Make temporary directory
        mkdir -p results/merge_sites/temp/

        # Run griffin_merge_sites_script when mappability_correction is False
        # whenn mappability_correction is False, GC_map_corrected_bw_path is set to none
        
        python3 /Griffin/scripts/griffin_merge_sites.py \
        --sample_name ~{sample_name} \
        --uncorrected_bw_path ~{uncorrected_bw} \
        --GC_corrected_bw_path ~{GC_corrected_bw} \
        --GC_map_corrected_bw_path none \
        --mappability_correction False \
        --tmp_dir results/merge_sites/temp \
        --results_dir results/merge_sites \
        --mappability_bw ~{mappability_bw} \
        --chrom_sizes_path ~{chrom_sizes} \
        --sites_yaml ~{sites_yaml} \
        --griffin_scripts /Griffin/scripts/ \
        --chrom_column ~{chrom_column} \
        --position_column ~{position_column} \
        --strand_column ~{strand_column} \
        --chroms ~{sep=" " chroms} \
        --norm_window ~{sep=" " norm_window} \
        --save_window ~{sep=" " save_window} \
        --center_window ~{sep=" " center_window} \
        --fft_window ~{sep=" " fft_window} \
        --fft_index ~{fft_index} \
        --smoothing_length ~{smoothing_length} \
        --exclude_paths ~{encode_exclude} ~{centromeres} ~{gaps} ~{patches} ~{alternative_haplotypes} \
        --step ~{step} \
        --CNA_normalization ~{CNA_normalization} \
        --individual ~{individual} \
        --smoothing ~{smoothing} \
        --exclude_outliers ~{exclude_outliers} \
        --exclude_zero_mappability ~{exclude_zero_mappability} \
        --number_of_sites ~{number_of_sites} \
        --sort_by ~{sort_by} \
        --ascending ~{ascending} \
        --CPU ~{cpu_num}

    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3
        }
    output {
        File uncorrected_cov = "results/merge_sites/~{sample_name}/~{sample_name}.uncorrected.coverage.tsv"
        File GC_corrected_cov = "results/merge_sites/~{sample_name}/~{sample_name}.GC_corrected.coverage.tsv"
        }
}

#This task generates plots from the merged sites
task generate_plots {
    input {
        String griffin_docker
        String sample_name
        File uncorrected_cov
        File GC_corrected_cov
        File GC_bias_file
        File bam_file
        Array[Int] save_window
        Int step
        Boolean individual
        # If cpu, mem, and disk size were not specified, use 8 cores, 10GB, and 100 GB as default
        Int cpu_num = 1
        Int mem_size = 4
        Int disk_size = 100
        }
    command <<<
        set -e
        # Make temporary directory and copy files to results directory
        mkdir -p results/~{sample_name}/
        cp ~{uncorrected_cov} results/~{sample_name}/
        cp ~{GC_corrected_cov} results/~{sample_name}/

        # Make config directory to store samples yaml file
        mkdir -p griffin_nucleosome_profiling_files/config
        mkdir -p results/plots

        # Create a sample yaml file from input sample_name and bam_file
        echo "samples:
                  ~{sample_name}:
                    bam: ~{bam_file}
                    GC_bias: ~{GC_bias_file}" > griffin_nucleosome_profiling_files/config/samples.GC.yaml

        # Run griffin_generate_plots_script
        
        python3 /Griffin/scripts/griffin_plot.py \
        --in_dir results \
        --samples_yaml griffin_nucleosome_profiling_files/config/samples.GC.yaml \
        --mappability_correction False \
        --save_window ~{sep=" " save_window} \
        --step ~{step} \
        --individual ~{individual} \
        --out_dir results

        tar -czvf ~{sample_name}.plots.tar.gz results/plots/*.pdf

    >>>
    runtime {
        docker: griffin_docker
        bootDiskSizeGb: 12
        cpu: cpu_num
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 2
        maxRetries: 3
        }
    output {
        File output_plots = "results/~{sample_name}.plots.tar.gz"
        File samples_yaml = "griffin_nucleosome_profiling_files/config/samples.GC.yaml"
        }
}
