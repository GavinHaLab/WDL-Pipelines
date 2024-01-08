# Griffin wdl
Griffin source code [here](https://github.com/GavinHaLab/Griffin)- this is where you can find complete info on Griffin.

Griffin has two steps that this WDL divides into the following tasks:
1) GC correction
  -  `GC_counts`
  -  `GC_bias`
2) Nucleosome profiling
  - `calc_cov`
  - `merge_sites`
  - `generate_plots`

### Inputs
For each sample:

"`sample_name`"  
"`bam_file`"- path to sample bam file  

For all samples:

"`griffin.griffin_docker`"- Default: "vortexing/griffin:v0.2.0"  
"`griffin.GC_counts.disk_size`"- Default: 120  
"`griffin.GC_counts.cpu_num`"- Default: 8  
"`griffin.GC_counts.mem_size`": Default: 8  

"`griffin.mappability_bw`"- Path to your file that looks like: k100.Umap.MultiTrackMappability.bw  
"`griffin.reference_genome`":- Path to your file that looks like: Homo_sapiens_assembly38.fasta  
"`griffin.gaps`"- Default: "/Griffin/Ref/hg38_gaps.bed"  
"`griffin.genome_GC_frequency`"- Default: "/Griffin/Ref/genome_GC_frequency"  
"`griffin.mappable_regions`"- Default: "/Griffin/Ref/k100_minus_exclusion_lists.mappable_regions.gh38.bed"  
"`griffin.chrom_sizes`"- Default: "/Griffin/Ref/hg38.standard.chrom.sizes"  
"`griffin.patches`"- Default: "/Griffin/Ref/hg38_fix_patches.bed"  
"`griffin.centromeres`"- Default: "/Griffin/Ref/hg38_centromeres.bed"  
"`griffin.alternative_haplotypes`"- Default: "/Griffin/Ref/hg38_alternative_haplotypes.bed"  
"`griffin.encode_exclude`"- Default: "/Griffin/Ref/encode_unified_GRCh38_exclusion_list.bed"  
"`griffin.map_quality`"- minimum mapping quality to keep a read. Default: 20  
"`griffin.GC_bias_size_range`"- i.e. [15,500]  

"`griffin.sites_yaml`"- Path to your sites yaml file. Used to run calc coverage step on all sites  
"`griffin.chroms`"- Default: ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]  
"`griffin.size_range`"- Range of fragment lengths to be used for analysis. 100 to 200 captures nucleosome sized fragments. 15 to 500 also okay  
"`griffin.norm_window`"- Window around each site for normalizing to 1 (-5000 5000 bp for typical TFBS WGS analysis)  
"`griffin.save_window`"- window around each site to save to outputs  
"`griffin.center_window`"- range of positions used to calculate the central coverage feature  
"`griffin.fft_window`"- array of two integers  
"`griffin.fft_index`"-  integer  
"`griffin.number_of_sites`"- how many sites to analyze. use 0 to analyze all sites  
"`griffin.sort_by`"- column to use for sorting sites. use "none" if analyzing all sites  
"`griffin.ascending`"- true/false whether to sort sites in ascending order (False if you want to keep sites with a large value in the sort_by column), use 'none' to analyze all sites  
"`griffin.step`"- integer  
"`griffin.chrom_column`"- Column containing the chromosome in your sites file. Default: "Chrom"  
"`griffin.strand_column`"- column for indicating site direction. If this column doesn't exist, the script will assume non-directional sites. Default: "Strand"  
"`griffin.position_column`"- Column containing the site position. Default: "position"  
"`griffin.CNA_normalization`"- true/false  
"`griffin.exclude_outliers`"- true/false  
"`griffin.exclude_zero_mappability`"- true/false  
"`griffin.individual`"- true/false  
"`griffin.smoothing`"- true/false  
"`griffin.smoothing_length`"- approximately the fragment length  
"`griffin.call_nucleosome_profiling`"- true/false

### Outputs
GC_counts_file: `results/GC_counts/~{sample_name}.GC_counts.txt`

GC_bias_file: `GC_bias.GC_bias_file #results/GC_bias/~{sample_name}.GC_bias.txt`

GC_plots: `results/GC_plots/~{sample_name}.GC_bias.summary.pdf`

plot: `results/GC_plots/~{sample_name}.GC_bias.pdf`

key_lengths: `results/GC_plots/~{sample_name}.GC_bias.key_lengths.pdf`

sites_yaml: `griffin_nucleosome_profiling_files/sites/sites.yaml`

uncorrected_bw; `results/calc_cov/temp/~{sample_name}/tmp_bigWig/~{sample_name}.uncorrected.bw`

GC_corrected_bw: `results/calc_cov/temp/~{sample_name}/tmp_bigWig/~{sample_name}.GC_corrected.bw`

uncorrected_cov: `results/merge_sites/~{sample_name}/~{sample_name}.uncorrected.coverage.tsv`

GC_corrected_cov: `results/merge_sites/~{sample_name}/~{sample_name}.GC_corrected.coverage.tsv`

output_plots: `glob("results/plots/*.pdf")`

samples_yaml: `griffin_nucleosome_profiling_files/config/samples.GC.yaml`


## Demo

Based on the [snakemake demo](https://github.com/adoebley/Griffin/wiki) of griffin.

Download the reference genome from the link below (if you don't have wget, you can open the link in a browser) and unzip it (the download may take ~10 minutes):\
`wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz`\
`gunzip hg38.fa.gz`\

Download the mappability track (1.2gb) from the link below (if you don't have wget, you can open the link in a browser)(the download may take ~15 minutes):
`wget https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k100.Umap.MultiTrackMappability.bw`

Convert the demo cram file to bam and create an index (takes ~1 minute):
`samtools view -b -T Ref/hg38.fa -o demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam`\
`demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.cram`\

`samtools index demo/bam/Healthy_GSM1833219_downsampled.sorted.mini.bam`\

Download the griffin wdl and griffin_inputs.json. 

In the json file:

Edit `"griffin.samples"` to include the samples and the path to the bam files. Edit `"griffin.mappability_bw"` and `"griffin.reference_genome"` so that your cromwell server can access it.

Any directories starting with `/Griffin/Ref/` can be replaced if necessary, but should be able to run on its own as those files are provided in the docker. 

Add `"griffin.sites_yaml"`.

If you don't want to do nucleosome profiling, set `"griffin.call_nucleosome_profiling"` to false.

Edit any of the other inputs as needed, but should be fine left alone.

### Running via command line
Log into your chosen node and set the current directory to where your griffin.wdl, griffin_inputs.json, and your optional workflow_options.json files are stored. Load the modules necessary:

`module load java`

`module load cromwell`

To execute the WDL, do the following (assuming the WDL is named griffin.wdl):

`java -jar path/to/cromwell.jar run griffin.wdl -i griffin_inputs.json`

For FH server:

`java -jar $EBROOTCROMWELL/cromwell.jar run griffin.wdl -i griffin_inputs.json`

To [configure the workflow](https://github.com/GavinHaLab/WDL_Pipelines/tree/main/workflow-options), add `--options workflow_options.json` to the line.

### Running via shiny app - for FH users

The DaSL has a helpful [tutorial on running workflows via the Shiny App](https://hutchdatascience.org/FH_WDL101_Cromwell/using-shiny-to-manage-workflows.html).

Open the [shiny app](https://cromwellapp.fredhutch.org/) and connect to your cromwell server. Under "Submit a Workflow", upload the WDL and griffin_inputs.json. If you wanted to configure the workflow, [take a workflow options file](https://github.com/GavinHaLab/WDL_Pipelines/tree/main/workflow-options) and upload it as well.
