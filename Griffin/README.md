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
Example files provided in `griffin_inputs.json`.\
If you wish not to complete the nucleosome profiling step, set `griffin.call_nucleosome_profiling` to false.
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
Based off of the [snakemake demo](https://github.com/adoebley/Griffin/wiki) of griffin.
