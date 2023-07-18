# ichorCNA-wdl
A revised version of [this WDL](https://github.com/GavinHaLab/wdl-repo/tree/main/processes/ichorCNA) so that it outputs plots.
[Here is info on a snakemake version](https://github.com/broadinstitute/ichorCNA/wiki/SnakeMake-pipeline-for-ichorCNA) for comparison. 

The ichorCNA.wdl workflow has two tasks:
- `read_counter`
- `run_ichorCNA`

The wdl calls a `scattter` and calls both tasks for each sample in `batchSamples`.

The workflow outputs 8 files:
- `~{sampleName}.bin~{binSize}.wig`
- `~{sampleId}.correctedDepth.txt`
- `~{sampleId}.cna.seg`
- `~{sampleId}.params.txt`
- `~{sampleId}.seg.txt`
- `~{sampleId}.seg`
- `~{sampleId}.RData`
- `~{sampleId}.plots.tar.gz`

## ichorCNA Wiki Page
For more info on the ichorCNA tool, visit the [Github Wiki page for ichorCNA](https://github.com/broadinstitute/ichorCNA/wiki) by the BROAD Institute.

## Configuration
### ichorCNA.wdl file
`taskDocker`: By default fredhutch/ichorcna:v0.5.0 is set, and doesn't need to be changed unless you wish to.\
`taskCPU`: CPU variable for run_ichorCNA. May have to edit `read_counter` CPU in its runtime section.\
`memory`: This may have to be edited in both `read_counter` and `run_ichorCNA` in their runtime sections.

### inputs.json file
`ichorCNA.batchSamples`: Each sample encased in `{}` contains a `sampleName`, `sex`, `tumorBam`, `tumorBai`, optional `normalBam`, optional `normalBai`, `normalPanel`, `genomeBuild`, and `genomeStyle`- which is NCBI or UCSC if hg38.

Everything else in the inputs.json file have example inputs that you can refer to. For the `ichorCNA.plotFileType` it can be either png or pdf.

## Outputs
A complete list of outputs can be found in [this Github wiki page](https://github.com/broadinstitute/ichorCNA/wiki/Output) along with parameter info.

## Instructions
You can run it if your workplace already has a Cromwell server configured or by other WDL execution tools.

This is just a rough explanation of what I've tested this with and how you can do it yourself, from a beginner's perspective.
You can run this on your own if you have [cromwell installed](https://github.com/broadinstitute/cromwell/releases/tag/85) and want to use the `run` or `server` mode.

### Run mode
To do this, enter this into your terminal:

    java -jar cromwell-XX.jar run ichorCNA.wdl -i inputs.json -o options.json

with `XX` being the version of cromwell you have. Make sure all of your files (WDL, input, options) are in the same folder as the cromwell .jar file.

### Server mode
[Here is a tutorial on how to run Cromwell's server mode](https://cromwell.readthedocs.io/en/stable/tutorials/ServerMode/). Skip the Five Minute Introduction if you've already downloaded Cromwell and familiar with it.
