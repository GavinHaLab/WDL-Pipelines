# ichorCNA-WDL
A revised version of [this WDL](https://github.com/GavinHaLab/wdl-repo/tree/main/processes/ichorCNA) so that it outputs plots.
Instructions on how to run a WDL on Cromwell [here](#instructions).

The ichorCNA.wdl workflow has two tasks:
- read_counter
- run_ichorCNA

and outputs 8 files:
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
Docker: By default fredhutch/ichorcna:v0.5.0 is set, and doesn't need to be changed unless you wish to.


## Outputs


## Instructions
You can run it if your workplace already has a Cromwell server configured or by other WDL execution tools.

This is just a rough explanation of what I've tested this with and how you can do it yourself, from a beginner's perspective.
You can run this on your own if you have [cromwell installed](https://github.com/broadinstitute/cromwell/releases/tag/85) and want to use the `run` or `server` mode.

### Run mode
To do this, enter this into your terminal:

    java -jar cromwell-XX.jar run ichorCNA.wdl -i inputs.json -o options.json

with `XX` being the version of cromwell you have. Make sure all of your files (WDL, input, options) in the same folder as the cromwell .jar file.

### Server mode
