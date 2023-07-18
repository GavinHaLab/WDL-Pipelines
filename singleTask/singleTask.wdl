version 1.0
# testing 123~!

workflow oneTask {
  input {
    Array[File] inputFileArray = ["fh/scratch/delete90/ha_g/file1.tsv","fh/scratch/delete90/ha_g/file2.tsv"]
  }

  ## Docker containers this workflow has been validated with
  String moduleToUse = "Python/3.9.6-GCCore-11.2.0"

scatter (thisone in inputFileArray) {

    call summarize {
      input: 
        fileInput = thisone,
        taskModule = moduleToUse
    }

 } # End Scatter
  # Outputs that will be retained when execution is complete
  output {
    Array[File] summarizeFiles = summarize.outputFile
    }
# End workflow
}

#### TASK DEFINITIONS

task summarize {
  input {
    File fileInput
    String taskModule
  }
    String fileName = basename(fileInput, ".tsv")
  command {
    set -eo pipefail
    python myscript.py --input ~{fileInput} --out ~{fileName}.summary.txt
  }
  runtime {
    cpu: 1
    modules: taskModule
  }
  output {
    File outputFile = "~{fileName}.summary.txt"
  }
}