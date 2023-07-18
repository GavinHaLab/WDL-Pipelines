# wdl-repo
A central repository for tools wrapped in WDL for the Ha lab.  

## TGR data management
This contains support documentation for working with data annotation and tracking processes supported by the Fred Hutch Translational Genomics Repository.  

## containers
For easier transitions between computing infrastructure and workflow sharing, ideally workflow tasks can be run in the simplest docker container feasible.  The documentation of custom docker containers used in workflows in this repo should be stored, container-by-container in folders in this folder.  Ideally a dockerfile for each would be included in these folders, and information in the README about where docker containers can be found would be provided. 


## workflows
This folder can house orchestrating workflows.  Perhaps these are WDLs that import commonly-occuring processes (subworkflows containing tasks that are always performed together) from the `processes` folder.  Otherwise they are workflows with defined tasks in the WDL files with either templates for inputs or parameter collections or the filled in input files for the group.  One folder per workflow is ideal, and in each folder one or more WDL files, one or more input json's, and potentially a README explaining the suggested use of the workflow would be considered fully documented.  

## processes
This folder contains commonly used processes that are themselves complete WDL workflow definitions that describe processes that link several tasks together to create an output.  Often these are small workflows which create intermediate files that are very unlikely to ever be used directly and thus not passed back the main workflows that call these subworkflows.  
> NOTE:  Some of these workflows include a custome runtime parameter `dockerSL`.  These workflows can be run on small (<20GB) input files with Cromwell in general, or on other platforms such as AWS, Azure or Terra, IF you exchange `docker` for `dockerSL` in the files.  For those running workflows using the Fred Hutch cluster and the Fred Hutch configuration of Cromwell, these files can be run as is if the Cromwell configuration [version v1.1](https://github.com/FredHutch/diy-cromwell-server/releases/tag/v1.1) is used.  Keep in mind, these workflows can be run with `docker` even on larger files, but they will take substantially longer.  There is no other difference in the results as the difference is solely a processing improvement not anything influencing the tasks themselves.  

# Docker Containers in Use


Beyond those in this folder, Docker containers and their documentation from external providers are listed here:

- GATK:  broadinstitute/gatk:4.2.2.0; https://github.com/broadinstitute/gatk/releases/tag/4.2.2.0

- GATK with BWA for piped BWA tasks: fredhutch/gatk-bwa:4.2.20-0.7.17, https://github.com/FredHutch/gatk-bwa


## To update a docker container:
- Clone this repository to your local computer where Docker is installed
- In Terminal change to the directory where the dockerfile you want to build is
- Do:
```
docker build . -t vortexing/ichorcna:v0.5.0
```
> Here, replace `vortexing` with your DockerHub username or the institution name you will be pushing to docker hub with, and `ichorcna` with the tool name that you're focusing on in this container.  Please also include the `:v0.5.0` with an appropriate version name for your software that will specify to users exactly what version to expect inside of the container. 

- Fix issues
- Push container to DockerHub using:
```
docker push vortexing/ichorcna:v0.5.0
```

Now you're ready to use the container in your runtime block of your WDLs!
