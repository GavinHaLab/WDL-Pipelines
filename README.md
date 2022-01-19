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
