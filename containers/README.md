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