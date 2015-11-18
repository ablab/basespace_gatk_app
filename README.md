# GATK3 Best Practices App for BaseSpace
Source code and [BaseSpace](https://basespace.illumina.com) configs for [GATK3 Best Practices](https://www.broadinstitute.org/gatk/guide/best-practices) wrapping App.

Docker container with this sources is at https://hub.docker.com/r/gurevich/gatk3_best_practices/. It is based on **ubuntu** Docker. Additional installations are (after `sudo apt-get update`):
* python v2.7 (`sudo apt-get install python`)
* java v1.7 (`sudo apt-get install openjdk-7-jre-headless`)
* samtools v.0.1.19 (`sudo apt-get install samtools`)
* tabix and bgzip (`sudo apt-get install tabix`)

Deployment protocol and workflow overview are [here](https://docs.google.com/document/d/1W0y86a8eVzYiCcyY_PZD2seJ3oX-mD3_iQLGUg5LNmM/edit#).
