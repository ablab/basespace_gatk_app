function launchSpec(dataProvider)
{
    return {
              commandLine:  ["python","/pipeline/pipeline.py"],
              containerImageId:"almiheenko/gatk3best",
              Options: ["bsfs.enabled=true"]
        }
}

