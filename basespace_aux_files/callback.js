function launchSpec(dataProvider)
{
    return {
              commandLine:  ["python","/pipeline/pipeline.py", "basespace"],
              containerImageId:"almiheenko/gatk3best",
              Options: ["bsfs.enabled=false"]
        }
}

