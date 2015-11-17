function launchSpec(dataProvider)
{
    return {
              commandLine:  ["python","/pipeline/pipeline.py", "basespace"],
              containerImageId:"gurevich/gatk3best",
              Options: ["bsfs.enabled=false"]
        }
}

