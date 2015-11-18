function launchSpec(dataProvider)
{
    return {
              commandLine:  ["python","/pipeline/pipeline.py", "basespace"],
              containerImageId:"gurevich/gatk3_best_practices",
              Options: ["bsfs.enabled=false"]
        }
}

