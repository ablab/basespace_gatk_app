function launchSpec(dataProvider)
{
    return {
              commandLine:  ["python","/pipeline/pipeline.py", "basespace"],
              containerImageId:"gurevich/gatk3_best_practices",
              Options: ["bsfs.enabled=false"]
        }
}

function formUpdates(dataProvider)
{
    var is_custom_ref = ('' + dataProvider.GetProperty("input.select-ref")) === "1";
    if (is_custom_ref) {
        dataProvider.AttributeUpdates.Remove({ ElementId: "checkbox-full", AttributeName: "checked" });
    }
}
