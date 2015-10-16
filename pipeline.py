#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################

import os
import shutil
import sys
import json

import config
references = ['/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
              '/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa']

def main(args):
    jsonfile = open('/data/input/AppSession.json')
    jsonObject = json.load(jsonfile)

    sampleIDs = []
    sample_names = []
    samples = []

    ref_num = 0
    project_id = 0
    ref_fpath = None

    config.temp_output_dir = config.OUT_ROOT

    if os.path.isdir(config.temp_output_dir):
        shutil.rmtree(config.temp_output_dir)
    os.makedirs(config.temp_output_dir)
    # parsing "Properties"
    for entry in jsonObject["Properties"]["Items"]:
        if "Name" not in entry:
            continue
        if entry["Name"] == "Input.project-id":
            project_id = entry["Content"]["Id"]
            config.basespace_output_dir = os.path.join(config.UROOT, project_id, 'Results')
            if os.path.isdir(config.basespace_output_dir):
                shutil.rmtree(config.basespace_output_dir)
            os.makedirs(config.basespace_output_dir)

        # sample properties
        elif entry['Name'] == 'Input.AppResults':
            for sample in range(len(entry['Items'])):
                sampleID = entry['Items'][sample]['Id']
                sampleDir = '/data/input/appresults/%s/' % sampleID
                for root, dirs, files in os.walk(str(sampleDir)):
                    for name in files:
                        if name.endswith('.bam'):
                            samples.append(os.path.join(root, name))
                            sampleIDs.append(entry['Items'][sample]['Id'])
                            sample_names.append(entry['Items'][sample]['Name'])
        elif entry['Name'] == 'Input.select-ref':
            ref_num = int(entry['Content'])
        elif entry['Name'] == 'Input.Files':
            for sample in range(len(entry['Items'])):
                fileID = entry['Items'][sample]['ParentAppResult']['Id']
                refDir = '/data/input/appresults/%s/' % fileID
                for root, dirs, files in os.walk(str(refDir)):
                    for name in files:
                        if name.endswith('.fasta') or name.endswith('.fa') or name.endswith('.fna'):
                            raw_ref_fpath = os.path.join(root, name)
                            ref_fpath = os.path.join(config.temp_output_dir, name)
                            shutil.copy(raw_ref_fpath, ref_fpath)
        elif entry['Name'] == 'Input.checkbox-pipeline':
            config.reduced_workflow = False

    if ref_num == 0 or not ref_fpath:
        ref_fpath = references[0]
    else:
        config.reduced_workflow = True

    config.threads = 24
    config.max_memory = 60

    # one sample per launch in multi-node mode
    import preprocessing
    bam_fpaths = preprocessing.do(ref_fpath, samples, sampleIDs, config.temp_output_dir, config.basespace_output_dir, project_id)
    import variant_calling
    variant_calling.do(ref_fpath, sampleIDs, bam_fpaths, config.temp_output_dir, config.basespace_output_dir, project_id, samples, sample_names)
    return

if __name__ == '__main__':
    main(sys.argv)