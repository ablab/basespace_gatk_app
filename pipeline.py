#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################

import fnmatch
import os
import shutil
import sys
import json

import config
import variant_calling
import preprocessing
references = ['Ecoli', 'Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa', 'Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa']

def main(args):
    jsonfile = open('/data/input/AppSession.json')
    jsonObject = json.load(jsonfile)

    sampleIDs = []
    sampleNames = []
    sampleHref = []
    sampleName = []
    sampleDir = []
    samples = []

    ref_num = 0
    project_id = 0
    is_human = False

    # parsing "Properties"
    for entry in jsonObject["Properties"]["Items"]:
        if "Name" not in entry:
            continue
        if entry["Name"] == "Input.project-id":
            project_id = entry["Content"]["Id"]
            config.basespace_output_dir = os.path.join(config.UROOT, project_id, 'Results')

        # sample properties
        elif entry['Name'] == 'Input.Files':
            for sample in range(len(entry['Items'])):
                sampleID = entry['Items'][sample]['ParentAppResult']['Id']
                sampleDir = '/data/input/appresults/%s/' % sampleID
                for root, dirs, files in os.walk(str(sampleDir)):
                    for name in files:
                        if name.endswith('.bam'):
                            samples.append(os.path.join(root, name))
                            sampleIDs.append(entry['Items'][sample]['ParentAppResult']['Id'])
                            sampleNames.append(entry['Items'][sample]['ParentAppResult']['Name'])
        elif entry['Name'] == 'Input.select-ref':
            ref_num = int(entry['Content'])
            if ref_num > 0:
                is_human = True

    config.temp_output_dir = config.OUT_ROOT
    config.threads = 32
    config.memory = 60

    if os.path.isdir(config.basespace_output_dir):
        shutil.rmtree(config.basespace_output_dir)

    os.makedirs(config.basespace_output_dir)

    if os.path.isdir(config.temp_output_dir):
        shutil.rmtree(config.temp_output_dir)
    os.makedirs(config.temp_output_dir)

    ref_fpath = os.path.join('/genomes', references[ref_num])
    # one sample per launch in multi-node mode

    bam_fpaths = preprocessing.do(ref_fpath, samples, sampleIDs, config.temp_output_dir, config.basespace_output_dir, is_human, project_id)
    variant_calling.do(ref_fpath, sampleIDs, bam_fpaths, config.temp_output_dir, config.basespace_output_dir, is_human, project_id)
    return

if __name__ == '__main__':
    main(sys.argv)
