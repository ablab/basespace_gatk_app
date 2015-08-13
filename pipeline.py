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
references = ['Ecoli', 'human']

def main(args):
    jsonfile = open('/data/input/AppSession.json')
    jsonObject = json.load(jsonfile)

    sampleID = []
    sampleHref = []
    sampleName = []
    sampleDir = []
    samples = []

    ref_num = 0
    is_human = False

    # parsing "Properties"
    for entry in jsonObject["Properties"]["Items"]:
        if "Name" not in entry:
            continue
        if entry["Name"] == "Input.project-id":
            project_id = entry["Content"]["Id"]
            config.basespace_output_dir = os.path.join(config.UROOT, project_id)

        # sample properties
        elif entry['Name'] == 'Input.Samples':
            for sample in range(len(entry['Items'])):
                sampleID.append(entry['Items'][sample]['Id'])
                sampleHref.append(entry['Items'][sample]['Href'])
                sampleName.append(entry['Items'][sample]['Name'])
                sampleDir = '/data/input/samples/%s/Data/Intensities/BaseCalls/' % sampleID[-1]
                if not os.path.exists(sampleDir):
                    sampleDir = '/data/input/samples/%s/' % sampleID[-1]
                for root, dirs, files in os.walk(str(sampleDir)):
                    R1files = fnmatch.filter(files,'*_R1_*')
                    R2files = fnmatch.filter(files,'*_R2_*')
                if len(R1files) != len(R2files):
                    print "number of R1 and R2 files do not match"
                    sys.exit()
                R1files.sort()
                R2files.sort()
                samples.append(map(lambda x, y: (sampleDir + x, sampleDir + y), R1files, R2files))
                config.basespace_output_dir = os.path.join(config.basespace_output_dir, entry['Items'][sample]['Name'])
        elif entry['Name'] == 'Input.select-ref':
            ref_num = int(entry['Content'])
            if ref_num == 1:
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

    ref_fpath = os.path.join(config.DIR_HOME, 'genomes', references[ref_num] + '.fa')
    # one sample per launch in multi-node mode
    bam_fpath = preprocessing.do(ref_fpath, samples[0], sampleID[0], config.temp_output_dir, config.basespace_output_dir, is_human)
    variant_calling.do(ref_fpath, sampleID[0], bam_fpath, config.temp_output_dir, config.basespace_output_dir, is_human)
    return

if __name__ == '__main__':
    main(sys.argv)