#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################
import getopt

import os
import shutil
import sys
import json

import config
import utils
from os.path import join

references = ['/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
              '/genomes/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa']

def parse_basespace_input(project_id):
    ref_num = 0
    ref_fpath = None

    samples = []
    sample_ids = []
    sample_names = []
    config.INPUT_DIR = config.INPUT_DIR_CLOUD
    json_fname = join(config.INPUT_DIR, 'AppSession.json')
    if not os.path.isfile(json_fname):
        config.INPUT_DIR = config.INPUT_DIR_LOCAL
        json_fname = join(config.INPUT_DIR, 'AppSession.json')
    json_file = open(json_fname)
    json_object = json.load(json_file)
    config.temp_output_dir = config.SCRATCH_DIR

    if os.path.isdir(config.temp_output_dir):
        shutil.rmtree(config.temp_output_dir)
    os.makedirs(config.temp_output_dir)
    # parsing "Properties"
    for entry in json_object["Properties"]["Items"]:
        if "Name" not in entry:
            continue
        if entry["Name"] == "Input.project-id":
            project_id = entry["Content"]["Id"]
            config.output_dir = join(config.RESULTS_DIR, project_id, 'Results')
            if os.path.isdir(config.output_dir):
                shutil.rmtree(config.output_dir)
            os.makedirs(config.output_dir)

        # sample properties
        elif entry['Name'] == 'Input.AppResults':
            for sample in range(len(entry['Items'])):
                sample_id = entry['Items'][sample]['Id']
                sample_dir = join(config.INPUT_DIR, 'appresults', sample_id)
                for root, dirs, files in os.walk(str(sample_dir)):
                    for name in files:
                        if name.endswith('.bam'):
                            samples.append(join(root, name))
                            sample_ids.append(sample_id)
                            sample_names.append(entry['Items'][sample]['Name'])
        elif entry['Name'] == 'Input.select-ref':
            ref_num = int(entry['Content'])
        elif entry['Name'] == 'Input.Files':
            for sample in range(len(entry['Items'])):
                file_id = entry['Items'][sample]['ParentAppResult']['Id']
                ref_dir = join(config.INPUT_DIR, 'appresults', file_id)
                for root, dirs, files in os.walk(str(ref_dir)):
                    for name in files:
                        if name.endswith('.fasta') or name.endswith('.fa') or name.endswith('.fna'):
                            raw_ref_fpath = join(root, name)
                            ref_fpath = join(config.temp_output_dir, name)
                            shutil.copy(raw_ref_fpath, ref_fpath)
        elif entry['Name'] == 'Input.checkbox-full':
            config.reduced_workflow = False
        elif entry['Name'] == 'Input.checkbox-lowemit':
            config.low_emit = True

    if ref_num == 0 or not ref_fpath:
        ref_fpath = references[0]
    else:
        config.reduced_workflow = True

    config.max_threads = 24
    config.max_memory = 60
    return ref_fpath, samples, sample_ids, sample_names, project_id

def main(args):
    project_id = 'project'
    if args[0] == 'basespace':
        ref_fpath, samples_fpaths, sample_ids, sample_names, project_id = parse_basespace_input(project_id)
    else:
        ref_fpath = None
        sample_ids = []
        sample_names = []
        try:
            options, samples_fpaths = getopt.gnu_getopt(args, config.short_options, config.long_options)
        except getopt.GetoptError:
            _, exc_value, _ = sys.exc_info()
            print >> sys.stderr, exc_value
            print >> sys.stderr
            sys.exit(2)

        for opt, arg in options[:]:
            if opt in ('-S', "-samples"):
                sample_ids = arg.split(',')
            elif opt in ('-R', "--reference"):
                ref_fpath = utils.assert_file_exists(arg, 'reference')
            elif opt in ('-f', "--full"):
                config.reduced_workflow = False
            elif opt in ('-l', "--low-emit"):
                config.low_emit = True
            elif opt in ('-t', "--threads"):
                config.max_threads = int(arg)
            elif opt in ('-m', "--memory"):
                config.max_memory = int(arg)

        if not ref_fpath:
            print 'Error! You should specify the reference file with option -R'
            sys.exit(1)

    config.max_single_gatk_mem = config.max_memory // config.max_threads
    if config.max_single_gatk_mem < 2:
        config.max_single_gatk_mem = 2
        config.max_threads = max(1, config.max_memory/2)
    config.max_gatk_threads = config.max_threads
    utils.check_gatk()
    utils.check_samtools()

    # one sample per launch in multi-node mode
    import preprocessing
    bam_fpaths = preprocessing.do(ref_fpath, samples_fpaths, sample_ids, config.temp_output_dir, config.output_dir)
    import variant_calling
    variant_calling.do(ref_fpath, sample_ids, bam_fpaths, config.temp_output_dir, config.output_dir, project_id, samples_fpaths, sample_names)
    return

if __name__ == '__main__':
    main(sys.argv[1:])