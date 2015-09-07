#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################

import os
import shlex
import config
import utils

picard_dirpath = os.path.join(config.LIBS_LOCATION, 'picard')
gatk_dirpath = os.path.join(config.LIBS_LOCATION, 'gatk')

def gatk_fpath():
    return os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

def picard_fpath():
    return os.path.join(picard_dirpath, 'picard.jar')


def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True


def process_single_sample(ref_fpath, sampleID, bam_fpath, scratch_dirpath, output_dirpath, is_human, project_id, num_threads):
    log_fpath = os.path.join(output_dirpath, sampleID + '.log')
    final_bam_fpath = os.path.join(output_dirpath, sampleID + '.bam')

    realignedbam_fpath = os.path.join(scratch_dirpath, sampleID + '.realigned.bam')
    recaltable_fpath = os.path.join(scratch_dirpath, sampleID + '.table')
    post_recaltable_fpath = os.path.join(scratch_dirpath, sampleID + '_post.table')

    csv_fpath = os.path.join(output_dirpath, sampleID + '_recalibration_plots.csv')

    if is_human:
        targetintervals_fpath = os.path.join(scratch_dirpath, sampleID + '_realignment_targets.list')
        known_fpath = os.path.join(config.DIR_HOME, 'genomes', 'gold_indels.vcf')
        dbsnp_fpath = '/genomes/Homo_sapiens/UCSC/hg19/Annotation/dbsnp_132.hg19.vcf'
        print 'Realign indels...'
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'RealignerTargetCreator', '-R', ref_fpath, '-nt', num_threads, 
                                '-I', bam_fpath, '-known', known_fpath, '-o', targetintervals_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'IndelRealigner', '-R', ref_fpath,
                                '-I', bam_fpath, '-targetIntervals',  targetintervals_fpath, '-known', known_fpath,
                              '-o', realignedbam_fpath], stderr=open(log_fpath, 'a'))
        print 'Recalibrate bases...'
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'BaseRecalibrator', '-R', ref_fpath, '-nct', num_threads, 
                                '-I', realignedbam_fpath, '-knownSites', dbsnp_fpath, '-knownSites', known_fpath,
                              '-o', recaltable_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'BaseRecalibrator', '-R', ref_fpath, '-nct', num_threads, 
                                '-I', realignedbam_fpath, '-knownSites',  dbsnp_fpath, '-knownSites', known_fpath,
                              '-BQSR', recaltable_fpath, '-o', post_recaltable_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'AnalyzeCovariates', '-R', ref_fpath,
                                '-before', recaltable_fpath, '-after',  post_recaltable_fpath,
                              '-csv', csv_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'PrintReads', '-R', ref_fpath, '-nct', num_threads,
                                '-I', realignedbam_fpath, '-BQSR', recaltable_fpath,
                              '-o', final_bam_fpath], stderr=open(log_fpath, 'a'))
    else:
        final_bam_fpath = bam_fpath
    print 'Building BAM index...'
    utils.call_subprocess(['java', '-jar', picard_fpath(), 'BuildBamIndex', 'INPUT=%s' % final_bam_fpath,
                           'OUTPUT=%s' % final_bam_fpath + '.bai'], stderr=open(log_fpath, 'a'))
    return final_bam_fpath


def do(ref_fpath, samples, sampleIDs, scratch_dirpath, output_dirpath, is_human, project_id):
    from libs.joblib import Parallel, delayed
    n_jobs = min(len(samples), config.threads)
    num_threads = min(config.max_gatk_threads, config.threads//n_jobs)
    final_bam_fpaths = Parallel(n_jobs=n_jobs)(delayed(process_single_sample)(ref_fpath, sampleIDs[i], samples[i], scratch_dirpath, output_dirpath, is_human, project_id, str(num_threads))
                                               for i in range(len(samples)))
    return final_bam_fpaths
