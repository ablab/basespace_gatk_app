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

bwa_dirpath = os.path.join(config.LIBS_LOCATION, 'bwa')
picard_dirpath = os.path.join(config.LIBS_LOCATION, 'picard')
samtools_dirpath = os.path.join(config.LIBS_LOCATION, 'samtools')
gatk_dirpath = os.path.join(config.LIBS_LOCATION, 'gatk')

def gatk_fpath():
    return os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

def samtools_fpath():
    return os.path.join(samtools_dirpath, 'samtools')

def bwa_fpath():
    return os.path.join(bwa_dirpath, 'bwa')

def picard_fpath():
    return os.path.join(picard_dirpath, 'picard.jar')


def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True


def prepare_reference(ref_fpath, output_dirpath):
    log_fpath = output_dirpath + 'ref.log'
    index_algorithm = 'is'
    if os.path.getsize(ref_fpath) > 1000000000:
        index_algorithm = 'bwtsw'
    utils.call_subprocess([bwa_fpath(), 'index', '-a', index_algorithm, ref_fpath], stderr=open(log_fpath, 'w'))
    if not os.path.exists(ref_fpath + '.fai'):
        utils.call_subprocess([samtools_fpath(), 'faidx', ref_fpath], stderr=open(log_fpath, 'a'))
    if not os.path.exists(os.path.splitext(ref_fpath)[0] + '.dict'):
        utils.call_subprocess(['java', '-jar', picard_fpath(), 'CreateSequenceDictionary', 'REFERENCE=%s' % ref_fpath,
                           'OUTPUT=%s' % os.path.splitext(ref_fpath)[0] + '.dict'], stderr=open(log_fpath, 'a'))
    return


def process_single_library(ref_fpath, sampleID, bwa_threads, reads_fpath, sam_fpath, bam_fpath, bam_rg_fpath, log_fpath):
    # mapping reads
    cmd = bwa_fpath() + ' mem -M -t %s %s %s' % (str(bwa_threads), ref_fpath, ' '.join(reads_fpath))
    utils.call_subprocess(shlex.split(cmd), stdout=open(sam_fpath, 'w'), stderr=open(log_fpath, 'w'))
    utils.call_subprocess(['java', '-jar', picard_fpath(), 'SortSam', 'INPUT=%s' % sam_fpath, 'OUTPUT=%s' % bam_fpath,
                           'SORT_ORDER=coordinate'], stderr=open(log_fpath, 'a'))
    file_name = os.path.splitext(os.path.basename(reads_fpath[0]))[0]
    rg_identifier = file_name.split('_')
    rg_lib = rg_identifier[-3]
    rg_unit = os.path.splitext(rg_identifier[-1])[0]
    utils.call_subprocess(['java', '-jar', picard_fpath(), 'AddOrReplaceReadGroups', 'INPUT=%s' % bam_fpath,
                            'OUTPUT=%s' % bam_rg_fpath, 'RGLB=%s' % rg_lib, 'RGPL=illumina',
                           'RGPU=%s' % rg_unit, 'RGSM=%s' % sampleID], stderr=open(log_fpath, 'a'))


def process_single_sample(ref_fpath, sampleID, reads_fpaths, scratch_dirpath, output_dirpath, is_human):
    res_fpath = os.path.join(output_dirpath, sampleID + '_alignment_stats.txt')
    log_fpath = os.path.join(output_dirpath, sampleID + '.log')


    bam_fpath = os.path.join(scratch_dirpath, sampleID + '.temp.bam')
    sam_fpaths = [os.path.join(scratch_dirpath, sampleID + str(i+1) + '.sam') for i in range(len(reads_fpaths))]
    bam_fpaths = [os.path.join(scratch_dirpath, sampleID + str(i+1) + '.temp.bam') for i in range(len(reads_fpaths))]
    bam_rg_fpaths = [os.path.join(scratch_dirpath, sampleID + str(i+1) + '.bam') for i in range(len(reads_fpaths))]
    dedup_fpath = os.path.join(scratch_dirpath, sampleID + '.dedup.bam')
    metrics_fpath = os.path.join(output_dirpath, sampleID + '_metrics.txt')

    realignedbam_fpath = os.path.join(scratch_dirpath, sampleID + '.realigned.bam')
    recaltable_fpath = os.path.join(scratch_dirpath, sampleID + '.table')
    post_recaltable_fpath = os.path.join(scratch_dirpath, sampleID + '_post.table')

    plots_fpath = os.path.join(output_dirpath, sampleID + '_recalibration_plots.table')
    finalbam_fpath = os.path.join(output_dirpath, sampleID + '.bam')

    n_jobs = min(config.threads, len(reads_fpaths))
    bwa_threads = max(1, config.threads // len(reads_fpaths))
    print 'Running BWA...'
    from libs.joblib import Parallel, delayed
    Parallel(n_jobs=n_jobs)(delayed(process_single_library)(ref_fpath, sampleID, bwa_threads, reads_fpath,
                                                sam_fpaths[index], bam_fpaths[index], bam_rg_fpaths[index], log_fpath) for index, reads_fpath in enumerate(reads_fpaths))

    if len(reads_fpaths) > 1:
        cmd = samtools_fpath() + ' merge %s %s' % (bam_fpath, ' '.join(bam_rg_fpaths))
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'w'))
    else:
        bam_fpath = bam_fpaths[0]
    utils.call_subprocess(['java', '-jar', picard_fpath(), 'MarkDuplicates', 'INPUT=%s' % bam_fpath,
                            'OUTPUT=%s' % (dedup_fpath if is_human else finalbam_fpath), 'METRICS_FILE=%s' % metrics_fpath], stderr=open(log_fpath, 'a'))

    print 'Alignment is finished.'
    if is_human:
        targetintervals_fpath = os.path.join(scratch_dirpath, sampleID + '_realignment_targets.list')
        known_fpath = os.path.join(config.DIR_HOME, 'genomes', 'gold_indels.vcf')
        dbsnp_fpath = os.path.join(config.DIR_HOME, 'genomes', 'dbsnp.vcf')
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'RealignerTargetCreator', '-R', ref_fpath,
                                '-I', dedup_fpath, '-known', known_fpath, '-o', targetintervals_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'IndelRealigner', '-R', ref_fpath,
                                '-I', dedup_fpath, '-targetIntervals',  targetintervals_fpath, '-known', known_fpath,
                              '-o', realignedbam_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'BaseRecalibrator', '-R', ref_fpath,
                                '-I', realignedbam_fpath, '-knownSites', dbsnp_fpath, '-knownSites', known_fpath,
                              '-o', recaltable_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'BaseRecalibrator', '-R', ref_fpath,
                                '-I', realignedbam_fpath, '-knownSites',  dbsnp_fpath, '-knownSites', known_fpath,
                              '-BQSR', recaltable_fpath, '-o', post_recaltable_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'AnalyzeCovariates', '-R', ref_fpath,
                                '-before', recaltable_fpath, '-after',  post_recaltable_fpath,
                              '-plots', plots_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'PrintReads', '-R', ref_fpath,
                                '-I', realignedbam_fpath, '-BQSR', recaltable_fpath,
                              '-o', finalbam_fpath], stderr=open(log_fpath, 'a'))

    utils.call_subprocess(['java', '-jar', picard_fpath(), 'BuildBamIndex', 'INPUT=%s' % finalbam_fpath,
                           'OUTPUT=%s' % finalbam_fpath + '.bai'], stderr=open(log_fpath, 'a'))

    utils.call_subprocess(['java', '-jar', picard_fpath(), 'CollectAlignmentSummaryMetrics', 'INPUT=%s' % finalbam_fpath,
                            'OUTPUT=%s' % res_fpath, 'R=%s' % ref_fpath], stderr=open(log_fpath, 'a'))
    return finalbam_fpath


def do(ref_fpath, reads_fpaths, sampleID, scratch_dirpath, output_dirpath, is_human):
    prepare_reference(ref_fpath, scratch_dirpath)
    bam_fpath = process_single_sample(ref_fpath, sampleID, reads_fpaths,scratch_dirpath, output_dirpath, is_human)
    return bam_fpath