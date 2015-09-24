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
ignored_errors_in_bam = ["INVALID_QUALITY_FORMAT","INVALID_FLAG_PROPER_PAIR","INVALID_FLAG_MATE_UNMAPPED",
                         "MISMATCH_FLAG_MATE_UNMAPPED","INVALID_FLAG_MATE_NEG_STRAND","MISMATCH_FLAG_MATE_NEG_STRAND",
                         "INVALID_FLAG_FIRST_OF_PAIR","INVALID_FLAG_SECOND_OF_PAIR","PAIRED_READ_NOT_MARKED_AS_FIRST_OR_SECOND",
                         "INVALID_FLAG_NOT_PRIM_ALIGNMENT","INVALID_FLAG_SUPPLEMENTARY_ALIGNMENT","INVALID_FLAG_READ_UNMAPPED",
                         "INVALID_INSERT_SIZE","INVALID_MAPPING_QUALITY","INVALID_CIGAR","ADJACENT_INDEL_IN_CIGAR","INVALID_MATE_REF_INDEX",
                         "MISMATCH_MATE_REF_INDEX","INVALID_REFERENCE_INDEX","INVALID_ALIGNMENT_START","MISMATCH_MATE_ALIGNMENT_START",
                         "MATE_FIELD_MISMATCH","INVALID_TAG_NM","MISSING_TAG_NM","MISSING_HEADER","MISSING_SEQUENCE_DICTIONARY","RECORD_OUT_OF_ORDER",
                         "RECORD_MISSING_READ_GROUP","INVALID_INDEXING_BIN","MISSING_VERSION_NUMBER","INVALID_VERSION_NUMBER","TRUNCATED_FILE",
                         "MISMATCH_READ_LENGTH_AND_QUALS_LENGTH","EMPTY_READ","CIGAR_MAPS_OFF_REFERENCE","MISMATCH_READ_LENGTH_AND_E2_LENGTH",
                         "MISMATCH_READ_LENGTH_AND_U2_LENGTH","E2_BASE_EQUALS_PRIMARY_BASE","BAM_FILE_MISSING_TERMINATOR_BLOCK","UNRECOGNIZED_HEADER_TYPE",
                         "POORLY_FORMATTED_HEADER_TAG","HEADER_TAG_MULTIPLY_DEFINED","HEADER_RECORD_MISSING_REQUIRED_TAG","INVALID_DATE_STRING",
                         "TAG_VALUE_TOO_LARGE","INVALID_INDEX_FILE_POINTER","INVALID_PREDICTED_MEDIAN_INSERT_SIZE","DUPLICATE_PROGRAM_GROUP_ID",
                         "MATE_NOT_FOUND","MATES_ARE_SAME_END","MISMATCH_MATE_CIGAR_STRING","MATE_CIGAR_STRING_INVALID_PRESENCE"]

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

    replace_rg_fpath = os.path.join(scratch_dirpath, sampleID + '.temp.bam')
    realignedbam_fpath = os.path.join(scratch_dirpath, sampleID + '.realigned.bam')
    recaltable_fpath = os.path.join(scratch_dirpath, sampleID + '.table')
    post_recaltable_fpath = os.path.join(scratch_dirpath, sampleID + '_post.table')

    csv_fpath = os.path.join(output_dirpath, sampleID + '_recalibration_plots.csv')

    if is_human:
        targetintervals_fpath = os.path.join(scratch_dirpath, sampleID + '_realignment_targets.list')
        known_fpath = os.path.join(config.DIR_HOME, 'genomes', 'gold_indels.vcf')
        dbsnp_fpath = '/genomes/Homo_sapiens/UCSC/hg19/Annotation/dbsnp_132.hg19.vcf'
        validate_log_fpath = os.path.join(output_dirpath, sampleID + '.validate.txt')
        ignored_errors = ['IGNORE=%s' % error for error in ignored_errors_in_bam]
        ignored_errors = ' '.join(ignored_errors)
        picard_path = picard_fpath()
        cmd = ('java -jar {picard_path} ValidateSamFile INPUT={bam_fpath} OUTPUT={validate_log_fpath} IGNORE_WARNINGS=true MAX_OUTPUT=1 '
               '{ignored_errors}').format(**locals())
        is_corrupted_file = utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))
        if is_corrupted_file:
            new_bam_fpath = os.path.join(scratch_dirpath, sampleID + '_sort.bam')
            utils.call_subprocess(['java', '-jar', picard_fpath(), 'SortSam', 'INPUT=%s' % bam_fpath,
                                   'OUTPUT=%s' % new_bam_fpath, 'SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=LENIENT',
                                   'CREATE_INDEX=true'], stderr=open(log_fpath, 'a'))
            utils.call_subprocess(['java', '-jar', picard_fpath(), 'AddOrReplaceReadGroups', 'INPUT=%s' % new_bam_fpath,
                               'OUTPUT=%s' % replace_rg_fpath, 'RGPL=illumina', 'RGSM=%s' % sampleID,
                                   'RGLB=lib', 'RGPU=adapter', 'VALIDATION_STRINGENCY=LENIENT',
                                   'CREATE_INDEX=true'], stderr=open(log_fpath, 'a'))
            bam_fpath = replace_rg_fpath
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
