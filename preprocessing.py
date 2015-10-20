#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################

import os
import utils
import config

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

def all_required_binaries_exist(bin_dirpath, binary):
    if not os.path.isfile(os.path.join(bin_dirpath, binary)):
        return False
    return True


def process_single_sample(ref_fpath, sample_id, bam_fpath, scratch_dirpath, output_dirpath, project_id, num_threads):
    log_fpath = os.path.join(output_dirpath, sample_id + '.log')
    final_bam_fpath = os.path.join(output_dirpath, sample_id + '.bam')

    replace_rg_fpath = os.path.join(scratch_dirpath, sample_id + '.temp.bam')
    bqsr_fpath = os.path.join(scratch_dirpath, sample_id + '.bqsr.bam')

    targetintervals_fpath = os.path.join(scratch_dirpath, sample_id + '_realignment_targets.list')
    validate_log_fpath = os.path.join(output_dirpath, sample_id + '.validate.txt')
    ignored_errors = ['IGNORE=%s' % error for error in ignored_errors_in_bam]
    cmd = ['java', '-jar', config.picard_fpath, 'ValidateSamFile', 'INPUT=%s' % bam_fpath, 'OUTPUT=%s' % validate_log_fpath, 'IGNORE_WARNINGS=true',
           'MAX_OUTPUT=1'] + ignored_errors
    is_corrupted_file = utils.call_subprocess(cmd, stderr=open(log_fpath, 'a'))
    if is_corrupted_file:
        new_bam_fpath = os.path.join(scratch_dirpath, sample_id + '_sort.bam')
        utils.call_subprocess(['java', '-jar', config.picard_fpath, 'SortSam', 'INPUT=%s' % bam_fpath,
                               'OUTPUT=%s' % new_bam_fpath, 'SORT_ORDER=coordinate', 'VALIDATION_STRINGENCY=LENIENT',
                               'CREATE_INDEX=true'], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', config.picard_fpath, 'AddOrReplaceReadGroups', 'INPUT=%s' % new_bam_fpath,
                              'OUTPUT=%s' % replace_rg_fpath, 'RGPL=illumina', 'RGSM=%s' % sample_id,
                               'RGLB=lib', 'RGPU=adapter', 'VALIDATION_STRINGENCY=LENIENT',
                               'CREATE_INDEX=true'], stderr=open(log_fpath, 'a'))
        bam_fpath = replace_rg_fpath
    if not os.path.exists(ref_fpath + '.fai'):
        print 'Preparing reference file...'
        samtools_fpath = os.path.join(config.samtools_dirpath, 'samtools')
        utils.call_subprocess([samtools_fpath, 'faidx', ref_fpath], stderr=open(log_fpath, 'a'))
        get_chr_lengths(ref_fpath)
        ref_fname, _ = os.path.splitext(ref_fpath)
        utils.call_subprocess(['java', '-jar', config.picard_fpath, 'CreateSequenceDictionary', 'R=%s' % ref_fpath,
                               'O=%s' % ref_fname + '.dict'], stderr=open(log_fpath, 'a'))

    print 'Realign indels...'
    cmd = ['java', '-jar', config.gatk_fpath, '-T', 'RealignerTargetCreator', '-R', ref_fpath, '-nt', num_threads,
                            '-I', bam_fpath, '-o', targetintervals_fpath]
    if not config.reduced_workflow:
        cmd += ['-known', config.known_fpath, '-known', config.tg_indels_fpath]
    utils.call_subprocess(cmd, stderr=open(log_fpath, 'a'))
    cmd = ['java', '-jar', config.gatk_fpath, '-T', 'IndelRealigner', '-R', ref_fpath,
                            '-I', bam_fpath, '-targetIntervals',  targetintervals_fpath, '-o', final_bam_fpath]
    if not config.reduced_workflow:
        cmd += ['-known', config.known_fpath, '-known', config.tg_indels_fpath]
    utils.call_subprocess(cmd, stderr=open(log_fpath, 'a'))

    if not config.reduced_workflow:
        print 'Recalibrate bases...'
        recaltable_fpath = os.path.join(scratch_dirpath, sample_id + '.table')
        utils.call_subprocess(['java', '-jar', config.gatk_fpath, '-T', 'BaseRecalibrator', '-R', ref_fpath, '-nct', num_threads,
                                '-I', final_bam_fpath, '-knownSites', config.dbsnp_fpath, '-dt', 'ALL_READS', '-dfrac', '0.10 ',
                               '-o', recaltable_fpath], stderr=open(log_fpath, 'a'))
        utils.call_subprocess(['java', '-jar', config.gatk_fpath, '-T', 'PrintReads', '-R', ref_fpath, '-nct', num_threads,
                                '-I', final_bam_fpath, '-BQSR', recaltable_fpath,
                                '-o', bqsr_fpath], stderr=open(log_fpath, 'a'))

    print 'Building BAM index...'
    utils.call_subprocess(['java', '-jar', config.picard_fpath, 'BuildBamIndex', 'INPUT=%s' % final_bam_fpath, 'VALIDATION_STRINGENCY=LENIENT',
                           'OUTPUT=%s' % final_bam_fpath + '.bai'], stderr=open(log_fpath, 'a'))
    return final_bam_fpath


def get_chr_lengths(ref_fpath):
    ref_index_file = open(ref_fpath + '.fai')
    config.chr_lengths = {}
    config.chr_names = []
    for line in ref_index_file.read().split('\n'):
        if line:
            line = line.split()
            config.chr_names.append(line[0])
            config.chr_lengths[line[0]] = line[1]

def do(ref_fpath, samples, sample_ids, scratch_dirpath, output_dirpath, project_id):
    from libs.joblib import Parallel, delayed
    n_jobs = min(len(samples), config.max_threads)
    num_threads = max(1, config.max_threads//n_jobs)
    final_bam_fpaths = Parallel(n_jobs=n_jobs)(delayed(process_single_sample)(ref_fpath, sample_ids[i], samples[i], scratch_dirpath, output_dirpath, project_id, str(num_threads))
                                               for i in range(len(samples)))
    return final_bam_fpaths
