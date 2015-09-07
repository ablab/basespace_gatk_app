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
from libs.joblib import Parallel, delayed
gatk_dirpath = os.path.join(config.LIBS_LOCATION, 'gatk')
chr_names_hg19 = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                  'chr21', 'chr22', 'chrX', 'chrY', 'chrM']

def gatk_fpath():
    return os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

def process_single_file(ref_fpath, sampleID, bam_fpath, output_dirpath, scratch_dirpath):
    log_fpath = os.path.join(output_dirpath, sampleID + '.log')
    n_jobs = min(len(chr_names_hg19), config.max_gatk_threads)
    raw_vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(process_single_chr)(ref_fpath, sampleID, bam_fpath, scratch_dirpath, log_fpath, chr)
                                               for chr in chr_names_hg19)
    return raw_vcf_fpaths

def merge_vcfs(output_dirpath, sampleID, raw_vcf_fpaths, ref_fpath):
    merge_vcf_fpath = os.path.join(output_dirpath, sampleID + '.g.vcf')
    log_fpath = os.path.join(output_dirpath, sampleID + '.log')
    variants = '--variant ' + ' --variant '.join(raw_vcf_fpaths)
    gatk_path = gatk_fpath()
    cmd = ("java -jar {gatk_path} -T CombineGVCFs -R {ref_fpath} {variants} "
               "-o {merge_vcf_fpath}").format(**locals())
    utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))
    return merge_vcf_fpath

def process_single_chr(ref_fpath, sampleID, bam_fpath, output_dirpath, log_fpath, chr):
    raw_g_vcf_fpath = os.path.join(output_dirpath, sampleID + chr + '.g.vcf')
    utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'HaplotypeCaller', '-R', ref_fpath, '-L', chr,
                            '-I', bam_fpath, '-stand_emit_conf', '10',
                          '-stand_call_conf', '30', '-ERC', 'GVCF', '-variant_index_type', 'LINEAR',
                           '-variant_index_parameter', '128000', '-o', raw_g_vcf_fpath], stderr=open(log_fpath, 'a'))
    return raw_g_vcf_fpath

def process_files(ref_fpath, sampleIDs, bam_fpaths, scratch_dirpath, output_dirpath, is_human, project_id):
    print 'Calling variants...'
    raw_vcf_fpaths = [process_single_file(ref_fpath, sampleIDs[i], bam_fpaths[i], output_dirpath, scratch_dirpath)
                                               for i in range(len(bam_fpaths))]
    n_jobs = min(len(raw_vcf_fpaths), config.max_gatk_threads)
    raw_vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(merge_vcfs)(output_dirpath, sampleIDs[i], raw_vcf_fpaths[i], ref_fpath)
                                               for i in range(len(raw_vcf_fpaths)))
    vcf_fpath = os.path.join(output_dirpath, project_id + '.vcf')
    if is_human:
        gatk_path = gatk_fpath()
        log_fpath = os.path.join(output_dirpath, project_id + '.log')
        raw_vcf_fpath = os.path.join(scratch_dirpath, project_id + '.vcf')
        hapmap_fpath = os.path.join(config.DIR_HOME, 'genomes', 'hapmap_3.3.vcf')
        omni_fpath = os.path.join(config.DIR_HOME, 'genomes', '1000G_omni.vcf')
        dbsnp_fpath = '/genomes/Homo_sapiens/UCSC/hg19/Annotation/dbsnp_132.hg19.vcf'
        mills_fpath = os.path.join(config.DIR_HOME, 'genomes', 'gold_indels.vcf')

        recal_fpath = os.path.join(scratch_dirpath, project_id + '_SNP.recal')
        tranches_fpath = os.path.join(scratch_dirpath, project_id + '_SNP.tranches')

        raw_indels_vcf_fpath = os.path.join(scratch_dirpath, project_id + '_raw_indels.vcf')
        recal_indel_fpath = os.path.join(scratch_dirpath, project_id + '_INDEL.recal')
        tranches_indel_fpath = os.path.join(scratch_dirpath, project_id + '_INDEL.tranches')
        variants = '--variant ' + ' --variant '.join(raw_vcf_fpaths)
        print 'Joint genotyping...'
        num_threads = str(config.max_gatk_threads)
        cmd = ("java -jar {gatk_path} -T GenotypeGVCFs -R {ref_fpath} {variants} -nt {num_threads}"
               "-o {raw_vcf_fpath}").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

        print 'Filtering variants...'
        # variant filtering
        cmd = ("java -jar {gatk_path} -T VariantRecalibrator -R {ref_fpath} -input {raw_vcf_fpath} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "
               " {hapmap_fpath} -resource:omni,known=false,training=true,truth=true,prior=12.0 {omni_fpath} "
               "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp_fpath} -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum "
               "-mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile {recal_fpath} "
               "-tranchesFile {tranches_fpath}").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

        cmd = ("java -jar {gatk_path} -T ApplyRecalibration -R {ref_fpath} -input {raw_vcf_fpath} -mode SNP --ts_filter_level 99.0 "
               "-recalFile  {recal_fpath} -tranchesFile {tranches_fpath} -o {raw_indels_vcf_fpath}").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

        cmd = ("java -jar {gatk_path} -T VariantRecalibrator -R {ref_fpath} -input {raw_indels_vcf_fpath} -resource:mills,known=true,training=true,truth=true,prior=12.0 "
               "{mills_fpath} -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0"
               " -tranche 90.0 --maxGaussians 4 -recalFile {recal_indel_fpath} -tranchesFile {tranches_indel_fpath} ").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))
        cmd = ("java -jar {gatk_path} -T ApplyRecalibration -R {ref_fpath} -input {raw_indels_vcf_fpath} -mode INDEL --ts_filter_level 99.0 -recalFile {recal_indel_fpath} "
               "-tranchesFile {tranches_indel_fpath} -o {vcf_fpath} ").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

    return vcf_fpath


def do(ref_fpath, sampleIDs, bam_fpaths, scratch_dirpath, output_dirpath, is_human, project_id):
    vcf_fpath = process_files(ref_fpath, sampleIDs, bam_fpaths, scratch_dirpath, output_dirpath, is_human, project_id)
    return vcf_fpath
