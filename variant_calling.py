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

gatk_dirpath = os.path.join(config.LIBS_LOCATION, 'gatk')

def gatk_fpath():
    return os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

def process_single_file(ref_fpath, sampleID, bam_fpath, scratch_dirpath, output_dirpath, is_human):

    vcf_fpath = os.path.join(output_dirpath, sampleID + '.vcf')
    raw_vcf_fpath = os.path.join(scratch_dirpath, sampleID + '.g.vcf')
    log_fpath = os.path.join(output_dirpath, sampleID + '.log')

    utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'HaplotypeCaller', '-R', ref_fpath,
                            '-I', bam_fpath, '-stand_emit_conf', '10',
                          '-stand_call_conf', '30', '-ERC', 'GVCF', '-variant_index_type', 'LINEAR',
                           '-variant_index_parameter', '128000', '-o', raw_vcf_fpath if is_human else vcf_fpath], stderr=open(log_fpath, 'a'))
    if is_human:
        hapmap_fpath = os.path.join(config.DIR_HOME, 'genomes', 'hapmap_3.3.vcf')
        omni_fpath = os.path.join(config.DIR_HOME, 'genomes', '1000G_omni.vcf')
        dbsnp_fpath = os.path.join(config.DIR_HOME, 'genomes', 'dbsnp.vcf')
        mills_fpath = os.path.join(config.DIR_HOME, 'genomes', 'gold_indels.vcf')

        recal_fpath = os.path.join(scratch_dirpath, sampleID + '_SNP.recal')
        tranches_fpath = os.path.join(scratch_dirpath, sampleID + '_SNP.tranches')

        raw_indels_vcf_fpath = os.path.join(scratch_dirpath, sampleID + '_raw_indels.vcf')
        recal_indel_fpath = os.path.join(scratch_dirpath, sampleID + '_INDEL.recal')
        tranches_indel_fpath = os.path.join(scratch_dirpath, sampleID + '_INDEL.tranches')
        # variant filtering
        cmd = ("java -jar {gatk_fpath()} -T VariantRecalibrator -R {ref_fpath} -input {raw_vcf_fpath} -resource:hapmap,known=false,training=true,truth=true,prior=15.0 "
               " {hapmap_fpath} -resource:omni,known=false,training=true,truth=true,prior=12.0 {omni_fpath} "
               "-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {dbsnp_fpath} -an DP -an QD -an FS -an SOR -an MQ-an MQRankSum -an ReadPosRankSum -an InbreedingCoeff"
               "-mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile {recal_fpath} "
               "-tranchesFile {tranches_fpath}").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

        cmd = ("java -jar {gatk_fpath()} -T ApplyRecalibration -R {ref_fpath} -input {raw_vcf_fpath} -mode SNP --ts_filter_level 99.0 "
               "-recalFile  {recal_fpath} -tranchesFile {tranches_fpath} -o {raw_indels_vcf_fpath}").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

        cmd = ("java -jar {gatk_fpath()} -T VariantRecalibrator -R {ref_fpath} -input {raw_indels_vcf_fpath} -resource:mills,known=true,training=true,truth=true,prior=12.0 "
               "{mills_fpath} -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff -mode INDEL -tranche 100.0 -tranche 99.9 -tranche 99.0"
               " -tranche 90.0 --maxGaussians 4 -recalFile {recal_indel_fpath} -tranchesFile {tranches_indel_fpath} ").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))
        cmd = ("java -jar {gatk_fpath()} -T ApplyRecalibration -R {ref_fpath} -input {raw_indels_vcf_fpath} -mode INDEL --ts_filter_level 99.0 -recalFile {recal_indel_fpath} "
               "-tranchesFile {tranches_indel_fpath} -o {vcf_fpath} ").format(**locals())
        utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

    return vcf_fpath


def do(ref_fpath, sampleID, bam_fpath, scratch_dirpath, output_dirpath, is_human):
    vcf_fpath = process_single_file(ref_fpath, sampleID, bam_fpath, scratch_dirpath, output_dirpath, is_human)
    return vcf_fpath