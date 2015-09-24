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
bgzip_fpath = os.path.join(config.LIBS_LOCATION, 'misc', 'bgzip')
tabix_fpath = os.path.join(config.LIBS_LOCATION, 'misc', 'tabix')

chr_names_hg19 = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                  'chr21', 'chr22', 'chrX', 'chrY']
col_names_count_vars = ["CountVariants","CompRod","EvalRod","Sample","nProcessedLoci","nCalledLoci","nRefLoci",
                        "nVariantLoci","variantRate","variantRatePerBp","nSNPs","nMNPs","nInsertions","nDeletions",
                        "nComplex","nSymbolic","nMixed","nNoCalls","nHets","nHomRef","nHomVar","nSingletons","nHomDerived",
                        "heterozygosity","heterozygosityPerBp","hetHomRatio","indelRate","indelRatePerBp","insertionDeletionRatio"]
col_names_titv = ["TiTvVariantEvaluator","CompRod","EvalRod","Sample","nTi","nTv","tiTvRatio","nTiInComp","nTvInComp","TiTvRatioStandard","nTiDerived","nTvDerived","tiTvDerivedRatio"]
report_col_names = ["Number","App Result Name","File"]
first_col = ["App Result"]
snp_col_values = ["nSNPs","nSingletons","hetHomRatio"]
snp_col_names = ["Total SNPs Count","Singletons Count","Het/Hom Ratio"]
indel_col_values = ["nInsertions","nDeletions","indelRatePerBp","insertionDeletionRatio"]
indel_col_names = ["Insertions Count","Deletions Count","Indel Rate Per Bp","Insertion/Deletion Ratio"]
tstv_col_values = ["nTi","nTv","tiTvRatio"]
tstv_col_names = ["Transitions Count","Transversions Count","Ts/Tv Ratio"]

def val_to_str(val):
    if val is None:
        return '-'
    else:
        return str(val)

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
    variants = '-V ' + ' -V '.join(raw_vcf_fpaths)
    gatk_path = gatk_fpath()
    cmd = ("java -cp {gatk_path} org.broadinstitute.gatk.tools.CatVariants -R {ref_fpath} {variants} -assumeSorted "
               "-out {merge_vcf_fpath}").format(**locals())
    utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))
    return merge_vcf_fpath

def process_single_chr(ref_fpath, sampleID, bam_fpath, output_dirpath, log_fpath, chr):
    raw_g_vcf_fpath = os.path.join(output_dirpath, sampleID + chr + '.g.vcf')
    utils.call_subprocess(['java', '-jar', gatk_fpath(), '-T', 'HaplotypeCaller', '-R', ref_fpath, '-L', chr,
                            '-I', bam_fpath, '-stand_emit_conf', '10',
                          '-stand_call_conf', '30', '-ERC', 'GVCF', '-variant_index_type', 'LINEAR',
                           '-variant_index_parameter', '128000', '-o', raw_g_vcf_fpath], stderr=open(log_fpath, 'a'))
    return raw_g_vcf_fpath

def process_files(ref_fpath, sampleIDs, bam_fpaths, scratch_dirpath, output_dirpath, is_human, project_id, sample_files, sample_names):
    log_fpath = os.path.join(output_dirpath, project_id + '.log')
    print 'Calling variants...'
    raw_vcf_fpaths = [process_single_file(ref_fpath, sampleIDs[i], bam_fpaths[i], output_dirpath, scratch_dirpath)
                                               for i in range(len(bam_fpaths))]
    n_jobs = min(len(raw_vcf_fpaths), config.max_gatk_threads)
    raw_vcf_fpaths = Parallel(n_jobs=n_jobs)(delayed(merge_vcfs)(output_dirpath, sampleIDs[i], raw_vcf_fpaths[i], ref_fpath)
                                               for i in range(len(raw_vcf_fpaths)))
    vcf_fpath = os.path.join(output_dirpath, project_id + '.vcf')
    if is_human:
        gatk_path = gatk_fpath()
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
        variants = ['-V %s' % raw_vcf_fpaths[i] for i in range(len(raw_vcf_fpaths))]
        variants = ' '.join(variants)

        print 'Joint genotyping...'
        num_threads = str(config.max_gatk_threads)
        cmd = ("java -jar {gatk_path} -T GenotypeGVCFs -R {ref_fpath} {variants} -nt {num_threads} "
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
    report_vars_fpath = os.path.join(scratch_dirpath, project_id + '.var.txt')
    cmd = ("java -jar {gatk_path} -T VariantEval -R {ref_fpath} -eval {vcf_fpath} "
               "-noST -noEV -EV CountVariants -ST Sample -o {report_vars_fpath}").format(**locals())
    utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))

    report_tstv_fpath = os.path.join(scratch_dirpath, project_id + '.tv.txt')
    cmd = ("java -jar {gatk_path} -T VariantEval -R {ref_fpath} -eval {vcf_fpath} "
               "-noST -noEV -EV TiTvVariantEvaluator -ST Sample -o {report_tstv_fpath}").format(**locals())
    utils.call_subprocess(shlex.split(cmd), stderr=open(log_fpath, 'a'))
    printReport(report_vars_fpath, report_tstv_fpath, sample_names, sampleIDs, sample_files, output_dirpath)

    utils.call_subprocess([bgzip_fpath, vcf_fpath], stderr=open(log_fpath, 'a'))
    utils.call_subprocess([tabix_fpath, '-p', 'vcf', vcf_fpath + '.gz'], stderr=open(log_fpath, 'a'))
    return vcf_fpath

def printReport(report_vars_fpath, report_tstv_fpath, sample_names, sample_IDs, sample_files, output_dirpath):
    all_values = {}
    samples = []
    if not os.path.exists(report_vars_fpath) or not os.path.exists(report_tstv_fpath):
        return
    report_vars = open(report_vars_fpath).read().split('\n')
    for line in report_vars:
        if not line:
            break
        if line[0] == "#":
            continue
        line = line.split()
        if line[3] in sample_names or line[3] in sample_IDs:
            if line[3] not in all_values:
                samples.append(line[3])
                all_values[line[3]] = []
            all_values[line[3]].extend(zip(col_names_count_vars, line))
    report_vars = open(report_tstv_fpath).read().split('\n')
    for line in report_vars:
        if not line:
            break
        if line[0] == "#":
            continue
        line = line.split()
        if line[3] in sample_names or line[3] in sample_IDs:
            if line[3] not in all_values:
                samples.append(line[3])
                all_values[line[3]] = []
            all_values[line[3]].extend(zip(col_names_titv, line))
    report_values = [[i+1, samples[i], sample_files[i]] for i in range(len(sample_names))]
    snp_values = [[all_values[id][row_num][1] for row_num in range(len(all_values[samples[0]])) if all_values[samples[0]][row_num][0] in snp_col_values] for id in samples]
    indel_values = [[all_values[id][row_num][1] for row_num in range(len(all_values[samples[0]])) if all_values[samples[0]][row_num][0] in indel_col_values] for id in samples]
    tstv_values = [[all_values[id][row_num][1] for row_num in range(len(all_values[samples[0]])) if all_values[samples[0]][row_num][0] in tstv_col_values] for id in samples]

    info_fpath = os.path.join(output_dirpath, "File_Info.tsv")
    saveTsv(info_fpath, report_col_names, report_values)
    snp_fpath = os.path.join(output_dirpath, "SNP_Info.tsv")
    saveTsv(snp_fpath, first_col + snp_col_names, snp_values, samples)
    indels_fpath = os.path.join(output_dirpath, "Indel_Info.tsv")
    saveTsv(indels_fpath, first_col + indel_col_names, indel_values, samples)
    tstv_fpath = os.path.join(output_dirpath, "TsTv_Info.tsv")
    saveTsv(tstv_fpath, first_col + tstv_col_names, tstv_values, samples)


def saveTsv(tsv_fpath, col_names, row_values, sample_names=None):
    tsv_file = open(tsv_fpath, 'w')
    print >>tsv_file, '\t'.join(col_names)
    for i, row in enumerate(row_values):
        if sample_names:
            print >>tsv_file, '\t'.join([sample_names[i]] + map(val_to_str, row))
        else:
            print >>tsv_file, '\t'.join(map(val_to_str, row))
    tsv_file.close()

def do(ref_fpath, sampleIDs, bam_fpaths, scratch_dirpath, output_dirpath, is_human, project_id, sample_files, sample_names):
    vcf_fpath = process_files(ref_fpath, sampleIDs, bam_fpaths, scratch_dirpath, output_dirpath, is_human, project_id, sample_files, sample_names)
    return vcf_fpath