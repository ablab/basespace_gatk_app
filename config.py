############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
############################################################################

import datetime
import os
import platform
import sys

DIR_HOME = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
LIBS_LOCATION = os.path.join(DIR_HOME, 'libs')
temp_output_dir = os.path.join(DIR_HOME, 'output')
output_dir = os.path.join(DIR_HOME, 'output')

max_threads = 4
max_single_gatk_mem = 2
max_gatk_threads = 4
max_memory = 8
reduced_workflow = True

low_emit = False
stand_call_conf = '30'
stand_emit_conf = '30'
low_call_conf = '25'
low_emit_conf = '10'

short_options = 'R:S:flt:m:'
long_options  = 'reference= samples= full low-emit threads= memory='

# for BaseSpace
RESULTS_DIR = "/data/output/appresults/"
INPUT_DIR_LOCAL = "/data/input_override/"
INPUT_DIR_CLOUD = "/data/input/"
INPUT_DIR = INPUT_DIR_CLOUD
SCRATCH_DIR = "/data/scratch/"

db_dirname = 'db'

gatk_dirpath = os.path.join(LIBS_LOCATION, 'GATK')
gatk_fpath = os.path.join(gatk_dirpath, 'GenomeAnalysisTK.jar')

picard_dirpath = os.path.join(LIBS_LOCATION, 'picard')
picard_fpath = os.path.join(picard_dirpath, 'picard.jar')
picard_tmp_dirpath = os.path.join(SCRATCH_DIR, 'tmp')

dbsnp_fpath = os.path.join(DIR_HOME, db_dirname, 'dbsnp.vcf')
gold_indels_fpath = os.path.join(DIR_HOME, db_dirname, 'gold_indels.vcf')
tg_indels_fpath = os.path.join(DIR_HOME, db_dirname, '1000G_phase1.indels.hg19.vcf')
hapmap_fpath = os.path.join(DIR_HOME, db_dirname, 'hapmap_3.3.vcf')
omni_fpath = os.path.join(DIR_HOME, db_dirname, '1000G_omni.vcf')
mills_fpath = os.path.join(DIR_HOME, db_dirname, 'gold_indels.vcf')

chr_lengths = {}
chr_names = []
hg19_chr_names = ['chrM', 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                  'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                  'chr21', 'chr22', 'chrX', 'chrY']

hg19_chr_lengths = {'chrM': 16571, 'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
                'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
                'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
                'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
                'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}
