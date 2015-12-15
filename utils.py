#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################

import os
import subprocess
import sys
import config


def name_from_fpath(fpath):
    return os.path.splitext(os.path.basename(fpath))[0]


def call_subprocess(args, stdin=None, stdout=None, stderr=None, env=None, debug=False):
    printed_args = args[:]
    if stdin:
        printed_args += ['<', stdin.name]
    if stdout:
        printed_args += ['>>' if stdout.mode in ['a', 'a+'] else '>', stdout.name]
    if stderr:
        printed_args += ['2>>' if stderr.mode in ['a', 'a+'] else '2>', stderr.name]

    if debug:
        print '    Starting:', ' '.join(printed_args)
    return_code = subprocess.call(args, stdin=stdin, stdout=stdout, stderr=stderr, env=env)
    if return_code != 0 and 'ValidateSamFile' not in args:
        print 'ERROR!'
        sys.exit(1)

    return return_code


def assert_file_exists(fpath, message=''):
    if not os.path.isfile(fpath):
        print ("File not found (%s): %s" % (message, fpath))
        sys.exit(2)
    return fpath


def get_path_to_program(program):
    """
    returns the path to an executable or None if it can't be found
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return exe_file
    return None


def check_gatk():
    if not os.path.exists(config.gatk_fpath):
        print 'ERROR! GATK was not found. Please move GenomeAnalysisTK.jar to ' + config.gatk_dirpath + ' and restart the application.'
        sys.exit(1)


def check_external_programs(names):
    for name in names:
        if not get_path_to_program(name):
            print 'ERROR! %s was not found in PATH, please install it manually and/or add to PATH variable' % name
            sys.exit(1)


def check_dbsnp():
    if not os.path.exists(config.dbsnp_fpath):
        print 'Warning: dbSNP was not found. Full pipeline requires dbSNP, thus reduced workflow will be used for now. ' \
              'If you want to use full pipeline, please rename the archive with dbSNP as dbsnp.vcf.gz, move it to ' + \
              config.db_dirpath + ' and restart the application.'
        config.reduced_workflow = True
        return False
    return True

def set_threads_and_mem():
    try:
        import multiprocessing
        config.max_threads = max(1, multiprocessing.cpu_count() / 2)
        config.max_memory = max(2, get_free_memory() / 2)

    except (ImportError, NotImplementedError):
        pass


def get_free_memory():
    from psutil import virtual_memory
    mem = virtual_memory()

    return mem.total / 10 ** 9


def set_reduced_pipeline():
    print 'Warning: full pipeline requires hg19 reference genome. Reduced workflow will be used for now.'
    config.reduced_workflow = True


def prepare_reference(ref_fpath, output_dirpath, silent=False):
    if not silent:
        print 'Preparing reference file...'
    log_fpath = os.path.join(output_dirpath, 'processing_reference.log')
    if not os.path.exists(ref_fpath + '.fai'):
        call_subprocess(['samtools', 'faidx', ref_fpath], stderr=open(log_fpath, 'a'))
    ref_fname, _ = os.path.splitext(ref_fpath)
    if not os.path.exists(ref_fname + '.dict'):
        call_subprocess(['java', '-jar', config.picard_fpath, 'CreateSequenceDictionary', 'R=%s' % ref_fpath,
                               'O=%s' % ref_fname + '.dict'], stderr=open(log_fpath, 'a'))
    get_chr_lengths(ref_fpath)
    if not config.is_run_on_cloud:
        if not config.reduced_workflow or silent:  # check for hg19 reference
            if config.chr_names != config.hg19_chr_names:
                if not silent:
                    set_reduced_pipeline()
                return False
            for chr in config.chr_lengths:
                if chr not in config.hg19_chr_lengths or config.chr_lengths[chr] != config.hg19_chr_lengths[chr]:
                    if not silent:
                        set_reduced_pipeline()
                    return False
    return True


def get_chr_lengths(ref_fpath):
    ref_index_file = open(ref_fpath + '.fai')
    config.chr_lengths = {}
    config.chr_names = []
    for line in ref_index_file.read().split('\n'):
        if line:
            line = line.split()
            config.chr_names.append(line[0])
            config.chr_lengths[line[0]] = int(line[1])
    ref_index_file.close()