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

def call_subprocess(args, stdin=None, stdout=None, stderr=None, env=None):
    printed_args = args[:]
    if stdin:
        printed_args += ['<', stdin.name]
    if stdout:
        printed_args += ['>>' if stdout.mode in ['a', 'a+'] else '>', stdout.name]
    if stderr:
        printed_args += ['2>>' if stderr.mode in ['a', 'a+'] else '2>', stderr.name]

    return_code = subprocess.call(args, stdin=stdin, stdout=stdout, stderr=stderr, env=env)
    if return_code != 0 and 'ValidateSamFile' not in args:
        print 'ERROR! See log files'
        sys.exit(0)

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

def check_samtools():
    if not os.path.exists(os.path.join(config.samtools_dirpath, 'samtools')):
        samtools_path = get_path_to_program('samtools')
        if not samtools_path:
            print 'Compiling SAMtools (details are in ' + os.path.join(config.samtools_dirpath, 'make.log') + ' and make.err)'
            return_code = call_subprocess(['make', '-C', config.samtools_dirpath], stdout=open(os.path.join(config.samtools_dirpath, 'make.log'), 'w'),
                                          stderr=open(os.path.join(config.samtools_dirpath, 'make.err'), 'w'), )

            if return_code != 0 or not os.path.exists(os.path.join(config.samtools_dirpath, 'samtools')):
                print 'Failed to compile SAMtools (' + config.samtools_dirpath + ')! Try to compile it manually. '
                sys.exit(1)
        else:
            config.samtools_dirpath = ''