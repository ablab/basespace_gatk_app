#!/usr/bin/env python

############################################################################
# Copyright (c) 2015 Saint Petersburg State University
# All Rights Reserved
# See file LICENSE for details.
#############################################################################

import os
from os.path import relpath
import subprocess


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
    if return_code != 0:
        print 'The tool returned non-zero.'

    return return_code