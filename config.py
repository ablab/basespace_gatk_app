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
threads = 4
memory = 2
basespace_output_dir = ''

# for BaseSpace
UROOT = "/data/output/appresults/"
DROOT = "/data/input/samples/"
OUT_ROOT = "/data/scratch/output/"

