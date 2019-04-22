"""
test_basic_plot.py
"""
import sys
import os
import pytest
import shutil

# define location of input files for testing
mydir = os.path.dirname(os.path.abspath(__file__))

# import functions to aid testing
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
from helper import *

from quanformer.basic_plot import *

import matplotlib as mpl
mpl.use("Agg")  # for Mac OS X error of NSInvalidArgumentException on Travis CI

# -----------------------


def test_basic_plot():
    basic_plot(os.path.join(mydir, 'data_tests', 'survey_confs', 'divrefine-220.sdf'),
        tag='QM Psi4 Final Opt. Energy (Har) b3lyp-d3mbj/def2-tzvp',
        style='line',
        take_relative=True,
        har_to_kcal=True)
    os.remove('output_0.png')
    os.remove('output_1.png')

def test_combine_files_plot():
    os.chdir(mydir)
    combine_files_plot(os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch_ene.in'),
        molname='Div_6',
        verbose=True,
        take_relative=True,
        har_to_kcal=True)
    os.remove('combined.dat')
    os.remove('combined.png')


# test manually without pytest
if 0:
    sys.path.insert(0, '/home/limvt/Documents/quanformer/quanformer')
    from basic_plot import *
    test_basic_plot()
