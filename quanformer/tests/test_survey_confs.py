"""
test_survey_confs.py
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

from quanformer.survey_confs import *

import matplotlib as mpl
#mpl.use("Agg")  # for Mac OS X error of NSInvalidArgumentException on Travis CI

# -----------------------

# TODO ADD FILE CHECKS!!

def test_survey_times():
    survey_confs(
        os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch_time.in'),
        analyze_energies=False,
        analyze_times=True,
        ref_index=None,
        plot_enes=False
    )
    os.remove('timeAvgs.txt')
    os.remove('rename_me.png')

def test_survey_times_wrongtag():
    with pytest.raises(NameError):
        survey_confs(
            os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch_wrongtag.in'),
            analyze_energies=False,
            analyze_times=True,
            ref_index=None,
            plot_enes=False
        )
    os.remove('timeAvgs.txt')
    assert True

def test_survey_energies():
    survey_confs(
        os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch_ene.in'),
        analyze_energies=True,
        analyze_times=False,
        ref_index=None,
        plot_enes=False
    )
    os.remove('relene.dat')

def test_survey_energies_reference():
    survey_confs(
        os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch_ene.in'),
        analyze_energies=True,
        analyze_times=False,
        ref_index=0,
        plot_enes=False
    )
    assert os.path.getsize('relene-rmsd.dat') == 951
    os.remove('relene-rmsd.dat')

def test_survey_energies_plot():
    survey_confs(
        os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch_ene.in'),
        analyze_energies=True,
        analyze_times=False,
        ref_index=0,
        plot_enes=True
    )
    os.remove('relene-rmsd.dat')
    os.remove('barchart.png')


# test manually without pytest
if 0:
    sys.path.insert(0, '/home/limvt/Documents/quanformer/quanformer')
    from survey_confs import *
    test_survey_times()
