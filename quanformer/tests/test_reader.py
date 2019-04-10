"""
test_reader.py
"""
import sys
import os
import pytest

# define location of input files for testing
mydir = os.path.dirname(os.path.abspath(__file__))

from quanformer.reader import *

# -----------------------

def test_read_mols():
    mols = read_mols(os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'))
    assert len(list(mols)) == 2

def test_read_mols_fail():
    with pytest.raises(FileNotFoundError):
        mols = read_mols(os.path.join(mydir, 'data_tests', 'blah.sdf'))
    assert True

#def test_read_text_input():
#    mols = read_mols(os.path.join(mydir, 'data_tests', 'survey_confs', 'stitch.in'))
#    assert wholedict[1]['theory'] == ''
#    assert wholedict[1]['fname'] == ''
#    assert wholedict[1]['tagkey'] == ''
#    assert wholedict[1]['label'] == ''
#    assert ref_index == 0


# test manually without pytest
if 0:
    sys.path.insert(0, '/home/limvt/Documents/quanformer/quanformer')
    from utils import *
    test_convert_extension()
