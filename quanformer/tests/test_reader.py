"""
test_reader.py
"""
import sys
import os
import pytest
import openeye.oechem as oechem

# define location of input files for testing
mydir = os.path.dirname(os.path.abspath(__file__))

from quanformer.reader import *

# -----------------------

def test_read_mols():
    mols = read_mols(os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'))
    list_from_gen = list(mols)
    assert len(list_from_gen) == 2
    #mol = list_from_gen[0]
    #conf = list(mol.GetConfs())[0]
    #assert oechem.OEHasSDData(conf, "MM Szybki SD Energy") == True

def test_read_mols_slice():
    mlist = read_mols(
        os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'),
        mol_slice=[0,1,1])
    assert len(mlist) == 1
    assert mlist[0].GetTitle() == 'AlkEthOH_c312'
    assert mlist[0].NumConfs() == 9
    conf = list(mlist[0].GetConfs())[0]
    assert oechem.OEHasSDData(conf, "MM Szybki SD Energy") == True

    #AlkEthOH_c1178, Div_2, Div_6, Div_9, Div_3b, Div_7b, Div_8b, AlkEthOH_r187
    mlist = read_mols(
        os.path.join(mydir, 'data_tests', 'eight_mols.sdf'),
        mol_slice=[2,8,2])
    assert len(mlist) == 3
    assert mlist[0].GetTitle() == 'Div_6'
    assert mlist[1].GetTitle() == 'Div_3b'
    assert mlist[2].GetTitle() == 'Div_8b'

def test_read_mols_slice_invalid_incr():
    with pytest.raises(ValueError):
        mlist = read_mols(
            os.path.join(mydir, 'data_tests', 'eight_mols.sdf'),
            mol_slice=[2,8,-2])
    assert True

def test_read_mols_slice_invalid_start():
    with pytest.raises(ValueError):
        mlist = read_mols(
            os.path.join(mydir, 'data_tests', 'eight_mols.sdf'),
            mol_slice=[8,2,2])
    assert True

def test_read_mols_slice_invalid_len():
    with pytest.raises(ValueError):
        mlist = read_mols(
            os.path.join(mydir, 'data_tests', 'eight_mols.sdf'),
            mol_slice=[8])
    assert True

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
