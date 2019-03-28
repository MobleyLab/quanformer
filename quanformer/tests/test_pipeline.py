"""
test_pipeline.py
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

#from quanformer.pipeline import *

# -----------------------

def test_name_manager_pass():
    curr_dir, checked_infile, prefix, ext, no_path_infile = name_manager(
        os.path.join(mydir, 'data_tests', 'two_alkanes.smi'))
    assert prefix == 'two_alkanes'
    assert ext == '.smi'
    assert no_path_infile == 'two_alkanes.smi'

def test_name_manager_fail():
    try:
        curr_dir, checked_infile, prefix, ext, no_path_infile = name_manager(
            os.path.join(mydir, 'data_tests', 'blah.sdf'))
    except FileNotFoundError:
        assert True

def test_setup_conformers():
    setup_conformers(os.path.join(mydir, 'data_tests', 'two_alkanes.smi'))

    # check pre-filtered file
    mols = read_mol(os.path.join(mydir, 'two_alkanes.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 9
    mol = next(mols)
    assert mol.NumConfs() == 3

    # check post-filtered file
    mols = read_mol(os.path.join(mydir, 'two_alkanes-200.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 3
    mol = next(mols)
    assert mol.NumConfs() == 3
    os.remove(os.path.join(mydir, 'two_alkanes.sdf'))
    os.remove(os.path.join(mydir, 'two_alkanes-200.sdf'))
    os.remove(os.path.join(mydir, 'numConfs.txt'))

def test_setup_calculations_true():
    setup_calculations(
        os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'),
        'mp2',
        'def2-sv(p)',
        'spe')
    # check file byte size (this line should be updated if confs_to_psi changes)
    assert os.path.getsize(os.path.join('methane', '1', 'input.dat')) == 370
    shutil.rmtree('methane')

def test_setup_calculations_false():
    try:
        setup_calculations(
            os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'),
            'mp2',
            'def2-sv(p)',
            'blah')
    except ValueError:
        assert True

# test manually without pytest
if 1:
    sys.path.insert(0, '/home/limvt/Documents/quanformer/quanformer')
    from pipeline import *
    test_name_manager_pass()
    test_name_manager_fail()
    test_setup_conformers()
    test_setup_calculations_true()
    test_setup_calculations_false()
