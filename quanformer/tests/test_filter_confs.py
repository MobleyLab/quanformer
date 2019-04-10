"""
test_filter_confs.py
"""
import sys
import os
import pytest

# define location of input files for testing
mydir = os.path.dirname(os.path.abspath(__file__))

# import functions to aid testing
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
from helper import *

from quanformer.filter_confs import *

# -----------------------

def test_identify_minima():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 9
    # use same params defined in filter_confs.py script
    assert identify_minima(mol, 'MM Szybki SD Energy', 5.E-4, 0.2) is True
    assert mol.NumConfs() == 3

def test_identify_minima_fail():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 9
    assert identify_minima(mol, 'blah', 5.E-4, 0.2) is False

def test_identify_minima_one():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 1
    assert identify_minima(mol, 'MM Szybki SD Energy', 5.E-4, 0.2) is True
    assert mol.NumConfs() == 1

def test_identify_minima_one_fail():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 1
    assert identify_minima(mol, 'blah', 5.E-4, 0.2) is False

def test_filter_confs():
    filter_confs(
        os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'), 'MM Szybki SD Energy',
        os.path.join(mydir, 'data_tests', 'output.sdf'))
    mols = read_mol(os.path.join(mydir, 'data_tests', 'output.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 3
    os.remove(os.path.join(mydir, 'data_tests', 'output.sdf'))
    os.remove(os.path.join(os.getcwd(), 'numConfs.txt')) # don't use mydir here

def test_filter_confs_exists():
    with pytest.raises(FileExistsError):
        filter_confs(
            os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'), 'MM Szybki SD Energy',
            os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'))
    assert True

def test_filter_confs_unfinished():
    filter_confs(
        os.path.join(mydir, 'data_tests', 'gbi_prefilt.sdf'), 'QM Psi4 Final Opt. Energy (Har) mp2/def2-SV(P)',
        os.path.join(mydir, 'data_tests', 'output.sdf'))
    mols = read_mol(os.path.join(mydir, 'data_tests', 'output.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 1
    os.remove(os.path.join(mydir, 'data_tests', 'output.sdf'))
    os.remove(os.path.join(os.getcwd(), 'numConfs.txt')) # don't use mydir here

# test manually without pytest
if 0:
    test_identify_minima()
    test_filter_confs()
