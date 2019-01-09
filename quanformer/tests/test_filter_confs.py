"""
test_filter_confs.py
"""
# local testing vs. travis testing
try:
    from quanformer.filter_confs import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, '/home/limvt/Documents/off_psi4/quanformer')
    from filter_confs import *

# define location of input files for testing
import os
mydir = os.path.dirname(os.path.abspath(__file__))

# -----------------------

import pytest
from helper import *


def test_identify_minima():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'gbi.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 36
    # use same params defined in filter_confs.py script
    assert identify_minima(mol, 'MM Szybki SD Energy', 5.E-4, 0.2) is True
    assert mol.NumConfs() == 5


def test_filter_confs():
    filter_confs(
        os.path.join(mydir, 'data_tests', 'gbi.sdf'), 'MM Szybki SD Energy',
        'output.sdf')
    mols = read_mol(os.path.join(os.getcwd(), 'output.sdf'), True)
    mol = next(mols)
    assert mol.NumConfs() == 5
    os.remove(os.path.join(os.getcwd(), 'output.sdf'))
    os.remove(os.path.join(os.getcwd(), 'numConfs.txt'))


# test manually without pytest
if 0:
    test_identify_minima()
    test_filter_confs()
