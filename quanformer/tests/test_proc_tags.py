"""
test_proc_tags.py
"""
# local testing vs. travis testing
try:
    from quanformer.proc_tags import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, '/home/limvt/Documents/off_psi4/quanformer')
    from proc_tags import *

# define location of input files for testing
import os
mydir = os.path.dirname(os.path.abspath(__file__))

# -----------------------

import pytest
from helper import *


def test_get_sd_list():
    # TODO
    pass


def test_set_sd_tags_hess():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'gbi_single.sdf'))
    props = {
        'method': 'test-m',
        'basis': 'test-b',
        'package': 'test-p',
        'time': -1
    }
    set_sd_tags(mol, props, 'hess')
    assert oechem.OEHasSDData(
        mol, "QM test-p Hessian Runtime (sec) test-m/test-b") == True
    assert oechem.OEGetSDData(
        mol, "QM test-p Hessian Runtime (sec) test-m/test-b") == '-1'


def test_set_sd_tags_spe_notfinish():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'gbi_single.sdf'))
    props = {
        'method': 'test-m',
        'basis': 'test-b',
        'package': 'test-p',
        'time': -1
    }
    set_sd_tags(mol, props, 'spe')
    assert oechem.OEHasSDData(
        mol, "QM test-p Single Pt. Runtime (sec) test-m/test-b") == True
    assert oechem.OEGetSDData(
        mol, "QM test-p Single Pt. Runtime (sec) test-m/test-b") == '-1'
    assert oechem.OEHasSDData(mol, "Note on Single Pt. test-m/test-b") == True
    assert oechem.OEGetSDData(
        mol, "Note on Single Pt. test-m/test-b") == "JOB DID NOT FINISH"


def test_set_sd_tags_spe_didfinish():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'gbi_single.sdf'))
    props = {
        'method': 'test-m',
        'basis': 'test-b',
        'package': 'test-p',
        'time': -1,
        'finalEnergy': -2
    }
    set_sd_tags(mol, props, 'spe')
    assert oechem.OEHasSDData(
        mol, "QM test-p Single Pt. Runtime (sec) test-m/test-b") == True
    assert oechem.OEGetSDData(
        mol, "QM test-p Single Pt. Runtime (sec) test-m/test-b") == '-1'
    assert oechem.OEHasSDData(
        mol, "QM test-p Final Single Pt. Energy (Har) test-m/test-b") == True
    assert oechem.OEGetSDData(
        mol, "QM test-p Final Single Pt. Energy (Har) test-m/test-b") == '-2'


def test_set_sd_tags_spe_scs():
    # TODO
    pass


def test_set_sd_tags_opt():
    # TODO
    pass


def test_set_sd_tags_opt_scs():
    # TODO
    pass


def test_delete_tag():
    # TODO
    pass


# test manually without pytest
if 0:
    test_set_sd_tags_hess()
    test_set_sd_tags_spe_notfinish()
    test_set_sd_tags_spe_didfinish()
