"""
test_proc_tags.py
"""
import sys
import os
import pytest

# define location of input files for testing
mydir = os.path.dirname(os.path.abspath(__file__))

# import functions to aid testing
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
from helper import *

from quanformer.proc_tags import *

# -----------------------

def test_define_tag():
    tag = define_tag("QM opt energy", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Final Opt. Energy (Har) mp2/def2-SV(P)"
    tag = define_tag("QM opt energy scs", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Final Opt. Energy (Har) SCS-mp2/def2-SV(P)"
    tag = define_tag("QM opt energy initial", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Initial Opt. Energy (Har) mp2/def2-SV(P)"
    tag = define_tag("QM spe", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Single Pt. Energy (Har) mp2/def2-SV(P)"
    tag = define_tag("QM spe scs", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Single Pt. Energy (Har) SCS-mp2/def2-SV(P)"
    tag = define_tag("MM opt energy", None, None, None)
    assert tag == "MM Szybki SD Energy"
    tag = define_tag("original index", None, None, None)
    assert tag == "Original omega conformer number"
    tag = define_tag("opt runtime", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Opt. Runtime (sec) mp2/def2-SV(P)"
    tag = define_tag("spe runtime", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Single Pt. Runtime (sec) mp2/def2-SV(P)"
    tag = define_tag("opt step", "Psi4", "mp2", "def2-SV(P)")
    assert tag == "QM Psi4 Opt. Steps mp2/def2-SV(P)"

def test_get_sd_list():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'), True)
    mol = next(mols)
    data = get_sd_list(mol, "MM opt energy")
    assert len(data) == 9
    assert data[3] == '6.736580988585'

def test_get_sd_list_fail():
    mols = read_mol(os.path.join(mydir, 'data_tests', 'two_alkanes_prefilt.sdf'), True)
    mol = next(mols)
    try:
        data = get_sd_list(mol, "blah")
    except ValueError:
        assert True

def test_set_sd_tags_hess():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
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
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
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
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
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
    mols = read_mol(os.path.join(mydir, 'data_tests', 'carbon-222.sdf'), True)
    mol = next(mols)
    conf = list(mol.GetConfs())[0]
    taglabel = 'MM Szybki SD Energy'
    assert oechem.OEHasSDData(conf, taglabel) == True
    delete_tag(mol, taglabel)
    assert oechem.OEHasSDData(conf, taglabel) == False


# test manually without pytest
if 1:
    test_get_sd_list()
