"""
test_confs_to_psi.py
"""
# local testing vs. travis testing
try:
    from quanformer.confs_to_psi import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, '/home/limvt/Documents/off_psi4/quanformer')
    from confs_to_psi import *

# define location of input files for testing
import os
mydir = os.path.dirname(os.path.abspath(__file__))

# -----------------------

import pytest
from helper import *


def test_make_hessian():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
    test_string = make_psi_input(mol, mol.GetTitle(), 'mp2', 'aug-cc-pVTZ',
                                 'hess')
    assert "H, wfn = hessian('mp2', return_wfn=True)" in test_string
    assert "wfn.hessian().print_out()" in test_string


def test_make_frozen():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'freeze.sdf'))
    test_string = make_psi_input(mol, mol.GetTitle(), 'mp2', 'aug-cc-pVTZ')
    assert "4 xyz" in test_string
    assert "1 xyz" in test_string
    assert "3 xyz" in test_string
    assert "12 xyz" in test_string


def test_make_dfmp2_dunning():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
    test_string = make_psi_input(mol, mol.GetTitle(), 'mp2', 'aug-cc-pVTZ')
    assert "df_basis_mp2" not in test_string


def test_make_dfmp2_qzvpd():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
    test_string = make_psi_input(mol, mol.GetTitle(), 'mp2', 'def2-qzvpd')
    assert "df_basis_mp2" not in test_string
    return


def test_make_dfmp2_svpp():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'))
    test_string = make_psi_input(mol, mol.GetTitle(), 'mp2', 'def2-sv(p)')
    assert "def2-sv_p_-ri" in test_string
    return


def test_confs_to_psi():
    confs_to_psi(
        os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'), 'mp2',
        'def2-sv(p)')
    # check file byte size (this line should be updated if confs_to_psi changes)
    assert os.path.getsize(os.path.join('methane', '1', 'input.dat')) == 358
    shutil.rmtree('methane')
    return


def test_confs_to_psi_json():
    confs_to_psi(
        os.path.join(mydir, 'data_tests', 'methane_c2p.sdf'),
        'mp2',
        'def2-sv(p)',
        calctype='spe',
        via_json=True)
    # check file byte size (this line should be updated if confs_to_psi changes)
    assert os.path.getsize(os.path.join('methane', '1', 'input.py')) == 1032
    shutil.rmtree('methane')
    return


# test manually without pytest
if 0:
    test_confs_to_psi()
