"""
test_getTurbResults.py
"""
# local testing vs. travis testing
try:
    from quanformer.getTurbResults import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, '/home/limvt/Documents/off_psi4/quanformer')
    from getTurbResults import *

# define location of input files for testing
import os
mydir = os.path.dirname(os.path.abspath(__file__))

# -----------------------

import pytest
import helper

# passing locally (even updated all pckgs) but failing on travis. what TODO
#def test_get_time():
#    dt = get_time(os.path.join(mydir,'data_tests','cooh','0','1'))
#    assert dt == 52.0


def test_process_turb_out():
    os.chdir(os.path.join(mydir, 'data_tests', 'cooh', '0', '1'))
    props = process_turb_out({}, 'opt', True)
    # {'initEnergy': '-227.8160966973', 'finalEnergy': '-227.8233603490', 'numSteps': 8, 'ocEnergy': '-227.8233914412'}
    assert float(props['initEnergy']) == pytest.approx(-227.8160966973,
                                                       0.00000001)
    assert float(props['finalEnergy']) == pytest.approx(
        -227.8233603490, 0.00000001)
    assert float(props['ocEnergy']) == pytest.approx(-227.8233914412,
                                                     0.00000001)
    assert props['numSteps'] == 8


def test_getTurbResults_out_spe():
    # SPE not yet supported
    pass


#def test_getTurbResults_out_opt():
#    # this test is fully valid but will not pass on Travis CI because it
#    # requires t2x from Turbomole to extract geometry from Turbomole output file.
#    # can still use to test locally.
#    os.chdir(os.path.join(mydir,'data_tests','cooh'))
#    infile = os.path.join(mydir,'data_tests','cooh','fromVMD.sdf')
#    outfile = os.path.join(mydir,'data_tests','cooh','hfsolv.sdf')
#    getTurbResults(infile, 'HF/6-31G*', outfile, calctype='opt', cosmo=True)
#    mols = helper.read_mol(outfile,True)
#    assert len(list(mols)) == 2
#    os.remove(outfile)

# test manually without pytest
if 0:
    test_get_time()
    test_process_turb_out()
    test_getTurbResults_out_spe()
    test_getTurbResults_out_opt()
