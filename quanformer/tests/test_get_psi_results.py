"""
test_get_psi_results.py
"""
# local testing vs. travis testing
try:
    from quanformer.get_psi_results import *
except ModuleNotFoundError:
    import sys
    sys.path.insert(0, '/home/limvt/Documents/off_psi4/quanformer')
    from get_psi_results import *

# define location of input files for testing
import os
mydir = os.path.dirname(os.path.abspath(__file__))

# -----------------------

import pytest
from helper import *


def test_initiate_dict():
    d = initiate_dict()
    assert d['package'] == 'Psi4'
    assert d['missing'] == False


def test_get_conf_data():
    # TODO
    pass


def test_set_conf_data():
    # TODO
    pass


def test_check_title():
    mol = read_mol(os.path.join(mydir, 'data_tests', 'methane_title-1.0.sdf'))
    mol = check_title(
        mol, os.path.join(mydir, 'data_tests', 'methane_title-1.0.sdf'))
    assert mol.GetTitle() == 'methane_title10'


def test_get_psi_time():
    time = get_psi_time(os.path.join(mydir, 'data_tests', 'timer.dat'))
    assert time == 847.00


def test_process_psi_out_spe():
    spe_dict = process_psi_out(
        os.path.join(mydir, 'data_tests', 'output_spe.dat'), {}, 'spe')
    # {'basis': 'def2-tzvp', 'method': 'b3lyp-d3mbj', 'finalEnergy': -466.37497616456733}
    # actual energy: -463.3016772141738784
    # scs final: -465.1942227292489633
    assert spe_dict['basis'] == 'cc-pVTZ'
    assert spe_dict['method'] == 'mp2'
    assert spe_dict['finalEnergy'] == pytest.approx(-463.3016772141738784,
                                                    0.000000000001)
    assert spe_dict['finalSCSEnergy'] == pytest.approx(-465.1942227292489633,
                                                       0.000000000001)


def test_process_psi_out_hess():
    hess_dict = process_psi_out(
        os.path.join(mydir, 'data_tests', 'output_hess.dat'), {}, 'hess')
    assert hess_dict['basis'] == 'def2-sv(p)'
    assert hess_dict['method'] == 'mp2'
    hess = hess_dict['hessian']
    assert hess.shape == (36, 36)
    assert hess[1, 1] == 0.77147127082666


def test_process_psi_out_opt():
    opt_dict = process_psi_out(
        os.path.join(mydir, 'data_tests', 'output_opt.dat'), {}, 'opt')
    # initial energy: -582.148922080397
    # final energy: -582.1568394053036
    # scs final: -582.0904465352028865
    assert opt_dict['basis'] == 'def2-SV(P)'
    assert opt_dict['method'] == 'mp2'
    assert opt_dict['numSteps'] == '8'
    assert opt_dict['initEnergy'] == pytest.approx(-582.148922080397,
                                                   0.000000000001)
    assert opt_dict['finalEnergy'] == pytest.approx(-582.1568394053036,
                                                    0.000000000001)
    assert opt_dict['finalSCSEnergy'] == pytest.approx(-582.0904465352028865,
                                                       0.000000000001)
    assert opt_dict['coords'][0] == -0.0238533448
    assert len(opt_dict['coords']) == 69


def test_process_psi_out_two():
    # TODO what happens if passed in psi4 output with opt-->hess, or opt-->spe
    pass


def test_getPsiOne():
    infile = os.path.join(mydir, 'data_tests', 'gbi_single.sdf')
    outfile = os.path.join(mydir, 'data_tests', 'gbi_single-210.sdf')
    psiout = os.path.join(mydir, 'data_tests', 'GBI', '1', 'output.dat')
    timeout = os.path.join(mydir, 'data_tests', 'GBI', '1', 'timer.dat')
    mol, props = getPsiOne(infile, outfile, 'opt', psiout, timeout)
    assert oechem.OEHasSDData(
        mol, "QM Psi4 Opt. Runtime (sec) mp2/def2-SV(P)") == True
    assert oechem.OEGetSDData(
        mol, "QM Psi4 Opt. Runtime (sec) mp2/def2-SV(P)") == '847.0'
    assert oechem.OEHasSDData(
        mol, "QM Psi4 Final Opt. Energy (Har) mp2/def2-SV(P)") == True
    assert oechem.OEGetSDData(
        mol, "QM Psi4 Final Opt. Energy (Har) mp2/def2-SV(P)"
    ) == '-582.1568394053036'
    assert oechem.OEHasSDData(
        mol, "QM Psi4 Initial Opt. Energy (Har) mp2/def2-SV(P)") == True
    assert oechem.OEGetSDData(
        mol, "QM Psi4 Initial Opt. Energy (Har) mp2/def2-SV(P)"
    ) == '-582.148922080397'
    assert oechem.OEHasSDData(mol, "QM Psi4 Opt. Steps mp2/def2-SV(P)") == True
    assert oechem.OEGetSDData(mol, "QM Psi4 Opt. Steps mp2/def2-SV(P)") == '8'
    os.remove(outfile)


# TODO
#def test_get_psi_results_none(capsys):
#    infile = os.path.join(mydir,'data_tests','gbi-200.sdf')
#    outfile = os.path.join(mydir,'data_tests','gbi-210.sdf')
#    with pytest.raises(SystemExit):
#        m, b = get_psi_results(infile, outfile, calctype='blah', psiout="output.dat", timeout="timer.dat")
#    out, err = capsys.readouterr()
#    assert err == "Specify a valid calculation type."
#    print(out, err)


def test_get_psi_results_spe():
    # TODO
    pass


def test_get_psi_results_hess():
    infile = os.path.join(mydir, 'data_tests', 'carbon-222.sdf')
    outfile = os.path.join(mydir, 'data_tests', 'carbon_hess-222.sdf')
    m, b = get_psi_results(infile, outfile, 'hess', "output.dat", "timer.dat")
    # check the method and basis set returned
    assert m == 'mp2'
    assert b == 'def2-tzvp'

    # check tag of Hessian calculation time for one conformer
    mols = read_mol(outfile, True)
    conf = list(next(mols).GetConfs())[0]  # only conf of the s1 mol
    assert oechem.OEGetSDData(
        conf, 'QM Psi4 Hessian Runtime (sec) mp2/def2-tzvp') == '253.0'

    # check a few values of the Hessian matrix for each conformer
    hpickle = os.path.join(mydir, 'data_tests', 'carbon_hess-222.hess.pickle')
    hdict = pickle.load(open(hpickle, 'rb'))
    # mol s1, conf 1, row 10, column 23
    assert hdict['s1'][1][9][22] == 0.00065859359798
    # mol t1, conf 1, row 3, column 9
    assert hdict['t1'][1][2][8] == -0.42670095298438

    # clean up
    os.remove(outfile)
    os.remove(hpickle)


def test_get_psi_results_opt():
    infile = os.path.join(mydir, 'data_tests', 'gbi-200.sdf')
    outfile = os.path.join(mydir, 'data_tests', 'gbi-210.sdf')
    m, b = get_psi_results(infile, outfile, 'opt', "output.dat", "timer.dat")
    assert m == 'mp2'
    assert b == 'def2-SV(P)'

    # read the single mol (from generator) and get its three conformers
    # confs 2 and 4 failed opt so we should have data from confs 1 3 5
    mol = read_mol(outfile, True)
    confs = list(next(mol).GetConfs())
    assert len(confs) == 3

    # get one conf and check its tag data
    conf = confs[1]
    # handy way to get all tags and data (.GetTag(), .GetValue())
    #print([x.GetValue() for x in list(oechem.OEGetSDDataPairs(conf))])
    assert oechem.OEGetSDData(
        conf, 'QM Psi4 Opt. Runtime (sec) mp2/def2-SV(P)') == '2635.0'
    assert oechem.OEGetSDData(
        conf, 'QM Psi4 Final Opt. Energy (Har) mp2/def2-SV(P)'
    ) == '-582.1570265488717'
    assert oechem.OEGetSDData(conf,
                              'QM Psi4 Opt. Steps mp2/def2-SV(P)') == '21'
    assert oechem.OEGetSDData(
        conf, 'QM Psi4 Initial Opt. Energy (Har) mp2/def2-SV(P)'
    ) == '-582.146838702714'
    os.remove(outfile)


# test manually without pytest
if 0:
    #test_getPsiOne()
    #test_get_psi_results_none()
    test_get_psi_results_hess()
