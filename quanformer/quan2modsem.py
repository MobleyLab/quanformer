#!/usr/bin/env python
"""
quan2modsem.py

Purpose:    This script interfaces the output of Psi4 Hessian calculations
            (organized by Quanformer) to the modified Seminario method for
            obtaining force constants.

            The modified Seminario method is applied from this code:
            https://github.com/aa840/ModSeminario_Py

Reference:  Allen et al., 10.1021/acs.jctc.7b00785

Usage:      [1] Obtain modified Seminario Python code (see reference below)
            [2] Alter modified_Seminario_method.py:
                [2.1] Comment out all lines with "input_data_processing". (should be 2)
                [2.2] Replace def line with the following:
                      def modified_Seminario_method(bond_list, angle_list, coords, N, hessian, atom_names, inputfilefolder, outputfilefolder, vibrational_scaling):
            [3] Comment out the two lines at the bottom: import sys, modified_Seminario_method(...)
            [4] python quan2modsem.py -i file.sdf -h file.hess.pickle

By:         Victoria T. Lim

Version:    Oct 22 2018

Assumptions:
 - Input SDF file should be from setting up Hessian calculations. All
   molecules/conformers should be represented here as in directory layout to
   ensure correct data directory when iterating over confs.

Notes:
 - This is the most straightforward way to change the modsem code. Other
   approaches would be (1) write my own input_data_processing, but this would
   entail repeatedly opening and closing the SDF file or changing parameters of
   both modsem scripts; (2) add the input data processing code directly into
   the main modsem script, which would be unwieldy.

"""

import os, sys
import numpy as np
import itertools
import pickle
import openeye.oechem as oechem

sys.path.insert(
    0,
    '/beegfs/DATA/mobley/limvt/openforcefield/hessian/modsem/Python_Modified_Seminario_Method'
)
import modified_Seminario_method_vtl2


def hbb_to_kaa(hessian):
    """
    Unit conversions on input Hessian matrix from (Hartrees/Bohr/Bohr)
    (kcal/mol/Angstrom/Angstrom).

    """
    hessian = (hessian * 627.509474) / (0.529177**2)
    return hessian


def prep_hess(mol, hessian):

    # extract bonds as list of indices
    bond_list = []
    for bond in mol.GetBonds():
        bond_list.append([bond.GetBgnIdx(), bond.GetEndIdx()])

    # obtain angles as list of indices
    angle_list = []
    atom_names = []
    for atom in mol.GetAtoms():
        ## get atom names
        #atom_names.append(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()))

        nbor_list = [n.GetIdx() for n in list(atom.GetAtoms())]

        # get all combinations of two atoms which are the outer angle atoms
        for outers in itertools.combinations(nbor_list, 2):
            this_angle = list(outers)
            this_angle.insert(1, atom.GetIdx())
            angle_list.append(this_angle)

    # format coordinates into Nx3 array
    coords = list(mol.GetCoords().values())
    coords = np.array(coords)

    # define variable for number of atoms
    N = mol.NumAtoms()

    # convert hessian units
    hessian = hbb_to_kaa(hessian)

    # the atom names are being returned as a list of string-formatted integers
    # string for modified Seminario code; ints for use with CCB ChemPer
    atom_names = [str(i) for i in range(0, N)]
    return bond_list, angle_list, coords, N, hessian, atom_names


def quan2modsem(infile, pfile):

    hdir, fname = os.path.split(infile)
    wdir = os.getcwd()

    # read in sdf file and distinguish each molecule's conformers
    ifs = oechem.oemolistream()
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifs.open(infile):
        sys.exit("Unable to open %s for reading" % infile)
    molecules = ifs.GetOEMols()

    # open quanformer-generated pickle file with dictionary of hessians
    hdict = pickle.load(open(pfile, 'rb'))

    for mol in molecules:
        print("===== %s =====" % (mol.GetTitle()))
        for j, conf in enumerate(mol.GetConfs()):

            # set file locations; dir for modsem needs / at end of string
            datadir = os.path.join(hdir, "%s/%s/" % (mol.GetTitle(), j + 1))

            # extract hessian from the quanformer-generated dictionary (get_psi_results)
            hessian = hdict[mol.GetTitle()][j + 1]

            # run modsem
            bond_list, angle_list, coords, N, hessian, atom_names = prep_hess(
                mol, hessian)
            modified_Seminario_method_vtl2.modified_Seminario_method(
                bond_list,
                angle_list,
                coords,
                N,
                hessian,
                atom_names,
                datadir,
                datadir,
                vibrational_scaling=1)

            # check to make sure files were generated, and note the ones that didn't work
            # TODO

    ifs.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", required=True,
        help="Input SDF file with all conformers of Hessian data")
    parser.add_argument("-p", "--pfile", required=True,
        help="Associated pickle file with dictionary of extracted Hessian matrices")

    args = parser.parse_args()
    opt = vars(args)

    quan2modsem(args.infile, args.pfile)
