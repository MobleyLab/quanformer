#!/usr/bin/env python
"""
initialize_confs.py

This script takes a list of smiles (or multi-molecule file) and for each:
 - generates conformers with Omega
 - resolves electrostatic clashes with Cartesian optimization using Szybki
 - performs quick optimization with steepest descent opt with Szybki
 - writes out multi-mol multi-conformer file

To use this script as an imported module:
- import initialize_confs
- initialize_confs.initialize_confs('filename', resolve_clash=True, do_opt=True)

By: Christopher I. Bayly, Victoria T. Lim

"""

import os, sys
import openeye.oechem as oechem
import openeye.oeomega as oeomega
import openeye.oeszybki as oeszybki

### ------------------- Functions -------------------


def generate_confs(mol):
    """
    Generate conformers of molecule from its SMILES string.

    Parameters
    ----------
    mol : OEChem molecule

    Returns
    -------
    molWithConfs : OEChem molecule with omega-generated conformers

    """
    molWithConfs = oechem.OEMol(mol)
    omega = oeomega.OEOmega()
    maxConfs = 0
    omega.SetMaxConfs(maxConfs)
    omega.SetStrictStereo(False)
    omega.SetSampleHydrogens(True)
    omega.SetEnumNitrogen(oeomega.OENitrogenEnumeration_All)
    if not omega(molWithConfs):
        print("omega failed on %s" % mol.GetTitle())
        return None
    else:
        return molWithConfs


def resolve_clashes(mol, clashfile):
    """
    Minimize conformers with severe steric interaction.

    Parameters
    ----------
    mol : single OEChem molecule (single conformer)
    clashfile : string
        name of file to write output

    Returns
    -------
    boolean
        True if completed successfully, False otherwise.

    """

    # set general energy options along with the single-point specification
    spSzybki = oeszybki.OESzybkiOptions()
    spSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    spSzybki.SetSolventModel(oeszybki.OESolventModel_Sheffield)
    spSzybki.SetRunType(oeszybki.OERunType_SinglePoint)
    # generate the szybki MMFF94 engine for single points
    szSP = oeszybki.OESzybki(spSzybki)
    # construct minimiz options from single-points options to get general optns
    optSzybki = oeszybki.OESzybkiOptions(spSzybki)
    # now reset the option for minimization
    optSzybki.SetRunType(oeszybki.OERunType_CartesiansOpt)
    # generate szybki MMFF94 engine for minimization
    szOpt = oeszybki.OESzybki(optSzybki)
    # add strong harmonic restraints to nonHs
    szOpt.SetHarmonicConstraints(10.0)
    # construct a results object to contain the results of a szybki calculation
    szResults = oeszybki.OESzybkiResults()
    # work on a copy of the molecule
    tmpmol = oechem.OEMol(mol)
    if not szSP(tmpmol, szResults):
        print('szybki run failed for %s' % tmpmol.GetTitle())
        return False
    Etotsp = szResults.GetTotalEnergy()
    Evdwsp = szResults.GetEnergyTerm(oeszybki.OEPotentialTerms_MMFFVdW)
    if Evdwsp > 35:
        if not szOpt(tmpmol, szResults):
            print('szybki run failed for %s' % tmpmol.GetTitle())
            return False
        Etot = szResults.GetTotalEnergy()
        Evdw = szResults.GetEnergyTerm(oeszybki.OEPotentialTerms_MMFFVdW)
        wfile = open(clashfile, 'a')
        wfile.write(
            '%s resolved bad clash: initial vdW: %.4f ; '
            'resolved EvdW: %.4f\n' % (tmpmol.GetTitle(), Evdwsp, Evdw))
        wfile.close()
        mol.SetCoords(tmpmol.GetCoords())
    oechem.OESetSDData(mol, oechem.OESDDataPair('MM Szybki Single Point Energy'\
, "%.12f" % szResults.GetTotalEnergy()))
    return True


def quick_opt(mol):
    """
    Fast MM optimization to whittle down number of conformers before QM.
    Default Szybki OEOptType type set to steepest descent (SD) based on
       preliminary comparisons.

    Parameters
    ----------
    mol : single OEChem molecule (single conformer)

    Returns
    -------
    boolean
        True if completed successfully, False otherwise.

    """
    # set general energy options along with the run type specification
    optSzybki = oeszybki.OESzybkiOptions()
    optSzybki.SetForceFieldType(oeszybki.OEForceFieldType_MMFF94S)
    optSzybki.SetSolventModel(oeszybki.OESolventModel_Sheffield)
    optSzybki.SetOptimizerType(oeszybki.OEOptType_SD)
    taglabel = 'MM Szybki SD Energy'

    # generate szybki MMFF94 engine for minimization
    szOpt = oeszybki.OESzybki(optSzybki)
    # construct a results object to contain the results of a szybki calculation
    szResults = oeszybki.OESzybkiResults()
    # work on a copy of the molecule
    tmpmol = oechem.OEMol(mol)
    if not szOpt(tmpmol, szResults):
        print('szybki run failed for %s' % tmpmol.GetTitle())
        return False
    mol.SetCoords(tmpmol.GetCoords())
    oechem.OESetSDData(mol, oechem.OESDDataPair(taglabel, "%.12f" \
        % szResults.GetTotalEnergy()))
    return True


### ------------------- Script -------------------


def initialize_confs(smiles, resolve_clash=True, do_opt=True):
    """
    From a file containing smiles strings, generate omega conformers,
       resolve steric clashes, do a quick MM opt, and write SDF output.

    Parameters
    ----------
    smiles : string
        Name of the input molecule file. Note: this isn't strictly limited
        to SMILES files, can also work with SDF, MOL2, etc. The file can be in
        a different location than the directory of which this function is
        called. The output SDF and txt files will be placed in the directory
        of which the function is called.
    resolve_clash : Boolean
        True to resolve steric clashes by geometry optimization using OESzybki
    do_opt : Boolean
        True to run quick geometry optimization using OESzybki with MMFF94S
        force field and Sheffield solvent model

    """
    base, extension = os.path.splitext(os.path.basename(smiles))
    if extension == '.sdf':
        base = base + '_quanformer'
    sdfout = base + '.sdf'

    ### Read in smiles file.
    ifs = oechem.oemolistream()
    if not ifs.open(smiles):
        oechem.OEThrow.Warning("Unable to open %s for reading" % smiles)

    ### Open output file to write molecules.
    ofs = oechem.oemolostream()
    if os.path.exists(sdfout):
        print(
            "Output .sdf file already exists. Exiting initialize_confs.\n{}\n".format(
                os.path.abspath(sdfout)))
        return
    if not ofs.open(sdfout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % sdfout)

    ### Output files detailing number of resolved clashes
    ###   and original number of conformers before MM opt.
    conffile = open('numConfs.txt', 'a')
    conffile.write("Number of original conformers\n")

    ### For each molecule: label atoms, generate confs, resolve clashes, optimize.
    for smimol in ifs.GetOEMols():
        oechem.OETriposAtomNames(smimol)
        oechem.OEAddExplicitHydrogens(smimol)
        mol = generate_confs(smimol)
        if mol is None:
            continue
        conffile.write("%s\t%s\n" % (mol.GetTitle(), mol.NumConfs()))

        for i, conf in enumerate(mol.GetConfs()):
            print(mol.GetTitle(), i + 1)
            ### Resolve bad clashes.
            if resolve_clash:
                print("Resolving bad clashes...")
                if not resolve_clashes(conf, "numClashes.txt"):
                    print('Resolving bad clashes failed for molecule %s \
conformer %d:' % (mol.GetTitle(), i + 1))
                    continue
            ### MM optimization.
            if do_opt:
                print("Doing a quick MM (SD) optimization...")
                if not quick_opt(conf):
                    print('Quick optimization failed for molecule %s \
conformer %d:' % (mol.GetTitle(), i + 1))
                    continue
        oechem.OEWriteConstMolecule(ofs, mol)

    ### Close files.
    ifs.close()
    ofs.close()
    conffile.close()


if __name__ == "__main__":
    initialize_confs(sys.argv[1], sys.argv[2], sys.argv[3])
