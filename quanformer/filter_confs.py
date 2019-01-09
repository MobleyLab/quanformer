#!/usr/bin/env python

## By: Victoria T. Lim

## This script takes an SDF file and, for each molecule:
##  - conformers are compared by energy
##  - conformers are compared by RMSD
## in order to roughly filter out duplicate minima and keep unique ones.
## Filtered conformers for all molecules are written out in SDF file.

## Import and call filter_confs.filter_confs(rmsdfile, tag, rmsdout)

import re
import os, sys, glob
import openeye.oechem as oechem

### ------------------- Functions -------------------


def identify_minima(mol, tag, ThresholdE, ThresholdRMSD):
    """
    For a molecule's set of conformers computed with some level of theory,
        whittle down unique conformers based on energy and RMSD.

    Parameters
    ----------
    mol           OEChem molecule with all of its conformers
    tag           string name of the SD tag in this molecule
    ThresholdE    float value for abs(E1-E2), below which 2 confs are "same"
        Units are hartrees (default output units of Psi4)
    ThresholdR    float value for RMSD, below which 2 confs are "same"
        Units are in Angstrom (Psi4 default)

    Returns
    -------
    boolean True if successful filter + delete. False if there's only
        one conf and it didn't optimize, or something else funky.

    """
    # Parameters for OpenEye RMSD calculation
    automorph = True
    heavyOnly = False
    overlay = True

    # declare variables for conformers to delete
    confsToDel = set()
    delCount = 0

    # check if SD tag exists for the case of single conformer
    if mol.NumConfs() == 1:
        testmol = mol.GetConfs().next()
        for x in oechem.OEGetSDDataPairs(mol):
            if tag.lower() in x.GetTag().lower():
                return True
            else:
                return False

    # Loop over conformers twice (NxN diagonal comparison of RMSDs)
    for confRef in mol.GetConfs():
        print(" ~ Reference: %s conformer %d" % (mol.GetTitle(),
                                                 confRef.GetIdx() + 1))

        # get real tag (correct for capitalization)
        for x in oechem.OEGetSDDataPairs(confRef):
            if tag.lower() in x.GetTag().lower():
                taglabel = x.GetTag()

        # check that there was a successfully found tag
        try:
            taglabel
        except NameError:
            print("Unable to filter bc missing SD data based on tag: {}".
                  format(tag))
            return False

        # delete cases that don't have energy (opt not converged; or other)
        if not oechem.OEHasSDData(confRef, taglabel):
            confsToDel.add(confRef.GetIdx())
            delCount += 1
            continue
        refE = float(oechem.OEGetSDData(confRef, taglabel))

        for confTest in mol.GetConfs():
            # upper right triangle comparison
            if confTest.GetIdx() <= confRef.GetIdx():
                continue
            # skip cases already set for removal
            if confTest.GetIdx() in confsToDel:
                continue
            # delete cases that don't have energy
            if not oechem.OEHasSDData(confTest, taglabel):
                confsToDel.add(confTest.GetIdx())
                continue

            testE = float(oechem.OEGetSDData(confTest, taglabel))
            # if MM (not Psi4) energies, convert absERel to Hartrees
            if 'mm' in taglabel.lower():
                absERel = abs(refE - testE) / 627.5095
            else:
                absERel = abs(refE - testE)
            # if energies are diff enough --> confs are diff --> keep & skip ahead
            if absERel > ThresholdE:
                continue
            # if energies are similar, see if they are diff by RMSD
            rmsd = oechem.OERMSD(confRef, confTest, automorph, heavyOnly,
                                 overlay)
            # if measured_RMSD < threshold_RMSD --> confs are same --> delete
            if rmsd < ThresholdRMSD:
                confsToDel.add(confTest.GetIdx())

    # for the same molecule, delete tagged conformers
    print("%s original number of conformers: %d" % (mol.GetTitle(),
                                                    mol.NumConfs()))
    if delCount == mol.NumConfs():
        # all conformers in this mol has been tagged for deletion
        return False
    for conf in mol.GetConfs():
        if conf.GetIdx() in confsToDel:
            print('Removing %s conformer index %d' % (mol.GetTitle(),
                                                      conf.GetIdx()))
            if not mol.DeleteConf(conf):
                oechem.OEThrow.Fatal("Unable to delete %s GetIdx() %d" \
                                  % (mol.GetTitle(), conf.GetIdx()))
    return True


### ------------------- Script -------------------


def filter_confs(rmsdfile, tag, rmsdout):
    """
    Read in OEMols (and each of their conformers) in 'rmsdfile'.
    For each molecule:
        rough filter conformers based on energy differences specified by 'tag',
        fine filter conformers based on RMSD values.

    Parameters
    ----------
    rmsdfile : str
        Name of SDF file with conformers to be filtered
    tag : str
        SD tag name with the energy value to roughly screen conformers before RMSD
        Screening works by removing conformers of very similar energies, where
        "similar" is defined by thresE parameter. Examples:
        - "QM Psi4 Final Opt. Energy (Har) mp2/def-sv(p)"
        - "QM Psi4 Final Single Pt. Energy (Har) mp2/def-sv(p)"
    rmsdout : str
        Name of the output file with filtered conformers

    """
    # Parameters for distinguishing cutoff of conformer similarity
    thresE = 5.E-4  # declare confs diff & skip RMSD comparison above this threshold
    thresRMSD = 0.2  # above this threshold (Angstrom), confs are "diff" minima

    wdir, fname = os.path.split(rmsdfile)
    numConfsF = open(os.path.join(os.getcwd(), "numConfs.txt"), 'a')
    numConfsF.write("\n{}\n".format(tag))

    # Open file to be processed.
    rmsd_ifs = oechem.oemolistream()
    if not rmsd_ifs.open(rmsdfile):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % rmsdfile)
    rmsd_ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    rmsd_molecules = rmsd_ifs.GetOEMols()

    # Open outstream file.
    rmsd_ofs = oechem.oemolostream()
    if os.path.exists(rmsdout):
        print("%s output file already exists in %s. Skip filtering.\n" %
              (rmsdout, os.getcwd()))
        return
    if not rmsd_ofs.open(rmsdout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % rmsdout)

    # Identify minima and write output file.
    for mol in rmsd_molecules:
        if identify_minima(mol, tag, thresE, thresRMSD):
            numConfsF.write("%s\t%s\n" % (mol.GetTitle(), mol.NumConfs()))
            oechem.OEWriteConstMolecule(rmsd_ofs, mol)
        else:
            numConfsF.write("%s\t0\n" % (mol.GetTitle()))
    rmsd_ifs.close()
    numConfsF.close()
    rmsd_ofs.close()

    print("Done filtering %s to %s.\n" % (fname, rmsdout))
