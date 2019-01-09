#!/usr/bin/env python
"""
proc_tags.py

This script parses output of Psi4 calculations and writes data in SD tags.

Usage:
- import proc_Tags as pt
- pt.set_sd_tags(args)

By: Victoria T. Lim

"""

import openeye.oechem as oechem
import sys


def get_sd_list(mol, datum, Package='Psi4', Method=None, Basisset=None):
    """
    Get list of specified SD tag for all confs in mol.

    Parameters
    ----------
    mol:        OEChem molecule with all of its conformers
    datum:       string description of property of interest
        options implemented: "QM opt energy" "MM opt energy"
    Package:    software package used for QM calculation. Psi4 or Turbomole.
    Method:     string, for specific properties. e.g. 'mp2'
    Basisset:   string, for specific properties. e.g. '6-31+G(d)'

    Returns
    -------
    sdlist: A 1D N-length list for N conformers with property from SDTag.
    """

    # TODO: dictionary
    if datum == "QM opt energy":
        taglabel = "QM %s Final Opt. Energy (Har) %s/%s" % (Package, Method,
                                                            Basisset)
    if datum == "QM opt energy scs":
        taglabel = "QM %s Final Opt. Energy (Har) SCS-%s/%s" % (
            Package, Method, Basisset)
    if datum == "QM opt energy initial":
        taglabel = "QM %s Initial Opt. Energy (Har) %s/%s" % (Package, Method,
                                                              Basisset)
    if datum == "QM spe":
        taglabel = "QM %s Single Pt. Energy (Har) %s/%s" % (Package, Method,
                                                            Basisset)
    if datum == "QM spe scs":
        taglabel = "QM %s Single Pt. Energy (Har) SCS-%s/%s" % (
            Package, Method, Basisset)
    if datum == "MM opt energy":
        taglabel = "MM Szybki Newton Energy"

    if datum == "original index":
        taglabel = "Original omega conformer number"

    if datum == "opt runtime":
        taglabel = "QM %s Opt. Runtime (sec) %s/%s" % (Package, Method,
                                                       Basisset)
    if datum == "spe runtime":
        taglabel = "QM %s Single Pt. Runtime (sec) %s/%s" % (Package, Method,
                                                             Basisset)
    if datum == "opt step":
        taglabel = "QM %s Opt. Steps %s/%s" % (Package, Method, Basisset)

    try:
        taglabel
    except UnboundLocalError as e:
        sys.exit("Error in input tag of extracting SD data.")

    SDList = []
    for j, conf in enumerate(mol.GetConfs()):
        for x in oechem.OEGetSDDataPairs(conf):
            # Case: opt did not finish --> append nan
            if "note on opt." in x.GetTag().lower(
            ) and "did not finish" in x.GetValue().lower():
                SDList.append('nan')
                break
            # Case: want energy value OR want original index number
            elif taglabel.lower() in x.GetTag().lower():
                SDList.append(x.GetValue())
                break
    return SDList


def set_sd_tags(Conf, Props, calctype):
    """
    For one particular conformer, set all available SD tags based on data
    in Props dictionary.

    Warning
    -------
    If the exact tag already exists, and you want to add a new one then there
    will be duplicate tags with maybe different data. (NOT recommended).
    Then the function to get SDList will only get one or the other;
    I think it just gets the first matching tag.

    TODO: maybe add some kind of checking to prevent duplicate tags added

    Parameters
    ----------
    Conf:       Single conformer from OEChem molecule
    Props:      Dictionary output from ProcessOutput function.
                Should contain the keys: basis, method, numSteps,
                initEnergy, finalEnergy, coords, time, pkg
    calctype: string; one of 'opt','spe','hess' for geometry optimization,
        single point energy calculation, or Hessian calculation

    """

    # get level of theory for setting SD tags
    method = Props['method']
    basisset = Props['basis']
    pkg = Props['package']

    # turn parameters into tag descriptions
    full_method = "{}/{}".format(method, basisset)
    cdict = {'spe': 'Single Pt.', 'opt': 'Opt.', 'hess': 'Hessian'}

    # time info can be set for all cases
    taglabel = "QM {} {} Runtime (sec) {}".format(pkg, cdict[calctype],
                                                  full_method)
    oechem.OEAddSDData(Conf, taglabel, str(Props['time']))

    # hessian has no other info for sd tag
    if calctype == 'hess':
        return

    # check that finalEnergy is there. if not, opt probably did not finish
    # make a note of that in SD tag then quit function
    if not 'finalEnergy' in Props:
        taglabel = "Note on {} {}".format(cdict[calctype], full_method)
        oechem.OEAddSDData(Conf, taglabel, "JOB DID NOT FINISH")
        return

    # Set new SD tag for conformer's final energy
    taglabel = "QM {} Final {} Energy (Har) {}".format(pkg, cdict[calctype],
                                                       full_method)
    oechem.OEAddSDData(Conf, taglabel, str(Props['finalEnergy']))

    # Set new SD tag for final SCS-MP2 energy if method is MP2
    if method.lower() == 'mp2':
        taglabel = "QM {} Final SCS-{} Energy (Har) {}".format(
            pkg, cdict[calctype], full_method)
        oechem.OEAddSDData(Conf, taglabel, str(Props['finalSCSEnergy']))

    # Add COSMO energy with outlying charge correction. Turbomole only!
    if 'ocEnergy' in Props:
        if calctype == 'spe':
            print(
                "Extraction of COSMO OC energy from Turbomole not yet supported for SPE calcns"
            )
        elif calctype == 'opt':
            taglabel = "QM {} Final {} Energy with OC correction (Har) {}".format(
                pkg, cdict[calctype], full_method)
            oechem.OEAddSDData(Conf, taglabel, str(Props['ocEnergy']))

    # spe has no other relevant info for sd tag
    if calctype == 'spe':
        return

    # Set new SD tag for original conformer number if not existing
    # !! Opt2 files should ALREADY have this !! Opt2 index is NOT orig index !!
    taglabel = "Original omega conformer number"
    if not oechem.OEHasSDData(Conf, taglabel):
        # if not working with confs, will have no GetIdx
        try:
            oechem.OEAddSDData(Conf, taglabel, str(Conf.GetIdx() + 1))
        except AttributeError as err:
            pass
    # if tag exists, append new conformer ID after the old one
    else:
        try:
            oldid = oechem.OEGetSDData(Conf, taglabel)
            newid = str(Conf.GetIdx() + 1)
            totid = "{}, {}".format(oldid, newid)
            oechem.OESetSDData(Conf, taglabel, totid)
        except AttributeError as err:
            pass

    # Set new SD tag for numSteps of geom. opt.
    taglabel = "QM {} {} Steps {}".format(pkg, cdict[calctype], full_method)
    oechem.OEAddSDData(Conf, taglabel, str(Props['numSteps']))

    # Set new SD tag for conformer's initial energy
    taglabel = "QM {} Initial {} Energy (Har) {}".format(
        pkg, cdict[calctype], full_method)
    oechem.OEAddSDData(Conf, taglabel, str(Props['initEnergy']))


def delete_tag(mol, tag):
    """
    Delete specified SD tag from all conformers of mol.

    Parameters
    ----------
    mol:        OEChem molecule with all of its conformers
    tag:        exact string label of the data to delete

    """
    for j, conf in enumerate(mol.GetConfs()):
        oechem.OEDeleteSDData(conf, tag)
