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


def define_tag(datum, package, method, basisset):
    # CASE SENSITIVE
    # this function is case-sensitive though get_sd_list should remove sensitivity

    tagdict = {
        "QM opt energy": "QM {} Final Opt. Energy (Har) {}/{}".format(
            package, method, basisset),
        "QM opt energy scs": "QM {} Final Opt. Energy (Har) SCS-{}/{}".format(
            package, method, basisset),
        "QM opt energy initial": "QM {} Initial Opt. Energy (Har) {}/{}".format(
            package, method, basisset),
        "QM spe": "QM {} Final Single Pt. Energy (Har) {}/{}".format(
            package, method, basisset),
        "QM spe scs": "QM {} Single Pt. Energy (Har) SCS-{}/{}".format(
            package, method, basisset),
        "MM opt energy": "MM Szybki SD Energy",
        "original index": "Original omega conformer number",
        "opt runtime": "QM {} Opt. Runtime (sec) {}/{}".format(
            package, method, basisset),
        "spe runtime": "QM {} Single Pt. Runtime (sec) {}/{}".format(
            package, method, basisset),
        "opt step": "QM {} Opt. Steps {}/{}".format(package, method, basisset),
    }
    try:
        taglabel = tagdict[datum]
    except KeyError as exc:
        raise NameError("Error in defining string to extract SD data for "
            f"[ {datum} ]. Please verify that the key for define_tag is one of: ",
            tagdict.keys()) from exc

    return taglabel

def get_sd_list(mol, datum, package='Psi4', method=None, basisset=None, taglabel=None):
    """
    Get list of specified SD tag for all confs in mol.

    Parameters
    ----------
    mol:        OEChem molecule with all of its conformers
    datum:       string description of property of interest
        options implemented: "QM opt energy" "MM opt energy"
    package:    software package used for QM calculation. Psi4 or Turbomole.
    method:     string, for specific properties. e.g. 'mp2'
    basisset:   string, for specific properties. e.g. '6-31+G(d)'
    taglabel : string
        exact tag string from which to extract SD data

    Returns
    -------
    sdlist: A 1D N-length list for N conformers with property from SDTag.
    """

    if taglabel is None:
        taglabel = define_tag(datum, package, method, basisset)

    sd_list = []
    for j, conf in enumerate(mol.GetConfs()):
        for x in oechem.OEGetSDDataPairs(conf):

            # Case: opt did not finish --> append nan
            if "note on opt." in x.GetTag().lower(
            ) and "did not finish" in x.GetValue().lower():
                sd_list.append('nan')
                break

            # Case: want energy value OR want original index number
            elif taglabel.lower() in x.GetTag().lower():
                sd_list.append(x.GetValue())
                break

    return sd_list


def set_sd_tags(Conf, Props, calctype):
    """
    For one particular conformer, set all available SD tags based on data
    in Props dictionary.

    Warning
    -------
    If the exact tag already exists, and you want to add a new one then there
    will be duplicate tags with maybe different data. (NOT recommended).
    Then the function to get sd_list will only get one or the other;
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

    Note: Multi-conformer molecule must be specified
    else will get AttributeError:
    'OEGraphMol' object has no attribute 'GetConfs'.

    Parameters
    ----------
    mol : multi-conformer OEChem molecule
    tag : string
        exact label of the data to delete

    """
    for j, conf in enumerate(mol.GetConfs()):
        oechem.OEDeleteSDData(conf, tag)
