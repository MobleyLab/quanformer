#!/usr/bin/env python
"""
reader.py

Functions to parse or read in files or OEMols.

Version:    Apr 2 2019
By:         Victoria T. Lim

"""

import os
import openeye.oechem as oechem
import collections  # ordered dictionary


def read_mols(infile):
    """
    Open a molecule file and return molecules and conformers as OEMols.

    Parameters
    ----------
    infile : string
        name of input file with molecules

    Returns
    -------
    mols : OEMols

    """

    ifs = oechem.oemolistream()
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifs.open(infile):
        raise FileNotFoundError(f"Unable to open {infile} for reading")
    mols = ifs.GetOEMols()

    return mols


def read_text_input(infile, reffile=None, ref_index=None):
    """
    Read input file into an ordered dictionary.
    http://stackoverflow.com/questions/25924244/creating-2d-dictionary-in-python

    Parameters
    ----------
    TODO

    """
    linecount = 0
    wholedict = collections.OrderedDict()
    with open(infile) as f:
        for line in f:
            if line.startswith('#'):
                continue
            dataline = [x.strip() for x in line.split(',')]
            wholedict[linecount] = {
                'theory': dataline[0],
                'fname': dataline[1],
                'tagkey': dataline[2],
                'label': dataline[3],
            }
            linecount += 1

    return wholedict
