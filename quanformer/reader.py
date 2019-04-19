#!/usr/bin/env python
"""
reader.py

Functions to parse or read in files or OEMols.

Version:    Apr 2 2019
By:         Victoria T. Lim

"""

import openeye.oechem as oechem
import collections  # ordered dictionary
import copy


def read_mols(infile, mol_slice=None):
    """
    Open a molecule file and return molecules and conformers as OEMols.

    Parameters
    ----------
    infile : string
        name of input file with molecules
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]

    Returns
    -------
    mols : OEMols

    """
    ifs = oechem.oemolistream()
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifs.open(infile):
        raise FileNotFoundError(f"Unable to open {infile} for reading")
    mols = ifs.GetOEMols()

    if mol_slice is not None:
        if len(mol_slice) != 3 or mol_slice[0] >= mol_slice[1] or mol_slice[2] <= 0:
            raise ValueError("Check input to mol_slice. Should have len 3, "
                "start value < stop value, step >= 1.")
        # TODO more efficient. can't itertools bc lost mol info (name, SD) after next()
        # adding copy/deepcopy doesnt work on generator objects
        # also doesn't work to convert generator to list then slice list
        #mols = itertools.islice(mols, mol_slice[0], mol_slice[1], mol_slice[2])
        #mlist = mlist[mol_slice[0]:mol_slice[1]:mol_slice[2]]

        def incrementer(count, mols, step):
            if step == 1:
                count += 1
                return count
            # use step-1 because for loop already increments once
            for j in range(step-1):
                count += 1
                next(mols)
            return count

        mlist = []
        count = 0
        for i, m in enumerate(mols):

            if count >= mol_slice[1]:
                return mlist
            elif count < mol_slice[0]:
                count += 1
                continue
            else:
                # important to append copy else still linked to orig generator
                mlist.append(copy.copy(m))
                try:
                    count = incrementer(count, mols, mol_slice[2])
                except StopIteration:
                    return mlist

        return mlist

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
