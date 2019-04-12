#!/usr/bin/env python

"""
cat_mols.py

Purpose: Combine molecules from separate files into single output file.
 - Can take multiple input molecule files
 - Can take input files of different extensions
 - Can ignore (remove) molecules in first file for provided list of molnames
 - Can be used with single input file if just want to remove some mols

Examples:
 - python cat_mols.py -i input1.mol2 input2.mol2 -o output.mol2
 - python cat_mols.py -i input1.sdf input2.mol2  -list exclude_from_input1.txt -o output.sdf
 - python cat_mols.py -i input1.sdf input2.mol2  -title molname_1 molname_2 -o output.sdf

Written by OpenEye, adapted by Victoria Lim
https://docs.eyesopen.com/toolkits/python/_downloads/catmols.py
 - Edited to remove the OEAddMols command to keep molecules separate
 - Edited to skip molecules with given titles from first file
 - Edited to sort final molecules by title from
   [https://docs.eyesopen.com/toolkits/cookbook/python/_downloads/ecdfcd69f00dc4e2d4e0f826e749b6b0/moldb_titlesort.py]

"""

#############################################################################
# Copyright (C) 2008-2015 OpenEye Scientific Software, Inc.
#############################################################################
# This program concatenates molecules into one file.
# It can be useful for generating ROCS queries or reattach ligands to an
# protein structure
#############################################################################
import re
import sys
from openeye.oechem import *

def CatMols(infnames, outfname,nameset):
    ofs = oemolostream()
    if not ofs.open(outfname):
        OEThrow.Fatal("Unable to open %s for writing" % outfname)

    for i, fname in enumerate(infnames):
        print(fname)
        ifs = oemolistream()
        if ifs.open(fname):
            for imol in ifs.GetOEGraphMols():
                if imol.GetTitle() in nameset and i==0:
                    continue
                else:
                    OEWriteMolecule(ofs, imol)
        else:
            OEThrow.Fatal("Unable to open %s for reading" % fname)

def natural_key(astr):
    """https://stackoverflow.com/questions/34518/natural-sorting-algorithm/34528"""
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', astr[0])]

def SortByTitle(ifs, ofs):
    moldb = OEMolDatabase(ifs)

    titles = [(t, i) for i, t in enumerate(moldb.GetTitles())]

    # VTL TODO: remove associated .idx file from sorting
    titles.sort(key=natural_key)

    indices = [i for t, i in titles]

    moldb.Order(indices)
    moldb.Save(ofs)


Interface = """
!BRIEF -i <infile1> [<infile2>...] -o <outfile>
!PARAMETER -i
  !ALIAS -in
  !TYPE string
  !LIST true
  !REQUIRED true
  !BRIEF input file name(s)
!END
!PARAMETER -o
  !ALIAS -out
  !TYPE string
  !REQUIRED true
  !BRIEF output file name
!END
!PARAMETER -title
  !ALIAS -t
  !TYPE string
  !LIST true
  !BRIEF Single or space-separated list of titles to exclude from parent file
!END
!PARAMETER -list
  !ALIAS -l
  !TYPE string
  !BRIEF List file of mol titles to exclude from parent file
!END
"""


def main(argv=[__name__]):
    itf = OEInterface(Interface, argv)

    # collect names
    nameset = set()
    if itf.HasString("-list"):
        try:
            lfs = open(itf.GetString("-list"))
        except IOError:
            OEThrow.Fatal("Unable to open %s for reading" % itf.GetString("-list"))
        for name in lfs.readlines():
            name = name.strip()
            nameset.add(name)
    elif itf.HasString("-title"):
        for t in itf.GetStringList("-title"):
            nameset.add(t)

    # perform the combination
    CatMols(itf.GetStringList("-i"), '____temp____.sdf', nameset)

    # sort molecules by title
    SortByTitle('____temp____.sdf', itf.GetString("-o"))

if __name__ == "__main__":
    sys.exit(main(sys.argv))
