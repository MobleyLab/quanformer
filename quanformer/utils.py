#!/usr/bin/env python
"""
utils.py

Simple tools for use with the quanformer package

Version:    Mar 29 2019
By:         Victoria T. Lim

"""
import openeye.oechem as oechem
import quanformer.reader as reader

def convert_extension(infile, outfile):

    # open input file
    mols = reader.read_mols(infile)

    # open output file
    ofs = oechem.oemolostream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    # write to output
    for mol in mols:
        for conf in mol.GetConfs():
            oechem.OEWriteConstMolecule(ofs, conf)

    # close filestreams
    ofs.close()
