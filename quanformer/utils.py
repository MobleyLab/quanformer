#!/usr/bin/env python
"""
utils.py

Simple tools for use with the quanformer package

Version:    Mar 29 2019
By:         Victoria T. Lim

"""
import openeye.oechem as oechem

def convert_extension(infile, outfile):

    # open input file
    ifs = oechem.oemolistream()
    if not ifs.open(infile):
        oechem.OEThrow.Fatal("Unable to open {} for reading".format(infile))
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())

    # open output file
    ofs = oechem.oemolostream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    # write to output
    for mol in ifs.GetOEMols():
        for conf in mol.GetConfs():
            oechem.OEWriteConstMolecule(ofs, conf)

    # close filestreams
    ifs.close()
    ofs.close()
