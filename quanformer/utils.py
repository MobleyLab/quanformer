#!/usr/bin/env python
"""
utils.py

Simple tools for use with the quanformer package

Version:    Mar 29 2019
By:         Victoria T. Lim

"""
import openeye.oechem as oechem
import quanformer.reader as reader

def convert_extension(infile, outfile, canonical=False):
    """
    Convert one molecule file format into another using OpenEye tools.
    The user may also assign canonical smiles as name before writing output.

    """
    # open input file
    mols = reader.read_mols(infile)

    # open output file
    ofs = oechem.oemolostream()
    if not ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    # write to output
    for mol in mols:
        if canonical:
            smi = oechem.OEMolToSmiles(mol)
        for conf in mol.GetConfs():
            if canonical:
                conf.SetTitle(smi)
            oechem.OEWriteConstMolecule(ofs, conf)

    # close filestreams
    ofs.close()

def convert_extension_separate(infile, presuffix, canonical=False, separate='mol'):
    """
    Convert one molecule file format into another using OpenEye tools.
    The user may also assign canonical smiles as name before writing output.
    Separate output into (each mol with all confs) or (each conf).

    presuffix : list
        first item contains prefix of output name, last item contains extension
        ex. ['alkyl', '.xyz']
    separate : string
        'mol' or 'conf'

    """
    # open input file
    mols = reader.read_mols(infile)

    # write to output
    for i, mol in enumerate(mols):

        # open output file
        if separate=='mol':
            ofs = oechem.oemolostream()
            if not ofs.open('{}_{}{}'.format(presuffix[0], str(i), presuffix[1])):
                oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)
        if canonical:
            smi = oechem.OEMolToSmiles(mol)

        for j, conf in enumerate(mol.GetConfs()):
            # open output file
            if separate=='conf':
                ofs = oechem.oemolostream()
                if not ofs.open('{}_{}_{}{}'.format(presuffix[0], str(i), str(j), presuffix[1])):
                    oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)
            if canonical:
                conf.SetTitle(smi)
            oechem.OEWriteConstMolecule(ofs, conf)

    # close filestreams
    ofs.close()
