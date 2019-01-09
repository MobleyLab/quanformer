"""
helper.py

Helper functions used in the testing scripts.

"""

import openeye.oechem as oechem


def read_mol(infile, many=False):
    # many for multiple conformers, set False for single conf/mol
    ifs = oechem.oemolistream()
    if not ifs.open(infile):
        oechem.OEThrow.Fatal("Unable to open {} for reading".format(infile))
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    if many:
        return ifs.GetOEMols()
    mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, mol)
    ifs.close()
    return mol
