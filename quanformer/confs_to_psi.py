#!/usr/bin/env python
"""
confs_to_psi.py

Purpose:    Generate Psi4 inputs for many molecules/conformers.
By:         Victoria T. Lim, Christopher I. Bayly
Version:    Nov 30 2018

"""

import os, sys
import openeye.oechem as oechem
import shutil
import json


def make_psi_input(mol, label, method, basisset, calctype='opt', mem=None):
    """
    Get coordinates from input mol, and generate/format input text for
    Psi4 calculation.

    Parameters
    ----------
    mol : OpenEye OEMol
        OEMol with coordinates
    label : string
        Name of molecule with integer identifier (for conformers).
    method: string
        Name of the method as understood by Psi4. Example: "mp2"
    basis : string
        Name of the basis set as understood by Psi4. Example: "def2-sv(p)"
    calctype : string
        What kind of Psi4 calculation to run. Supported inputs are:
        'opt' for geometry optimization,
        'spe' for single point energy calculation, and
        'hess' for Hessian calculation.
    memory : string
        How much memory each Psi4 job should take. If not specified, the
        default in Psi4 is 500 Mb. Examples: "2000 MB" "1.5 GB"
        http://www.psicode.org/psi4manual/master/psithoninput.html

    Returns
    -------
    inputstring : string
        Contents of Psi4 input file to be written out

    """

    # check that specified calctype is valid
    if calctype not in {'opt', 'spe', 'hess'}:
        sys.exit("Specify a valid calculation type.")

    inputstring = ""

    # specify memory requirements, if defined
    if mem != None:
        inputstring += "memory %s\n" % mem
    inputstring += ('molecule %s {\n' % label)

    # charge and multiplicity; multiplicity hardwired to singlet (usually is)
    netCharge = oechem.OENetCharge(mol)
    inputstring += ('  %s 1' % netCharge)

    # get atomic symbol and coordinates of each atom
    xyz = oechem.OEFloatArray(3)
    for atom in mol.GetAtoms():
        mol.GetCoords(atom, xyz)
        inputstring+=( '\n  %s %10.4f %10.4f  %10.4f' \
                       %(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()),
                       xyz[0], xyz[1], xyz[2]) )
    inputstring += ('\n  units angstrom\n}')

    # check if mol has a "freeze" tag
    for x in oechem.OEGetSDDataPairs(mol):
        if calctype == "opt" and "atoms to freeze" in x.GetTag():
            b = x.GetValue()
            y = b.replace("[", "")
            z = y.replace("]", "")
            a = z.replace(" ", "")
            freeze_list = a.split(",")
            inputstring += (
                "\n\nfreeze_list = \"\"\"\n  {} xyz\n  {} xyz\n  {} "
                "xyz\n  {} xyz\n\"\"\"".format(freeze_list[0], freeze_list[1],
                                               freeze_list[2], freeze_list[3]))
            inputstring += "\nset optking frozen_cartesian $freeze_list"
            inputstring += (
                "\nset optking dynamic_level = 1\nset optking "
                "consecutive_backsteps = 2\nset optking intrafrag_step_limit = "
                "0.1\nset optking interfrag_step_limit = 0.1\n")

    # best practices for scf calculations
    # http://www.psicode.org/psi4manual/master/scf.html#recommendations
    # http://www.psicode.org/psi4manual/master/dft.html#recommendations
    inputstring += '\n\nset scf_type df'
    inputstring += '\nset guess sad'

    # explicitly specify MP2 RI-auxiliary basis for [Ahlrichs] basis set
    # http://www.psicode.org/psi4manual/master/basissets_byfamily.html
    # DFMP2 *should* get MP2 aux sets fine for [Pople and Dunning] sets
    # http://www.psicode.org/psi4manual/master/dfmp2.html
    if method.lower() == 'mp2' and 'def2' in basisset:
        if basisset.lower() == 'def2-sv(p)':
            inputstring += ('\nset df_basis_mp2 def2-sv_p_-ri')
        elif basisset.lower() != 'def2-qzvpd':  # no aux set for qzvpd 10-6-18
            inputstring += ('\nset df_basis_mp2 %s-ri' % (basisset))

    inputstring += ('\n\nset basis %s' % (basisset))
    inputstring += ('\nset freeze_core True')
    # specify command for type of calculation
    if calctype == 'opt':
        inputstring += ('\noptimize(\'%s\')\n\n' % (method))
    elif calctype == 'spe':
        inputstring += ('\nenergy(\'%s\')\n\n' % (method))
    elif calctype == 'hess':
        inputstring += (
            '\nH, wfn = hessian(\'%s\', return_wfn=True)\nwfn.hessian().print_out()\n\n'
            % (method))

    return inputstring


def make_psi_json(mol, label, method, basisset, calctype='opt', mem=None):
    """
    THIS FUNCTION IS A WORK IN PROGRESS.

    Get coordinates from input mol, and generate/format input text for
    Psi4 calculation via JSON wrapper.

    Parameters
    ----------
    mol : OpenEye OEMol
        OEMol with coordinates
    label : string
        Name of molecule with integer identifier (for conformers).
    method: string
        Name of the method as understood by Psi4. Example: "mp2"
    basis : string
        Name of the basis set as understood by Psi4. Example: "def2-sv(p)"
    calctype : string
        What kind of Psi4 calculation to run. Supported inputs are:
        'opt' for geometry optimization,
        'spe' for single point energy calculation, and
        'hess' for Hessian calculation.
    memory : string
        How much memory each Psi4 job should take. If not specified, the
        default in Psi4 is 500 Mb. Examples: "2000 MB" "1.5 GB"
        http://www.psicode.org/psi4manual/master/psithoninput.html

    Returns
    -------
    inputstring : string
        Contents of Psi4 input file to be written out

    """
    # check that specified calctype is valid
    if calctype not in {'opt', 'spe', 'hess'}:
        sys.exit("Specify a valid calculation type.")

    inputdict = {}
    moldict = {}
    modeldict = {}
    keydict = {}
    inputdict["schema_name"] = "qc_schema_input"
    inputdict["schema_version"] = 1

    # specify memory requirements, if defined
    if mem != None:
        inputdict["memory"] = mem

    #TODO -- json version
    # charge and multiplicity; multiplicity hardwired to singlet (usually is)
    #inputdict["charge"] = oechem.OENetCharge( mol)
    #inputdict["multiplicity"] = 1

    # get atomic symbol and coordinates of each atom
    geom_list = []
    elem_list = []
    xyz = oechem.OEFloatArray(3)
    for atom in mol.GetAtoms():
        mol.GetCoords(atom, xyz)
        geom_list.append(xyz[0])
        geom_list.append(xyz[1])
        geom_list.append(xyz[2])
        elem_list.append(oechem.OEGetAtomicSymbol(atom.GetAtomicNum()))
    moldict["geometry"] = geom_list
    moldict["symbols"] = elem_list
    inputdict["molecule"] = moldict

    #TODO -- json version
    ## check if mol has a "freeze" tag
    #for x in oechem.OEGetSDDataPairs(mol):
    #    if calctype=="opt" and "atoms to freeze" in x.GetTag():
    #        b = x.GetValue()
    #        y = b.replace("[", "")
    #        z = y.replace("]", "")
    #        a = z.replace(" ", "")
    #        freeze_list = a.split(",")
    #        inputstring += ("\n\nfreeze_list = \"\"\"\n  {} xyz\n  {} xyz\n  {} "
    #                       "xyz\n  {} xyz\n\"\"\"".format(freeze_list[0],
    #                       freeze_list[1], freeze_list[2], freeze_list[3]))
    #        inputstring += "\nset optking frozen_cartesian $freeze_list"
    #        inputstring += ("\nset optking dynamic_level = 1\nset optking "
    #            "consecutive_backsteps = 2\nset optking intrafrag_step_limit = "
    #            "0.1\nset optking interfrag_step_limit = 0.1\n")

    #TODO -- json version
    ## explicitly specify MP2 RI-auxiliary basis for Ahlrichs basis set
    ## http://www.psicode.org/psi4manual/master/basissets_byfamily.html
    ## DFMP2 *should* get MP2 aux sets fine for Pople/Dunning
    ## http://www.psicode.org/psi4manual/master/dfmp2.html
    #if method.lower()=='mp2' and 'def2' in basisset:
    #    if basisset.lower()=='def2-sv(p)':
    #        inputstring+=('\nset df_basis_mp2 def2-sv_p_-ri')
    #    elif basisset.lower()!='def2-qzvpd':  # no aux set for qzvpd 10-6-18
    #        inputstring+=('\nset df_basis_mp2 %s-ri' % (basisset))

    modeldict["basis"] = basisset
    modeldict["method"] = method
    inputdict["model"] = modeldict

    #inputstring+=('\nset freeze_core True')
    # specify command for type of calculation
    if calctype == 'opt':
        # TODO
        pass
    elif calctype == 'spe':
        inputdict["driver"] = 'energy'
    elif calctype == 'hess':
        inputdict["driver"] = 'hessian'
        keydict["return_wfn"] = True
    inputdict["keywords"] = keydict

    return inputdict


def confs_to_psi(insdf,
                 method,
                 basis,
                 calctype='opt',
                 memory=None,
                 via_json=False):
    """
    Read in molecule(s) (and conformers, if present) in insdf file. Create
    Psi4 input calculations for each structure.

    Parameters
    ----------
    insdf: string
        Name of the molecule file for which to create Psi4 input file.
        SDF format can contain multiple molecules and multiple conformers per
        molecule in a single file.
    method: string
        Name of the method as understood by Psi4. Example: "mp2"
    basis : string
        Name of the basis set as understood by Psi4. Example: "def2-sv(p)"
    calctype : string
        What kind of Psi4 calculation to run. Supported inputs are:
        'opt' for geometry optimization,
        'spe' for single point energy calculation, and
        'hess' for Hessian calculation.
    memory : string
        How much memory each Psi4 job should take. If not specified, the
        default in Psi4 is 500 Mb. Examples: "2000 MB" "1.5 GB"
        http://www.psicode.org/psi4manual/master/psithoninput.html
    via_json : Boolean
        If True, use JSON wrapper for Psi4 input and output.
        - Psi4 input would be in "input.py", called with python
        - Psi4 output would be in "output.json"
        If False, use normal text files for Psi4 input and output.
        - Psi4 input would be in "input.dat"
        - Psi4 output would be in "output.dat"
    """
    wdir = os.getcwd()

    ### Read in .sdf file and distinguish each molecule's conformers
    ifs = oechem.oemolistream()
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifs.open(insdf):
        oechem.OEThrow.Warning("Unable to open %s for reading" % insdf)
        return

    ### For each molecule: for each conf, generate input
    for mol in ifs.GetOEMols():
        print(mol.GetTitle(), mol.NumConfs())
        if not mol.GetTitle():
            sys.exit("ERROR: OEMol must have title assigned! Exiting.")
        for i, conf in enumerate(mol.GetConfs()):
            # change into subdirectory ./mol/conf/
            subdir = os.path.join(wdir, "%s/%s" % (mol.GetTitle(), i + 1))
            if not os.path.isdir(subdir):
                os.makedirs(subdir)
            if os.path.exists(os.path.join(subdir, 'input.dat')):
                print("Input file already exists. Skipping.\n{}\n".format(
                    os.path.join(subdir, 'input.dat')))
                continue
            label = mol.GetTitle() + '_' + str(i + 1)
            if via_json:
                ofile = open(os.path.join(subdir, 'input.py'), 'w')
                ofile.write("# molecule {}\n\nimport numpy as np\nimport psi4"
                            "\nimport json\n\njson_data = ".format(label))
                json.dump(
                    make_psi_json(conf, label, method, basis, calctype,
                                  memory),
                    ofile,
                    indent=4,
                    separators=(',', ': '))
                ofile.write(
                    "\njson_ret = psi4.json_wrapper.run_json(json_data)\n\n")
                ofile.write("with open(\"output.json\", \"w\") as ofile:\n\t"
                            "json.dump(json_ret, ofile, indent=2)\n\n")
            else:
                ofile = open(os.path.join(subdir, 'input.dat'), 'w')
                ofile.write(
                    make_psi_input(conf, label, method, basis, calctype,
                                   memory))
            ofile.close()
    ifs.close()
