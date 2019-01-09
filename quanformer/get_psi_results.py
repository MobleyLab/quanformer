#!/usr/bin/env python
"""
get_psi_results.py

Purpose:    Process results from a Psi4 QM calculation.
By:         Victoria T. Lim
Version:    Oct 12 2018

NOTE:       Psi4 results are obtained by parsing output text file.
            Not sure if this is most efficient.
            JSON wrapper has some limitations though (as of Dec 2018).
            Also, not sure best way of parsing: via iterator (current
            approach) or by using some sort of find or regex command?

"""

import re
import os, sys, glob
import pickle
import openeye.oechem as oechem

# local testing vs. travis testing
try:
    import quanformer.proc_tags as pt
except ModuleNotFoundError:
    import proc_tags as pt  # VTL temporary bc travis fails to import

### ------------------- Functions -------------------


def initiate_dict():
    """
    Initiate conformer's data dict and set initial details

    Returns
    -------
    dictionary with keys for 'package' and 'missing'
        'package' is the name of the QM software package
        'missing' is whether the output file can be found

    """

    props = {}
    props['package'] = "Psi4"
    props['missing'] = False

    return props


def get_conf_data(props, calctype, timeout, psiout):
    """
    Call the relevant functions to process the output
    files for one calculation.

    Parameters
    ----------
    """

    # get wall clock time of the job
    props['time'] = get_psi_time(timeout)

    # get calculation details
    props = process_psi_out(psiout, props, calctype)

    return props


def set_conf_data(mol, props, calctype):

    # Set last coordinates from optimization. skip if missing.
    if 'coords' in props and len(props['coords']) != 0:
        mol.SetCoords(oechem.OEFloatArray(props['coords']))

    # Set SD tags for this molecule
    pt.set_sd_tags(mol, props, calctype)

    return mol


def check_title(mol, infile):
    """
    Check if mol has title which is required for setting up new subdirectories.
    If not title listed, rename molecule based on filename. Process filename
    to remove dashes and dots due to file name limitations of pipeline.

    Parameters
    ----------
    mol : OpenEye OEMol
    infile : string, name of the file from which the mol came

    Returns
    -------
    OpenEye OEMol

    """
    extract_fname = os.path.splitext(os.path.basename(infile))[0]
    extract_fname = extract_fname.replace("-", "")
    extract_fname = extract_fname.replace(".", "")
    if mol.GetTitle() == "":
        mol.SetTitle(extract_fname)

    return mol


def get_psi_time(filename):
    """
    Get wall-clock time from Psi4 file. If multiple times are present,
        the average will be taken. Used in CompareTimes(...) function.

    Parameters
    ----------
    filename: string name of the timefile. E.g. "timer.dat"

    Returns
    -------
    time: float of the average wall-clock time of a single timefile
        OR string of error message if file not found

    """

    # check whether file exists
    if not os.path.isfile(filename):
        print("*** ERROR: timer file not found: {} ***".format(filename))
        return ("Timer output file not found")

    # read file and extract time
    with open(filename) as fname:
        times = []
        for line in fname:
            if "Wall Time:" in line:
                times.append(float(line.split()[2]))
    time = sum(times) / float(len(times))

    return time


def get_scs_mp2(lines):
    """
    Get final SCS MP2 energy from output file. This is handled a little
    differently, outside of iterator because (1) have to do for both SPE
    and OPT, (2) SCS-MP2 may be listed multiple times, not just at end of file,
    so must be able to get the last one directly ideally without storing them
    all in a list first, (3) checking the qm method for each line seems silly,
    since SCS-MP2 energies would only be sought when method is mp2.

    This function called by the process_psi_out function.

    Note: MP2 Hessian output also has this but not sure if should be extracted
    since am not extracting optimized energies.

    Parameters
    ----------
    lines : list
        list of lines of Psi4 output file from Python readlines function

    Returns
    -------
    scs_ene : float
        value of the final SCS-MP2 total energy in original units of Psi4 (Har)

    """
    matching = [s for s in lines if "SCS Total Energy" in s]
    ene = float(matching[-1].split()[4])
    return ene


def process_psi_out(filename, properties, calctype='opt'):
    """
    Go through output file and get level of theory (method and basis set),
        number of optimization steps, initial and final energies, and
        optimized coordinates. Returns all this information in a dictionary
        that was passed to this function.

    Parameters
    ----------
    filename: string name of the output file. E.g. "output.dat"
    properties: dictionary where all the data will go. Can be empty or not.
    calctype: string; one of 'opt','spe','hess' for geometry optimization,
        single point energy calculation, or Hessian calculation

    Returns
    -------
    properties: dictionary with summarized data from output file.
        spe keys:  basis, method, finalEnergy
        opt keys:  basis, method, numSteps, initEnergy, finalEnergy, coords
        hess keys: basis, method, hessian

    """

    # check whether output file exists
    if not os.path.isfile(filename):
        print("*** ERROR: Output file not found: {} ***".format(filename))
        properties['missing'] = True
        return properties

    # open and read file
    f = open(filename, "r")
    lines = f.readlines()
    it = iter(lines)

    # process results for single point energy calculation
    if calctype == 'spe':
        for line in it:
            if "set basis" in line:
                properties['basis'] = line.split()[2]
            if "energy(" in line:
                properties['method'] = line.split('\'')[1]
            # this line should only show up once for spe
            if "Total Energy =" in line:
                properties['finalEnergy'] = float(line.split()[3])

        # check for scs-mp2 energy if applicable
        if properties['method'].lower() == 'mp2':
            scs_ene = get_scs_mp2(lines)
            properties['finalSCSEnergy'] = scs_ene

        return properties

    # process results for Hessian calculation
    elif calctype == 'hess':
        import math
        import numpy as np
        hess = []
        for line in it:
            if "set basis" in line:
                properties['basis'] = line.split()[2]
            if "=hessian(" in ''.join(line.split()):  # rm spaces around eql
                properties['method'] = line.split('\'')[1]
            if "## Hessian" in line:
                line = next(it)  # "Irrep:"
                line = next(it)  # blank line
                line = next(it)  # column index labels
                line = next(it)  # blank line

                # convert the rest of file/iterator to list
                rough = list(it)
                # remove last 5 lines: four blank, one exit message
                rough = rough[:-5]
                # convert strings to floats
                rough = [[float(i) for i in line.split()] for line in rough]
                # find blank sublists
                empty_indices = [i for i, x in enumerate(rough) if not x]
                # blank sublists are organized in groups of three lines:
                # (1) blank (2) column labels (3) blank. empty_ind has (1) and (3)
                # delete in reverse to maintain index consistency
                num_chunks = 1  # how many chunks psi4 printed the Hessian into
                for i in reversed(range(0, len(empty_indices), 2)):
                    del rough[empty_indices[i + 1]]  # (3)
                    del rough[empty_indices[i] + 1]  # (2)
                    del rough[empty_indices[i]]  # (1)
                    num_chunks += 1
                # now remove first element of every list of row label
                _ = [l.pop(0) for l in rough]
                # get dimension of Hessian (3N for 3Nx3N matrix)
                three_n = int(len(rough) / num_chunks)
                # concatenate the chunks
                hess = np.array([]).reshape(three_n, 0)
                for i in range(num_chunks):
                    beg_ind = i * three_n
                    end_ind = (i + 1) * three_n
                    chunk_i = np.array(rough[beg_ind:end_ind])
                    hess = np.concatenate((hess, chunk_i), axis=1)
                # check that final matrix is symmetric
                if not np.allclose(hess, hess.T, atol=0):
                    print("ERROR: Quanformer did not read symmetric Hessian "
                          "from Psi4 output file")
                    properties['hessian'] = "Hessian not found in output file"
                    return properties
                properties['hessian'] = hess
            # the Hessian pattern was never found
            else:
                properties['hessian'] = "Hessian not found in output file"

        return properties

    # process results for geometry optimization
    # loop through file to get method, basis set, numSteps, energies, coords
    rough = []
    coords = []
    for line in it:
        if "set basis" in line:
            properties['basis'] = line.split()[2]
        if "optimize(" in line:
            properties['method'] = line.split('\'')[1]
        if "Optimization is complete" in line:
            properties['numSteps'] = line.strip().split(' ')[5]
            for _ in range(8):
                line = next(it)
            properties['initEnergy'] = float(line.split()[1])
        if "Final energy" in line:
            properties['finalEnergy'] = float(line.split()[3])
            line = next(it)  # "Final (previous) structure:"
            line = next(it)  # "Cartesian Geometry (in Angstrom)"
            line = next(it)  # Start of optimized geometry
            while "Saving final" not in line:
                rough.append(line.split()[1:4])
                line = next(it)

    # check for scs-mp2 energy if applicable
    if properties['method'].lower() == 'mp2':
        scs_ene = get_scs_mp2(lines)
        properties['finalSCSEnergy'] = scs_ene

    # convert the 2D (3xN) coordinates to 1D list of length 3N (N atoms)
    for atomi in rough:
        coords += [float(i) for i in atomi]
    properties['coords'] = coords
    f.close()
    return properties


### ------------------- Script -------------------


def get_psi_results(origsdf,
                    finsdf,
                    calctype='opt',
                    psiout="output.dat",
                    timeout="timer.dat"):
    """
    Read in OEMols (and each of their conformers) in origsdf file,
        get results from Psi4 calculations in the same directory as origsdf,
        and write out results into finsdf file.
    Directory layout is .../maindir/molName/confNumber/outputfiles .
    Both origsdf and finsdf are located in maindir.

    Parameters
    ----------
    origsdf:  string - original SDF file of input structures of QM calculation
    finsdf:   string - full name of final SDF file with optimized results.
    calctype: string; one of 'opt','spe','hess' for geometry optimization,
        single point energy calculation, or Hessian calculation
    psiout:   string - name of the Psi4 output files. Default is "output.dat"
    timeout: string - name of the Psi4 timer files. Default is "timer.dat"

    Returns
    -------
    OpenEye OEMol with data in SD tags

    None is returned if the function returns early (e.g., if output file
       already exists) or if there is KeyError from processing last
       iteration of output file (last conf of last mol).

    """

    hdir, fname = os.path.split(origsdf)
    wdir = os.getcwd()

    # check that specified calctype is valid
    if calctype not in {'opt', 'spe', 'hess'}:
        sys.exit("Specify a valid calculation type.")

    # read in sdf file and distinguish each molecule's conformers
    ifs = oechem.oemolistream()
    ifs.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifs.open(origsdf):
        sys.exit("Unable to open %s for reading" % origsdf)
    molecules = ifs.GetOEMols()

    # open outstream file
    writeout = os.path.join(wdir, finsdf)
    write_ofs = oechem.oemolostream()
    if os.path.exists(writeout):
        print("File already exists: %s. Skip getting results.\n" % (finsdf))
        return (None, None)
    if not write_ofs.open(writeout):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % writeout)

    # Hessian dictionary, where hdict['molTitle']['confIndex'] has np array
    if calctype == 'hess':
        hdict = {}

    # for each conformer, process output file and write new data to SDF file
    for mol in molecules:
        print("===== %s =====" % (mol.GetTitle()))
        if calctype == 'hess':
            hdict[mol.GetTitle()] = {}

        for j, conf in enumerate(mol.GetConfs()):

            props = initiate_dict()

            # set file locations
            timef = os.path.join(hdir,
                                 "%s/%s/%s" % (mol.GetTitle(), j + 1, timeout))
            outf = os.path.join(hdir,
                                "%s/%s/%s" % (mol.GetTitle(), j + 1, psiout))

            # process output and get dictionary results
            props = get_conf_data(props, calctype, timef, outf)

            # if output was missing or are missing calculation details
            # move on to next conformer
            if props['missing'] or (calctype == 'opt' and not all(
                    key in props
                    for key in ['numSteps', 'finalEnergy', 'coords'])):
                print(
                    "ERROR reading {}\nEither Psi4 job was incomplete OR wrong calctype specified\n"
                    .format(outf))
                continue

            # add data to oemol
            conf = set_conf_data(conf, props, calctype)

            # if hessian, append to dict bc does not go to SD tag
            if calctype == 'hess':
                hdict[mol.GetTitle()][j + 1] = props['hessian']

            # check mol title
            conf = check_title(conf, origsdf)

            # write output file
            oechem.OEWriteConstMolecule(write_ofs, conf)

    # if hessian, write hdict out to separate file
    if calctype == 'hess':
        hfile = os.path.join(wdir,
                             os.path.splitext(finsdf)[0] + '.hess.pickle')
        pickle.dump(hdict, open(hfile, 'wb'))

    # close file streams
    ifs.close()
    write_ofs.close()
    try:
        return props['method'], props['basis']
    except KeyError:
        return None, None


def getPsiOne(infile,
              outfile,
              calctype='opt',
              psiout="output.dat",
              timeout="timer.dat"):
    """
    Write out Psi4 optimized geometry details into a new OEMol.

    Parameters
    ----------
    infile : string
        Name of input geometry file THAT WAS USED TO WRITE THE PSI INPUT.
        To ensure that the atom orderings remain constant.
    outfile : string
        Name of output geometry file with optimized results.
    calctype: string; one of 'opt','spe','hess' for geometry optimization,
        single point energy calculation, or Hessian calculation
    psiout : string
        Name of the Psi4 output files. Default is "output.dat"
    timeout : string
        Name of the Psi4 timer files. Default is "timer.dat"

    Returns
    -------
    method: string - QM method from Psi4 calculations
    basisset: string - QM basis set from Psi4 calculations

    None is returned if the function returns early (e.g., if output file
       already exists)

    """
    # check that specified calctype is valid
    if calctype not in {'opt', 'spe', 'hess'}:
        sys.exit("Specify a valid calculation type.")

    # Read in SINGLE MOLECULE .sdf file
    ifs = oechem.oemolistream()
    mol = oechem.OEGraphMol()
    if not ifs.open(infile):
        sys.exit("Unable to open %s for reading" % origsdf)
    oechem.OEReadMolecule(ifs, mol)

    # Open outstream file.
    write_ofs = oechem.oemolostream()
    if os.path.exists(outfile):
        print("File already exists: %s. Skip getting results.\n" % (outfile))
        return (None, None)
    if not write_ofs.open(outfile):
        oechem.OEThrow.Fatal("Unable to open %s for writing" % outfile)

    props = initiate_dict()

    # process output and get dictionary results
    props = get_conf_data(props, calctype, timeout, psiout)

    # if output was missing or are missing calculation details
    # move on to next conformer
    if props['missing'] or (calctype == 'opt' and not all(
            key in props for key in ['numSteps', 'finalEnergy', 'coords'])):
        sys.exit(
            "ERROR reading {}\nEither Psi4 job was incomplete OR wrong calctype specified\n"
            .format(psiout))

    # add data to oemol
    mol = set_conf_data(mol, props, calctype)

    # if hessian, write hdict out to separate file
    if calctype == 'hess':
        hfile = os.path.join(os.path.splitext(outfile)[0] + '.hess.pickle')
        pickle.dump(props['hessian'], open(hfile, 'wb'))

    # check mol title
    mol = check_title(mol, infile)

    # write output file
    oechem.OEWriteConstMolecule(write_ofs, mol)
    write_ofs.close()

    return mol, props


if __name__ == "__main__":
    get_psi_results(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                    sys.argv[5])
