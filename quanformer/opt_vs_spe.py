#!/usr/bin/env python
"""
opt_vs_spe.py

By: Jessica Maat, Victoria T. Lim

TODO: double check that this works with only single input file

WARNING: If comparing many files, they should all be analogous!
 E.g., same number of conformers without filtering, so that
 reference conformer subtracted among all will be the same reference conf.

"""

import os
import sys
import openeye.oechem as oechem
import numpy as np
import proc_tags as pt
import math
from collections import defaultdict

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

### ------------------- Functions -------------------


def plotRMSD(rmsdict):
    """
    reference for 3d bar plot on stackoverflow- https://tinyurl.com/y9butx3l
    """
    zs = []
    # plot dimensions: x is mol, y is file, z is rmsd
    for keyFile, valueDict in rmsdict.items():
        fileData = sorted(valueDict.items())
        fileRMSDs = [y for _, y in fileData]
        zs.append(fileRMSDs)  # data should be all xs for y0, then for y1, etc

    # get list lengths and tick labels
    row_names = rmsdict.keys()  # file names, not sorted
    column_names = [x for x, _ in fileData]  # mol names, sorted
    numFiles = len(rmsdict)
    numMols = np.max([len(v) for v in rmsdict.values()])

    # set up mesh of positions
    xpos = np.arange(0, numMols, 1)
    ypos = np.arange(0, numFiles, 1)
    xpos, ypos = np.meshgrid(xpos + 0.25, ypos + 0.25)

    # convert positions to 1D array
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(numMols * numFiles)

    # lengths of bar widths, uniform in xy but diff in z
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = np.asarray(zs).flatten()  # plotted by (x,y0,z), (x,y1,z), ...

    # for each file being unique color
    values = np.linspace(0, 1., numFiles)
    colors = mpl.cm.tab20(values)
    colors = np.repeat(colors, numMols, axis=0)

    # for each mol being unique color
#    values = np.linspace(0, 1., numMols)
#    colors = mpl.cm.tab20(values)
#    colors = np.tile(colors,(numFiles,1))

    # draw figure
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color=colors)

    # format ticks, add labels
    ticksx = np.arange(0.5, numMols, 1)
    plt.xticks(ticksx, column_names, rotation=40, ha='right')
    ticksy = np.arange(0.5, numFiles, 1)
    plt.yticks(ticksy, row_names, rotation=-30, ha='left')
    for t in ax.zaxis.get_major_ticks():
        t.label.set_fontsize(14)

    ax.set_xlabel('Molecule')
#    ax.set_ylabel('QM theory')
    ax.set_zlabel('RMSD (kcal/mol)', fontsize=14)
    plt.savefig('barplot3d.png', bbox_inches='tight')
    plt.show()


def getRMSD(sdfRef, theory, rmsdict, package='Psi4'):
    """
    Perform RMSD calculation from an SDF file for molecule and its conformers.

    sdfRef: string, pathname of the SDF file with energies of opt 1 and opt 2
    theory: string, level of theory in format of mp2/6-31G*
    rmsdict: dictionary (can be empty) which will be populated in form of
             rmsdict[theory][molName] = 0.000  if the RMSD of before/after energies are 0.000
    package: string, name of software package used for QM calculation. only Psi4 currently supported

    """

    method, basis = theory.split('/')[0].strip(), theory.split('/')[1].strip()

    # create a molecule read in stream
    print("Opening SDF file %s" % sdfRef)
    ifsRef = oechem.oemolistream()
    ifsRef.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifsRef.open(sdfRef):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % sdfRef)
    # store all molecules in molsRef
    molsRef = ifsRef.GetOEMols()

    # create file object for output RMSD calculation
    RMSD = open("RMSD.txt", 'a')
    RMSD.write(
        "\nAnalyzing file: %s\n# Level of theory: %s\n" % (sdfRef, theory))

    # create file object for initial and final energies
    energies = open("energies_breakdown.txt", 'a')
    maximum = open("maxenergies.txt", "a")

    # Grab energies, perform RMSD calculation, write data to txt files.
    for rmol in molsRef:
        molName = rmol.GetTitle()
        tmol = np.asarray(
            pt.get_sd_list(rmol, 'QM opt energy', 'Psi4', method, basis),
            dtype=float)
        imol = np.asarray(
            pt.get_sd_list(rmol, 'QM opt energy initial', 'Psi4', method,
                           basis),
            dtype=float)
        final = tmol.copy()
        initial = imol.copy()

        # subtract conformer[0] energies from all conformers
        try:
            tmol -= tmol[0]
        except IndexError as e:
            sys.exit("No energies found for {} {}/{}! Check that data is \
stored in tags. Exiting.".format(rmol.GetTitle(), method, basis))
        imol -= imol[0]

        #subtracts initial minus final and sqaures all values
        fmol = np.subtract(tmol, imol)
        fmol = fmol[~np.isnan(fmol)]
        fmol = np.square(fmol)

        #sums all energies of conformers for given rmol and then takes average with respect to n-1 number of conformers
        tot = 0
        for n in fmol:
            tot += n
        average = math.sqrt(tot / (fmol.size - 1))

        #convert average from Hartree to Kcal/mol
        average = average * 627.5095

        # puts RMSD values into .txt file, and store in dict for plotting.
        RMSD.write("#%s\t%.5f RMSD(Kcal/mol)\n" % (molName, average))
        rmsdict[theory][molName] = average

        # store energies of initial and final for molecules conformers in energies.txt
        energies.write(
            "\n#%s\n#%s\n#RMSD = %.5f(y)\t\t(x=Hartree, y=kcal/mol)\n#conf. init. Energy(x)  \t final Energy(x) \t diff.(x)\tdiff. (y) \n"
            % (theory, molName, average))

        # get list of conformer indices to identify high RMSD ones
        conflist = pt.get_sd_list(rmol, "original index", package, method,
                                  basis)
        conformer = []
        for item in conflist:
            conformer.append(item.split(',')[0])  # append orig conf
        conformer = np.asarray(conformer, dtype=int)
        difference = np.array([])
        for i in range(len(tmol)):
            energies.write(
                "%r \t %5.9f \t %5.9f \t %5.9f\t%5.9f \n" %
                (conformer[i], initial[i], final[i], final[i] - initial[i],
                 (final[i] - initial[i]) * 627.5095))
            difference = np.append(difference,
                                   [(final[i] - initial[i]) * 627.5095])

        # find max 3 confs with highest RMSDs
        try:
            difference = np.absolute(difference)
            confmax1 = (np.nanargmax(difference))
            # set max conf to zero to find next highest
            difference[confmax1] = 0

            difference = np.absolute(difference)
            confmax2 = (np.nanargmax(difference))
            difference[confmax2] = 0

            difference = np.absolute(difference)
            confmax3 = (np.nanargmax(difference))
            difference[confmax3] = 0

            max1 = conformer[confmax1]
            max2 = conformer[confmax2]
            max3 = conformer[confmax3]
        except ValueError as e:
            #print("ValueError: {}".format(e))
            # TODO don't plot this mol for all nan's
            print("All RMSDs in list for file {} mol {} are nan!!!".format(
                sdfRef, molName))
            max1 = max2 = max3 = -1

        energies.write(
            "#*** Max energy differences are conformers (hi-->low): %r, %r, %r ***\n\n"
            % (max1, max2, max3))
        maximum.write(
            "%s, %s : %r, %r, %r\n" % (theory, molName, max1, max2, max3))

    maximum.close()
    RMSD.close()
    energies.close()

    return rmsdict


### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
        help="Required argument on name of text file with information on\
              levels of theory and associated filenames to process.\
              Format example: MP2/6-31G*, file1.sdf\
              See README file or examples for more details.")

    args = parser.parse_args()
    opt = vars(args)
    if not os.path.exists(opt['input']):
        raise parser.error("Input file %s does not exist." % opt['filename'])

    # Read input file and store each file's information in two lists.
    sdfList = []
    thryList = []
    with open(opt['input']) as f:
        for line in f:
            if line.startswith('#'):
                continue
            dataline = [x.strip() for x in line.split(',')]
            if dataline == ['']:
                continue
            thryList.append(dataline[0])
            sdfList.append(dataline[1])

    # Check that each file exists before starting
    while True:
        list1 = []
        for f in sdfList:
            list1.append(os.path.isfile(f))
        if all(list1):
            break
        else:
            raise parser.error("One or more analysis files are missing!!")

    xdict = lambda: defaultdict(xdict)
    rmsdict = xdict()
    for fname, theory in zip(sdfList, thryList):
        rmsdict = getRMSD(fname, theory, rmsdict)

    # make sure all mols exist before plotting
    mostMols = max(rmsdict.values(), key=len)  # method with most mols
    mols2del = []  # mols that are missing from other files
    for keyfile in rmsdict:  # make sure each file has mols of mostMols
        for keymol in mostMols:
            if type(rmsdict[keyfile][keymol]) == defaultdict:
                print("Skip plotting {} bc missing data from {}".format(
                    keymol, keyfile))
                mols2del.append(keymol)
    # delete missing mols from main dict
    if len(mols2del) != 0:
        rmsdict_orig = rmsdict.copy()  # just in case will need later
        for delmol in mols2del:
            for keyfile in rmsdict:
                rmsdict[keyfile].pop(delmol, None)
    plotRMSD(rmsdict)
