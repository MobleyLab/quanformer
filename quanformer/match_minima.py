#!/usr/bin/env python
"""
match_minima.py

By:     Victoria T. Lim
"""

import os
import sys
import openeye.oechem as oechem
import numpy as np
from numpy import nan
import pickle
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
import proc_tags as pt  # for get_sd_list

### ------------------- Functions -------------------


def compare_two_mols(rmol, qmol):
    """
    For two identical molecules, with varying conformers,
        make an M by N comparison to match the M minima of
        rmol to the N minima of qmol. Match is declared
        for lowest RMSD between the two conformers and
        if the RMSD is below 0.5 Angstrom.

    Parameters
    ----------
    rmol:       reference OEChem molecule with all its filtered conformers
    qmol:       query OEChem molecule with all its filtered conformers

    Returns
    -------
    molIndices: 1D list of qmol conformer indices that correspond to rmol confs

    """

    automorph = True  # take into acct symmetry related transformations
    heavyOnly = False  # do consider hydrogen atoms for automorphisms
    overlay = True  # find the lowest possible RMSD

    molIndices = []  # 1D list, stores indices of matched qmol confs wrt rmol

    for Rconf in rmol.GetConfs():
        print(">>>> Matching %s conformers to minima: %d <<<<"\
            % (qmol.GetTitle(),Rconf.GetIdx()+1))

        # for this Rconf, calculate/store RMSDs with all of qmol's conformers
        rsublist = []
        for Qconf in qmol.GetConfs():
            rms = oechem.OERMSD(Rconf, Qconf, automorph, heavyOnly, overlay)
            rsublist.append(rms)

        # for this Rconf, get qmol conformer index for minimum RMSD
        thisMin = [i for i, j in enumerate(rsublist) if j == min(rsublist)][0]
        if rsublist[thisMin] <= 0.5:
            molIndices.append(thisMin)
        else:
            print('no match bc rmsd is ', rsublist[thisMin])
            molIndices.append(None)

    return molIndices


def plot_mol_minima(molName, minimaE, xticklabels, selected=None, stag=False):
    '''

    selected: list of indexes for methods to be plotted. 0-based index.
              consider also renaming figname to not overwrite existing plots.
    stag: Boolean, True to stagger plots to better see line trends.
          Works best with few (<4?) numFiles.
          Be very wary of using this to not confuse energy distributions.

    # minimaE for ONE molecule
    '''

    refNumConfs = len(minimaE[0])
    refFile = xticklabels[0]
    numFiles = len(minimaE)

    ### Flatten this 2D list into a 1D to find min and max for plot
    flatten = [item for sublist in minimaE for item in sublist]
    floor = min(flatten)
    ceiling = max(flatten)
    if (ceiling - floor) > 4.0:
        ystep = (ceiling - floor) / 9  # have 10 increments of y-axis
        ystep = round(ystep * 2) / 2  # round the step to nearest 0.5
    else:
        ystep = (ceiling - floor)

    ### Stagger each of the component files of minimaE for ease of viewing.
    if stag == True:
        tempMinimaE = []
        for i, fileE in enumerate(minimaE):
            tempMinimaE.append([x + i / 2. for x in fileE])
        minimaE = tempMinimaE
        ceiling = ceiling + numFiles

    ### Figure labels.
    plttitle = "Relative Energies of %s Minima" % (molName)
    plttitle += "\nby Reference File %s" % (refFile)
    ylabel = "Relative energy (kcal/mol)"
    figname = "minimaE_%s.png" % (molName)

    ### Set x-axis values and labels. Can either be letter or number.

    # LETTER LABELS
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'  # label x-axis by letter
    rpt = int((len(minimaE[0]) / 26) + 1)
    xlabs = [''.join(i)
             for i in itertools.product(letters, repeat=rpt)][:refNumConfs]
    # OR, NUMBER LABELS
    #xlabs = range(len(minimaE[0]))

    ### Set this up for grid.
    fig = plt.figure(figsize=(20, 10))
    ax = fig.gca()
    ax.set_xticks(np.arange(-1, refNumConfs + 1, 2))

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle, fontsize=16)
    plt.ylabel(ylabel, fontsize=20)
    plt.xlabel("conformer minimum", fontsize=20)
    plt.xticks(list(range(refNumConfs)), xlabs, fontsize=18)
    plt.yticks(fontsize=18)

    ### Plot the data.
    colors = mpl.cm.rainbow(np.linspace(0, 1, numFiles))
    markers = [
        "x", "^", "8", "d", "o", "s", "*", "p", "v", "<", "D", "+", ">", "."
    ] * 10
    #for i in reversed(range(numFiles)):
    for i, FileE in enumerate(minimaE):
        if selected is not None and i not in selected:
            continue
        xi = list(range(refNumConfs))
        #yi = [item[i] for item in minimaE]
        yi = FileE
        plt.plot(xi,yi,color=colors[i],label=xticklabels[i],\
            marker=markers[i],markersize=9)

    ### Add legend and set plot limits.
    plt.legend(bbox_to_anchor=(0.96, 1), loc=2, prop={'size': 18})
    plt.xlim(-1, refNumConfs + 1)
    # y axis limits: min, max, step
    ax.set_yticks(
        np.arange(int(round(floor)) - 2,
                  int(round(ceiling)) + 2, ystep))
    plt.grid()

    plt.savefig(figname, bbox_inches='tight')
    #    plt.show()
    plt.clf()


def plot_avg_times(molName, avgTimes, sdTimes, xticklabels):
    plttitle = "Conformer-Averaged Wall Times\nfor %s" % (molName)
    plttitle += "\nGeometry Optimization in Psi4"
    ylabel = "time (s)"
    figname = "timebars_%s.png" % molName
    x = list(range(len(avgTimes)))

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle, fontsize=20)
    plt.ylabel(ylabel, fontsize=18)
    plt.xticks(x, xticklabels, fontsize=14, rotation=-30, ha='left')
    plt.yticks(fontsize=14)

    ### Plot the data.
    colors = mpl.cm.rainbow(np.linspace(0, 1, len(x)))
    plt.bar(
        x, avgTimes, color=colors, align='center', yerr=sdTimes, ecolor='k')
    plt.savefig(figname, bbox_inches='tight')
    #    plt.show()
    plt.clf()


def match_minima(sdfList, thryList, *tags):
    """
    For list of SDF files, match the conformer minima to those of the reference
       SDF file. Ex. Conf G of reference file matches with conf R of file3.

    Parameters
    ----------
    sdfList: str list - list of the SDF file names to be analyzed.
          This list should include reference SDF file (sdfRef) as first element.
    thryList: str list - list of levels of theory corresponding to the files
          in sdfList. E.g., ['MP2/def2-TZVP','B3LYP-D3MBJ/6-311++G**']
    tags: variable number of arguments for the SD tag to get during matching

    Returns
    -------
    moldict: nested dictionary. example format:
        moldict['molName']['QM opt energy']
          = 2d list, [[file1: conf1 conf2] [file2: conf1 conf2]]
        Second-level keys are from tags variable.

    """

    def load_file(fname):
        ifs = oechem.oemolistream()
        ifs.SetConfTest(oechem.OEAbsoluteConfTest())
        if not ifs.open(fname):
            oechem.OEThrow.Fatal("Unable to open %s for reading" % fname)
        mols = ifs.GetOEMols()
        return mols

    sdfRef = sdfList[0]
    allIndices = []  # for M mols, N reference minima of each mol, P matching indices for each ref minimia
    moldict = {}  # nested dictionary with 1st layer of mol names, 2nd layer of properties (energies, opt times, etc)

    for i, sdfQuery in enumerate(sdfList):
        qthry = thryList[i]
        qmethod = qthry.split('/')[0].strip()
        qbasis = qthry.split('/')[1].strip()

        print("\n\nOpening reference file %s" % sdfRef)
        molsRef = load_file(sdfRef)

        print("Opening query file %s, and using [ %s ] energies" % (sdfQuery,
                                                                    qthry))
        molsQuery = load_file(sdfQuery)

        # loop over each molecule in reference file and in query file
        for rmol in molsRef:
            name = rmol.GetTitle()
            noMatch = True
            for qmol in molsQuery:
                if rmol.GetTitle() == qmol.GetTitle():
                    noMatch = False  # only continue match for same mols
                    break
            # before reaching verdict of match, check tag holders in moldict
            if name not in moldict: moldict[name] = {}
            for t in tags:
                if t not in moldict[name]: moldict[name][t] = []
            if 'indices' not in moldict[name]: moldict[name]['indices'] = []
            if 'refNumConfs' not in moldict[name]:
                moldict[name]['refNumConfs'] = []
            # no match was found; don't continue match function
            if noMatch:
                moldict[name]['indices'].append([-2] * rmol.NumConfs())
                for t in tags:
                    moldict[name][t].append([nan] * rmol.NumConfs())
                print(
                    "No %s molecule found in %s" % (rmol.GetTitle(), sdfQuery))
                # reset molsQuery generator
                molsQuery = load_file(sdfQuery)
                continue

            # get data from tags
            for t in tags:
                if t not in moldict[name]: moldict[name][t] = []
                moldict[name][t].append(
                    list(
                        map(float,
                            pt.get_sd_list(qmol, t, 'Psi4', qmethod, qbasis))))

            # Skip minmatch if this query file is same as reference file;
            #    before skip, get data for elists, refNumConfs, allIndices.
            if sdfQuery == sdfRef:
                print("\nSkipping comparison against self.")
                moldict[name]['refNumConfs'].append(rmol.NumConfs())
                moldict[name]['indices'].append([-1] * rmol.NumConfs())
                continue

            # get indices of qmol conformers that match rmol conformers
            molIndices = compare_two_mols(rmol, qmol)
            moldict[name]['indices'].append(molIndices)

    return moldict


def calc_rms_error(trimE, zeroes):
    """
    From relative energies with respect to some conformer from calc_rel_ene,
       calculate the root mean square error with respect to the relative
       conformer energies of the first (reference) file.

    Parameters
    ----------
    trimE: 3D list of energies, where trimE[i][j][k] represents the
      ith molecule, jth file, kth conformer rel energy
    zeroes: a 1D list of index of the reference conformer per each mol

    Returns
    -------
    relByFile: a 1D list of RMS errors for each file with reference to
      first input file

    """
    relByFile = []
    for i, molist in enumerate(trimE):
        molEnes = []
        for j, filelist in enumerate(molist):
            errs = np.asarray(filelist) - np.asarray(
                molist[0])  # subtract ref file
            sqrs = errs**2.  # squared
            sqrs = np.delete(
                sqrs, zeroes[i])  # delete reference conformer (zero rel ene)
            sqrs = sqrs[
                ~np.isnan(sqrs)]  # delete nan values to get an rmse********
            mse = np.mean(sqrs)
            rmse = np.sqrt(mse)
            molEnes.append(rmse)
        relByFile.append(molEnes)

    return relByFile


def get_ratio_times(allMolTimes, zeroes):
    """
    From all molecule times, calculate relative time ratios for matched minima.
       If a conf has nan or is not matched, that time is not considered.
       After dividing by reference time, the matched conformer files are
       averaged for a particular file opt.

    Parameters
    ----------
    allMolTimes: 3D list of times, where allMolTimes[i][j][k] represents the
      ith molecule, jth file, kth conformer time (sec)

    Returns
    -------
    relByFile: a 1D list of times ratios for each file with reference to
      first input file
    sdByFile: a 1D list of standard deviation of conformer-averaged times
      relative to first input file

    """

    relByFile = []
    sdByFile = []
    for i, molist in enumerate(allMolTimes):
        molTimes = []
        molStds = []
        for j, filelist in enumerate(molist):
            rels = np.asarray(filelist) / np.asarray(molist[0])
            rels = rels[
                ~np.isnan(rels)]  # delete nan values to get avg********
            avg = np.mean(rels)
            sd = np.std(rels)
            molTimes.append(avg)
            molStds.append(sd)
        relByFile.append(molTimes)
        sdByFile.append(molStds)
    return relByFile, sdByFile


def calc_rel_ene(minimaE):
    """
    Calculate the relative energy. For each file, take conformer energy
       relative to minimum X. The conformer minimum is chosen from
       the first conformer for which all files have an energy value.
       Note that relative energies are taken with a file's conformers,
       not subtracting one file from another (see calc_rms_error).

    Parameters
    ----------
    minimaE: 3D list of energies, where minimaE[i][j][k] represents the
      ith molecule, jth file, kth minima of that ith molecule

    Returns
    -------
    trimE: 3D list of energies as above except with relative energies
      in kcal/mol (instead of absolute energies in Hartrees). For mols
      with a single conformer it doesn't make sense to calculate rel
      energies. These mols are deleted from minimaE.
    zeroes: a 1D list of index of the reference conformer per each mol

    """

    zeroes = []
    mols2del = []
    for i, molist in enumerate(minimaE):  # loop over ith molecule
        refCount = len(molist[0])  # number of confs in ref file

        # find first conformer with least nan's.
        nanCnt = []
        for j in range(refCount):  # loop over confs of ref file
            nanCnt.append(sum(np.isnan([item[j] for item in molist])))
        print("mol {} nanCnt: {}".format(i, nanCnt))

        # get indices of confs with least nan's. no argmin bc want all idx
        where0 = np.empty(0)
        counter = 0
        while where0.size == 0 and counter < refCount:
            where0 = np.where(np.asarray(nanCnt) == counter)[0]
            counter += 1
        if where0.size > 0:  # if there ARE confs with 0 nan's, find min E
            leastNanEs = np.take(molist[0], where0)  # energies of the where0s
            winningconfIdx = where0[np.argmin(leastNanEs)]  # idx of lowest E
            zeroes.append(winningconfIdx)
        else:
            zeroes.append(nanCnt.index(min(nanCnt)))

    # calc relative energies, and convert Hartrees to kcal/mol.
    trimE = []
    for z, molE in zip(zeroes, minimaE):
        temp = []  # temp list for this mol's relative energies
        for fileE in molE:
            temp.append(
                [627.5095 * (fileE[i] - fileE[z]) for i in range(len(fileE))])
        trimE.append(temp)

    return trimE, zeroes


def write_rel_ene(molName, rmse, relEnes, zero, thryList, prefix='relene'):
    """

    """

    compF = open(prefix + '_' + molName + '.dat', 'w')
    compF.write("# Molecule %s\n" % molName)
    compF.write("# Energies (kcal/mol) for each matched conformer relative to "
                + "conformer " + str(zero) + " across each column.\n")
    compF.write("# Rows represent conformers of this molecule; columns " +
                "represent some calculations from a particular file.\n")
    compF.write("# Columns are ordered by conformer index, then the " +
                "following levels of theory:")

    # write methods, RMSEs, integer column header
    rmsheader = "\n# "
    colheader = "\n\n# "
    for i, t in enumerate(thryList):
        compF.write("\n# %d %s" % ((i + 1), t))
        rmsheader += '\t%.4f' % rmse[i]
        colheader += '\t' + str(i + 1)

    compF.write("\n\n# RMS errors by level of theory, with respect to the " +
                "first level of theory listed:")
    compF.write(rmsheader)
    compF.write(colheader)

    # write each opt's relative energies
    for i in range(len(relEnes[0])):
        compF.write('\n' + str(i) + '\t')
        thisline = [x[i] for x in relEnes]
        thisline = ['%.4f' % elem for elem in thisline]
        thisline = '\t'.join(map(str, thisline))
        compF.write(thisline)
    compF.close()


def extract_matches(moldict, *tags):
    """
    This function checks if minima is matched, using indices lists inside dict.
        If there is no match, the value (energy, time, etc) listed in the dict
        is not used, and nan is added as placeholder. After matching from
        indices, the new list is stored under a new key with the original tag
        label plus '_matched' suffix.

    TODO - remove need for tags var since tags of interest are already in dct


    Parameter
    ---------
    moldict

    Returns
    -------
    same as the parameter with new key/value pairs for matched data
    """

    for m in moldict:
        fileIndices = moldict[m]['indices']

        for t in tags:
            tagData = moldict[m][
                t]  # 2D list, 1st layer is file, 2nd layer is ene/time/etc
            updated = []

            for i, f in enumerate(
                    fileIndices):  # 1D list of indices for this file
                fileData = []

                for j, confNum in enumerate(f):
                    if confNum == None:
                        print(
                            "No matching conformer within 0.5 A RMSD for {}th\
 conf of {} mol in {}th file.".format(j, m, i))
                        fileData.append(nan)
                    elif confNum == -2:
                        # only print this warning message once per mol
                        if j == 0:
                            print("!!!! The entire {} mol is not found\
 in {}th file. !!!!".format(m, i))
                        fileData.append(nan)
                    elif len(tagData[i]) == 0:
                        print("!!!! Mol {} was found and confs were matched by\
 RMSD but there are no energies of {}th method. !!!!".format(m, i))
                        fileData.append(nan)
                    elif confNum == -1:  # -1 signifies reference theory
                        fileData.append(float(tagData[i][j]))
                    else:
                        #print(theArray[i][j])
                        fileData.append(float(tagData[i][confNum]))
                updated.append(fileData)
            moldict[m][t + '_matched'] = updated
    return moldict


### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--input",
        help="Required argument on name of text file with information on\
              file(s) and levels of theory to process.\
              See README file or examples for more details. TODO")

    parser.add_argument("--readpickle", action="store_true", default=False,
        help="If specified, read in data from pickle files from each \
              directory. Input file can be same as for heat plot inputs, \
              and pickle files will be read from same directory as \
              specified output files.")

    parser.add_argument("--verbose", action="store_true", default=False,
        help="If specified, write out relative energies in kcal/mol for \
              all conformers of all mols for all files. If in doubt, \
              do specify this option.")

    parser.add_argument("--eplot",action="store_true", default=False,
        help="Generate line plots for every molecule with relative energies.")

    parser.add_argument("--tplot", action="store_true", default=False,
        help="Generate bar plots of conformer-averaged time per each \
              optimization. One plot generated per molecule.")

    args = parser.parse_args()
    opt = vars(args)
    if not os.path.exists(opt['input']):
        raise parser.error("Input file %s does not exist." % opt['filename'])
    sys.stdout.flush()

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

    # Check that each file exists before starting -_-
    while True:
        list1 = []
        for f in sdfList:
            list1.append(os.path.isfile(f))

        if all(list1):
            # All elements are True. Therefore all the files exist. Run %run commands
            break
        else:
            raise parser.error("One or more input files are missing!!")

    # run the workhorse, unless reading in from pickle file
    if not opt['readpickle']:
        moldict = match_minima(sdfList, thryList, 'QM opt energy',
                               'opt runtime', 'opt step')
        pickle.dump(moldict, open('match.pickle', 'wb'))
    else:
        moldict = pickle.load(open('match.pickle', 'rb'))

    # process dictionary to match the data values by RMSD-matched conformers
    numMols = len(moldict)
    #moldict = extract_matches(moldict, 'QM opt energy', 'opt runtime', 'opt step')
    moldict = extract_matches(moldict, 'QM opt energy', 'opt runtime')
    minimaE = []
    for m in moldict:
        minimaE.append(moldict[m]['QM opt energy_matched'])

    # with matched energies, calculate relative values and RMS error
    trimE, zeroes = calc_rel_ene(minimaE)
    rmselist = calc_rms_error(trimE, zeroes)

    # =========================================================================

    molNames = moldict.keys()  # temp workaround since adding dictionary (TODO)
    if opt['verbose']:
        for i, mn in enumerate(molNames):
            try:
                write_rel_ene(mn, rmselist[i], trimE[i], zeroes[i], thryList)
            except IndexError:
                zeroes.append(nan)
                write_rel_ene(mn, [nan] * len(thryList), elists[i], zeroes[i],
                              thryList)

    if opt['tplot']:
        timesByMol = []
        for m in moldict:
            timesByMol.append(moldict[m]['opt runtime_matched'])
        relTimes, sdTimes = get_ratio_times(timesByMol, zeroes)

        allFileTimes = [[] for i in range(numMols)]
        allFileStds = [[] for i in range(numMols)]

        for i, molTimes in enumerate(timesByMol):
            # loop by each file to remove nans
            for j, molFileList in enumerate(molTimes):
                shortmf = np.asarray(molFileList)
                shortmf = shortmf[~np.isnan(shortmf)]
                allFileTimes[i].append(np.mean(shortmf))
                allFileStds[i].append(np.std(shortmf))

        # bar plot of average times with stdevs
        for name, fileTimes, stdevs in zip(molNames, allFileTimes,
                                           allFileStds):
            to_exclude = {}  # declare empty set if want to use all levels of thry

            fileTimes_i = [
                element for i, element in enumerate(fileTimes)
                if i not in to_exclude
            ]
            stdevs_i = [
                element for i, element in enumerate(stdevs)
                if i not in to_exclude
            ]
            thryList_i = [
                element for i, element in enumerate(thryList)
                if i not in to_exclude
            ]
            plot_avg_times(name, fileTimes_i, stdevs_i, thryList_i)

        if opt['verbose']:  # append time to relative energies file
            for i, name in enumerate(molNames):
                compF = open('relene_' + name + '.dat', 'a')
                compF.write("\n\n# Four rows: (1) avg times, (2) time stdevs,\
 (3) avg time ratios wrt ref method, (4) stdevs of time ratios:")
                avgline = stdline = a2line = s2line = "\n# "
                for j, t in enumerate(thryList):
                    avgline += ' %.4f' % allFileTimes[i][j]
                    stdline += ' %.4f' % allFileStds[i][j]
                    a2line += ' %.4f' % relTimes[i][j]
                    s2line += ' %.4f' % sdTimes[i][j]

                compF.write(avgline)
                compF.write(stdline)
                compF.write(a2line)
                compF.write(s2line)
                compF.close()

    if opt['eplot']:
        for i, m in enumerate(moldict):
#            if m != 'AlkEthOH_c1178': continue
            plot_mol_minima(m, trimE[i], thryList)
#            plot_mol_minima(m, trimE[i], thryList, selected=[0]) # zero based index
