#!/usr/bin/env python
"""
survey_conformers.py

Purpose:    Compare energies and calculation times of the same molecule set
            across different QM methods.

Version:    Apr 24 2019
By:         Victoria T. Lim

"""

import os
import numpy as np
import operator as o
from scipy import stats
import matplotlib.pyplot as plt

import quanformer.proc_tags as pt
import quanformer.reader as reader



def avg_mol_time(titles, infile, method, basis, tag, mol_slice=[]):
    """
    For an SDF file with all confs of all mols, get the average runtime
        of all conformers for each molecule.
    The input dictionary may or not be empty. If it is, append the avg/stdev
        time of the calculation from this infile to existing molecule's value.

    Parameters
    ----------
    titles : dictionary
        Keys are molecule names.
        Values are [[qm1_avg, qm1_std], [qm2_avg, qm2_std], ... ]
        where the index refers to a particular level of theory.
        Dictionary may or may not be empty.
    infile : string
        name of the SDF file from which to extract time data from SD tag
    method : string
        QM method
    basis : string
        QM basis set
    tag : string
        datum of interest, e.g., "QM opt energy"
        See keys in the define_tag function of proc_tags module.
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]

    Returns
    -------
    titles : dictionary
        dictionary with extracted data from SDF file; keys are molnames,
        values are lists of list of [avg_time, stdev_time] for many QM methods

    """

    # Open molecule file.
    if len(mol_slice) == 3:
        mols = reader.read_mols(infile, mol_slice)
    else:
        mols = reader.read_mols(infile)


    # Prepare text file to write extracted data.
    timeF = open("timeAvgs.txt", 'a')
    timeF.write("\nFile: {}\n".format(infile))
    timeF.write("Average [{}/{}] [{}s] over all confs for each molecule\n".format(method, basis, tag))

    for mol_i in mols:

        # get array of all conformer data of this mol
        try:
            time_array = np.fromiter(pt.get_sd_list(mol_i, tag, 'Psi4', method, basis), dtype=np.float64)
        except ValueError:
            time_array = np.asarray([np.nan])*mol_i.NumConfs()

        # exclude conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(time_array))
        for i in reversed(nanIndices):  # loop in reverse to delete correctly
            time_array = np.delete(time_array, i)
        meantime = np.mean(time_array)
        stdtime = np.std(time_array)

        # write out data to file and store in dictionary
        timeF.write("  %s\t%d confs\t\t%.3f +- %.3f\n" % (mol_i.GetTitle(), time_array.size, meantime, stdtime))
        name = mol_i.GetTitle()
        if name not in titles:
            titles[name] = []
        titles[name].append([meantime, stdtime])

    timeF.close()
    return titles


def plot_groupedbar(ax, labels_data_std, errbar=False):
    """
    Function modified from Peter Kerpedjiev.
    http://emptypipes.org/2013/11/09/matplotlib-multicategory-barchart/

    Create a barchart for data across different x_labels with
    multiple color_labels for each x_label.

    Parameters
    ----------
    ax: The plotting axes from matplotlib.
    labels_data_std: The data set as an (n, 3) numpy array

    """
    # Aggregate the color_labels and the x_labels according to their mean values
    color_labels = [(c, np.mean(labels_data_std[labels_data_std[:, 0] == c][:, 2].astype(float))) for c in np.unique(labels_data_std[:, 0])]
    x_labels = [(c, np.mean(labels_data_std[labels_data_std[:, 1] == c][:, 2].astype(float))) for c in np.unique(labels_data_std[:, 1])]

    # sort the color_labels, x_labels and data so that the bars in
    # the plot will be ordered by x_category and color_label
    color_labels = [c[0] for c in sorted(color_labels, key=o.itemgetter(1))]
    x_labels = [c[0] for c in sorted(x_labels, key=o.itemgetter(1))]
    labels_data_std = np.array(sorted(labels_data_std, key=lambda x: x_labels.index(x[1])))

    # the space between each set of bars
    space = 0.3
    n = len(color_labels)
    width = (1 - space) / (len(color_labels))
    indices = range(len(x_labels))

    # set color cycle
    if n < 9:
        ax.set_prop_cycle(plt.cycler('color', plt.cm.Accent(np.linspace(0, 1, n))))
    elif n < 11:
        ax.set_prop_cycle(plt.cycler('color', plt.cm.tab10(np.linspace(0, 1, n))))
    elif n < 13:
        ax.set_prop_cycle(plt.cycler('color', plt.cm.Paired(np.linspace(0, 1, n))))
    else:
        ax.set_prop_cycle(plt.cycler('color', plt.cm.tab20(np.linspace(0, 1, n))))

    # Create a set of bars at each position
    for i, cond in enumerate(color_labels):
        vals = labels_data_std[labels_data_std[:, 0] == cond][:, 2].astype(np.float)
        pos = [j - (1 - space) / 2. + i * width for j in indices]
        if errbar:
            stds = labels_data_std[labels_data_std[:, 0] == cond][:, 3].astype(np.float)
            ax.bar(pos, vals, yerr=stds, width=width, label=cond)
        else:
            ax.bar(pos, vals, width=width, label=cond)

    # Set the x-axis tick labels to be equal to the x_labels
    if len(indices) > 7: # adjust int as needed for better placement of many labels
        indices = [x - 1 for x in indices]
    ax.set_xticks(indices)
    ax.set_xticklabels(x_labels)
    plt.setp(plt.xticks()[1], rotation=50)

    # Add a legend
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper left')


def arrange_and_plot(wholedict, ptitle):
    """
    Organize and plot data.

    Parameters
    ----------
    wholedict : OrderedDict
        dictionary of information to be plotted, ordered by input file
        required keys are 'label' 'titleMols' 'rmsds'
    ptitle : string
        Title on plot

    """
    # fill in dictionary for ref file for list comprehension later
    wholedict[0]['titleMols'] = np.full(len(wholedict[1]['titleMols']), np.nan)
    wholedict[0]['rmsds'] = np.full(len(wholedict[1]['rmsds']), np.nan)

    # extract part of dictionary using list comprehension
    subset = np.array([[(wholedict[fd][key]) for key in\
('label','titleMols','rmsds')] for fd in list(wholedict.keys())[1:]], dtype=object).T

    # build plot list
    plotlist = []
    for m in range(len(wholedict[1]['titleMols'])):
        for f in range(len(wholedict) - 1):
            temp = []
            temp.append(subset[0][f])
            temp.append(subset[1][f][m])
            temp.append(subset[2][f][m])
            plotlist.append(temp)
    plotlist = np.array(plotlist)

    # generate plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_groupedbar(ax, plotlist)
    plt.title(ptitle, fontsize=14)
    plt.savefig('barchart.png', bbox_inches='tight')
    plt.show()


def extract_enes(dict1, mol_slice=[]):
    """
    From files in input dictionaries, read in molecules, extract information
    from SD tags for conformer energies and indices.

    Parameters
    ----------
    dict1 : dict
        dictionary of input files and information to extract from SD tags
        keys are: 'theory' 'fname' 'tagkey' 'label'
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]

    Returns
    -------
    titleMols : list of strings
        names of all molecules in the SDF file
    confNums : list of ints
        conformer index numbers
    enes : list of numpy arrays
        conformer energies of the compared file (kcal/mol)
    confNans : list of numpy arrays
        indices of enes where the values are nans

    """

    # Open molecule file.
    if len(mol_slice) == 3:
        mols = reader.read_mols(dict1['fname'], mol_slice)
    else:
        mols = reader.read_mols(dict1['fname'])

    short_tag = dict1['tagkey']
    qmethod, qbasis = reader.separated_theory(dict1['theory'])

    titleMols = []
    confNums = []
    enes = []
    confNans = []

    for imol in mols:

        # Get absolute energies from the SD tags
        iabs = np.array(list(map(float, pt.get_sd_list(imol, short_tag, 'Psi4', qmethod, qbasis))))

        # Get omega conformer number of first, for reference info
        # whole list can be used for matching purposes
        indices_orig = pt.get_sd_list(imol, "original index")

        # find conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(iabs))

        # convert energies from Hartrees to kcal/mol
        iabs = 627.5095 * iabs

        titleMols.append(imol.GetTitle())
        confNans.append(nanIndices)
        confNums.append(indices_orig)
        enes.append(iabs)

    return titleMols, confNums, enes, confNans


def remove_dead_conformers(enelist, idxlist, nanlist):
    # enelist[i][j] is for file i, mol j

    # reshape list for that i goes to mol and j goes to file
    enelist = np.array(enelist).T
    idxlist = np.array(idxlist).T
    nanlist = np.array(nanlist).T

    for i, mols_nan in enumerate(nanlist):
        for file_mols_nan in mols_nan:
            mols_ene = enelist[i]
            mols_idx = idxlist[i]
            for j, (file_mols_ene, file_mols_idx) in enumerate(zip(mols_ene, mols_idx)):
                enelist[i][j] = np.delete(file_mols_ene, file_mols_nan)
                idxlist[i][j] = np.delete(file_mols_idx, file_mols_nan)

    # transpose it back to return
    return idxlist.T, enelist.T


def relative_energies(enelist, method, ref_conf_int=0):
    """
    Compute relative energies by either subtracting or dividing
        by a reference conformer's energy value. This energy value
        is different for every mol in every file.

    Parameters
    ----------
    enelist : list of lists
        enelist[i][j] is the energy for file i, mol j
    method : string
        either 'subtract' or 'divide' to take all values of mol
        and subtract by first conf or divide by first conf
    ref_conf_int : integer
        0-based index of the integer for which to use as reference
        for taking relative energies

    Returns
    --------
    rel_enes : list of lists
        same structure as input enelist but with relative energies

    """

    rel_enes = []
    for i, file_ene in enumerate(enelist):
        rel_enes.append([])
        for mol_ene in file_ene:
            if method == 'subtract':
                rel = mol_ene - mol_ene[ref_conf_int]
            elif method == 'divide':
                rel = mol_ene/mol_ene[ref_conf_int]
            rel_enes[i].append(rel)
    return rel_enes


def avg_coeffvar(enelist, ref_conf_int=0):
    """
    For the ABSOLUTE energies of a single conformer, compute
        coefficient of variation (CV) across methods, then
        average the CVs across all conformers for each molecule.

    Note: Input energies can be raw energies, or scaled energies (such that
        energies are still on ratio scale).
    Note: It is incorrect to compute the coefficient of variation on relative
        values. See wikipedia for an example on temperature.
    Note: The units of the input shouldn't matter since the CV is a ratio,
        meaning that units should cancel.

    Parameters
    ----------
    enelist : list of lists
        enelist[i][j] is absolute energy of file i, mol j

    Returns
    -------
    spread_by_mol : list
        len(spread_by_mol) == number of molecules
        spread_by_mol[i] == avg(CVs(conformer enes from diff methods) of mol i
    cvlist_all_mols : list of lists
        len(cvlist_all_mols) == number of molecules
        len(cvlist_all_mols[i]) == number of conformers for molecule i minus 1

    Notes
    -----
    * This takes population stdevs, as opposed to sample stdevs.
    * Validated for single conformer, and for three conformers.

    """

    # 1D list for average of CVs of each mol
    spread_by_mol = []
    cvlist_all_mols = []

    # reshape list so that enelist[i][j] is mol i, file j
    enelist = np.array(enelist).T

    # here, mols_ene[x][y] is ith molecule, xth file, and yth conformer
    for i, mols_ene in enumerate(enelist):
        cvlist = []

        # check that at least two conformers exist
        num_confs = mols_ene[0].shape[0]
        if num_confs <= 1:
            print(f"skipping molecule {i} since only 0 or 1 conformer")
            spread_by_mol.append(np.nan)
            cvlist_all_mols.append(num_confs*[np.nan])
            continue

        # number of files is mols_ene.shape[0]; should be same for all mols
        for j in range(num_confs):

            # for single conformer j, get its energies from all k files
            confs_ene = [mols_ene[k][j] for k in range(mols_ene.shape[0])]

            # compute coeff of variation
            with np.errstate(divide='ignore', invalid='ignore'):
                cv = np.std(confs_ene) / np.average(confs_ene)
            cvlist.append(cv)

        # remove ref conf from avg calc bc energies scaled by this conf
        cvlist.pop(ref_conf_int)

        # average the coeffs of variation
        spread_by_mol.append(np.average(cvlist))
        cvlist_all_mols.append(cvlist)

    return spread_by_mol, cvlist_all_mols


def normalized_deviation(enelist):
    """
    For the relative energies of each conformer of some molecule,
    subtract the average energy for the different methods, then
    take the standard deviation of all conformers of all methods.

    TODO -- should the stdev be computed for all the conformers
        of the mol (currently coded) or across the methods
        (so add the line under np.average)? Would the stdevs of
        each conformer then be averaged?
    AKA should you be taking the stdev of the normalized data (current)
        or take avg of the stdev?

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO

    """

    # enelist[i][j] is for file i, mol j
    # spread_by_mol is 1D list of len(num_mols)
    spread_by_mol = []

    # reshape list for that i goes to mol and j goes to file
    enelist = np.array(enelist).T
    for i, mols_ene in enumerate(enelist):
        mollist = []

        # number of files is mols_ene.shape[0]
        for j in range(1, mols_ene[0].shape[0]):
            confs_ene = [mols_ene[k][j] for k in range(mols_ene.shape[0])]

            # subtract the average of methods for the conformer
            # note: not subtracting average conformer energy for each method
            # bc that would give us estimate of spread conformer energies as
            # opposed to estimating spread of method energies
            confs_ene = confs_ene - np.average(confs_ene)
            mollist.append(confs_ene)

        spread_by_mol.append(np.std(mollist))

    return spread_by_mol


def rmsd_two_files(dict1, dict2, mol_slice=[]):
    """
    From files in input dictionaries, read in molecules, extract information
    from SD tags, compute relative conformer energies wrt to first conformer
    of each molecule, and calculate RMSD of energies wrt reference data.

    Parameters
    ----------
    dict1 : dict
        REFERENCE dictionary from which to calculate energy RMSDs
        dictionary of input files and information to extract from SD tags
        keys are: 'theory' 'fname' 'tagkey' 'label'
    dict2 : dict
        dictionary of input files and information to extract from SD tags
        keys are: 'theory' 'fname' 'tagkey' 'label'
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]

    Returns
    -------
    titleMols : list of strings
        names of all molecules in the SDF file
    rmsds : list of floats
    confNums : list of ints
        conformer index numbers, excluding nans of both dict1 and dict2
        e.g.,  if dict1 mol1 has confs of 0 1 2 3 5
              and dict2 mol1 has confs of 0 2 3 4 5
                     then output would be 0 2 3 5
    reference_enes : list of numpy arrays
        relative conformer energies of the reference file (kcal/mol)
    compared_enes : list of numpy arrays
        relative conformer energies of the compared file (kcal/mol)

    """

    titleMols = []
    rmsds = []
    confNums = []
    reference_enes = []
    compared_enes = []

    # load molecules and extract data
    mols1_titles, mols1_indices, mols1_enes, mols1_nans = extract_enes(dict1, mol_slice)
    mols2_titles, mols2_indices, mols2_enes, mols2_nans = extract_enes(dict2, mol_slice)

    # check that number of molecules exist in both files
    if len(mols1_titles) != len(mols2_titles):
        raise AssertionError("The number of molecules do not match for files:"
                             "\n{}\n{}".format(dict1['fname'], dict2['fname']))

    for i in range(len(mols1_titles)):
        title1 = mols1_titles[i]
        title2 = mols2_titles[i]
        indices1_orig = mols1_indices[i]
        indices2_orig = mols2_indices[i]
        enes1 = mols1_enes[i]
        enes2 = mols2_enes[i]
        nans1 = mols1_nans[i]
        nans2 = mols2_nans[i]

        # check that the names and conformers match for both mols
        if title1 != title2 or (indices1_orig != indices2_orig):
            raise AssertionError("ERROR: Either names or number of conformers "
                                 "differ for mol {}/{} in:\n{}\n{}".format(
                                 title1, title2, dict1['fname'], dict2['fname']))

        # exclude conformers for which job did not finish (nan)
        indices1_prune1, enes1_prune1 = remove_dead_conformers(enes1, indices1_orig, nans1)
        indices2_prune1, enes2_prune1 = remove_dead_conformers(enes2, indices2_orig, nans1)
        indices1_prune2, enes1_prune2 = remove_dead_conformers(enes1_prune1, indices1_prune1, nans2)
        indices2_prune2, enes2_prune2 = remove_dead_conformers(enes2_prune1, indices2_prune1, nans2)

        # make sure conformer indices match after pruning
        if not (indices1_prune2 == indices2_prune2).all():
            raise AssertionError("List of conformer indices do not match "
                                 "after pruning conformers for {}: in files:\n{}\n{}\n".format(
                                     title1, dict1['fname'], dict2['fname']))

        # check that there exists at least two conformers
        num_confs = indices1_prune2.shape[0]
        if num_confs <= 1:
            print("skipping molecule {} since only 0 or 1 conformer".format(title1))
            continue

        # take relative energy to first conf
        irel = enes1 - enes1[0]
        jrel = enes2 - enes2[0]

        # take RMSD of conformer energies for this particular mol
        dev = irel - jrel
        sqd = np.square(dev)
        mn = np.sum(sqd) / (np.shape(sqd)[0] - 1)
        rt = np.sqrt(mn)

        # store data
        titleMols.append(title1)
        rmsds.append(rt)
        confNums.append(indices1_prune2)
        reference_enes.append(irel)
        compared_enes.append(jrel)

    return titleMols, rmsds, confNums, reference_enes, compared_enes


### ------------------- Script -------------------


def survey_energies_ref(wholedict, ref_index, mol_slice=[], outfn='relene-rmsd.dat'):
    """
    Compute RMSD of relative conformer energies with respect to the data
    in ref_index spot. The relative energies by conformer are computed first,
    then the RMSD is calculated with respect to reference.

    Parameters
    ----------
    wholedict : OrderedDict
        ordered dictionary of input files and information to extract from SD tags
        keys are: 'theory' 'fname' 'tagkey' 'label'
    ref_index : int
        integer to specify reference file of wholedict[ref_index]
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]
    outfn : str
        name of the output file with RMSDs of energies

    Returns
    -------
    wholedict : OrderedDict
        same as input wholedict with additional keys:
        'titleMols' 'rmsds' 'confNums' 'reference_enes' 'compared_enes'

    """

    description = 'relative energies (kcal/mol)'
    infile = wholedict[0]['fname']
    print("Using reference file: %s " % infile)

    # Write description in output file.
    compF = open(outfn, 'w')
    compF.write(f"# RMSD of {description} wrt file 0\n")

    for i, d in enumerate(wholedict.values()):
        compF.write("# File %d: %s\n" % (i, d['fname']))
        compF.write("#   calc=%s, %s\n" % (d['tagkey'], d['theory']))

    for i in range(1, len(wholedict)):
        print("\nStarting comparison on file: %s" % wholedict[i]['fname'])

        # each of the four returned vars is (file) list of (mols) lists
        titleMols, rmsds, confNums, reference_enes, compared_enes =\
            rmsd_two_files(wholedict[ref_index], wholedict[i], mol_slice)
        wholedict[i]['titleMols'] = titleMols
        wholedict[i]['rmsds'] = rmsds
        wholedict[i]['confNums'] = confNums
        wholedict[i]['reference_enes'] = reference_enes
        wholedict[i]['compared_enes'] = compared_enes

    # loop over each mol and write energies from wholedict by column
    for m in range(len(wholedict[1]['titleMols'])):
        compF.write('\n\n# Mol ' + wholedict[1]['titleMols'][m])

        # for this mol, write the rmsd from each file side by side
        line = ' RMSDs: '
        for i in range(1, len(wholedict)):
            line += (str(wholedict[i]['rmsds'][m]) + '\t')
        compF.write(line)
        compF.write('\n# ==================================================================')
        compF.write(f'\n# conf\t{description} in column order by file (listed at top)')

        # for this mol, write the compared_enes from each file by columns
        for c in range(len(wholedict[1]['confNums'][m])):
            line = '\n' + str(wholedict[1]['confNums'][m][c]) + '\t'
            line += str(wholedict[1]['reference_enes'][m][c]) + '\t'
            for i in range(1, len(wholedict)):
                line += (str(wholedict[i]['compared_enes'][m][c]) + '\t')
            compF.write(line)

    compF.close()

    return wholedict


def survey_energies(wholedict, mol_slice=[], outfn='relene.dat', ref_conf=0):
    """
    Compute the spread of the conformer energies for each molecule in each file.
    The spread can be computed as the coefficient of variation of all methods,
    computed for each conformer, which is averaged over all conformers. Another
    approach to computing spread in this function is simply normalizing the
    energies by remove the average, then taking the stdev of the normalized data.
     ^ TODO reconsider this, and add variable in def if wanting to include
       TODO also add tests for these spread functions

    Conformers that are missing for one calculation/file are taken out of
    the analysis for all other calculations/files.

    Parameters
    ----------
    wholedict : OrderedDict
        ordered dictionary of input files and information to extract from SD tags
        keys are: 'theory' 'fname' 'tagkey' 'label'
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]
    outfn : str
        name of the output file with RMSDs of energies

    Returns
    -------
    wholedict : OrderedDict
        same as input wholedict with additional keys:
        'titleMols' 'confNums' 'compared_enes'

    """

    num_files = len(wholedict)

    # Write description in output file.
    compF = open(outfn, 'w')
    compF.write(f"# Variability in energies from diff QM methods averaged over conformers\n")
    for i, d in enumerate(wholedict.values()):
        compF.write("# File %d: %s\n" % (i+1, d['fname']))
        compF.write("#   calc=%s, %s\n" % (d['tagkey'], d['theory']))

    enelist = []
    idxlist = []
    nanlist = []
    for i in range(num_files):
        print("Extracting data for file: %s" % wholedict[i]['fname'])
        titleMols, confNums, compared_enes, confNans = extract_enes(wholedict[i], mol_slice)
        enelist.append(compared_enes)
        idxlist.append(confNums)
        nanlist.append(confNans)
        wholedict[i]['titleMols'] = titleMols

    print("Removing mols missing from 1+ files...")
    # find molecules common in all files
    molname_sets = [set(wholedict[x]['titleMols']) for x in range(num_files)]
    molname_common = set.intersection(*molname_sets)

    for i in range(num_files):
        # find the extraneous mols in this file
        mols_in_file = wholedict[i]['titleMols']
        mols_extraneous = list(set(mols_in_file) - molname_common)

        # get indices of mols to be removed for this file
        removal_idx = []
        for molname in mols_extraneous:
            removal_idx.append(mols_in_file.index(molname))
        for index in sorted(removal_idx, reverse=True):
            print(f"  Removing {wholedict[i]['titleMols'][index]} from {wholedict[i]['fname']}")
            del wholedict[i]['titleMols'][index]
            del enelist[i][index]
            del idxlist[i][index]
            del nanlist[i][index]

    print("Removing un-finished conformers and computing relative energies...")
    idxlist, enelist = remove_dead_conformers(enelist, idxlist, nanlist)

    # scale energies: scale energies by subtracting or dividing conf_j
    # by conf_i for all j confs in all mols for all files
    # (1) subtract first conf, or (2) divide first conf
    if isinstance(ref_conf, int):
        ref_conf_int = ref_conf
    rel_enelist = relative_energies(enelist, 'subtract', ref_conf_int)
    rat_enelist = relative_energies(enelist, 'divide', ref_conf_int)

    for i in range(num_files):
        wholedict[i]['compared_enes'] = rel_enelist[i]
        wholedict[i]['confNums'] = idxlist[i]

    # check that entire array is not zero, if ALL mols have one conf
    all_zeros = np.all(np.asarray(rel_enelist)==0)
    if all_zeros:
        print("WARNING: All relative energy values are zero. "
            "Cannot calculate energy spread.")
        spreadlist = [-1]*len(rel_enelist[0])
    else:
        # estimate spread of data on energies
        spreadlist, cvlist_all_mols = avg_coeffvar(rat_enelist, ref_conf_int)
        #spreadlist = normalized_deviation(rel_enelist)

    # loop over each mol and write energies from wholedict by column
    for m in range(len(wholedict[1]['titleMols'])):
        compF.write('\n\n# Mol ' + wholedict[1]['titleMols'][m])

        # for this mol, write the rmsd from each file side by side
        line = ' {:.2f} ppm averaged CV'.format(spreadlist[m]*1000000)
        compF.write(line)
        compF.write('\n# ==================================================================')
        compF.write(f'\n# rel enes (kcal/mol), column i = file i, row j = conformer j')

        # for this mol, write the compared_enes from each file by columns
        this_mol_conf_inds = wholedict[1]['confNums'][m]
        for c in range(len(this_mol_conf_inds)):
            line = '\n' + str(this_mol_conf_inds[c]) + '\t'
            for i in range(num_files):
                line += '{:.4f}\t'.format(wholedict[i]['compared_enes'][m][c])

            # print overall contribution to spread, skip ref conf since scaled
            # include second condition in case ref_conf_int is negative
            if c == ref_conf_int or c == len(this_mol_conf_inds)+ref_conf_int:
                line += '#  nan ppm CV'
            else:
                line += '# {:.2f} ppm CV'.format(cvlist_all_mols[m][c-1]*1000000)
            compF.write(line)

    compF.close()

    return wholedict


def survey_times(wholedict, mol_slice=[]):
    """
    Average all conformers' calculations times for each molecule in each file.
    Generate grouped bar plot where x=QM method, y=time(s). Molecules are
    grouped together by method and color coded by molecule identity.

    Parameters
    ----------
    wholedict : OrderedDict
        ordered dictionary of input files and information to extract from SD tags
        keys are: 'theory' 'fname' 'tagkey' 'label'
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]

    """

    # create storage space
    titles = {}
    timeplot = []  # timeplot[i][j] == avg of conformer times for mol i, file j
    stdplot = []  # stdplot[i][j] == stdev of conformer times for mol i, file j

    # loop over each files and compute avg conformer times for all mols
    for i in wholedict:
        qmethod, qbasis = reader.separated_theory(wholedict[i]['theory'])
        qtag = wholedict[i]['tagkey']
        titles = avg_mol_time(titles, wholedict[i]['fname'], qmethod, qbasis, qtag, mol_slice)

    # extract data in dictionary to arrays
    for mol in titles.values():
        timeplot.append([k[0] for k in mol])
        stdplot.append([k[1] for k in mol])

    # calculate number of datafiles for each mol
    lens = [len(x) for x in timeplot]

    # obtain mode representing most common number of files
    # first [0] to get modes not counts, second [0] to get highest mode
    m = stats.mode(lens)[0][0]

    # delete sublists (mols) with missing data; asarray takes uniform length sublists
    # should only happen if 1+ method fails job for ALL confs of some mol
    # TODO: this might be addressed in a better way
    tracker = []
    for i, t in enumerate(timeplot):
        if len(t) != m:
            tracker.append(i)
    for index in sorted(tracker, reverse=True):
        mol_name = list(titles.keys())[index]
        print(f"removing {mol_name} from time analysis")
        del timeplot[index]
        del stdplot[index]
        del titles[mol_name]

    # prepare plot data
    timeplot = np.asarray(timeplot)
    stdplot = np.asarray(stdplot)
    timearray = timeplot.flatten('F')
    stdarray = stdplot.flatten('F')

    # generate mol labels (colors) from dict keys, repeated by number of files
    # e.g., ['a', 'b', 'c', 'a', 'b', 'c']
    mol_names = list(titles.keys()) * len(wholedict)

    # generate method labels (x-axis) from dict value of 'theory', repeated
    # by number of mols. e.g., ['a', 'a', 'b', 'b', 'c', 'c']
    # - remove duplicated spaces from method listing in input file
    # - get the first word in tagkey which should be 'opt' or 'spe'
    # https://stackoverflow.com/questions/2449077/duplicate-each-member-in-a-list-python
    method_labels_short = [
        " ".join(wholedict[item]['theory'].split()) + ' ' +
        wholedict[item]['tagkey'].split(" ")[0] for item in wholedict]
    method_labels = [val for val in method_labels_short for _ in range(len(titles))]

    # combine color labels, x labels, and plot data
    plotlist = np.column_stack((mol_names, method_labels,
        timearray.astype(np.object), stdarray.astype(np.object)))
    print(plotlist)

    # generate plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plot_groupedbar(ax, plotlist)

    plt.title("conformer-averaged wall-clock times", fontsize=16)
    plt.ylabel("time (s)", fontsize=14)
    plt.yticks(fontsize=12)

    plt.savefig('rename_me.png', bbox_inches='tight')
    plt.show()


def survey_confs(infile, analyze_energies=False, analyze_times=False, ref_index=None, plot_enes=True, mol_slice=[]):
    """
    Main interface to survey_times, survey_energies, or survey_energies_ref.

    Parameters
    ----------
    infile : string
        Name of text file with information on file(s) and SD tag info on
        the level of theory to process.
    analyze_energies : Boolean
        analyze energies for the provided input files
    analyze_times : Boolean
        analyze calculation times for the provided input files
    ref_index : integer
        provide the integer of the line number in the input file to treat
        as reference file (such as to compute RMSD of energies);
        line numbers start at ZERO and do not count commented lines
    plot_enes : Boolean
        generate plot for energies analysis
    mol_slice : list
        list of indices from which to slice mols generator for read_mols
        [start, stop, step]

    """

    if not os.path.exists(infile):
        raise FileNotFoundError(f"Input file {infile} does not exist.")

    wholedict = reader.read_text_input(infile)

    if (analyze_energies == analyze_times):
        raise ValueError("Please specify ONE of analyze_energies or analyze_times")

    if analyze_times:
        survey_times(wholedict, mol_slice)

    if analyze_energies:
        if ref_index is None:
            survey_energies(wholedict, mol_slice)
        else:
            wholedict = survey_energies_ref(wholedict, ref_index, mol_slice)
            # TODO move this part inside the survey_energies fx
            if plot_enes:
                arrange_and_plot(wholedict, "RMSDs of relative conformer energies")
