#!/usr/bin/env python
"""
stitch_spe.py

Purpose:    Compare energies from different methods of single point energy
            calculations on the same molecule set. If a reference file is
            provided, then the RMSD in relative energies is computed for each
            comparison file. If a reference file is NOT provided, then [TODO]

Version:    Nov 16 2018
By:         Victoria T. Lim

"""

## TODO: Histogram the deviations of relative energies.
## TODO: Add line plotting functionality for specified molecules.

import os
import openeye.oechem as oechem
import numpy as np
import argparse
import proc_tags as pt
import collections  # ordered dictionary

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import operator as o


def barplot(ax, dpoints, ptitle=""):
    """
    Function modified from Peter Kerpedjiev.
    http://emptypipes.org/2013/11/09/matplotlib-multicategory-barchart/

    Create a barchart for data across different categories with
    multiple conditions for each category.

    Parameters
    ----------
    ax: The plotting axes from matplotlib.
    dpoints: The data set as an (n, 3) numpy array
    ptitle: String title for plot

    """

    # Aggregate the conditions and the categories according to their
    # mean values
    conditions = [(c, np.mean(dpoints[dpoints[:, 0] == c][:, 2].astype(float)))
                  for c in np.unique(dpoints[:, 0])]
    categories = [(c, np.mean(dpoints[dpoints[:, 1] == c][:, 2].astype(float)))
                  for c in np.unique(dpoints[:, 1])]

    # sort the conditions, categories and data so that the bars in
    # the plot will be ordered by category and condition
    conditions = [c[0] for c in sorted(conditions, key=o.itemgetter(1))]
    categories = [c[0] for c in sorted(categories, key=o.itemgetter(1))]

    dpoints = np.array(sorted(dpoints, key=lambda x: categories.index(x[1])))

    # the space between each set of bars
    space = 0.3
    n = len(conditions)
    width = (1 - space) / (len(conditions))

    # Create a set of bars at each position
    for i, cond in enumerate(conditions):
        indeces = range(len(categories))
        #indeces = range(1, len(categories)+1)
        vals = dpoints[dpoints[:, 0] == cond][:, 2].astype(np.float)
        pos = [j - (1 - space) / 2. + i * width for j in indeces]
        ax.bar(
            pos, vals, width=width, label=cond, color=cm.Accent(float(i) / n))

    # Set the x-axis tick labels to be equal to the categories
    ax.set_xticks(indeces)
    ax.set_xticklabels(categories)
    plt.setp(plt.xticks()[1], rotation=50)

    # Add the axis labels
    ax.set_ylabel("energy (kcal/mol)")
    ax.set_title(ptitle)

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
        required keys are 'ftitle' 'titleMols' 'rmsds'
    ptitle : string
        Title on plot

    """
    # fill in dictionary for ref file for list comprehension later
    wholedict[0]['titleMols'] = np.full(len(wholedict[1]['titleMols']), np.nan)
    wholedict[0]['rmsds'] = np.full(len(wholedict[1]['rmsds']), np.nan)

    # extract part of dictionary using list comprehension
    subset = np.array([[(wholedict[fd][key]) for key in\
('ftitle','titleMols','rmsds')] for fd in list(wholedict.keys())[1:]], dtype=object).T

    # build plot list
    plotlist = []
    for m in range(len(wholedict[1]['titleMols'])):
        for f in range(len(wholedict) - 1):
            temp = []
            #temp.append(subset[0][f].split('/')[-1].split('.')[0])
            temp.append(subset[0][f])
            temp.append(subset[1][f][m])
            temp.append(subset[2][f][m])
            plotlist.append(temp)
    plotlist = np.array(plotlist)

    # generate plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    barplot(ax, plotlist, ptitle)
    plt.savefig('barchart.png', bbox_inches='tight')
    plt.show()


def read_mols_tag(insdf, calctype):

    if calctype not in {'opt', 'spe'}:
        sys.exit("Specify a valid calculation type for {}.".format(insdf))

    # Open file.
    ifs1 = oechem.oemolistream()
    ifs1.SetConfTest(oechem.OEAbsoluteConfTest())
    if not ifs1.open(insdf):
        oechem.OEThrow.Fatal("Unable to open %s for reading" % insdf)
    mols = ifs1.GetOEMols()

    # Determine SD tag from which to obtain energy.
    if calctype.lower() == 'spe':
        tagword = "QM spe"
    else:
        tagword = "QM opt energy"
    return mols, tagword


def extract_enes(dict1):
    """

    From the files in input dictionaries, read in molecules, extract information
    from SD tags for conformer energies and indices.

    Parameters
    ----------
    dict1 : dict
        dictionary of input files and information to extract from SD tags
        keys are: 'ftitle' 'fname' 'calctype' 'method' 'basisset'

    Returns
    -------
    titleMols : list of strings
        names of all molecules in the SDF file
    confNums : list of ints
        conformer index numbers
    compEnes : list of numpy arrays
        conformer energies of the compared file (kcal/mol)

    """

    mols1, tag1 = read_mols_tag(dict1['fname'], dict1['calctype'])
    titleMols = []
    rmsds = []
    confNums = []
    compEnes = []
    confNans = []

    for imol in mols1:

        # Get absolute energies from the SD tags
        #print(dict2['calctype'],tag2, dict2['method'],dict2['basisset']) # for debugging
        #print(pt.get_sd_list(jmol, tag2, 'Psi4', dict2['method'],dict2['basisset'])) # for debugging
        iabs = np.array(
            list(
                map(
                    float,
                    pt.get_sd_list(imol, tag1, 'Psi4', dict1['method'],
                                 dict1['basisset']))))

        # Get omega conformer number of first, for reference info
        # whole list can be used for matching purposes
        originum = pt.get_sd_list(imol, "original index")

        # find conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(iabs))

        # convert energies from Hartrees to kcal/mol
        iabs = 627.5095 * iabs

        titleMols.append(imol.GetTitle())
        confNans.append(nanIndices)
        confNums.append(originum)
        compEnes.append(iabs)

    return titleMols, confNums, compEnes, confNans


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
            for j, (file_mols_ene, file_mols_idx) in enumerate(
                    zip(mols_ene, mols_idx)):
                enelist[i][j] = np.delete(file_mols_ene, file_mols_nan)
                idxlist[i][j] = np.delete(file_mols_idx, file_mols_nan)

    return idxlist.T, enelist.T


def relative_energies(enelist):
    # enelist[i][j] is for file i, mol j

    rel_enes = []
    for i, file_ene in enumerate(enelist):
        rel_enes.append([])
        for mol_ene in file_ene:
            rel = mol_ene - mol_ene[0]
            rel_enes[i].append(rel)
    return rel_enes


def avg_of_stdev(enelist):
    """
    For the relative energies of a single conformer, compute coefficient of
    variation (cv) across methods. Then average the CVs across all conformers
    for each molecule.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO

    Notes
    -----
    * This takes population stdevs, as opposed to sample stdevs.
    * Validated for single conformer, and for three conformers.

    """
    # enelist[i][j] is for file i, mol j
    # spread_by_mol is 1D list of len(num_mols)
    spread_by_mol = []

    # reshape list for that i goes to mol and j goes to file
    enelist = np.array(enelist).T
    for i, mols_ene in enumerate(enelist):
        sdlist = []

        # number of files is mols_ene.shape[0]
        for j in range(mols_ene[0].shape[0]):
            confs_ene = [mols_ene[k][j] for k in range(mols_ene.shape[0])]
            sd = np.std(confs_ene) / np.average(confs_ene)
            #sd = np.std(confs_ene)
            sdlist.append(sd)

        ## check that first element is 0. bc energies relative to first conf
        #if sdlist[0] != 0:
        #    print("WARNING: Something wrong with input avg_of_stdev function. "
        #          "First conformers of each mol should have zero relative "
        #          "energies and thus zero standard deviation across files.")
        #    return

        sdlist = np.delete(sdlist, 0)
        spread_by_mol.append(np.average(sdlist))

    return spread_by_mol


def norm_and_dev(enelist):
    """
    For the relative energies of each conformer of some molecule,
    subtract the average energy for the different methods, then
    take the standard deviation of all conformers of all methods.

    Parameters
    ----------
    TODO

    Returns
    -------
    TODO

    Notes
    -----
    * This takes population stdevs, as opposed to sample stdevs.
    * Validated for single conformer, and for three conformers.
    * sd and sdlist would be more accurately named cv and cvlist

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
            # note: not subtracting the average conformer energy for each method
            # bc that would give us estimate of spread conformer energies as
            # opposed to estimating spread of method energies
            confs_ene = confs_ene - np.average(confs_ene)
            mollist.append(confs_ene)

        spread_by_mol.append(np.std(mollist))

    return spread_by_mol


def extract_and_rmsd(dict1, dict2):
    """
    TODO: this can probably be rewritten or removed now to use extract_enes, remove_dead_conformers, relative_energies

    From the files in input dictionaries, read in molecules, extract information
    from SD tags, compute relative conformer energies wrt to the first conformer
    per molecule, then calculate RMSD of energies wrt reference data.

    Parameters
    ----------
    dict1 : dict
        REFERENCE dictionary from which to calculate energy RMSDs
        dictionary of input files and information to extract from SD tags
        keys are: 'ftitle' 'fname' 'calctype' 'method' 'basisset'
    dict2 : dict
        dictionary of input files and information to extract from SD tags
        keys are: 'ftitle' 'fname' 'calctype' 'method' 'basisset'

    Returns
    -------
    titleMols : list of strings
        names of all molecules in the SDF file
    rmsds : list of floats
    confNums : list of ints
        conformer index numbers, not including those with nans in either dict1/dict2 comparison
        e.g.,  if dict1 mol1 has confs of 0 1 2 3 5
              and dict2 mol1 has confs of 0 2 3 4 5
                     then output would be 0 2 3 5
    refEnes : list of numpy arrays
        relative conformer energies of the reference file (kcal/mol)
    compEnes : list of numpy arrays
        relative conformer energies of the compared file (kcal/mol)

    """

    mols1, tag1 = read_mols_tag(dict1['fname'], dict1['calctype'])
    mols2, tag2 = read_mols_tag(dict2['fname'], dict2['calctype'])
    titleMols = []
    rmsds = []
    confNums = []
    refEnes = []
    compEnes = []

    for imol in mols1:
        jmol = next(mols2)

        # check that both mols match by comparing titles and numconfs
        if (imol.NumConfs() != jmol.NumConfs()) or (imol.GetTitle() !=
                                                    jmol.GetTitle()):
            sys.exit(
                "ERROR: Either titles or number of conformers differ for mol {} in the files: "
                "{}\n{}".format(imol.GetTitle(), dict1['fname'],
                                dict2['fname']))

        # Get absolute energies from the SD tags
        #print(dict2['calctype'],tag2, dict2['method'],dict2['basisset']) # for debugging
        #print(pt.get_sd_list(jmol, tag2, 'Psi4', dict2['method'],dict2['basisset'])) # for debugging
        iabs = np.array(
            list(
                map(
                    float,
                    pt.get_sd_list(imol, tag1, 'Psi4', dict1['method'],
                                 dict1['basisset']))))
        jabs = np.array(
            list(
                map(
                    float,
                    pt.get_sd_list(jmol, tag2, 'Psi4', dict2['method'],
                                 dict2['basisset']))))

        # Get omega conformer number of first, for reference info
        # whole list can be used for matching purposes
        originum = pt.get_sd_list(imol, "original index")
        origjnum = pt.get_sd_list(jmol, "original index")

        # exclude conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(iabs))
        iabs = np.delete(iabs, nanIndices)
        jabs = np.delete(jabs, nanIndices)
        originum = np.delete(np.asarray(originum), nanIndices)
        origjnum = np.delete(np.asarray(origjnum), nanIndices)
        # same check for file2
        nanIndices = np.argwhere(np.isnan(jabs))
        iabs = np.delete(iabs, nanIndices)
        jabs = np.delete(jabs, nanIndices)
        originum = np.delete(np.asarray(originum), nanIndices)
        origjnum = np.delete(np.asarray(origjnum), nanIndices)

        # Take relative energy to first conf
        irel = iabs - iabs[0]
        jrel = jabs - jabs[0]

        # take RMSD of conformer energies for this particular mol
        dev = irel - jrel
        sqd = np.square(dev)
        mn = np.sum(sqd) / (np.shape(sqd)[0] - 1)
        rt = 627.5095 * np.sqrt(mn)

        # convert relative energies from Hartrees to kcal/mol
        irel = 627.5095 * irel
        jrel = 627.5095 * jrel

        titleMols.append(imol.GetTitle())
        rmsds.append(rt)
        confNums.append(originum)
        refEnes.append(irel)
        compEnes.append(jrel)

    return titleMols, rmsds, confNums, refEnes, compEnes


### ------------------- Script -------------------


def stitch_with_ref(wholedict, ref_index, outfn='relene-rmsd.dat'):
    """
    Compute RMSD of relative conformer energies with respect to the data
    in ref_index spot. The relative energies by conformer are computed first,
    then the RMSD is calculated with respect to reference.

    Parameters
    ----------
    wholedict : OrderedDict
        ordered dictionary of input files and information to extract from SD tags
        keys are: 'ftitle' 'fname' 'calctype' 'method' 'basisset'
    ref_index : int
        integer to specify reference file of wholedict[ref_index]
    outfn : str
        name of the output file with RMSDs of energies

    Returns
    -------
    wholedict : OrderedDict
        same as input wholedict with additional keys:
        'titleMols' 'rmsds' 'confNums' 'refEnes' 'compEnes'

    """

    sdfRef = wholedict[0]['fname']
    print("Using reference file: %s " % sdfRef)

    # Write description in output file.
    compF = open(outfn, 'w')
    compF.write("# RMSD of relative energies (kcal/mol) wrt file 0\n")

    for i, d in enumerate(wholedict.values()):
        compF.write("# File %d: %s\n" % (i, d['fname']))
        compF.write("#   calc=%s, %s/%s\n" % (d['calctype'], d['method'],
                                              d['basisset']))

    for i in range(1, len(wholedict)):
        print("Starting comparison on file: %s" % wholedict[i]['fname'])

        # each of the four returned vars is (file) list of (mols) lists
        titleMols, rmsds,confNums, refEnes, compEnes =\
            extract_and_rmsd(wholedict[ref_index], wholedict[i])
        wholedict[i]['titleMols'] = titleMols
        wholedict[i]['rmsds'] = rmsds
        wholedict[i]['confNums'] = confNums
        wholedict[i]['refEnes'] = refEnes
        wholedict[i]['compEnes'] = compEnes

    # loop over each mol and write energies from wholedict by column
    for m in range(len(wholedict[1]['titleMols'])):
        compF.write('\n\n# Mol ' + wholedict[1]['titleMols'][m])

        # for this mol, write the rmsd from each file side by side
        line = ' RMSDs: '
        for i in range(1, len(wholedict)):
            line += (str(wholedict[i]['rmsds'][m]) + '\t')
        compF.write(line)
        compF.write(
            '\n# =================================================================='
        )
        compF.write(
            '\n# conf\trel. enes. in column order by file (listed at top)')

        # for this mol, write the compEnes from each file by columns
        for c in range(len(wholedict[1]['confNums'][m])):
            line = '\n' + str(wholedict[1]['confNums'][m][c]) + '\t'
            line += str(wholedict[1]['refEnes'][m][c]) + '\t'
            for i in range(1, len(wholedict)):
                line += (str(wholedict[i]['compEnes'][m][c]) + '\t')
            compF.write(line)

    compF.close()

    return wholedict


def stitch_spe(wholedict, outfn='relene.dat'):
    """
    TODO

    Parameters
    ----------
    wholedict : OrderedDict
        ordered dictionary of input files and information to extract from SD tags
        keys are: 'ftitle' 'fname' 'calctype' 'method' 'basisset'
    outfn : str
        name of the output file with RMSDs of energies

    Returns
    -------
    wholedict : OrderedDict
        same as input wholedict with additional keys:
        'titleMols' 'stdevs' 'confNums' 'compEnes'

    """

    # Write description in output file.
    compF = open(outfn, 'w')
    compF.write(
        "# Spread of relative energies (kcal/mol) over diff QM methods averaged over conformers\n"
    )
    for i, d in enumerate(wholedict.values()):
        compF.write("# File %d: %s\n" % (i, d['fname']))
        compF.write("#   calc=%s, %s/%s\n" % (d['calctype'], d['method'],
                                              d['basisset']))

    enelist = []
    idxlist = []
    nanlist = []
    for i in range(len(wholedict)):

        print("Extracting data for file: %s" % wholedict[i]['fname'])
        titleMols, confNums, compEnes, confNans = extract_enes(wholedict[i])
        enelist.append(compEnes)
        idxlist.append(confNums)
        nanlist.append(confNans)
        wholedict[i]['titleMols'] = titleMols

    print("Removing un-finished conformers and computing relative energies...")
    idxlist, enelist = remove_dead_conformers(enelist, idxlist, nanlist)
    relenelist = relative_energies(enelist)
    for i in range(len(wholedict)):
        wholedict[i]['compEnes'] = relenelist[i]
        wholedict[i]['confNums'] = idxlist[i]

    # estimate spread of data
    spreadlist = avg_of_stdev(relenelist)
    #spreadlist = norm_and_dev(relenelist)

    # loop over each mol and write energies from wholedict by column
    for m in range(len(wholedict[1]['titleMols'])):
        compF.write('\n\n# Mol ' + wholedict[1]['titleMols'][m])

        # for this mol, write the rmsd from each file side by side
        line = ' Spread: {}'.format(spreadlist[m])
        compF.write(line)
        compF.write(
            '\n# =================================================================='
        )
        compF.write(
            '\n# conf\trel. enes. in column order by file (listed at top)')

        # for this mol, write the compEnes from each file by columns
        for c in range(len(wholedict[1]['confNums'][m])):
            line = '\n' + str(wholedict[1]['confNums'][m][c]) + '\t'
            for i in range(len(wholedict)):
                line += (str(wholedict[i]['compEnes'][m][c]) + '\t')
            compF.write(line)

    compF.close()

    return wholedict


### ------------------- Parser -------------------

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", required=True,
        help="Name of text file with information on file(s) and SD tag info on "
             "the level of theory to process.")

    parser.add_argument("--reffile", default=None,
        help="If specified, treat the specified file as a reference file from which "
             "to compute RMSD of energies. If not, write out avgs/stdevs for "
             "each conformer's relative energies in kcal/mol.")

    #    parser.add_argument("-l", "--lineplot",
    #        help=("Optional, list of molecule name(s) for which line plots will be "
    #              "generated. Each molecule's line plots shows relative energies "
    #              "for each file.")

    parser.add_argument("--plotbars", action="store_true", default=False,
        help="If specified, generate bar plots of each comparison file wrt "
             "reference file (first entry of input file).")

    args = parser.parse_args()
    if not os.path.exists(args.infile):
        raise parser.error("Input file %s does not exist." % args.infile)

    # Read input file into an ordered dictionary.
    # http://stackoverflow.com/questions/25924244/creating-2d-dictionary-in-python
    linecount = 0
    wholedict = collections.OrderedDict()
    with open(args.infile) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if args.reffile is not None and args.reffile in line:
                ref_index = linecount
            dataline = [x.strip() for x in line.split(',')]
            wholedict[linecount] = {
                'ftitle': dataline[0],
                'fname': dataline[1],
                'calctype': dataline[2],
                'method': dataline[3],
                'basisset': dataline[4]
            }
            linecount += 1

    if args.reffile is not None:
        wholedict = stitch_with_ref(wholedict, ref_index)
        if args.plotbars:
            arrange_and_plot(wholedict, "RMSDs of relative conformer energies")
    else:
        wholedict = stitch_spe(wholedict)
        if args.plotbars:
            print('still working on this')
        #    arrange_and_plot(wholedict, "RMSDs of relative conformer energies")
