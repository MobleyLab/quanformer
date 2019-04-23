#!/usr/bin/env python
"""
basic_plot.py

Basic plots (lines or bars) of information in SD tags.

Version:    Apr 12 2019
By:         Victoria T. Lim


TODO:
- plot diff tags in same file

"""
import numpy as np
import itertools

import matplotlib.pyplot as plt

import quanformer.proc_tags as pt
import quanformer.reader as reader


def basic_plot(infile, tag, style, molname=None, take_relative=False, har_to_kcal=False):
    """
    TODO

    Parameters
    ----------
    infile : string
        Name of SDF file with information in SD tags.
    tag : string
        Full tag string directly as listed in the SD file.
    style : string
        plot style. can be 'scatter', 'line', or 'bar'
        TODO
    take_relative : Boolean
        subtract lowest value
    har_to_kcal : Boolean
        multiply data in Hartrees by 627.5095 to yield kcal/mol

    """
    # Open molecule file.
    mols = reader.read_mols(infile)

    for i, mol_i in enumerate(mols):
        if molname is not None and mol_i.GetTitle() != molname:
            continue

        # get array of all conformer data of this mol
        try:
            data_array = np.fromiter(pt.get_sd_list(mol_i, datum='', taglabel=tag), dtype=np.float64)
        except ValueError:
            data_array = np.asarray([np.nan])*mol_i.NumConfs()

        # exclude conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(data_array))
        for i in reversed(nanIndices):  # loop in reverse to delete correctly
            data_array = np.delete(data_array, i)

        if take_relative:
            data_array = data_array - np.amin(data_array)
        if har_to_kcal:
            data_array = 627.5095*data_array

        # generate plot
        plt.plot(data_array)
        plt.grid()
        plt.title(mol_i.GetTitle()+'\n'+tag, fontsize=14)
        plt.savefig(f'output_{i}.png', bbox_inches='tight')
        plt.show()


def combine_files_plot(infile, figname='combined.png', molname=None, verbose=False, take_relative=False, har_to_kcal=False):
    """
    TODO

    This only supports plotting of ONE specified molecule across different files.

    Note on take_relative:
        [1] Subtracting global minimum (single value) from all energies
        doesn't work since everything is still on different scale.
    subtract: (1) first conformer of each?, (2) global minimum?, (3) minimum of each?

    Parameters
    ----------
    infile : str
        Filename with information on the files to read in, and
        the SDF tags to be extracted from each. Columns are:
        (1) QM method/basis, (2) sdf file, (3) tag key in sdf (like 'QM spe'),
        (4) arbitrary label for plotting. Separate columns by comma.
    molname
    verbose

    """
    wholedict = reader.read_text_input(infile)

    numFiles = len(wholedict)
    xarray = []
    yarray = []
    labels = []
    titles = []
    for i in wholedict:
        print("Reading molecule(s) from file: ", wholedict[i]['fname'])
        mols = reader.read_mols(wholedict[i]['fname'])
        qmethod, qbasis = reader.separated_theory(wholedict[i]['theory'])
        short_tag = wholedict[i]['tagkey']

        for j, mol_j in enumerate(mols):
            if molname is not None and mol_j.GetTitle() != molname:
                continue
            data_array = np.array(list(map(float,
                pt.get_sd_list(mol_j, short_tag, 'Psi4', qmethod, qbasis))))

        if take_relative:
            data_array = data_array - data_array[0]
            #data_array = data_array/data_array[0]
        if har_to_kcal:
            data_array = 627.5095*data_array

        titles.append(mol_j.GetTitle())
        labels.append(wholedict[i]['label'])
        yarray.append(data_array)
        xarray.append(range(len(data_array)))

    if verbose:
        header = '{}\n'.format(molname)
        for l in labels:
            header += ("%s\n" % l)
        xydata = np.vstack((xarray[0], yarray)).T
        np.savetxt('combined.dat', xydata, delimiter='\t', header=header,
            fmt=' '.join(['%i'] + ['%10.4f']*numFiles))

    # letter labels for x-axis
    num_confs = len(xarray[0])
    letters = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    rpt = int((num_confs / 26) + 1)
    xlabs = [''.join(i)
             for i in itertools.product(letters, repeat=rpt)][:num_confs]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    xlabel='conformer'
    ylabel="energy"

    # vtl print max range of relative energies
    conf_then_file = np.array(yarray).T
    ranges = []
    for c in conf_then_file:
        c_spread = max(c)-min(c)
        ranges.append(c_spread)
    print(f'mol {molname} max range: {max(ranges)}')

    ax.set_prop_cycle(plt.cycler('color', plt.cm.rainbow(np.linspace(0, 1, len(yarray)))))
    for i, (xs, ys) in enumerate(zip(xarray,yarray)):
        plt.plot(xs, ys, '-o', lw=0.8, label=labels[i])

    # publication view
#    plt.ylabel(ylabel,fontsize=8)
#    plt.xlabel(xlabel,fontsize=8)
#    plt.legend(bbox_to_anchor=(0.08,1.05),loc=3,fontsize=8)
#    fig.set_size_inches(3.37,1.7)

    # standard view
    plt.ylabel(ylabel,fontsize=14)
    plt.xlabel(xlabel,fontsize=14)
    plt.xticks(list(range(num_confs)), xlabs)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

    plt.title(molname)
    plt.grid()
    plt.savefig(figname, bbox_inches='tight',dpi=300)
    plt.show()
