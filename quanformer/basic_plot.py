#!/usr/bin/env python
"""
basic_plot.py

Basic plots (lines or bars) of information in SD tags.

Version:    Apr 12 2019
By:         Victoria T. Lim


TODO:
- plot same tag from multiple files (plotMultSDF.py) (http://tinyurl.com/yy62nqqt)
- plot diff tags in same file

"""
import numpy as np

import matplotlib.pyplot as plt

import quanformer.proc_tags as pt
import quanformer.reader as reader


def basic_plot(infile, tag, style, take_relative=False, har_to_kcal=False):
    """
    TODO

    Parameters
    ----------
    infile : string
        Name of SDF file with information in SD tags.
    tag : string

    style : string
        plot style. can be 'scatter', 'line', or 'bar'
        TODO

    """
    # Open molecule file.
    mols = reader.read_mols(infile)

    for i, mol_i in enumerate(mols):

        # get array of all conformer data of this mol
        try:
            data_array = np.fromiter(pt.get_sd_list(mol_i, datum='', taglabel=tag), dtype=np.float64)
        except ValueError:
            data_array = np.asarray([np.nan])*mol_i.NumConfs()

        # exclude conformers for which job did not finish (nan)
        nanIndices = np.argwhere(np.isnan(data_array))
        for i in reversed(nanIndices):  # loop in reverse to delete correctly
            data_array = np.delete(data_array, i)

        # take relative values by subtracting minimum
        if take_relative:
            data_array = data_array - np.amin(data_array)

        # generate plot
        plt.plot(data_array)
        plt.title(mol_i.GetTitle()+'\n'+tag, fontsize=14)
        plt.savefig(f'output_{i}.png', bbox_inches='tight')
        plt.show()

