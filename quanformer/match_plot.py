#!/usr/bin/env python
"""
match_plot.py

Purpose:    Generate heat plots, scatter plots, and 3D plots from output of match_minima.py.
            This extends analysis capabilities of match_minima.py, which
            generates simple bar plots (for average compute times) and
            line plots (for relative conformer energies).

By:         Victoria T. Lim
Version:    Nov 30 2018

"""

import os
import sys
import numpy as np
import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl


def shift_array(rmsArray):
    """
    Place the first element of array in diagonal spot.
    Order of everything before and everything after is retained.
    Example,
      0 4 6                0 4 6
      0 2 7  --becomes-->  2 0 7
      0 1 3                1 3 0

    """
    for i, sublist in enumerate(rmsArray):
        sublist.insert(i, sublist.pop(0))
    return rmsArray


def plot_heat_rmse(molName,
                   rmsArray,
                   ticklabels,
                   ptitle='RMS error (kcal/mol)',
                   fprefix='rmse',
                   colors='PRGn_r'):
    """
    """
    plttitle = "%s\n%s" % (ptitle, molName)
    figname = "%s_%s.png" % (fprefix, molName)
    x = list(range(len(rmsArray)))
    y = list(range(len(rmsArray)))

    plt.figure(figsize=(10, 5))

    ### Tranpose and plot data - imshow swaps x and y
    plt.imshow(np.asarray(rmsArray).T, cmap=colors, origin='lower')
    cbar = plt.colorbar()

    ### Label figure. Label xticks before plot for better spacing.
    #    plt.title(plttitle,fontsize=20)
    plt.xticks(x, ticklabels, fontsize=12, rotation=-40, ha='left')
    plt.yticks(y, ticklabels, fontsize=14)
    plt.xlabel("reference", fontsize=14)
    plt.ylabel("compared", fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    ### Save/show plot.
    plt.savefig(figname, bbox_inches='tight')
    #    plt.show()
    plt.clf()


def plot_rmse_time(molName, eneArray, timeArray, ticklabels,
                   fprefix='scatter'):
    """
    """
    plttitle = "RMS error vs. log ratio of wall time\n%s" % molName
    figname = "%s_%s.png" % (fprefix, molName)
    colors = mpl.cm.rainbow(np.linspace(0, 1, len(eneArray)))
    markers = [
        "x", "^", "8", "d", "o", "s", "*", "p", "v", "<", "D", "+", ">", "."
    ] * 10

    # use plt.plot instead of scatter to label each point
    for i, (x, y) in enumerate(zip(eneArray, timeArray)):
        plt.scatter(x, y, c=colors[i], marker=markers[i], label=ticklabels[i])

    ### Label figure. Label xticks before plot for better spacing.
    plt.title(plttitle, fontsize=20)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
    #plt.xticks(x,ticklabels,fontsize=12,rotation=-20, ha='left')
    #plt.yticks(y,ticklabels,fontsize=12)
    plt.xlabel("RMS error (kcal/mol)", fontsize=14)
    plt.ylabel("log ratio of wall time", fontsize=14)

    ### Edit legend colors. All is one color since each sublist
    # colored by spectrum.
    ax = plt.gca()
    leg = ax.get_legend()
    for i in range(len(eneArray)):
        leg.legendHandles[i].set_color(colors[i])

    ### Save/show plot.
    plt.gcf().set_size_inches(8, 6)
    plt.savefig(figname, bbox_inches='tight')
    #    plt.show()
    plt.clf()


def match_plot(args):

    # Read input file and store each file's information in two lists.
    dat_list = []
    thry_list = []
    with open(args.infile) as f:
        for line in f:
            # ignore commented lines
            if line.startswith('#'):
                continue
            dataline = [x.strip() for x in line.split(',')]
            # ignore blank lines
            if dataline == ['']:
                continue
            thry_list.append(dataline[0])
            dat_list.append(dataline[1])

    # HEAT PLOT OF RMSE
    if args.eheatplot:
        rmsArray = []
        for datfile in dat_list:
            with open(datfile) as f:
                for line in f:
                    if "RMS error" in line:
                        rmse = next(itertools.islice(f, 1))
                        rmse = [float(s) for s in rmse.split()[1:]]
                        rmsArray.append(rmse)
                        break
        rmsArray = shift_array(rmsArray)

        # adjust data when dat files have more columns than # of method files were input, identify which index first
        #[l.pop(-1) for l in rmsArray]

        plot_heat_rmse(args.title, rmsArray, thry_list)

    # HEAT PLOT OF TIMES
    if args.theatplot:
        rmsArray = []
        for datfile in dat_list:
            with open(datfile) as f:
                for line in f:
                    if "avg time" in line:
                        ravgs = itertools.islice(f, 3)  # ratio line via iter
                        for j in ravgs:  # get last item of iterator
                            pass
                        ravgs = [float(s) for s in j.split()[1:]]
                        rmsArray.append(ravgs)
                        break
        rmsArray = shift_array(rmsArray)

        # adjust data when dat files have more columns than # of method files were input, identify which index first
        #[l.pop(-1) for l in rmsArray]

        # plot log ratio of wall times
        plot_heat_rmse(
            args.title,
            np.log10(rmsArray),
            thry_list,
            ptitle='log ratio of wall times',
            fprefix='times',
            colors='seismic')
        # plot direct ratio of wall times
        #plot_heat_rmse(args.title,rmsArray,thry_list, ptitle='ratio of wall times',fprefix='times')

    # SCATTER PLOT OF TIME VS RMSE
    if args.etscatter:
        eArray = []
        tArray = []
        for datfile in dat_list:
            with open(datfile) as f:
                for line in f:
                    if "RMS error" in line:
                        rmse = next(itertools.islice(f, 1))
                        rmse = [float(s) for s in rmse.split()[1:]]
                        eArray.append(rmse)
                    if "avg time" in line:
                        ravgs = itertools.islice(f, 3)  # ratio line via iter
                        for j in ravgs:  # get last item of iterator
                            pass
                        ravgs = [float(s) for s in j.split()[1:]]
                        tArray.append(ravgs)
                        break
        eArray = shift_array(eArray)
        tArray = shift_array(tArray)

        # adjust data when dat files have more columns than # of method files were input, identify which index first
        #[l.pop(-1) for l in eArray]
        #[l.pop(-1) for l in tArray]

        for i in range(len(dat_list)):
            #            if i<13: continue
            plot_rmse_time(
                args.title,
                eArray[i],
                np.log10(tArray[i]),
                thry_list,
                fprefix='scatter' + str(i + 1))


def single_scatter(args):
    # SINGLE SCATTER PLOT
    eArray = []
    tArray = []
    thry_list = []
    read_thry_list = True
    with open(args.infile) as f:
        for line in f:
            try:
                index1 = line.split()[1]
                if read_thry_list and index1.isdigit():
                    thry_list.append(line.split()[2])
            except IndexError:
                read_thry_list = False
                pass
            if "RMS error" in line:
                rmse = next(itertools.islice(f, 1))
                rmse = [float(s) for s in rmse.split()[1:]]
                eArray.append(rmse)
            if "avg time" in line:
                ravgs = itertools.islice(f, 3)  # ratio line via iterator
                for j in ravgs:  # get last item of iterator
                    pass
                ravgs = [float(s) for s in j.split()[1:]]
                tArray.append(ravgs)
                break
    eArray = shift_array(eArray)
    tArray = shift_array(tArray)

    # delete nan indices; bc must also know which methods to remove from thry_list
    # take [0] bc *Array is list of lists but only has one elem in single_scatter fx
    nanIndices = np.argwhere(np.isnan(eArray[0]))
    eArray = np.delete(eArray[0], nanIndices)
    tArray = np.delete(tArray[0], nanIndices)
    thry_list = np.delete(thry_list, nanIndices)

    plot_rmse_time(
        args.title, eArray, np.log10(tArray), thry_list, fprefix='scatter')


### ------------------- Parser -------------------

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("-i", "--infile", required=True,
        help=("Input text file with info on dat files (from match_minima.py) "
              "and associated levels of theory. This file is analogous to the "
              "input file from match_minima.py except that the files should "
              "point to the relative energies .dat files instead of to the "
              "molecule SDF files. Also don't need to have True/False for "
              "SPE/OPT. *** The number of entries in this input file MUST MATCH "
              "the number of columns in each dat file. Else the heat plots "
              "won't be square/sensible. However, if you only want to generate "
              "a single scatter rmse vs. time plot, then use the numscatter input flag."))

    parser.add_argument("-t", "--title", default="",
        help="Name/identifier of molecule")

    parser.add_argument("--eheatplot", action="store_true", default=False,
        help=("Generate heat plot of root mean square error of energies by methods. "
              "Looks for line in .dat file with the following: "
              "\"RMS errors by level of theory\""))

    parser.add_argument("--theatplot", action="store_true", default=False,
        help=("Generate heat plot of relative optimization times by methods. "
              "Looks for line in .dat file with the following: "
              "\"Four rows: (1) avg times,\""))

    parser.add_argument("--etscatter", action="store_true", default=False,
        help="Generate scatter plots of log ratio of opt times vs. RMSE.")

    parser.add_argument("--onescatter", action="store_true", default=False,
        help=("Generate single scatter plot from one of the .dat files of "
              "match_minima.py. This may be desired when a mol of some method "
              "failed to optimize. In that case, match_minima will have a column "
              "of nan's for the RMS errors of that method. Since that method "
              "does not have a .dat file for that mol, heat plots will not "
              "be aligned, and scatter plots will not be accurate unless that "
              "column is deleted in all other dat files (i.e., see lines w/pop). "
              "This flag allows user to specify one .dat file AS THE INPUT file of "
              "this match_plot.py script to generate a time vs rmse scatter plot."))

    args = parser.parse_args()
    if args.onescatter:
        single_scatter(args)
    else:
        match_plot(args)
