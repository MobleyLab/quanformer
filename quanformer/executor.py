#!/usr/bin/env python
"""
executor.py

Main interface to Quanformer pipeline.

Version:    Jan 6 2019
By:         Victoria T. Lim

"""

import os, sys
import argparse

import initialize_confs
import filter_confs
import confs_to_psi
import get_psi_results


def name_manager(infile):
    curr_dir = os.getcwd()

    # if infile does not contain full path, then its path is curr dir
    checked_infile = os.path.abspath(infile)

    # get base name without suffix and without extension
    inpath, no_path_infile = os.path.split(checked_infile)

    # get extension of .sdf, .smi, etc.
    all_but_ext, ext = os.path.splitext(checked_infile)

    # replace - with # and split by # to get basename without suffix/extension
    prefix = os.path.basename(all_but_ext).replace('-', '#').split('#')[0]

    return curr_dir, checked_infile, prefix, ext, no_path_infile


def main(**kwargs):

    curr_dir, checked_infile, prefix, ext, no_path_infile = name_manager(
        opt['filename'])

    if opt['setup']:

        # default of pipeline uses '200' suffix for MM opt/filtering output
        if opt['suffix'] is None:
            suffix = '200'
        else:
            suffix = opt['suffix'][0]

        # if input is SMILES, do MM opt then filter
        if ext == '.smi':
            print("\nGenerating and filtering conformers for %s" %
                  opt['filename'])
            initialize_confs.initialize_confs(checked_infile)
            pre_filt = os.path.join(curr_dir, prefix + '.sdf')
            # prefix is base filename with same dir as input file
            post_filt = os.path.join(curr_dir, "{}-{}.sdf".format(
                prefix, suffix))
            filter_confs.filter_confs(pre_filt, "MM Szybki SD Energy",
                                      post_filt)
        # if input is SDF, don't generate confs/filter, generate QM inputs
        else:
            post_filt = checked_infile

        # generate Psi4 inputs
        print("\nCreating Psi4 input files for %s..." % prefix)
        confs_to_psi.confs_to_psi(post_filt, opt['method'], opt['basisset'],
                                  opt['calctype'], opt['mem'])

    else:  # ========== AFTER QM =========== #

        # default of pipeline goes '200' --> '210'/'220' --> '221/'222'
        if opt['suffix'] is None:
            if '-200.sdf' in no_path_infile:
                out_results = os.path.join(curr_dir, prefix + '-210.sdf')
                out_filter = os.path.join(curr_dir, prefix + '-220.sdf')
            elif '-220.sdf' in no_path_infile:
                out_results = os.path.join(curr_dir, prefix + '-221.sdf')
                out_filter = os.path.join(curr_dir, prefix + '-222.sdf')
            else:
                sys.exit(
                    "ERROR: Input file does not have usual 200-series "
                    "suffixes (see README for details).\nPlease specify "
                    "one suffix if calling setup or if extracting Hessian "
                    "data, or specify two suffixes if extracting/filtering "
                    "optimization results.\nSee usage in argparse.")
        else:
            out_results = os.path.join(
                curr_dir, "{}-{}.sdf".format(prefix, opt['suffix'][0]))
            if opt['calctype'] == 'opt':
                out_filter = os.path.join(
                    curr_dir, "{}-{}.sdf".format(prefix, opt['suffix'][1]))

        # get psi4 results
        print("Getting Psi4 results for %s ..." % (checked_infile))
        method, basisset = get_psi_results.get_psi_results(
            checked_infile, out_results, calctype=opt['calctype'])

        # only filter structures after opts; spe/hess should not change geoms
        if opt['calctype'] == 'opt':

            # if get_psi_results exits early (output file already exists),
            # look for method from command line call for filtering
            if None in [method, basisset]:
                if None in [opt['method'], opt['basisset']]:
                    print(
                        "\nERROR: no results obtained and no conformers "
                        "filtered. If you want to filter an already-existing "
                        "output file, specify method and basis set in command "
                        "line call with -m [method] -b [basis]\n")
                    return
                else:
                    method = opt['method']
                    basisset = opt['basisset']

            tag = "QM Psi4 Final Opt. Energy (Har) %s/%s" % (method, basisset)
            print("Filtering Psi4 results for %s ..." % (out_results))
            filter_confs.filter_confs(out_results, tag, out_filter)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", "--filename", required=True,
        help="SDF file (with FULL path) to be set up or processed.")

    # setup calculations or process results for specified calctype
    parser.add_argument("--setup", action="store_true", default=False,
        help="If True (default=False), generate and filter conformers (if "
              "starting from *.smi), then generate Psi4 input files.")
    parser.add_argument("--results", action="store_true", default=False,
        help="If True (default=False), process Psi4 output files and filter "
             "conformers.")
    parser.add_argument("-t", "--calctype", default="opt",
        help="Specify either 'opt' for geometry optimizations, 'spe' for "
             "single point energy calculations, or 'hess' for Hessian "
             "calculations. Default is 'opt'.")

    # qm job detail
    parser.add_argument("-m", "--method",
        help="Name of QM method. Put this in 'quotes'.")
    parser.add_argument("-b", "--basisset",
        help="Name of QM basis set. Put this in 'quotes'.")
    parser.add_argument("--mem", default="5.0 Gb",
        help="Memory specification for each Psi4 calculation.")

    # custom suffixes for pipeline outputs
    parser.add_argument("--suffix", nargs='+',
        help="For custom naming of results and filtered files throughout "
             "pipeline. If called with --setup option, include ONE suffix "
             "for output of filtered conformers. If suffix is called with "
             "--results option, include TWO suffixes for (1) output of QM "
             "calculations and (2) the filtered file of (1)."
             "Examples: --suffix 'filt'; --suffix 'qm' 'qmfilt' ")

    args = parser.parse_args()
    opt = vars(args)

    # check that both 'setup' and 'spe' are not both true or both false
    if opt['setup'] == opt['results']:
        raise parser.error(
            "Specify exactly one of either --setup or --results.")

    # check that input file exists
    if not os.path.exists(opt['filename']):
        raise parser.error(
            "Input file %s does not exist. Try again." % opt['filename'])

    # check that specified calctype is valid
    if opt['calctype'] not in {'opt', 'spe', 'hess'}:
        raise parser.error("Specify a valid calculation type.")

    main(**opt)
