#!/usr/bin/env python
"""
pipeline.py

Main interface to Quanformer pipeline.

Version:    Mar 28 2019
By:         Victoria T. Lim

"""

import os

import quanformer.initialize_confs as initialize_confs
import quanformer.filter_confs     as filter_confs
import quanformer.confs_to_psi     as confs_to_psi
import quanformer.get_psi_results  as get_psi_results


def name_manager(infile):
    """
    File-checking and parsing for internal use.

    Parameters
    ----------
    infile : str
        Ex. "/path/to/file.sdf"

    Returns
    -------
    curr_dir : str
        Ex. "/path/to"
    checked_infile : str
        Ex. "/path/to/file.sdf"
    prefix : str
        Ex. "file"
    ext : str
        Ex. ".sdf"
    no_path_infile : str
        Ex. "file.sdf"

    """
    curr_dir = os.getcwd()

    if not os.path.exists(infile):
        raise FileNotFoundError("No such file: {}".format(infile))

    # if infile does not contain full path, then its path is curr dir
    checked_infile = os.path.abspath(infile)

    # get base name without suffix and without extension
    inpath, no_path_infile = os.path.split(checked_infile)

    # get extension of .sdf, .smi, etc.
    all_but_ext, ext = os.path.splitext(checked_infile)

    # replace - with # and split by # to get basename without suffix/extension
    prefix = os.path.basename(all_but_ext).replace('-', '#').split('#')[0]

    return curr_dir, checked_infile, prefix, ext, no_path_infile


def setup_conformers(infile):
    """
    Generate and filter conformers for input list of SMILES strings.
    Output is saved in SDF file.

    Parameters
    ----------
    infile : str
        filename of the SMILES input data

    """
    curr_dir, checked_infile, prefix, ext, no_path_infile = name_manager(
        infile)

    # append value of suffix to filename after MM opt/filtering
    suffix = '200'

    if ext != '.smi':
        raise ValueError("Input should be a SMILES file with extension of .smi")
    print("\nGenerating and filtering conformers for {}".format(infile))

    # generate conformers and MM opt
    initialize_confs.initialize_confs(checked_infile)

    # set filenames; prefix is same filename (w/o extension) of input file
    pre_filt = os.path.join(curr_dir, prefix + '.sdf')
    post_filt = os.path.join(curr_dir, "{}-{}.sdf".format(
        prefix, suffix))

    # filter conformers
    filter_confs.filter_confs(pre_filt, "MM Szybki SD Energy",
                              post_filt)


def setup_calculations(infile, method, basisset, calctype='opt', mem='5.0 Gb'):
    """
    Write input files for Psi4 calculations.

    Parameters
    ----------
    infile : str
        filename of the SDF input molecules
    method : str
        name of QM method
    basisset : str
        name of QM basis set
    calctype : str
        'opt' for geometry optimizations,
        'spe' for single point energy calculations,
        'hess' for Hessian calculation
        default in pipeline is 'opt'
    mem : str
        allotted memory for each Psi4 calculation, default in pipeline is 5 Gb


    """
    curr_dir, checked_infile, prefix, ext, no_path_infile = name_manager(
        infile)

    # check that specified calctype is valid
    if calctype not in {'opt', 'spe', 'hess'}:
        raise ValueError("Specify a valid calculation type.")

    print("\nCreating Psi4 input files for %s..." % prefix)
    confs_to_psi.confs_to_psi(checked_infile, method, basisset, calctype, mem)


def process_results(infile, calctype='opt', suffix=[], psiout='output.dat', timeout='timer.dat'):
    """
    Process Psi4 output files and filter conformers.

    Parameters
    ----------
    infile : str
        filename of the SDF input molecules
    suffix : list
        suffix of filename if not following numbering convention of pipeline
        ('200,'210','220','221','222'). this list should contain one string
        for appending to filename with extracted QM results. if calctype is
        'opt', the list should also contain a second string for extracted
        and filtered QM results. Ex. ['qm','qmfilt']

    """
    curr_dir, checked_infile, prefix, ext, no_path_infile = name_manager(
        infile)

    # default of pipeline goes '200' --> '210'/'220' --> '221/'222'
    if len(suffix) > 0:
        out_results = os.path.join(
            curr_dir, "{}-{}.sdf".format(prefix, suffix[0]))
        if len(suffix) == 2:
            out_filter = os.path.join(
                curr_dir, "{}-{}.sdf".format(prefix, suffix[1]))
    else:
        if '-200.sdf' in no_path_infile:
            out_results = os.path.join(curr_dir, prefix + '-210.sdf')
            out_filter = os.path.join(curr_dir, prefix + '-220.sdf')
        elif '-220.sdf' in no_path_infile:
            out_results = os.path.join(curr_dir, prefix + '-221.sdf')
            out_filter = os.path.join(curr_dir, prefix + '-222.sdf')
        else:
            raise ValueError(
                "ERROR: Input file does not have usual 200-series "
                "suffixes (see README).\nPlease specify suffix(es) in a "
                "list in accordance with documentation.")

    # get psi4 results
    print("Getting Psi4 results for %s ..." % (checked_infile))
    method, basisset = get_psi_results.get_psi_results(
        checked_infile, out_results, calctype=calctype, psiout=psiout, timeout=timeout)

    # only filter structures after opts; spe/hess should not change geoms
    if calctype == 'opt' and None not in [method, basisset]:
        filter_results(out_results, out_filter, method, basisset)

def filter_results(infile, outfile, method, basisset):
    # may call this function directly to filter and not extract data

    tag = "QM Psi4 Final Opt. Energy (Har) %s/%s" % (method, basisset)
    print("Filtering Psi4 results for %s ..." % (outfile))
    filter_confs.filter_confs(infile, tag, outfile)

