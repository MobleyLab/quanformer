.. quanformer documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

quanformer: quantum mechanical analysis of molecular conformers
=========================================================

**quanformer** is a Python-based pipeline for generating conformers, preparing quantum mechanical (QM) calculations, and processing QM results for a set of input molecules. 

For each molecule, conformers are generated and optimized using a molecular mechanics force field.
Input files for QM calculations are then prepared for geometry optimizations, single point energy (SPE) calculations, or Hessian calculations.
The user specifies any QM method and basis set that is supported in the QM software package (Psi4 or Turbomole).
With completed calculations, the pipeline extracts final energies and geometries and collects job-related details such as calculation time and number of optimization steps.
Analysis scripts are provided for comparing different methods in the lens of conformer energies and calculation times via easily-interpreted plots.
This pipeline was tested to robustly process hundreds of conformers per molecule and hundreds of molecules, though it can likely handle more.

This code is open source on `Github <https://github.com/MobleyLab/quanformer>`_.

.. toctree::
    :maxdepth: 1
    :caption: User Documentation

    install
    api
    examples
