Installing quanformer
====================

Anaconda
--------

1. Install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_ (or Anaconda) with Python 3.6 (version 3.7 is not yet tested).

2. Create a conda environment for **quanformer**::

    conda create --name quanformer python=3.6 matplotlib scipy
    conda activate quanformer

3. Obtain OpenEye, Psi4, and Psi4 dependencies::

    conda install -c openeye openeye-toolkits
    conda install -c psi4 psi4 dftd3 gcp [gpu_dfcc]

4. Check installations.

    a. Psi4::

        psi4 --test

    b. OpenEye (in interactive Python session)::

        import openeye.oechem as oechem
        mol = oechem.OEMol()


Quanformer
----------
 
Install **quanformer** from source by the following:


1. Clone the repository from `Github <https://github.com/MobleyLab/quanformer>`_::

    git clone https://github.com/MobleyLab/quanformer
    cd quanformer

2. Install::

    python setup.py install

   or use ``pip`` for a local install (editable installation)::

    pip install -e .

