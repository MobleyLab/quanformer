Installing quanformer
====================

Install **quanformer** from source by the following:


1. Clone the repository from `Github <https://github.com/MobleyLab/quanformer>`_::

    git clone https://github.com/MobleyLab/quanformer
    cd quanformer

2. If you wish to install this package into a new conda environment::

    conda env create -f devtools/conda-envs/quanformer.yaml python=3.6 
    conda activate quanformer

3. Install::

    python setup.py install

or use ``pip`` for a local install (editable installation)::

    pip install -e .

