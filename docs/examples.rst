Examples
========

Example files are provided on `Github <https://github.com/MobleyLab/quanformer/tree/master/examples/main>`_.

Remember to activate your conda environment containing **quanformer**::

    conda activate quanformer

Basic usage
-----------
    
Set up structures and prepare calculations:

.. code-block:: python

    import quanformer.pipeline as qp

    qp.setup_conformers('two_alkanes.smi')
    qp.setup_calculations('two_alkanes-200.sdf','mp2','def2-sv(p)')

Your directory structure should look something like this::
    .
    ├── AlkEthOH_c1008
    │   ├── 1
    │   │   └── input.dat
    │   ├── 2
    │   │   └── input.dat
    │   └── 3
    │       └── input.dat
    ├── AlkEthOH_c312
    │   ├── 1
    │   │   └── input.dat
    │   ├── 2
    │   │   └── input.dat
    │   └── 3
    │       └── input.dat
    ├── numConfs.txt
    ├── two_alkanes-200.sdf
    └── two_alkanes.sdf
    
After running calculations from provided input files, collect results:

.. code-block:: python

    qp.process_results('two_alkanes-200.sdf')

