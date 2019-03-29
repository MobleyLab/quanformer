"""
simple_use_case.py
--------------------------
Final directory structure:

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

"""
import quanformer.pipeline as qp

qp.setup_conformers('two_alkanes.smi')
qp.setup_calculations('two_alkanes-200.sdf','mp2','def2-sv(p)')
#qf.pipeline.process_results('two_alkanes-200.sdf')
