"""
simple_use_case.py
"""
import sys
import os
import pytest

os.system("python ../executor.py -f data_tests/two_alkanes.smi --setup -m 'mp2' -b 'def2-sv(p)'")

# TODO: add checks for these?
# files that should be created
"""
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
