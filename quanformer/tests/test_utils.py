"""
test_utils.py
"""
import sys
import os
import pytest
import shutil

# define location of input files for testing
mydir = os.path.dirname(os.path.abspath(__file__))

# import functions to aid testing
sys.path.append(os.path.join(os.path.dirname(__file__), 'helpers'))
from helper import *

from quanformer.utils import *

# -----------------------

def test_convert_extension():
    convert_extension(
        os.path.join(mydir, 'data_tests', 'carbon-222.sdf'),
        os.path.join(mydir, 'data_tests', 'carbon-222.mol2'))
    assert os.path.getsize(os.path.join(mydir, 'data_tests', 'carbon-222.mol2')) == 1320
    os.remove(os.path.join(mydir, 'data_tests', 'carbon-222.mol2'))

# test manually without pytest
if 0:
    sys.path.insert(0, '/home/limvt/Documents/quanformer/quanformer')
    from utils import *
    test_convert_extension()
