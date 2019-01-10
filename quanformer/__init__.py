"""
Quanformer
Quanformer is a Python-based pipeline for generating conformers, preparing quantum mechanical (QM) calculations, and processing QM results for a set of molecules and their conformers.
"""

# Make Python 2 and 3 imports work the same
# Safe to remove with Python 3-only code
from __future__ import absolute_import

# Add imports here
from .initialize_confs import *  # lgtm [py/polluting-import]

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
