"""
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.

Python bindings for the Patches class.
"""
from ._patchyAnalysisTools import Patches as Patches_cpp

class Patches(Patches_cpp):
    """
    Python bindings for the Patches class.
    """
    pass