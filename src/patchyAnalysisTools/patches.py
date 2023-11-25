"""
This file holds the python bindings of the patches class.
"""
from ._patchyAnalysisTools import Patches as Patches_cpp

class Patches(Patches_cpp):
    """
    Python bindings for the Patches class.
    """
    pass