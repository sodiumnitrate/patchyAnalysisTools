"""
patchyAnalysisTools
by Irem Altan

See LICENSE AND README.md in https://github.com/sodiumnitrate/patchyAnalysisTools.

Python bindings for the Clusters class.
"""
from ._patchyAnalysisTools import Clusters as Clusters_cpp

class Clusters(Clusters_cpp):
    """
    Python bindings for the Clusters class.
    """
    pass