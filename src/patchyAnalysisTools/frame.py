"""
This file holds the python bindings of the frame class.
"""
from ._patchyAnalysisTools import Frame as Frame_cpp

class Frame(Frame_cpp):
    """
    Python bindings for the Frame class.
    """
    pass