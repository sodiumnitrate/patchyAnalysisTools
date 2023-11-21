"""
This file holds the python bindings of the frame class.
"""
from .patchyAnalysisTools_cpp import Frame as Frame_cpp

class Frame(Frame_cpp):
    """
    Python bindings for the Frame class.
    """
    pass