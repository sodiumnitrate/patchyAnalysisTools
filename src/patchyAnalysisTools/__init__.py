"""
__init__ file for the patchyAnalysiosTools module
"""
from __future__ import annotations
from .patchyAnalysisTools_cpp import __doc__, __version__

from .frame import Frame
from .trajectory import Trajectory

__all__ =["__doc__", "__version__", "Frame", "Trajectory"]