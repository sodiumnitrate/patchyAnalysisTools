"""
This file holds the python bindings of the Trajectory class.
"""
from .patchyAnalysisTools_cpp import Trajectory as Trajectory_cpp

class Trajectory(Trajectory_cpp):
    """
    Python bindings for the Trajectory class.
    """
    pass