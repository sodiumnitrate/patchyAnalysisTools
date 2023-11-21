"""
Unit tests for Frame.
"""

from patchyAnalysisTools import Frame

class TestFrame:
    def test_init(self):
        f = Frame([[0,0,0],[1,1,1]], [[0,0,0],[0,0,0]])