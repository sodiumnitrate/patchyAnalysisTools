"""
Unit tests for Frame.
"""

import tempfile

from patchyAnalysisTools import Frame

class TestFrame:
    def test_init(self):
        f = Frame([[0,0,0],[1,1,1]], [[0,0,0],[0,0,0]], [1, 1, 1])
        assert f.get_cell() == [1, 1, 1]

    def test_N_frame_num(self):
        f = Frame([[1,1,1],[2,2,2]], [[0,0,0], [0,0,0]], [2, 2, 2])
        f.set_frame_num(10)
        assert f.get_frame_num() == 10
        assert f.get_N() == 2
        assert f.frame_num == 10

    def test_get_coords_orients(self):
        f = Frame([[1,1,1],[2,2,2]], [[0,0,0], [0,0,0]], [2, 2, 2])
        coords = f.get_coordinates()
        orients = f.get_orientations()
        assert coords == [[1,1,1],[2,2,2]]
        assert orients == [[0,0,0], [0,0,0]]

    def test_write_xyz(self):
        f = Frame([[1,1,1],[2,2,2]], [[0,0,0], [0,0,0]], [2, 2, 2])
        with tempfile.TemporaryDirectory() as tmpdirname:
            f.write_xyz(f"{tmpdirname}/test.xyz")

    def test_particle_disp(self):
        f = Frame([[0,0,0], [1,1,1],[2,2,2]], [[0,0,0], [0,0,0]], [2, 2, 2])
        disp = f.get_displacement(0,2)
        assert disp == [0,0,0]