from patchyAnalysisTools import Patches

class TestPatches:
    def test_init_from_file(self):
        patches = Patches()
        patches.read_from_file("aux_files/patches.dat")

        assert patches.get_n_patch() == 5
        assert patches.get_max_lambda() == 1.1

    def test_init_add_patch(self):
        patches = Patches()
        
        assert patches.get_n_patch() == 0
        assert patches.get_max_lambda() == 0

        patches.add_patch(1, 1.5, 0.92, [0, 1, 0], 0)

        assert patches.get_n_patch() == 1
        assert patches.get_max_lambda() == 1.5