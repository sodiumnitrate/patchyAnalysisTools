from patchyAnalysisTools import Frame, Patches, Clusters

class TestFrame:
    def test_init(self):
        coords = [[1, 1, 1], [2, 2, 2]]
        orients = [[0,0,0],[0,0,0]]
        frame = Frame(coords, orients, [3, 3, 3])

        assert frame.get_N() == 2
        assert frame.get_time_stamp() == 0
        assert len(frame.get_coordinates()) == 2
        assert len(frame.get_orientations_as_angles()) == 2

    def test_set_types(self):
        coords = [[1, 1, 1], [2, 2, 2]]
        orients = [[0,0,0],[0,0,0]]
        frame = Frame(coords, orients, [3, 3, 3])

        types = frame.get_types()
        assert types[0] == 0
        assert types[1] == 0

        types_2 = [1, 5]
        frame.set_types(types_2)
        assert frame.get_types()[0] == 1
        assert frame.get_types()[1] == 5

        frame.set_time_stamp(5)
        assert frame.get_time_stamp() == 5

    def test_interactions(self):
        coords = []
        orients = []
        sep = 1.05
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    coords.append([i * sep, j * sep, k * sep])
                    orients.append([0,0,0])
        cell = [sep * 3, sep * 3, sep * 3]

        patches = Patches()
        patches.add_patch(1, 1.1, 0.92, [1,0,0], 0)
        patches.add_patch(1, 1.1, 0.92, [-1,0,0], 0)
        patches.add_patch(1, 1.1, 0.92, [0,1,0], 0)
        patches.add_patch(1, 1.1, 0.92, [0,-1,0], 0)
        patches.add_patch(1, 1.1, 0.92, [0,0,1], 0)
        patches.add_patch(1, 1.1, 0.92, [0,0,-1], 0)
        patches.make_all_patches_adjacent()

        frame = Frame(coords, orients, cell)
        frame.determine_bond_list(patches)
        bond_list = frame.get_bond_list()
        assert len(bond_list) == len(coords) * 3

    def test_clusters(self):
        coords = []
        orients = []
        sep = 1.05
        for i in range(3):
            coords.append([i * sep, 0, 0])
            orients.append([0,0,0])
        cell = [sep * 3, sep * 3, sep * 3]

        patches = Patches()
        patches.add_patch(1, 1.1, 0.92, [1,0,0], 0)
        patches.add_patch(1, 1.1, 0.92, [-1,0,0], 0)
        patches.make_all_patches_adjacent()

        frame = Frame(coords, orients, cell)
        frame.determine_bond_list(patches)
        frame.determine_clusters()
        frame.determine_percolation(patches)
        clusters = frame.get_clusters()
        clusters_list = clusters.get_clusters()
        assert len(clusters_list) == 1
        assert len(clusters_list[0]) == 3