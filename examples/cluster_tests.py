"""
This file contains some tests, focusing mostly on cluster properties.
"""

import patchyAnalysisTools.trajectory as trj
import matplotlib.pyplot as plt
import pdb

# load trajectory
traj = trj.trajectory(file_name="last_frame_clust.conf")
# read patch file
traj.set_patch_info("patches_2.dat")

# get the last frame
last_frame = traj.get_last_frame()
last_frame.get_list_of_interacting_pairs()
last_frame.check_calculated_bonds_against_bond_numbers()

# get cluster info
last_frame.set_cluster_info()
print("There are %d clusters!" % len(last_frame.cluster_info.clusters))

# find percolating clusters and check if system is percolated
last_frame.find_percolating_clusters()
last_frame.is_system_percolated()

if last_frame.percolated:
    print("system is percolated!")
else:
    print("system is not percolated!")

last_frame.cluster_info.get_cluster_props()

# write frame to xyz file with cluster info
last_frame.write_frame_with_cluster_info("last_frame_clusters_clust.xyz")

# write the biggest cluster
last_frame.write_biggest_cluster_xyz("biggest_cluster_clust.xyz")

# calculate cluster radii of gyration
last_frame.get_cluster_rg()

# get loops and write loops of the first cluster to file
loops = last_frame.cluster_info.loops
selected = []
sel_types = []
atom = 0
for l in loops[0]:
    # each loop with a different atom
    atom = (atom + 1) % 5
    for p in l:
        sel_types.append(atom)
        selected.append(p)
traj.write_xyz("first_cluster_loops.xyz",selected=selected, sel_types=sel_types)

# get chains and write longest chain for each cluster to file
chains = last_frame.cluster_info.chains
selected = []
sel_types = []
atom = 0
for c in chains:
    # each chain with a different atom
    atom = (atom + 1) % 5
    for p in c:
        sel_types.append(atom)
        selected.append(p)
traj.write_xyz("longest_chains.xyz",selected=selected, sel_types=sel_types)
