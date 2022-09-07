import matplotlib.pyplot as plt
import pdb
import patchyAnalysisTools.trajectory as trj
import patchyAnalysisTools.utils as utils

# load trajectory
traj = trj.trajectory(file_name="trajectory.conf")
# read patch file
traj.set_patch_info("patches.dat")
print("number of patches: ", traj.patches.n_patch)


# get the first ten frames as a subtrajectory
subtraj = traj.get_frame(0,10)

# write the subtrajectory to a file
subtraj.write_xyz("frames_0_10.xyz")

# get the first frame
one_frame = traj.get_frame(0)

# get bond probability
pb = one_frame.get_bond_probability()
print("bond probability: %lf"%pb)

# create a sliced image from the frame
#slice = utils.get_frame_slice(one_frame,grid_spacing=0.02,rad=25,start=0,end=5)
#plt.imshow(slice)
#plt.show()

# get the last frame
last_frame = traj.get_last_frame()
last_frame.get_list_of_interacting_pairs()
last_frame.check_calculated_bonds_against_bond_numbers()

last_frame.set_cluster_info()
print("There are %d clusters!"%len(last_frame.cluster_info.clusters))
#trajectory.plot_clusters(last_frame,clusters)

# find percolating clusters and check if the system is percolated
last_frame.find_percolating_clusters()
last_frame.is_system_percolated()
if last_frame.percolated:
    print("system is percolated!")
else:
    print("system is not percolated!")

# write frame to xyz file with cluster info
last_frame.write_frame_with_cluster_info("last_frame_clusters.xyz")
last_frame.write_frame_pdb("last_frame.pdb",write_bond_info=True)

# calculate rgs
last_frame.get_cluster_rg()
plt.hist(last_frame.cluster_info.R_g)
plt.title("histogram of cluster R_g")
plt.show()

# get number of loops
plt.hist(last_frame.cluster_info.N_loop)
plt.title("histogram of cluster N_loops")
plt.show()

# print length of chains
plt.hist(last_frame.cluster_info.L_chain)
plt.title("histogram of cluster L_chain")
plt.show()

# write the biggest cluster to file
last_frame.write_biggest_cluster_xyz("biggest_cluster.xyz")

print("energy: ", last_frame.energy[0])
print("number of bonds: ", len(last_frame.bonds_calculated[0]))