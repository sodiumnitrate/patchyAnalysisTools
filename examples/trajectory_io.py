from more_itertools import last
from patchyAnalysisTools.trajectory import trajectory
import matplotlib.pyplot as plt
import pdb

# load trajectory
traj = trajectory.read_trajectory("trajectory.conf")

# get the first ten frames as a subtrajectory
subtraj = traj.get_frame(0,10)

# write the subtrajectory to a file
subtraj.write_xyz("frames_0_10.xyz")

# get the first frame
one_frame = traj.get_frame(0)

# read patch file
patches_object = trajectory.read_patch_info("patches.dat")
print("number of patches: ", patches_object.n_patch)

pb = trajectory.get_bond_probability(one_frame,patches_object)
print("bond probability: %lf"%pb)

# create a sliced image from the frame
#slice = trajectory.get_frame_slice(one_frame,grid_spacing=0.02,rad=25,start=0,end=5)
#plt.imshow(slice)
#plt.show()

# get the last frame
last_frame = traj.get_last_frame()
bonds = last_frame.get_list_of_interacting_pairs(patches_object)
check = trajectory.check_calculated_bonds_against_bond_numbers(last_frame,bonds)
if check:
    print("bond counts match!")

clusters = trajectory.find_all_clusters(bonds)
print("There are %d clusters!"%len(clusters))
trajectory.plot_clusters(last_frame,clusters)

percolated_clusters = trajectory.find_percolating_clusters(last_frame, clusters)
if trajectory.is_system_percolated(last_frame,clusters):
    print("system is percolated!")
else:
    print("system is not percolated!")