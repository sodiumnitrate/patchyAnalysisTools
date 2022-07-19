from patchyAnalysisTools.trajectory import trajectory
import matplotlib.pyplot as plt

# load trajectory
traj = trajectory.read_trajectory("trajectory.conf")

# get the first ten frames as a subtrajectory
subtraj = traj.get_frame(0,10)

# write the subtrajectory to a file
subtraj.write_xyz("frames_0_10.xyz")

# get the first frame
one_frame = traj.get_frame(0)

# create a sliced image from the frame
slice = trajectory.get_frame_slice(one_frame,grid_spacing=0.02,rad=25,start=0,end=5)
plt.imshow(slice)
plt.show()
