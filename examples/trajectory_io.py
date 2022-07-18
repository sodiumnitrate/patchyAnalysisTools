from patchyAnalysisTools.trajectory import trajectory

traj = trajectory.read_trajectory("trajectory.conf")

subtraj = traj.get_frame(0,10)

subtraj.write_xyz("frames_0_10.xyz")