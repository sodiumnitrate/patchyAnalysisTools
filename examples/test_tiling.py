"""
This file tests the stitch_snapshots function in utils.
"""

from patchyAnalysisTools.utils import stitch_snapshots
import pdb
import matplotlib.pyplot as plt
import numpy as np

# as a test, just tile the same box in a 3 by 3 array
list_of_frames = []
for i in range(27):
    list_of_frames.append("trajectory.conf")

# stitch snapshots together
frame = stitch_snapshots(list_of_frames)

# write to .xyz file to visualize with vmd/ovito
frame.write_xyz("tiled_3_by_3.xyz")

# write pdb file to be able to calculate g(r) with vmd
frame.write_frame_pdb("tiled_3_by_3.pdb")

# calculate gr
r,gr = frame.calculate_rdf()

# load VMD data
data = np.loadtxt("gr_vmd_tiled_3_by_3.dat")

# plot
plt.plot(r,gr,label='g(r)')
plt.plot(data[:,0],data[:,1],label='g(r) from vmd')
plt.legend()
plt.yscale('log')
plt.show()