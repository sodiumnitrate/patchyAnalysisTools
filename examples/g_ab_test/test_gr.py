"""
This file exists to test all the g(r) results against those 
determined with VMD.
"""

import patchyAnalysisTools.trajectory as trj
import matplotlib.pyplot as plt
import pdb
import numpy as np

# load trajectory
traj = trj.trajectory(file_name="final.conf")

# get the last frame
last_frame = traj.get_last_frame()

# get g(r)
r,gr = last_frame.calculate_rdf()

# get g_00(r)
r,gr_00 = last_frame.calculate_rdf(selection=(0,0))

# get g_11(r)
r,gr_11 = last_frame.calculate_rdf(selection=(1,1))

# get g_01(r)
r,gr_01 = last_frame.calculate_rdf(selection=(0,1))

# get g_10(r)
r,gr_10 = last_frame.calculate_rdf(selection=(1,0))

# load VMD data
data0 = np.loadtxt('g00.dat')
data1 = np.loadtxt('g11.dat')
data01 = np.loadtxt('g01.dat')

# plot
plt.plot(r,gr,label="g(r)")
plt.plot(r,gr_00,label="g_00")
plt.plot(data0[:,0],data0[:,1],label="g_00 from vmd")
plt.plot(data1[:,0],data1[:,1],label="g_11 from vmd")
plt.plot(data01[:,0],data01[:,1],label='g_01 from vmd')
plt.plot(r,gr_11,label="g_11")
plt.plot(r,gr_01,label="g_01")
plt.plot(r,gr_10,label="g_10")
plt.legend()
plt.yscale('log')
plt.show()