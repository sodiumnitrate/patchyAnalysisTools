# patchyAnalysisTools

A set of tools to analyze the output of patchy particle simulations. Works with a very specific format of trajectory files, but the I/O methods can be extended for the analysis tools to work on arbitrary trajectory file formats. See examples for details.

## Installation

	git clone https://github.com/sodiumnitrate/patchyAnalysisTools.git
	cd patchyAnalysisTools
	python setup.py install

## Examples

(See the `examples` folder for more details.)

### Loading trajectories

To load a trajectory:

	traj = trj.trajectory(file_name="traj.conf")

where `traj.conf` is a trajectory file with specific formatting requirements (see examples). Patch information is then read from a file and set for the `traj` object by

	traj.set_patch_info("patches.dat")

where `patches.dat` is again a special file (see examples). _Some_ analysis can be done without reading and setting patch information, but most require this for determining bonds, etc.

### Selecting frames or subtrajectories

To select the last frame:

	last_frame = traj.get_last_frame()

To select an arbitrary frame:

	frame = traj.get_frame(10)

To select a range of frames:

	subtraj = traj.get_frame(10,20)

### Bonding and cluster info

Determination of bonding and cluster information requires the patches to have been set in the trajectory before a frame is selected (see above).

First, we determine bonds:

	frame.get_list_of_interacting_pairs()

If the trajectory file has the number of bonds per particle written, we can double check the calculated bonds against these values:

	frame.check_calculated_bonds_against_bond_numbers()

We then determine the clusters by

	frame.set_cluster_info()

To find the percolating clusters:

	frame.find_percolating_clusters()

To check if system is percolated:

	frame.is_system_percolated()

which returns a boolean.

## What's missing
The main thing that is incomplete in this package is support for the output of Gibbs-Ensemble simulations (which contain two boxes). Some of the tools in this package support it, some don't. I may complete it at one point, but for now, everything is only tested for regular simulation outputs.

It is also important to note that not all the tools in this package are thoroughly tested for simulation output with non-cubic boxes.