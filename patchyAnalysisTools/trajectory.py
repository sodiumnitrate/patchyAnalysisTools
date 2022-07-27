import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt
import random
import networkx as nx
from . import utils
from . import rdf_sq
from . import clusters
from . import patches

class frame():
    '''

    The frame class holds all info related to a snapshot in a simulation trajectory.
    
    '''
    def __init__(self, n_particles, n_frame, coordinates, cell, time_stamp, orientation=None, bonding=None, type=None):
        self.n_particles = n_particles
        self.n_frame = n_frame
        self.coordinates = coordinates
        self.orientation = orientation
        self.bonding = bonding
        self.type = type
        self.cell = cell
        self.time_stamp = time_stamp

        # the following are optional and can be set using methods below
        self.patches = None
        self.cluster_info = None
        self.bonds_calculated = None
        self.percolated = None

        # check that the provided data types make sense
        self.check_data()

    def set_cluster_info(self):
        # creates a cluster_info object and populates the necessary info
        if self.bonds_calculated is None:
            self.get_list_of_interacting_pairs()

        cl = clusters.cluster_info(self.bonds_calculated)
        self.cluster_info = cl

    def set_patch_info(self,file_name):
        # reads patch info from file and creates patches object 
        self.patches = patches.patches(file_name)

    def check_data(self):
        # do some checking to see if all data is the correct size and type
        if not isinstance(self.n_particles, int):
            print("WARNING: non-integer particle number")
        elif self.n_particles <= 0:
            print("WARNING: particle number <= 0")

        if not isinstance(self.n_frame, int):
            print("WARNING: non-integer frame number")
        elif self.n_frame < 0:
            print("WARNING: frame number < 0")

        if type(self.coordinates).__module__ != np.__name__:
            print("WARNING: coordinates must be contained in a numpy array")
        else:
            m, n = self.coordinates.shape
            if m != self.n_particles:
                print(
                    "WARNING: number of coordinates do not match the number of particles")
            if n != 3:
                print("WARNING: invalid number of columns for coordinates")

        if type(self.cell).__module__ != np.__name__:
            print("WARNING: cell info must be contained in a numpy array")
        else:
            m = self.cell.shape[0]
            if m != 6:
                print("WARNING: invalid cell array shape")

        if self.orientation is not None:
            if type(self.orientation).__module__ != np.__name__:
                print("WARNING: orientations must be contained in a numpy array")
            else:
                m, n = self.orientation.shape
                if m != self.n_particles:
                    print(
                        "WARNING: number of orientations do not match the number of particles")
                if n != 3:
                    print("WARNING: invalid number of columns for orientations")

        if self.bonding is not None:
            if type(self.orientation).__module__ != np.__name__:
                print("WARNING: bonding info must be contained in a numpy array")
            else:
                m = self.bonding.shape[0]
                if m != self.n_particles:
                    print(
                        "WARNING: number of bonding info do not match the number of particles")

        if self.type is not None:
            if type(self.type).__module__ != np.__name__:
                print("WARNING: type info must be contained in a numpy array")
            else:
                m = self.type.shape[0]
                if m != self.n_particles:
                    print(
                        "WARNING: number of bonditypeng info do not match the number of particles")

    def get_list_of_interacting_pairs(self):
        # returns a list of bonds. Requires the patches object to be initialized.
        if self.patches is None:
            sys.exit("ERROR: need to set patch info via frame.set_patch_info(file_name)")

        patch_obj = self.patches
        # TODO: make this more efficient with neighbor lists
        cell = self.cell

        # get the square of the maximum lambda value to use as a threshold distance
        max_lambda_sq = (np.max(patch_obj.lambda_vals))**2

        bonds = []
        # loop over pairs of particles (i,j) s.t. j>i
        for i in range(self.n_particles):
            for j in range(i+1, self.n_particles):
                pos_i = self.coordinates[i, :]
                pos_j = self.coordinates[j, :]

                # rij vector
                dist = pos_i - pos_j

                # periodic boundary conditions
                dist = utils.nearest_image(dist, cell)

                # squared distance between particles
                d2 = dist[0]**2 + dist[1]**2 + dist[2]**2

                # check if particles are interacting
                interacting = False
                # check if distance is less than the threshold value chosen above
                if d2 <= max_lambda_sq:
                    # types of particles
                    type_i = self.type[i]
                    type_j = self.type[j]

                    # loop over all patch pairs
                    # TODO: you can create this list above and avoid using the continue statements
                    for pi in range(patch_obj.n_patch):
                        if patch_obj.types[pi] != type_i:
                            continue
                        for pj in range(patch_obj.n_patch):
                            if patch_obj.types[pj] != type_j:
                                continue

                            # get orientations of particles
                            angles_i = self.orientation[i, :]
                            angles_j = self.orientation[j, :]
                            # get vector of the two patches
                            pi_v = patch_obj.patch_vectors[pi]
                            pj_v = patch_obj.patch_vectors[pj]
                            # get cos_delta values for the two patches
                            pi_cosdelta = patch_obj.cos_delta_vals[pi]
                            pj_cosdelta = patch_obj.cos_delta_vals[pj]
                            # check if the orientation is correct for interaction
                            correct_orientation = utils.check_reciprocal_interaction(
                                angles_i, angles_j, pi_v, pj_v, pi_cosdelta, pj_cosdelta, dist)
                            # if the orientation is correct
                            if correct_orientation:
                                # make sure the interacting patches have the same lambda value
                                assert(
                                    patch_obj.lambda_vals[pi] == patch_obj.lambda_vals[pj])
                                d = np.sqrt(d2)
                                # final check of distance
                                if d <= patch_obj.lambda_vals[pi]:
                                    # the particles are interacting
                                    interacting = True
                if interacting:
                    bonds.append((i, j))

        self.bonds_calculated = bonds

    def check_calculated_bonds_against_bond_numbers(self):
        # make sure that the calculated number of bonds per particle matches what's written in the trajectory file
        if self.bonds_calculated is None:
            self.get_list_of_interacting_pairs()
        if self.bonding is None:
            print("bonding info not available")
            return
        bond_numbers = np.zeros(self.n_particles)
        for bond in self.bonds_calculated:
            bond_numbers[bond[0]] += 1
            bond_numbers[bond[1]] += 1
            
        correct = True
        for i in range(self.n_particles):
            if bond_numbers[i] != self.bonding[i]:
                correct = False
                print("Particle i has %d bonds from calc, but traj file says %d" %
                  (bond_numbers[i], self.bonding[i]))

        if not correct:
            sys.exit("There's a problem with the calculated bond numbers")
        else:
            print("Bond counts match!")
        
    
    def calculate_rdf(self,binsize=0.1):
        # calculate the radial distribution function
        cell=self.cell
        # TODO: (low-prio) generalize the rdf function to a rectangular box

        # make sure we have a cubic box
        np.testing.assert_allclose(cell[:3],cell[0],rtol=1e-5)
        r,gr = rdf_sq.calculate_rdf(self.coordinates,self.n_particles,cell[0],binsize=binsize)
        return r,gr

    def calculate_sq(self,g=30):
        # calculate the structure factors
        cell = self.cell
        # make sure we have a cubic box
        # TODO: generalize to a rectangular box
        np.testing.assert_allclose(cell[:3],cell[0],rtol=1e-5)
        k,sk = rdf_sq.calculate_sq(self.n_particles,self.coordinates,cell[0],g=g)
        return k,sk
        
    def get_bond_probability(self):
        # calcluate the bond probability
        assert(self.patches is not None)
        patch_object = self.patches
        unique_types = list(set(patch_object.types))
        n_patches_per_type = np.zeros(len(unique_types))
        for utype, i in enumerate(unique_types):
            for type in patch_object.types:
                if utype == type:
                    n_patches_per_type[i] += 1
                    
        max_bonds = 0
        curr_bonds = 0
        for i in range(self.n_particles):
            type = self.type[i]
            bonds = self.bonding[i]
            max_bonds += n_patches_per_type[type]
            curr_bonds += bonds
            
        return curr_bonds/max_bonds

    def find_percolating_clusters(self,max_lambda=1.1):
        # find out which clusters percolate
        cell = self.cell
        grid_spacing = 0.1
        clusters = self.cluster_info.clusters

        percolated_clusters = np.zeros(len(clusters))
        
        nx = int(np.round(cell[0] / grid_spacing))
        ny = int(np.round(cell[1] / grid_spacing))
        nz = int(np.round(cell[2] / grid_spacing))
        
        grid_spacing_x = cell[0] / nx
        grid_spacing_y = cell[1] / ny
        grid_spacing_z = cell[2] / nz
        
        rad = int(np.round((max_lambda/2)/grid_spacing_x))
        
        valid_triples = []
        for ix in range(-rad, rad):
            for iy in range(-rad, rad):
                for iz in range(-rad, rad):
                    if ix**2 + iy**2 + iz**2 <= rad**2:
                        valid_triples.append((ix, iy, iz))
                        
        for ci,cluster in enumerate(clusters):
            grid = np.zeros((nx, ny, nz))
            for p in cluster:
                x, y, z = self.coordinates[p,0], self.coordinates[p,1], self.coordinates[p,2]
                sx, sy, sz = int(np.round(x / grid_spacing_x)), int(np.round(y / grid_spacing_y)), int(np.round(z / grid_spacing_z))
                
                for triple in valid_triples:
                    ixx = (sx + triple[0]) % nx
                    iyy = (sy + triple[1]) % ny
                    izz = (sz + triple[2]) % nz
                    
                    grid[ixx, iyy, izz] = 1
                    
                    
                xy_plane = np.sum(grid, axis=2)
                yz_plane = np.sum(grid, axis=0)
                
                # reduce dimensions to project connectivity onto axes and bin the data to overcome discretization errors
                binsize = 0.2
                nxs = int(cell[0] / binsize) + 1
                nys = int(cell[1] / binsize) + 1
                nzs = int(cell[2] / binsize) + 1
                
                x_axis = np.zeros(nxs)
                y_axis = np.zeros(nys)
                z_axis = np.zeros(nzs)
                for i in range(nx):
                    bin_ind = int(i*grid_spacing_x/binsize)
                if np.sum(xy_plane[i, :]) >= 1:
                    x_axis[bin_ind] += 1
                    
                for i in range(ny):
                    bin_ind = int(i*grid_spacing_y/binsize)
                if np.sum(xy_plane[:, i]) >= 1:
                    y_axis[bin_ind] += 1
                    
                for i in range(nz):
                    bin_ind = int(i*grid_spacing_z/binsize)
                if np.sum(yz_plane[:, i]) >= 1:
                    z_axis[bin_ind] += 1
                    
                # percolated if all values >= 1 in along at least one dimension
                if all(i >= 1 for i in x_axis) or all(i >= 1 for i in y_axis) or all(i >= 1 for i in z_axis):
                    percolated_clusters[ci] = 1

        self.cluster_info.percolated_clusters = percolated_clusters
        
    def is_system_percolated(self, max_lambda=1.1):
        # determine if the system percolates
        if self.cluster_info.percolated_clusters is None:
            self.find_percolating_clusters(self,max_lambda=max_lambda)
            
        if np.sum(self.cluster_info.percolated_clusters) >= 1:
            self.percolated = True
        else:
            self.percolated = False

    def write_frame_with_cluster_info(self,file_name):
        # write an xyz file such that each cluster is denoted by one of the following random atom types
        types = ['N','C','O','H','S','P','Na', 'Mg', 'Cl', 'Si']
        
        f = open(file_name, 'w')
        f.write("%d\n" % self.n_particles)
        f.write("Lx= %lf, Ly= %lf, Lz= %lf\n" % (self.cell[0], self.cell[1], self.cell[2]))
        for cluster in self.cluster_info.clusters:
            type = random.choice(types)
            for p in cluster:
                f.write(type)
                f.write(" %lf %lf %lf\n" % (self.coordinates[p, 0], self.coordinates[p, 1], self.coordinates[p, 2]))
        f.close()

    def get_cluster_rg(self):
        # calculate the R_g of clusters
        # (may want to exclude percolating ones when calculating stats)
        cell = self.cell
        rgs = []
        for i,cluster in enumerate(self.cluster_info.clusters):
            relevant_bonds = self.cluster_info.relevant_bonds[i]
            positions = utils.make_molecule_whole(self.coordinates,cell,relevant_bonds)
            rg = utils.calculate_rg(positions)
            rgs.append(rg)

        self.cluster_info.R_g = rgs

    def write_biggest_cluster_xyz(self,file_name):
        # get the biggest cluster and write it to an .xyz file
        max_size = 0
        relevant_bonds = []
        for i,cluster in enumerate(self.cluster_info.clusters):
            if len(cluster) > max_size:
                max_size = len(cluster)
                relevant_bonds = self.cluster_info.relevant_bonds[i]

        xyz = utils.make_molecule_whole(self.coordinates,self.cell,relevant_bonds)
        
        f = open(file_name, 'w')
        f.write("%d\n"%max_size)
        f.write("Lx= %lf, Ly= %lf, Lz= %lf\n" % (self.cell[0], self.cell[1], self.cell[2]))
        for pos in xyz:
            f.write('N ')
            f.write(" %lf %lf %lf\n" % (pos[0],pos[1],pos[2]))
        f.close()


class trajectory():
    def __init__(self, file_name=None, list_of_frames=None):
        if file_name is not None:
            self.read_trajectory(file_name)
        else:
            assert(list_of_frames is not None)
            self.list_of_frames = list_of_frames
        
        self.n_frames = len(self.list_of_frames)
        self.n_particles = self.list_of_frames[0].n_particles

        self.patches = None

        self.check_data()

    def check_data(self):
        if type(self.list_of_frames) is not list:
            sys.exit("ERROR: trajectory requires a list of frames")
        for snapshot in self.list_of_frames:
            if not isinstance(snapshot, frame):
                sys.exit("ERROR: trajectory must consist of frames")

    def get_frame(self, start, end=None):
        # get a portion of the trajectory
        # TODO: make sure start and end are within range
        if end is None or start == end:
            # TODO: should we return a trajectory with a single frame instead of a frame object?
            fr =  self.list_of_frames[start]
            fr.patches = self.patches
            return fr
        else:
            if start > end:
                dummy = end
                end = start
                start = dummy

        trj = trajectory(file_name=None, list_of_frames=self.list_of_frames[start:end])
        trj.patches = self.patches
        return trj

    def get_last_frame(self):
        # get the last frame
        frame =  self.list_of_frames[-1]
        frame.patches = self.patches
        return frame

    def write_xyz(self, file_name):
        # up to 5 different types of particles are currently supported
        types = {0: 'N', 1: 'C', 2: 'O', 3: 'H', 4: 'S'}

        # TODO: add a start and end option
        f = open(file_name, 'w')
        for frame_obj in self.list_of_frames:
            f.write("%d\n" % self.n_particles)
            f.write("frame = %d" % frame_obj.n_frame)
            f.write(", Lx= %lf, Ly= %lf, Lz= %lf\n" %
                    (frame_obj.cell[0], frame_obj.cell[1], frame_obj.cell[2]))

            for i in range(self.n_particles):
                f.write(types[frame_obj.type[i]])
                f.write(" %lf %lf %lf\n" % (
                    frame_obj.coordinates[i, 0], frame_obj.coordinates[i, 1], frame_obj.coordinates[i, 2]))

        f.close()

    def set_patch_info(self,file_name):
        # set patch info
        self.patches = patches.patches(file_name)

    def position_autocorrelation(self,start=None,dt=500,dx=1):
        # calculate autocorrelation function for positions
        if start is None:
            # if no start is provided, start from 1/3rd of the traj
            start = int(self.n_frames/3)

        times = []
        corr_vals = []
        for i,frame in enumerate(self.list_of_frames):
            xyz = frame.coordinates
            if i == start:
                time = 0
                xyz_first = np.copy(xyz)
            # skip frames before start
            elif i < start:
                continue

            time += dt
            corr = 0

            for j, atom in enumerate(xyz):
                # check how much the particle moved
                dist = xyz[j] - xyz_first[j]
                # take into account pbc
                dist = utils.nearest_image(dist,frame.cell)
                # calculate distance
                d = np.linalg.norm(dist)

                # if the distance is smaller than a cutoff distance
                # (typically chosen as 1 particle diameter)
                if d < dx:
                    # then decide that the particle has not moved away sufficiently
                    corr += 1.

            # get the fraction of immobile particles
            corr /= frame.n_particles

            times.append(time)
            corr_vals.append(corr)
        return times, corr_vals

    def read_trajectory(self,file_name):
        # read trajectory from file 

        # TODO: modify this to allow for reading on the last frame without loading the whole file
 
        # read first 10 lines of file to determine what info is contained
        with open(file_name) as f:
            head = [next(f) for x in range(10)]
            
        col_nums = [len(line.split()) for line in head]
        col_parts = col_nums[-1]
        
        if col_parts == 3:
            orientation = None
            bonding = None
            type = None
        elif col_parts == 6:
            bonding = None
            type = None
            orientation = []
        elif col_parts > 8:
            sys.exit(
            "Error reading trajectory file. Particle info has more than 8 columns")    
        elif col_parts < 3:
            sys.exit(
            "Error reading trajectory file. Particle info has fewer than 3 columns.")
        elif col_parts == 7:
            sys.exit("Error reading trajectory file. File has 7 columns.")
        else:
            bonding = []
            type = []
            orientation = []
            
        if col_nums[2] != 2:
            sys.exit(
            "Error reading trajectory file. The third line should have 2 columns.")
            
        # read the whole file
        with open(file_name, 'r') as f:
            lines = f.readlines()
            
        # get number of frames and particles
        frames = 0
        particles = 0
        for line in lines:
            num = len(line.split())
            if num == 2:
                frames += 1
            elif num == col_parts:
                particles += 1
                
        particles = int(particles/frames)
        frames = int(frames)
        
        # number of lines of info for each frame
        T = particles + 3
        
        # loop over frames
        traj = []
        for fr in range(frames):
            # index where coords start
            ind1 = T*fr + 3
            # index where coordinates end
            ind2 = T*(fr+1) - 1
            
            # get cell
            newCellLine = lines[T*fr].split()
            cell = np.array([float(newCellLine[0]), float(newCellLine[1]), float(newCellLine[2]), 90., 90., 90.])
            
            # get time_stamp
            time_stamp_line = lines[T*fr+1].split()
            time_stamp = int(time_stamp_line[2])
            
            data = []
            # loop through the lines that contain particle coordinate data
            for i in range(ind1, ind2+1):
                # format and add the coordinates into the data list
                data.append([float(l) for l in lines[i].split()])
                
            # convert data list to numpy array for easier manipulation
            data = np.array(data)
            
            # extract particle coordinates
            xyz = data[:, :3]
            
            # TODO: check column numbers instead and remove the list initializations above
            if orientation is not None:
                orientation = data[:, 3:6]
                
            if bonding is not None:
                bonding = data[:, 7].astype(int)
                
            if type is not None:
                type = data[:, 6].astype(int)
                
            frame_obj = frame(particles, fr, xyz, cell,
                          time_stamp, orientation, bonding, type)
            traj.append(frame_obj)

        self.list_of_frames = traj
