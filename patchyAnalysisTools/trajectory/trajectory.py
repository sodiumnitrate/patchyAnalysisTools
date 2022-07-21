import numpy as np
import sys
import pdb
import matplotlib.pyplot as plt
import random


class frame():
    def __init__(self, particles, frame, coordinates, cell, time_stamp, orientation=None, bonding=None, type=None):
        # TODO: change var names that are scalars to contain the word "number" or letter "n"
        self.particle = particles
        self.frame = frame
        self.coordinates = coordinates
        self.orientation = orientation
        self.bonding = bonding
        self.type = type
        self.cell = cell
        self.time_stamp = time_stamp

        # check that the provided data types make sense
        self.check_data()

    def check_data(self):
        if not isinstance(self.particle, int):
            print("WARNING: non-integer particle number")
        elif self.particle <= 0:
            print("WARNING: particle number <= 0")

        if not isinstance(self.frame, int):
            print("WARNING: non-integer frame number")
        elif self.frame < 0:
            print("WARNING: frame number < 0")

        if type(self.coordinates).__module__ != np.__name__:
            print("WARNING: coordinates must be contained in a numpy array")
        else:
            m, n = self.coordinates.shape
            if m != self.particle:
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
                if m != self.particle:
                    print(
                        "WARNING: number of orientations do not match the number of particles")
                if n != 3:
                    print("WARNING: invalid number of columns for orientations")

        if self.bonding is not None:
            if type(self.orientation).__module__ != np.__name__:
                print("WARNING: bonding info must be contained in a numpy array")
            else:
                m = self.bonding.shape[0]
                if m != self.particle:
                    print(
                        "WARNING: number of bonding info do not match the number of particles")

        if self.type is not None:
            if type(self.type).__module__ != np.__name__:
                print("WARNING: type info must be contained in a numpy array")
            else:
                m = self.type.shape[0]
                if m != self.particle:
                    print(
                        "WARNING: number of bonditypeng info do not match the number of particles")

    def get_list_of_interacting_pairs(self, patch_obj):
        # TODO: make this more efficient with neighbor lists
        cell = self.cell

        # get the square of the maximum lambda value to use as a threshold distance
        max_lambda_sq = (np.max(patch_obj.lambda_vals))**2

        bonds = []
        # loop over pairs of particles (i,j) s.t. j>i
        for i in range(self.particle):
            for j in range(i+1, self.particle):
                pos_i = self.coordinates[i, :]
                pos_j = self.coordinates[j, :]

                # rij vector
                dist = pos_i - pos_j

                # periodic boundary conditions
                dist = nearest_image(dist, cell)

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
                            correct_orientation = check_reciprocal_interaction(
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

        return bonds


def nearest_image(d, cell):
    for i in range(3):
        hbox = cell[i] / 2
        if d[i] > hbox:
            d[i] -= cell[i]
        elif d[i] < -hbox:
            d[i] += cell[i]

    return d


def check_calculated_bonds_against_bond_numbers(frame, bonds):
    bond_numbers = np.zeros(frame.particle)
    for bond in bonds:
        bond_numbers[bond[0]] += 1
        bond_numbers[bond[1]] += 1

    correct = True
    for i in range(frame.particle):
        if bond_numbers[i] != frame.bonding[i]:
            correct = False
            print("Particle i has %d bonds from calc, but traj file says %d" %
                  (bond_numbers[i], frame.bonding[i]))
    return correct


def get_bond_probability(frame, patch_object):
    unique_types = list(set(patch_object.types))
    n_patches_per_type = np.zeros(len(unique_types))
    for utype, i in enumerate(unique_types):
        for type in patch_object.types:
            if utype == type:
                n_patches_per_type[i] += 1

    max_bonds = 0
    curr_bonds = 0
    for i in range(frame.particle):
        type = frame.type[i]
        bonds = frame.bonding[i]
        max_bonds += n_patches_per_type[type]
        curr_bonds += bonds

    return curr_bonds/max_bonds


def check_reciprocal_interaction(angles_i, angles_j, pi_v, pj_v, pi_cosdelta, pj_cosdelta, rji):
    pi_v = rotate_vector(angles_i, pi_v)
    pj_v = rotate_vector(angles_j, pj_v)

    omega_i = pi_v.dot(-rji)/np.linalg.norm(rji)
    omega_j = pj_v.dot(rji)/np.linalg.norm(rji)

    if omega_i >= pi_cosdelta and omega_j >= pj_cosdelta:
        return True
    else:
        return False


def select_cluster(node, bonds):
    Q = [node]
    selected = []
    while len(Q) > 0:
        curr = Q[0]
        Q = Q[1:]
        selected.append(curr)
        for (i, j) in bonds:
            if curr == i and j not in selected:
                Q.append(j)
            elif curr == j and i not in selected:
                Q.append(i)
    return selected


def find_all_clusters(bonds):
    list_of_particles_in_bonds = []
    for (i, j) in bonds:
        list_of_particles_in_bonds.append(i)
        list_of_particles_in_bonds.append(j)

    list_of_particles_in_bonds = list(set(list_of_particles_in_bonds))

    clusters = []
    while len(list_of_particles_in_bonds) > 0:
        node = list_of_particles_in_bonds[0]
        selected = select_cluster(node, bonds)
        selected = list(set(selected))
        clusters.append(selected)
        for p in selected:
            list_of_particles_in_bonds.remove(p)
    return clusters


def plot_clusters(frame, clusters):
    n_clusters = len(clusters)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for cluster in clusters:
        # generate random color for each cluster
        color = "#"+''.join([random.choice('0123456789ABCDEF')
                            for i in range(6)])
        for i in cluster:
            x = frame.coordinates[i, 0]
            y = frame.coordinates[i, 1]
            z = frame.coordinates[i, 2]
            ax.scatter(x, y, z, c=color)

    plt.show()


def find_percolating_clusters(frame_obj, clusters, max_lambda=1.1):
    cell = frame_obj.cell
    grid_spacing = 0.1

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
            x, y, z = frame_obj.coordinates[p,0], frame_obj.coordinates[p,1], frame_obj.coordinates[p,2]
            sx, sy, sz = int(np.round(x / grid_spacing_x)), int(np.round(y / grid_spacing_y)), int(np.round(z / grid_spacing_z))

            for triple in valid_triples:
                ixx = (sx + triple[0]) % nx
                iyy = (sy + triple[1]) % ny
                izz = (sz + triple[2]) % nz

                grid[ixx, iyy, izz] = 1

        xy_plane = np.sum(grid, axis=2)
        yz_plane = np.sum(grid, axis=0)
        #xz_plane = np.sum(grid, axis=1)

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

    return percolated_clusters

def is_system_percolated(frame_obj, clusters, max_lambda=1.1):
    percolated_clusters = find_percolating_clusters(frame_obj,clusters,max_lambda=max_lambda)
    if np.sum(percolated_clusters) >= 1:
        return True
    else:
        return False


def rotate_vector(angles, vector):
    phi = angles[0]
    theta = angles[1]
    psi = angles[2]
    R = np.zeros((3, 3))
    R[0][0] = np.cos(phi) * np.cos(psi) - np.cos(theta) * \
        np.sin(phi) * np.sin(psi)
    R[0][1] = -np.sin(phi) * np.cos(theta) * np.cos(psi) - \
        np.cos(phi) * np.sin(psi)
    R[0][2] = np.sin(phi) * np.sin(theta)
    R[1][0] = np.cos(psi) * np.sin(phi) + np.cos(phi) * \
        np.cos(theta) * np.sin(psi)
    R[1][1] = np.cos(phi) * np.cos(theta) * np.cos(psi) - \
        np.sin(phi) * np.sin(psi)
    R[1][2] = -np.cos(phi) * np.sin(theta)
    R[2][0] = np.sin(psi) * np.sin(theta)
    R[2][1] = np.cos(psi) * np.sin(theta)
    R[2][2] = np.cos(theta)

    return R.dot(vector)


class trajectory():
    def __init__(self, list_of_frames, particles=None, frames=None):
        self.list_of_frames = list_of_frames

        self.check_data()

        if particles is None:
            self.particles = self.list_of_frames[0].particle
        else:
            self.particles = particles
            if self.particles != self.list_of_frames[0].particle:
                sys.exit("ERROR: particle number mismatch")

        if frames is None:
            self.frames = len(self.list_of_frames)
        else:
            self.frames = frames
            if self.frames != len(self.list_of_frames):
                sys.exit("ERROR: number of frames mismatch")

    def check_data(self):
        if type(self.list_of_frames) is not list:
            sys.exit("ERROR: trajectory requires a list of frames")
        for snapshot in self.list_of_frames:
            if not isinstance(snapshot, frame):
                sys.exit("ERROR: trajectory must consist of frames")

    def get_frame(self, start, end=None):
        # TODO: make sure start and end are within range
        if end is None or start == end:
            # TODO: should we return a trajectory with a single frame instead of a frame object?
            return self.list_of_frames[start]
        else:
            if start > end:
                dummy = end
                end = start
                start = dummy

        return trajectory(self.list_of_frames[start:end])

    def get_last_frame(self):
        return self.list_of_frames[-1]

    def write_xyz(self, file_name):
        # up to 5 different types of particles are currently supported
        types = {0: 'N', 1: 'C', 2: 'O', 3: 'H', 4: 'S'}

        # TODO: add a start and end option
        f = open(file_name, 'w')
        for frame_obj in self.list_of_frames:
            f.write("%d\n" % self.particles)
            f.write("frame = %d" % frame_obj.frame)
            f.write(", Lx= %lf, Ly= %lf, Lz= %lf" %
                    (frame_obj.cell[0], frame_obj.cell[1], frame_obj.cell[2]))

            for i in range(self.particles):
                f.write(types[frame_obj.type[i]])
                f.write(" %lf %lf %lf\n" % (
                    frame_obj.coordinates[i, 0], frame_obj.coordinates[i, 1], frame_obj.coordinates[i, 2]))


def read_trajectory(file_name):
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
        cell = np.array([float(newCellLine[0]), float(
            newCellLine[1]), float(newCellLine[2]), 90., 90., 90.])

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

    traj_obj = trajectory(traj, particles, frames)

    return traj_obj


def get_frame_slice(frame_object, grid_spacing=0.02, rad=25, start=0, end=2):
    cell = frame_object.cell

    if start < 0 or start > cell[0] or end < 0 or end > cell[0]:
        sys.exit(
            "Error: slice start=%d and end=%d values are problematic" % (start, end))

    if start > end:
        dummy = start
        start = end
        end = dummy

    nx = int(np.round(cell[0] / grid_spacing))
    ny = int(np.round(cell[1] / grid_spacing))
    nz = int(np.round(cell[2] / grid_spacing))

    grid_spacing_x = cell[0] / nx
    grid_spacing_y = cell[1] / ny
    grid_spacing_z = cell[2] / nz

    start /= grid_spacing_z
    end /= grid_spacing_z
    start = int(start)
    end = int(end)

    if start == end:
        end += 1

    if end > nz:
        end = nz

    grid = np.zeros((nx, ny, nz))

    valid_triples = []
    for ix in range(-rad, rad):
        for iy in range(-rad, rad):
            for iz in range(-rad, rad):
                if ix**2 + iy**2 + iz**2 <= rad**2:
                    valid_triples.append((ix, iy, iz))

    for particle in frame_object.coordinates:
        x, y, z = particle[0], particle[1], particle[2]
        sx, sy, sz = int(np.round(x/grid_spacing_x)), int(np.round(y /
                                                                   grid_spacing_y)), int(np.round(z/grid_spacing_z))

        for triple in valid_triples:
            ixx = (sx + triple[0]) % nx
            iyy = (sy + triple[1]) % ny
            izz = (sz + triple[2]) % nz

            grid[ixx, iyy, izz] = 1

    slice = np.sum(grid[:, :, start:end], axis=2)
    slice[slice > 1] = 1

    return slice


class patches():
    # this assumes all patches interact with every other patch
    def __init__(self, n_patch, eps_vals, lambda_vals, cos_delta_vals, patch_vectors, types):
        self.n_patch = n_patch
        self.eps_vals = eps_vals
        self.lambda_vals = lambda_vals
        self.cos_delta_vals = cos_delta_vals
        self.patch_vectors = patch_vectors
        self.types = types


def read_patch_info(patch_file_name):
    f = open(patch_file_name, 'r')
    line = f.readline()             # npatch diameter pm_switch (labels)
    line = f.readline().split()     # npatch diameter pm_switch
    n_patch = int(line[0])
    line = f.readline()             # ts  z[0]    z[1]    z[2] (labels)
    line = f.readline()  # ts  z[0]    z[1]    z[2]

    eps_vals = []
    lambda_vals = []
    cos_delta_vals = []
    patch_vectors = []
    types = []
    for i in range(n_patch):
        labels = f.readline()
        line = f.readline().split()
        eps_vals.append(float(line[0]))
        lambda_vals.append(float(line[1]))
        cos_delta_vals.append(float(line[2]))

        labels = f.readline()
        line = f.readline().split()
        vector = [float(line[i]) for i in range(3)]
        patch_vectors.append(np.array(vector))
        types.append(int(line[3]))
        labels = f.readline()
        labels = f.readline()

    return patches(n_patch, eps_vals, lambda_vals, cos_delta_vals, patch_vectors, types)
