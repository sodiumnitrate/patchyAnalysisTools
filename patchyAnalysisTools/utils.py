"""
This file holds a bunch of independent utility functions.
"""
from . import trajectory
import numpy as np
import matplotlib.pyplot as plt
import sys
import random
import networkx as nx

def get_cluster_rg(frame_obj,clusters):
    # calculate the radius of gyration of clusters
    cell = frame_obj.cell
    rgs = []
    for cluster in clusters:
        positions = []
        for i in cluster:
            x = frame_obj.coordinates[i,0]
            y = frame_obj.coordinates[i,1]
            z = frame_obj.coordinates[i,2]

            positions.append([x,y,z])

        positions = np.array(positions)
        positions = make_molecule_whole(positions,cell)
        rg = calculate_rg(positions)
        rgs.append(rg)

    return rgs

def rotate_vector(angles, vector):
    # rotate a vector given Euler angles
    # ZXZ counterclockwise (?)
    # (note that this mirrors the rotation scheme in the simulation code)
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

def check_reciprocal_interaction(angles_i, angles_j, pi_v, pj_v, pi_cosdelta, pj_cosdelta, rji,tol=1e-6):
    # function to check the (angular) reciprocal interaction of two given particles
    pi_v = rotate_vector(angles_i, pi_v)
    pj_v = rotate_vector(angles_j, pj_v)

    omega_i = pi_v.dot(-rji)/np.linalg.norm(rji)
    omega_j = pj_v.dot(rji)/np.linalg.norm(rji)

    if omega_i+tol >= pi_cosdelta and omega_j+tol >= pj_cosdelta:
        return True
    else:
        return False


def plot_clusters(frame):
    # function to plot clusters in 3d, using matplotlib
    # NOTE: very resource intensive. Use the print .xyz methods and use 
    # a visulization program like VMD or OVITO instead.

    assert(frame.cluster_info is not None)
    clusters = frame.cluster_info.clusters

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


def calculate_rg(xyz):
    # calculate the radius of gyration of the given coordinates
    N,d = xyz.shape
    com = np.sum(xyz,axis=0)
    com /= N
    xyz -= com
    
    rg = 0.
    for pos in xyz:
        rg += pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2

    rg /= N
    return np.sqrt(rg)

def nearest_image(d, cell):
    # apply periodic boundary conditions
    for i in range(3):
        hbox = cell[i] / 2

        while d[i] > hbox:
            d[i] -= cell[i]

        while d[i] < -hbox:
            d[i] += cell[i]

    return d

def make_molecule_whole(xyz,cell,relevant_bonds):
    # unfragment a cluster that has been broken up due to periodic boundary conditions

    G = nx.Graph()
    G.add_edges_from(relevant_bonds)
    pos = dict(enumerate(xyz))
    
    dx = (cell[0] / 2) * any( 2 * abs(pos[u][0] - pos[v][0]) > cell[0] for u,v in relevant_bonds)
    dy = (cell[1] / 2) * any( 2 * abs(pos[u][1] - pos[v][1]) > cell[1] for u,v in relevant_bonds)
    dz = (cell[2] / 2) * any( 2 * abs(pos[u][2] - pos[v][2]) > cell[2] for u,v in relevant_bonds)

    new_pos = [[(x+dx)%cell[0], (y+dy)%cell[1], (z+dz)%cell[2]] for v,(x,y,z) in pos.items() if v in list(G.nodes)]

    return np.array(new_pos)


def autocorrelation(start,n_frames,A):
    '''
    Calculates autocorrelation for an observable A
    '''

    mean_A = np.mean(A[start:])
    autocorr = np.zeros(n_frames-start)
    for f in range(1,n_frames-start,1):
        integrand = 0
        count = 0
        # loop over all possible starting points
        for ff in range(start+1,n_frames-f):
            count += 1
            # calculate delta A at the starting point
            dAt = A[ff] - mean_A

            # 
            dAtt = A[ff:ff+f+1] - mean_A

            # sum over the product
            integrand += np.sum(dAt*dAtt)

        if count != 0:
            # normalize so that we have a proper average
            integrand /= count * f
            autocorr[f] = integrand

    if autocorr[1] == 0:
        return None
    else:
        return autocorr[1:] / autocorr[1]


def get_frame_slice(frame_object, grid_spacing=0.02, rad=25, start=0, end=2):
    # function to select a chosen slab of the particles
    # (for visualization purposes)

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

def generate_valid_triples(rad,grid_spacing):
    # select grid points within +/- radius (of a chosen point)
    radius = int(np.round((rad/2)/grid_spacing))
    valid_triples = []
    for ix in range(-radius,radius):
        for iy in range(-radius,radius):
            for iz in range(-radius,radius):
                if ix**2 + iy**2 + iz**2 <= radius**2:
                    valid_triples.append((ix,iy,iz))
    return valid_triples
