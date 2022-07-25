from . import trajectory
import numpy as np

def get_cluster_rg(frame_obj,clusters):
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


def make_molecule_whole(xyz,cell):
    reference_particle = xyz[0,:]
    new_xyz = np.zeros(xyz.shape)
    for i,pos in enumerate(xyz):
        dist = pos - reference_particle
        dist = nearest_image(dist,cell)
        pos = dist + reference_particle
        new_xyz[i,:] = pos

    return new_xyz


def calculate_rg(xyz):
    N,d = xyz.shape
    com = np.sum(xyz,axis=0)
    com /= N
    xyz -= com
    
    rg = 0.
    for pos in xyz:
        rg += pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2

    rg /= N
    return rg

def nearest_image(d, cell):
    for i in range(3):
        hbox = cell[i] / 2
        if d[i] > hbox:
            d[i] -= cell[i]
        elif d[i] < -hbox:
            d[i] += cell[i]

    return d

def make_molecule_whole(xyz,cell):
    # TODO: there's a bug in this function somewhere -- fix it
    reference_particle = xyz[0,:]
    new_xyz = np.zeros(xyz.shape)
    for i,pos in enumerate(xyz):
        dist = pos - reference_particle
        dist = nearest_image(dist,cell)
        pos = dist + reference_particle
        new_xyz[i,:] = pos

    return new_xyz

def calculate_rg(xyz):
    N,_ = xyz.shape
    com = np.sum(xyz,axis=0)
    com /= N
    xyz -= com
    
    rg = 0.
    for pos in xyz:
        rg += pos[0] ** 2 + pos[1] ** 2 + pos[2] ** 2

    rg /= N
    return rg