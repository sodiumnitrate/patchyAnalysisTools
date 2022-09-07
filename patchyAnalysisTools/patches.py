from . import trajectory
from . import utils
import numpy as np
import pdb

class patches():
    # this assumes all patches interact with every other patch
    def __init__(self, file_name):
        self.n_patch = None
        self.eps_vals = None
        self.lambda_vals = None
        self.cos_delta_vals = None
        self.patch_vectors = None
        self.types = None

        self.adjacency = None

        self.read_patch_info(file_name)

    def read_patch_info(self, file_name):
        f = open(file_name, 'r')
        line = f.readline()         # npatch diameter pm_switch (labels)
        line = f.readline().split() # npatch diameter pm_switch 
        n_patch = int(line[0])
        pm_switch = int(line[-1])
        line = f.readline()         # ts  z[0]  z[1]  z[2]  (labels)
        line = f.readline()         # ts  z[0]  z[1]  z[2]  

        adjacency = np.zeros((n_patch,n_patch))
        if pm_switch == 0:
            adjacency += 1

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
            if len(line) >= 4:
                types.append(int(line[3]))
            else:
                types = None
            labels = f.readline()
            int_line = f.readline().strip()
            interacts_with = int(int_line)
            adjacency[i,interacts_with] = 1



        self.n_patch = n_patch
        self.eps_vals = eps_vals
        self.lambda_vals = lambda_vals
        self.cos_delta_vals = cos_delta_vals
        self.patch_vectors = patch_vectors
        self.types = types
        self.adjacency = adjacency