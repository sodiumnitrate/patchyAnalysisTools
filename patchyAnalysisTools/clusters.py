import copy
import networkx as nx
import pdb
import numpy as np

'''

This files defines a cluster_info class to hold cluster-related info for a given frame.
(This is then used as an attribute for the frame class)

'''

class cluster_info():
    # initialize. List of bonds need to be already calculated (see the frame class).
    def __init__(self,bonds):
        self.clusters = None
        self.bonds = bonds

        self.relevant_bonds = None

        self.N_loop = None
        self.loops = None
        self.L_chain = None
        self.chains = None
        self.R_g = None
        self.S_cluster = None

        self.percolated_clusters = None

        # populates the clusters attribute
        self.find_all_clusters()
        print("all clusters found")
        self.get_number_of_cycles() # sets N_loop
        print("N_loop calculated")
        self.get_chain_lengths()    # sets L_chain
        print("L_chain calculated")
        self.get_cluster_size()     # sets S_cluster
        print("S_cluster calculated")

    def find_all_clusters(self):
        # function to find the list of all clusters

        list_of_particles_in_bonds = []
        all_bonds = []
        for k in range(len(self.bonds)):
            for (i, j) in self.bonds[k]:
                list_of_particles_in_bonds.append(i)
                list_of_particles_in_bonds.append(j)
                all_bonds.append((i,j))

        list_of_particles_in_bonds = list(set(list_of_particles_in_bonds))
        all_bonds = list(set(all_bonds))

        clusters = []
        while len(list_of_particles_in_bonds) > 0:
            node = list_of_particles_in_bonds[0]
            selected = select_cluster(node, all_bonds)
            selected = list(set(selected))
            clusters.append(selected)
            for p in selected:
                list_of_particles_in_bonds.remove(p)

        self.clusters = clusters

    def get_relevant_bonds(self):
        # finds a list of bonds for each cluster
        relevant_bonds = []
        for cluster in self.clusters:
            bondlist = []
            for particle in cluster:
                for k in range(len(self.bonds)):
                    for bond in self.bonds[k]:
                        if particle in bond:
                            bondlist.append(bond)

            bondlist = list(set(bondlist))
            relevant_bonds.append(copy.copy(bondlist))

        self.relevant_bonds = relevant_bonds

    def get_cluster_size(self):
        # function to populate S_cluster
        S_cluster = []
        for cluster in self.clusters:
            S_cluster.append(len(cluster))

        self.S_cluster = S_cluster

    def get_number_of_cycles(self):
        # finds the number of (basis) cycles in each cluster
        if self.relevant_bonds is None:
            self.get_relevant_bonds()
        N_loop = []
        loops = []
        for i,cluster in enumerate(self.clusters):
            G = nx.Graph()
            G.add_edges_from(self.relevant_bonds[i])
            cycles = list(nx.cycle_basis(G))
            N_loop.append(len(cycles))
            loops.append(cycles)

        self.N_loop = N_loop
        self.loops = loops

    def get_chain_lengths(self):
        # calculates the chain length for each cluster
        # (defined as the longest path in a cluster)
        if self.relevant_bonds is None:
            self.get_relevant_bonds()

        L_chain = []
        chain = []
        for i,cluster in enumerate(self.clusters):
            G = nx.Graph()
            G.add_edges_from(self.relevant_bonds[i])
            
            p = nx.shortest_path(G)
            paths = []
            chains = []
            for i in cluster:
                for j in cluster:
                    length = len(p[i][j])
                    paths.append(length)
                    chains.append(p[i][j])

            ind = np.argmax(paths)
            L_chain.append(paths[ind])
            chain.append(chains[ind])

        self.L_chain = L_chain
        self.chains = chain


def select_cluster(node, bonds):
    # function to select clusters, given a list of bonds and a node
    # (flood-fill algorithm)
    Q = [node]
    selected = []
    while len(Q) > 0:
        curr = Q[0]
        Q = Q[1:]
        selected.append(curr)
        # TODO: use nx.adj here instead for a possible performance improvement
        for (i, j) in bonds:
            if curr == i and j not in selected:
                Q.append(j)
            elif curr == j and i not in selected:
                Q.append(i)
    return selected
