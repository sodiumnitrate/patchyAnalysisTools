import numpy as np
from . import utils
from . import trajectory
import sq
import RDF as rdfcpp

def calculate_rdf(points, N, L, binsize=0.1):
    '''
    Function to calculate the radial distribution function, given
    N - number of particles
    L - (cubic) box size
    binsize - bin size (default=0.1)
    '''
    # distance cutoff for calculating r(g)
    cutoff = np.sqrt(3)/2 * L

    # number of bins
    nbins = int(cutoff/binsize)+1
    # re-adjust bin size based on nbins
    u = np.arange(0,nbins,1)
    # r vector of g(r)
    r = binsize*(u+0.5)
    # calculate density
    rho = float(N)/(L**3)

    # x, y, and z coordinates
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]

    # calculate gr
    gr = rdfcpp.RDF_virtualcopies(x, y, z, L, cutoff, binsize, nbins, rho, N)
    return r[:-1], gr[:-1]

def calculate_sq(n_particles, points, L, g=30):
    '''
    Function to calculate the 1d structure factor, given
    n_particles - number of particles
    points - coordinates
    g - max num of bins (default = 30)
    '''
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]

    points = np.array([x,y,z])

    # get bin size
    binsize = 2*np.pi/L
    # max k vector
    kmax = binsize * np.sqrt(3) * g
    # number of kbins
    kbins = int(kmax/binsize)+1
    # create k vector
    k = np.linspace(0+0.5*binsize,kbins*binsize+0.5*binsize,kbins)
    # calculate s(k)
    sk = sq.sq(g,L,points,n_particles)

    return k[1:], sk[1:kbins]