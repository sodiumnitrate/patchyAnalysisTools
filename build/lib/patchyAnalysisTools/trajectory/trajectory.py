import numpy as np
import sys 
import pdb

def test_add(x,y):
    return x+y

class frame():
    def __init__(self,particles,frame,coordinates,cell,orientation=None,bonding=None,type=None):
        #TODO: change var names that are scalars to contain the word "number" or letter "n"
        self.particle = particles
        self.frame = frame
        self.coordinates = coordinates
        self.orientation = orientation
        self.bonding = bonding
        self.type = type
        self.cell = cell

        # check that the provided data types make sense
        self.check_data()

    def check_data(self):
        if not isinstance(self.particle,int): 
            print("WARNING: non-integer particle number")
        elif self.particle <= 0:
            print("WARNING: particle number <= 0")

        if not isinstance(self.frame,int):
            print("WARNING: non-integer frame number")
        elif self.frame < 0:
            print("WARNING: frame number < 0")
        
        if type(self.coordinates).__module__ != np.__name__:
            print("WARNING: coordinates must be contained in a numpy array")
        else:
            m,n = self.coordinates.shape
            if m != self.particle:
                print("WARNING: number of coordinates do not match the number of particles")
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
                m,n = self.orientation.shape
                if m != self.particle:
                    print("WARNING: number of orientations do not match the number of particles")
                if n != 3:
                    print("WARNING: invalid number of columns for orientations")

        if self.bonding is not None:
            if type(self.orientation).__module__ != np.__name__:
                print("WARNING: bonding info must be contained in a numpy array")
            else:
                m = self.bonding.shape[0]
                if m != self.particle:
                    print("WARNING: number of bonding info do not match the number of particles")

        if self.type is not None:
            if type(self.type).__module__ != np.__name__:
                print("WARNING: type info must be contained in a numpy array")
            else:
                m = self.type.shape[0]
                if m != self.particle:
                    print("WARNING: number of bonditypeng info do not match the number of particles")

class trajectory():
    def __init__(self,list_of_frames,particles=None,frames=None):
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
            if not isinstance(snapshot,frame):
                sys.exit("ERROR: trajectory must consist of frames")

    def get_frame(self,start,end=None):
        #TODO: make sure start and end are within range
        if end is None or start == end:
            #TODO: should we return a trajectory with a single frame instead of a frame object?
            return self.list_of_frames[start]
        else:
            if start > end:
                dummy = end
                end = start
                start = dummy

        return trajectory(self.list_of_frames[start:end])

    def write_xyz(self,file_name):
        # up to 5 different types of particles are currently supported
        types = {0:'N', 1:'C', 2:'O', 3:'H', 4:'S'}

        #TODO: add a start and end option
        f = open(file_name,'w')
        for frame_obj in self.list_of_frames:
            f.write("%d\n" % self.particles)
            f.write("frame = %d" % frame_obj.frame)
            f.write(", Lx= %lf, Ly= %lf, Lz= %lf" % (frame_obj.cell[0],frame_obj.cell[1],frame_obj.cell[2]))

            for i in range(self.particles):
                f.write(types[frame_obj.type[i]])
                f.write(" %lf %lf %lf\n"%(frame_obj.coordinates[i,0],frame_obj.coordinates[i,1],frame_obj.coordinates[i,2]))



def read_trajectory(file_name):
    #TODO: modify this to allow for reading on the last frame without loading the whole file

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
        sys.exit("Error reading trajectory file. Particle info has more than 8 columns")
    elif col_parts < 3:
        sys.exit("Error reading trajectory file. Particle info has fewer than 3 columns.")
    elif col_parts == 7:
        sys.exit("Error reading trajectory file. File has 7 columns.")
    else:
        bonding = []
        type = []
        orientation = []

    if col_nums[2] != 2:
        sys.exit("Error reading trajectory file. The third line should have 2 columns.")

    # read the whole file
    with open(file_name,'r') as f:
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

        data = []

        # loop through the lines that contain particle coordinate data
        for i in range(ind1, ind2+1):
            # format and add the coordinates into the data list
            data.append([float(l) for l in lines[i].split()])

        # convert data list to numpy array for easier manipulation
        data = np.array(data)

        # extract particle coordinates
        xyz = data[:,:3]

        #TODO: check column numbers instead and remove the list initializations above
        if orientation is not None:
            orientation = data[:,3:6]

        if bonding is not None:
            bonding = data[:,6]

        if type is not None:
            type = data[:,7]

        frame_obj = frame(particles,fr,xyz,cell,orientation,bonding,type)
        traj.append(frame_obj)

    traj_obj = trajectory(traj,particles,frames)

    return traj_obj