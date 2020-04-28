import numpy as np
import scipy as sp
import copy

# timestep vector (for experimental purposes for now) 48 hours? minutes? seconds? tacos? im hungry for tacos
time = np.arange(48)

# begin class code
class Th17cell:
    # some random sh!t that may come into use later, still tryna figure out the role of variables defined outside of
    # initialization method
    growth_factor = 1
    radius = 1
    # initialize variables, initializtion method, whatever u wanna call it. takes in position of form [x,y,z]
    def __init__(self, pos):
        self.kmil17 = 100 # michaelis menten half concentration to max rate blah blah blah constant (units arbitrary rn)
        self.maxil17 = 200 # michaelis menten max rate. this is a made up value, here to hold up the skeleton of a model
        self.pos = pos # position within the 3d box specified in cellmats 3rd,4th,and 5th arrays. eek!
        self.dil17 = 0 # nitialize rate of il17, hence the d
        self.dgmcsf = 0 # initialize rate of gmcsf
        self.xpos = int(self.pos[0]) # making locations into ints for use in matrix indexing
        self.ypos = int(self.pos[1])
        self.zpos = int(self.pos[2])
        # below is legit work with the 5D array
        # the MAGICAL 6d cell matrix, location is within a 2x2x2 cube space just for shits, irrelevant tbh
        # backbone of the model, everything comes back to this 5dimensional thing which i can barely conceptualize
        # (cytokine, concentration, x, y, z, time). yuh
        self.cellmat = np.zeros((4,len(time),2,2,2,len(time)))
        # this series of for loops assigns the time vector to be the "dimension of the 6D array. 6d means 6 dimensions
        # thats a lot of dimensions
        for cytokine in range(4):
            for concentration in range(len(time)):
                for xloc in range(2):
                    for yloc in range(2):
                        for zloc in range(2):
                            self.cellmat[cytokine,concentration,xloc,yloc,zloc,:] = time
    # sense cytokines, assign to position in array. for this we assume inital conditions, t=0
    def sense(self, il6, il1b):
        # throughout this code, i refer to the position input, because that is the location within the box thingy
        # where the cell exists and senses and secretes cytokines
        self.il1b = float(il1b)
        self.il6 = float(il6)
        self.cellmat[0, :, self.xpos, self.ypos, self.zpos, 0] = self.il1b * np.ones(len(time))
        self.cellmat[1, :,  self.xpos, self.ypos, self.zpos, 0] = self.il6 * np.ones(len(time))
    # secrete cytokines, assign to values in matrix. matrixes are cool
    # future work: defining net rate of il1b and il6. will require programming of fibrolast class
    # will incorporate a shifting time value based on when the cell senses the cytokines, but not now. It's 2:18am
    # and i need to sleep. Also i just turned 20 like 2 hours ago! yay
    def secrete(self):
        # cytokine rate according to michaelis menten kinetics with switching functions, needs work.
        # currently both cytokines would have the same rate, will look into this later.
        # we are going for a spooky scary skeleton of a model rn
        self.dil17 = self.maxil17*(self.il6/(self.kmil17 + self.il6)) * self.maxil17*(self.il1b/(self.kmil17 + self.il1b))
        self.dgmcsf = self.maxil17*(self.il6/(self.kmil17 + self.il6)) * self.maxil17*(self.il1b/(self.kmil17 + self.il1b))
        # this section exists for graphing purposes, storing the progressive il17/gmcsf  values in their
        # respective vectors. by vectors i mean a vector within the 5D array that changes with the time values
        # python calls them lists. WE the people of matlab will never refer to a VECTOR as a list.
        # it is so much more than a list
        for x in range(len(time)):
            self.cellmat[2, :, self.xpos, self.ypos, self.zpos, x] = self.dil17*x*np.ones(len(time))  # storing il17 values at given times in the matrix
            self.cellmat[3,:, self.xpos, self.ypos, self.zpos, x] = self.dgmcsf*x*np.ones(len(time))   # same, but for gmc-sf

# our first th17! so cute
cellt = Th17cell([0,0,0])
# cells sense at time 0 and secrete cytokines from time 1-48
cellt.sense(4.0,4.0)
cellt.secrete()

# for now, cytokine concentrations at a given time are stored in a length 48 vector wtih identical entries.
# Concentration changes wieh time:
# This would print Il17 concentration at timestep 47
print(cellt.cellmat[2,:,0,0,0,47])

#csmajor
#pythonfun
#chemewho?
#ihopedaveapproves
#hashtag