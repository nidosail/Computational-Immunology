import numpy as np
import scipy as sp
import copy
import random

# begin class code
class Th17cell:
    # some random sh!t that may come into use later, still tryna figure out the role of variables defined outside of
    # initialization method
    Thresh2 = 70  # threshold of il2 required for cell prolif
    Thresh6 = 100  # threshold of il6 reuired for cell proliferation
    Thresh7 = 100  # threshold of il23 required for cell prolif
    v23 = -3  # value that determines maximum enhancement to speed of prolif in the presence of  IL-23
    k23 = 100  # vale that determines sensitivite to IL-23
    size = 1  # number of voxels that this cell occupies

    # initialize variables, initializtion method, whatever u wanna call it. takes in position of form [x,y,z]
    def __init__(self, pos=[0, 0, 0], k17=100, v17=200, kgm=100, vgm=200, dil17=np.zeros(4), dgmcsf=np.zeros(4),
                 dblTmr=12, actTmr=0, dieTmr=36, divNum=0):
        self.k17 = k17  # michaelis menten half concentration to max rate blah blah blah constant (units arbitrary rn)
        self.v17 = v17  # michaelis menten max rate. this is a made up value, here to hold up the skeleton of a model
        self.kgm = kgm
        self.vgm = vgm
        self.pos = pos  # position within the 3d box specified in cellmats 3rd,4th,and 5th arrays. eek!
        self.dil17 = dil17  # initialize rate of il17 production for the current and next 3 timesteps. ie, dil17[0] would be the current rate of production, dil17[1] would be the rate of production one timestep in the future.
        self.dgmcsf = dgmcsf  # initialize rate of gmcsf production for the current and next 3 timesteps
        self.pos = pos  # store the cells position. I think this wont be necessary.
        self.dblTmr = dblTmr
        self.actTmr = actTmr
        self.dieTmr = dieTmr
        self.divNum = divNum

    def secrete(self, il6, il1b):
        # cytokine rate according to michaelis menten kinetics with switching functions, needs work.
        # currently both cytokines would have the same rate, will look into this later.

        # we are going for a spooky scary skeleton of a model rn
        self.il6 = il6
        self.il1b = il1b
        dil17_0 = self.dil17[0]  # This pulls the values of the rate of IL17 that will be secreted at the current timestep
        dgmcsf_0 = self.dgmcsf[0]  # This pulls the value of rate of GMCSF that will be secreted at the current timestep

        # this code updates the values stored for future secretion. I think we need to think more about this and maybe change the functions.
        if self.actTmr > 0:
            self.dil17 = self.dil17[1:] + [
                self.v17 * (self.il6 / (self.k17 + self.il6)) * self.v17 * (self.il1b / (self.k17 + self.il1b))]
            self.dgmcsf = self.dgmcsf[1:] + [
                self.vgm * (self.il6 / (self.kgm + self.il6)) * self.vgm * (self.il1b / (self.kgm + self.il1b))]
            self.actTmr = self.actTmr - 1
        else:
            self.dil17 = self.dil17[1:] + [0]
            self.dgmcsf = self.dil17[1:] + [0]
        # return the values of cytokines to be secreted
        return [dil17_0, dgmcsf_0]

    def dblOrDie(self, il6, il23=0, il2=0, il7=0):
        if self.dblTmr <= 0 and self.actTmr > 0 and divNum < 6 and (il6 >= self
        .Thresh6 or il2 >= self.Thresh2 or il7 >= self.Thresh7):
            self.dblTmr = 12 + self.v23 * il23 / (il23 + self.k23)
            self.dieTmr = 36
            self.divNum = self.divNum + 1
            return 'dbl'
        elif dieTmr == 0:
            return 'die'

# %%
# This will end up being the function that runs the actual simulation and kinda just wraps everything together.
def main():
    # First, greet the user. I want to make this a sort of user friendly interface to some degree,
    # so I'm starting out with that in mind. But, I am giving a hidden "skip" option in the beginning
    # for when we want to run a bunch of simulations, and it will basically skip straight to have you state how many
    # sims you want to run, and uploading a .txt file or something of the sort that contains your starting parameters

    if not input(
            'Hello, and welcome to the Agent-Based Cytokine-Driven model of Rheumatoid Arthritis, or "The ABCDs of RA". Press enter to begin.\n').lower() == 'skip':

        # This next section is again in place for user friendliness, but will mostly be in there for users who are want to be walked through how the model is initialized.
        if input(
                'Would you like to use the prompt interface to set initialization parameters, or would you prefer to use a .txt file to load parameters. For "prompt" enter "P" and for .txt enter "T"\n').lower() == 'p':
            # Load the default matrix size for now
            l = 100

            # Prompt asking if you want to change from the default lattice size
            if input(
                    'The default size is a square lattice of 100 voxel sides (1 mm^3). Would you like to change the default size? Y/N\n').lower() == 'y':
                l = int(input('Please put the integer number that represents a side of the lattice.\n'))

            # Prompts to get the starting number of various cell types.
            Th17_0 = int(input('Please input the number of Th17 cells that you would like to initialize:\n'))
            FLS_0 = int(input('Please input the number of FLS cells that you would like to initialize:\n'))

            # Piece to calculate the total number of cells that will be placed
            cellCount = Th17_0 + FLS_0

            # Here we prompt the user for the geometry of how they want the cells placed. We will likely mostly use either the None, Transwell, or Physiological
            distType = input('Please specifiy how you would like the cells to be distributed in space. Options are: None (N), Random (R), Transwell (T), or Physiological (P), Co-culture:\n')

            if distType.lower() == 't':
                Th17place = [1] * Th17_0 + [0] * ((l ** 2) - Th17_0)
                random.shuffle(Th17place) # this creates a vector of size 1000 populated with zeroes and ones randomly
                FLSplace = [1] * Th17_0 + [0] * ((l ** 2) - FLS_0)
                random.shuffle(FLSplace) # this creates a vector of size 1000 populated with zeroes and ones randomly
                Th17cellmat = np.zeros((l,l), dtype=Th17cell)
                FLScellmat = np.zeros((l,l), dtype=Th17cell)

                for i in range(int(l ** (1/2))):
                    Th17cellmat[i, :] = Th17place[int(i*l):int(i*l + (l))]
                    FLScellmat[i, :] = FLSplace[int(i*l):int(i*l + (l))]
                # for row in range():
                    for element in range(len(Th17cellmat[i, :])):
                        if Th17cellmat[i, element] == 1:
                            Th17cellmat[i, element] = Th17cell(pos = [i,element,0])
                        else:
                            Th17cellmat[i, element] = 0
                print(Th17cellmat[1, :])



#csmajor
#pythonfun
#chemewho?
#ihopedaveapproves
#hashtag
