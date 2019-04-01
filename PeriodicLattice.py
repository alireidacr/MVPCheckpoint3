# class to represent a lattice with periodic boundary conditions

import numpy as np
import random as rnd

class PeriodicLattice(object):

    def __init__(self, dimensions, NNType, initGrid):
        self.grid = initGrid
        self.dimensions = dimensions
        self.NNType = NNType # stores whether to use 4 or 8 nearest neighbours
        if (NNType == 4):
            self.NNFunc = self.getNNs4
        else:
            self.NNFunc = self.getNNs8



    def __getitem__(self, pos):
        return self.grid[pos]

    def __setitem__(self, pos, val):
        self.grid[pos] = val

    def size(self):
        return (self.dimensions)

    def getLattice(self):
        return self.grid

    def setLattice(self, newLattice):
        self.grid = newLattice

    def getNearestNeighbours(self, pos):
        NNs = self.NNFunc(pos)

        for i, pair in enumerate(NNs):
            for j, var in enumerate(pair):
                if (var == -1):
                    NNs[i][j] = self.dimensions -1
                elif (var == self.dimensions):
                    NNs[i][j] = 0

        # iterate over NNs to convert pairs to immutable tuples
        returnTup = [tuple(NN) for NN in NNs]
        return returnTup

    def getNNs8(self, pos):

        NNs = [] # list to contian lists of nearest neighbours coordinates

        NNs.append([pos[0], pos[1] +1])
        NNs.append([pos[0] +1, pos[1] +1])
        NNs.append([pos[0] +1, pos[1]])
        NNs.append([pos[0] +1, pos[1] -1])
        NNs.append([pos[0], pos[1] -1])
        NNs.append([pos[0] -1, pos[1] -1])
        NNs.append([pos[0] -1, pos[1]])
        NNs.append([pos[0] -1, pos[1] +1])

        return NNs

    def getNNs4(self, pos):
        NNs = [] # list to contain tuples of nearest neighbours coordinates

        NNs.append([pos[0], pos[1] +1])
        NNs.append([pos[0] +1, pos[1]])
        NNs.append([pos[0], pos[1] -1])
        NNs.append([pos[0] -1, pos[1]])

        return NNs
