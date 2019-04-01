# python script to simulate cahn-hilliard equation on periodic lattice

from PeriodicLattice import PeriodicLattice
import sys
import numpy as np
import random as rnd
import matplotlib.pyplot as plt

def getInputParams():
    dimensions = int(sys.argv[1])
    mobility = float(sys.argv[2])
    kappa = float(sys.argv[3])
    aConst = float(sys.argv[4])
    spatDisc = float(sys.argv[5])
    tempDisc = float(sys.argv[6])
    initVal = float(sys.argv[7])

    return dimensions, mobility, kappa, aConst, spatDisc, tempDisc, initVal

def generateInitState(dimensions, initVal):
    grid = np.zeros((dimensions, dimensions))
    for pos, state in np.ndenumerate(grid):
        grid[pos] = initVal + rnd.uniform(-0.1, 0.1)

    return grid

def findChemPotential(lattice, kappa, aConst, pos, spatDisc):
    NNs = lattice.getNearestNeighbours(pos)

    chemPotential = -(aConst * lattice[pos]) + (aConst * lattice[pos]**3)
    kappaCont = sum([lattice[NN] for NN in NNs]) - (4 * lattice[pos])

    chemPotential -= (kappa / spatDisc**2) * kappaCont
    return chemPotential



def updateLattice(lattice, mobility, kappa, aConst, spatDisc, tempDisc):
    returnGrid = np.zeros((lattice.size(), lattice.size()))

    for pos, state in np.ndenumerate(lattice.getLattice()):
        chemPot = findChemPotential(lattice, kappa, aConst, pos, spatDisc)

        state = lattice[pos]
        NNs = lattice.getNearestNeighbours(pos)
        mobCont =  sum([findChemPotential(lattice, kappa, aConst, NN, spatDisc) for NN in NNs]) - (4 * chemPot)
        returnGrid[pos] = state + ((mobility * tempDisc)/(spatDisc**2) * mobCont)

    lattice.setLattice(returnGrid);


def main():
    dimensions, mobility, kappa, aConst, spatDisc, tempDisc, initVal = getInputParams()

    initgrid = generateInitState(dimensions, initVal)
    lattice = PeriodicLattice(dimensions, 8, initgrid)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    im = ax.imshow(lattice.getLattice(), cmap='cool', vmin=-1, vmax=1)
    cbar = plt.colorbar(im, ticks=[-1, 0, 1])
    #cbar.ax.set_yticklabels(['-1', '0', '1'])
    plt.show(block=False)

    for sweep in range(1000):
        updateLattice(lattice, mobility, kappa, aConst, spatDisc, tempDisc)

        im.set_array(lattice.getLattice())
        fig.canvas.draw()
        print(lattice.getLattice())

main()
