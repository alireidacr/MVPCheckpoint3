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

    print("dimensions: " + str(dimensions))
    print("mobility: " + str(mobility))
    print("kappa: " + str(kappa))
    print("aConst: " + str(aConst))
    print("spatDisc:" + str(spatDisc))
    print("tempDisc: " + str(tempDisc))
    print("initVal: " + str(initVal))
    return dimensions, mobility, kappa, aConst, spatDisc, tempDisc, initVal

def generateInitState(dimensions, initVal):
    grid = np.zeros((dimensions, dimensions))
    for pos, state in np.ndenumerate(grid):
        grid[pos] = initVal + (rnd.random() - 0.5) * 2 * 0.01

    return grid

def findChemPotential(lattice, kappa, aConst, pos, spatDisc):
    NNs = lattice.getNearestNeighbours(pos)

    intermediate1 = aConst * (lattice[pos]**3 - lattice[pos])
    intermediate2 = sum([lattice[NN] for NN in NNs]) - (4.0 * lattice[pos])
    intermediate3 = ((kappa / spatDisc**2) * intermediate2)

    chemPotential = intermediate1 + intermediate3
    return chemPotential


def updateLattice(lattice, mobility, kappa, aConst, spatDisc, tempDisc):
    returnGrid = np.zeros((lattice.size(), lattice.size()))

    for pos, state in np.ndenumerate(lattice.getLattice()):
        chemPot = findChemPotential(lattice, kappa, aConst, pos, spatDisc)
        NNs = lattice.getNearestNeighbours(pos)

        mobCont = sum([findChemPotential(lattice, kappa, aConst, NN, spatDisc) for NN in NNs]) - (4 * chemPot)
        factor = (mobility * tempDisc)/(spatDisc**2)
        returnGrid[pos] = state + (factor * mobCont)

    lattice.setLattice(returnGrid)


def freeEnergy(lattice, aConst, kappa, spatDisc):
    intermediate1 = 0
    intermediate2 = 0
    intermediate3 = 0

    for pos, state in np.ndenumerate(lattice.getLattice()):
         intermediate1 += state**2
         intermediate2 += state**4
         NNs = lattice.getNearestNeighbours(pos)
         intermediate3 += (NNs[0]-NNs[2]+NNs[1]-NNs[3])

    freeEn = (aConst/2) * intermediate1 + (aConst/4) * intermdiate2 + (kappa/(4*spatDist)) * intermediate3

    return freeEn


def main():
    dimensions, mobility, kappa, aConst, spatDisc, tempDisc, initVal = getInputParams()

    initgrid = generateInitState(dimensions, initVal)
    lattice = PeriodicLattice(dimensions, 4, initgrid)

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
        #print(lattice.getLattice())

main()
