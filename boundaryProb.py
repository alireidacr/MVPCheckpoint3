# python script to solve Poisson equation for electrostatics

import sys
import random as rnd
import numpy as np

def getInputParams():
    dimensions = int(sys.argv[1])

    return dimensions

def generateInitState(dimensions):
    lattice = np.zeros((dimensions, dimensions, dimensions))
    for pos, state in np.ndenumerate(lattice):
            lattice[pos] = rnd.random()



    return lattice


def jacobiUpdate(pos, potentialLat, charge):
    intermediate = (potentialLat[(pos[0]+1, pos[1], pos[2])] + potentialLat[(pos[0]-1, pos[1], pos[2])]
                  + potentialLat[(pos[0], pos[1]+1, pos[2])] + potentialLat[(pos[0], pos[1]-1, pos[2])]
                  + potentialLat[(pos[0], pos[1], pos[2]+1)] + potentialLat[(pos[0], pos[1], pos[2]-1)])

    return (1/6)*(intermediate + charge)


def updateLattice(chargeLat, potentialLat, dimensions):
    returnLat = np.zeros((dimensions, dimensions, dimensions))

    for pos, state in np.ndenumerate(potentialLat):
        if (dimensions -1 in pos or 0 in pos):
            returnLat[pos] = 0
        else:
            returnLat[pos] = jacobiUpdate(pos, potentialLat, chargeLat[pos])

    potentialLat = returnLat

def writeElectricField(potentialLat, dimensions):
    eFile = open("electricField.txt", "a+")

    for pos, state in np.ndenumerate(potentialLat):
        if (0 in pos or dimensions-1 in pos):
            pass
            #eFile.write(str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(0) + " " + str(0) + " " + str(0) + "\n")
        else:
            xComp = (potentialLat[(pos[0]+1, pos[1], pos[2])] - potentialLat[(pos[0]-1, pos[1], pos[2])])/2
            yComp = (potentialLat[(pos[0], pos[1]+1,pos[2])] - potentialLat[(pos[0], pos[1]-1, pos[2])])/2
            zComp = (potentialLat[pos[0], pos[1], pos[2]+1] - potentialLat[pos[0], pos[1], pos[2]-1])/2

            eFile.write(str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(-xComp) + " " + str(-yComp) + " " + str(-zComp) + "\n")

    eFile.close()

def main():
    dimensions = getInputParams()

    chargeLat = np.zeros((dimensions, dimensions, dimensions))
    centre = int(dimensions/2)
    chargeLat[(centre, centre, centre)] = 1 # place a point charge at the centre of the box
    potentialLat = generateInitState(dimensions)

    for sweep in range(10000):
        updateLattice(chargeLat, potentialLat, dimensions)
        print(sweep)

    writeElectricField(potentialLat, dimensions)

main()
