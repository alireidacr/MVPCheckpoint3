# python script to solve Poisson equation for electrostatics

import sys
import random as rnd
import numpy as np
import math

def getInputParams():
    dimensions = int(sys.argv[1])
    accuracy = float(sys.argv[2])
    type = sys.argv[3]

    return dimensions, accuracy, type

def generateInitState(dimensions):
    lattice = np.zeros((dimensions, dimensions, dimensions))
    for pos, state in np.ndenumerate(lattice):
            lattice[pos] = rnd.random()



    return lattice


def jacobiHelper(pos, potentialLat, charge):
    intermediate = (potentialLat[(pos[0]+1, pos[1], pos[2])] + potentialLat[(pos[0]-1, pos[1], pos[2])]
                  + potentialLat[(pos[0], pos[1]+1, pos[2])] + potentialLat[(pos[0], pos[1]-1, pos[2])]
                  + potentialLat[(pos[0], pos[1], pos[2]+1)] + potentialLat[(pos[0], pos[1], pos[2]-1)])

    return (1/6)*(intermediate + charge)


def jacobiUpdate(chargeLat, potentialLat, dimensions):
    returnLat = np.zeros((dimensions, dimensions, dimensions))

    for pos, state in np.ndenumerate(potentialLat):
        if (dimensions -1 in pos or 0 in pos):
            returnLat[pos] = 0
        else:
            returnLat[pos] = jacobiHelper(pos, potentialLat, chargeLat[pos])

    return returnLat


def GSUpdate(chargeLat, potentialLat, dimensions):
    returnLat = np.copy(potentialLat)

    for pos, state in np.ndenumerate(returnLat):
        if (dimensions -1 in pos or 0 in pos):
            returnLat[pos] = 0
        else:
            returnLat[pos] = jacobiHelper(pos, potentialLat, chargeLat[pos])

    return returnLat

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

def writePotential(potentialLat, centre):
    potFile = open("Potential.txt", "a+")
    potSlice = open("PotentialSlice.txt", "a+")

    for pos, state in np.ndenumerate(potentialLat):
        potFile.write(str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(state) + "\n")

        if (pos[2] == centre):
            potSlice.write(str(pos[0]) + " " + str(pos[1]) + " " + str(state) + "\n")

    potFile.close()
    potSlice.close()

def writeEFieldSlice(potentialLat, dimensions, centre):
    eSliceFile = open("ESlice.txt", "a+")
    eMagFile = open("EMag.txt", "a+")

    for pos, state in np.ndenumerate(potentialLat):
        if (pos[1] == centre and not 0 in pos and not dimensions-1 in pos):
            xComp = (potentialLat[(pos[0]+1, pos[1], pos[2])] - potentialLat[(pos[0]-1, pos[1], pos[2])])/2
            zComp = (potentialLat[pos[0], pos[1], pos[2]+1] - potentialLat[pos[0], pos[1], pos[2]-1])/2
            eSliceFile.write(str(pos[0]) + " " + str(pos[2]) + " " + str(-xComp) + " " + str(-zComp) + "\n")

            rad = math.sqrt((pos[0] - centre)**2 + (pos[2] - centre)**2)
            eMag = math.sqrt(xComp**2 + zComp**2)
            eMagFile.write(str(rad) + " " + str(eMag) + "\n")

    eSliceFile.close()
    eMagFile.close()


def findMaxDif(oldLat, newLat):
    maxDiff = 0
    for pos, state in np.ndenumerate(oldLat):
        if (abs(state - newLat[pos]) > maxDiff):
            maxDiff = abs(state - newLat[pos])

    return maxDiff

def writeMagField(potentialLat, dimensions, centre):
    magFile = open("magField.txt", "a+")
    magSliceFile = open("magFieldSlice.txt", "a+")

    for pos, state in np.ndenumerate(potentialLat):
        if (not 0 in pos and not dimensions-1 in pos):
            xComp = (potentialLat[(pos[0], pos[1]+1, pos[2])] - potentialLat[(pos[0], pos[1]-1, pos[2])])
            yComp = (potentialLat[(pos[0]+1, pos[1], pos[2])] - potentialLat[(pos[0]-1, pos[1], pos[2])])
            zComp = 0

            magFile.write(str(pos[0]) + " " + str(pos[1]) + " " + str(pos[2]) + " " + str(xComp) + " " + str(yComp) + " " + str(zComp) + "\n")
            if (pos[2] == centre):
                magSliceFile.write(str(pos[0]) + " " + str(pos[1]) + " " + str(xComp) + " " + str(yComp) + "\n")

    magFile.close()
    magSliceFile.close()


def main():
    dimensions, accuracy, type = getInputParams()

    if (type.lower() == "jacobi"):
        updateLattice = jacobiUpdate
    else:
        updateLattice = GSUpdate

    chargeLat = np.zeros((dimensions, dimensions, dimensions))
    centre = int(dimensions/2)
    #chargeLat[(centre, centre, centre)] = 10 # place a point charge at the centre of the box
    for i in range(dimensions):
        pos = (centre, centre, i)
        chargeLat[pos] = 10

    potentialLat = generateInitState(dimensions)

    maxDiff = accuracy + 1
    sweep = 0
    while maxDiff > accuracy:

        returnLat = updateLattice(chargeLat, potentialLat, dimensions)
        maxDiff = findMaxDif(potentialLat, returnLat)
        potentialLat = returnLat
        sweep += 1
        print(sweep)
        print("max Diff: " + str(maxDiff))

    #writeElectricField(potentialLat, dimensions)
    #writeEFieldSlice(potentialLat, dimensions, centre)
    writePotential(potentialLat, centre)
    writeMagField(potentialLat, dimensions, centre)

main()
