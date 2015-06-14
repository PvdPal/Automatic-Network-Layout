## Thesis: AUTOMATIC NETWORK LAYOUT USING HILBERT AND MOORE CURVES
## NAME: Patrick van der Pal
## STUDENT ID: 6223931 / 10013555

## HilbertAlgorithms.py
## This file contains all the functions used during the thesis project to
## build the design algorithms, based on building blocks made from Hilbert
## curves.

import math
import matplotlib.pyplot as plt
from random import randint

import Curves
import Measurements

"""
The Uneven Building Blocks
This function sets up the right amount of building blocks, given the number
of nodes and the largest iteration. This function allows the number of
building blocks to be uneven. If this is the case the middle building block
is a Moore curve on top of the two stacks of Hilbert curve.
"""
def setupBuildingBlocksUneven(nrNodes, iteration):
    listOfNodes = []
    nrBuildingBlocks = math.ceil(nrNodes / math.pow(4, iteration))
    for i in range(int(math.floor(nrBuildingBlocks / 2))):
        listOfNodes.extend(Curves.hilbertCurve(iteration, 'leftReverse', 0, i * math.pow(2, iteration)))
    if(nrBuildingBlocks % 2 != 0):
        listOfNodes.extend(Curves.mooreCurve(iteration, 'top', 0.5 * math.pow(2, iteration), math.floor(nrBuildingBlocks / 2) * math.pow(2, iteration)))
    for i in range(int(math.floor(nrBuildingBlocks / 2))):
        listOfNodes.extend(Curves.hilbertCurve(iteration, 'rightReverse', math.pow(2, iteration), (math.floor(nrBuildingBlocks / 2) - i - 1) * math.pow(2, iteration)))
    return listOfNodes

"""
The Even Building Blocks
This function sets up the right amount of building blocks, given the number
of nodes and the largest iteration. This function does not allow the number
of building blocks to be uneven. Therefore, this returns two stacks of
Hilbert curves.
"""
def setupBuildingBlocksEven(nrNodes, iteration):
    listOfNodes = []
    nrBuildingBlocks = math.ceil(nrNodes / math.pow(4, iteration))
    for i in range(int(math.ceil(nrBuildingBlocks / 2))):
        listOfNodes.extend(Curves.hilbertCurve(iteration, 'leftReverse', 0, i * math.pow(2, iteration)))
    for i in range(int(math.ceil(nrBuildingBlocks / 2))):
        listOfNodes.extend(Curves.hilbertCurve(iteration, 'rightReverse', math.pow(2, iteration), (math.ceil(nrBuildingBlocks / 2) - i - 1) * math.pow(2, iteration)))
    return listOfNodes

""" This function cuts off nodes from the listOfNodes randomly. """
def cuttingOffNodes(listOfNodes):
    del listOfNodes[randint(0, len(listOfNodes) - 1)]

"""
This function cuts off the single nodes from the given listOfNodes that are
farthest from the calculated midpoint.
"""
def cuttingOffNodesOptimal(listOfNodes):
    listOfNominees = []
    averageX, averageY = Measurements.midpoint(listOfNodes)
    for i in range(len(listOfNodes)):
        newX = listOfNodes[i][0]
        newY = listOfNodes[i][1]
        distance = math.sqrt(math.pow(newX - averageX, 2) + math.pow(newY - averageY, 2))
        listOfNominees.append((i, distance))
    index = zip(*listOfNominees)[1].index(max(*zip(*listOfNominees)[1]))
    del listOfNodes[listOfNominees[index][0]]

"""
This function cuts off pairs of nodes from the given listOfNodes. This uses
the 'Cutting Mechanism' described in the report.
"""
def cuttingOffPairs(listOfNodes):
    for i in range(len(listOfNodes) - 3):
        distance = math.sqrt((listOfNodes[i][0] - listOfNodes[i + 3][0])**2 + (listOfNodes[i][1] - listOfNodes[i + 3][1])**2)
        if(distance == 1.0):
            del listOfNodes[i + 1]
            del listOfNodes[i + 1]
            break

"""
This function cuts off pairs of nodes from the given listOfNodes. The pair
that is farthest from the calculated midpoint will be removed. This uses the
'Cutting Mechanism' described in the report.
"""
def cuttingOffPairsOptimal(listOfNodes):
    listOfNominees = []
    averageX, averageY = Measurements.midpoint(listOfNodes)
    for i in range(len(listOfNodes) - 3):
        distance = math.sqrt((listOfNodes[i][0] - listOfNodes[i + 3][0])**2 + (listOfNodes[i][1] - listOfNodes[i + 3][1])**2)
        if(distance == 1.0):
            newX = (listOfNodes[i + 1][0] + listOfNodes[i + 2][0]) / 2
            newY = (listOfNodes[i + 1][1] + listOfNodes[i + 2][1]) / 2
            distance = math.sqrt(math.pow(newX - averageX, 2) + math.pow(newY - averageY, 2))
            listOfNominees.append((i + 1, i + 2, distance))
    index = zip(*listOfNominees)[2].index(max(*zip(*listOfNominees)[2]))
    del listOfNodes[listOfNominees[index][0]]
    del listOfNodes[listOfNominees[index][0]]

def addingPairs(listOfNodes):
    pass

"""
These functions are the several combinations of the functions above,
corresponding with the numbering of the algorithms in the thesis.
"""
def algorithm1(nodes, iteration):
    listOfNodes = setupBuildingBlocksUneven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffPairs(listOfNodes)
    return listOfNodes

def algorithm2(nodes, iteration):
    listOfNodes = setupBuildingBlocksUneven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffPairsOptimal(listOfNodes)
    return listOfNodes

def algorithm3(nodes, iteration):
    listOfNodes = setupBuildingBlocksEven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffPairs(listOfNodes)
    return listOfNodes

def algorithm4(nodes, iteration):
    listOfNodes = setupBuildingBlocksEven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffPairsOptimal(listOfNodes)
    return listOfNodes

def algorithm5(nodes, iteration):
    listOfNodes = setupBuildingBlocksUneven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffNodes(listOfNodes)
    return listOfNodes

def algorithm6(nodes, iteration):
    listOfNodes = setupBuildingBlocksUneven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffNodesOptimal(listOfNodes)
    return listOfNodes

def algorithm7(nodes, iteration):
    listOfNodes = setupBuildingBlocksEven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffNodes(listOfNodes)
    return listOfNodes

def algorithm8(nodes, iteration):
    listOfNodes = setupBuildingBlocksEven(nodes, iteration)
    while(nodes < len(listOfNodes)):
        cuttingOffNodesOptimal(listOfNodes)
    return listOfNodes
