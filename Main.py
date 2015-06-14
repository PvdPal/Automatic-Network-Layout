## Thesis: AUTOMATIC NETWORK LAYOUT USING HILBERT AND MOORE CURVES
## NAME: Patrick van der Pal
## STUDENT ID: 6223931 / 10013555

## Main.py
## This file contains the 'networkLayout()'-function that, given a number of
## nodes and a requested variant of the design algorithm, runs the algorithm
## and collects the statistics. The only way to control which statistics are
## collected and which are not, is to use '#' to comment out the appropriate
## lines.

## This is the numbering of the different algorithms. These correspond with
## the numbering in the thesis.
##    01 = rounding up to an uneven amount of building blocks, based on
##         Hilbert curves, and cutting off pairs of nodes (default)
##    02 = rounding up to an uneven amount of building blocks, based on
##         Hilbert curves, and cutting off pairs of nodes (optimized)
##    03 = rounding up to an even amount of building blocks, based on Hilbert
##         curves, and cutting off pairs of nodes (default)
##    04 = rounding up to an even amount of building blocks, based on Hilbert
##         curves, and cutting off pairs of nodes (optimized)
##    05 = rounding up to an uneven amount of building blocks, based on
##         Hilbert curves, and cutting off single nodes (random)
##    06 = rounding up to an uneven amount of building blocks, based on
##         Hilbert curves, and cutting off single nodes (optimized)
##    07 = rounding up to an even amount of building blocks, based on Hilbert
##         curves, and cutting off single nodes (random)
##    08 = rounding up to an even amount of building blocks, based on Hilbert
##         curves, and cutting off single nodes (optimized)
##    09 = rounding up to an uneven amount of building blocks, based on
##         Z-order curves, and cutting off pairs of nodes (default)
##    10 = rounding up to an uneven amount of building blocks, based on
##         Z-order curves, and cutting off pairs of nodes (optimized)
##    11 = rounding up to an even amount of building blocks, based on Z-order
##         curves, and cutting off pairs of nodes (default)
##    12 = rounding up to an even amount of building blocks, based on Z-order
##         curves, and cutting off pairs of nodes (optimized)
##    13 = rounding up to an uneven amount of building blocks, based on
##         Z-order curves, and cutting off single nodes (random)
##    14 = rounding up to an uneven amount of building blocks, based on
##         Z-order curves, and cutting off single nodes (optimized)
##    15 = rounding up to an even amount of building blocks, based on Z-order
##         curves, and cutting off single nodes (random)
##    16 = rounding up to an even amount of building blocks, based on Z-order
##         curves, and cutting off single nodes (optimized)

import numpy as np
from scipy.spatial import ConvexHull

import math
import matplotlib.pyplot as plt
import time

import HilbertAlgorithms
import Z_OrderAlgorithms
import Curves
import Measurements

def convexHull(listOfNodes):
    iteration = math.floor(math.log(len(listOfNodes), 4))
    fig = plt.figure()
    points = np.array(listOfNodes)
    hull = ConvexHull(points)
    for simplex in hull.simplices:
        plt.plot(points[simplex, 0], points[simplex, 1], 'r--', figure=fig)
    plt.axis([-1, 2 * math.pow(2, iteration), -1, math.ceil(math.ceil(len(listOfNodes) / math.pow(4, iteration)) / 2) * math.pow(2, iteration)])
    plt.show()
    plt.close(fig)
    Measurements.areaConvexHull(points)

def setupDesign(nrNodes, algorithm):
    iteration = math.floor(math.log(nrNodes, 4))
    if algorithm == 1:
        return HilbertAlgorithms.algorithm1(nrNodes, iteration)
    elif algorithm == 2:
        return HilbertAlgorithms.algorithm2(nrNodes, iteration)
    elif algorithm == 3:
        return HilbertAlgorithms.algorithm3(nrNodes, iteration)
    elif algorithm == 4:
        return HilbertAlgorithms.algorithm4(nrNodes, iteration)
    elif algorithm == 5:
        return HilbertAlgorithms.algorithm5(nrNodes, iteration)
    elif algorithm == 6:
        return HilbertAlgorithms.algorithm6(nrNodes, iteration)
    elif algorithm == 7:
        return HilbertAlgorithms.algorithm7(nrNodes, iteration)
    elif algorithm == 8:
        return HilbertAlgorithms.algorithm8(nrNodes, iteration)

    elif algorithm == 9:
        return Z_OrderAlgorithms.algorithm1(nrNodes, iteration)
    elif algorithm == 10:
        return Z_OrderAlgorithms.algorithm2(nrNodes, iteration)
    elif algorithm == 11:
        return Z_OrderAlgorithms.algorithm3(nrNodes, iteration)
    elif algorithm == 12:
        return Z_OrderAlgorithms.algorithm4(nrNodes, iteration)
    elif algorithm == 13:
        return Z_OrderAlgorithms.algorithm5(nrNodes, iteration)
    elif algorithm == 14:
        return Z_OrderAlgorithms.algorithm6(nrNodes, iteration)
    elif algorithm == 15:
        return Z_OrderAlgorithms.algorithm7(nrNodes, iteration)
    elif algorithm == 16:
        return Z_OrderAlgorithms.algorithm8(nrNodes, iteration)

def createDesign(nrNodes, algorithm):
    listOfNodes = setupDesign(nrNodes, algorithm)
    nrIterations = math.floor(math.log(len(listOfNodes), 4))
    
    fig = plt.figure()
    listOfNodes = [(x, y) for x, y in listOfNodes]
    plt.plot(*zip(*listOfNodes), figure=fig)
    plt.axis([-1, 2 * math.pow(2, nrIterations), -1, math.ceil(math.ceil(nrNodes / math.pow(4, nrIterations)) / 2) * math.pow(2, nrIterations)])
    plt.show()
    plt.close(fig)
    return listOfNodes

def createDesigns(nrNodesMin, nrNodesMax, algorithm):
    for nrNodes in range(nrNodesMin, nrNodesMax, 2):
        createDesign(nrNodes, algorithm)

def analyzeDesign(nrNodes, algorithm):
    listOfNodes = setupDesign(nrNodes, algorithm)

    convexHull(listOfNodes)
    Measurements.overview(listOfNodes)

def analyzeDesigns(nrNodesMin, nrNodesMax, algorithm):
    for nrNodes in range(nrNodesMin, nrNodesMax, 2):
        analyzeDesign(nrNodes, algorithm)

def analysisOverview(nrNodesMin, nrNodesMax, algorithm):
    Measurements.inTheGrandSchemeOfThings(nrNodesMin, nrNodesMax, algorithm)

def analyzeAlgorithm(nrNodesMin, nrNodesMax, algorithm, level):
    listOfAverage = []
    listOfMaximum = []
    listOfMedian = []
    for nrNodes in range(nrNodesMin, nrNodesMax, 2):
        listOfNodes = setupDesign(nrNodes, algorithm)
        average = Measurements.averageDistance(listOfNodes, level)
        listOfAverage.append([nrNodes, average])
        maximum = Measurements.maximumDistance(listOfNodes, level)
        listOfMaximum.append([nrNodes, maximum])
        median = Measurements.distanceMedian(listOfNodes, level)
        listOfMedian.append([nrNodes, median])
    plt.plot(*zip(*listOfAverage), label = 'Average')
    plt.plot(*zip(*listOfMaximum), label = 'Maximum')
    plt.plot(*zip(*listOfMedian), label = 'Median')
    plt.axis([nrNodesMin, nrNodesMax, 0, 10])
    plt.legend()
    plt.show()
        
def algorithmOverview(nrNodesMin, nrNodesMax, category):
    Measurements.inTheGrandestSchemeOfThings(nrNodesMin, nrNodesMax, category)

def designOverview(nrNodes, algorithm):
    createDesign(nrNodes, algorithm)
    analyzeDesign(nrNodes, algorithm)

def designsOverview(nrNodesMin, nrNodesMax, algorithm):
    for nrNodes in range(nrNodesMin, nrNodesMax, 2):
        createDesign(nrNodes, algorithm)
        analyzeDesign(nrNodes, algorithm)

if __name__ == '__main__':
    import os, sys, argparse
    #createDesign(40, 8)
    #analyzeDesign(40, 8)
    algorithmOverview(6, 66, 1)
