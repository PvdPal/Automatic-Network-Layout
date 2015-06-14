## Thesis: AUTOMATIC NETWORK LAYOUT USING HILBERT AND MOORE CURVES
## NAME: Patrick van der Pal
## STUDENT ID: 6223931 / 10013555

## Measurements.py
## This file contains all the functions used during the thesis project to
## help express the performance of the several algorithms.

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay
from matplotlib import cm
import numpy as np

from Main import setupDesign
import Measurements

"""
These first three functions are helper functions for the other functions
following these three. The first function calculates the mean value of a
list. The second one determines the median value of a list, while the third
function calculates the area of a triangle given the three corner coordinates.
"""
def mean(list):
    return math.fsum(list) / len(list)

def median(list):
    list = sorted(list)
    if len(list) < 1:
        return None
    if len(list) %2 == 1:
        return list[((len(list)+1)/2)-1]
    else:
        return float(sum(list[(len(list)/2)-1:(len(list)/2)+1]))/2.0

def triangleArea(triangle):
    a, b, c = triangle
    first = a[0] * (b[1] - c[1]) + b[0] * (c[1] - a[1]) + c[0] * (a[1] - b[1])
    return math.fabs(first / 2.0)

"""
This function calculates the midpoint of the given listOfNodes.
"""
def midpoint(listOfNodes):
    sumX = sumY = 0
    for x, y in listOfNodes:
        sumX += x
        sumY += y
    return sumX / len(listOfNodes), sumY / len(listOfNodes)

"""
This function calculates the average distance between every node, given the
listOfNodes and the 'level'. This 'level'-variable is explained in the report.
"""
def averageDistance(listOfNodes, level):
    length = len(listOfNodes)
    j = distance = 0
    for i in range(length - level):
        distance += math.sqrt((listOfNodes[i][0] - listOfNodes[i + level][0])**2 +
                                  (listOfNodes[i][1] - listOfNodes[i + level][1])**2)
        j += 1
    print "Average Distance on Level", level, ":", distance / j
    return distance / j

"""
This function calculates the maximum distance between every node, given the
listOfNodes and 'level'. This 'level'-variable is explained in the report.
"""
def maximumDistance(listOfNodes, level):
    length = len(listOfNodes)
    j = distance = maximumDistance= 0
    for i in range(length - level):
        distance = math.sqrt((listOfNodes[i][0] - listOfNodes[i + level][0])**2 +
                             (listOfNodes[i][1] - listOfNodes[i + level][1])**2)
        if distance > maximumDistance:
            maximumDistance = distance
    print "Maximum Distance on Level", level, ":", maximumDistance

    return maximumDistance

"""
This function creates a histogram of all the distances between the nodes,
given the listOfNodes and the 'level'. The 'level'-variable is explained in
the report.
"""
def distanceHistogram(listOfNodes, level):
    listOfDistances = []
    length = len(listOfNodes)
    for i in range(length - level):
        listOfDistances.append(math.sqrt((listOfNodes[i][0] - listOfNodes[i + level][0])**2 +
                                  (listOfNodes[i][1] - listOfNodes[i + level][1])**2))

    plt.hist(listOfDistances, 50, normed=0, facecolor='g')
    plt.grid(True)
    plt.show()

"""
This function calculates the standard deviation of all the distances between
the nodes, given the listOfNodes and the 'level'. This 'level'-variable is
explained in the report.
"""
def distanceStdDev(listOfNodes, level):
    listOfDistances = []
    length = len(listOfNodes)
    for i in range(length - level):
        listOfDistances.append(math.sqrt((listOfNodes[i][0] - listOfNodes[i + level][0])**2 +
                                  (listOfNodes[i][1] - listOfNodes[i + level][1])**2))
    
    x = 0
    n = len(listOfDistances)
    sum = math.fsum(listOfDistances)
    average = mean(listOfDistances)

    for i in listOfDistances:
        x += math.pow(i - average, 2)

    stddev = math.sqrt(x / n)

    print "The standard deviation of the distance on Level", level, ":", stddev

"""
This function calculates the median value of all the distances between the
nodes, given the listOfNodes and the 'level'. This 'level'-variable is
explained in the report.
"""
def distanceMedian(listOfNodes, level):
    listOfDistances = []
    length = len(listOfNodes)
    for i in range(length - level):
        listOfDistances.append(math.sqrt((listOfNodes[i][0] - listOfNodes[i + level][0])**2 +
                                  (listOfNodes[i][1] - listOfNodes[i + level][1])**2))

    print "The median of the distances on Level", level, ":", median(listOfDistances)
    return median(listOfDistances)

"""
This function calculates the area of the convex hull of the entire design.
This uses the Delaunay function from the 'scipy.spatial'-library.
"""
def areaConvexHull(points):
    dt = Delaunay(points)
    tets = dt.points[dt.simplices]

    polygonArea = 0
    for triangle in tets:
        polygonArea += triangleArea(triangle)
    print "The area of the polygon", polygonArea

"""
This function creates an overview plot of the average distance, the maximum
distance and the median value of all the distances. It uses the appropriate
functions, written above.
"""
def overview(listOfNodes):
    listOfMedian = []
    listOfMaximum = []
    listOfAverage = []
    maxLevel = len(listOfNodes) - 1

    for level in range(1, maxLevel + 1):
        listOfMedian.append([level, distanceMedian(listOfNodes, level)])
        listOfMaximum.append([level, maximumDistance(listOfNodes, level)])
        listOfAverage.append([level, averageDistance(listOfNodes, level)])

    plt.plot(*zip(*listOfMedian), label='Median')
    plt.plot(*zip(*listOfMaximum), label='Maximum')
    plt.plot(*zip(*listOfAverage), label='Average')
    t = range(0, maxLevel + 1, 1)
    plt.plot(t, np.sqrt(t), 'k:', label='Square Root')
    plt.legend()
    plt.xlabel('Level')
    plt.ylabel('Distance')
    plt.xticks(range(0, 42, 2))
    plt.yticks(range(0, 42, 2))
    plt.axis([0, 40, 0, 40])
    plt.show()
    return listOfAverage, listOfMaximum, listOfMedian

"""
'inTheGrandSchemeOfThings' produces a graph where the average distance of
each level is compared with the maximum distance of each level in a
one-dimensional network, i.e. the maximum distance of level 10 in a
one-dimensional network is 10. The same is done with the maximum distance and
the median value of all the distances. This comparison is then displayed as a
percentage.
"""
def inTheGrandSchemeOfThings(nrNodesMin, nrNodesMax, algorithm):
    listOfAverage = [[0.0], []]
    listOfMaximum = [[0.0], []]
    listOfMedian = [[0.0], []]

    for nrNodes in range(nrNodesMin, nrNodesMax, 2):
        listOfNodes = setupDesign(nrNodes, algorithm)
        maxLevel = len(listOfNodes) - 1
        while len(listOfAverage) <= maxLevel and len(listOfMaximum) <= maxLevel and len(listOfMedian) <= maxLevel:
            listOfAverage.append([])
            listOfMaximum.append([])
            listOfMedian.append([])
        for level in range(1, maxLevel + 1):
            average = Measurements.averageDistance(listOfNodes, level)
            listOfAverage[level].append(average / level)
            maximum = Measurements.maximumDistance(listOfNodes, level)
            listOfMaximum[level].append(maximum / level)
            median = Measurements.distanceMedian(listOfNodes, level)
            listOfMedian[level].append(median / level)
    for i in range(1, len(listOfAverage)):
        listOfAverage[i] = math.fsum(listOfAverage[i]) / len(listOfAverage[i])
        listOfMaximum[i] = math.fsum(listOfMaximum[i]) / len(listOfMaximum[i])
        listOfMedian[i] = math.fsum(listOfMedian[i]) / len(listOfMedian[i])
    del listOfAverage[0]
    del listOfMaximum[0]
    del listOfMedian[0]
    plt.plot(range(1, 1 + len(listOfAverage)), listOfAverage, label = 'Average')
    plt.plot(range(1, 1 + len(listOfMaximum)), listOfMaximum, label = 'Maximum')
    plt.plot(range(1, 1 + len(listOfMedian)), listOfMedian, label = 'Median')
    plt.legend(loc=3)
    plt.axis([0, maxLevel, 0, 1.1])
    plt.show()

"""
'inTheGrandestSchemeOfThings' does the same as 'inTheGrandSchemeOfThings',
with the differences that the plot of the maximum distance and the median
value of all the distances are removed. In their place, the percentages,
based on the average distances are compared between the four algorithms in
the chosen category. Below is stated which category holds which algorithms.
"""
def inTheGrandestSchemeOfThings(nrNodesMin, nrNodesMax, category):
    listOfAlgorithm1 = [[0.0], []]
    listOfAlgorithm2 = [[0.0], []]
    listOfAlgorithm3 = [[0.0], []]
    listOfAlgorithm4 = [[0.0], []]

    for nrNodes in range(nrNodesMin, nrNodesMax, 2):
        if category == 1:
            """A design based on Hilbert curves with an even amount of nodes."""
            listOfNodes1 = setupDesign(nrNodes, 1)
            listOfNodes2 = setupDesign(nrNodes, 2)
            listOfNodes3 = setupDesign(nrNodes, 3)
            listOfNodes4 = setupDesign(nrNodes, 4)
        elif category == 2:
            """A design based on Hilbert curves with an uneven amount of nodes."""
            listOfNodes1 = setupDesign(nrNodes, 5)
            listOfNodes2 = setupDesign(nrNodes, 6)
            listOfNodes3 = setupDesign(nrNodes, 7)
            listOfNodes4 = setupDesign(nrNodes, 8)
        elif category == 3:
            """A design based on Z-order curves with an even amount of nodes."""
            listOfNodes1 = setupDesign(nrNodes, 9)
            listOfNodes2 = setupDesign(nrNodes, 10)
            listOfNodes3 = setupDesign(nrNodes, 11)
            listOfNodes4 = setupDesign(nrNodes, 12)
        elif category == 4:
            """A design based on Z-order curves with an even amount of nodes."""
            listOfNodes1 = setupDesign(nrNodes, 13)
            listOfNodes2 = setupDesign(nrNodes, 14)
            listOfNodes3 = setupDesign(nrNodes, 15)
            listOfNodes4 = setupDesign(nrNodes, 16)

        maxLevel = len(listOfNodes2) - 1
        while len(listOfAlgorithm1) <= maxLevel and len(listOfAlgorithm2) <= maxLevel and len(listOfAlgorithm3) <= maxLevel and len(listOfAlgorithm4) <= maxLevel:
            listOfAlgorithm1.append([])
            listOfAlgorithm2.append([])
            listOfAlgorithm3.append([])
            listOfAlgorithm4.append([])
        for level in range(1, maxLevel + 1):
            average = Measurements.averageDistance(listOfNodes1, level)
            listOfAlgorithm1[level].append(average / level)
            average = Measurements.averageDistance(listOfNodes2, level)
            listOfAlgorithm2[level].append(average / level)
            average = Measurements.averageDistance(listOfNodes3, level)
            listOfAlgorithm3[level].append(average / level)
            average = Measurements.averageDistance(listOfNodes4, level)
            listOfAlgorithm4[level].append(average / level)
    for i in range(1, 1 + maxLevel):
        listOfAlgorithm1[i] = math.fsum(listOfAlgorithm1[i]) / len(listOfAlgorithm1[i])
        listOfAlgorithm2[i] = math.fsum(listOfAlgorithm2[i]) / len(listOfAlgorithm2[i])
        listOfAlgorithm3[i] = math.fsum(listOfAlgorithm3[i]) / len(listOfAlgorithm3[i])
        listOfAlgorithm4[i] = math.fsum(listOfAlgorithm4[i]) / len(listOfAlgorithm4[i])
    del listOfAlgorithm1[0]
    del listOfAlgorithm2[0]
    del listOfAlgorithm3[0]
    del listOfAlgorithm4[0]
    plt.plot(range(1, 1 + len(listOfAlgorithm1)), listOfAlgorithm1, label = 'Algorithm1')
    plt.plot(range(1, 1 + len(listOfAlgorithm2)), listOfAlgorithm2, label = 'Algorithm2')
    plt.plot(range(1, 1 + len(listOfAlgorithm3)), listOfAlgorithm3, label = 'Algorithm3')
    plt.plot(range(1, 1 + len(listOfAlgorithm4)), listOfAlgorithm4, label = 'Algorithm4')
    plt.legend(loc=3)
    plt.axis([0, maxLevel, 0, 1.1])
    plt.show()
