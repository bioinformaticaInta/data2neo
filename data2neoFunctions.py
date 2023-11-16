#!/usr/bin/python

import py4j
from filesplit.split import Split
from math import ceil
import os
import sh

def splitWithHeaderProcs(file, headerChar, nprocs, manFileName = "manifest"):
    Nheader, Nlines = 0, 0
    chroms = set()
    with open(file) as allFile:
        for vcfLine in allFile:
            if vcfLine.startswith(headerChar):
                Nheader += 1
            else:
                Nlines += 1
                chroms.add(vcfLine.split()[0])
    linesForFile = ceil((Nlines+nprocs*Nheader)/(nprocs))
    allFile = Split(file,".")
    allFile.manfilename = manFileName
    allFile.bylinecount(linesForFile, Nheader)
    return chroms

def splitWithHeaderLines(file, headerChar, chunkLines, manFileName = "manifest"):
    Nheader = 0
    allFile = open(file)
    while allFile.readline().startswith(headerChar):
        Nheader += 1
    allFile.close()
    allFile = Split(file,".")
    allFile.manfilename = manFileName
    allFile.bylinecount(chunkLines, Nheader)

def removeSplitFiles(manFileName = "manifest"):
    with open(manFileName) as mfFile:
        next(mfFile)
        for mfLine in mfFile:
            mfLine = mfLine[:-1].split(",")
            os.remove(mfLine[0])
    os.remove(manFileName)

def createVertexIndexes(graph, dataSet):
    for labelName,propName in dataSet:
        graph.createNodeIndex(labelName, propName)

def loadVCFDataToGraph(graphData, listOfOuts):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])

    edges, dataForIndexes = set(), set()
    print("Grouping edges...")
    for output in listOfOuts:
        edges |= output[0]
        dataForIndexes |= output[1]

    createVertexIndexes(graph, dataForIndexes)

    print("Loading edges...")
    graph.merge(py4j.SubGraph(set(),edges))
    print("Loaded edges!!!")

def csvRelsHeaderFormating(csvFile):
    csv = open(csvFile)
    firstLineList = csv.readline()[:-1].split(",")
    csv.close()
    firstLineNew = "start:"+firstLineList[0]
    pos = 1
    for col in firstLineList[1:]:
        begin = "rel:"
        if pos == len(firstLineList)-1:
            begin = "end:"
        firstLineNew += ","+begin+col
        pos += 1
    sh.sed("-i", "1s/.*/" + firstLineNew + "/", csvFile)

def csvRelsHeaderReFormating(csvFile):
    csv = open(csvFile)
    firstLineList = csv.readline()[:-1].split(",")
    csv.close()
    pos = 0
    firstLineNew = ":".join(firstLineList[0].split(":")[1:])
    for col in firstLineList[1:]:
        firstLineNew += ","+":".join(col.split(":")[1:])
    sh.sed("-i", "1s/.*/" + firstLineNew + "/", csvFile)

def firstVertexInSet(vertexs, labels=None, props=None):
    labels = labels if labels else ()
    props = props if props else {}
    vertex = None
    for actual in vertexs:
        if actual.isNode(*labels, **props):
            vertex = actual
            break
    return vertex

def firstEdgeInSet(edges, startLabels=None, startProps=None, relType=None, relProps=None, endLabels=None, endProps=None):
    startLabels = startLabels if startLabels else ()
    endLabels = endLabels if endLabels else ()
    startProps = startProps if startProps else {}
    endProps = endProps if endProps else {}
    relProps = relProps if relProps else {}
    edge = None
    for actual in edges:
        if actual.startNode().isNode(*startLabels, **startProps) and actual.endNode().isNode(*endLabels, **endProps) and actual.isEdge(relType, **relProps):
            edge = actual
            break
    return edge

def allVertexsInSet(vertexs, labels=None, props=None):
    labels = labels if labels else ()
    props = props if props else {}
    vertexs = set()
    for actual in vertexs:
        if actual.isNode(*labels, **props):
            vertexs.add(actual)
    return vertex

def allEdgesInSet(edges, startLabels=None, startProps=None, relType=None, relProps=None, endLabels=None, endProps=None):
    startLabels = startLabels if startLabels else ()
    endLabels = endLabels if endLabels else ()
    startProps = startProps if startProps else {}
    endProps = endProps if endProps else {}
    relProps = relProps if relProps else {}
    edges = set()
    for actual in edges:
        if actual.startNode().isNode(*startLabels, **startProps) and actual.endNode().isNode(*endLabels, **endProps) and actual.isEdge(relType, **relProps):
            edges.add(actual)
    return edge
