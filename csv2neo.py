#!/usr/bin/python

import py4j
import pandas as pd
import data2neoFunctions
import copy as cp
import datetime

###############################################################################
##############################AUXILIAR FUNCTIONS###############################
###############################################################################

def findOrCreateVertex(graph, label, propName, propValue, createNodes):
    vertex = graph.nodeMatch(label, **{propName:propValue}).first()
    if not vertex and createNodes:
        vertex = py4j.Node(label, **{propName:propValue})
    return vertex

def propsFromRow(rowData, colNames, propF):
    propDict = {}
    namePos = 0 if len(colNames[0].split(":")) == 2 else 1
    for col in colNames:
        propData = col.split(":")
        propValue = rowData[col]
        if str(propValue) != "nan":
            propDict[propData[namePos]] = propF[col](propValue)
    return propDict
        

###############################################################################
###############################################################################
###############################################################################

#header: prop1Name:type , prop2Name:type , ...... , propLName:type
#lines:  prop1Value , prop2Value , ...... , propLValue
def csvNodes2neo(graphData, labelsList):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd']) 
    csv = pd.read_csv(graphData['file'], on_bad_lines='skip')

    dicF = {"date":datetime.date.fromisoformat , "bool":bool , "int":int , "str":str , "float":float}
    colNames = tuple(csv.columns.values)
    propNames = tuple(col.split(":")[1] for col in colNames)
    propF = {col : dicF[col.split(":")[1]] for col in colNames}
    labelsList = tuple(x.title() for x in labelsList)

    nodes = set()
    for rowNumber in csv.index:
        props = propsFromRow(csv.iloc[rowNumber], colNames, propF)
        vertex = py4j.Node(*labelsList, **props)
        nodes.add(vertex)
   
    graph.createNodeIndex(labelsList[0], propNames[0])
    graph.merge(py4j.SubGraph(nodes, set()))

#header: start:propStartName:type , rel:proprel1Name:type, rel:propRel2Name:type , ...... , rel:propRelLName:type , end:propEndName:type
#lines:  propStartValue , proprel1Value, propRel2Value , ...... , propRelLValue , propEndValue
def csvRels2neo(graphData, startLabel, endLabel, relType, createNodes):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])
    csv = pd.read_csv(graphData['file'], on_bad_lines='skip')
    
    dicF = {"date":datetime.date.fromisoformat , "bool":bool , "int":int , "str":str , "float":float}
    colNames = tuple(csv.columns.values)
    propNames = tuple(col.split(":")[1] for col in colNames)
    propF = {col : dicF[col.split(":")[2]] for col in colNames}

    relType = relType.upper()
    startLabel = startLabel.title()
    endLabel = endLabel.title()

    edges = set()
    for rowNumber in csv.index:
        rowData = csv.iloc[rowNumber]
        startVertex = findOrCreateVertex(graph, startLabel, propNames[0], propF[colNames[0]](rowData[colNames[0]]), createNodes)
        endVertex = findOrCreateVertex(graph, endLabel, propNames[-1], propF[colNames[-1]](rowData[colNames[-1]]), createNodes)
        if startVertex and endVertex:
            props = propsFromRow(rowData, colNames[1:-1], propF)
            edge = py4j.Relationship(startVertex, endVertex, relType, **props)
            edges.add(edge)

    graph.createNodeIndex(startLabel, propNames[0])
    graph.createNodeIndex(endLabel, propNames[-1])
    graph.merge(py4j.SubGraph(set(), edges))
    
if __name__ == '__main__':
    csvNodes2neo()
