#!/usr/bin/python

import py4j
import pandas as pd
import data2neoFunctions

###############################################################################
#############################AUXILIAR FUNCTIONS################################
###############################################################################

def getHeaderLen(gafFile, headerChar):
    gaf = open(gafFile)
    headerLen = 0
    line = gaf.readline()
    while line.startswith(headerChar):
        headerLen += 1
        line = gaf.readline()
    gaf.close()
    return headerLen

def insertGAFRelation(gafRow, nameFieldNumber, edges, graph):
    gene = graph.nodeMatch("Gene", name=gafRow[nameFieldNumber]).first()
    goTerm = None
    if ":" in gafRow[4]:
        goTerm = graph.nodeMatch("GO_Term", accession=gafRow[4].split(":")[1]).first()
    if gene and goTerm:
        relType = gafRow[3].upper() if str(gafRow[3]) != "nan" else "ANNOTATED_AS" 
        relType = relType.replace("NOT|","NOT_")
        geneToTerm = py4j.Relationship(gene, goTerm, relType, evidenceCode = gafRow[6], assignedBy = gafRow[14])
        edges.add(geneToTerm)
    return gene!=None and goTerm!=None

###############################################################################
###############################################################################
###############################################################################
# GAF 2.1/2.2/2.0
# 1  ->  Gene name
# 3  ->  Qualifier (connection type)
# 4  ->  GO ID
# 6  ->  Evidence code (annotation type)
# 14 ->  Source of annotation

def gaf2neo(graphData, gafFile, nameFieldNumber):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])
    
    headerLen = getHeaderLen(gafFile, "!")
    insertedRels, totalRels = 0, 0
    edges = set()

    gaf = pd.read_table(gafFile, header=None, delimiter='\t', on_bad_lines='skip', usecols=[nameFieldNumber,3,4,6,14], skiprows=headerLen)
    for rowNumber in gaf.index:
        totalRels += 1
        insertedRels += int(insertGAFRelation(gaf.iloc[rowNumber], nameFieldNumber, edges, graph))

    #VER QUE PASA CON EL REL TYPE QUE TIENE UN PIPE CON "NOT"
    graph.merge(py4j.SubGraph(set(), edges))
    print("Total relationships in file:", totalRels)
    print("Inserted relationships:", insertedRels)
    
if __name__ == '__main__':
    gaf2neo()
