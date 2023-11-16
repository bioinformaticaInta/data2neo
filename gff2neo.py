#!/usr/bin/python

from BCBio import GFF
import py4j
import warnings
import data2neoFunctions

#######################################################################
###########################AUXILIAR FUNCTIONS##########################
#######################################################################

def getStrandChar(number):
    strand = "+"
    if number == -1:
        strand = "-"
    return strand

def vertexFromFeature(feature, level = 1):
    #Extraer y capitalizar el tipo para usarlo como label del vertice
    labels = [feature.type.title().replace("_","")]
    name = feature.id
    if name:
        if ":" in name:
            name = name.split(":")[1]
    else:
        name = feature.qualifiers["Name"][0]
    if level == 2: labels.append("Transcript")
    return py4j.Node(*labels, name = name, source = feature.qualifiers["source"][0])

def addVertexIndexToSet(vertex, indexProp, dataSet):
    for label in vertex.labels():
        dataSet.add((label, indexProp))

def insertProtein(edges, mrna):
    prot = py4j.Node("Protein", name = mrna["name"])
    protToVertex = py4j.Relationship(prot , mrna, "TRANSLATED_FROM")
    edges.add(protToVertex)
    return prot

def getChromChain(gffRec):
    resolution = 200
    chromChainVertexs = set()
    chromChainEdges = set()
    namePrefix = str(gffRec.id)+":"
    chromLen = len(gffRec.seq)
    
    pre = None
    start = 0
    while start <= chromLen:
        end = start+resolution-1 if start+resolution-1<=chromLen else chromLen
        chainName = namePrefix+str(start)+"-"+str(end)
        actual = py4j.Node("ChromosomeChain", name=chainName, start=start, end=end)
        chromChainVertexs.add(actual)
        print(actual)
        if pre:
            chromChainEdges.add(py4j.Relationship(pre, actual, "NEXT"))
        pre = actual
        start += resolution
    return chromChainVertexs, chromChainEdges        

#######################################################################
#######################################################################
#######################################################################

#genes ->  level 1
#RNAs  ->  level 2
#exons ->  level 3
def insertFeature(feature, edges, chrom, vertexIndexesData, levelMax, levelAct = 1, parent = None, exonSet = set()):
    vertex = vertexFromFeature(feature, levelAct)
    vertexInGraph = None
    #Si es un exon, controlamos que no este como vertice en el grafo
    if levelAct == 3:
        vertexInGraph = tuple(filter(lambda x: (x.startNode().hasLabel(vertex.labels()[0]) and x["start"]==int(feature.location.start) and x["end"]==int(feature.location.end) and x["strand"]==getStrandChar(feature.location.strand) and x.endNode()==chrom), exonSet))
        #Si el exon ya existe, se usa el vertice existente
        if vertexInGraph: vertex = vertexInGraph[0].startNode()
    #Se inserta el vertice si no existe en el grafo
    if not vertexInGraph:
        addVertexIndexToSet(vertex, "name", vertexIndexesData)
        #Se inserta una arista con la localizacion en el cromosoma si es gen o exon
        if levelAct in (1,3):
            #Se crea la arista
            vertexToChrom = py4j.Relationship(vertex, chrom, "LOCATED_AT",  start = int(feature.location.start), end = int(feature.location.end), strand = getStrandChar(feature.location.strand))
            edges.add(vertexToChrom)
            #Se almacena la arista-exonToChrom en la lista de exones del gen
            if levelAct == 1:
                exonSet = set()
            else:
                exonSet.add(vertexToChrom)
        #Si es un mensajero no se localiza pero se crea una proteina
        if str(feature.type) == "mRNA":
            protVertex = insertProtein(edges, vertex)
            addVertexIndexToSet(protVertex, "name", vertexIndexesData)
        #Se crea la arista de tipo PART_OF si tiene un feature parent
        if parent:
            relLabel = "PART_OF"
            if levelAct == 2: relLabel = "TRANSCRIBED_FROM" #RNAs
            vertexToParent = py4j.Relationship(vertex, parent, relLabel)
            edges.add(vertexToParent)
        #Se sigue con la insercion de los subfeatures, por ejemplo: mRNAs de genes o exones de mRNAs
        if levelAct < levelMax:
            if feature.sub_features:
                for subFeature in feature.sub_features:
                    if (levelAct == 1) or (levelAct == 2 and str(subFeature.type) == "exon"): #En nivel 3 solo se insertan exones
                        insertFeature(subFeature, edges, chrom, vertexIndexesData, levelMax, levelAct+1, vertex, exonSet)

#######################################################################
#######################################################################
#######################################################################

def gff2neo(graphData, gffFile, structured, loadExons, chromFeature):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])

    if not structured and loadExons:
        warnings.warn("The option -e is ignored for unstructured gff files")

    maxLevel = 1
    if structured:
        maxLevel = 3 if loadExons else 2
    
    edges, vertexIndexesData = set(), set()
    gff = open(gffFile)
    for rec in GFF.parse(gff):
        print(rec)
        print(len(rec.seq))
        chromChain = graph.nodeMatch("ChromosomeChain", chromName=rec.id).all()
        if not chromChain: chromChain = getChromChain(rec)
        #recFeatures = rec.features[1:] if chromFeature else rec.features
        #print("Insert features from:", rec.id)
        #for feature in recFeatures:
        #    insertFeature(feature, edges, chrom, vertexIndexesData, maxLevel)
    gff.close()

    #data2neoFunctions.createVertexIndexes(graph, vertexIndexesData)
    #graph.merge(py4j.SubGraph(set(), edges))

if __name__ == '__main__':
    gff2neo()
