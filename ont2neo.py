#!/usr/bin/python

import pronto
import py4j
import data2neoFunctions
import time

###############################################################################
###############################AUXILIAR FUNCTIONS##############################
###############################################################################

def insertSuperTerms(ontology, term, preTerms, preRels, termVertex, insertedTerms, edges):
    for superTerm in term.superclasses(1,False):
        superVertex = insertTerm(ontology, superTerm, preTerms, preRels, insertedTerms, edges)
        if superVertex and not data2neoFunctions.firstEdgeInSet(preRels, startProps={"accession":termVertex["accession"]}, relType="instanceOf", endProps={"accession":superVertex["accession"]}):
            termToSuper = py4j.Relationship(termVertex , superVertex, "IS_A")
            edges.add(termToSuper)

def insertRelatedTerms(ontology, term, preTerms, preRels, termVertex, insertedTerms, edges):
    analogTypes = {"PART_OF":"componentOf" , "HAS_PART":"hasPart"}
    for relationType in term.relationships:
        for relatedTerm in term.relationships[ontology.get_relationship(relationType.id)]:
            ##Terms with a specific relation type
            relatedVertex = insertTerm(ontology, relatedTerm, preTerms, preRels, insertedTerms, edges)
            relTypeID = str(relationType.id).upper()
            relToSearch = analogTypes[relTypeID] if relTypeID in ("HAS_PART","PART_OF") else relTypeID
            if relatedVertex and not data2neoFunctions.firstEdgeInSet(preRels, startProps={"accession":termVertex["accession"]}, relType=relToSearch, endProps={"accession":relatedVertex["accession"]}):
                termToRelated = py4j.Relationship(termVertex, relatedVertex, relTypeID)
                edges.add(termToRelated)

def insertTerm(ontology, term, preTerms, preRels, insertedTerms, edges):
    ontName = ontology.metadata.ontology.upper()
    termAcc = str(term.id).split(":")[1]
    termVertex = data2neoFunctions.firstVertexInSet(insertedTerms, props={"accession":termAcc})
    if not termVertex and not term.obsolete:
        ##Term insert
        termVertex = data2neoFunctions.firstVertexInSet(preTerms, props={"accession":termAcc})
        if not termVertex:
            labels = (ontName+"_"+str(term.namespace).title().replace("_",""), ontName+"_Term")
            termVertex = py4j.Node(*labels, accession = termAcc, databaseName = ontName, definition = str(term.definition), displayName = str(term.name), name = str(term.name))
        insertedTerms.add(termVertex)
        ##superClasses insert with "is_a" connections (fathers)
        insertSuperTerms(ontology, term, preTerms, preRels, termVertex, insertedTerms, edges)
        ##Relationships insert with "part_of", "has_part" or other type of connections
        insertRelatedTerms(ontology, term, preTerms, preRels, termVertex, insertedTerms, edges)
    return termVertex

###############################################################################
###############################################################################
###############################################################################

def ont2neo(graphData, ontFile):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])
    ontology = pronto.Ontology(ontFile)
    ontName = ontology.metadata.ontology.upper()
    insertedTerms, edges = set(), set()
    preTerms = set(graph.nodeMatch(ontName+"_Term").all())
    preRels = set(py4j.Match(graph.run("MATCH (s:%s)-[l]-(e:%s) RETURN s,l,e" % (ontName+"_Term",ontName+"_Term"), "l")).all())
    print("Start ontology route...")
    startTime=time.time()
    for term in ontology.terms():
        insertTerm(ontology, term, preTerms, preRels, insertedTerms, edges)
    stepOne = time.time()
    print("End ontology route:", stepOne-startTime)
    print("Start DB Load...")
    ontologyName = str(term.id).split(":")[0]
    graph.createNodeIndex(ontName+"_Term", "accession")
    graph.createNodeIndex(ontName+"_Term", "name")
    graph.merge(py4j.SubGraph(set(), edges))
    stepTwo = time.time()
    print("End DB load:", stepTwo-stepOne)

if __name__ == '__main__':
    ont2neo()
