#!/usr/bin/python

import vcf as VCF
import multiprocessing as mp
import py4j
import data2neoFunctions

###############################################################################
###############################AUXILIAR FUNCTIONS##############################
###############################################################################

def insertVariant(feature, chrom, edges, indexDataSet):
    #Extraer y capitalizar el tipo para usarlo como label del vertice
    labels = (feature.var_type.title().replace('_',''), "Variant")
    #Create variant vertex / TODO: Agregar mas datos???
    name = feature.ID
    if not name:
        name = str(feature.CHROM) + "_" + str(feature.POS)
    vertex = py4j.Node(*labels, name = name)
    #Insert edge with genomic localization
    edges.add(py4j.Relationship(vertex , chrom, "LOCATED_AT", position=int(feature.POS)))
    
    indexDataSet.add((labels[0], 'name'))

    #Insert alleles
    alleles = {}
    allAlleles = [feature.REF] + feature.ALT
    for number,var in enumerate(allAlleles):
        allele = py4j.Node("Allele", name = name+"_"+str(number), reference = (number==0), variation = str(var))
        #Connect allele to variant
        edges.add(py4j.Relationship(allele , vertex, "VARIANT_OF"))
        alleles[number] = allele
    return alleles

def insertGenotypeAlleles(genotypeData, genotype, edges, alleles):
    #Connect genotype to alleles
    alleleToGenotypes = {}
    for alleleNumber in genotypeData.split('|'):
        try:
            alleleNumber = int(alleleNumber)
            if alleleNumber in alleleToGenotypes:
                alleleToGenotypes[alleleNumber] = py4j.Relationship(alleles[alleleNumber], genotype, "PRESENT_IN", homozygous=True)
            else:
                alleleToGenotypes[alleleNumber] = py4j.Relationship(alleles[alleleNumber], genotype, "PRESENT_IN", homozygous=False)
        except:
            pass
    edges |= set(alleleToGenotypes.values())

###############################################################################
###############################################################################
###############################################################################

def vcf2neo(graphData, vcfFile):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])
    chroms = set(graph.nodeMatch("Chromosome").all())
    genotypes = set(graph.nodeMatch("Genotype").all())
    
    vcf = VCF.Reader(open(vcfFile))
    dataForIndexes = {("Variant","name")}
    edges = set()
    chrom = None
    for feature in vcf:
        #Se inserta variante y retorna diccionario de vertices de alelos
        if not chrom or chrom["name"] != feature.CHROM:
            chrom = data2neoFunctions.firstVertexInSet(chroms, props={"name":feature.CHROM})
        alleles = insertVariant(feature, chrom, edges, dataForIndexes)
        #Se inserta informacion de genotipos
        for genotypeData in feature.samples:
            genotype = data2neoFunctions.firstVertexInSet(genotypes, props={"name":genotypeData.sample})
            insertGenotypeAlleles(genotypeData['GT'], genotype, edges, alleles)
    return edges, dataForIndexes

def createGenotypes(graphData, vcfFile):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])
    vcf = VCF.Reader(open(vcfFile))
    genotypes = set()
    for sample in vcf.samples:
        genotypes.add(py4j.Node("Genotype", name = sample))
    graph.createNodeIndex("Genotype", "name") 
    graph.merge(py4j.SubGraph(genotypes,set()))

def insertChroms(graphData, chromsSet, createChrom=True):
    graph = py4j.Graph("neo4j://"+graphData['ip']+":"+graphData['port'], graphData['dbname'], graphData['user'], graphData['passwd'])

    chroms = set()
    for chromName in chromsSet:
        chrom = graph.nodeMatch("Chromosome", name=chromName).first()
        if not chrom:
            if createChrom:
                chroms.add(py4j.Node("Chromosome", name = chromName))
            else:
                raise Exception("The chromosome %s is not in the database" % (chromName))
    graph.createNodeIndex("Chromosome", "name") 
    graph.create(py4j.SubGraph(chroms,set()))

if __name__ == '__main__':
    vcf2neo()
