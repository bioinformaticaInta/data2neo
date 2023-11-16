#!/usr/bin/python

import click
import data2neoFunctions
import multiprocessing as mp
from math import ceil
from itertools import repeat
import socket as sck
import gff2neo
import vcf2neo
import gaf2neo
import ont2neo
import csv2neo

###############################################################################
###############################################################################

@click.group()
@click.option('-d', '--dbname', type=str, required=True, help='Neo4j database name')
@click.option('-u', '--user', type=str, required=True, help='Neo4j user name')
@click.option('-P', '--passwd', type=str, required=True, help='Neo4j user name password')
@click.option('-p', '--port', type=str, required=False, default="7687", help='Neo4j port (default: 7687)')
@click.option('-i', '--ip', type=str, required=False, default=None, help='Neo4j server IP or hostname (default: Local hostname)')
@click.pass_context
def data2neo(ctx, dbname, user, passwd, port, ip):
    '''Load data files to neo4j database.'''
    ctx.obj = dict()
    ctx.obj['dbname'] = dbname
    ctx.obj['user'] = user
    ctx.obj['passwd'] = passwd
    ctx.obj['port'] = port
    if not ip:
        ip = sck.gethostname()
    ctx.obj['ip'] = ip

@data2neo.command()
@click.pass_context
@click.option('-f', '--filename', type=str, required=True, help='GFF3 file name')
@click.option('-s', '--structured', is_flag=True, help='GFF3 file has structure with levels')
@click.option('-e', '--loadexons', is_flag=True, help='Load exons as features. Only for structured GFF3 files')
@click.option('-c', '--chromfeature', is_flag=True, help='Chromosome present as feature in GFF3 file')
def gff(ctx, filename, structured, loadexons, chromfeature):
    '''Load GFF3 files'''
    gff2neo.gff2neo(ctx.obj, filename, structured, loadexons, chromfeature)

@data2neo.command()
@click.pass_context
@click.option('-f', '--filename', type=str, required=True, help='VCF file name')
@click.option('-n', '--nprocs', type=int, required=False, default=None, help='Number of processes per chunk (default: ALL)')
@click.option('-l', '--nlines', type=int, required=False, default=100000, help='Number of lines per chunk (default: 100000)')
@click.option('-c', '--createchroms', is_flag=True, help='Create unknown chromosome nodes')
def vcf(ctx, filename, nprocs, nlines, createchroms):
    '''Load VCF files'''
    nprocsMax = mp.cpu_count()
    if not nprocs or nprocs > nprocsMax:
        nprocs = nprocsMax
    #Insert genotypes
    vcf2neo.createGenotypes(ctx.obj, filename)
    data2neoFunctions.splitWithHeaderLines(filename, "#", nlines, "allManifest")
    try:
        with open("allManifest") as mfFileAll:
            next(mfFileAll)
            nPool = 1
            for mfLineAll in mfFileAll:
                mfLineAll = mfLineAll[:-1].split(",")
                chroms = data2neoFunctions.splitWithHeaderProcs(mfLineAll[0], "#", nprocs, "manifest"+str(nPool))
                #Insert Chromosomes
                vcf2neo.insertChroms(ctx.obj, chroms, createchroms)
                mfFile = open("manifest"+str(nPool))
                mfFileNames = (mfLine.split(',')[0] for mfLine in mfFile.readlines()[1:])
                mfFile.close()
                print("Running pool number:",nPool)
                poolOfProcs = mp.Pool(nprocs)
                listOfOuts = poolOfProcs.starmap(vcf2neo.vcf2neo, zip(repeat(ctx.obj),mfFileNames))
                poolOfProcs.close()
                poolOfProcs.join()
                data2neoFunctions.loadVCFDataToGraph(ctx.obj, listOfOuts)
                nPool += 1
    except Exception as e:
        print("Error in proccessing")
        print(str(e))
    finally:
        for n in range(1,nPool): data2neoFunctions.removeSplitFiles("manifest"+str(n))
        data2neoFunctions.removeSplitFiles("allManifest")

@data2neo.command()
@click.pass_context
@click.option('-f', '--filename', type=str, required=True, help='Ontology file name or URL')
def ont(ctx, filename):
    '''Load ontologies (OBO files)'''
    ont2neo.ont2neo(ctx.obj, filename)

@data2neo.command()
@click.pass_context
@click.option('-f', '--filename', type=str, required=True, help='GAF file name')
@click.option('-s', '--objectsymbol', is_flag=True, help='Use DB object symbol(column 3) instead DB object ID (column 2) as gene number')
def gaf(ctx, filename, objectsymbol):
    '''Load functional annotation (GAF files)'''
    nameFieldNumber = 1
    if objectsymbol: nameFieldNumber = 2
    gaf2neo.gaf2neo(ctx.obj, filename, nameFieldNumber)

@data2neo.group()
@click.pass_context
@click.option('-f', '--filename', type=str, required=True, help='CSV file name')
def csv(ctx, filename):
    '''Load nodes or relationships from CSV files'''
    ctx.obj['file'] = filename

@csv.command()
@click.pass_context
@click.option('-l', '--label', type=str, required=True, default=None, help='Node label name (coma separated for multiple labels')
def nodes(ctx, label):
    '''Load nodes from CSV files\n
    Existing nodes will be modified\n\n
    CSV file format:\n
    header: prop1Name:type , prop2Name:type , ...... , propLName:type\n
    line1:  prop1Value , prop2Value , ...... , propLValue\n
                           .........\n
    lineM:  prop1Value , prop2Value , ...... , propLValue\n
    prop1 (first column) in used as primary property and index creation\n
    Accepted types: date(YYYY-MM-DD), bool, int, float and str
    '''
    csv2neo.csvNodes2neo(ctx.obj, label.split(','))

@csv.command()
@click.pass_context
@click.option('-s', '--labelstart', type=str, required=True, help='Label of start node')
@click.option('-e', '--labelend', type=str, required=True, help='Label of end node')
@click.option('-t', '--typerel', type=str, required=True, help='Relationship type name')
@click.option('-c', '--createnodes', is_flag=True, help='Create unknown nodes')
def rels(ctx, labelstart, labelend, typerel, createnodes):
    '''Load relationships from CSV files\n\n
    CSV file format:\n
    header: propStartName:type , proprel1Name:type, propRel2Name:type , ...... , propRelLName:type , propEndName:type\n
    line1:  propStartValue , proprel1Value, propRel2Value , ...... , propRelLValue , propEndValue\n
                             .........\n
    lineM:  propStartValue , proprel1Value, propRel2Value , ...... , propRelLValue , propEndValue\n
    Accepted types: date(YYYY-MM-DD), bool, int, float and str
    '''
    data2neoFunctions.csvRelsHeaderFormating(ctx.obj['file'])
    try:
        csv2neo.csvRels2neo(ctx.obj, labelstart, labelend, typerel, createnodes)
    except Exception as e:
        print("Error in file, data not loaded to graph database")
        print(str(e))
    finally:
        data2neoFunctions.csvRelsHeaderReFormating(ctx.obj['file'])

if __name__ == '__main__':
    data2neo()
