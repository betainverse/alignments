#!/usr/bin/python
"""
gi2taxid.py takes an integer (provided on the command-line as a string),
and searches gi_taxid_prot.dmp for that gi identifier, to retrieve the
corresponding taxonomic ID.
"""

import sys
import MySQLdb
from Bio import Entrez
import xml.etree.ElementTree as ET

# def slowgi2taxid(gi,database):
#     gistr = str(gi)
#     taxid = -1
#     for line in database:
#         columns = line.split()
#         if gistr == columns[0]:
#             return columns[1]
#     return taxid

#seek min line length is ~7, max is 15
#59630862 lines
#987086046 characters
#1000 lines in test
#9835 characters

def blast2fasta(XMLblastfilename,fastafilename):
    tree = ET.parse(XMLblastfilename)
    root = tree.getroot()
    iterations = root.find('BlastOutput_iterations')
    Iteration = iterations.find('Iteration')
    Iteration_hits = Iteration.find('Iteration_hits')
    hits = Iteration_hits.findall('Hit')
    giList = []
    for hit in hits:
        geneInfoID = None
        accept = True
        hitnum = hit.find('Hit_num').text
        title = hit.find('Hit_def').text
        hitids = hit.find('Hit_id').text.split('|')
        if hitids[0] == 'gi':
            geneInfoID = hitids[1]
            giList.append(geneInfoID)
    FastaText = ''
    fastafile = open(fastafilename,'w')
    for gi in giList:
        sys.stderr.write(gi+'\n')
        fastafile.write(genInfo2fasta(gi))

def taxid2name(taxid):
    # http://zetcode.com/databases/mysqlpythontutorial/
    scientificName = 'Error: %d'%taxid
    connection = None
    connection = MySQLdb.connect('localhost', 'root', '','bioseqdb');
    with connection:
        cursor = connection.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute("select taxon_id from taxon where ncbi_taxon_id = %d"%taxid)
        tabletaxid = cursor.fetchone()["taxon_id"]
        cursor.execute("select * from taxon_name where taxon_id = %s and name_class = \"scientific name\""%tabletaxid)
        row = cursor.fetchone()
        scientificName = row["name"]
    return scientificName

def getparenttaxid(taxid):
    scientificName = 'Error: %d'%taxid
    connection = None
    connection = MySQLdb.connect('localhost', 'root', '','bioseqdb');
    with connection:
        cursor = connection.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute("select parent_taxon_id from taxon where ncbi_taxon_id = %d"%taxid)
        tabletaxid = cursor.fetchone()["parent_taxon_id"]
        cursor.execute("select ncbi_taxon_id from taxon where taxon_id = %s"%tabletaxid)
        row = cursor.fetchone()
        parenttaxid = int(row["ncbi_taxon_id"])
    return parenttaxid

def genInfo2fasta(gi):
    Entrez.email = "kedmonds@bu.edu"   
    net_handle = Entrez.efetch(db="nucleotide", id=gi,rettype="fasta", retmode="text")
    return net_handle.read()

def gi2taxid(gi):
    minidx = 0
    maxidx = 987086046 #total chars in file
    dbfile = '/home/edmonds/Documents/Alignments/biosql/gi_taxid_prot.dmp'
    filehandle = open(dbfile,'r')
    taxid = searchfile(minidx,maxidx,filehandle,gi)
    filehandle.close()
    return taxid

def linearsearch(minidx,filehandle,quarry):
    startidx = max(minidx-15,0)
    filehandle.seek(startidx)
    filehandle.readline()
    line = filehandle.readline()
    #print minidx,line[0:-1]
    thisID = int(line.split()[0])
    while thisID < quarry:
        line = filehandle.readline()
        #print line[0:-1]
        thisID = int(line.split()[0])
    if quarry == thisID:
        return int(line.split()[1])
    else:
        return -1
    
def searchfile(minidx,maxidx,filehandle,quarry):
    if maxidx-minidx < 1000:
        return linearsearch(minidx,filehandle,quarry)
    thisidx = minidx+(maxidx-minidx)/2
    filehandle.seek(thisidx)
    filehandle.readline()
    line = filehandle.readline()
    #print minidx,thisidx,maxidx,line[0:-1]
    thisID = int(line.split()[0])
    if quarry < thisID:
        return searchfile(minidx,thisidx,filehandle,quarry)
    elif quarry > thisID:
        return searchfile(thisidx,maxidx,filehandle,quarry)
    elif quarry == thisID:
        return int(line.split()[1])
    

def main():
    #if len(sys.argv) < 3:
    #    dbfile = '/home/edmonds/Documents/Alignments/biosql/gi_taxid_prot.dmp'
        #dbfile = '/home/edmonds/Documents/Alignments/biosql/test.dmp'
    #else:
    #    dbfile = sys.argv[2]
    if len(sys.argv) < 2:
        gi = 288
        gi = 301769165
    else:
        gi = int(sys.argv[1])
    #database = open(dbfile,'r')
    taxid = gi2taxid(gi)
    #taxid = searchfile(0,987086046,database,gi)   
    #database.close()
    print taxid

#main()
