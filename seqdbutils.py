#!/usr/bin/python
"""
seqdbutils.py contains a variety of utility functions for explore taxa,
genInfo ids, and blast. 

 - gi2taxid.py takes an integer (provided on the command-line as a string),
and searches gi_taxid_prot.dmp for that gi identifier, to retrieve the
corresponding taxonomic ID, as an integer. It uses a combination of
_binarysearch() and _linearsearch() because the file is nearly 1Gb, so
a strictly linear search would be painfully slow. 

 - taxid2name is a function to search the BioSQL database for the name of a taxon,
given its integer identifier (integer -> string).

 - genInfo2fasta is a function that retrieves the fasta sequence for any integer gi number (integer -> text).

 - getparenttaxid retrieves the taxon id of the parent of the provided taxon, to enable traversing the taxonomic tree (integer -> integer). 

 - blast2fasta extracts gi numbers from XML blast output, and retrieves the fasta-formatted sequence info for each gi number (XML file -> FASTA file). 

 - taxid2rank retrieves the rank of a taxon (eg. kingdom, order, genus, species, superfamily, parvorder, subspieces, etc.) (integer -> string)

 - taxid2cladename(taxid,level) traverses the taxonomic tree to find the taxonomic name of a particular rank level. (int,string -> string)
    taxid2clade(9606,'genus') -> Homo
    taxid2cladename(9606,'superkingdom') -> Eukaryota

 - browsetaxonomy(taxid) traverses the taxonomic tree to find the names of all the clades containing this taxon (int -> string)

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
            giList.append(int(geneInfoID))
    FastaText = ''
    fastafile = open(fastafilename,'w')
    for gi in giList:
        sys.stderr.write('%d\n'%gi)
        fastafile.write(genInfo2fasta(gi))

def blast2superkingdom(XMLblastfilename):
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
            giList.append(int(geneInfoID))
    for gi in giList:
        #sys.stderr.write('%d\n'%gi2taxid(gi))
        superkingdom = getsuperkingdom(gi2taxid(gi))
        #if superkingdom == 'Eukaryota':
        #    sys.stderr.write('.')
        #else:
        #    sys.stderr.write('\n%s\n'%superkingdom)
        sys.stderr.write('\n%s\n'%superkingdom)
   

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

def taxid2rank(taxid):
    # http://zetcode.com/databases/mysqlpythontutorial/
    rank = 'Error: %d'%taxid
    connection = None
    connection = MySQLdb.connect('localhost', 'root', '','bioseqdb');
    with connection:
        cursor = connection.cursor(MySQLdb.cursors.DictCursor)
        cursor.execute("select taxon_id from taxon where ncbi_taxon_id = %d"%taxid)
        tabletaxid = cursor.fetchone()["taxon_id"]
        cursor.execute("select node_rank from taxon where taxon_id = %s"%tabletaxid)
        row = cursor.fetchone()
        rank = row["node_rank"]
    return rank   

def taxid2cladename(taxid,level):
    """
    Traverse up the tree of life until the desired level is reached, and return the taxonomic name. Examples: 
    taxid2clade(9606,'genus') -> Homo
    taxid2cladename(9606,'superkingdom') -> Eukaryota
    taxid2cladename(71280,'superkingdom') -> Archaea
    taxid2cladename(562,'superkingdom') -> Bacteria
    """
    if taxid2rank(taxid) == level:
        return taxid2name(taxid)
    elif taxid == 1: # taxid 1 corresponds to root
        return 'Error: %s not found'%level
    else:
        return taxid2cladename(getparenttaxid(taxid),level)

def browsetaxonomy(taxid):
    name = taxid2name(taxid)
    rank = taxid2rank(taxid)
    print taxid,':',name,':',rank
    if name != 'root':
        browsetaxonomy(getparenttaxid(taxid))
                    

def getsuperkingdom(taxid):
    return taxid2cladename(taxid,'superkingdom')

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
    maxidx = 987086046 #total chars in file, may be updated in the future.
    dbfile = '/home/edmonds/Documents/Alignments/biosql/gi_taxid_prot.dmp'
    filehandle = open(dbfile,'r')
    taxid = _binarysearch(minidx,maxidx,filehandle,int(gi))
    filehandle.close()
    return int(taxid)

def _linearsearch(minidx,filehandle,quarry):
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
    
def _binarysearch(minidx,maxidx,filehandle,quarry):
    if maxidx-minidx < 1000:
        return _linearsearch(minidx,filehandle,quarry)
    thisidx = minidx+(maxidx-minidx)/2
    filehandle.seek(thisidx)
    filehandle.readline()
    line = filehandle.readline()
    #print minidx,thisidx,maxidx,line[0:-1]
    thisID = int(line.split()[0])
    if quarry < thisID:
        return _binarysearch(minidx,thisidx,filehandle,quarry)
    elif quarry > thisID:
        return _binarysearch(thisidx,maxidx,filehandle,quarry)
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
    #taxid = _binarysearch(0,987086046,database,gi)   
    #database.close()
    print taxid

#main()
