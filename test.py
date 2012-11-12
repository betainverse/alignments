#!/usr/bin/python
#from BioSQL import BioSeqDatabase
#import MySQLdb as mdb
from seqdbutils import gi2taxid,taxid2name,genInfo2fasta,getparenttaxid,blast2fasta,taxid2rank,taxid2cladename,blast2superkingdom,browsetaxonomy
#server = BioSeqDatabase.open_database(driver="MySQLdb", user="root",
#                     passwd = "", host = "localhost", db="bioinfo")

#gilist = [84043963,355751524,355565935,387763418,114579120,73969383,40788346,
#          332251518,5002645,301769165,410352203,168267520]

#taxidlist = [gi2taxid(x) for x in gilist]
#print taxidlist

#namelist = [taxid2name(x) for x in taxidlist]
#print namelist

#parentlist = [taxid2name(getparenttaxid(getparenttaxid(getparenttaxid(getparenttaxid(x))))) for x in taxidlist]
#print parentlist

#ranklist = [taxid2rank(x) for x in taxidlist]
#print ranklist

#browsetaxonomy(9606)
browsetaxonomy(32630)
#print taxid2cladename(32630,'superkingdom') # Error - synthetic

#blast2superkingdom('/home/edmonds/Documents/5B/blast/automated/2759.xml')

# print taxid2cladename(71280,'superkingdom')
# print taxid2cladename(562,'superkingdom')
# print taxid2cladename(9606,'superkingdom')
# print taxid2cladename(9606,'genus')
# print getsuperkingdom(562)

#inputfile = '5BhitsculledbygenusUnder80.xml'
#outputfile = '5BhitsculledbygenusUnder80.fasta'
#blast2fasta(inputfile,outputfile)

#outputfile = open('test.out','w')

#for gi in gilist:
#    outputfile.write(genInfo2fasta(gi))

#outputfile.close()

# con = None

# con = mdb.connect('localhost', 'root', '','bioseqdb');
# with con: 
#     cur = con.cursor()
#     cur.execute("select * from taxon_name where taxon_id = 4092 and name_class = \"scientific name\"")
#     #cur.execute("select * from taxon_name where taxon_id = 4092")
#     rows = cur.fetchall()
    
#     for row in rows:
#         print row[1]

    #data = cur.fetchone()
    

# for identifier in ['6273291', '6273290', '6273289'] :
#     seq_record = db.lookup(gi=identifier)
#     print seq_record.id, seq_record.description[:50] + "..."
#     print "Sequence length %i," % len(seq_record.seq)


