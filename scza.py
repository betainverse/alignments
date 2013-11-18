#!/usr/bin/env python
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from seqdbutils import filterblast2fasta,simplifyFASTAtitles

#record = SeqIO.read(open("pneumo-scza.fasta"), format="fasta")
#result_handle = NCBIWWW.qblast("blastp", "refseq_protein", record.seq,hitlist_size=500)
#save_file = open("scza.xml", "w")
#save_file.write(result_handle.read())
#save_file.close()
#result_handle.close()

XMLblastfilename = 'scza.xml'
fastafilename = 'scza.fasta'
filterblast2fasta(XMLblastfilename,fastafilename)

simplifyFASTAtitles(fastafilename,'outputscza.fasta')

