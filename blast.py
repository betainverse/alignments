#!/usr/bin/python
from Bio.Blast import NCBIWWW
from Bio import SeqIO
record = SeqIO.read(open("delta5B.fasta"), format="fasta")
result_handle = NCBIWWW.qblast("blastp", "nr", record.seq,entrez_query="txid2759[orgn]",hitlist_size=500)
#result_handle = NCBIWWW.qblast("blastp", "nr", record.seq)
#result_handle = NCBIWWW.qblast("blastp", "gpipe/2759/ref_contig", record.seq) ##hangs##
save_file = open("2759.xml", "w")
save_file.write(result_handle.read())
save_file.close()
result_handle.close()
#entrez_query="Mus musculus[orgn]",record.format("fasta")
#entrez_query="txid10090[orgn]" is equivalent
# set db as gpipe/10090/ref_contig, for mouse
# eukaryota = taxid:2759
# see http://www.dalkescientific.com/writings/NBN/blast_searching.html
