#!/usr/bin/python
from Bio import AlignIO
from seqdbutils import filterFastaColumnsByFirstSeq
#alignment = AlignIO.read("5BhitsculledbygenusUnder80clustalo.fst", "fasta")
#print _findFastaPosWOdashes(alignment[0].seq,626)
#print alignment[:20,2342:2352]
#print alignment[:20,2342:2984]

filterFastaColumnsByFirstSeq("5BhitsculledbygenusUnder80clustalo.fst","5BhitsculledbygenusUnder80clustalo627-1220.fst",627,1220)
filterFastaColumnsByFirstSeq("5BhitsculledbygenusUnder80clustalo.fst","5BhitsculledbygenusUnder80clustalo951-1085.fst",951,1085)
