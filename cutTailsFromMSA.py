#!/usr/bin/python

from seqdbutils import filterFastaColumnsByFirstSeq

infile = "5BhitsculledbygenusUnder80linsi-trimmed.fasta"
outfile = "5BhitsculledbygenusUnder80linsi-trimmed-tailless.fasta"
StartResidue = 627 #First sequence
EndResidue = 1221 #First sequence
filterFastaColumnsByFirstSeq(infile,outfile,StartResidue,EndResidue)
