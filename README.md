alignments
==========

This repository contains a collection of python scripts for managing the
results of Blast searches. 

blast.py is a simple script for running blast searches. There may be problems with limiting a search by taxonomy id. 

seqdbutils.py contains several handy functions:

 - gi2taxid is function to search the file gi_taxid_prot.dmp for a gi identifier, 
and find the taxonomic ID.

 - taxid2name is a function to search the BioSQL database for the name of a taxon,
given its identifier

 - genInfo2fasta is a function that retrieves the fasta sequence for any gi number.

 - getparenttaxid retrieves the taxon id of the parent of the provided taxon, to enable traversing the taxonomic tree. 

Useful links:
[1] ftp://ftp.ncbi.nih.gov/pub/taxonomy/

http://www.mathworks.com/help/bioinfo/examples/performing-a-metagenomic-analysis-of-a-sargasso-sea-sample.html