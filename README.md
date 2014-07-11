alignments
==========

This repository contains a collection of python scripts for managing the
results of Blast searches. 

blast.py is a simple script for running blast searches. There may be problems with limiting a search by taxonomy id. 

seqdbutils.py contains several handy functions:

 - gi2taxid is function to search the file gi_taxid_prot.dmp for a gi identifier, 
and find the taxonomic ID. (integer -> integer)

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

 - filterFastaColumnsByFirstSeq(infile,outfile,startres,endres) removes columns from a multiple sequence alignment, preserving only those corresponding to a given residue range in the first sequence. 

To do:

Parse CSV of blast output
Calculate %coverage from q.start and q.end
fields: query id, subject ids, % identity, % positives, alignment length,
mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
somehow parse subject IDs to determine whether you've seen that species
before. (gi2taxid helps) If so, compare coverage to prior sequence from same species. If
coverage is more, replace prior sequence in dictionary keyed by species. 

Useful links:
[1] ftp://ftp.ncbi.nih.gov/pub/taxonomy/

http://www.mathworks.com/help/bioinfo/examples/performing-a-metagenomic-analysis-of-a-sargasso-sea-sample.html
