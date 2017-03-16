MLSTtrimm
=========

Python Scripts for inferring Allele trimming patterns for MLST alleles 

The goal of this scripts is to provide internal trimming patterns as IUPAC codes to MLST alleles given a fasta file of alleles for an MLST locus. The user as an option to select a minimmum threshold value to which a base in a given position should be considered and the number of bases for both the 5' and 3' trimming patterns

#Usage#

    MLST_trimpattern.py <fasta_filename> <threshold (0-1)> <pattern size> 

#Description#


Reads a fasta file and outputs the trimming pattern sequences in IUPAC code by the selected size and threshold."

##Developers##

João André Carriço

Jorge Diamantino Miranda 

jcarrico@fm.ul.pt, jmiranda@medicina.ulisboa.pt  30/05/2014"
