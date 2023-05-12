# COMSV

### Notes
Please rename the chromosome id in the alignment to be integer before running COMSV. For example, chromosome X to be chromosome 23 and chromosome Y to be chromosome Y. 

COMSV support SV detection for other genomes, but remember to change chromosome id to be integer.

COMSV may have problem to call SVs for the alignments with extremely long molecules/contig (e.g., >100Mbp) or with extremely high depth (e.g., >1000x). Current version is not optimized for these extreme cases as it will take much more memory and running time. Please feel free to leave your comments if you face this issue, and we can specifically modified the program for your needs.
