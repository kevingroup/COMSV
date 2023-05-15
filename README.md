# COMSV
COMSV (Cancer Optical Mapping based Structural Variation detection) is a pipeline for SV detection based on nanochannel optical mapping data.

Copyright &copy; 2023  Le Li (<ll863@cornell.edu>)
## License
See LICENSE.

## Installation
No installation is needed. 
The pipeline is constructed with c++. 
Just download the whole package and compile all c++ source files with the script `bash compile_source.sh`.
Remember to give **execute** permission to all compiled c++ programs (`chmod +x`).

## One-line command
We provide a one-line command for a quick start with example data (alignments of sample 261T on chromosome 1).
In the path of the package, simply type `bash run_all.sh`. After progrem stops, check the output folder `COMSV_SVs_261T` for the output SV lists.

## Different SV detection components
There are several components to call SVs in our pipeline:

    dfdfdfdf

## Notes
Please rename the chromosome IDs in the alignment to be positive integers before running COMSV. For example, chromosome X to be chromosome 23 and chromosome Y to be chromosome Y. 

COMSV support SV detection for other genomes, but remember to change chromosome IDs to be integers.

COMSV may have problem to call SVs for the alignments with extremely long molecules/contig (e.g., >100Mbp) or with extremely high depth (e.g., >1000x). Current version is not optimized for these extreme cases as it will take much more memory and running time. Please feel free to leave your comments if you face this issue, and we can specifically modified the program for your needs.
