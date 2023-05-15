# COMSV
COMSV (Cancer Optical Mapping based Structural Variation detection) is a pipeline for SV detection based on nanochannel optical mapping data. COMSV accepts .oma (OM Alignment) format files as alignment results and .cmap (consensus map) file as refrence genome maps.

Copyright &copy; 2023  Le Li (<ll863@cornell.edu>)
## License
See LICENSE.

## Installation
No installation is needed. 
The pipeline is constructed with c++. 
Just download the whole package and compile all c++ source files with the script 

    bash compile_source.sh
    
Remember to give **execute** permission to all compiled c++ programs (`chmod +x`).

## One-line command
We provide a one-line command for a quick start with example data (alignments of sample 261T on chromosome 1).
In the path of the package, simply type 

    bash run_all.sh
    
After progrem stops, check the output folder `COMSV_SVs_261T` for the output SV lists.

## Different SV detection components
There are several components to call SVs in our pipeline:

    run_Comsv_indel.sh: this component call indels based on molecule alignments with c++ program `COMSV_indel`. In the header of this script, you can change the sample name, label, alignment filename, and reference filename for your own use. A simple use of `COMSV_indel` with default parameters is included. Please use `./COMSV_indel -help` for help of setting parameters of `COMSV_indel`.
    
    run_Comsv_contig.sh: this component call indels based on contig alignments with c++ program `COMSV_indel_contig`. In the header of this script, you can change the sample name, label, alignment filename, and reference filename for your own use. A simple use of `COMSV_indel_contig` with default parameters is included. Please use `./COMSV_indel_contig -help` for help of setting parameters of `COMSV_indel_contig`.
    
    run_Comsv_Complex.sh: this component call complex SVs based on molecule alignments, contig alignments, as well as two-round split alignment with c++ program `COMSV_complex`. This script outputs three lists of complex SVs, which can be integrated into one list by using `assemble_SVs.sh`. In the header of this script, you can change the sample name, label, alignment filename, and reference filename for your own use. A simple use of `COMSV_complex` with default parameters is included. Please use `./COMSV_complex -help` for help of setting parameters of `COMSV_complex`.

    assemble_SVs.sh: this script does postprocessing and integrates the SVs called by the previous components into three lists of SVs: indel list, invdup list, and translocation list. Remember to change the sample name `smp` to the one used in SV detections.
    
## Tips
Please rename the chromosome IDs in the alignment to be **positive integers** before running COMSV. For example, chromosome X to be chromosome 23 and chromosome Y to be chromosome Y. 

COMSV support SV detection for other genomes, but remember to change chromosome IDs to be integers.

COMSV may have problem to call SVs for the alignments with extremely long molecules/contig (e.g., >100Mbp) or with extremely high depth (e.g., >1000x). Current version is not optimized for these extreme cases as it will take much more memory and running time. Please feel free to leave your comments if you face this issue, and we can specifically modified the program for your needs.

## Citation

    Le Li*, Chenyang Hong*, Jie Xu*, Claire Yik-Lok Chung, Alden King-Yung Leung, Delbert Almerick T. Boncan, Lixin Cheng, Kwok-Wai Lo, Paul B. S. Lai, John Wong, Jingying Zhou, Alfred Sze-Lok Cheng, Ting-Fung Chan, Feng Yue# and Kevin Y. Yip#. Accurate identication of structural variations from cancer samples. Under review, 2023.
