# COMSV
COMSV (Cancer Optical Mapping based Structural Variation detection) is a pipeline for SV detection based on nanochannel optical mapping data. COMSV accepts .oma (OM Alignment) format files as alignment results and .cmap (consensus map) file as refrence genome maps.

Copyright &copy; 2023  Le Li (<ll863@cornell.edu>)

![alt text](https://github.com/kevingroup/COMSV/blob/main/FigPipelines.png?raw=true)

The COMSV pipelines, which include the indel pipeline (**a-d**) and the complex SV pipeline (**e-f**). **a** Observed-to-expected distance ratios between two labels can deviate from 1 due to scaling or alignment errors (first row). After normalization, for alignments that still have these abnormal ratios (second row, dashed rectangles), the mutually reinforcing ones are considered indel candidates while the others are discarded (third row). Isolated label pairs with abnormal distances are instead corrected (fourth row, green bar). **b** Adjacent label pairs are collected from individual molecules to form candidate indel regions based on their overlaps, which enables the collection of distances between non-adjacent labels from a molecule if they are adjacent in another molecule (shown in green). **c** For each label pair, the collected distances are clustered to identify the number of distinct alleles and the number of molecules that support each of them. **d** For each potential indel region, the clustering result of each label pair suggests an initial SV type, and the final SV type of the whole region is determined by considering these suggestions jointly. **e** Different types of SVs that can be identified from split-alignments. An SV break point is called if a split alignment suggests an SV event, while a complete SV is called if the full span of the SV can be determined. **f** An example of calling SV break points and complete SVs. The reference is repeatedly displayed multiple times to show the split-alignment of each molecule clearly. Vertical dotted lines mark key aligned labels that define the boundary of each aligned segment on the reference. The corresponding labels on the reference are shown in red. Break point information is first collected from individual molecules and then considered together to refine the break point locations and determine whether complete SVs can be called.

## License
See LICENSE.

## Installation
No installation is needed. 
The pipeline is constructed with c++. 
Just download the whole package and compile all c++ source files with the script 

    bash compile_source.sh
    
Remember to give **execute** permission to all compiled c++ programs (`chmod +x`).

## One-line command
We provide a one-line command for a quick start with example data (alignments of sample 261T on chromosome 1). Please unzip the data first before running the pipeline:

    gunzip 261T*

After that, simply type 

    bash run_all.sh
    
After progrem stops, check the output folder `COMSV_SVs_261T` for the output SV lists (indels_dedup.osv, invdup_dedup.osv, and translocation_dedup.osv).

## Different SV detection components
There are several components to call SVs in our pipeline:

    run_Comsv_indel.sh: this component call indels based on molecule alignments with c++ program `COMSV_indel`. In the header of this script, you can change the sample name, label, alignment filename, and reference filename for your own use. A simple use of `COMSV_indel` with default parameters is included. Please use `./COMSV_indel -help` for help of setting parameters of `COMSV_indel`.
    
    run_Comsv_contig.sh: this component call indels based on contig alignments with c++ program `COMSV_indel_contig`. In the header of this script, you can change the sample name, label, alignment filename, and reference filename for your own use. A simple use of `COMSV_indel_contig` with default parameters is included. Please use `./COMSV_indel_contig -help` for help of setting parameters of `COMSV_indel_contig`.
    
    run_Comsv_Complex.sh: this component call complex SVs based on molecule alignments, contig alignments, as well as two-round split alignment with c++ program `COMSV_complex`. This script outputs three lists of complex SVs, which can be integrated into one list by using `assemble_SVs.sh`. In the header of this script, you can change the sample name, label, alignment filename, and reference filename for your own use. A simple use of `COMSV_complex` with default parameters is included. Please use `./COMSV_complex -help` for help of setting parameters of `COMSV_complex`.

    assemble_SVs.sh: this script does postprocessing and integrates the SVs called by the previous components into three lists of SVs: indel list, invdup list, and translocation list. Remember to change the sample name `smp` to the one used in SV detections.

## Predicting scores of indels called from molecule alignments
`score_pred.py` provides score prediction for indels. It requires the python package `sklearn`.

    python score_pred.py INPUTFILE OUTPUTFILE
    
**INPUTFILE**: the name of tab-temilited input file, which contains three columns: Size of SVs, Support rate of SVs, and Coverage of SVs. It must contain one header line with the exact names: "Size", "Support", and "Coverage". Size should be the absolute value of Size1.

**OUTPUTFILE**: the output file contains Score1 (score of SV without considering zygosity) and Score2 (score of SV considering zygosity) for each SV.
 
## Tips
Please rename the chromosome IDs in the alignment to be **positive integers** before running COMSV. For example, chromosome X to be chromosome 23 and chromosome Y to be chromosome Y. 

COMSV support SV detection for other genomes, but remember to change chromosome IDs to be integers.

COMSV may have problem to call SVs for the alignments with extremely long molecules/contig (e.g., >100Mbp) or with extremely high depth (e.g., >1000x). Current version is not optimized for these extreme cases as it will take much more memory and running time. Please feel free to leave your comments if you face this issue, and we can specifically modified the program for your needs.

## Data availability
All simulation data and real cancer data used in the following paper are available upon request.

## Citation

    Le Li*, Chenyang Hong*, Jie Xu*, Claire Yik-Lok Chung, Alden King-Yung Leung, Delbert Almerick T. Boncan, Lixin Cheng, Kwok-Wai Lo, Paul B. S. Lai, John Wong, Jingying Zhou, Alfred Sze-Lok Cheng, Ting-Fung Chan, Feng Yue# and Kevin Y. Yip#. Accurate identication of structural variations from cancer samples. Under review, 2023.
