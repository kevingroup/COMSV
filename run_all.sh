mkdir temp

bash compile_source.sh
bash run_Comsv_indel.sh
bash run_Comsv_indel_contig.sh
bash run_Comsv_Complex.sh
bash assemble_SVs.sh
rm -r temp
