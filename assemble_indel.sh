smp="261T"
outputFolder="COMSV_SVs_$smp"

sed '/#/d' $outputFolder/Complex_contig.rosv | sed '/###/d' | awk '$16="Contig-Alignment"' OFS='\t' > $outputFolder/Complex_All_${smp}.rosv
sed '/#/d' $outputFolder/Complex.rosv | sed '/###/d' | sed '/Duplication/!d' | awk '$16="Single-Molecule-Alignment"' OFS='\t' >> $outputFolder/Complex_All_${smp}.rosv
sed '/#/d' $outputFolder/Complex_2step.rosv | sed '/###/d' | sed '/-Split/d' | sed '/Large-/d' | awk '$16="Two-round-split-alignment"' OFS='\t' >> $outputFolder/Complex_All_${smp}.rosv
bash rosv_to_osv.sh $outputFolder/Complex_All_${smp}.rosv tp1_${smp}.osv $outputFolder/translocation_${smp}.rosv
sort -n -k1,1 -k2,2 -k3,3 tp1_${smp}.osv > tp2_${smp}.osv
./Deduplicate_osv_complex -inputSV tp2_${smp}.osv -outputFile compiled_SVs_${smp}.osv &> /dev/null

sed '/#/d' $outputFolder/translocation_${smp}.rosv | cut -f2- | sed 's/\.\.\./\./g' | sed 's/\/Deletion//g' > translocation_${smp}.txt
./Deduplicate_osv_trans -inputSV translocation_${smp}.txt -outputFile $outputFolder/translocation_dedup.osv
#rm translocation_${smp}_dedup_unfilter.osv


    ctof=30
    rco=5
    sed '/#/d' $outputFolder/Indel_All.osv | awk '$10="Single-Molecule-Alignment"' OFS='\t' | awk -v var=$rco '$13>=var' > Indels_${smp}.osv
##### here remove assembly indels
    sed '/#/d' $outputFolder/Indel_Contig_All.osv | awk '$10="Contig-Alignment"' OFS='\t' >> indels_${smp}.osv

    # Complement Indels from complex SV pipeline
    awk '$4=="Insertion"' compiled_SVs_${smp}.osv | awk '$10="Complex-Complement"' OFS='\t' >> indels_${smp}.osv
    awk '$4=="Deletion"' compiled_SVs_${smp}.osv | sed 's/sizeChange=/sizeChange=-/g' | awk '$10="Complex-Complement"' OFS='\t' >> indels_${smp}.osv
    sort -n -k1,1 -k2,2 -k3,3 indels_${smp}.osv > indels_added_${smp}.osv
    # Detect inversions and complement indels from inversion rescue
    ./infer_inv -inputSVFile indels_added_${smp}.osv -outputFile extra_inv_${smp}.bed &> /dev/null
    sed '/#/d' extra_inv_${smp}.bed | awk -F"\t" '$6="-"' OFS="\t" | awk -F"\t" '$7="-"' OFS="\t" | awk '$10="Indels-Inferred"' OFS='\t' >> compiled_SVs_${smp}.osv
    sort -n -k1,1 -k2,2 -k3,3 compiled_SVs_${smp}.osv > tp2_${smp}.osv
    ./Deduplicate_osv_complex -inputSV tp2_${smp}.osv -outputFile compiled_SVs_${smp}.osv &> /dev/null

    # Filter inversions with duplications
    grep "Duplication" compiled_SVs_${smp}.osv > duplication_${smp}.osv
    ./compBeds_full_new -inputFSFile duplication_${smp}.osv -inputSVFile compiled_SVs_${smp}.osv -outputFile anti_dup_val_${smp}.bed > /dev/null
    sed '/#/d' anti_dup_val_${smp}.bed | cut -f12 > anti_dup_val_${smp}.txt
    paste compiled_SVs_${smp}.osv anti_dup_val_${smp}.txt > compiled_SVs_val_${smp}.txt
    grep '#' compiled_SVs_${smp}.osv > compiled_SVs_${smp}.txt
    awk '$(NF)==0' OFS='\t' compiled_SVs_val_${smp}.txt | rev | cut -f2- | rev >> compiled_SVs_${smp}.txt
    mv compiled_SVs_${smp}.txt compiled_SVs_${smp}.osv
    ./Deduplicate_osv_complex -inputSV compiled_SVs_${smp}.osv -outputFile $outputFolder/compiled_SVs_dedup.osv -mode 3 &> /dev/null
    rm compiled_SVs_dedup_${smp}_unfilter.osv

    # Filter indels with complex SV callset
    sed 's/sizeChange=//g' compiled_SVs_${smp}.osv | grep 'Inversion\|Duplication' | awk '$8<500000' | awk '$8="sizeChange="$8' OFS='\t' | sort -nk1,1 -k2,2 -k3,3 > invdup_filter_${smp}.osv
    ./compBeds_full_new -inputFSFile invdup_filter_${smp}.osv -inputSVFile indels_added_${smp}.osv -outputFile indels_added_val_${smp}.bed > /dev/null
    sed '/#/d' indels_added_val_${smp}.bed | cut -f12 > indel_val_${smp}.txt
    paste indels_added_${smp}.osv indel_val_${smp}.txt > indels_added_val_${smp}.txt
    awk '$(NF)==0' OFS='\t' indels_added_val_${smp}.txt | rev | cut -f2- | rev > indels_val_${smp}.txt
    ./Deduplicate_osv_indel -inputSV indels_val_${smp}.txt -outputFile indels_${smp}.osv &> /dev/null
    ./Deduplicate_osv_indel -inputSV indels_val_${smp}.txt -outputFile $outputFolder/indels_dedup.osv -mode 3 &> /dev/null

    rm *${smp}*.osv
    rm *${smp}*.txt
    rm *${smp}*.rosv
