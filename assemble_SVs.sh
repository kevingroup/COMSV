smp="261T"
outputFolder="COMSV_SVs_$smp"

rm -f $outputFolder/Complex_All_${smp}.rosv
if [ -f $outputFolder/Complex_contig.rosv ]
then
    sed '/#/d' $outputFolder/Complex_contig.rosv | sed '/###/d' | awk '$16="Contig-Alignment"' OFS='\t' > $outputFolder/Complex_All_${smp}.rosv
fi
if [ -f $outputFolder/Complex.rosv ]
then
    if [ -f $outputFolder/Complex_All_${smp}.rosv ]
    then
        sed '/#/d' $outputFolder/Complex.rosv | sed '/###/d' | sed '/Duplication/!d' | awk '$16="Single-Molecule-Alignment"' OFS='\t' >> $outputFolder/Complex_All_${smp}.rosv
    else
        sed '/#/d' $outputFolder/Complex.rosv | sed '/###/d' | sed '/Duplication/!d' | awk '$16="Single-Molecule-Alignment"' OFS='\t' > $outputFolder/Complex_All_${smp}.rosv
    fi
fi
if [ -f $outputFolder/Complex_2step.rosv ]
then
    if [ -f $outputFolder/Complex_All_${smp}.rosv ]
    then
        sed '/#/d' $outputFolder/Complex_2step.rosv | sed '/###/d' | sed '/-Split/d' | sed '/Large-/d' | awk '$16="Two-round-split-alignment"' OFS='\t' >> $outputFolder/Complex_All_${smp}.rosv
    else
        sed '/#/d' $outputFolder/Complex_2step.rosv | sed '/###/d' | sed '/-Split/d' | sed '/Large-/d' | awk '$16="Two-round-split-alignment"' OFS='\t' > $outputFolder/Complex_All_${smp}.rosv
    fi
fi
bash rosv_to_osv.sh $outputFolder/Complex_All_${smp}.rosv tp1_${smp}.osv $outputFolder/translocation_${smp}.rosv
sort -n -k1,1 -k2,2 -k3,3 tp1_${smp}.osv > tp2_${smp}.osv
./Deduplicate_osv_complex -inputSV tp2_${smp}.osv -outputFile compiled_SVs_${smp}.osv &> /dev/null

sed '/#/d' $outputFolder/translocation_${smp}.rosv | cut -f2- | sed 's/\.\.\./\./g' | sed 's/\/Deletion//g' > translocation_${smp}.txt
./Deduplicate_osv_trans -inputSV translocation_${smp}.txt -outputFile $outputFolder/translocation_dedup.osv

rco=5
rm -f indels_${smp}.osv
if [ -f $outputFolder/Indel_All.osv ]
then
    sed '/#/d' $outputFolder/Indel_All.osv | awk '$10="Single-Molecule-Alignment"' OFS='\t' | awk -v var=$rco '$13>=var' > indels_${smp}.osv
fi
if [ -f $outputFolder/Indel_Contig_All.osv ]
then
    if [ -f indels_${smp}.osv ]
    then
        sed '/#/d' $outputFolder/Indel_Contig_All.osv | awk '$10="Contig-Alignment"' OFS='\t' >> indels_${smp}.osv
    else
        sed '/#/d' $outputFolder/Indel_Contig_All.osv | awk '$10="Contig-Alignment"' OFS='\t' > indels_${smp}.osv
    fi
fi

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

# Filter indels with complex SV callset
sed 's/sizeChange=//g' compiled_SVs_${smp}.osv | grep 'Inversion\|Duplication' | awk '$8<500000' | awk '$8="sizeChange="$8' OFS='\t' | sort -nk1,1 -k2,2 -k3,3 > invdup_filter_${smp}.osv
./compBeds_full_new -inputFSFile invdup_filter_${smp}.osv -inputSVFile indels_added_${smp}.osv -outputFile indels_added_val_${smp}.bed > /dev/null
sed '/#/d' indels_added_val_${smp}.bed | cut -f12 > indel_val_${smp}.txt
paste indels_added_${smp}.osv indel_val_${smp}.txt > indels_added_val_${smp}.txt
awk '$(NF)==0' OFS='\t' indels_added_val_${smp}.txt | rev | cut -f2- | rev > indels_val_${smp}.txt
./Deduplicate_osv_indel -inputSV indels_val_${smp}.txt -outputFile indels_${smp}.osv &> /dev/null
./Deduplicate_osv_indel -inputSV indels_val_${smp}.txt -outputFile $outputFolder/indels_dedup.osv -mode 3 &> /dev/null



sed -i 's/\/Deletion//g' $outputFolder/translocation_dedup.osv
./compBeds_3col_trans -inputSVFile $outputFolder/translocation_dedup.osv -inputFSFile density_regions_sorted.bed -outputFile translocation_${smp}_dedup_rem.bed
sed '/compBeds_3col/d' translocation_${smp}_dedup_rem.bed | cut -f5 > val_remove_$smp.txt
paste $outputFolder/translocation_dedup.osv val_remove_$smp.txt | awk 'NR==1 || $NF==0' | cut -f1-10 > $outputFolder/translocation_dedup_filter.osv && mv $outputFolder/translocation_dedup_filter.osv $outputFolder/translocation_dedup.osv
#Use the following command to filter SVs on Y-chromosome
## awk '$1<24 && $3<24 && $NF==0' $outputFolder/translocation_dedup.osv | sed 's/\/Deletion//g' > translocation_${smp}_dedup.osv && mv translocation_${smp}_dedup.osv $outputFolder/translocation_dedup.osv


./compBeds_3col -inputSVFile $outputFolder/indels_dedup.osv -inputFSFile density_regions_sorted.bed -outputFile indels_${smp}_rem.bed
sed '/compBeds_3col/d' indels_${smp}_rem.bed | cut -f4 > val_remove_$smp.txt
paste $outputFolder/indels_dedup.osv val_remove_$smp.txt | awk 'NR==1 || $NF==0' | cut -f1-4,6,8,9,11-14 | sed 's/OMBlast/Single-Molecule-Alignment/g' | sed 's/Insertions/Multi-allelic-double-insertion/g' | sed 's/Deletions/Multi-allelic-double-deletion/g' | sed 's/Both\|Both-ND/Multi-allelic-insertion-deletion/g' > $outputFolder/indels_dedup_filter.osv && mv $outputFolder/indels_dedup_filter.osv $outputFolder/indels_dedup.osv
#Use the following command to filter SVs on Y-chromosome
## awk '$NF==0 && $1<24' $outputFolder/indels_dedup.osv > indels_dedup_${smp}.osv && mv indels_dedup_${smp}.osv $outputFolder/indels_dedup.osv

./compBeds_3col -inputSVFile $outputFolder/compiled_SVs_dedup.osv -inputFSFile density_regions_sorted.bed -outputFile compiled_SVs_${smp}_rem.bed
sed '/compBeds_3col/d' compiled_SVs_${smp}_rem.bed | cut -f4 > val_remove_$smp.txt
paste $outputFolder/compiled_SVs_dedup.osv val_remove_$smp.txt | awk 'NR==1 || $NF==0' | cut -f1-8,13,16 | sed '/Deletion\|Insertion/d' | sed 's/-Split/-realigned/g' > $outputFolder/compiled_SVs_dedup_filter.osv && mv $outputFolder/compiled_SVs_dedup_filter.osv $outputFolder/compiled_SVs_dedup.osv
#Use the following command to filter SVs on Y-chromosome
#awk '$1<24 && $NF==0' $outputFolder/compiled_SVs_dedup.osv > compiled_SVs_dedup_${smp}.osv && mv compiled_SVs_dedup_${smp}.osv $outputFolder/compiled_SVs_dedup.osv


rm -f *${smp}*.osv
rm -f *${smp}*.bed
rm -f *${smp}*.txt
cd $outputFolder
ls -1 | grep -v "dedup" | xargs rm -f
cd ..
