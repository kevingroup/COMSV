
smp="261T"
outputFolder="COMSV_complex_$smp"

sed '/#/d' $outputFolder/Label301_${conf}.rosv | sed '/###/d' | awk '$16="Contig-Alignment"' OFS='\t' > $outputFolder/Complex_All_${smp}.rosv
sed '/#/d' $outputFolder/Label401_${conf}.rosv | sed '/###/d' | sed '/Duplication/!d' | awk '$16="Single-Molecule-Alignment"' OFS='\t' >> $outputFolder/Complex_All_${smp}.rosv
sed '/#/d' $outputFolder/Label501_${conf}.rosv | sed '/###/d' | sed '/-Split/d' | sed '/Large-/d' | awk '$16="Two-round-split-alignment"' OFS='\t' >> $outputFolder/Complex_All_${smp}.rosv
bash rosv_to_osv.sh $outputFolder/Complex_All_${smp}.rosv $outputFolder/tp1_${smp}.osv $outputFolder/translocation_${smp}.rosv
sort -n -k1,1 -k2,2 -k3,3 $outputFolder/tp1_${smp}.osv > $outputFolder/tp2_${smp}.osv
./Deduplicate_osv_complex -inputSV $outputFolder/tp2_${smp}.osv -outputFile $outputFolder/compiled_SVs_${smp}.osv &> /dev/null

sed '/#/d' $outputFolder/translocation_${smp}.rosv | cut -f2- | sed 's/\.\.\./\./g' | sed 's/\/Deletion//g' > $outputFolder/translocation_${smp}.txt
./Deduplicate_osv_trans -inputSV $outputFolder/translocation_${smp}.txt -outputFile $outputFolder/translocation_${smp}_dedup.osv
rm $outputFolder/translocation_${smp}_dedup_unfilter.osv

