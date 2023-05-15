##convert rosv to osv##

sed '/#/d' $1 | sed '/-TL/!d' | sort -n -k2,2 -k4,4 -k3,3 -k5,5 > $3
sed '/#/d' $1 | sed '/-TL/d' | sed 's/Invert-//g' | cut -f2,3,5-12,16 | awk '$2>$3{temp=$3;$3=$2;$2=temp}1' OFS='\t' | awk '$7="["$7"]-["$8"]"' OFS='\t' | awk '{if($4~/^Inversion/ || $9==0) $8="sizeChange="$3-$2; else $8="sizeChange="sqrt($9^2);}1' OFS='\t' | awk '$9="None\t"$11"\t-\t-"' OFS='\t' | sort -n -k1,1 -k2,2 -k3,3 | sed 's/Large-//g' | cut -f-13 > $2
