for sp in -8 29 -1 10 3 4 -4 5 -3 7 6 -5 2 -7 -6 9 -2 8 1 20 21 19 17 15 11 18 13 14 16 12 24 23 22
#for sp in -1 10 -4 5 -3 7 -5 2 -2 8
#for sp in 27
#for((sp=-8;sp<=26;sp+=1))
do
if [ $sp -eq 0 ]
then
    continue
fi
if [ $sp -eq -8 ]
then
       smp1="A549_BSPQI"
       smp="A549"
       CTP="LUNG_BRONCHUS"
       STP="A549"
elif [ $sp -eq -7 ]
then
       smp1="PANC-1"
       smp="PANC-1"
       CTP="PANCREAS"
       STP="PANC-1"
elif [ $sp -eq -6 ]
then
       smp1="SK-N-MC"
       smp="SK-N-MC"
       CTP="BRAIN"
       STP="SK-N-MC"
elif [ $sp -eq -5 ]
then
       smp1="H460_BspQI"
       smp="NCI-H460 Nt.BspQI"
       CTP="LUNG_BRONCHUS"
       STP="NCI-H460"
elif [ $sp -eq -4 ]
then
       smp1="K562_BspQI"
       smp="K562 Nt.BspQI"
       CTP="BLOOD_BONE_MARROW_HEMATOPOIETIC_SYS"
       STP="K562"
elif [ $sp -eq -3 ]
then
       smp1="LNCaP_BspQI"
       smp="LNCaP Nt.BspQI"
       CTP="PROSTATE_GLAND"
       STP="LNCAP"
elif [ $sp -eq -2 ]
then
       smp1="T47D_BspQI"
       smp="T47D Nt.BspQI"
       CTP="BREAST"
       STP="T47D"
elif [ $sp -eq -1 ]
then
       smp1="Caki2_BspQI"
       smp="Caki2 Nt.BspQI"
       CTP="KIDNEY"
       STP="Caki2"
elif [ $sp -eq 1 ]
then
       smp="261T"
       smp1="261T"
       CTP="LIVER"
       STP="-"
#       rm validation_databases/validation_results_all_novel.txt
elif [ $sp -eq 2 ]
then
       smp="NCI-H460 DLE-1"
       smp1="H460_DLE-1"
       CTP="LUNG_BRONCHUS"
       STP="NCI-H460"
elif [ $sp -eq 3 ]
then
#       continue
       smp="HCC827_DLE1"
       CTP="-"
       STP="-"
elif [ $sp -eq 4 ]
then
       smp="Huh7"
       smp1="Huh7"
       CTP="LIVER"
       STP="-"
elif [ $sp -eq 5 ]
then
       smp="K562 DLE-1"
       smp1="K562_DLE-1"
       CTP="BLOOD_BONE_MARROW_HEMATOPOIETIC_SYS"
       STP="K562"
elif [ $sp -eq 6 ]
then
#       continue
       smp="LO2"
       smp1="LO2"
       CTP="LIVER"
       STP="-"
elif [ $sp -eq 7 ]
then
       smp="LNCaP DLE-1"
       smp1="LNCaP_DLE-1"
       CTP="PROSTATE_GLAND"
       STP="LNCAP"
elif [ $sp -eq 8 ]
then
       smp="T47D DLE-1"
       smp1="T47D_DLE-1"
       CTP="BREAST"
       STP="T47D"
elif [ $sp -eq 9 ]
then
       smp="SK-BR-3"
       smp1="SK-BR-3"
       CTP="BREAST"
       STP="SKBR3"
elif [ $sp -eq 10 ]
then
       smp="Caki2 DLE-1"
       smp1="Caki2_DLE-1"
       CTP="KIDNEY"
       STP="Caki2"
elif [ $sp -eq 11 ]
then
#       continue
       smp="7403 Non-cancer"
       smp1="7403_Non-cancer"
       CTP="BLOOD_BONE_MARROW_HEMATOPOIETIC_SYS"
       STP="-"
elif [ $sp -eq 12 ]
then
#       continue
       smp="7622 Non-cancer"
       smp1="7622_Non-cancer"
       CTP="BLOOD_BONE_MARROW_HEMATOPOIETIC_SYS"
       STP="-"
elif [ $sp -eq 13 ]
then
       smp="7528 Cancer"
       smp1="7528_Cancer"
       CTP="MOUTH"
       STP="-"
elif [ $sp -eq 14 ]
then
#       continue
       smp="7528 Non-cancer"
       smp1="7528_Non-cancer"
       CTP="BLOOD_BONE_MARROW_HEMATOPOIETIC_SYS"
       STP="-"
elif [ $sp -eq 15 ]
then
       smp="7403 Cancer"
       smp1="7403_Cancer"
       CTP="MOUTH"
       STP="-"
elif [ $sp -eq 16 ]
then
       smp1="7622_Cancer"
       smp="7622 Cancer"
       CTP="MOUTH"
       STP="-"
elif [ $sp -eq 17 ]
then
#       continue
       smp1="7052_Non-cancer"
       smp="7052 Non-cancer"
       CTP="BLOOD_BONE_MARROW_HEMATOPOIETIC_SYS"
       STP="-"
elif [ $sp -eq 18 ]
then
       smp="7518"
       smp1="7518"
       CTP="THYROID_GLAND"
       STP="-"
elif [ $sp -eq 19 ]
then
       smp1="7052_Cancer"
       smp="7052 Cancer"
       CTP="MOUTH"
       STP="-"
elif [ $sp -eq 20 ]
then
       smp1="3096"
       smp="3096"
       CTP="URINARY_BLADDER"
       STP="-"
elif [ $sp -eq 21 ]
then
       smp1="3717"
       smp="3717"
       CTP="THYROID_GLAND"
       STP="-"
elif [ $sp -eq 22 ]
then
       smp1="14369"
       smp="14369"
       CTP="LUNG_BRONCHUS"
       STP="-"
elif [ $sp -eq 23 ]
then
       smp1="10974"
       smp="10974"
       CTP="LIVER"
       STP="-"
elif [ $sp -eq 24 ]
then
       smp1="7708"
       smp="7708"
       CTP="THYROID_GLAND"
       STP="-"
elif [ $sp -eq 25 ]
then
       bash generate_shuffle_indels.sh
       smp="Shuffle DLE-1"
       smp1="Shuffle_DLE-1"
       CTP="-"
       STP="-"
elif [ $sp -eq 26 ]
then
       bash generate_shuffle_indels_BSPQI.sh
       smp="Shuffle Nt.BspQI"
       smp1="Shuffle_BSPQI"
       CTP="-"
       STP="-"
elif [ $sp -eq 27 ]
then
       smp="Total DLE-1"
       smp1="Total_DLE-1"
       CTP="-"
       STP="-"
elif [ $sp -eq 28 ]
then
       smp="Total Nt.BspQI"
       smp1="Total_BSPQI"
       CTP="-"
       STP="-"
elif [ $sp -eq 29 ]
then
       smp1="C666-1"
       smp="C666-1"
       CTP="ESOPHAGUS"
       STP="C666-1"
fi

    if ! [ -f indels_dedup_${sp}_unfilter.osv ]
    then
        cp indels_dedup_${sp}.osv indels_dedup_${sp}_unfilter.osv
        ./compBeds_3col -inputSVFile indels_dedup_${sp}.osv -inputFSFile density_regions_sorted.bed -outputFile indels_${sp}_rem.bed
        sed '/compBeds_3col/d' indels_${sp}_rem.bed | cut -f4 > val_remove.txt
        paste indels_dedup_${sp}.osv val_remove.txt > indels_dedup_${sp}_filter.osv
        awk '$NF==0 && $1<24' indels_dedup_${sp}_filter.osv > indels_dedup_${sp}.osv
    fi
done
