#!/bin/bash
####  Created by Le Li, May 11, 2023
####  Script of running indel caller of COMSV pipeline

    if [ "$1" == "help" ] || [ "$1" == "HELP" ] || [ "$1" == "--help" ] || [ "$1" == "--HELP" ]
    then
        ./COMSV_indel -help
        exit
    fi

    # the input data for the sample (sample name, label, alignment file, reference file)
    smp="261T"
    label="101"
    inputAlign="${smp}_chr1.oma"
    ref="${smp}_r_chr1.cmap"
    sed -i '/Unmap/d' $inputAlign # make sure Unmap entries are not included


    # default settings
    siz=2000
    minInsClusterSize=3
    minDelClusterSize=3
    numberOfSupport=20
    tpFolder="temp/"
    outputFolder="COMSV_SVs_$smp"
    reapF="density_PAR_region_DLE1.bed"
    mkdir -p $outputFolder
    rm -f $outputFolder/*

    # running initialization with chromosomeID=-1
    ./COMSV_indel -inputRepeatFileName $reapF -inputLabel  $label  -outputFolder  $outputFolder  -chrMapFile $ref -optTempFolder $tpFolder -minIndelSize $siz -chromosomeID -1 -optAlignFile $inputAlign
    # running COMSV in parallel for all chromosomes (please change the numbers 1-24 based on your own chromosome IDs)
#    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
    for chr in 1
    do
        if [ $chr -eq 1 ]
        then
            ./COMSV_indel -inputRepeatFileName $reapF -inputLabel  $label  -outputFolder  $outputFolder  -chrMapFile $ref -optTempFolder $tpFolder -minInsClusterSize $minInsClusterSize -minDelClusterSize $minDelClusterSize -numberOfSupportIndelMolecule $numberOfSupport -minIndelSize $siz -chromosomeID $chr -optAlignFile $inputAlign && pkill sleep && echo "The longest chromosome has been finished, wake up!" &
        else
            ./COMSV_indel -inputRepeatFileName $reapF -inputLabel  $label  -outputFolder  $outputFolder  -chrMapFile $ref -optTempFolder $tpFolder -minInsClusterSize $minInsClusterSize -minDelClusterSize $minDelClusterSize -numberOfSupportIndelMolecule $numberOfSupport -minIndelSize $siz -chromosomeID $chr -optAlignFile $inputAlign || echo "error in $chr" &
        fi
    done

    stp=0
    while [ $stp -eq 0 ]
    do
        if ps ax | grep -v grep | grep "optAlignFile $inputAlign" > /dev/null
        then
            echo "Not finished yet, snooze two minutes (Zzz..)"
            sleep 5m
        else
            stp=1
        fi
    done

    # integrate the SVs on all chromosomes
#    for chr in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
    # call SVs on chromosome 1
    for chr in 1
    do
        if [ ! -f "$outputFolder/Indel_All.osv" ]
        then
            cp $outputFolder/Label${label}_Chr${chr}_${siz}.osv $outputFolder/Indel_All.osv
        else
            sed '/#/d' $outputFolder/Label${label}_Chr${chr}_${siz}.osv >> $outputFolder/Indel_All.osv
        fi
    done
