#!/bin/bash
####  Created by Le Li, May 11, 2023
####  Script of running complex SV caller of COMSV pipeline

        smp="261T"
        conf=10
        minSegLen=10000
        minSites=3
        outputFolder="COMSV_SVs_$smp"
        inputAlign="${smp}_chr1.oma"
        ref="${smp}_r_chr1.cmap"
        inputAlignContig="${smp}_contig_chr1.oma"
        refContig="${smp}_contig_r_chr1.cmap"
        inputAlign2round="${smp}_2step_chr1.oma"
        ref2round="${smp}_2step_r_chr1.cmap"

        ./COMSV_complex -SVoutputFile_CR $outputFolder/Complex_contig.rosv -minS 1 -optCRAlignFile $inputAlignContig -highDensityFile density_PAR_region_DLE1.bed -minScore $conf -refMapFile $refContig -minSegLen $minSegLen -minSites $minSites
        ./COMSV_complex -SVoutputFile_CR $outputFolder/Complex.rosv -minS 2 -optCRAlignFile $inputAlign -highDensityFile density_PAR_region_DLE1.bed -minScore $conf -refMapFile $ref -minSegLen $minSegLen -minSites $minSites
        ./COMSV_complex -SVoutputFile_CR $outputFolder/Complex_2step.rosv -minS 3 -optCRAlignFile $inputAlign2round -highDensityFile density_PAR_region_DLE1.bed -minScore $conf -refMapFile $ref2round -minSegLen $minSegLen -minSites $minSites
