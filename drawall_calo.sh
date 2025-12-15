#!/bin/bash

IT=0
NEVT=0


#while [ $IT -lt $(( `cat chi2filescosmic.txt | wc -l` + 1 )) ]; do
#    UP=$(( $IT + 10 ))
#    root -b -q -l 'drawcalo.C('${IT}','${UP}',0,"chi2filescosmic.txt",0,-1,-1,0,1)'
#    IT=$(( $IT + 10 ))
#done

rm ranges.txt
MAX=$(( `cat chi2filesblair.txt | wc -l` + 1 ))
while [ $IT -lt $MAX ]; do
    UP=$(( $IT + 10 ))
    #root -b -q -l "drawcalo.C(${IT},${UP},1)"
    #NEVT=$(( $NEVT + $? ))
    
    #root -b -q -l "compare_events.C(${IT},${UP})" 
#    root -b -q -l "timing_hists.C(${IT},${UP},0)"
    #root -b -q -l "drawf.C(${IT},${UP})"
    echo "${IT} ${UP}" >> ranges.txt
    IT=$(( $IT + 10 ))

done

condor_submit condraw.job
#IT=0
#while [ $IT -lt $(( `cat chi2files50.txt | wc -l` + 1 )) ]; do
#    UP=$(( $IT + 100 ))
#    root -b -q -l "timing_hists.C(${IT},${UP},1)"
#    root -b -q -l "timing_hists.C(${IT},${UP},3)"
#    root -b -q -l "timing_hists.C(${IT},${UP},2)"
#
#    root -b -q -l "cuteff.C(${IT},${UP},1)"
#    root -b -q -l "cuteff.C(${IT},${UP},3)"
#    root -b -q -l "cuteff.C(${IT},${UP},2)"
#    IT=$(( $IT + 100 ))
#done


