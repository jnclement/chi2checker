#!/bin/bash

IT=0
NEVT=0

while [ $IT -lt $(( `cat chi2filesdat.txt | wc -l` + 1 )) ]; do
    UP=$(( $IT + 250 ))
    root -b -q -l "drawcalo.C(${IT},${UP},1)"
    NEVT=$(( $NEVT + $? ))
    #root -b -q -l "timing_hists.C(${IT},${UP},0)"
    #root -b -q -l "drawf.C(${IT},${UP})"
    #root -b -q -l "drawcalo.C(${IT},${UP},1)"
    IT=$(( $IT + 250 ))
    echo $IT
done
echo $NEVT
IT=0
#while [ $IT -lt $(( `cat chi2files50.txt | wc -l` + 1 )) ]; do
#    UP=$(( $IT + 10 ))
#    root -b -q -l "timing_hists.C(${IT},${UP},1)"
#    root -b -q -l "timing_hists.C(${IT},${UP},2)"
#    IT=$(( $IT + 10 ))
#done


