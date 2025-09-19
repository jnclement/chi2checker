#!/bin/bash

IT=0
NEVT=0
while [ $IT -lt $(( `cat chi2files.txt | wc -l` + 1 )) ]; do
    UP=$(( $IT + 100 ))
    root -b -q -l "drawcalo.C(${IT},${UP},1)"
    NEVT=$(( $NEVT + $? ))
    #root -b -q -l "drawf.C(${IT},${UP})"
    #root -b -q -l "drawcalo.C(${IT},${UP},1)"
    IT=$(( $IT + 100 ))
    echo $IT
done
echo $NEVT
