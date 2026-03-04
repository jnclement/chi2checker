#!/bin/bash



#while [ $IT -lt $(( `cat chi2filescosmic.txt | wc -l` + 1 )) ]; do
#    UP=$(( $IT + 10 ))
#    root -b -q -l 'drawcalo.C('${IT}','${UP}',0,"chi2filescosmic.txt",0,-1,-1,0,1)'
#    IT=$(( $IT + 10 ))
#done
rm ranges.txt
for sname in jet30; do
    IT=0
    NEVT=0
    MAX=$(( `cat chi2files${sname}.txt | wc -l` + 1 ))
    while [ $IT -lt $MAX ]; do
	UP=$(( $IT + 20 ))
	echo "${sname} ${IT} ${UP}" >> ranges.txt
	IT=$(( $IT + 20 ))
    done
done
condor_submit condraw.job
