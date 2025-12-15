#!/bin/bash

rm chi2filesdat.txt
#rm chi2filesjet30.txt
#rm chi2filesjet50.txt
#rm chi2filesjet10.txt
#rm wavefiles.txt
#rm chi2files50.txt
#rm chi2files70.txt
#rm chi2files30.txt
#rm chi2files20.txt
#rm chi2files10.txt
#rm chi2filesbak.txt
for rn in `ls /sphenix/tg/tg01/jets/jocl/chi2`; do
    echo $rn
    ls /sphenix/tg/tg01/jets/jocl/chi2/${rn}/*20251211_dat*chi2file* >> chi2filesdat.txt
    #ls /sphenix/tg/tg01/jets/jocl/chi2/${rn}/*20251106jet10*chi2file* >> chi2filesjet10.txt
    #ls /sphenix/tg/tg01/jets/jocl/chi2/${rn}/*20251106jet30*chi2file* >> chi2filesjet30.txt
    #ls /sphenix/tg/tg01/jets/jocl/chi2/${rn}/*20251106jet50*chi2file* >> chi2filesjet50.txt
done
#for rn in {0..9999}; do #`ls /sphenix/tg/tg01/jets/jocl/chi2`; do
    #ls ../chi2/${rn}/*20250925dat*chi2* >> chi2filesdat.txt
    #ls ../chi2/${rn}/*20250925dat*wave* >> wavefiles.txt
#    if [ $(( $rn % 100 )) -eq 0 ]; then
#	echo $rn
#    fi
#    ls ../chi2/${rn}/*20250925sim10*chi2* >> chi2files10.txt
#    ls ../chi2/${rn}/*20250925sim20*chi2* >> chi2files20.txt
#    ls ../chi2/${rn}/*20250925sim30*chi2* >> chi2files30.txt
#    ls ../chi2/${rn}/*20250925sim50*chi2* >> chi2files50.txt
#    ls ../chi2/${rn}/*20250925sim70*chi2* >> chi2files70.txt
#done
#rm chi2filescosmic.txt
#for rn in {0..499}; do
#    ls ../chi2/${rn}/*20250924cosmic1000to2000* >> chi2filescosmic.txt
#done
#for rn in `ls ../chi2_bak`; do
#    ls ../chi2_bak/$rn/*20250918*chi2* >> chi2filesbak.txt
#done
