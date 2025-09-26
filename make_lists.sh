#!/bin/bash

#rm chi2filesdat.txt
#rm wavefiles.txt
rm chi2files50.txt
rm chi2files70.txt
rm chi2files30.txt
#rm chi2filesbak.txt
for rn in `ls ../chi2`; do
    #ls ../chi2/${rn}/*20250919dat*chi2* >> chi2filesdat.txt
    #ls ../chi2/${rn}/*20250919dat*wave* >> wavefiles.txt
    ls ../chi2/${rn}/*20250925sim50*chi2* >> chi2files50.txt
    ls ../chi2/${rn}/*20250925sim70*chi2* >> chi2files70.txt
    ls ../chi2/${rn}/*20250925sim30*chi2* >> chi2files30.txt
done

#for rn in `ls ../chi2_bak`; do
#    ls ../chi2_bak/$rn/*20250918*chi2* >> chi2filesbak.txt
#done
