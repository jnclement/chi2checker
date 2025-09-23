#!/bin/bash

rm chi2filesdat.txt
rm wavefiles.txt
rm chi2files50.txt
rm chi2files70.txt

for rn in `ls ../chi2`; do
    ls ../chi2/${rn}/*20250919dat*chi2* >> chi2filesdat.txt
    ls ../chi2/${rn}/*20250919dat*wave* >> wavefiles.txt
    ls ../chi2/${rn}/*20250919sim50*chi2* >> chi2files50.txt
    ls ../chi2/${rn}/*20250919sim70*chi2* >> chi2files70.txt
done
