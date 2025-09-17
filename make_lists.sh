#!/bin/bash

rm chi2files.txt
rm wavefiles.txt

for rn in `ls ../chi2`; do
    ls ../chi2/${rn}/*testnewprod*chi2* >> chi2files.txt
    ls ../chi2/${rn}/*testnewprod*wave* >> wavefiles.txt
done
