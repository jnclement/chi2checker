#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n new
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
export TESTINSTALL=/sphenix/user/jocl/projects/testinstall
if [[ ! -z "$_CONDOR_SCRATCH_DIR" && -d $_CONDOR_SCRATCH_DIR ]]; then
    cd $_CONDOR_SCRATCH_DIR
else
    echo condor scratch NOT set
    exit -1
fi

cp /sphenix/user/jocl/projects/chi2checker/src/drawcalo.C .

START=$2
END=$3

mkdir -p inputfiles/
mkdir -p savedhists/
sed -n "$((START+1)),$((END))p" /sphenix/user/jocl/projects/chi2checker/src/chi2files${1}.txt > thelist.list

while IFS= read -r infile; do
    # Skip empty lines
    [[ -z "$infile" ]] && continue

    # Extract base filename
    fname=$(basename "$infile")

    # Copy with dd
    echo "Copying $infile to $fname"
    dd if="$infile" of="inputfiles/$fname" bs=64M
done < thelist.list

ls -v inputfiles/* > inroots.list

mkdir images/

#root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesjet10.txt",0,-1,-1,0,0)'
#root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesjet30.txt",0,-1,-1,0,0)'
#root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesjet50.txt",0,-1,-1,0,0)'
root -b -l -q 'drawcalo.C('${2}','${3}',1,"'${1}'")'
cp savedhists/* /sphenix/user/jocl/projects/chi2checker/savedhists/
cp images/* /sphenix/user/jocl/projects/chi2checker/images/disp/
