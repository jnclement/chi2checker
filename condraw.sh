#!/bin/bash

source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.515
source /opt/sphenix/core/bin/setup_local.sh "/sphenix/user/jocl/projects/testinstall"
export HOME=/sphenix/u/jocl
export TESTINSTALL=/sphenix/user/jocl/projects/testinstall

#root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesjet10.txt",0,-1,-1,0,0)'
#root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesjet30.txt",0,-1,-1,0,0)'
#root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesjet50.txt",0,-1,-1,0,0)'
root -b -l -q 'drawcalo.C('${1}','${2}',1,"chi2filesblair.txt")'
