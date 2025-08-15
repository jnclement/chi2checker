#!/bin/bash

IT=0

while [ $IT -lt 25000 ]; do
    UP=$(( $IT + 1000 ))
    root -b -q -l "drawcalo.C(${IT},${UP},0)"
    root -b -q -l "drawf.C(${IT},${UP})"
    #root -b -q -l "drawcalo.C(${IT},${UP},1)"
    IT=$(( $IT + 1000 ))
    echo $IT
done
