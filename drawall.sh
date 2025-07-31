#!/bin/bash

IT=0

while [ $IT -lt 50000 ]; do
    UP=$(( $IT + 100 ))
    root -b -q -l "drawcalo.C(${IT},${UP},0)"
    root -b -q -l "drawcalo.C(${IT},${UP},1)"
    IT=$(( $IT + 100 ))
    echo $IT
done
