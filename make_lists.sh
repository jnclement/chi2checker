#!/bin/bash

find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260216*chi2file*" | sort -V > chi2filesdat.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet5_*chi2file*" | sort -V > chi2filesjet5.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet12*chi2file*" | sort -V > chi2filesjet12.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet20*chi2file*" | sort -V > chi2filesjet20.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet30*chi2file*" | sort -V > chi2filesjet30.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet40*chi2file*" | sort -V > chi2filesjet40.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet50*chi2file*" | sort -V > chi2filesjet50.txt
find /sphenix/tg/tg01/jets/jocl/chi2/ -type f -name "*20260219_jet60*chi2file*" | sort -V > chi2filesjet60.txt
