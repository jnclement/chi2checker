#!/bin/bash

root -b -q -l "draw_timingcut.C(0)"
root -b -q -l "draw_timingcut.C(1)"
root -b -q -l "draw_spec.C(31,45)"
root -b -q -l "draw_spec.C(46,55)"
root -b -l -q "draw_dijet_eff.C"
