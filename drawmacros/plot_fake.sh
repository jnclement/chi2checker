#!/bin/bash

root -b -q -l "draw_spec_fake.C(0,0,31,40)"

#root -b -q -l "draw_spec_fake.C(0,21,40)"
#root -b -q -l "draw_spec_fake.C(1,21,40)"
#root -b -q -l "draw_spec_fake.C(2,21,40)"
#root -b -q -l "draw_spec_fake.C(0,41,60)"
#root -b -q -l "draw_spec_fake.C(1,41,60)"
#root -b -q -l "draw_spec_fake.C(2,41,60)"
root -b -q -l "draw_spec_fake_multihist.C"
root -b -q -l 'draw_spec_fake_multihist.C("2")'
root -b -q -l "draw_fake.C(0,1)"
