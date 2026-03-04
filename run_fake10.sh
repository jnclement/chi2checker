#!/bin/bash


#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_jet10.root",true,false,"jet10")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_jet30.root",true,false,"jet30")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_jet50.root",true,false,"jet50")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_dat.root",false,false,"dat")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_dat.root",false,false,"dat",true)'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_dat.root",false,false,"dat",true,true)'


root -l -b -q 'fake10.C("./savedhistsdat.list",false,false,"dat",true,false,"20250220")'
root -l -b -q 'fake10.C("./savedhistsdat.list",false,true,"dat",true,false,"20250220")'
root -l -b -q 'fake10.C("./savedhistsjet30.list",true,true,"jet30",true,false,"20250220")'
root -l -b -q 'fake10.C("./savedhistsjet30.list",true,true,"jet30",false,false,"20250220")'
root -l -b -q 'fake10.C("./savedhistsjet30.list",true,false,"jet30",true,false,"20250220")'
root -l -b -q 'fake10.C("./savedhistsjet30.list",true,false,"jet30",false,false,"20250220")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20250208_dat.root",false,false,"dat",true,false,"20250208")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251211_dat.root",false,false,"dat",true,true,"20251211")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251211_dat.root",false,false,"dat",true,false,"20251211")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_dat.root",false,false,"dat",true,true)'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_jet10.root",true,true,"jet10")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_jet30.root",true,true,"jet30")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_jet50.root",true,true,"jet50")'
#root -l -b -q 'fake10.C("../savedhists/savedhists_20251108_dat.root",false,true,"dat",true)'
