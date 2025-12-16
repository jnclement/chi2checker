#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam.root", "hptdtemfracdat", {{21,100}}, {{71,130}}, {2}, {kAzure}, {20}, {"#splitline{|#Delta-t|<3 ns && -8 ns<t_{lead}<4 ns}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "emfrac_timingcuts.pdf", "Normalized Counts")'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam.root", "hptdtemfracnotdat", {{21,100}}, {{91,120}}, {1}, {kAzure}, {20}, {"#splitline{0.9<E^{EMCal}/E_{jet}^{uncalib}}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "dt_90emfrac.pdf", "Normalized Counts")'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam_slewed.root", "hptdtemfracnotdat", {{16,30}}, {{96,120}}, {1}, {kAzure}, {20}, {"#splitline{0.95<E^{EMCal}/E_{jet}^{uncalib}}{&& 15 GeV<p_{T}^{uncalib}<30 GeV}"}, "dt_95emfrac.pdf", "Normalized Counts")'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_bdatsam.root", "hptdtemfracnotbdat", {{16,30}}, {{96,120}}, {1}, {kAzure}, {20}, {"#splitline{0.95<E^{EMCal}/E_{jet}^{uncalib}}{&& 15 GeV<p_{T}^{uncalib}<30 GeV}"}, "dt_95emfrac.pdf", "Normalized Counts")'

root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam20251211.root", "hpttemfracdat", {{16,30},{16,30},{16,30},{16,30}}, {{1,20},{41,60},{61,80},{101,120}}, {1,1,1,1}, {kAzure,kGreen,kRed,kViolet}, {20,21,71,72}, {"E^{EMCal}/E_{jet}^{uncalib}<0.1","0.3<E^{EMCal}/E_{jet}^{uncalib}<0.5","0.5<E^{EMCal}/E_{jet}^{uncalib}<0.7","E^{EMCal}/E_{jet}^{uncalib}>0.9"}, "../../images/efd/t_emfrac_tem.pdf", "Normalized Counts")'

root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam20251211.root", "hpttohfracdat", {{16,30},{16,30},{16,30},{16,30}}, {{1,20},{41,60},{61,80},{101,120}}, {1,1,1,1}, {kAzure,kGreen,kRed,kViolet}, {20,21,71,72}, {"E^{EMCal}/E_{jet}^{uncalib}<0.1","0.3<E^{EMCal}/E_{jet}^{uncalib}<0.5","0.5<E^{EMCal}/E_{jet}^{uncalib}<0.7","E^{EMCal}/E_{jet}^{uncalib}>0.9"}, "../../images/efd/t_emfrac_toh.pdf", "Normalized Counts")'

root -b -q -l 'draw_spec_fake_calotime.C("xy",1,20)'
root -b -q -l 'draw_spec_fake_calotime.C("xy",41,60)'
root -b -q -l 'draw_spec_fake_calotime.C("xy",61,80)'
root -b -q -l 'draw_spec_fake_calotime.C("xy",101,120)'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam.root", "hptdtemfracnotdat", {{21,100}}, {{100,120}}, {1}, {kAzure}, {20}, {"#splitline{0.99<E^{EMCal}/E_{jet}^{uncalib}}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "dt_99emfrac.pdf", "Normalized Counts")'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_datsam.root", "hptdtemfracnotdat", {{21,100}}, {{1,90}}, {1}, {kAzure}, {20}, {"#splitline{0.9>E^{EMCal}/E_{jet}^{uncalib}}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "dt_l90emfrac.pdf", "Normalized Counts")'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_jet30.root", "hptdtemfracnotjet30", {{21,100}}, {{91,120}}, {1}, {kAzure}, {20}, {"#splitline{0.9<E^{EMCal}/E_{jet}^{uncalib}}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "j30dt_90emfrac.pdf", "Normalized Counts",0)'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_jet10.root", "hptdtemfracnotjet10", {{16,30}}, {{96,120}}, {1}, {kAzure}, {20}, {"#splitline{0.95<E^{EMCal}/E_{jet}^{uncalib}}{&& 15 GeV<p_{T}^{uncalib}<30 GeV}"}, "j10dt_95emfrac.pdf", "Normalized Counts",0)'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_jet30.root", "hptdtemfracnotjet30", {{21,100}}, {{100,120}}, {1}, {kAzure}, {20}, {"#splitline{0.99<E^{EMCal}/E_{jet}^{uncalib}}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "j30dt_99emfrac.pdf", "Normalized Counts",0)'

#root -l -b -q 'drawprettyeff.C("../hists_mbdtimereq_out_jet30.root", "hptdtemfracnotjet30", {{21,100}}, {{1,90}}, {1}, {kAzure}, {20}, {"#splitline{0.9>E^{EMCal}/E_{jet}^{uncalib}}{&& 20 GeV<p_{T}^{uncalib}<100 GeV}"}, "j30dt_l90emfrac.pdf", "Normalized Counts",0)'
