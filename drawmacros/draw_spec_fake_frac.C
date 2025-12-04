#include <../dlUtility.h>
int draw_spec_fake_frac(string tt = "dt", string ap = "yz", int lo = 30, int hi = 45)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* inf = TFile::Open("../hists_mbdtimereq_out_datsam.root","READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get(("hpt"+tt+"fracdat").c_str());

  TCanvas* can = new TCanvas("","",1500,1500);
  
  can->cd();
  gPad->SetTopMargin(0.15);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);

  //gPad->SetLogz();
  TBox* box = new TBox(-3,-3,3,3);
  box->SetFillStyle(0);
  box->SetLineColor(kBlue);
  box->SetLineWidth(3);
  if(ap=="xy") h3_pt_lem_loh->GetZaxis()->SetRange(lo,hi);
  else if(ap=="xz") h3_pt_lem_loh->GetYaxis()->SetRange(lo,hi);
  else if(ap=="yz") h3_pt_lem_loh->GetXaxis()->SetRange(lo,hi);
  TH2D* h2_t_dt = (TH2D*)(h3_pt_lem_loh->Project3D(ap.c_str()));
  //h2_t_dt->GetXaxis()->SetTitle("MBD - t_{lead} [ns]");
  //h2_t_dt->GetZaxis()->SetTitle("N_{jet}");
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  /*
  h2_t_dt->Rebin2D(4,4);
  h2_t_dt->GetXaxis()->SetRangeUser(-50,50);
  h2_t_dt->GetYaxis()->SetRangeUser(-25,25);
  h2_t_dt->GetYaxis()->SetTitleSize(0.05);
  h2_t_dt->GetXaxis()->SetTitleSize(0.045);
  */
  h2_t_dt->Draw("COLZ");
  //box->Draw();
  maintexts(0.7,0.25,0,0.04,1,0);
  stringstream tdss;
  string axstr;
  float lobound, hibound;
  string unitstr;
  if(ap=="yz")
    {
      axstr = "<p_{T}^{uncalib}<";
      lobound = h3_pt_lem_loh->GetXaxis()->GetBinLowEdge(lo);
      hibound = h3_pt_lem_loh->GetXaxis()->GetBinLowEdge(hi+1);
      unitstr = " GeV";
    }
  else if(ap=="xz")
    {
      axstr = tt=="dt"?"<#Delta-t<":"<t_{lead}<";
      lobound = h3_pt_lem_loh->GetYaxis()->GetBinLowEdge(lo);
      hibound = h3_pt_lem_loh->GetYaxis()->GetBinLowEdge(hi+1);
      unitstr = " ns";
    }
  else if(ap=="xy")
    {
      axstr = "<E^{jet}_{OHCal}/E^{jet}<";
      lobound = h3_pt_lem_loh->GetZaxis()->GetBinLowEdge(lo);
      hibound = h3_pt_lem_loh->GetZaxis()->GetBinLowEdge(hi+1);
      unitstr = "";
    }
  tdss << std::fixed << std::setprecision(1) << "Leading jets " << lobound << axstr << hibound << unitstr;
  drawText(tdss.str().c_str(),0.25,0.75,0,kBlack,0.04);
  drawText("No reconstructed z_{vtx} requirement",0.25,0.55,0,kBlack,0.04);
  can->SaveAs(("../../images/efd/data_"+tt+"_"+ap+"_proj2d_"+to_string(lo)+"-"+to_string(hi)+".pdf").c_str());
  gPad->SetLogz();
  can->SaveAs(("../../images/efd/data_"+tt+"_"+ap+"_proj2d_"+to_string(lo)+"-"+to_string(hi)+"_log.pdf").c_str());
  return 0;
}
