#include <../dlUtility.h>
int isdat = 0;
string stype = "jet10";
int draw_spec_fake_frac(string tt = "dt", string ap = "yz", int lo = 30, int hi = 45)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* inf = TFile::Open(("../hists_mbdtimereq_out_"+stype+(isdat?"sam":"")+".root").c_str(),"READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get(("hpt"+tt+"frac"+stype).c_str());

  TCanvas* can = new TCanvas("","",1500,1500);
  
  can->cd();
  gPad->SetTopMargin(0.05);
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
  maintexts(0.89,0.25,0,0.04,isdat,0);
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
  TLine* line0 = new TLine(-0.1,-3,1.1,-3);
  TLine* line1 = new TLine(-0.1,3,1.1,3);
  TLine* line2 = new TLine(-0.1,-4,1.1,-4);
  TLine* line3 = new TLine(-0.1,4,1.1,4);
  TProfile* prof = h2_t_dt->ProfileX();
  TLegend* tleg = new TLegend(0.4,0.12,0.9,0.25);
  tleg->SetFillStyle(0);
  tleg->SetBorderSize(0);
  tleg->SetFillColor(0);
  if(ap=="yz")
    {
      line0->SetLineColor(kRed);
      line2->SetLineColor(kAzure+2);
      line1->SetLineColor(kRed);
      line3->SetLineColor(kAzure+2);
      line0->SetLineWidth(2);
      line1->SetLineWidth(2);
      line2->SetLineWidth(2);
      line3->SetLineWidth(2);
      tleg->AddEntry(line0,"Nominal cut bounds","l");
      tleg->AddEntry(line2,"Cut variation bounds","l");
      line0->Draw();
      line1->Draw();
      line2->Draw();
      line3->Draw();
      tleg->Draw();
      prof->SetMarkerStyle(20);
      prof->SetMarkerSize(2);
      prof->SetMarkerColor(kGreen+2);
      prof->SetLineColor(kGreen+2);
      prof->Draw("SAME PE");
    }
      
  tdss << std::fixed << std::setprecision(1) << "Leading jets " << lobound << axstr << hibound << unitstr;
  drawText(tdss.str().c_str(),0.25,0.76,0,kBlack,0.04);
  drawText("No reconstructed z_{vtx} requirement",0.25,0.71,0,kBlack,0.04);
  if(!isdat) drawText(("PYTHIA "+stype+" sample").c_str(),0.25,0.66,0,kBlack,0.04);
  can->SaveAs(("../../images/efd/"+stype+"_"+tt+"_"+ap+"_proj2d_"+to_string(lo)+"-"+to_string(hi)+".pdf").c_str());
  gPad->SetLogz();
  can->SaveAs(("../../images/efd/"+stype+"_"+tt+"_"+ap+"_proj2d_"+to_string(lo)+"-"+to_string(hi)+"_log.pdf").c_str());
  return 0;
}
