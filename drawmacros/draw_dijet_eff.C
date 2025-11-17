#include <../dlUtility.h>

int draw_dijet_eff()
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inf = TFile::Open("../../cuteff/final_cuteff.root","READ");

  TH1D* h1_pt = (TH1D*)inf->Get("h3_pt_lem_loh");
  TH1D* h1_pt_dijet = (TH1D*)inf->Get("h3_pt_lem_loh_dijet");

  //TH1D* h1_pt_dijet = h3_pt_dijet->ProjectionX("projx_dijet");
  //TH1D* h1_pt = h3_pt->ProjectionX("projx");

  h1_pt->Rebin(5);
  h1_pt_dijet->Rebin(5);

  TH1D* eff = (TH1D*)h1_pt->Clone("eff");

  eff->Divide(h1_pt_dijet,h1_pt,1,1,"B");

  TCanvas* can = new TCanvas("","",1500,1500);
  ratioPanelCanvas(can,0.3);
  can->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.05);
  gPad->SetLogy();
  TLegend* leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetNColumns(2);

  h1_pt->SetMarkerColor(kBlack);
  h1_pt->SetLineColor(kBlack);
  h1_pt->SetLineWidth(2);
  h1_pt->SetMarkerStyle(20);
  h1_pt->SetMarkerSize(2);
  leg->AddEntry(h1_pt,"No cuts","p");
  h1_pt_dijet->SetMarkerColor(kRed);
  h1_pt_dijet->SetLineColor(kRed);
  h1_pt_dijet->SetLineWidth(2);
  h1_pt_dijet->SetMarkerStyle(71);
  h1_pt_dijet->SetMarkerSize(2);
  leg->AddEntry(h1_pt_dijet,"Dijet cut","p");
  
  eff->SetMarkerColor(kRed);
  eff->SetLineColor(kRed);
  eff->SetLineWidth(2);
  eff->SetMarkerStyle(20);
  eff->SetMarkerSize(2);

  h1_pt->GetYaxis()->SetTitleSize(0.06);
  h1_pt->Draw("PE");
  h1_pt->GetXaxis()->SetRangeUser(15,70);
  h1_pt->GetYaxis()->SetTitle("Scaled Counts");
  h1_pt->GetXaxis()->SetTitle("Uncalibrated p_{T}^{jet}");
  h1_pt_dijet->Draw("SAME PE");
  
  leg->Draw();
  can->cd(2);
  eff->GetYaxis()->SetRangeUser(0,1.199);
  eff->GetXaxis()->SetRangeUser(15,70);
  eff->GetYaxis()->SetTitle("Efficiency");
  eff->GetXaxis()->SetLabelSize(0.1);
  eff->GetYaxis()->SetLabelSize(0.1);
  eff->GetXaxis()->SetTitleSize(0.1);
  eff->GetYaxis()->SetTitleSize(0.1);
  eff->Draw("PE");
  TLine* line = new TLine(15,1,70,1);
  line->SetLineStyle(9);
  line->Draw();
  can->cd(0);
  
  maintexts(0.48,0.15,0,0.03,0,0);
  drawText("Truth-reco matched jets",0.15,0.375,0,kBlack,0.03);
  drawText("No reconstructed z_{vtx} requirement",0.15,0.34,0,kBlack,0.03);

  can->SaveAs("../../images/dnp/sim_dijet_eff.pdf");

  return 0;
}
