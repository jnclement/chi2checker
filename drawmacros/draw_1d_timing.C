#include "dlUtility.h"

int format_hist(TH1D* hist, int color = kBlack, int markerstyle = 20, float markersize = 1.5)
{
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerstyle);
  hist->SetMarkerSize(markersize);
  hist->SetLineColor(color);
  return 0;
}

int draw_1d_timing()
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inf = TFile::Open("../../timing/hadded_timing.root","READ");

  TH1D* data_dt = (TH1D*)(((TH3D*)inf->Get("datah3_pt_dt_lt_both"))->ProjectionY("dt_good",31,46));
  TH1D* data_lt = (TH1D*)(((TH3D*)inf->Get("datah3_pt_dt_lt_both"))->ProjectionZ("lt_good",31,46));
  TH1D* sim_dt = (TH1D*)(((TH3D*)inf->Get("simh3_pt_dt_lt_both"))->ProjectionY("dt_good",31,46));
  TH1D* sim_lt = (TH1D*)(((TH3D*)inf->Get("simh3_pt_dt_lt_both"))->ProjectionZ("lt_good",31,46));
  
  TH1D* dt_rat = (TH1D*)data_dt->Clone("dt_rat");
  TH1D* lt_rat = (TH1D*)data_lt->Clone("lt_rat");

  data_dt->Scale(1./data_dt->Integral());
  data_lt->Scale(1./data_lt->Integral());
  sim_dt->Scale(1./sim_dt->Integral());
  sim_lt->Scale(1./sim_lt->Integral());

  dt_rat->Divide(data_dt,sim_dt);
  lt_rat->Divide(data_lt,sim_lt);

  data_dt->Fit("gaus","IM");
  data_lt->Fit("gaus","IM");
  sim_dt->Fit("gaus","IM");
  sim_dt->Fit("gaus","IM");

  format_hist(data_dt);
  format_hist(dt_rat);
  format_hist(data_lt);
  format_hist(lt_rat);
  format_hist(sim_dt,kRed);
  format_hist(sim_lt,kRed);
  
  data_dt->GetFunction("gaus")->SetLineColor(kBlack);
  data_lt->GetFunction("gaus")->SetLineColor(kBlack);

  TLegend* leg = new TLegend(0.6,0.6,0.9,0.8);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->AddEntry(data_dt,"Data","p");
  leg->AddEntry(sim_dt,"Simulation","p");
  

  float w_dt_data = data_dt->GetFunction("gaus")->GetParameter(2);
  float m_dt_data = data_dt->GetFunction("gaus")->GetParameter(1);
  float w_lt_data = data_lt->GetFunction("gaus")->GetParameter(2);
  float m_lt_data = data_lt->GetFunction("gaus")->GetParameter(1);

  float w_dt_sim = sim_dt->GetFunction("gaus")->GetParameter(2);
  float m_dt_sim = sim_dt->GetFunction("gaus")->GetParameter(1);
  float w_lt_sim = sim_lt->GetFunction("gaus")->GetParameter(2);
  float m_lt_sim = sim_lt->GetFunction("gaus")->GetParameter(1);

  std::stringstream dt_mw, lt_mw;

  dt_mw << "#mu_{data}=" << std::fixed << std::setprecision(3) << m_dt_data << "#pm" << w_dt_data << " ns; #mu_{sim}=" << m_dt_sim << "#pm" << w_dt_sim;

  lt_mw << "#mu_{data}=" << std::fixed << std::setprecision(3) << m_lt_data << "#pm" << w_lt_data << " ns; #mu_{sim}=" << m_lt_sim << "#pm" << w_lt_sim;
  
  TCanvas* c = new TCanvas("","",1000,1500);

  string dt_mw_str = dt_mw.str();
  string lt_mw_str = lt_mw.str();

  string cut = "Leading jets for 30 GeV<E_{T}^{lead}<45 GeV, pass E-frac & dijet cuts";

  float md = data_dt->GetMaximum()>data_lt->GetMaximum()?data_dt->GetMaximum():data_lt->GetMaximum();
  float ms = sim_dt->GetMaximum()>sim_lt->GetMaximum()?sim_dt->GetMaximum():sim_lt->GetMaximum();
  
  float max = md>ms?md:ms;

  data_dt->GetYaxis()->SetRangeUser(0,1.1*max);
  data_lt->GetYaxis()->SetRangeUser(0,1.1*max);
  
  ratioPanelCanvas(c,0.3);
  c->cd(1);
  gPad->SetTopMargin(0.2);
  data_dt->Draw("PE");
  sim_dt->Draw("SAME PE");
  leg->Draw();
  d->cd(2);
  dt_rat->Draw("PE");
  c->cd(0);
  maintexts(0.98,0.6,0,0.02);
  drawText(dt_mw_str.c_str(),0.6,0.94,0,kBlack,0.02);
  drawText(cut.c_str(),0.3,0.94,0,kBlack,0.02);
  c->SaveAs("../../images/timing/dt_1d_simdat.png");
  c->Clear();

  ratioPanelCanvas(c,0.3);
  c->cd(1);
  gPad->SetTopMargin(0.2);
  data_lt->Draw("PE");
  sim_lt->Draw("SAME PE");
  leg->Draw();
  d->cd(2);
  lt_rat->Draw("PE");
  c->cd(0);
  maintexts(0.98,0.6,0,0.02);
  drawText(lt_mw_str.c_str(),0.6,0.94,0,kBlack,0.02);
  drawText(cut.c_str(),0.3,0.94,0,kBlack,0.02);
  c->SaveAs("../../images/timing/lt_1d_simdat.png");
  c->Clear();

  return 0;
}
