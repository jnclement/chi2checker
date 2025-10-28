#include <dlUtility.h>

float max(float a, float b)
{
  return a>b?a:b;
}

float min(float a, float b)
{
  return a<b?a:b;
}

int format_hist(TH1D* hist, int color = kBlack, int markerstyle = 20, float markersize = 1.5)
{
  hist->SetMarkerColor(color);
  hist->SetMarkerStyle(markerstyle);
  hist->SetMarkerSize(markersize);
  hist->SetLineColor(color);
  return 0;
}

int drawtiming(int lobinx, int hibinx, int lobiny, int hibiny, int lobinz, int hibinz, int isdat = 1)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  string typestr = (isdat?"dat":"sim");
  TFile* inf = TFile::Open(("../timing/finaltiming_"+typestr+".root").c_str(),"READ");
  //TFile* datf = TFile::Open("../timing/finaltiming_dat.root","READ");

  //const int nfile = 2;
  
  //TFile* files[nfile] = {simf,datf};
  
  const int ncut = 4;
  string cuttype[ncut] = {"both","dijet","frac","all"};

  TH3D* h3_pt_adt_t[ncut][2];
  TH3D* h3_pt_dt_t[ncut][2];
  TH2D* h2_pt_t[ncut];
  TH2D* h2_lt_st[ncut];
  const int projtype = 3;
  const int rattype = 2;
  string proj3dtypes[projtype] = {"yx","zy","xz"};

  TH2D* h3_proj3d[ncut][projtype][rattype];
  TH1D* h3_projx[ncut][projtype][rattype];

  TH2D* h3_aproj3d[ncut][projtype][rattype];
  TH1D* h3_aprojx[ncut][projtype][rattype];
  
  for(int i=0; i<ncut; ++i)
    {
      h3_pt_dt_t[i][0] = (TH3D*)inf->Get(("h3_pt_dt_t_"+cuttype[i]).c_str());
      h3_pt_dt_t[i][0]->GetXaxis()->SetRange(31,100);
      h3_pt_dt_t[i][1] = (TH3D*)h3_pt_dt_t[i][0]->Clone((string(h3_pt_dt_t[i][0]->GetName())+"clone").c_str());

      h3_pt_adt_t[i][0] = (TH3D*)inf->Get(("h3_pt_adt_t_"+cuttype[i]).c_str());
      h3_pt_adt_t[i][0]->GetXaxis()->SetRange(31,100);
      h3_pt_adt_t[i][1] = (TH3D*)h3_pt_adt_t[i][0]->Clone((string(h3_pt_adt_t[i][0]->GetName())+"clone").c_str());
      
      h2_pt_t[i] = (TH2D*)inf->Get(("h2_pt_t_"+cuttype[i]).c_str());
      h2_lt_st[i] = (TH2D*)inf->Get(("h2_lt_st_"+cuttype[i]).c_str());
      for(int j=0; j<projtype; ++j)
	{
	  h3_proj3d[i][j][0] = (TH2D*)h3_pt_dt_t[i][0]->Project3D(proj3dtypes[j].c_str());
	  h3_aproj3d[i][j][0] = (TH2D*)h3_pt_adt_t[i][0]->Project3D(proj3dtypes[j].c_str());
	  h3_projx[i][j][0] = h3_proj3d[i][j][0]->ProjectionX();
	  h3_aprojx[i][j][0] = h3_aproj3d[i][j][0]->ProjectionX();
	}
      //h3_pt_dt_t[i][1]->GetXaxis()->SetRange(lobinx,hibinx);
      h3_pt_dt_t[i][1]->GetYaxis()->SetRange(lobiny,hibiny);
      h3_pt_dt_t[i][1]->GetZaxis()->SetRange(lobinz,hibinz);

      //h3_pt_adt_t[i][1]->GetYaxis()->SetRange(0,hibiny);
      //h3_pt_adt_t[i][1]->GetZaxis()->SetRange(lobinz,hibinz);
      for(int j=0; j<projtype; ++j)
	{	  
	  h3_proj3d[i][j][1] = (TH2D*)h3_pt_dt_t[i][1]->Project3D(proj3dtypes[j].c_str());
	  h3_projx[i][j][1] = h3_proj3d[i][j][1]->ProjectionX();

	  h3_aproj3d[i][j][1] = (TH2D*)h3_pt_adt_t[i][1]->Project3D(proj3dtypes[j].c_str());
	  h3_aprojx[i][j][1] = h3_aproj3d[i][j][1]->ProjectionX();
	}
    }

  TH1D* thehists[4] = {(TH1D*)h3_projx[1][0][0]->Clone((string(h3_projx[1][0][0]->GetName())+"clone").c_str()),(TH1D*)h3_projx[1][0][1]->Clone((string(h3_projx[1][0][1]->GetName())+"clone").c_str()),(TH1D*)h3_projx[3][0][0]->Clone((string(h3_projx[3][0][0]->GetName())+"clone").c_str()),(TH1D*)h3_projx[3][0][1]->Clone((string(h3_projx[3][0][1]->GetName())+"clone").c_str())};

  TH1D* h_ratio_frac = (TH1D*)h3_projx[2][0][0]->Clone("h_ratio_frac");
  h_ratio_frac->Divide(h3_projx[2][0][1],h3_projx[2][0][0],1,1,"B");

  TH1D* h_ratio_both = (TH1D*)h3_projx[0][0][0]->Clone("h_ratio_both");
  h_ratio_both->Divide(h3_projx[0][0][1],h3_projx[0][0][0],1,1,"B");
  
  h_ratio_frac->SetMarkerStyle(20);
  h_ratio_frac->SetMarkerColor(kAzure+2);
  h_ratio_frac->SetLineColor(kAzure+2);
  h_ratio_frac->SetMarkerSize(1.5);

  TH1D* h_ratio_dijetcut = (TH1D*)thehists[3]->Clone("h_ratio_dijetcut");
  TH1D* h_ratio_all = (TH1D*)thehists[3]->Clone("h_ratio_all");
  h_ratio_dijetcut->Divide(thehists[1],thehists[0],1,1,"B");
  h_ratio_all->Divide(thehists[3],thehists[2],1,1,"B");

  TH1D* h_spectrum_dijetcut = h3_projx[1][0][0];
  TH1D* h_spectrum_tdijetcut = h3_projx[1][0][1];
  TH1D* h_spectrum_all = h3_projx[3][0][0];
  TH1D* h_spectrum_tall = h3_projx[3][0][1];
  TH1D* h_spectrum_both = h3_projx[0][0][0];
  TH1D* h_spectrum_tboth = h3_projx[0][0][1];
  TH1D* h_spectrum_frac = h3_projx[2][0][0];
  TH1D* h_spectrum_tfrac = h3_projx[2][0][1];

  format_hist(h_spectrum_both,kAzure+2);
  format_hist(h_spectrum_tboth,kRed+2);
  format_hist(h_ratio_both,kRed+2);
  format_hist(h_spectrum_frac,kAzure+2);
  format_hist(h_spectrum_tfrac,kRed+2);
  format_hist(h_ratio_frac,kRed+2);
  
  TH2D* h2_dt_t_proj_all = h3_proj3d[3][1][0];
  TH2D* h2_dt_t_proj_dijet = h3_proj3d[1][1][0];

  TCanvas* c = new TCanvas("","",1000,1000);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);
  c->SetTopMargin(0.15);
  c->SetLogz();
  c->cd();

  h_ratio_frac->GetYaxis()->SetRangeUser(0,1.1);
  h_ratio_frac->Draw("PE");
  c->SaveAs(("../images/timing/"+typestr+"_ratio_frac.png").c_str());
  
  h2_dt_t_proj_all->GetZaxis()->SetRangeUser(isdat?0.5:1e-12,isdat?1e4:1);//h2_dt_t_proj_all->GetMinimum()/2,h2_dt_t_proj_all->GetMaximum()*2);
  h2_dt_t_proj_all->Draw("COLZ");


  vector<string> towrite = {"All jets with E_{T} > 30 GeV in events passing dijet or fraction cut OR with a leading E_{T}^{jet} > 50 GeV","Leading and subleading jets with E_{T} > 30 GeV in events passing dijet cut","All jets with E_{T} > 30 GeV in events passing fraction cut","All jets with E_{T} > 30 GeV in events passing both cuts"};

  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[0].c_str(),0.1,0.9,0,kBlack,0.02);

  c->SaveAs(("../images/timing/"+typestr+"_h2_dt_t_proj_all.png").c_str());
  h2_dt_t_proj_dijet->GetZaxis()->SetRangeUser(isdat?0.5:1e-12,isdat?1e4:1);//(h2_dt_t_proj_dijet->GetMinimum()/2,h2_dt_t_proj_dijet->GetMaximum()*2);
  h2_dt_t_proj_dijet->Draw("COLZ");
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[1].c_str(),0.1,0.9,0,kBlack,0.02);
  c->SaveAs(("../images/timing/"+typestr+"_h2_dt_t_proj_dijet.png").c_str());

  
  h3_aproj3d[1][0][0]->GetZaxis()->SetRangeUser(isdat?0.5:1e-12,isdat?1e4:1);
  h3_aproj3d[1][0][0]->GetXaxis()->SetRangeUser(30,100);
  h3_aproj3d[1][0][0]->Draw("COLZ");
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[1].c_str(),0.1,0.9,0,kBlack,0.02);
  c->SaveAs(("../images/timing/"+typestr+"_h3_pt_adt_dijet.png").c_str());


  h3_aproj3d[3][0][0]->GetZaxis()->SetRangeUser(isdat?0.5:1e-12,isdat?1e4:1);
  h3_aproj3d[3][0][0]->GetXaxis()->SetRangeUser(30,100);
  h3_aproj3d[3][0][0]->Draw("COLZ");
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[0].c_str(),0.1,0.9,0,kBlack,0.02);
  c->SaveAs(("../images/timing/"+typestr+"_h3_pt_adt_all.png").c_str());
  
  /*
  if(!isdat)
    {
      for(int i=42; i<71; ++i)
	{
	  h_spectrum_tdijetcut->SetBinContent(i,h_spectrum_tdijetcut->GetBinContent(i)*15);
	  h_spectrum_tall->SetBinContent(i,h_spectrum_tall->GetBinContent(i)*15);
	  h_spectrum_dijetcut->SetBinContent(i,h_spectrum_dijetcut->GetBinContent(i)*15);
	  h_spectrum_all->SetBinContent(i,h_spectrum_all->GetBinContent(i)*15);
	}
    }
  */
  TCanvas* d = new TCanvas("","",1500,1500);
  ratioPanelCanvas(d, 0.3);

  TLegend* leg = new TLegend(0.5,0.6,0.9,0.8);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  d->cd(1);
  gPad->SetTopMargin(0.2);
  h_ratio_dijetcut->SetMarkerColor(kRed+2);
  h_ratio_dijetcut->SetMarkerStyle(20);
  h_ratio_dijetcut->SetMarkerSize(2);
  h_ratio_dijetcut->SetLineColor(kRed+2);
  h_ratio_all->SetMarkerColor(kRed+2);
  h_ratio_all->SetMarkerStyle(20);
  h_ratio_all->SetMarkerSize(2);
  h_ratio_all->SetLineColor(kRed+2);
  h_spectrum_tdijetcut->SetMarkerColor(kRed+2);
  h_spectrum_tdijetcut->SetMarkerStyle(20);
  h_spectrum_tdijetcut->SetMarkerSize(2);
  h_spectrum_tdijetcut->SetLineColor(kRed+2);
  h_spectrum_tall->SetMarkerColor(kRed+2);
  h_spectrum_tall->SetMarkerStyle(20);
  h_spectrum_tall->SetMarkerSize(2);
  h_spectrum_tall->SetLineColor(kRed+2);
  
  h_spectrum_all->SetMarkerColor(kAzure+2);
  h_spectrum_all->SetMarkerStyle(20);
  h_spectrum_all->SetMarkerSize(2);
  h_spectrum_all->SetLineColor(kAzure+2);
  h_spectrum_dijetcut->SetMarkerColor(kAzure+2);
  h_spectrum_dijetcut->SetMarkerStyle(20);
  h_spectrum_dijetcut->SetMarkerSize(2);
  h_spectrum_dijetcut->SetLineColor(kAzure+2);

  leg->AddEntry(h_spectrum_tdijetcut,"Timing cut","p");
  leg->AddEntry(h_spectrum_dijetcut,"No timing cut","p");
  
  d->cd();
  d->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetLogy();
  float maxval = isdat?1e4:1;//2*max(h_spectrum_tdijetcut->GetMaximum(),h_spectrum_dijetcut->GetMaximum());
  float minval = isdat?0.5:1e-12;//2*min(h_spectrum_tdijetcut->GetMinimum(),h_spectrum_dijetcut->GetMinimum());
  h_spectrum_dijetcut->GetYaxis()->SetRangeUser(minval,maxval);
  h_spectrum_dijetcut->Draw("PE");
  h_spectrum_tdijetcut->Draw("SAME PE");
  leg->Draw();
  d->cd(2);
  h_ratio_dijetcut->Draw("PE");
  d->cd(0);
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[1].c_str(),0.1,0.9,0,kBlack,0.02);
  //d->Update();
  d->SaveAs(("../images/timing/"+typestr+"_ratio_dijet.png").c_str());



  d->Clear();
  ratioPanelCanvas(d, 0.3);
  d->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetLogy();
  h_spectrum_all->GetYaxis()->SetRangeUser(minval,maxval);
  h_spectrum_all->Draw("PE");
  h_spectrum_tall->Draw("SAME PE");
  leg->Draw();
  d->cd(2);
  
  h_ratio_all->Draw("PE");
  d->cd(0);
  
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[0].c_str(),0.1,0.9,0,kBlack,0.02);
  d->SaveAs(("../images/timing/"+typestr+"_ratio_all.png").c_str());

  d->Clear();
  ratioPanelCanvas(d, 0.3);
  d->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetLogy();
  h_spectrum_frac->GetYaxis()->SetRangeUser(minval,maxval);
  h_spectrum_frac->Draw("PE");
  h_spectrum_tfrac->Draw("SAME PE");
  leg->Draw();
  d->cd(2);
  
  h_ratio_frac->Draw("PE");
  d->cd(0);
  
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[2].c_str(),0.1,0.9,0,kBlack,0.02);
  d->SaveAs(("../images/timing/"+typestr+"_ratio_frac.png").c_str());



  d->Clear();
  ratioPanelCanvas(d, 0.3);
  d->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetLogy();
  h_spectrum_both->GetYaxis()->SetRangeUser(minval,maxval);
  h_spectrum_both->Draw("PE");
  h_spectrum_tboth->Draw("SAME PE");
  leg->Draw();
  d->cd(2);
  h_ratio_both->Draw("PE");
  d->cd(0);
  
  sphenixtext(0.7,0.98,0,0.02);
  sqrt_s_text(0.7,0.96,0,0.02);
  antikt_text(0.4,0.7,0.94,0,0.02);
  drawText((isdat?"Data":"Simulation"),0.7,0.92,0,kBlack,0.02);
  drawText(towrite[3].c_str(),0.1,0.9,0,kBlack,0.02);
  d->SaveAs(("../images/timing/"+typestr+"_ratio_both.png").c_str());

  return 0;
}
