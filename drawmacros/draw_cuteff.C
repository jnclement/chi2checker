#include <../dlUtility.h>

int drawprettyeff(TH3D* hist3, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, int axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, string title)
{

  if(ybounds.size() != zbounds.size())
    {
      cout << "error! bounds must be same size vectors!" << endl;
      return 1;
    }

  std::vector<TH1D*> nums1 = {};
  TH1D* den1;

  if(axis==0) den1 = hist3->ProjectionX((string(hist3->GetName())+"_den").c_str(),0,-1,0,-1,"e");
  if(axis==1) den1 = hist3->ProjectionX((string(hist3->GetName())+"_den").c_str(),0,-1,0,-1,"e");
  if(axis==2) den1 = hist3->ProjectionX((string(hist3->GetName())+"_den").c_str(),0,-1,0,-1,"e");


  
  for(int i=0; i<ybounds.size(); ++i)
    {
      int ylo = ybounds.at(i).at(0);
      int yhi = ybounds.at(i).at(1);
      int zlo = zbounds.at(i).at(0);
      int zhi = zbounds.at(i).at(1);
      if(axis==0) nums1.push_back(hist3->ProjectionX((string(hist3->GetName())+"_proj_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis==1) nums1.push_back(hist3->ProjectionY((string(hist3->GetName())+"_proj_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis==2) nums1.push_back(hist3->ProjectionZ((string(hist3->GetName())+"_proj_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
    }

  std::vector<TH1D*> effs1 = {};

  for(int i=0; i<nums1.size(); ++i)
    {
      effs1.push_back((TH1D*)nums1.at(i)->Clone((string(nums1.at(i)->GetName())+"_eff").c_str()));
      effs1.at(i)->Divide(nums1.at(i),den1,1,1,"B");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  ratioPanelCanvas(can,0.3);
  can->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.05);

  TLegend* leg = new TLegend(0.6,0.5,0.9,0.8);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetNColumns(2);

  den1->SetMarkerColor(kBlack);
  den1->SetLineColor(kBlack);
  den1->SetMarkerStyle(20);
  den1->SetMarkerSize(2);
  leg->AddEntry(den1,"No cuts","p");

  
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->SetMarkerColor(colors.at(i));
      nums1.at(i)->SetLineColor(colors.at(i));
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      leg->AddEntry(nums1.at(i),numlabels.at(i).c_str(),"p");

      effs1.at(i)->SetMarkerColor(colors.at(i));
      effs1.at(i)->SetLineColor(colors.at(i));
      effs1.at(i)->SetMarkerStyle(markers.at(i));
      effs1.at(i)->SetMarkerSize(2);
    }

  den1->Draw("PE");
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->Draw("SAME PE");
    }

  can->cd(2);
  for(int i=0; i<effs1.size(); ++i)
    {
      if(i==0) effs1.at(i)->Draw("PE");
      else effs1.at(i)->Draw("SAME PE");
    }

  leg->Draw();
  can->cd(0);

  maintexts(0.98,0.6,0,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;

  return 0;
}


int draw_cuteff()
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inf = TFile::Open("../../cuteff/hadded_cuteff.root","READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get("h3_pt_lem_loh");
  TH3D* h3_tpt_lem_loh = (TH3D*)inf->Get("h3_tpt_lem_loh");

  std::vector<vector<int>> ybounds = {{1,100,1,120},{21,120,1,120},{1,120,1,120},{1,120,1,120}};
  std::vector<vector<int>> zbounds = {{1,120,1,120},{1,120,1,120},{1,100,1,120},{21,120,1,120}};
  int axis = 0;
  std::vector<int> colors = {kSpring, kAzure, kViolet, kOrange};
  std::vector<int> markers = {20, 21, 71, 72};
  std::vector<string> numlabels = {"EM fraction < 0.9 only","EM fraction > 0.1 only","OH fraction < 0.9 only","OH Fraction > 0.1 only"};

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"h3_pt_lem_loh_effs.png");
  
  return 0;
}
