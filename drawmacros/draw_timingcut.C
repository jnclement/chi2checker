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

  den1->Rebin(2);
  
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
      nums1.at(i)->Rebin(2);
      effs1.push_back((TH1D*)nums1.at(i)->Clone((string(nums1.at(i)->GetName())+"_eff").c_str()));
      effs1.at(i)->Divide(nums1.at(i),den1,1,1,"B");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  ratioPanelCanvas(can,0.3);
  can->cd(1);
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.05);
  gPad->SetLogy();
  TLegend* leg = new TLegend(0.3,0.6,0.9,0.8);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetNColumns(2);

  den1->SetMarkerColor(kBlack);
  den1->SetLineColor(kBlack);
  den1->SetLineWidth(2);
  den1->SetMarkerStyle(20);
  den1->SetMarkerSize(2);
  leg->AddEntry(den1,"Dijet cut only","p");

  
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->SetMarkerColor(colors.at(i)+2);
      nums1.at(i)->SetLineColor(colors.at(i)+2);
      nums1.at(i)->SetLineWidth(2);
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      leg->AddEntry(nums1.at(i),numlabels.at(i).c_str(),"p");

      effs1.at(i)->SetMarkerColor(colors.at(i)+2);
      effs1.at(i)->SetLineColor(colors.at(i)+2);
      effs1.at(i)->SetLineWidth(2);
      effs1.at(i)->SetMarkerStyle(markers.at(i));
      effs1.at(i)->SetMarkerSize(2);
    }

  den1->Draw("PE");
  den1->GetXaxis()->SetRangeUser(30,100);
  den1->GetYaxis()->SetTitle("Counts");
  den1->GetYaxis()->SetRangeUser(0.5,den1->GetMaximum()*1.1);
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->Draw("SAME PE");
    }
  leg->Draw();
  can->cd(2);
  for(int i=0; i<effs1.size(); ++i)
    {
      effs1.at(i)->GetYaxis()->SetRangeUser(0,1.199);
      effs1.at(i)->GetXaxis()->SetRangeUser(30,100);
      effs1.at(i)->GetXaxis()->SetLabelSize(0.07);
      effs1.at(i)->GetYaxis()->SetLabelSize(0.07);
      effs1.at(i)->GetXaxis()->SetTitleSize(0.07);
      effs1.at(i)->GetYaxis()->SetTitleSize(0.07);
      if(i==0) effs1.at(i)->Draw("PE");
      else effs1.at(i)->Draw("SAME PE");
    }


  can->cd(0);

  maintexts(0.98,0.6,0,0.03);
  //drawText("Jet30,50,70 PYTHIA",0.6,0.87,0,kBlack,0.03);
  drawText("No z_{vtx} cut (non-reconstructed included)",0.05,0.87,0,kBlack,0.03);
  //drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;

  return 0;
}


int draw_timingcut()
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inf = TFile::Open("../../timing/hadded_timing.root","READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get("dath3_apt_dt_t_dijet");
  TH3D* h3_tpt_lem_loh = (TH3D*)inf->Get("simh3_apt_dtem_dtoh_both");

  std::vector<vector<int>> ybounds = {{90,109},{92,107},{88,111}};
  std::vector<vector<int>> zbounds = {{72,111},{78,107},{66,117}};
  int axis = 0;
  std::vector<int> colors = {kAzure, kSpring, kViolet};
  std::vector<int> markers = {21, 20, 71};
  std::vector<string> numlabels = {"-7 ns<t_{lead}<3 ns && |#Delta t_{lead}|<2.5 ns","-5.5 ns<t_{lead}<1.5 ns && |#Delta t_{lead}|<2 ns","-8.5 ns<t_{lead}<4.5 ns && |#Delta t_{lead}|<3 ns"};

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"timing_cut.png");
  
  return 0;
}
