#include <../dlUtility.h>

int drawprettyeff(TH3D* hist3, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, std::vector<int> axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, string title)
{

  if(ybounds.size() != zbounds.size())
    {
      cout << "error! bounds must be same size vectors!" << endl;
      return 1;
    }

  std::vector<TH1D*> nums1 = {};
  
  for(int i=0; i<ybounds.size(); ++i)
    {
      int ylo = ybounds.at(i).at(0);
      int yhi = ybounds.at(i).at(1);
      int zlo = zbounds.at(i).at(0);
      int zhi = zbounds.at(i).at(1);
      if(axis.at(i)==0) nums1.push_back(hist3->ProjectionX((string(hist3->GetName())+"_projx_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis.at(i)==1) nums1.push_back(hist3->ProjectionY((string(hist3->GetName())+"_projy_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis.at(i)==2) nums1.push_back(hist3->ProjectionZ((string(hist3->GetName())+"_projz_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis.at(i) != 0)
	{
	  nums1.at(i)->Scale(1./nums1.at(i)->Integral());
	  nums1.at(i)->GetXaxis()->SetTitle("Energy Fraction");
	  nums1.at(i)->GetYaxis()->SetTitle("Self Normalized Counts");
	}
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  can->cd();
  gPad->SetTopMargin(0.2);
  gPad->SetRightMargin(0.05);
  //gPad->SetLogy();
  TLegend* leg = new TLegend(0.3,0.5,0.6,0.75);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);

  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->SetMarkerColor(colors.at(i)+2);
      nums1.at(i)->SetLineColor(colors.at(i)+2);
      nums1.at(i)->SetLineWidth(2);
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      leg->AddEntry(nums1.at(i),numlabels.at(i).c_str(),"p");
    }

  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->Draw((i==0?"PE":"SAME PE"));
    }
  leg->Draw();

  maintexts(0.98,0.6,0,0.03);
  drawText("Jet30,50,70 PYTHIA",0.6,0.87,0,kBlack,0.03);
  drawText("No z_{vtx} cut (non-reconstructed included)",0.05,0.87,0,kBlack,0.03);
  drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;

  return 0;
}


int draw_fracs()
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inf = TFile::Open("../../cuteff/hadded_cuteff.root","READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get("h3_pt_lem_loh");
  TH3D* h3_tpt_lem_loh = (TH3D*)inf->Get("h3_tpt_lem_loh");

  std::vector<vector<int>> ybounds = {{61,70},{61,70}};
  std::vector<vector<int>> zbounds = {{1,120},{1,120}};
  std::vector<int> axis = {1,2};
  std::vector<int> colors = {kSpring, kAzure};//, kViolet, kOrange, kGray};
  std::vector<int> markers = {20, 21};//, 71, 72, 88};
  std::vector<string> numlabels = {"E_{lead} EMCal Fraction","E_{lead} OHCal Fraction"};//"EM fraction < 0.9 only","EM fraction > 0.1 only","OH fraction < 0.9 only","OH Fraction > 0.1 only","EM frac <0.9 && OH frac > 0.1"};

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"h3_pt_lem_loh_fracs.png");
  
  return 0;
}
