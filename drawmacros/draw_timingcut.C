#include <../dlUtility.h>
int _singlespec = 0;
int drawprettyeff(TH3D* hist3, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, int axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, string title)
{

  if(ybounds.size() != zbounds.size())
    {
      cout << "error! bounds must be same size vectors!" << endl;
      return 1;
    }

  std::vector<TH1D*> nums1 = {};
  std::vector<TH1D*> outs1 = {};
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

      int ybins, zbins;

      if(axis==0) ybins = hist3->GetNbinsY();
      if(axis==1) ybins = hist3->GetNbinsX();
      if(axis==2) ybins = hist3->GetNbinsX();
      if(axis==0) ybins = hist3->GetNbinsZ();
      if(axis==1) ybins = hist3->GetNbinsZ();
      if(axis==2) ybins = hist3->GetNbinsY();

      if(axis==0) outs1.push_back(hist3->ProjectionX((string(hist3->GetName())+"_outside").c_str(),0,-1,0,-1,"e"));
      if(axis==1) outs1.push_back(hist3->ProjectionY((string(hist3->GetName())+"_outside").c_str(),0,-1,0,-1,"e"));
      if(axis==2) outs1.push_back(hist3->ProjectionZ((string(hist3->GetName())+"_outside").c_str(),0,-1,0,-1,"e"));

      outs1.at(i)->Add(nums1.at(i),-1);
      outs1.at(i)->Scale(1./53.5);
    }

  std::vector<TH1D*> effs1 = {};

  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->Rebin(2);
      outs1.at(i)->Rebin(2);
      effs1.push_back((TH1D*)nums1.at(i)->Clone((string(nums1.at(i)->GetName())+"_eff").c_str()));
      effs1.at(i)->Divide(nums1.at(i),den1,1,1,"B");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  if(!_singlespec) ratioPanelCanvas(can,0.3);
  can->cd(_singlespec?1:0);
  gPad->SetTopMargin(0.1);
  gPad->SetRightMargin(0.05);
  gPad->SetLogy();
  TLegend* leg = new TLegend(_singlespec?0.5:0.25,_singlespec?0.35:0.7,_singlespec?0.83:0.93,_singlespec?0.7:0.9);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  leg->SetNColumns(2);

  den1->SetMarkerColor(kBlack);
  den1->SetLineColor(kBlack);
  den1->SetLineWidth(2);
  den1->SetMarkerStyle(20);
  den1->SetMarkerSize(2);
  leg->AddEntry(den1,"#splitline{Dijet cut only}{(jets + background)}","p");
  den1->GetXaxis()->SetTitle("Uncalibrated p_{T}^{jet} [GeV]");
  //den1->Scale(1./100);
  std::vector<TH1D*> subs1 = {};
  for(int i=0; i<nums1.size(); ++i)
    {
      //nums1.at(i)->Scale(1./100);
      //outs1.at(i)->Scale(1./100);
      nums1.at(i)->SetMarkerColor(colors.at(i)+2);
      nums1.at(i)->SetLineColor(colors.at(i)+2);
      nums1.at(i)->SetLineWidth(2);
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      subs1.push_back((TH1D*)(nums1.at(i)->Clone((nums1.at(i)->GetName()+to_string(i)+"sub").c_str())));
      subs1.at(i)->Add(outs1.at(i),-1);
      subs1.at(i)->SetMarkerColor(kMagenta+2);
      subs1.at(i)->SetLineColor(kMagenta+2);
      subs1.at(i)->SetMarkerStyle(71);
      if(!_singlespec)
	{
	  leg->AddEntry(nums1.at(i),numlabels.at(i).c_str(),"p");
	  leg->AddEntry(subs1.at(i),"Signal - surviving bg.","p");
	  //leg->AddEntry((TObject*)0,"","");
	}
      outs1.at(i)->SetMarkerColor(kOrange+2);
      outs1.at(i)->SetLineColor(kOrange+2);
      outs1.at(i)->SetLineWidth(2);
      outs1.at(i)->SetMarkerStyle(markers.at(i));
      outs1.at(i)->SetMarkerSize(2);
      if(!_singlespec) leg->AddEntry(outs1.at(i),"#splitline{Surviving background}{estimate}","p");
      
      effs1.at(i)->SetMarkerColor(colors.at(i)+2);
      effs1.at(i)->SetLineColor(colors.at(i)+2);
      effs1.at(i)->SetLineWidth(2);
      effs1.at(i)->SetMarkerStyle(markers.at(i));
      effs1.at(i)->SetMarkerSize(2);
    }

  den1->GetYaxis()->SetTitle("Counts");
  den1->GetYaxis()->SetTitleOffset(_singlespec?1.2:1.5);
  effs1.at(0)->GetXaxis()->SetTitle("Uncalibrated p_{T}^{jet} [GeV]");
  den1->GetXaxis()->SetRangeUser(10,100);
  den1->GetYaxis()->SetRangeUser(_singlespec?0.02:0.05001,den1->GetMaximum()*2);
  den1->Draw("PE");
  leg->Draw();
  if(!_singlespec)
    {
      for(int i=0; i<nums1.size(); ++i)
	{
	  nums1.at(i)->Draw("SAME PE");
	  outs1.at(i)->Draw("SAME PE");
	  subs1.at(i)->Draw("SAME PE");
	}

      /*
      can->cd(2);
      for(int i=0; i<effs1.size(); ++i)
	{
	  effs1.at(i)->GetYaxis()->SetRangeUser(0,1.199);
	  effs1.at(i)->GetXaxis()->SetRangeUser(30,100);
	  effs1.at(i)->GetXaxis()->SetLabelSize(0.1);
	  effs1.at(i)->GetYaxis()->SetLabelSize(0.1);
	  effs1.at(i)->GetXaxis()->SetTitleSize(0.1);
	  effs1.at(i)->GetYaxis()->SetTitleSize(0.1);
	  if(i==0) effs1.at(i)->Draw("PE");
	  else effs1.at(i)->Draw("SAME PE");
	}
      */
      
    }
  can->cd(0);

  maintexts(_singlespec?0.75:0.67,0.5,0,0.03,1,0);
  //drawText("Jet30,50,70 PYTHIA",0.6,0.87,0,kBlack,0.03);
  drawText("No reconstructed z_{vtx} requirement",0.5,_singlespec?0.65:0.57,0,kBlack,0.03);
  //drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;

  return 0;
}


int draw_timingcut(int singlespec = 0)
{
  _singlespec = singlespec;
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* inf = TFile::Open("../hists_mbdtimereq_out_datsam.root","READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get("hpttdtdat");
  //TH3D* h3_tpt_lem_loh = (TH3D*)inf->Get("simh3_apt_dtem_dtoh_both");

  std::vector<vector<int>> ybounds = {{111,170}};
  std::vector<vector<int>> zbounds = {{136,165}};
  int axis = 0;
  std::vector<int> colors = {kAzure};
  std::vector<int> markers = {20};
  std::vector<string> numlabels = {"-8 ns<t_{lead}<4 ns && |#Delta t|<3 ns"};

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/dnp/timing_cut_"+to_string(singlespec)+".pdf");
  
  return 0;
}
