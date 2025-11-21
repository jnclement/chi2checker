#include <../dlUtility.h>

int drawprettyeff(std::vector<TH3D*> hist3, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, std::vector<int> axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, string title, string ytitle, std::vector<float> linex)
{

  if(ybounds.size() != zbounds.size())
    {
      cout << "error! bounds must be same size vectors!" << endl;
      return 1;
    }

  std::vector<TH1D*> nums1 = {};
  for(int j=0; j<hist3.size(); ++j)
    {
  for(int i=0; i<ybounds.size(); ++i)
    {
      int ylo = ybounds.at(i).at(0);
      int yhi = ybounds.at(i).at(1);
      int zlo = zbounds.at(i).at(0);
      int zhi = zbounds.at(i).at(1);
      cout << "test1" << endl;
      if(axis.at(i)==0) nums1.push_back(hist3.at(j)->ProjectionX((string(hist3.at(j)->GetName())+"_projx_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis.at(i)==1) nums1.push_back(hist3.at(j)->ProjectionY((string(hist3.at(j)->GetName())+"_projy_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      if(axis.at(i)==2) nums1.push_back(hist3.at(j)->ProjectionZ((string(hist3.at(j)->GetName())+"_projz_"+to_string(ylo)+"_"+to_string(yhi)+"_"+to_string(zlo)+"_"+to_string(yhi)).c_str(),ylo,yhi,zlo,zhi,"e"));
      cout << "test2" << endl;
      /*
      if(axis.at(i) != 0)
	{
	  //nums1.at(i)->Scale(1./nums1.at(i)->Integral());
	  nums1.at(i)->GetXaxis()->SetTitle("Energy Fraction");
	  nums1.at(i)->GetYaxis()->SetTitle("Self Normalized Counts");
	}
      */
      /*
      nums1.at(i)->Rebin(2);
      if(i==1) nums1.at(i)->Rebin(4);
      */
      while(nums1.at(j)->GetMaximum() < 100 && nums1.at(j)->GetNbinsX()%2==0)
	{
	  nums1.at(j)->Rebin(2);
	}
      
      
      //if(i==2) nums1.at(2)->Rebin(5);
      nums1.at(j)->Scale(1./nums1.at(j)->Integral("WIDTH"));
      if(ytitle!="")nums1.at(j)->GetYaxis()->SetTitle(ytitle.c_str());
      else nums1.at(j)->GetYaxis()->SetTitle("Integral Normalized Counts");
    }
    }
  TCanvas* can = new TCanvas("","",1500,1500);
  can->cd();
  gPad->SetTopMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  cout << "before tleg" << endl;
  TLegend* leg = new TLegend(0.62,0.35,0.92,0.6);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  float max = 0;
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->SetMarkerColor(colors.at(i)+2);
      nums1.at(i)->SetLineColor(colors.at(i)+2);
      nums1.at(i)->SetLineWidth(2);
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      nums1.at(i)->GetXaxis()->SetRangeUser(-25,25);
      if(nums1.at(i)->GetMaximum() > max) max = nums1.at(i)->GetMaximum();
      //nums1.at(i)->GetYaxis()->SetRangeUser(0.5,1e4);
      leg->AddEntry(nums1.at(i),numlabels.at(i).c_str(),"p");
    }
  max *= 1.2;
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->GetYaxis()->SetRangeUser(0,max);
      nums1.at(i)->Draw((i==0?"PE":"SAME PE"));
    }

  TLine* line1 = new TLine(linex.at(0),0,linex.at(0),max);
  TLine* line2 = new TLine(linex.at(1),0,linex.at(1),max);
  TBox* box1 = new TBox(linex.at(0),0,linex.at(1),max);
  box1->SetFillColorAlpha(kGreen+1,0.1);
  box1->SetLineWidth(0);
  leg->AddEntry(line1,"Cut bounds","l");
  leg->Draw();
  line1->SetLineStyle(9);
  line1->SetLineColor(kGreen+1);
  line1->SetLineWidth(3);
  line2->SetLineStyle(9);
  line2->SetLineColor(kGreen+1);
  line2->SetLineWidth(3);

  box1->Draw();
  line1->Draw();
  line2->Draw();

  maintexts(0.8,0.625,0,0.03,1,0);
  drawText("#Delta t \& t_{lead} cuts applied",0.625,0.69,0,kBlack,0.03);
  drawText((to_string(ybounds.at(0).at(0)-1)+" GeV<p_{T}<"+to_string(ybounds.at(0).at(1))+" GeV").c_str(),0.625,0.65,0,kBlack,0.03);
  //drawText("z_{vtx} requirement",0.625,0.61,0,kBlack,0.03);
  //drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;
  delete box1;
  delete line1;
  delete line2;

  return 0;
}


int draw_spec_fake_multihist(string num = "", int samint = 1, int mbdint = 0)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  string samstr = "";
  if(samint) samstr = "sam";
  string mbdstr = "";
  if(mbdint==1) mbdstr = "_mbdboth";
  if(mbdint==2) mbdstr = "_mbdeither";
  TFile* inf = TFile::Open(("../hists_mbdtimereq_out_dat"+samstr+".root").c_str(),"READ");


  string firsthist = "";
  if(num == "2")
    {
      firsthist = "hptdtmbdt2dat";
    }
  else
    {
      firsthist = "hptdtmbdt_ltcdat";
    }
  vector<TH3D*> h3_pt_lem_loh = {(TH3D*)inf->Get(firsthist.c_str()),(TH3D*)inf->Get(("hptdtmbdt_ltc_mbdboth"+num+"dat").c_str()),(TH3D*)inf->Get(("hptdtmbdt_ltc_mbdeither"+num+"dat").c_str())};

  //h3_pt_lem_loh->GetZaxis()->SetTitle("t_{lead} [ns]");

  std::vector<vector<int>> ybounds = {{21,40}};
  std::vector<vector<int>> zbounds = {{126,175}};
  std::vector<int> axis = {1};
  std::vector<int> colors = {kAzure,kOrange,kMagenta};//, kAzure};//, kViolet, kOrange, kGray};
  std::vector<int> markers = {21,20,71};//, 21};//, 71, 72, 88};
  std::vector<string> numlabels = {"Either or both have time","Both sides have time","Only one side has time"};//"EM fraction < 0.9 only","EM fraction > 0.1 only","OH fraction < 0.9 only","OH Fraction > 0.1 only","EM frac <0.9 && OH frac > 0.1"};

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/data_dt_proj1d"+mbdstr+samstr+"_overlay"+num+".pdf","",{-3,3});//,"dN_{jet}/d#Delta t_{l,sl} [ns^{-1}]");

  axis.at(0) = 2;

  zbounds.at(0).at(0) = 136;
  zbounds.at(0).at(1) = 165;

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/data_mbdt_proj1d"+mbdstr+samstr+to_string(ybounds.at(0).at(0)-1)+"-"+to_string(ybounds.at(0).at(1))+"_overlay"+num+".pdf","",{-3,3});

  ybounds.at(0).at(0) = 11;
  ybounds.at(0).at(1) = 20;
  
  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/data_mbdt_proj1d"+mbdstr+samstr+to_string(ybounds.at(0).at(0)-1)+"-"+to_string(ybounds.at(0).at(1))+"_overlay"+num+".pdf","",{-3,3});

  ybounds.at(0).at(0) = 31;
  ybounds.at(0).at(1) = 50;
  
  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/data_mbdt_proj1d"+mbdstr+samstr+to_string(ybounds.at(0).at(0)-1)+"-"+to_string(ybounds.at(0).at(1))+"_overlay"+num+".pdf","",{-3,3});

  ybounds.at(0).at(0) = 41;
  ybounds.at(0).at(1) = 100;

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/data_mbdt_proj1d"+mbdstr+samstr+to_string(ybounds.at(0).at(0)-1)+"-"+to_string(ybounds.at(0).at(1))+"_overlay"+num+".pdf","",{-3,3});

  TCanvas* c2 = new TCanvas("","",3000,1000);
  TH2D* hrnto = (TH2D*)(inf->Get("hrntodat"));
  TH2D* hrntj = (TH2D*)(inf->Get("hrntjdat"));

  hrnto->GetYaxis()->SetRangeUser(-15,15);
  hrntj->GetYaxis()->SetRangeUser(-15,15);
  
  TLine* lineo = new TLine(hrnto->GetXaxis()->GetBinLowEdge(1),hrnto->GetMean(2),hrnto->GetXaxis()->GetBinLowEdge(hrnto->GetNbinsX()),hrnto->GetMean(2));
  lineo->SetLineColor(kRed+1);
  lineo->SetLineWidth(2);

  cout << hrnto->GetMean(2) << endl;

  TLine* linej = new TLine(hrntj->GetXaxis()->GetBinLowEdge(1),hrntj->GetMean(2),hrntj->GetXaxis()->GetBinLowEdge(hrntj->GetNbinsX()),hrntj->GetMean(2));
  linej->SetLineColor(kRed+1);
  linej->SetLineWidth(2);

  hrnto->Draw("COLZ");
  lineo->Draw();
  c2->SaveAs("hrnto.pdf");

  hrntj->Draw("COLZ");
  linej->Draw();
  c2->SaveAs("hrntj.pdf");


  return 0;
}
