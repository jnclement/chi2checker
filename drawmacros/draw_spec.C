#include <../dlUtility.h>

int drawprettyeff(TH3D* hist3, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, std::vector<int> axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, string title, string ytitle, std::vector<float> linex)
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
      /*
      if(axis.at(i) != 0)
	{
	  //nums1.at(i)->Scale(1./nums1.at(i)->Integral());
	  nums1.at(i)->GetXaxis()->SetTitle("Energy Fraction");
	  nums1.at(i)->GetYaxis()->SetTitle("Self Normalized Counts");
	}
      */
      nums1.at(i)->Rebin(2);
      if(i==1) nums1.at(i)->Rebin(4);
      while(nums1.at(i)->GetMaximum() < 8 && nums1.at(i)->GetNbinsX()%2==0)
	{
	  nums1.at(i)->Rebin(2);
	}
      nums1.at(i)->Scale(1./nums1.at(i)->Integral("WIDTH"));
      if(ytitle!="")nums1.at(i)->GetYaxis()->SetTitle(ytitle.c_str());
      else nums1.at(i)->GetYaxis()->SetTitle("Integral Normalized Counts");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  can->cd();
  gPad->SetTopMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  TLegend* leg = new TLegend(0.62,0.55,0.92,0.8);
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

  maintexts(0.98,0.6,0,0.03);
  //drawText("Jet30,50,70 PYTHIA",0.6,0.87,0,kBlack,0.03);
  drawText("No z_{vtx} cut (non-reconstructed included)",0.05,0.87,0,kBlack,0.03);
  //drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;
  delete box1;
  delete line1;
  delete line2;

  return 0;
}


int draw_spec(int lo = 56, int hi = 70)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile* inf = TFile::Open("../../timing/hadded_timing_dat.root","READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get("dath3_apt_dt_t_dijet");

  h3_pt_lem_loh->GetZaxis()->SetTitle("t_{lead} [ns]");

  std::vector<vector<int>> ybounds = {{31,45},{46,55},{56,70}};
  std::vector<vector<int>> zbounds = {{0,-1},{0,-1},{0,-1}};
  std::vector<int> axis = {1,1,1};
  std::vector<int> colors = {kAzure,kOrange,kMagenta};//, kAzure};//, kViolet, kOrange, kGray};
  std::vector<int> markers = {21,20,71};//, 21};//, 71, 72, 88};
  std::vector<string> numlabels = {"Jets 30<p_{T}^{uncalib}<45 GeV","Jets 45<p_{T}^{uncalib}<55 GeV","Jets 55<p_{T}^{uncalib}<70 GeV"};//"EM fraction < 0.9 only","EM fraction > 0.1 only","OH fraction < 0.9 only","OH Fraction > 0.1 only","EM frac <0.9 && OH frac > 0.1"};

  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/dnp/data_dt_proj1d.png","",{-3,3});//,"dN_{jet}/d#Delta t_{l,sl} [ns^{-1}]");

  axis.at(0) = 2;
  axis.at(1) = 2;
  axis.at(2) = 2;
  
  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/dnp/data_t_proj1d.png","",{-8,4});//,"dN_{jet}/dt_{jet} [ns^{-1}]");


  TCanvas* can = new TCanvas("","",1500,750);
  
  can->cd();
  gPad->SetTopMargin(0.16);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);

  //gPad->SetLogz();
  TBox* box = new TBox(-8,-3,4,3);
  box->SetFillStyle(0);
  box->SetLineColor(kBlue);
  box->SetLineWidth(5);
  h3_pt_lem_loh->GetXaxis()->SetRange(lo,hi);
  TH2D* h2_t_dt = (TH2D*)(h3_pt_lem_loh->Project3D("yz"));
  h2_t_dt->GetXaxis()->SetTitle("t_{lead} [ns]");
  h2_t_dt->GetZaxis()->SetTitle("N_{jet}");
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  h2_t_dt->Rebin2D(4,4);
  h2_t_dt->GetXaxis()->SetRangeUser(-50,50);
  h2_t_dt->GetYaxis()->SetRangeUser(-25,25);
  h2_t_dt->Draw("COLZ");
  box->Draw();
  maintexts(0.98,0.6,0,0.03);
  drawText(("Jets "+to_string(lo-1)+"<p_{T}^{uncalib}<"+to_string(hi)+" GeV").c_str(),0.6,0.87,0,kBlack,0.03);
  drawText("No z_{vtx} cut (non-reconstructed included)",0.05,0.87,0,kBlack,0.03);
  can->SaveAs(("../../images/dnp/data_t_dt_proj2d_"+to_string(lo-1)+"-"+to_string(hi)+".png").c_str());
  return 0;
}
