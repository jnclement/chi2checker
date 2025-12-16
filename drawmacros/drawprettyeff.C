#include <../dlUtility.h>
int drawprettyeff(string filename, string histname, std::vector<vector<int>> ybounds, std::vector<vector<int>> zbounds, std::vector<int> axis, std::vector<int> colors, std::vector<int> markers, std::vector<string> numlabels, string title, string ytitle, int isdat=1)
{
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  TFile* inf = TFile::Open(filename.c_str(),"READ");

  TH3D* hist3 = (TH3D*)inf->Get(histname.c_str());
  if(ybounds.size() != zbounds.size())
    {
      cout << "error! bounds must be same size vectors!" << endl;
      return 1;
    }

  std::vector<TH1D*> nums1 = {};
  
  for(int i=0; i<ybounds.size(); ++i)
    {
      if(i>3) break;
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
      /*
      nums1.at(i)->Rebin(2);
      if(i==1) nums1.at(i)->Rebin(4);
      while(nums1.at(i)->GetMaximum() < 8 && nums1.at(i)->GetNbinsX()%2==0)
	{
	  nums1.at(i)->Rebin(2);
	}
      */
      //nums1.at(i)->Rebin(5);
      nums1.at(i)->Scale(1./nums1.at(i)->Integral("WIDTH"));
      if(ytitle!="")nums1.at(i)->GetYaxis()->SetTitle(ytitle.c_str());
      else nums1.at(i)->GetYaxis()->SetTitle("Integral Normalized Counts");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  can->cd();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  TLegend* leg = new TLegend(0.64,0.35,0.92,0.55);
  leg->SetFillStyle(0);
  leg->SetLineWidth(0);
  float max = 0;
  TF1* thef[4];
  thef[0] = new TF1("thef0","gaus",-10,10);
  thef[1] = new TF1("thef1","gaus",-10,10);
  thef[2] = new TF1("thef2","gaus",-10,10);
  thef[3] = new TF1("thef3","gaus",-10,10);

  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->SetMarkerColor(colors.at(i)+2);
      nums1.at(i)->SetLineColor(colors.at(i)+2);
      nums1.at(i)->SetLineWidth(2);
      nums1.at(i)->SetMarkerStyle(markers.at(i));
      nums1.at(i)->SetMarkerSize(2);
      //nums1.at(i)->GetXaxis()->SetRangeUser(-25,25);
      nums1.at(i)->Fit("gaus","IWL","",isdat?-5:-4,isdat?1.5:2.5);
      nums1.at(i)->GetFunction("gaus")->SetLineColor(colors.at(i));
      thef[i]->SetParameter(0,nums1.at(i)->GetFunction("gaus")->GetParameter(0));
      thef[i]->SetParameter(1,nums1.at(i)->GetFunction("gaus")->GetParameter(1));
      thef[i]->SetParameter(2,nums1.at(i)->GetFunction("gaus")->GetParameter(2));
      thef[i]->SetLineColor(colors.at(i));
      thef[i]->SetLineStyle(2);
      
      if(nums1.at(i)->GetMaximum() > max) max = nums1.at(i)->GetMaximum();
      //nums1.at(i)->GetYaxis()->SetRangeUser(0.5,1e4);
      leg->AddEntry(nums1.at(i),numlabels.at(i).c_str(),"p");
    }
  max *= 1.2;
  for(int i=0; i<nums1.size(); ++i)
    {
      nums1.at(i)->GetYaxis()->SetRangeUser(0,max);
      nums1.at(i)->Draw((i==0?"PE":"SAME PE"));
      thef[i]->Draw("SAME");
    }
  

  leg->Draw();
  maintexts(0.8,0.64,0,0.03,isdat,0);
  stringstream sstr;
  sstr << std::setprecision(0) << ybounds.at(0).at(0)-1 << " GeV<p_{T}^{uncalib}<" << ybounds.at(0).at(1) << " GeV";
  string ptstr = sstr.str();
  //drawText("#Delta t \& t_{lead} cuts applied",0.625,0.69,0,kBlack,0.03);
  drawText("No reconstructed",0.64,0.69,0,kBlack,0.03);
  drawText("z_{vtx} requirement",0.64,0.65,0,kBlack,0.03);
  drawText("Dijet pair required",0.64,0.61,0,kBlack,0.03);
  drawText(ptstr.c_str(),0.64,0.57,0,kBlack,0.03);
  //drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;
  return 0;
}
