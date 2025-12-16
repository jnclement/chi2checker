#include <../dlUtility.h>
int globalusefrac = 0;
int isdat = 0;
string stype = "jet30";
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
      if(i>2) break;
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
      int nbinx = nums1.at(i)->GetNbinsX();
      while(nums1.at(i)->GetMaximum() < 8 && (nbinx%2==0 || nbinx%3 == 0 || nbinx%5 == 0))
	{
	  if(nbinx%2==0) nums1.at(i)->Rebin(2);
	  else if(nbinx%5==0) nums1.at(i)->Rebin(5);
	  else if(nbinx%3==0) nums1.at(i)->Rebin(3);
	  nbinx = nums1.at(i)->GetNbinsX();
	}
      
      //nums1.at(i)->Rebin(5);
      nums1.at(i)->Scale(1./nums1.at(i)->Integral("WIDTH"));
      if(ytitle!="")nums1.at(i)->GetYaxis()->SetTitle(ytitle.c_str());
      else nums1.at(i)->GetYaxis()->SetTitle("Integral Normalized Counts");
    }

  TCanvas* can = new TCanvas("","",1500,1500);
  can->cd();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.03);
  gPad->SetLeftMargin(0.15);
  //gPad->SetLogy();
  TLegend* leg = new TLegend(0.73,0.35,0.94,0.6);
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
      if(axis.at(i)!=1) nums1.at(i)->GetXaxis()->SetRangeUser(-10,10);
      else nums1.at(i)->GetXaxis()->SetRangeUser(-15,15);
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

  maintexts(0.9,0.73,0,0.03,isdat,0);
  if(globalusefrac) drawText("Frac cut applied",0.73,0.71,0,kBlack,0.03);
  //drawText("#Delta t \& t_{lead} cuts applied",0.625,0.69,0,kBlack,0.03);
  drawText("No reconstructed",0.73,0.79,0,kBlack,0.03);
  drawText("z_{vtx} requirement",0.73,0.75,0,kBlack,0.03);
  if(!isdat) drawText(("PYTHIA "+stype+"").c_str(),0.73,0.71,0,kBlack,0.03);
  //drawText("Truth-reco matched jets",0.05,0.91,0,kBlack,0.03);

  can->SaveAs(title.c_str());

  delete leg;
  delete can;
  delete box1;
  delete line1;
  delete line2;

  return 0;
}


int draw_spec_fake(int usefrac = 0, int mbdint = 0, int lo = 46, int hi = 60)
{
  globalusefrac = usefrac;
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  string mbdstr = "";
  if(mbdint==1) mbdstr = "_mbdboth";
  if(mbdint==2) mbdstr = "_mbdeither";
  string fracstr = usefrac?"frac":"";
  TFile* inf = TFile::Open(("../hists_mbdtimereq_out_"+stype+fracstr+(isdat?"sam":"")+".root").c_str(),"READ");

  TH3D* h3_pt_lem_loh = (TH3D*)inf->Get(("hpttmbdt_dtc"+mbdstr+stype+fracstr).c_str());

  //h3_pt_lem_loh->GetZaxis()->SetTitle("t_{lead} [ns]");
  //136,165
  std::vector<vector<int>> ybounds = {{11,20},{21,30},{31,60}};
  if(stype=="jet30" || true)
    {
      ybounds[0][0]=31;
      ybounds[0][1]=45;
      ybounds[1][0]=46;
      ybounds[1][1]=55;
      ybounds[2][0]=56;
      ybounds[2][1]=70;
    }
  std::vector<vector<int>> zbounds = {{251,350},{251,350},{251,350}};
  std::vector<int> axis = {1,1,1};
  std::vector<int> colors = {kAzure,kOrange,kMagenta};//, kAzure};//, kViolet, kOrange, kGray};
  std::vector<int> markers = {21,20,71};//, 21};//, 71, 72, 88};
  std::vector<string> numlabels = {"Jets 10<p_{T}^{uncalib}<20 GeV","Jets 20<p_{T}^{uncalib}<30 GeV","Jets 30<p_{T}^{uncalib}<45 GeV"};//"EM fraction < 0.9 only","EM fraction > 0.1 only","OH fraction < 0.9 only","OH Fraction > 0.1 only","EM frac <0.9 && OH frac > 0.1"};
  numlabels[0] = "Jets 30 GeV<p_{T}^{uncalib}<45 GeV";
  numlabels[1] = "Jets 45 GeV<p_{T}^{uncalib}<55 GeV";
  numlabels[2] = "Jets 55 GeV<p_{T}^{uncalib}<70 GeV";



  axis.at(0) = 2;
  axis.at(1) = 2;
  axis.at(2) = 2;
  
  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/"+stype+"_mbdt_proj1d"+mbdstr+".pdf","",{-10,10});


  cout << fracstr << endl;
  h3_pt_lem_loh = (TH3D*)inf->Get(("hpttdt"+stype+fracstr).c_str());

  axis.at(0) = 1;
  axis.at(1) = 1;
  axis.at(2) = 1;
  
  drawprettyeff(h3_pt_lem_loh,ybounds,zbounds,axis,colors,markers,numlabels,"../../images/mbd/"+stype+"_t_proj1d"+mbdstr+".pdf","",{-6,6});//,"dN_{jet}/d#Delta t_{l,sl} [ns^{-1}]");
  std::vector<vector<int>> ybounds2 = {{11,20},{21,30},{31,45},{46,70}};
    if(stype=="jet30" || true)
    {
      ybounds2[0][0]=31;
      ybounds2[0][1]=45;
      ybounds2[1][0]=46;
      ybounds2[1][1]=55;
      ybounds2[2][0]=56;
      ybounds2[2][1]=70;
    }

  std::vector<vector<int>> zbounds2 = {{0,-1},{0,-1},{0,-1},{0,-1}};
  std::vector<int> axis2 = {2,2,2,2};
  std::vector<int> colors2 = {kAzure,kOrange,kMagenta,kSpring};//, kViolet, kOrange, kGray};
  std::vector<int> markers2 = {21,20,71,72};
  std::vector<string> numlabels2 = {"Jets 10<p_{T}^{uncalib}<20 GeV","Jets 20<p_{T}^{uncalib}<30 GeV","Jets 30<p_{T}^{uncalib}<45 GeV","Jets 45<p_{T}^{uncalib}<70 GeV"};//"EM fraction < 0.9 only","EM fraction > 0.1 only","OH fraction < 0.9 only","OH Fraction > 0.1 only","EM frac <0.9 && OH frac > 0.1"};
  numlabels2[0] = "Jets 30 GeV<p_{T}^{uncalib}<45 GeV";
  numlabels2[1] = "Jets 45 GeV<p_{T}^{uncalib}<55 GeV";
  numlabels2[2] = "Jets 55 GeV<p_{T}^{uncalib}<70 GeV";

  drawprettyeff(h3_pt_lem_loh,ybounds2,zbounds2,axis2,colors2,markers2,numlabels2,"../../images/mbd/"+stype+"_dt_proj1d"+fracstr+".pdf","",{-3,3});

  TCanvas* can = new TCanvas("","",1500,1500);
  
  can->cd();
  gPad->SetTopMargin(0.05);
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.15);

  //gPad->SetLogz();
  TBox* box = new TBox(-3,-6,3,6);
  box->SetFillStyle(0);
  box->SetLineColor(kBlue);
  box->SetLineWidth(3);
  h3_pt_lem_loh->GetXaxis()->SetRange(lo,hi);
  TH2D* h2_t_dt = (TH2D*)(h3_pt_lem_loh->Project3D("yz"));
  //h2_t_dt->GetXaxis()->SetTitle("MBD - t_{lead} [ns]");
  h2_t_dt->GetZaxis()->SetTitle("N_{jet}");
  h2_t_dt->Rebin2D(5,5);
  gStyle->SetPalette(kInvertedDarkBodyRadiator);
  /*
  h2_t_dt->Rebin2D(4,4);
  h2_t_dt->GetXaxis()->SetRangeUser(-50,50);
  h2_t_dt->GetYaxis()->SetRangeUser(-25,25);
  h2_t_dt->GetYaxis()->SetTitleSize(0.05);
  h2_t_dt->GetXaxis()->SetTitleSize(0.045);
  */
  h2_t_dt->Draw("COLZ");
  box->Draw();
  maintexts(0.9,0.3,0,0.04,isdat,0);
  drawText(("Jets "+to_string(lo-1)+"<p_{T}^{uncalib}<"+to_string(hi)+" GeV").c_str(),0.3,0.7,0,kBlack,0.04);
  drawText("No reconstructed z_{vtx} requirement",0.3,0.75,0,kBlack,0.04);
  //drawText("PYTHIA Jet30 Sample",0.3,0.5,0,kBlack,0.04);
  can->SaveAs(("../../images/mbd/"+stype+"_t_dt_proj2d_"+to_string(lo-1)+"-"+to_string(hi)+mbdstr+".pdf").c_str());
  return 0;
}
