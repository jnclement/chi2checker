#include <cmath>

void Pal6()
{   // for "magenta pixels"
  const Int_t nstp = 2;
  const Int_t ncol = 2;
  Double_t stp[nstp] = { 0.0000, 1.0000 };
  Double_t red[nstp] = { 0., 0. };
  Double_t grn[nstp] = { 1., 1. };
  Double_t blu[nstp] = { 0., 0. };
  static Int_t colors[ncol];
  static Bool_t initialized = kFALSE;
  if (!initialized) {
    Int_t FI = TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
    for (int i = 0; i < ncol; i++)
      colors[i] = FI;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(ncol, colors);
}

void Pal5()
{   // for "magenta pixels"
  const Int_t nstp = 2;
  const Int_t ncol = 2;
  Double_t stp[nstp] = { 0.0000, 1.0000 };
  Double_t red[nstp] = { 255./255., 1. };
  Double_t grn[nstp] = { 0./255., 0. };
  Double_t blu[nstp] = { 255./255., 1. };
  static Int_t colors[ncol];
  static Bool_t initialized = kFALSE;
  if (!initialized) {
    Int_t FI = TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
    for (int i = 0; i < ncol; i++)
      colors[i] = FI;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(ncol, colors);
}

void Pal4()
{   // for "brown pixels"
  const Int_t nstp = 2;
  const Int_t ncol = 2;
  Double_t stp[nstp] = { 0.0000, 1.0000 };
  Double_t red[nstp] = { 120./255., 120./255. };
  Double_t grn[nstp] = { 60./255., 60./255. };
  Double_t blu[nstp] = { 0./255., 0./255. };
  static Int_t colors[ncol];
  static Bool_t initialized = kFALSE;
  if (!initialized) {
    Int_t FI = TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
    for (int i = 0; i < ncol; i++)
      colors[i] = FI;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(ncol, colors);
}
bool _dosave;
void Pal2()
{   // for "black pixels"
  const Int_t nstp = 2;
  const Int_t ncol = 3;
  Double_t stp[nstp] = { 0.0000, 1.0000 };
  Double_t red[nstp] = { 0./255., 1./255. };
  Double_t grn[nstp] = { 0./255., 1./255. };
  Double_t blu[nstp] = { 0./255., 1./255. };
  static Int_t colors[ncol];
  static Bool_t initialized = kFALSE;
  if (!initialized) {
    Int_t FI = TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
    for (int i = 0; i < ncol; i++)
      colors[i] = FI + i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(ncol, colors);
}

void Pal1()
{
  /*
  const int ncol = 50;
  const int nstp = 3;
  double red[nstp] = {0.,1.,1.};
  double grn[nstp] = {0.,1.,0.};
  double blu[nstp] = {1.,1.,0.};
  double stp[nstp] = {0,2./27,1};
  */
  const Int_t nColors = 255;
  static Int_t colors[nColors];
  const Int_t nStops = 5;
   
  // Number of colors in the final palette (for a smooth gradient)

  
  // Define the color stop positions (from 0.0 to 1.0)
   // 0.0 is the low end of the data, 1.0 is the high end
  Double_t stops[nStops] = {0.00, 0.25, 0.50, 0.75, 1.00};
  Double_t red[nStops]   = {0.9, 1.00, 1.00, 0.50, 0.0};
  Double_t green[nStops] = {1.00, 1.00, 0.00, 0.00, 0.0};
  Double_t blue[nStops]  = {0.9, 0.00, 1.00, 1.00, 0.50};
  static Bool_t initialized = kFALSE;
  if (!initialized) {
    Int_t FI = TColor::CreateGradientColorTable(nStops, stops, red, green, blue, nColors);
    for (int i = 0; i < nColors; i++)
      colors[i] = FI + i;
    initialized = kTRUE;
    return;
  }
  gStyle->SetPalette(nColors, colors);
}

void Pal3()
{
  gStyle->SetPalette(kRainbow);
}
/*
void drawText(const char *text, float xp, float yp, bool isRightAlign=0, int textColor=kBlack, double textSize=0.04, int textFont = 42, bool isNDC=true){
  // when textfont 42, textSize=0.04                                                                                                                         
  // when textfont 43, textSize=18                                                                                                                           
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  //   if(bold)tex->SetTextFont(43);                                                                                                                         
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  if(isNDC) tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
  gPad->GetListOfPrimitives()->Add(tex);
  //delete tex;
}

void sqrt_snn_text(float xp = 0.7, float yp = 0.8, bool isRightAlign=0, double textsize = 0.04)
{
  drawText("#sqrt{S_{NN}} = 200 GeV",xp,yp,isRightAlign,kBlack,textsize);
}

void sphenixtext(float xpos = 0.7, float ypos = 0.96, int ra = 0, float textsize = 0.04)
{
  drawText("#bf{#it{sPHENIX}} Internal", xpos, ypos, ra, kBlack, textsize);
}
*/
#include "dlUtility.h"
int cancount = 0;
void drawCalo(float towersem[96][256], float towersih[24][64], float towersoh[24][64], float* jet_pt, float* jet_et, float* jet_ph, float* jet_t, int jet_n, float zvtx, int failscut, int runnum, int evtnum, float* frcoh, float* frcem, float* jet_e, int isbadem[96][256], int isbadih[24][64], int isbadoh[24][64], int ishotem[96][256], int ishotih[24][64], int ishotoh[24][64], int nocalem[96][256], int nocalih[24][64], int nocaloh[24][64], float jconem[24][64], float jconih[24][64], float jconoh[24][64], int isblt, float jetcut, float chi2em[96][256], float chi2ih[24][64], float chi2oh[24][64], float emt[96][256], float iht[24][64], float oht[24][64], bool rainbow = false, int isdat = 1, int iscosmic = 0)
{

  int maxindex = 0;
  float maxE = 0;
  float slE = 0;
  float lt, st;
  int slindex = 0;
  for(int i=0; i<jet_n; ++i)
    {
      if(jet_pt[i] > maxE)
	{
	  slindex = maxindex;
	  slE = maxE;
	  st = lt;
	  maxindex = i;
	  maxE = jet_pt[i];
	  lt = jet_t[i];
	}
      else if(jet_pt[i] > slE)
	{
	  slindex = i;
	  slE = jet_pt[i];
	  st = jet_t[i];
	}
    }
  if(maxE < 50) return;
  std::stringstream full_stream;
  full_stream << std::fixed << std::setprecision(3) << "Run " << runnum << ", event " << evtnum << ", Leading EM fraction: " << frcem[maxindex] << ", OH fraction: " << frcoh[maxindex] << ", z_{vtx} = " <<  std::setprecision(1) << (zvtx==0?NAN:zvtx) << std::setprecision(3) << " E_{sl}/E_{lead}=" << jet_e[slindex]/jet_e[maxindex] << ". Tower energy scale maxes out at 25, but actual energies may be higher. Values on plots are p_{T}^{jet}. ";
  std::string full_string = full_stream.str();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetNumberContours(50);
  int ncircle = 64;
  const int ncol = 50;
  const int nstp = 3;
  double red[nstp] = {0.,1.,1.};// = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[nstp] = {0.,1.,0.};// = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[nstp] = {1.,1.,0.};// = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[nstp] = {0,2./27,1};;// = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  

  TCanvas* c = new TCanvas("","",3000,1500);
  c->Divide(4,1,0,0.1);
  
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,-0.5,23.5,64,-0.5,63.5);
  TH2D* event_disrt[3];
  
  //TH2D* times[3];
  TH2D* deads[3];
  TH2D* bchi2[3];
  TH2D* nocal[3];
  //TH2D* jcons[3];
  TH2D* chi2s[3];
  
  for(int i=0; i<3; ++i)
    {
      int nbinx = (i==0?96:24);
      int nbiny = (i==0?256:64);
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      
      deads[i] = new TH2D(("deads"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      //times[i] = new TH2D(("times"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      bchi2[i] = new TH2D(("bchi2"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      nocal[i] = new TH2D(("nocal"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      //jcons[i] = new TH2D(("jcons"+to_string(i)).c_str(),"",24,-0.5,23.5,64,-0.5,63.5);
      //chi2s[i] = new TH2D(("chi2s"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      
    }
  //TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
  event_sum->Reset();
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->Reset();
      event_disrt[i]->GetZaxis()->SetTitleOffset(0.97);
      event_disrt[i]->GetYaxis()->SetTitleOffset(1.2);
      event_disrt[i]->GetZaxis()->SetTitle("Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(-2,25);

      event_disrt[i]->GetXaxis()->SetTitleSize(0.075);
      event_disrt[i]->GetYaxis()->SetTitleSize(0.075);
      event_disrt[i]->GetZaxis()->SetTitleSize(0.075);
      event_disrt[i]->GetXaxis()->SetTitleOffset(0.65);
      event_disrt[i]->GetXaxis()->SetLabelOffset(0.01);
      event_disrt[i]->GetXaxis()->SetBinLabel(1,"-1.1");
      event_disrt[i]->GetXaxis()->SetBinLabel(i==0?48:12,"0");
      event_disrt[i]->GetXaxis()->SetBinLabel(i==0?96:24,"1.1");
      event_disrt[i]->GetYaxis()->SetBinLabel(i==0?128:32,"0");
      event_disrt[i]->GetYaxis()->SetBinLabel(i==0?192:48,"#pi/2");
      event_disrt[i]->GetYaxis()->SetBinLabel(i==0?256:64,"#pi");
      event_disrt[i]->GetYaxis()->SetBinLabel(i==0?64:16,"-#pi/2");
      event_disrt[i]->GetYaxis()->SetBinLabel(1,"-#pi");
      event_disrt[i]->GetXaxis()->LabelsOption("h");
      event_disrt[i]->GetXaxis()->SetLabelSize(0.1);
      event_disrt[i]->GetYaxis()->SetLabelSize(0.1);
      event_disrt[i]->GetZaxis()->SetLabelSize(0.075);
      event_disrt[i]->GetXaxis()->SetNdivisions(4,0,0,kFALSE);
      event_disrt[i]->GetYaxis()->SetTickLength(0);
      event_disrt[i]->GetXaxis()->SetTickLength(0);
    }

  event_sum->GetXaxis()->SetBinLabel(1,"-1.1");
  event_sum->GetXaxis()->SetBinLabel(12,"0");
  event_sum->GetXaxis()->SetBinLabel(24,"1.1");
  event_sum->GetYaxis()->SetBinLabel(32,"0");
  event_sum->GetYaxis()->SetBinLabel(48,"#pi/2");
  event_sum->GetYaxis()->SetBinLabel(64,"#pi");
  event_sum->GetYaxis()->SetBinLabel(16,"-#pi/2");
  event_sum->GetYaxis()->SetBinLabel(1,"-#pi");
  event_sum->GetXaxis()->LabelsOption("h");
  event_sum->GetXaxis()->SetTitle("Calo Sum #eta");
  event_sum->GetYaxis()->SetTitle("Calo Sum #phi");
  event_sum->GetZaxis()->SetTitleOffset(0.97);
  event_sum->GetYaxis()->SetTitleOffset(1.2);
  event_sum->GetZaxis()->SetTitle("Tower Energy [GeV]");
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  
  event_sum->GetXaxis()->SetTitleSize(0.075);
  event_sum->GetYaxis()->SetTitleSize(0.075);
  event_sum->GetZaxis()->SetTitleSize(0.075);
  event_sum->GetXaxis()->SetTitleOffset(0.65);
  event_sum->GetXaxis()->SetLabelSize(0.1);
  event_sum->GetYaxis()->SetLabelSize(0.1);
  event_sum->GetZaxis()->SetLabelSize(0.075);
  event_sum->GetXaxis()->SetLabelOffset(0.01);
  
  event_sum->GetXaxis()->SetNdivisions(4,0,0,kFALSE);
  event_sum->GetYaxis()->SetTickLength(0);
  event_sum->GetXaxis()->SetTickLength(0);
  TExec* ex1 = new TExec("ex1","Pal1();");
  TExec* ex2 = new TExec("ex2","Pal2();");
  TExec* ex3 = new TExec("ex3","Pal3();");
  TExec* ex4 = new TExec("ex4","Pal4();");
  TExec* ex5 = new TExec("ex5","Pal5();");
  
  for(int j=0; j<3; ++j)
    {

      for(int eta=0; eta<(j==0?96:24); ++eta)
	{
	  for(int phi = 0; phi<(j==0?256:64); ++phi)
	    {
	      float energy = 0;
	      if(j==0) energy = towersem[eta][phi];
	      if(j==1) energy = towersih[eta][phi];
	      if(j==2) energy = towersoh[eta][phi];
	      event_disrt[j]->Fill(eta,phi,energy);
	      if(j==0)event_sum->Fill(eta/4,phi/4,energy);
	      else event_sum->Fill(eta,phi,energy);
	      
	      if(j==0) deads[j]->Fill(eta,phi,10*ishotem[eta][phi]);
	      if(j==1) deads[j]->Fill(eta,phi,10*ishotih[eta][phi]);
	      if(j==2) deads[j]->Fill(eta,phi,10*ishotoh[eta][phi]);
	      /*
	      if(j==0 && !std::isnan(emt[eta][phi])) times[j]->Fill(eta,phi,emt[eta][phi]);
	      else times[j]->Fill(eta,phi,0);
	      if(j==1 && !std::isnan(iht[eta][phi])) times[j]->Fill(eta,phi,iht[eta][phi]);
	      else times[j]->Fill(eta,phi,0);
	      if(j==2 && !std::isnan(oht[eta][phi])) times[j]->Fill(eta,phi,oht[eta][phi]);
	      else times[j]->Fill(eta,phi,0);
	      */
	      if(j==0) bchi2[j]->Fill(eta,phi,10*isbadem[eta][phi]);
	      if(j==1) bchi2[j]->Fill(eta,phi,10*isbadih[eta][phi]);
	      if(j==2) bchi2[j]->Fill(eta,phi,10*isbadoh[eta][phi]);

	      if(j==0) nocal[j]->Fill(eta,phi,10*nocalem[eta][phi]);
	      if(j==1) nocal[j]->Fill(eta,phi,10*nocalih[eta][phi]);
	      if(j==2) nocal[j]->Fill(eta,phi,10*nocaloh[eta][phi]);
	      
	      //if(j==0 && eta < 24 && phi < 64 && jconem[eta][phi] > jetcut) jcons[j]->Fill(eta,phi,10*jconem[eta][phi]);
	      //if(j==1 && jconih[eta][phi] > jetcut) jcons[j]->Fill(eta,phi,10*jconih[eta][phi]);
	      //if(j==2 && jconoh[eta][phi] > jetcut) jcons[j]->Fill(eta,phi,10*jconoh[eta][phi]);
	      /*
	      if(j==0 && !std::isnan(chi2em[eta][phi])) chi2s[j]->Fill(eta,phi,chi2em[eta][phi]);
	      if(j==1 && !std::isnan(chi2ih[eta][phi])) chi2s[j]->Fill(eta,phi,chi2ih[eta][phi]);
	      if(j==2 && !std::isnan(chi2oh[eta][phi])) chi2s[j]->Fill(eta,phi,chi2oh[eta][phi]);
	      */
	    }
	}
      
      deads[j]->SetMaximum(2);
      bchi2[j]->SetMaximum(2);
      nocal[j]->SetMaximum(2);
      
      c->cd(j+1);
      gPad->SetLogz(0);
      gPad->SetRightMargin(j==0?0.22:0.2);
      gPad->SetLeftMargin(0.2);
      gPad->SetTopMargin(0.05);
      
      event_disrt[j]->SetContour(ncol);
      event_disrt[j]->Draw("COLZ");
      rainbow?ex3->Draw("same"):ex1->Draw("same");
      event_disrt[j]->Draw("colz same");

      /*
      
      bchi2[j]->Draw("col same0");
      ex5->Draw();
      bchi2[j]->Draw("col same0");

      nocal[j]->Draw("col same0");
      ex4->Draw();
      nocal[j]->Draw("col same0");

      deads[j]->Draw("col same0");
      ex2->Draw();
      deads[j]->Draw("col same0");
      */
    }
  
  c->cd(4);

  gPad->SetLogz(0);            
  gPad->SetTopMargin(0.05);                                                       
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  event_sum->SetContour(ncol);
  rainbow?ex3->Draw("same"):ex1->Draw("same");
  event_sum->Draw("COLZ");
  rainbow?ex3->Draw("same"):ex1->Draw("same");
  event_sum->Draw("colz same");
  c->cd(0);
  //sphenixtext(0.96,0.96,1,0.04);
  maintexts(0.97,0.5,0,0.04,isdat,1);
  std::string fails = "";
  std::string whichcut = "";
  if(failscut == 0)
    {
      fails = "Passes dijet cut";
      whichcut = "dijet";
    }
  else if(failscut == 1)
    {
      fails = "Passes frac. cut";
      whichcut = "frac";
    }
  else if(failscut==2)
    {
      fails = "Passes both";
      whichcut = "both";
    }
  else
    {
      fails = "Passes neither";
      whichcut = "fail";
    }
  full_string += fails;
  
  //drawText(full_string.c_str(),0.05,0.925,0,kBlack,0.02);

  std::stringstream secondss;
  secondss << std::fixed << std::setprecision(2) << "E_{lead}="<<jet_e[maxindex]<<", E_{sl}="<<jet_e[slindex]<<", p_{T}^{lead}=" << maxE << ", p_{T}^{sl}=" << slE << ". " << (isblt?"Bad livetime region.":"Good livetime region.") << "  t_{lead}="<<lt*17.6 << " ns, t_{sub}="<<st*17.6 << " ns.";

  std::string secondstring = secondss.str();
  
  //drawText(secondstring.c_str(),0.05,0.9,0,kBlack,0.02);
  
  c->cd(4);
  float maxJetE = 0;
  for(int k=0; k<jet_n; ++k)
    {
      if(maxJetE < jet_pt[k])
	{
	  maxJetE = jet_pt[k];
	}
      if(jet_pt[k] < 15) continue;
      std::stringstream e_stream;
      e_stream << std::fixed << std::setprecision(2) << jet_pt[k];
      std::string e_string = e_stream.str();
      float drawbin = (((jet_ph[k]+(jet_ph[k]<0?2*M_PI:0))/(2*M_PI))*64);
      drawbin += (drawbin>55?-5:5);
      drawText((e_string+" GeV").c_str(),12,drawbin,0,kBlack,0.06,42,false);

    }
  //c->Update();
  string dirstring = "";
  for(int i=60; i<131; i+=10)
    {
      if(maxJetE < i)
	{
	  dirstring = to_string(i-10)+"to"+to_string(i);
	  break;
	}
    }
  if(maxJetE > 130) dirstring = "gr130";
  float minjetdraw = 50;
  bool past = !(abs(lt-st)*17.6<2.5 && abs(lt*17.6+2.0) < 5.0);
  //if(maxJetE > minjetdraw) c->SaveAs(("../images/disp/"+to_string(runnum)+"_"+to_string(evtnum)+"_allcalo_"+dirstring+"_"+whichcut+"_"+(isblt?"blt":"glt")+"_"+(rainbow?"rainbow":"normal")+(past?"_failtime":"_passtime")+".pdf").c_str());


  for(int i=0; i<3; ++i)
    {
      c->cd(i+1);
      gPad->SetLogz();
      event_disrt[i]->GetZaxis()->SetRangeUser(0.05,25);
      gPad->Update();
    }
  c->cd(4);
  gPad->SetLogz();
  event_sum->GetZaxis()->SetRangeUser(0.05,25);
  /*
  event_sum->Draw("COLZ");
  rainbow?ex3->Draw("same"):ex1->Draw("same");
  event_sum->Draw("colz same");
  */
  gPad->Update();
  if(maxJetE > minjetdraw) c->SaveAs(("../images/disp/"+to_string(runnum)+"_"+to_string(evtnum)+"_allcalo_"+dirstring+"_"+whichcut+"_"+(isblt?"blt":"glt")+"_"+(rainbow?"rainbow":"normal")+(past?"_failtime":"_passtime")+"_log.pdf").c_str());
  /*
  for(int i=0; i<3; ++i)
    {
      c->cd(i+1);
      chi2s[i]->GetYaxis()->SetTitle("Tower #phi");
      chi2s[i]->GetXaxis()->SetTitle("Tower #eta");
      chi2s[i]->GetZaxis()->SetTitle("#chi^{2} Value");
      chi2s[i]->GetZaxis()->SetRangeUser(0.1,1e8);
      chi2s[i]->Draw("COLZ");
    }
  c->cd(4);
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  gPad->SetLogz(0);
  if(maxJetE > minjetdraw && !rainbow) c->SaveAs(("../images/disp/"+to_string(runnum)+"_"+to_string(evtnum)+"_allcalo_"+dirstring+"_"+whichcut+"_"+(isblt?"blt":"glt")+(past?"_failtime":"_passtime")+"_chi2.png").c_str());

  
  //gStyle->SetPalette(2,kBird);
  for(int i=0; i<3; ++i)
    {
      
      c->cd(i+1);
      gPad->SetLogz(0);
      jcons[i]->GetYaxis()->SetTitle("Jet Constituent #phi");
      jcons[i]->GetXaxis()->SetTitle("Jet Constituent #eta");
      jcons[i]->GetZaxis()->SetTitle("Jet Constituent T/F Value");
      jcons[i]->GetZaxis()->SetRangeUser(0,1);
      jcons[i]->Draw("COLZ");
    }
  c->cd(4);
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  gPad->SetLogz(0);
  if(maxJetE > minjetdraw && !rainbow) c->SaveAs(("../images/disp/"+to_string(runnum)+"_"+to_string(evtnum)+"_allcalo_"+dirstring+"_"+whichcut+"_"+(isblt?"blt":"glt")+(past?"_failtime":"_passtime")+"_jetcon.png").c_str());

  gStyle->SetPalette(kBird);
  for(int i=0; i<3; ++i)
    { 
      c->cd(i+1);
      gPad->SetLogz(0);
      times[i]->GetYaxis()->SetTitle("Tower #phi");
      times[i]->GetXaxis()->SetTitle("Tower #eta");
      times[i]->GetZaxis()->SetTitle("Tower Fitted Time");
      times[i]->Draw("COLZ");
    }
  if(maxJetE > minjetdraw && !rainbow) c->SaveAs(("../images/disp/"+to_string(runnum)+"_"+to_string(evtnum)+"_allcalo_"+dirstring+"_"+whichcut+"_"+(isblt?"blt":"glt")+(past?"_failtime":"_passtime")+"_times.png").c_str());
  */
  event_sum->SetName(("event_sum_"+to_string(runnum)+"_"+to_string(evtnum)).c_str());
  event_disrt[0]->SetName(("emcal_"+to_string(runnum)+"_"+to_string(evtnum)).c_str());
  event_disrt[1]->SetName(("ihcal_"+to_string(runnum)+"_"+to_string(evtnum)).c_str());
  event_disrt[2]->SetName(("ohcal_"+to_string(runnum)+"_"+to_string(evtnum)).c_str());
  
  if(_dosave)
    {
      event_sum->Write();
      for(int i=0; i<3; ++i) event_disrt[i]->Write();
    }
  
  if(c) delete c;
  if(event_disrt[0]) delete event_disrt[0];
  if(event_disrt[1]) delete event_disrt[1];
  if(event_disrt[2]) delete event_disrt[2];
  if(event_sum) delete event_sum;
  
  if(deads[0]) delete deads[0];
  if(deads[1]) delete deads[1];
  if(deads[2]) delete deads[2];
  if(bchi2[0]) delete bchi2[0];
  if(bchi2[1]) delete bchi2[1];
  if(bchi2[2]) delete bchi2[2];
  if(nocal[0]) delete nocal[0];
  if(nocal[1]) delete nocal[1];
  if(nocal[2]) delete nocal[2];
  /*
  if(jcons[0]) delete jcons[0];
  if(jcons[1]) delete jcons[1];
  if(jcons[2]) delete jcons[2];
  for(int i=0; i<3; ++i)
    {
      if(chi2s[i]) delete chi2s[i];
      if(times[i]) delete times[i];
    }
  */
  if(ex1) delete ex1;
  if(ex2) delete ex2;
  if(ex3) delete ex3;
  if(ex4) delete ex4;
  if(ex5) delete ex5;
  
}


int drawcalo(int lo, int hi, int dosave = 1, string listfilename = "chi2filesdat.txt", int rainbow = 0, int rundraw = -1, int evtdraw = -1, int isdat = 1, int iscosmic = 0)
{
  cancount = lo;
  //TFile* evtfile = TFile::Open("../chi2/hadded_chi2file_20250902.root","READ");
  _dosave = dosave;
  TFile* outf;
  string simtype = "jet10";
  if(listfilename.find("jet30") != string::npos)
    {
      simtype = "jet30";
    }
  else if(listfilename.find("jet50") != string::npos)
    {
      simtype = "jet50";
    }
  if(dosave) outf = TFile::Open(("../savedhists/savedhists_blairtest_"+to_string(lo)+"_"+to_string(hi)+"_"+(isdat?"dat":simtype)+".root").c_str(),"RECREATE");
  TTree* outt;
  if(dosave) outt = new TTree("outt","an output tree");
  if(dosave) outt->SetDirectory(outf);
  int passfrac, passdijet, outevt, outrun, outnjet;
  unsigned int outmbdhit[2];
  float jetkin[100][7];
  float outz;
  float outtime[2];
  float frac[100][2];
  if(dosave)
    {
      outt->Branch("passfrac",&passfrac,"passfrac/I");
      outt->Branch("passdijet",&passdijet,"passdijet/I");
      outt->Branch("evt",&outevt,"evt/I");
      outt->Branch("run",&outrun,"run/I");
      outt->Branch("njet",&outnjet,"njet/I");
      outt->Branch("jetkin",jetkin,"jetkin[njet][7]/F");
      outt->Branch("zvtx",&outz,"zvtx/F");
      outt->Branch("mbdhit",outmbdhit,"mbdhit[2]/i");
      outt->Branch("avgt",&outtime,"avgt[2]/F");
      outt->Branch("frac",&frac,"frac[njet][2]/F");
    }
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  int jet_n, runnum, evtnum, failscut;
  float frcem[100];
  float frcoh[100];
  float jet_e[100];
  float jet_pt[100];
  float jet_eta[100];
  float jet_phi[100];
  float jet_t[100];
  float jtem[100];
  float jtoh[100];
  float emtow[96][256];
  float ihtow[24][64];
  float ohtow[24][64];
  float zvtx;
  int isbadem[96][256] = {0};
  int isbadih[24][64] = {0};
  int isbadoh[24][64] = {0};
  int ishotem[96][256] = {0};
  int ishotih[24][64] = {0};
  int ishotoh[24][64] = {0};
  int nocalem[96][256] = {0};
  int nocalih[24][64] = {0};
  int nocaloh[24][64] = {0};
  float emt[96][256] = {0};
  float iht[24][64] = {0};
  float oht[24][64] = {0};
  float jconem[24][64];
  float jconih[24][64];
  float jconoh[24][64];
  int isblt;
  float chi2em[96][256];
  float chi2ih[24][64];
  float chi2oh[24][64];
  
  float avgt[2];
  unsigned int mbdhit[2];
  int rnwf, enwf;
  
  
  TChain* jet_tree = new TChain("jet_tree");
  //TChain* wft = new TChain("wft");
  
  string tempinfilename;
  ifstream chi2list(listfilename);

  for(int i=0; i<lo; ++i)
    {
      std::getline(chi2list,tempinfilename);
      cout << i << endl;
    }
  int counter = lo;
  while(std::getline(chi2list,tempinfilename))
    {
      //cout << "Read: " << tempinfilename << endl;
      TFile* f = TFile::Open(tempinfilename.c_str(),"READ");
      if(f && !f->IsZombie())
	{
	  jet_tree->Add(tempinfilename.c_str());
	  f->Close();
	  delete f;
	  cout << "File " << tempinfilename << " is good" << endl;
	}
      else
	{
	  if(f)
	    {
	      f->Close();
	      delete f;
	      cout << "File " << tempinfilename << " is bad" << endl;
	    }
	  else
	    {
	      cout << "File " << tempinfilename << " does not exist???" << endl;
	    }
	}
      ++counter;
      cout << counter << " " << jet_tree->GetEntries() << endl;
      if(counter >= hi)
	{
	  break;
	}
    }


  /*
  ifstream wflist("wavefiles.txt");
  for(int i=0; i<lo; ++i)
    {
      std::getline(wflist,tempinfilename);
    }
  
  int wfcounter = lo;
  while(std::getline(wflist,tempinfilename))
    {
      wft->Add(tempinfilename.c_str());
      ++wfcounter;
      if(wfcounter >= hi) break;
    }
  */
  
  jet_tree->SetBranchStatus("*",0);
  jet_tree->SetBranchStatus("jet_n",1);
  jet_tree->SetBranchStatus("runnum",1);
  jet_tree->SetBranchStatus("evtnum",1);
  jet_tree->SetBranchStatus("failscut",1);
  jet_tree->SetBranchStatus("alljetfrcem",1);
  jet_tree->SetBranchStatus("alljetfrcoh",1);
  jet_tree->SetBranchStatus("jet_et",1);
  jet_tree->SetBranchStatus("jet_pt",1);
  jet_tree->SetBranchStatus("jet_t",1);
  jet_tree->SetBranchStatus("jet_eta",1);
  jet_tree->SetBranchStatus("jet_phi",1);
  jet_tree->SetBranchStatus("emtow",1);
  jet_tree->SetBranchStatus("ihtow",1);
  jet_tree->SetBranchStatus("ohtow",1);
  jet_tree->SetBranchStatus("zvtx",1);
  jet_tree->SetBranchStatus("jet_t_oh",1);
  jet_tree->SetBranchStatus("jet_t_em",1);
  
  /*
  jet_tree->SetBranchStatus("isbadem",1);
  jet_tree->SetBranchStatus("isbadih",1);
  jet_tree->SetBranchStatus("isbadoh",1);
  jet_tree->SetBranchStatus("nocalem",1);
  jet_tree->SetBranchStatus("nocalih",1);
  jet_tree->SetBranchStatus("nocaloh",1);
  jet_tree->SetBranchStatus("ishotem",1);
  jet_tree->SetBranchStatus("ishotih",1);
  jet_tree->SetBranchStatus("ishotoh",1);
  */
  
  jet_tree->SetBranchAddress("jet_n",&jet_n);
  jet_tree->SetBranchAddress("runnum",&runnum);
  jet_tree->SetBranchAddress("evtnum",&evtnum);
  jet_tree->SetBranchAddress("failscut",&failscut);
  jet_tree->SetBranchAddress("alljetfrcem",frcem);
  jet_tree->SetBranchAddress("alljetfrcoh",frcoh);
  jet_tree->SetBranchAddress("jet_et",jet_e);
  jet_tree->SetBranchAddress("jet_pt",jet_pt);
  jet_tree->SetBranchAddress("jet_t",jet_t);
  jet_tree->SetBranchAddress("jet_eta",jet_eta);
  jet_tree->SetBranchAddress("jet_phi",jet_phi);
  jet_tree->SetBranchAddress("emtow",emtow);
  jet_tree->SetBranchAddress("ihtow",ihtow);
  jet_tree->SetBranchAddress("ohtow",ohtow);
  jet_tree->SetBranchAddress("zvtx",&zvtx);
  jet_tree->SetBranchAddress("jet_t_em",jtem);
  jet_tree->SetBranchAddress("jet_t_oh",jtoh);
  /*
  jet_tree->SetBranchAddress("isbadem",isbadem);
  jet_tree->SetBranchAddress("isbadih",isbadih);
  jet_tree->SetBranchAddress("isbadoh",isbadoh);
  jet_tree->SetBranchAddress("nocalem",nocalem);
  jet_tree->SetBranchAddress("nocalih",nocalih);
  jet_tree->SetBranchAddress("nocaloh",nocaloh);
  jet_tree->SetBranchAddress("ishotem",ishotem);
  jet_tree->SetBranchAddress("ishotih",ishotih);
  jet_tree->SetBranchAddress("ishotoh",ishotoh);
  /*
  jet_tree->SetBranchAddress("jconem",jconem);
  jet_tree->SetBranchAddress("jconih",jconih);
  jet_tree->SetBranchAddress("jconoh",jconoh);
  jet_tree->SetBranchAddress("isblt",&isblt);
  jet_tree->SetBranchAddress("chi2em",chi2em);
  jet_tree->SetBranchAddress("chi2ih",chi2ih);
  jet_tree->SetBranchAddress("chi2oh",chi2oh);
  */
  
  jet_tree->SetBranchAddress("mbdavgt",avgt);
  jet_tree->SetBranchAddress("mbdhit",mbdhit);
  /*
  wft->SetBranchAddress("runnum",&rnwf);
  wft->SetBranchAddress("evtnum",&enwf);
  wft->SetBranchAddress("emt",emt);
  wft->SetBranchAddress("iht",iht);
  wft->SetBranchAddress("oht",oht);
  */
  int nhist = 0;
  
  float jetcut = 4;
  if(dosave) outf->cd();
  int wfte = 0;
  for(int i=0; i<jet_tree->GetEntries(); ++i)
    {
      if(runnum == 51161) continue;
      if(!(i%100)) cout << i << endl;
      jet_tree->GetEntry(i);
      //wft->GetEntry(wfte);
      int flag = 0;
      /*
      while(rnwf != runnum || enwf != evtnum)
	{
	  cout << "no match" << endl;
	  ++wfte;
	  if(wfte >= wft->GetEntries())
	    {
	      flag = 1;
	      wfte = 0;
	      break;
	    }
	  wft->GetEntry(wfte);
	}
      if(flag) continue;
      ++wfte;
      */
      //if((failscut > 2 || failscut < 0)) continue; // && i % 100 != 0 && evtdraw < 0) continue;
      ++nhist;
      if(dosave)
	{
	  outmbdhit[0] = mbdhit[0];
	  outmbdhit[1] = mbdhit[1];
	  outtime[0] = avgt[0];
	  outtime[1] = avgt[1];
	  outevt = evtnum;
	  outrun = runnum;
	  outnjet = jet_n;
	  outz = zvtx;
	  if(failscut == 2 || failscut == 0)
	    {
	      passdijet = 1;
	    }
	  else
	    {
	      passdijet = 0;
	    }
	  if(failscut > 0)
	    {
	      passfrac = 1;
	    }
	  else
	    {
	      passfrac = 0;
	    }
	  float radem = 93.5;
	  float radoh = 225.87;
	  for(int j=0; j<jet_n; ++j)
	    {
	      jetkin[j][0] = jet_pt[j];
	      float radius;
	      frac[j][0] = frcem[j];
	      frac[j][1] = frcoh[j];
	      if(frcem[j] > 0.6) radius = radem;
	      else if(frcoh[j] > 0.6) radius = radoh;
	      else radius = (radem + radoh)/2;
	      float jetz = radius/(tan(2*atan(exp(-jet_eta[j]))));
	      float newz = jetz + zvtx;
	      float newtheta = atan2(radius,newz);
	      float neweta = -log(tan(0.5*newtheta));
	      jetkin[j][1] = jet_eta[j];
	      jetkin[j][2] = jet_phi[j];
	      jetkin[j][3] = jet_e[j];
	      jetkin[j][4] = jet_t[j];
	      jetkin[j][5] = jtem[j];
	      jetkin[j][6] = jtoh[j];
	    }
	  outt->Fill();
	}
      //if((evtnum != evtdraw || runnum != rundraw) && evtnum > -1 && rundraw > -1) continue;
      //drawCalo(emtow,ihtow,ohtow,jet_pt,jet_eta,jet_phi,jet_t,jet_n,zvtx,failscut,runnum,evtnum,frcoh,frcem,jet_e,isbadem,isbadih,isbadoh,ishotem,ishotih,ishotoh,nocalem,nocalih,nocaloh,jconem,jconih,jconoh,isblt,jetcut,chi2em,chi2ih,chi2oh,emt,iht,oht,rainbow?true:false,isdat,iscosmic);
    }
  if(dosave) outf->Write();
  
  return nhist;
}
      
