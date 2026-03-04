#include <cmath>

float get_eta_true(float ztrue, float eta0)
{
  float R = 100;
  float theta0 = 2*atan(exp(-eta0));
  float num = R*tan(theta0);
  float den = R-ztrue*tan(theta0);
  float inner = 0.5*(atan2(num,den)+eta0<0?M_PI:0);
  return -log(tan(inner));
}
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
  drawText("#bf{#it{sPHENIX}} Simulation", xpos, ypos, ra, kBlack, textsize);
}
int cancount = 0;
void drawCalo(float towersem[96][256], float towersih[24][64], float towersoh[24][64], float* jet_pt, float* jet_et, float* jet_ph, float* jet_e, int jet_n, float zvtx, float tzvtx, float* frcem, float* frcoh)
{

  int maxindex = 0;
  float maxE = 0;
  float slE = 0;
  int slindex = 0;
  float maxEtrue = 0;
  for(int i=0; i<jet_n; ++i)
    {
      if(jet_pt[i] > maxE)
	{
	  slindex = maxindex;
	  slE = maxE;
	  maxindex = i;
	  maxE = jet_pt[i];
	  if(zvtx==0)
	    {
	      maxEtrue = maxE*cosh(get_eta_true(tzvtx,jet_et[i]))/cosh(jet_et[i]);
	    }
	}
      else if(jet_pt[i] > slE)
	{
	  slindex = i;
	  slE = jet_pt[i];
	}
    }
  if(maxE < 55) return;
  std::stringstream full_stream;
  full_stream << std::fixed << std::setprecision(3) << "Leading EM fraction: " << frcem[maxindex] << ", OH fraction: " << frcoh[maxindex] << ", z_{MBD} = " <<  std::setprecision(1) << (zvtx==0?NAN:zvtx)  << ", z_{truth} = " << tzvtx << ". ";

  if(zvtx==0) full_stream << "p_{T,lead}^{uncalib} with #eta from z_{truth}: " << maxEtrue << ".";
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


  TCanvas* c = new TCanvas("","",1900,1000);
  c->Divide(4,1,0,0.1);

  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,-0.5,23.5,64,-0.5,63.5);
  TH2D* event_disrt[3];


  for(int i=0; i<3; ++i)
    {
      int nbinx = (i==0?96:24);
      int nbiny = (i==0?256:64);
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
    }
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(0.95);
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
  event_sum->GetZaxis()->SetTitleOffset(0.95);
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
  event_sum->Reset();

  TExec* ex1 = new TExec("ex1","Pal1();");
  TExec* ex2 = new TExec("ex2","Pal2();");
  TExec* ex3 = new TExec("ex3","Pal3();");
  TExec* ex4 = new TExec("ex4","Pal4();");
  TExec* ex5 = new TExec("ex5","Pal5();");

  for(int j=0; j<3; ++j)
    {
      event_disrt[j]->Reset();
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
	    }
	}

      c->cd(j+1);
      gPad->SetLogz(0);
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.2);
      gPad->SetTopMargin(0.05);

      event_disrt[j]->SetContour(ncol);
      event_disrt[j]->Draw("COLZ");
      ex1->Draw("same");
      event_disrt[j]->Draw("colz same");

    }

  c->cd(4);

  gPad->SetLogz(0);            
  gPad->SetTopMargin(0.05);                                                       
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  event_sum->SetContour(ncol);
  ex1->Draw("same");
  event_sum->Draw("COLZ");
  ex1->Draw("same");
  event_sum->Draw("colz same");
  c->cd(0);
  sphenixtext(0.96,0.96,1,0.04);
  c->cd(4);
  float maxJetE = 0;
  for(int k=0; k<jet_n; ++k)
    {
      if(maxJetE < jet_pt[k])
	{
	  maxJetE = jet_pt[k];
	}
      if(jet_pt[k] < 10) continue;
      std::stringstream pt_stream;
      pt_stream << std::fixed << std::setprecision(2) << jet_pt[k];
      std::string pt_string = pt_stream.str();
      drawText(("p_{T}="+pt_string+" GeV").c_str(),12,((jet_ph[k]+(jet_ph[k]<0?2*M_PI:0))/(2*M_PI))*64,0,kBlack,0.06,42,false);

      std::stringstream e_stream;
      e_stream <<std::fixed << std::setprecision(2) << jet_e[k];
      std::string e_string = e_stream.str();
      float drawbin = (((jet_ph[k]+(jet_ph[k]<0?2*M_PI:0))/(2*M_PI))*64);
      drawbin += (drawbin>55?-5:5);
      drawText(("E="+e_string+" GeV").c_str(),12,drawbin,0,kBlack,0.06,42,false);

    }
  //c->Update();
  float minjetdraw = 50;
  c->cd(0);
  drawText(full_string.c_str(),0.05,0.92);
  c->SaveAs(("./images/allcalo_"+to_string(maxE)+".png").c_str());
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
  gPad->Update();
  c->SaveAs(("./images/allcalo_"+to_string(maxE)+"_log.png").c_str());
  
  if(c) delete c;
  if(event_disrt[0]) delete event_disrt[0];
  if(event_disrt[1]) delete event_disrt[1];
  if(event_disrt[2]) delete event_disrt[2];
  if(event_sum) delete event_sum;
  if(ex1) delete ex1;
  if(ex2) delete ex2;
  if(ex3) delete ex3;
  if(ex4) delete ex4;
  if(ex5) delete ex5;

}

int drawcalo(int lo, int hi, int dosave = 1, string simtype = "dat")
{
  string listfilename = "inroots.list";
  TFile* outf;
  if(dosave) outf = TFile::Open(("./savedhists/savedhists_20260219_"+to_string(lo)+"_"+to_string(hi)+"_"+simtype+".root").c_str(),"RECREATE");
  TTree* outt;
  if(dosave) outt = new TTree("outt","an output tree");
  if(dosave) outt->SetDirectory(outf);
  int passfrac, passdijet, outevt, outrun, outnjet;
  unsigned int outmbdhit[2];
  float jetkin[100][5];
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
      outt->Branch("jetkin",jetkin,"jetkin[njet][5]/F");
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
  float zvtx, tzvtx;
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
  string tempinfilename;
  ifstream chi2list(listfilename);
  /*
  for(int i=0; i<lo; ++i)
    {
      std::getline(chi2list,tempinfilename);
      cout << i << endl;
    }
  */
  while(std::getline(chi2list,tempinfilename))
    {
      jet_tree->Add(tempinfilename.c_str());
    }
  cout << "finished adding to tchain" << endl;
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
  jet_tree->SetBranchStatus("zvtx",1);
  jet_tree->SetBranchStatus("tzvtx",1);
  jet_tree->SetBranchStatus("mbdavgt",1);
  jet_tree->SetBranchStatus("mbdhit",1);
  jet_tree->SetBranchStatus("emtow",1);
  jet_tree->SetBranchStatus("ihtow",1);
  jet_tree->SetBranchStatus("ohtow",1);

  jet_tree->SetBranchAddress("emtow",emtow);
  jet_tree->SetBranchAddress("ihtow",ihtow);
  jet_tree->SetBranchAddress("ohtow",ohtow);
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
  jet_tree->SetBranchAddress("zvtx",&zvtx);
  jet_tree->SetBranchAddress("tzvtx",&tzvtx);
  
  jet_tree->SetBranchAddress("mbdavgt",avgt);
  jet_tree->SetBranchAddress("mbdhit",mbdhit);
  int nhist = 0;
  
  float jetcut = 4;
  if(dosave) outf->cd();
  cout << "cded to outf" << endl;
  int wfte = 0;
  int jentries = jet_tree->GetEntries();
  for(int i=0; i<jentries; ++i)
    {
      jet_tree->GetEntry(i);
      if(jet_n < 0 || jet_n > 99)
	{
	  cout << "garbage jet_n for run " << runnum << " event " << evtnum << endl;
	  continue;
	}
      if(runnum == 51161) continue;
      if(!(i%100)) cout << i << endl;
      int flag = 0;
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
	    }
	  outt->Fill();
	}
      drawCalo(emtow,ihtow,ohtow,jet_pt,jet_eta,jet_phi,jet_e,jet_n,zvtx,tzvtx,frcem,frcoh);
    }
  if(dosave) outf->Write();
  
  return nhist;
}
      
