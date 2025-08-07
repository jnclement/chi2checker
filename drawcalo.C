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
  const int ncol = 50;
  const int nstp = 3;
  double red[nstp] = {0.,1.,1.};
  double grn[nstp] = {0.,1.,0.};
  double blu[nstp] = {1.,1.,0.};
  double stp[nstp] = {0,2./27,1};
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
  drawText("#bf{#it{sPHENIX}} Internal", xpos, ypos, ra, kBlack, textsize);
}
int cancount = 0;
void drawCalo(float towersem[96][256], float towersih[24][64], float towersoh[24][64], float* jet_pt, float* jet_et, float* jet_ph, int jet_n, float zvtx, int failscut, int runnum, int evtnum, float* frcoh, float* frcem, float* jet_e, int isbadem[96][256], int isbadih[24][64], int isbadoh[24][64], int ishotem[96][256], int ishotih[24][64], int ishotoh[24][64], int nocalem[96][256], int nocalih[24][64], int nocaloh[24][64], float jconem[24][64], float jconih[24][64], float jconoh[24][64], int isblt, float jetcut, float chi2em[96][256], float chi2ih[24][64], float chi2oh[24][64],  bool rainbow = false)
{

  int maxindex = 0;
  float maxE = 0;
  float slE = 0;
  int slindex = 0;
  for(int i=0; i<jet_n; ++i)
    {
      if(jet_pt[i] > maxE)
	{
	  slindex = maxindex;
	  slE = maxE;
	  maxindex = i;
	  maxE = jet_pt[i];
	}
      else if(jet_pt[i] > slE)
	{
	  slindex = i;
	  slE = jet_pt[i];
	}
    }
  
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
  

  
  TCanvas* c = new TCanvas("","",1900,1000);
  c->Divide(4,1,0,0.1);
  
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,-0.5,23.5,64,-0.5,63.5);
  TH2D* event_disrt[3];
  TH2D* deads[3];
  TH2D* bchi2[3];
  TH2D* nocal[3];
  TH2D* jcons[3];
  TH2D* chi2s[3];
  for(int i=0; i<3; ++i)
    {
      int nbinx = (i==0?96:24);
      int nbiny = (i==0?256:64);
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      deads[i] = new TH2D(("deads"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      bchi2[i] = new TH2D(("bchi2"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      nocal[i] = new TH2D(("nocal"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
      jcons[i] = new TH2D(("jcons"+to_string(i)).c_str(),"",24,-0.5,23.5,64,-0.5,63.5);
      chi2s[i] = new TH2D(("chi2s"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);
    }
  //TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta Bin");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi Bin");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta Bin");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi Bin");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta Bin");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi Bin");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(1);
      event_disrt[i]->GetYaxis()->SetTitleOffset(1);
      event_disrt[i]->GetZaxis()->SetTitle("Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(-2,25);
      //event_disrt[i]->GetXaxis()->SetNdivisions(4,kFALSE);
      event_disrt[i]->GetXaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetYaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetZaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetXaxis()->SetTitleOffset(1);
      event_disrt[i]->GetXaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetYaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetZaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetXaxis()->SetLabelOffset(0.02);
    }
  event_sum->GetXaxis()->SetTitle("Calo Sum #eta Bin");
  event_sum->GetYaxis()->SetTitle("Calo Sum #phi Bin");
  event_sum->GetZaxis()->SetTitleOffset(1);
  event_sum->GetYaxis()->SetTitleOffset(1);
  event_sum->GetZaxis()->SetTitle("Tower Energy [GeV]");
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  //event_sum->GetXaxis()->SetNdivisions(4,kFALSE);
  event_sum->GetXaxis()->SetTitleSize(0.04);
  event_sum->GetYaxis()->SetTitleSize(0.04);
  event_sum->GetZaxis()->SetTitleSize(0.04);
  event_sum->GetXaxis()->SetTitleOffset(1);
  event_sum->GetXaxis()->SetLabelSize(0.04);
  event_sum->GetYaxis()->SetLabelSize(0.04);
  event_sum->GetZaxis()->SetLabelSize(0.04);
  event_sum->GetXaxis()->SetLabelOffset(0.02);
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
	      if(j==0) deads[j]->Fill(eta,phi,10*ishotem[eta][phi]);
	      if(j==1) deads[j]->Fill(eta,phi,10*ishotih[eta][phi]);
	      if(j==2) deads[j]->Fill(eta,phi,10*ishotoh[eta][phi]);

	      if(j==0) bchi2[j]->Fill(eta,phi,10*isbadem[eta][phi]);
	      if(j==1) bchi2[j]->Fill(eta,phi,10*isbadih[eta][phi]);
	      if(j==2) bchi2[j]->Fill(eta,phi,10*isbadoh[eta][phi]);

	      if(j==0) nocal[j]->Fill(eta,phi,10*nocalem[eta][phi]);
	      if(j==1) nocal[j]->Fill(eta,phi,10*nocalih[eta][phi]);
	      if(j==2) nocal[j]->Fill(eta,phi,10*nocaloh[eta][phi]);
	      
	      if(j==0 && eta < 24 && phi < 64 && jconem[eta][phi] > jetcut) jcons[j]->Fill(eta,phi,10*jconem[eta][phi]);
	      if(j==1 && jconih[eta][phi] > jetcut) jcons[j]->Fill(eta,phi,10*jconih[eta][phi]);
	      if(j==2 && jconoh[eta][phi] > jetcut) jcons[j]->Fill(eta,phi,10*jconoh[eta][phi]);

	      if(j==0) chi2s[j]->Fill(eta,phi,chi2em[eta][phi]);
	      if(j==1) chi2s[j]->Fill(eta,phi,chi2ih[eta][phi]);
	      if(j==2) chi2s[j]->Fill(eta,phi,chi2oh[eta][phi]);
	    }
	}
      deads[j]->SetMaximum(2);
      bchi2[j]->SetMaximum(2);
      nocal[j]->SetMaximum(2);
      c->cd(j+1);
      gPad->SetLogz(0);
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.2);
      gPad->SetTopMargin(0.05);
      event_disrt[j]->SetContour(ncol);
      event_disrt[j]->Draw("COLZ");
      rainbow?ex3->Draw("same"):ex1->Draw("same");
      event_disrt[j]->Draw("colz same");
      
      bchi2[j]->Draw("col same0");
      ex5->Draw();
      bchi2[j]->Draw("col same0");

      nocal[j]->Draw("col same0");
      ex4->Draw();
      nocal[j]->Draw("col same0");

      deads[j]->Draw("col same0");
      ex2->Draw();
      deads[j]->Draw("col same0");
    }

  c->cd(4);

  gPad->SetLogz(0);            
  gPad->SetTopMargin(0.05);                                                       
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  event_sum->SetContour(ncol);
  
  event_sum->Draw("COLZ");
  rainbow?ex3->Draw("same"):ex1->Draw("same");
  event_sum->Draw("colz same");
  c->cd(0);
  sphenixtext(0.96,0.96,1,0.04);
  
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
  
  drawText(full_string.c_str(),0.05,0.925,0,kBlack,0.02);

  std::stringstream secondss;
  secondss << std::fixed << std::setprecision(2) << "E_{lead}="<<jet_e[maxindex]<<", E_{sl}="<<jet_e[slindex]<<", p_{T}^{lead}=" << maxE << ", p_{T}^{sl}=" << slE << ". " << (isblt?"Bad livetime region.":"Good livetime region.");

  std::string secondstring = secondss.str();
  
  drawText(secondstring.c_str(),0.05,0.9,0,kBlack,0.02);
  
  c->cd(4);
  float maxJetE = 0;
  for(int k=0; k<jet_n; ++k)
    {
      if(maxJetE < jet_pt[k])
	{
	  maxJetE = jet_pt[k];
	}
      if(jet_pt[k] < jetcut) continue;
      std::stringstream e_stream;
      e_stream << std::fixed << std::setprecision(2) << jet_pt[k];
      std::string e_string = e_stream.str();
      drawText((e_string+" GeV").c_str(),12,((jet_ph[k]+(jet_ph[k]<0?2*M_PI:0))/(2*M_PI))*64,0,kBlack,0.04,42,false);

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
      
  if(maxJetE > 60) c->SaveAs(("../images/candidate_"+dirstring+"_"+to_string(runnum)+"_"+whichcut+"_"+to_string(evtnum)+"_"+(isblt?"blt":"glt")+"_"+(rainbow?"rainbow":"normal")+".png").c_str());


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
  if(maxJetE > 60) c->SaveAs(("../images/candidate_"+dirstring+"_"+to_string(runnum)+"_"+whichcut+"_"+to_string(evtnum)+"_"+(isblt?"blt":"glt")+"_"+(rainbow?"rainbow":"normal")+"_log.png").c_str());
  ++cancount;
  cout << "Saved" << endl;

  for(int i=0; i<3; ++i)
    {
      c->cd(i+1);
      chi2s[i]->GetYaxis()->SetTitle("Tower #phi Bin");
      chi2s[i]->GetXaxis()->SetTitle("Tower #eta Bin");
      chi2s[i]->GetZaxis()->SetTitle("#chi^{2}");
      chi2s[i]->GetZaxis()->SetRangeUser(0.1,1e8);
      chi2s[i]->Draw("COLZ");
    }
  c->cd(4);
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  gPad->SetLogz(0);
  if(maxJetE > 60 && !rainbow) c->SaveAs(("../images/candidate_"+dirstring+"_"+to_string(runnum)+"_"+whichcut+"_"+to_string(evtnum)+"_"+(isblt?"blt":"glt")+"_chi2.png").c_str());

  
  //gStyle->SetPalette(2,kBird);
  for(int i=0; i<3; ++i)
    {
      
      c->cd(i+1);
      gPad->SetLogz(0);
      jcons[i]->GetYaxis()->SetTitle("Jet Constituent #phi Bin");
      jcons[i]->GetXaxis()->SetTitle("Jet Constituent #eta Bin");
      jcons[i]->GetZaxis()->SetTitle("Jet Constituent T/F Value");
      jcons[i]->GetZaxis()->SetRangeUser(0,1);
      jcons[i]->Draw("COLZ");
    }
  c->cd(4);
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  gPad->SetLogz(0);
  if(maxJetE > 60 && !rainbow) c->SaveAs(("../images/candidate_"+dirstring+"_"+to_string(runnum)+"_"+whichcut+"_"+to_string(evtnum)+"_"+(isblt?"blt":"glt")+"_jetcon.png").c_str());
  
  if(c) delete c;
  if(event_disrt[0]) delete event_disrt[0];
  if(event_disrt[1]) delete event_disrt[1];
  if(event_disrt[2]) delete event_disrt[2];
  if(deads[0]) delete deads[0];
  if(deads[1]) delete deads[1];
  if(deads[2]) delete deads[2];
  if(bchi2[0]) delete bchi2[0];
  if(bchi2[1]) delete bchi2[1];
  if(bchi2[2]) delete bchi2[2];
  if(nocal[0]) delete nocal[0];
  if(nocal[1]) delete nocal[1];
  if(nocal[2]) delete nocal[2];
  if(jcons[0]) delete jcons[0];
  if(jcons[1]) delete jcons[1];
  if(jcons[2]) delete jcons[2];
  for(int i=0; i<3; ++i)
    {
      if(chi2s[i]) delete chi2s[i];
    }
  if(event_sum) delete event_sum;
  if(ex1) delete ex1;
  if(ex2) delete ex2;
  if(ex3) delete ex3;
  if(ex4) delete ex4;
  if(ex5) delete ex5;
}


int drawcalo(int lo, int hi, int rainbow = 0)
{
  cancount = lo;
  TFile* evtfile = TFile::Open("../events/allevents.root","READ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  int jet_n, runnum, evtnum, failscut;
  float frcem[100];
  float frcoh[100];
  float jet_e[100];
  float jet_pt[100];
  float jet_eta[100];
  float jet_phi[100];
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
  float jconem[24][64];
  float jconih[24][64];
  float jconoh[24][64];
  int isblt;
  float chi2em[96][256];
  float chi2ih[24][64];
  float chi2oh[24][64];
  
  TTree* jet_tree = (TTree*)evtfile->Get("jet_tree");

  jet_tree->SetBranchAddress("jet_n",&jet_n);
  jet_tree->SetBranchAddress("runnum",&runnum);
  jet_tree->SetBranchAddress("evtnum",&evtnum);
  jet_tree->SetBranchAddress("failscut",&failscut);
  jet_tree->SetBranchAddress("alljetfrcem",frcem);
  jet_tree->SetBranchAddress("alljetfrcoh",frcoh);
  jet_tree->SetBranchAddress("jet_et",jet_e);
  jet_tree->SetBranchAddress("jet_pt",jet_pt);
  jet_tree->SetBranchAddress("jet_eta",jet_eta);
  jet_tree->SetBranchAddress("jet_phi",jet_phi);
  jet_tree->SetBranchAddress("emtow",emtow);
  jet_tree->SetBranchAddress("ihtow",ihtow);
  jet_tree->SetBranchAddress("ohtow",ohtow);
  jet_tree->SetBranchAddress("zvtx",&zvtx);
  jet_tree->SetBranchAddress("isbadem",isbadem);
  jet_tree->SetBranchAddress("isbadih",isbadih);
  jet_tree->SetBranchAddress("isbadoh",isbadoh);
  jet_tree->SetBranchAddress("nocalem",nocalem);
  jet_tree->SetBranchAddress("nocalih",nocalih);
  jet_tree->SetBranchAddress("nocaloh",nocaloh);
  jet_tree->SetBranchAddress("ishotem",ishotem);
  jet_tree->SetBranchAddress("ishotih",ishotih);
  jet_tree->SetBranchAddress("ishotoh",ishotoh);
  jet_tree->SetBranchAddress("jconem",jconem);
  jet_tree->SetBranchAddress("jconih",jconih);
  jet_tree->SetBranchAddress("jconoh",jconoh);
  jet_tree->SetBranchAddress("isblt",&isblt);
  jet_tree->SetBranchAddress("chi2em",chi2em);
  jet_tree->SetBranchAddress("chi2ih",chi2ih);
  jet_tree->SetBranchAddress("chi2oh",chi2oh);
  float jetcut = 4;
  for(int i=lo; i<(hi>jet_tree->GetEntries()?jet_tree->GetEntries():hi); ++i)
    {
      jet_tree->GetEntry(i);
      if((failscut > 2 || failscut < 0) && i % 50 != 0) continue;
      drawCalo(emtow,ihtow,ohtow,jet_pt,jet_eta,jet_phi,jet_n,zvtx,failscut,runnum,evtnum,frcoh,frcem,jet_e,isbadem,isbadih,isbadoh,ishotem,ishotih,ishotoh,nocalem,nocalih,nocaloh,jconem,jconih,jconoh,isblt,jetcut,chi2em,chi2ih,chi2oh,rainbow?true:false);
    }

  return 0;
}
      
