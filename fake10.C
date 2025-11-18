int fake10(string filename, bool issim = false, bool dofrac = false, string simstr = "dat") {

  if(dofrac) simstr += "frac";
  
  gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetOptDate(111);
  
   TFile *f  = TFile::Open(filename.c_str(),"READ");
   TTree *outt = static_cast <TTree *> (f->Get("outt"));

   Int_t           passfrac;
   Int_t           passdijet;
   Int_t           evt;
   Int_t           run;
   Int_t           njet;
   Float_t         jetkin[60][5];
   Float_t         zvtx;
   Int_t           mbdhit[2];
   Float_t         avgt[2];
   Float_t         frac[60][2];   
   
   // Set branch addresses.
   outt->SetBranchAddress("passfrac",&passfrac);
   outt->SetBranchAddress("passdijet",&passdijet);
   outt->SetBranchAddress("evt",&evt);
   outt->SetBranchAddress("run",&run);
   outt->SetBranchAddress("njet",&njet);
   outt->SetBranchAddress("jetkin",jetkin);
   outt->SetBranchAddress("zvtx",&zvtx);
   outt->SetBranchAddress("mbdhit",&mbdhit);
   outt->SetBranchAddress("avgt",&avgt);
   outt->SetBranchAddress("frac",&frac);   

   //=======================
   // DECLARE HISTOGRAMS
   //=======================
   
   int nbins = 35;

   float binb[36];
   for (int i=0;i<25;i++) binb[i] = 0.0 + 1.0 * (float) i;
   for (int i=0;i<7;i++) binb[i+25] = 25.0 + 2.5 * (float) i;
   for (int i=0;i<4;i++) binb[i+26+6] = 45.0 + 5.0 * (float) i;      

   for (int i=0;i<36;i++) cout << binb[i] << endl;
   TH1D *hspectra[20];
   char fooname[100];
   for (int i=0;i<20;i++) {
     snprintf(fooname,100,"hspectra%d",i);
     hspectra[i]  = new TH1D(("hspectra"+to_string(i)+simstr).c_str(),fooname,nbins,binb);
   }

   int npx = 0;
   int ncx = 0;
   int npd = 0;
   // no other cuts...
   TH1D *hleadtimeNOMBD = new TH1D(("hleadtimeNOMBD"+simstr).c_str(),"hleadtimeNOMBD",300,-30.0,30.0);
   TH1D *hleadtimeYESMBD = new TH1D(("hleadtimeYESMBD"+simstr).c_str(),"hleadtimeYESMBD",300,-30.0,30.0);   
   TH1D *hleadtimeNOMBDwdijet = new TH1D(("hleadtimeNOMBDwdijet"+simstr).c_str(),"hleadtimeNOMBDwdijet",300,-30.0,30.0);
   TH1D *hleadtimeYESMBDwdijet = new TH1D(("hleadtimeYESMBDwdijet"+simstr).c_str(),"hleadtimeYESMBDwdijet",300,-30.0,30.0);   
   TH1D *hleadtimeNOMBDwdijetP = new TH1D(("hleadtimeNOMBDwdijetP"+simstr).c_str(),"hleadtimeNOMBDwdijetP",300,-30.0,30.0);
   TH1D *hleadtimeYESMBDwdijetP = new TH1D(("hleadtimeYESMBDwdijetP"+simstr).c_str(),"hleadtimeYESMBDwdijetP",300,-30.0,30.0);
   TH3D* htdtmbdt = new TH3D(("htdtmbdt"+simstr).c_str(),";Jet t [ns];#Delta t [ns]; MBD - t_{lead} [ns];Counts",300,-30,30,300,-30,30,300,-30,30);
   TH3D* hpttmbdt = new TH3D(("hpttmbdt"+simstr).c_str(),";p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,300,-30,30,300,-30,30);
   TH3D* hptdtmbdt = new TH3D(("hptdtmbdt"+simstr).c_str(),";p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,300,-30,30,300,-30,30);

   TH3D* hpttmbdt_dtc = new TH3D(("hpttmbdt_dtc"+simstr).c_str(),";p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,300,-30,30,300,-30,30);
   TH3D* hptdtmbdt_ltc = new TH3D(("hptdtmbdt_ltc"+simstr).c_str(),";p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,300,-30,30,300,-30,30);

   TH3D* hpttmbdt_dtc_mbdboth = new TH3D(("hpttmbdt_dtc_mbdboth"+simstr).c_str(),";p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,300,-30,30,300,-30,30);
   TH3D* hptdtmbdt_ltc_mbdboth = new TH3D(("hptdtmbdt_ltc_mbdboth"+simstr).c_str(),";p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,300,-30,30,300,-30,30);
   
   TH1D *hfracgood0 = new TH1D(("hfracgood0"+simstr).c_str(),"hfracgood0",100,-0.5,1.5);
   TH1D *hfracbad0 = new TH1D(("hfracbad0"+simstr).c_str(),"hfracbad0",100,-0.5,1.5);   
   
   //=========================================================
   // LOOP OVER EVENTS
   //=========================================================
   
   Long64_t nentries = outt->GetEntries();
   Long64_t nbytes = 0;
   for (Long64_t i=0; i<nentries;i++) {
      nbytes += outt->GetEntry(i);

      if (run == 51161) continue;  // somehow a bad apple

      // spectra of leading jets or all jets (does this make any difference?)
      // 1. find the leading jet index (kin index 0,1,2,3,4 = jet pt,eta,phi,E,t)
      //    - does it pass the timing cuts (regular, tight) 
      // 2. find the subleading jet index
      // 3. 

      // (( 1 ))
	
      int leadingjetindex = 0;
      for (int ijet=0;ijet<njet;ijet++) {
	double jetpt = jetkin[ijet][0];
	if (jetpt > jetkin[leadingjetindex][0]) {
	  leadingjetindex = ijet;
	}
      }
      // minimum jet pt > 10 GeV, and |eta| < 0.7
      //      if (! (jetkin[leadingjetindex][0] >= 10.0 && TMath::Abs(jetkin[leadingjetindex][1]) <=  0.7) ) continue;
      if (! (jetkin[leadingjetindex][0] >= 10.0) ) continue;      

      // (( 2 ))

      // check on (a) leading time status
      // check on (b) MBD hit status

      bool PassLeadTime = true;
      bool PassLeadTimeONLY = false;      
      bool MBDboth = true;
      bool MBDeither = false;

      bool PassMinimalFrac = true;
      if (frac[leadingjetindex][0] < 0.10 || frac[leadingjetindex][0]>0.9 ||
	  frac[leadingjetindex][1] < 0.10 || frac[leadingjetindex][1]>0.9) PassMinimalFrac = false;

      double jetleadtime = 17.6*(jetkin[leadingjetindex][4]);
      //if (TMath::Abs(jetleadtime + 2.0) < 6.0 || issim) PassLeadTime = true;
      if (TMath::Abs(jetleadtime + 2.0) < 6.0 || issim) PassLeadTimeONLY = true;      
      hleadtimeNOMBD->Fill(jetleadtime);
      
      // now check if MBD time is available, and if so, apply the stricter cut
      double jetleadtimeMBD = -9999.0;
      
      // Joey has NaN entries for avgt is not available (bad practice)
      double mbdtime = -9999.0;
      double mbdoffset = 0.0;
      if (std::isnan(avgt[0])) avgt[0] = -9999.0;
      if (std::isnan(avgt[1])) avgt[1] = -9999.0;	  
      
      mbdoffset = -2.07;
      if (run > 49375 && run < 49700) mbdoffset += +2.5;
      if (run > 48600 && run < 49375) mbdoffset += -7.5;
      if (run > 48180 && run < 48270) mbdoffset += -7.5;
      if (run > 48050 && run < 48090) mbdoffset += -7.5;

      if (avgt[0] > -99 && avgt[1] < -99) mbdtime = avgt[0];
      if (avgt[0] < -99 && avgt[1] > -99) mbdtime = avgt[1];
      if (avgt[0] > -99 || avgt[1] > -99) MBDeither = true;

      if (avgt[0] > -99 && avgt[1] > -99) {
	mbdtime = 0.5*(avgt[0]+avgt[1]);
	//	MBDboth = true;
      }

      jetleadtimeMBD = mbdoffset + mbdtime - 17.6*jetkin[leadingjetindex][4];
      bool PassLeadTimeWIDE = true;
      if (TMath::Abs(jetleadtimeMBD) > 3.0 && !issim /*&& mbdtime > -99*/) PassLeadTime = false;
      if (TMath::Abs(jetleadtimeMBD) > 10.0 && !issim /*&& mbdtime > -99*/) PassLeadTimeWIDE = false;
      hleadtimeYESMBD->Fill(jetleadtimeMBD);

      if (zvtx == 0) MBDboth = false;

      // jet pT spectra - no cuts, with time requirement, with MBD requirements
      //hspectra[0]->Fill(jetkin[leadingjetindex][0]);   //  raw
      //if (PassLeadTime) hspectra[1]->Fill(jetkin[leadingjetindex][0]);   //  lead time cut
      //      if (PassLeadTime && MBDeither) hspectra[2]->Fill(jetkin[leadingjetindex][0]);   //  lead time cut
      //if (PassLeadTime && MBDboth) hspectra[3]->Fill(jetkin[leadingjetindex][0]);   //  lead time cut            

      // now is there a dijet partner (what requirements do we want to apply....)

      int subleadingjetindex = 0;
      for (int ijet=0;ijet<njet;ijet++) {
	double jetpt = jetkin[ijet][0];
	if (jetpt < jetkin[leadingjetindex][0] && jetpt >= jetkin[subleadingjetindex][0]) {
	  subleadingjetindex = ijet;
	}
      }
      // minimum subleading jet, xJ > 0.3 (energy)
      bool DijetPartner = false;
      bool DijetAndt = false;
      if (jetkin[subleadingjetindex][3]/jetkin[leadingjetindex][3] > 0.3)
	{
	  DijetPartner = true;
	  DijetAndt = true;
	}

      if (DijetPartner) {
	++npx;
	// Returns the difference in the range [-pi, pi]
	double phi1 = jetkin[leadingjetindex][2];
	double phi2 = jetkin[subleadingjetindex][2];
	const double two_pi = 2.0 * TMath::Pi();
	double diff = fmod(phi2 - phi1 + TMath::Pi(), two_pi);
	if (diff < 0) diff += two_pi;
	double wrap_radians = diff - TMath::Pi();
	//if (!(wrap_radians > (3.0/4.0)*TMath::Pi())) DijetPartner = false;   // fails if not on opposite side 3pi/4
	float dphi = fabs(phi1-phi2);
	if(dphi > M_PI) dphi = 2*M_PI-dphi;
	if(dphi < 3*M_PI/4)
	  {
	    DijetAndt = false;
	    DijetPartner = false;
	  }
	if(DijetPartner)
	  {
	    ++npd;
	  }
	else
	  {
	    ++ncx;
	  }
      }


      // dijet timing cut
      double dijetTimediff = 17.6*jetkin[leadingjetindex][4] - 17.6*jetkin[subleadingjetindex][4];
      if(!issim)
	{

	  if (TMath::Abs(dijetTimediff) > 3.0) DijetAndt = false;
	}

      htdtmbdt->Fill(jetleadtime,dijetTimediff,jetleadtimeMBD);
      hpttmbdt->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetAndt) hpttmbdt_dtc->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetAndt & MBDboth) hpttmbdt_dtc_mbdboth->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetPartner)
	{
	  hptdtmbdt->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	  if(PassLeadTimeONLY) hptdtmbdt_ltc->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	  if(PassLeadTimeONLY && MBDboth) hptdtmbdt_ltc_mbdboth->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	}
      //========================
      // IF DijetPartner = true (then passed xJ cut, delta-phi cut, delta-t cut on pair!
      //========================      

      if (DijetAndt) {
	hleadtimeNOMBDwdijet->Fill(jetleadtime);
	hleadtimeYESMBDwdijet->Fill(jetleadtimeMBD);
	if (PassLeadTimeONLY) hleadtimeNOMBDwdijetP->Fill(jetleadtime);
	if (MBDeither && PassLeadTime) hleadtimeYESMBDwdijetP->Fill(jetleadtimeMBD);			
      }

      if(!dofrac || (dofrac && PassMinimalFrac))
	{
	  
      if (MBDboth && PassLeadTime && DijetPartner && PassMinimalFrac) hspectra[15]->Fill(jetkin[leadingjetindex][0]);
      
      if (PassLeadTimeONLY && DijetAndt && PassLeadTimeWIDE) hspectra[5]->Fill(jetkin[leadingjetindex][0]);
      if (PassLeadTimeONLY && DijetAndt && PassLeadTime) hspectra[1]->Fill(jetkin[leadingjetindex][0]);
      if (PassLeadTimeONLY && DijetAndt) hspectra[0]->Fill(jetkin[leadingjetindex][0]);
      
      if (MBDboth && PassLeadTimeWIDE && DijetAndt && PassLeadTimeONLY) hspectra[4]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTime && DijetAndt && PassLeadTimeONLY) hspectra[3]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTimeONLY && DijetAndt) hspectra[2]->Fill(jetkin[leadingjetindex][0]);



      if (MBDboth && PassLeadTime && DijetAndt && PassLeadTimeONLY) hspectra[7]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTime) hspectra[6]->Fill(jetkin[leadingjetindex][0]);      

      if (MBDboth && PassLeadTime && PassMinimalFrac) hspectra[16]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTime && PassLeadTimeONLY && PassMinimalFrac && DijetAndt) hspectra[17]->Fill(jetkin[leadingjetindex][0]);      
      
      if (DijetAndt && PassLeadTimeONLY) {
	hspectra[8]->Fill(jetkin[leadingjetindex][0]);
	if (MBDboth) {
	  hspectra[9]->Fill(jetkin[leadingjetindex][0]);
	  hfracgood0->Fill(frac[leadingjetindex][0]);
	} else {
	  hfracbad0->Fill(frac[leadingjetindex][0]);
	}
      }

      if(PassLeadTimeONLY)
	{
	  hspectra[10]->Fill(jetkin[leadingjetindex][0]);
	  if(MBDboth)
	    {
	      hspectra[11]->Fill(jetkin[leadingjetindex][0]);
	    }
	  if(DijetAndt)
	    {
	      hspectra[12]->Fill(jetkin[leadingjetindex][0]);
	    }
	}

      if(PassLeadTimeONLY && PassMinimalFrac)
	{
	  if(MBDboth)
	    {
	      hspectra[13]->Fill(jetkin[leadingjetindex][0]);
	    }
	  if(DijetAndt)
	    {
	      hspectra[14]->Fill(jetkin[leadingjetindex][0]);
	    }
	}
      
      if (DijetAndt && PassLeadTimeONLY && PassMinimalFrac) {
	hspectra[18]->Fill(jetkin[leadingjetindex][0]);
	if (MBDboth) {
	  hspectra[19]->Fill(jetkin[leadingjetindex][0]);
	} else {
	}
      }
	}
      
   } // end loop over events

   //==========================================================
   // Ended loop over events
   //==========================================================
   
   TCanvas *cspectra = new TCanvas("cspectra","cspectra",10,10,2500,1000);
   cspectra->Divide(3,1);
   cspectra->cd(1);
   gPad->SetLogy(1);
   gPad->SetTicks(1);

   int icolor[10] = {1,4,6,2,1,8,5,1,1,1};
   
   for (int i=0;i<20;i++) {
     hspectra[i]->SetLineColor(icolor[i]);
     hspectra[i]->SetLineWidth(3);     
     if (i==0) hspectra[0]->Draw();
     hspectra[i]->GetXaxis()->SetRangeUser(0.0,100.0);
     hspectra[i]->Draw("same");
   }

   TLegend *tleg = new TLegend(0.5,0.6,0.85,0.85,"sPHENIX Internal pp 200 GeV","brNDC");
   tleg->Draw("same");

   // ratios, ratios, ratios...
   cspectra->cd(2);
   gPad->SetTicks(1);

   TH1D *hratio[10];
   for (int i=0;i<10;i++) {
     snprintf(fooname,100,"hratio%d",i);
     hratio[i]  = new TH1D(("hratio"+to_string(i)+simstr).c_str(),fooname,nbins,binb);
     hratio[i]->SetLineColor(icolor[i]);
     hratio[i]->SetLineWidth(3);
   }
   // ratio of 7 to 6 and 9 to 8
   for (int j=1;j<=nbins;j++) {
     if (hspectra[6]->GetBinContent(j)>0)
       hratio[0]->SetBinContent(j, hspectra[7]->GetBinContent(j)/hspectra[6]->GetBinContent(j));
     if (hspectra[8]->GetBinContent(j)>0)
       hratio[1]->SetBinContent(j, hspectra[9]->GetBinContent(j)/hspectra[8]->GetBinContent(j));

     if (hspectra[16]->GetBinContent(j)>0)
       hratio[2]->SetBinContent(j, hspectra[17]->GetBinContent(j)/hspectra[16]->GetBinContent(j));
     if (hspectra[18]->GetBinContent(j)>0)
       hratio[3]->SetBinContent(j, hspectra[19]->GetBinContent(j)/hspectra[18]->GetBinContent(j));

   }

   hratio[0]->GetXaxis()->SetTitle("Jet Reco p_{T} [GeV]");
   hratio[0]->GetYaxis()->SetTitle("Ratio");
   hratio[0]->GetYaxis()->SetRangeUser(0.0,1.05);
   hratio[0]->GetXaxis()->SetRangeUser(0.0,65.0);
   hratio[0]->SetLineColor(kRed);
   hratio[0]->Draw();
   hratio[1]->SetLineColor(kBlue);
   hratio[1]->Draw("same");
   hratio[2]->SetLineColor(kRed);
   hratio[2]->SetLineStyle(9);
   hratio[3]->SetLineColor(kBlue);
   hratio[3]->SetLineStyle(9);
   hratio[2]->Draw("same");
   hratio[3]->Draw("same");
   tleg->AddEntry(hratio[0],"Dijet \& MBD / MBD","l");
   tleg->AddEntry(hratio[2],"Dijet \& MBD / MBD [w/ tight E-frac]","l");
   tleg->AddEntry(hratio[1],"Dijet \& MBD / Dijet","l");   
   tleg->AddEntry(hratio[3],"Dijet \& MBD / Dijet [w/ tight E-frac]","l");   
   tleg->Draw("same");
   
   cspectra->cd(3);
   gPad->SetTicks(1);
   hfracgood0->Draw();
   hfracbad0->SetLineColor(kRed);
   hfracbad0->Draw("same");
   
   
   cspectra->Update();

    // Define the function for the second axis (e.g., a linear scaling)
    // Here, we'll make the second axis go from 0 to 100 as the primary goes from 0 to 1.9
    // So, if primary_x = 0, secondary_x = 0. If primary_x = 1.9, secondary_x = 100.
    // secondary_x = (primary_x / 1.9) * 100
    TF1 *f_second_axis = new TF1("f_second_axis", "x/0.75", 0, 200.0);

    // Create the TGaxis for the second x-axis
    // Arguments: xmin, ymin, xmax, ymax, function_name, ndiv, options
    // We place it at the top of the pad (gPad->GetUymax())
    TGaxis *axis = new TGaxis(gPad->GetUxmin(), gPad->GetUymin()-0.2,
                              gPad->GetUxmax(), gPad->GetUymin()-0.2,
                              "f_second_axis", 510, "+L"); // "+L" option for labels on top and left
    axis->SetTitle("Secondary X-axis");
    axis->SetLabelOffset(0.01); // Adjust label offset if needed
    axis->SetTitleOffset(1.2);  // Adjust title offset if needed
    //    axis->Draw();

    TCanvas *ctime = new TCanvas("ctime","ctime",10,10,900,400);
    ctime->Divide(2,1);
    ctime->cd(1);
    gPad->SetTicks(1);
    hleadtimeNOMBD->SetXTitle("Jet Lead Time [ns]");
    hleadtimeNOMBD->Draw();    
    hleadtimeNOMBDwdijet->SetFillColorAlpha(kBlack,0.3);
    hleadtimeNOMBDwdijet->Draw("same");    
    hleadtimeNOMBDwdijetP->SetFillColorAlpha(kRed,0.6);
    hleadtimeNOMBDwdijetP->Draw("same");    
    
    ctime->cd(2);
    gPad->SetTicks(1);
    hleadtimeYESMBD->SetXTitle("Jet Lead Time - MBD Time [ns]");
    hleadtimeYESMBD->Draw();
    hleadtimeYESMBDwdijet->SetFillColorAlpha(kBlack,0.3);
    hleadtimeYESMBDwdijet->Draw("same");    
    hleadtimeYESMBDwdijetP->SetFillColorAlpha(kRed,0.6);
    hleadtimeYESMBDwdijetP->Draw("same");    


    TFile* outf = TFile::Open(("hists_mbdtimereq_out_"+simstr+".root").c_str(),"RECREATE");

    outf->cd();

    hleadtimeNOMBD->Write();
    hleadtimeYESMBD->Write();
    hleadtimeNOMBDwdijet->Write();
    hleadtimeYESMBDwdijet->Write();
    hleadtimeNOMBDwdijetP->Write();
    hleadtimeYESMBDwdijetP->Write();
    htdtmbdt->Write();
    hpttmbdt->Write();
    hptdtmbdt->Write();

    hpttmbdt_dtc->Write();
    hptdtmbdt_ltc->Write();

    hpttmbdt_dtc_mbdboth->Write();
    hptdtmbdt_ltc_mbdboth->Write();
    for(int i=0; i<20; ++i)
      {
	hspectra[i]->Write();
	if(i<10) hratio[i]->Write();
      }

    outf->Write();
    outf->Close();
    cout << npx << " " << ncx << " " << npd << endl;
    return 0;
}
