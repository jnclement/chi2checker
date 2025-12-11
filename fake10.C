void print_bounds_and_efficiency(string cutname, float* bounds, long long unsigned int* numerator, long long unsigned int* denominator)
{
  float eff = ((float)numerator[0])/denominator[0];
  cout << "Efficiency of cut " << cutname << " in range " << bounds[0] << "-" << bounds[3] << " GeV: " << eff << endl;
  for(int i=1; i<4; ++i)
    {
      eff = ((float)numerator[i])/denominator[i];
      cout << "Efficiency of cut " << cutname << " in range " << bounds[i-1] << "-" << bounds[i] << " GeV: " << eff << endl;
    }
}


int fake10(string filename, bool issim = false, bool dofrac = false, string simstr = "dat", bool samtime = false, bool slew = false) {

  if(dofrac) simstr += "frac";
  gErrorIgnoreLevel = kError;
  gROOT->SetStyle("Plain");
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   gStyle->SetOptDate(111);
   std::map<int, float> rntmap = {};
   TF1* slewcor;
   TFile* slewfile;
   if(slew)
     {
       slewfile = TFile::Open(("slewcorfit_t_"+simstr+".root").c_str(),"READ");
       slewcor = (TF1*)slewfile->Get("pol2");
     }

   ifstream file = ifstream("/sphenix/user/samfred/projects/mbdt0/histmaking/MbdPmt.corr");
   string line;
   while (getline(file,line)) {
     int irunnum;
     float t0;
     istringstream iss(line);
     iss >> irunnum >> t0;
     rntmap[irunnum] = t0;
   }
   TFile *f  = TFile::Open(filename.c_str(),"READ");
   TTree *outt = static_cast <TTree *> (f->Get("outt"));

   Int_t           passfrac;
   Int_t           passdijet;
   Int_t           evt;
   Int_t           run;
   Int_t           njet;
   Float_t         jetkin[60][5];
   Float_t         zvtx;
   UInt_t           mbdhit[2];
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
   outt->SetBranchAddress("mbdhit",mbdhit);
   outt->SetBranchAddress("avgt",avgt);
   outt->SetBranchAddress("frac",frac);   

   //=======================
   // DECLARE HISTOGRAMS
   //=======================
   
   int nbins = 35;

   float binb[36];
   for (int i=0;i<25;i++) binb[i] = 0.0 + 1.0 * (float) i;
   for (int i=0;i<7;i++) binb[i+25] = 25.0 + 2.5 * (float) i;
   for (int i=0;i<4;i++) binb[i+26+6] = 45.0 + 5.0 * (float) i;      

   //for (int i=0;i<36;i++) cout << binb[i] << endl;
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
   TH1D *hleadtimeNOMBD = new TH1D(("hleadtimeNOMBD"+simstr).c_str(),"hleadtimeNOMBD",600,-30.0,30.0);
   TH1D *hleadtimeYESMBD = new TH1D(("hleadtimeYESMBD"+simstr).c_str(),"hleadtimeYESMBD",600,-30.0,30.0);   
   TH1D *hleadtimeNOMBDwdijet = new TH1D(("hleadtimeNOMBDwdijet"+simstr).c_str(),"hleadtimeNOMBDwdijet",600,-30.0,30.0);
   TH1D *hleadtimeYESMBDwdijet = new TH1D(("hleadtimeYESMBDwdijet"+simstr).c_str(),"hleadtimeYESMBDwdijet",600,-30.0,30.0);   
   TH1D *hleadtimeNOMBDwdijetP = new TH1D(("hleadtimeNOMBDwdijetP"+simstr).c_str(),"hleadtimeNOMBDwdijetP",600,-30.0,30.0);
   TH1D *hleadtimeYESMBDwdijetP = new TH1D(("hleadtimeYESMBDwdijetP"+simstr).c_str(),"hleadtimeYESMBDwdijetP",600,-30.0,30.0);
   TH3D* htdtmbdt = new TH3D(("htdtmbdt"+simstr).c_str(),";Jet t [ns];#Delta t [ns]; MBD - t_{lead} [ns]",300,-30,30,300,-30,30,300,-30,30);

   TH3D* hptdtfrac = new TH3D(("hptdtfrac"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet} [GeV];#Delta-t [ns];Lead Jet Energy Fraction in OHCal",100,0,100,200,-10,10,120,-0.1,1.1);
   TH3D* hpttfrac = new TH3D(("hpttfrac"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet} [GeV];t_{lead} [ns];Lead Jet Energy Fraction in OHCal",100,0,100,200,-10,10,120,-0.1,1.1);
   TH3D* hptdtemfrac = new TH3D(("hptdtemfrac"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet};#Delta-t [ns];Lead Jet Energy Fraction in EMCal",100,0,100,200,-10,10,120,-0.1,1.1);
   TH3D* hptdtemfracnot = new TH3D(("hptdtemfracnot"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet};#Delta-t [ns];Lead Jet Energy Fraction in EMCal",100,0,100,200,-10,10,120,-0.1,1.1);
   TH3D* hptdtohfrac = new TH3D(("hptdtohfrac"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet};#Delta-t [ns];Lead Jet Energy Fraction in OHCal",100,0,100,200,-10,10,120,-0.1,1.1);

   TH3D* hptdtemfrac_all = new TH3D(("hptdtemfrac_all"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet};#Delta-t [ns];Lead Jet Energy Fraction in EMCal",100,0,100,200,-10,10,120,-0.1,1.1);
   TH3D* hptdtohfrac_all = new TH3D(("hptdtohfrac_all"+simstr).c_str(),";Uncalibrated p_{T}^{leadjet};#Delta-t [ns];Lead Jet Energy Fraction in OHCal",100,0,100,200,-10,10,120,-0.1,1.1);
   
   TH3D* hpttdt = new TH3D(("hpttdt"+simstr).c_str(),";Uncalibrated p_{T}^{jet} [GeV];Jet t [ns];#Delta t [ns]",100,0,100,600,-30,30,600,-30,30);
   TH3D* hpttmbdt = new TH3D(("hpttmbdt"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt = new TH3D(("hptdtmbdt"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt2 = new TH3D(("hptdtmbdt2"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD t (offset corrected) [ns]",60,0,60,600,-30,30,600,-30,30);

   TH3D* hpttmbdt_dtc = new TH3D(("hpttmbdt_dtc"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt_ltc = new TH3D(("hptdtmbdt_ltc"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);

   TH3D* hpttmbdt_dtc_mbdboth = new TH3D(("hpttmbdt_dtc_mbdboth"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt_ltc_mbdboth = new TH3D(("hptdtmbdt_ltc_mbdboth"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt_ltc_mbdboth2 = new TH3D(("hptdtmbdt_ltc_mbdboth2"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD t (offset corrected) [ns]",60,0,60,600,-30,30,600,-30,30);

   TH3D* hpttmbdt_dtc_mbdeither = new TH3D(("hpttmbdt_dtc_mbdeither"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];t_{lead} [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt_ltc_mbdeither = new TH3D(("hptdtmbdt_ltc_mbdeither"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD - t_{lead} [ns]",60,0,60,600,-30,30,600,-30,30);
   TH3D* hptdtmbdt_ltc_mbdeither2 = new TH3D(("hptdtmbdt_ltc_mbdeither2"+simstr).c_str(),";Uncalibrated p_{T}^{lead} [GeV];#Delta t [ns];MBD t (offset corrected) [ns]",60,0,60,600,-30,30,600,-30,30);
   
   TH1D *hfracgood0 = new TH1D(("hfracgood0"+simstr).c_str(),"hfracgood0",6576,47288.5,53864.5);
   TH1D *hfracbad0 = new TH1D(("hfracbad0"+simstr).c_str(),"hfracbad0",6576,47288.5,53864.5);

   TH2D *hrnt = new TH2D(("hrnt"+simstr).c_str(),";Run Number;MBD Time (No Offset Correction) [ns]",6576,47288.5,53864.5,600,-30,30);
   TH2D *hrnto = new TH2D(("hrnto"+simstr).c_str(),";Run Number;MBD Time (Offset Corrected) [ns]",6576,47288.5,53864.5,600,-30,30);
   TH2D *hrntj = new TH2D(("hrntj"+simstr).c_str(),";Run Number;MBD - t_{lead} [ns]",6576,47288.5,53864.5,600,-30,30);
   
   //=========================================================
   // LOOP OVER EVENTS
   //=========================================================


   long long unsigned int nDijetPartner[4] = {0};
   long long unsigned int nPassLeadTime[4] = {0};
   long long unsigned int nPassDeltat[4] = {0};
   long long unsigned int nPassBothTime[4] = {0};
   float bounds[4] = {30,35,40,45};
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
      if(rntmap.count(run) && samtime) mbdoffset = -rntmap[run];
      else if(samtime && !rntmap.count(run)) continue;

      if (avgt[0] > -99 && avgt[1] < -99) mbdtime = avgt[0];
      if (avgt[0] < -99 && avgt[1] > -99) mbdtime = avgt[1];
      if (avgt[0] > -99 || avgt[1] > -99) MBDeither = true;
      if (zvtx == 0 || abs(zvtx) > 999 || isnan(zvtx)) MBDboth = false;
      if (avgt[0] > -99 && avgt[1] > -99) {
	mbdtime = 0.5*(avgt[0]+avgt[1]);
	//	MBDboth = true;
      }

      int subleadingjetindex = 0;
      for (int ijet=0;ijet<njet;ijet++) {
	double jetpt = jetkin[ijet][0];
	if (jetpt < jetkin[leadingjetindex][0] && jetpt >= jetkin[subleadingjetindex][0]) {
	  subleadingjetindex = ijet;
	}
      }

      float leadslew = slew?slewcor->Eval(frac[leadingjetindex][1]):0;
      float subslew = slew?slewcor->Eval(frac[subleadingjetindex][1]):0;

      float totalslew = leadslew-subslew;

      jetleadtimeMBD = mbdoffset + mbdtime - 17.6*jetkin[leadingjetindex][4] + leadslew;
      jetleadtime -= leadslew;
      bool PassLeadTimeWIDE = true;
      bool inTailLeadTime = false;
      if (TMath::Abs(jetleadtimeMBD) > 3.0 && !issim /*&& mbdtime > -99*/) PassLeadTime = false;
      if (jetleadtimeMBD < -3.0 && jetleadtimeMBD > -30.0) inTailLeadTime = true;
      if (TMath::Abs(jetleadtimeMBD) > 10.0 && !issim /*&& mbdtime > -99*/) PassLeadTimeWIDE = false;
      hleadtimeYESMBD->Fill(jetleadtimeMBD);


      

      // jet pT spectra - no cuts, with time requirement, with MBD requirements
      //hspectra[0]->Fill(jetkin[leadingjetindex][0]);   //  raw
      //if (PassLeadTime) hspectra[1]->Fill(jetkin[leadingjetindex][0]);   //  lead time cut
      //      if (PassLeadTime && MBDeither) hspectra[2]->Fill(jetkin[leadingjetindex][0]);   //  lead time cut
      //if (PassLeadTime && MBDboth) hspectra[3]->Fill(jetkin[leadingjetindex][0]);   //  lead time cut            

      // now is there a dijet partner (what requirements do we want to apply....)

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
      double dijetTimediff = 17.6*jetkin[leadingjetindex][4] - 17.6*jetkin[subleadingjetindex][4] - totalslew;
      if(!issim)
	{

	  if (TMath::Abs(dijetTimediff) > 3.0) DijetAndt = false;
	}

      
      if(DijetPartner && (PassMinimalFrac || !dofrac))
	{
	  for(int j=0; j<njet; ++j)
	    {
	      if(jetkin[j][0] > bounds[0] && jetkin[j][0] < bounds[3])
		{
		  ++nDijetPartner[0];
		  if(PassLeadTimeONLY)
		    {
		      ++nPassLeadTime[0];
		    }
		  if(DijetAndt)
		    {
		      ++nPassDeltat[0];
		    }
		  if(DijetAndt && PassLeadTimeONLY)
		    {
		      ++nPassBothTime[0];
		    }
		  for(int k=1; k<4; ++k)
		    {
		      if(jetkin[j][0] < bounds[k])
			{
			  ++nDijetPartner[k];
			  if(PassLeadTimeONLY)
			    {
			      ++nPassLeadTime[k];
			    }
			  if(DijetAndt)
			    {
			      ++nPassDeltat[k];
			    }
			  if(DijetAndt && PassLeadTimeONLY)
			    {
			      ++nPassBothTime[k];
			    }
			  break;
			}
		    }
		}
	    }
	}

      if(DijetPartner)
	{
	  hptdtfrac->Fill(jetkin[leadingjetindex][0],dijetTimediff,frac[leadingjetindex][1]);
	  hpttfrac->Fill(jetkin[leadingjetindex][0],jetleadtime,frac[leadingjetindex][1]);
	  hptdtemfracnot->Fill(jetkin[leadingjetindex][0],dijetTimediff,frac[leadingjetindex][0]);
	  if(jetleadtime < 4 && jetleadtime > -8)
	    {
	      hptdtemfrac->Fill(jetkin[leadingjetindex][0],dijetTimediff,frac[leadingjetindex][0]);
	      hptdtohfrac->Fill(jetkin[leadingjetindex][0],dijetTimediff,frac[leadingjetindex][1]);
	    }
	  for(int j=0; j<njet; ++j)
	    {
	      if(jetleadtime < 4 && jetleadtime > -8)
		{
		  hptdtemfrac_all->Fill(jetkin[j][0],dijetTimediff,frac[j][0]);
		  hptdtohfrac_all->Fill(jetkin[j][0],dijetTimediff,frac[j][1]);
		}
	    }
	}
      
      if(!dofrac || (dofrac && PassMinimalFrac))
	{

      
      htdtmbdt->Fill(jetleadtime,dijetTimediff,jetleadtimeMBD);
      hpttmbdt->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetAndt) hpttmbdt_dtc->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetAndt && MBDboth) hpttmbdt_dtc_mbdboth->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetAndt && !MBDboth && MBDeither) hpttmbdt_dtc_mbdeither->Fill(jetkin[leadingjetindex][0],jetleadtime,jetleadtimeMBD);
      if(DijetPartner)
	{
	  hptdtmbdt->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	  if(PassLeadTimeONLY)
	    {
	      hptdtmbdt_ltc->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	      hptdtmbdt2->Fill(jetkin[leadingjetindex][0],dijetTimediff,mbdtime+mbdoffset);
	    }
	  if(PassLeadTimeONLY && MBDboth)
	    {
	      hptdtmbdt_ltc_mbdboth->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	      hptdtmbdt_ltc_mbdboth2->Fill(jetkin[leadingjetindex][0],dijetTimediff,mbdtime+mbdoffset);
	    }
	  if(PassLeadTimeONLY && !MBDboth && MBDeither)
	    {
	      hptdtmbdt_ltc_mbdeither->Fill(jetkin[leadingjetindex][0],dijetTimediff,jetleadtimeMBD);
	      hptdtmbdt_ltc_mbdeither2->Fill(jetkin[leadingjetindex][0],dijetTimediff,mbdtime+mbdoffset);
	    }
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


	  if(DijetPartner)
	    {
	      for(int j=0; j<njet; ++j)
		{
		  hpttdt->Fill(jetkin[j][0],jetleadtime,dijetTimediff);
		}
	    }

      if (MBDboth && PassLeadTime && DijetPartner && PassMinimalFrac) hspectra[15]->Fill(jetkin[leadingjetindex][0]);
      
      if (MBDeither && !MBDboth && PassLeadTimeONLY && DijetAndt && PassLeadTimeWIDE) hspectra[5]->Fill(jetkin[leadingjetindex][0]);
      if (MBDeither && !MBDboth && PassLeadTimeONLY && DijetAndt && PassLeadTime) hspectra[1]->Fill(jetkin[leadingjetindex][0]);
      if (MBDeither && !MBDboth && PassLeadTimeONLY && DijetAndt) hspectra[0]->Fill(jetkin[leadingjetindex][0]);
      
      if (MBDboth && PassLeadTimeWIDE && DijetAndt && PassLeadTimeONLY) hspectra[4]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTime && DijetAndt && PassLeadTimeONLY) hspectra[3]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTimeONLY && DijetAndt) hspectra[2]->Fill(jetkin[leadingjetindex][0]);



      if (MBDboth && PassLeadTime && PassLeadTimeONLY && DijetPartner) hspectra[7]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTime && PassLeadTimeONLY) hspectra[6]->Fill(jetkin[leadingjetindex][0]);      

      if (MBDboth && PassLeadTime && PassMinimalFrac) hspectra[16]->Fill(jetkin[leadingjetindex][0]);
      if (MBDboth && PassLeadTime && PassLeadTimeONLY && PassMinimalFrac && DijetAndt) hspectra[17]->Fill(jetkin[leadingjetindex][0]);      
      
      if (DijetAndt && PassLeadTimeONLY) {
	hspectra[8]->Fill(jetkin[leadingjetindex][0]);
	hfracgood0->Fill(run);
	hrnt->Fill(run,mbdtime);
	hrnto->Fill(run,mbdtime+mbdoffset);
	hrntj->Fill(run,jetleadtimeMBD);
	if(inTailLeadTime) hfracbad0->Fill(run);
	if (MBDboth) {
	  hspectra[9]->Fill(jetkin[leadingjetindex][0]);
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

   hratio[0]->GetXaxis()->SetTitle("Jet Reco Uncalibrated p_{T} [GeV]");
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


    TFile* outf = TFile::Open(("hists_mbdtimereq_out_"+simstr+(samtime?"sam":"")+(slew?"_slewed":"")+".root").c_str(),"RECREATE");

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
    hptdtmbdt2->Write();

    hptdtfrac->Write();
    hpttfrac->Write();
    hptdtemfrac->Write();
    hptdtohfrac->Write();

    hptdtemfrac_all->Write();
    hptdtohfrac_all->Write();
    
    hpttdt->Write();
    hpttmbdt_dtc->Write();
    hptdtmbdt_ltc->Write();

    hpttmbdt_dtc_mbdboth->Write();
    hptdtmbdt_ltc_mbdboth->Write();
    hptdtmbdt_ltc_mbdboth2->Write();
    hpttmbdt_dtc_mbdeither->Write();
    hptdtmbdt_ltc_mbdeither->Write();
    hptdtmbdt_ltc_mbdeither2->Write();
    hfracgood0->Write();
    hfracbad0->Write();
    hrnt->Write();
    hrnto->Write();
    hrntj->Write();
    hptdtemfracnot->Write();
    for(int i=0; i<20; ++i)
      {
	hspectra[i]->Write();
	if(i<10) hratio[i]->Write();
      }

    outf->Write();
    outf->Close();

    print_bounds_and_efficiency("lead time",bounds,nPassLeadTime,nDijetPartner);
    print_bounds_and_efficiency("delta t",bounds,nPassDeltat,nDijetPartner);
    print_bounds_and_efficiency("lead time && delta t",bounds,nPassBothTime,nDijetPartner);
    
    //cout << npx << " " << ncx << " " << npd << endl;
    return 0;
}
