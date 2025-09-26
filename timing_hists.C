#include <cmath>

int timing_hists(int lo, int hi, int type)
{
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  int jet_n, runnum, evtnum, failscut, tjet_n;
  float frcem[100];
  float frcoh[100];
  float jet_e[100];
  float jet_pt[100];
  float tjet_pt[100];
  float jet_eta[100];
  float jet_phi[100];
  float jet_t[100];
  float jet_t_em[100];
  float jet_t_oh[100];
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

  float scalefactor = 1;
  if(type==1) scalefactor = 7.2695e-9;
  if(type==2) scalefactor = 1.034e-11;
  if(type==3) scalefactor = 2.502e-6;
  TChain* jet_tree = new TChain("jet_tree");
  TChain* wft = new TChain("wft");
  
  string tempinfilename;
  string nametype = type?"sim":"dat";
  string samplename = "dat";
  string listname = "chi2files";
  if(type==1) samplename = "50";
  if(type==2) samplename = "70";
  if(type==3) samplename = "30";
  listname += samplename+".txt";
  
  ifstream chi2list(listname);

  for(int i=0; i<lo; ++i)
    {
      std::getline(chi2list,tempinfilename);
    }
  int counter = lo;
  cout << hi << endl;
  while(std::getline(chi2list,tempinfilename))
    {
      //cout << "Read: " << tempinfilename << endl;
      jet_tree->Add(tempinfilename.c_str());
      ++counter;
      if(counter >= hi)
	{
	  cout << counter << endl;
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
      cout << "Read: " << tempinfilename << endl;
      wft->Add(tempinfilename.c_str());
      ++wfcounter;
      if(wfcounter >= hi) break;
    }
  */
  //TTree* jet_tree = (TTree*)evtfile->Get("jet_tree");
  jet_tree->SetBranchStatus("*",0);
  jet_tree->SetBranchStatus("jet_n",1);
  jet_tree->SetBranchStatus("failscut",1);
  jet_tree->SetBranchStatus("alljetfrcem",1);
  jet_tree->SetBranchStatus("alljetfrcoh",1);
  jet_tree->SetBranchStatus("jet_et",1);
  jet_tree->SetBranchStatus("jet_pt",1);
  //jet_tree->SetBranchStatus("jet_eta",1);
  //jet_tree->SetBranchStatus("jet_phi",1);
  //jet_tree->SetBranchStatus("zvtx",1);
  jet_tree->SetBranchStatus("jet_t",1);
  jet_tree->SetBranchStatus("jet_t_em",1);
  jet_tree->SetBranchStatus("jet_t_oh",1);
  //jet_tree->SetBranchStatus("emtow",1);
  //jet_tree->SetBranchStatus("ihtow",1);
  //jet_tree->SetBranchStatus("ohtow",1);
  if(type != 0)
    {
      jet_tree->SetBranchStatus("tjet_pt",1);
      jet_tree->SetBranchAddress("tjet_pt",tjet_pt);
      jet_tree->SetBranchStatus("tjet_n",1);
      jet_tree->SetBranchAddress("tjet_n",&tjet_n);
    }
  jet_tree->SetBranchAddress("jet_n",&jet_n);
  jet_tree->SetBranchAddress("failscut",&failscut);
  jet_tree->SetBranchAddress("alljetfrcem",frcem);
  jet_tree->SetBranchAddress("alljetfrcoh",frcoh);
  jet_tree->SetBranchAddress("jet_et",jet_e);
  jet_tree->SetBranchAddress("jet_pt",jet_pt);
  //jet_tree->SetBranchAddress("jet_eta",jet_eta);
  //jet_tree->SetBranchAddress("jet_phi",jet_phi);
  //jet_tree->SetBranchAddress("zvtx",&zvtx);
  jet_tree->SetBranchAddress("jet_t",jet_t);
  jet_tree->SetBranchAddress("jet_t_em",jet_t_em);
  jet_tree->SetBranchAddress("jet_t_oh",jet_t_oh);
  //jet_tree->SetBranchAddress("emtow",emtow);
  //jet_tree->SetBranchAddress("ihtow",ihtow);
  //jet_tree->SetBranchAddress("ohtow",ohtow);
  
  const int ncut = 4;
  TH3D* h3_pt_dt_t[ncut];
  TH3D* h3_pt_dtem_tem[ncut];
  TH3D* h3_pt_dtoh_toh[ncut];
  TH3D* h3_apt_dt_t[ncut];
  TH3D* h3_pt_dt_lt[ncut];
  TH3D* h3_pt_adt_t[ncut];
  TH2D* h2_pt_t[ncut];
  TH2D* h2_lt_st[ncut];
  string cuttype[ncut] = {"all","dijet","frac","both"};
  for(int i=0; i<ncut; ++i)
    {
      h3_pt_dt_t[i] = new TH3D((nametype+"h3_pt_dt_t_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];#Delta t_{l,sl} [ns];t_{jet} [ns]",200,0,200,200,-25,25,200,-25,25);
      h3_pt_dtem_tem[i] = new TH3D((nametype+"h3_pt_dtem_tem_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];#Delta t_{EM} [ns];t_{jet,EM} [ns]",200,0,200,200,-25,25,200,-25,25);
      h3_pt_dtoh_toh[i] = new TH3D((nametype+"h3_pt_dtoh_toh_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];#Delta t_{OH} [ns];t_{jet,OH} [ns]",200,0,200,200,-25,25,200,-25,25);
      h3_apt_dt_t[i] = new TH3D((nametype+"h3_apt_dt_t_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];#Delta t_{l,sl} [ns];t_{jet} [ns]",200,0,200,200,-25,25,200,-25,25);
      h3_pt_dt_lt[i] = new TH3D((nametype+"h3_pt_dt_lt_"+cuttype[i]).c_str(),";p_{T}^{lead} [GeV];#Delta t_{l,sl} [ns];t_{lead} [ns]",200,0,200,200,-25,25,200,-25,25);
      h3_pt_adt_t[i] = new TH3D((nametype+"h3_pt_adt_t_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];#Delta t_{l,sl} [ns];t_{jet} [ns]",200,0,200,50,0,5,200,-25,25);
      h2_pt_t[i] = new TH2D((nametype+"h2_pt_t_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];t_{jet} [ns]",200,0,200,200,-25,25);
      h2_lt_st[i] = new TH2D((nametype+"h2_lt_st_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];t_{jet} [ns]",200,-25,25,200,-25,25);
    }
  
  int wfte = 0;
  for(int i=0; i<jet_tree->GetEntries(); ++i)
    {
      if(!(i%100)) cout << i << endl;
      jet_tree->GetEntry(i);
      /*
      float calosum = 0;
      for(int j=0; j<96; ++j)
	{
	  for(int k=0; k<256; ++k)
	    {
	      calosum += emtow[j][k];
	      if(j<24 && k<64)
		{
		  calosum += ihtow[j][k];
		  calosum += ohtow[j][k];
		}
	    }
	}
      if(calosum > 180) continue;
      */
      /*
      wft->GetEntry(wfte);
      int flag = 0;
      
      cout << rnwf << " " << runnum << " " << enwf << " " << evtnum << endl;
      
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
      int lindex = -1;
      int sindex = -1;
      float lpt = 0;
      float spt = 0;
      float lt = 10;
      float st = -10;
      float sumjetpt = 0;
      for(int j=0; j<jet_n; ++j)
	{
	  sumjetpt += jet_pt[j];
	  if(jet_pt[j] > lpt)
	    {
	      st = lt;
	      spt = lpt;
	      sindex = lindex;
	      lt = jet_t[j];
	      lpt = jet_pt[j];
	      lindex = j;
	    }
	  else if(jet_pt[j] > spt)
	    {
	      spt = jet_pt[j];
	      st = jet_t[j];
	      sindex = j;
	    }
	}
      if(type != 0)
	{
	  float ltpt = 0;
	  for(int j=0; j<tjet_n; ++j)
	    {
	      if(tjet_pt[j] > ltpt) ltpt = tjet_pt[j];
	    }
	  if(type == 1 && (ltpt > 71 || ltpt < 52)) continue;
	  if(type == 2 && ltpt < 71) continue;
	  if(type == 3 && ltpt > 52) continue;
	}
      //if(sumjetpt > 175) continue;
      float stons = 17.6;
      lt*=stons;
      st*=stons;
      float ltem = 0;
      if(frcem[lindex]*jet_e[lindex] > 5) ltem = jet_t_em[lindex]*stons;
      float stem = 0;
      if(frcem[sindex]*jet_e[sindex] > 5) stem = jet_t_em[sindex]*stons;
      float ltoh = 0;
      if(frcoh[lindex]*jet_e[lindex] > 5) ltoh = jet_t_oh[lindex]*stons;
      float stoh = 0;
      if(frcoh[sindex]*jet_e[sindex] > 5) stoh = jet_t_oh[sindex]*stons;
      if(sindex > -1)
	{
	  h2_lt_st[0]->Fill(lt,st,scalefactor);
	  if(!std::isnan(ltem) && !std::isnan(stem) && ltem != 0 && stem != 0)
	    {
	      h3_pt_dtem_tem[0]->Fill(lpt,ltem-stem,ltem,scalefactor);
	      h3_pt_dtem_tem[0]->Fill(spt,ltem-stem,stem,scalefactor);
	    }
	  if(!std::isnan(ltoh) && !std::isnan(stoh) && ltoh != 0 && stoh != 0)
	    {
	      h3_pt_dtoh_toh[0]->Fill(lpt,ltoh-stoh,ltoh,scalefactor);
	      h3_pt_dtoh_toh[0]->Fill(spt,ltoh-stoh,stoh,scalefactor);
	    }
	  if(failscut==0)
	    {
	      if(!std::isnan(ltem) && !std::isnan(stem) && ltem != 0 && stem != 0)
		{
		  h3_pt_dtem_tem[1]->Fill(lpt,ltem-stem,ltem,scalefactor);
		  h3_pt_dtem_tem[1]->Fill(spt,ltem-stem,stem,scalefactor);
		}
	      if(!std::isnan(ltoh) && !std::isnan(stoh) && ltoh != 0 && stoh != 0)
		{
		  h3_pt_dtoh_toh[1]->Fill(lpt,ltoh-stoh,ltoh,scalefactor);
		  h3_pt_dtoh_toh[1]->Fill(spt,ltoh-stoh,stoh,scalefactor);
		}
	      h2_lt_st[1]->Fill(lt,st,scalefactor);
	    }
	  if(failscut==1)
	    {
	      if(!std::isnan(ltem) && !std::isnan(stem) && ltem != 0 && stem != 0)
		{
		  h3_pt_dtem_tem[2]->Fill(lpt,ltem-stem,ltem,scalefactor);
		  h3_pt_dtem_tem[2]->Fill(spt,ltem-stem,stem,scalefactor);
		}
	      if(!std::isnan(ltoh) && !std::isnan(stoh) && ltoh != 0 && stoh != 0)
		{
		  h3_pt_dtoh_toh[2]->Fill(lpt,ltoh-stoh,ltoh,scalefactor);
		  h3_pt_dtoh_toh[2]->Fill(spt,ltoh-stoh,stoh,scalefactor);
		}
	      h2_lt_st[2]->Fill(lt,st,scalefactor);
	    }
	  if(failscut==2)
	    {
	      for(int k=1;k<4;++k)
		{
		  if(!std::isnan(ltem) && !std::isnan(stem) && ltem != 0 && stem != 0)
		    {
		      h3_pt_dtem_tem[k]->Fill(lpt,ltem-stem,ltem,scalefactor);
		      h3_pt_dtem_tem[k]->Fill(spt,ltem-stem,stem,scalefactor);
		    }
		  if(!std::isnan(ltoh) && !std::isnan(stoh) && ltoh != 0 && stoh != 0)
		    {
		      h3_pt_dtoh_toh[k]->Fill(lpt,ltoh-stoh,ltoh,scalefactor);
		      h3_pt_dtoh_toh[k]->Fill(spt,ltoh-stoh,stoh,scalefactor);
		    }
		  h2_lt_st[k]->Fill(lt,st,scalefactor);
		}
	    }
	}

      if(sindex > -1 && st != 0)
	{
	  h3_pt_dt_t[0]->Fill(lpt,lt-st,lt,scalefactor);
	  h3_pt_dt_t[0]->Fill(spt,lt-st,st,scalefactor);
	  h3_pt_dt_lt[0]->Fill(lpt,lt-st,lt,scalefactor);
	  h3_pt_adt_t[0]->Fill(lpt,abs(lt-st),lt,scalefactor);
	  h3_pt_adt_t[0]->Fill(spt,abs(lt-st),st,scalefactor);
	  
	  if(failscut==0)
	    {
	      h3_pt_dt_t[1]->Fill(lpt,lt-st,lt,scalefactor);
	      h3_pt_dt_t[1]->Fill(spt,lt-st,st,scalefactor);
	      h3_pt_dt_lt[1]->Fill(lpt,lt-st,lt,scalefactor);
	      h3_pt_adt_t[1]->Fill(lpt,abs(lt-st),lt,scalefactor);
	      h3_pt_adt_t[1]->Fill(spt,abs(lt-st),st,scalefactor);
	    }
	  if(failscut==1)
	    {
	      h3_pt_dt_t[2]->Fill(lpt,lt-st,lt,scalefactor);
	      h3_pt_dt_t[2]->Fill(spt,lt-st,st,scalefactor);
	      h3_pt_dt_lt[2]->Fill(lpt,lt-st,lt,scalefactor);
	      h3_pt_adt_t[2]->Fill(lpt,abs(lt-st),lt,scalefactor);
	      h3_pt_adt_t[2]->Fill(spt,abs(lt-st),st,scalefactor);
	    }
	  if(failscut==2)
	    {
	      for(int k=1;k<4;++k)
		{
		  h3_pt_dt_t[k]->Fill(lpt,lt-st,lt,scalefactor);
		  h3_pt_dt_t[k]->Fill(spt,lt-st,st,scalefactor);
		  h3_pt_dt_lt[k]->Fill(lpt,lt-st,lt,scalefactor);
		  h3_pt_adt_t[k]->Fill(lpt,abs(lt-st),lt,scalefactor);
		  h3_pt_adt_t[k]->Fill(spt,abs(lt-st),st,scalefactor);
		}
	    }
	}
      
      
      for(int j=0; j<jet_n; ++j)
	{
	  if(jet_t[j] == 0) continue;
	  float tns = jet_t[j]*stons;
	  if(sindex > -1 && st != 0)
	    {
	      h3_apt_dt_t[0]->Fill(jet_pt[j],lt-st,tns,scalefactor);
	      if(failscut == 0)
		{
		  h3_apt_dt_t[1]->Fill(jet_pt[j],lt-st,tns,scalefactor);
		}
	      if(failscut == 1)
		{
		  h3_apt_dt_t[2]->Fill(jet_pt[j],lt-st,tns,scalefactor);
		}
	      if(failscut == 2)
		{
		  for(int k=1; k<4; ++k) h3_apt_dt_t[k]->Fill(jet_pt[j],lt-st,tns,scalefactor);
		}
	    }
	  h2_pt_t[0]->Fill(jet_pt[j],tns,scalefactor);
	  if(failscut==0) h2_pt_t[1]->Fill(jet_pt[j],tns,scalefactor);
	  if(failscut==1) h2_pt_t[2]->Fill(jet_pt[j],tns,scalefactor);
	  if(failscut==2)
	    {
	      for(int k=1;k<4;++k) h2_pt_t[k]->Fill(jet_pt[j],tns,scalefactor);
	    }
	}
    }

  TFile* outf = TFile::Open(("../timing/timinghists"+to_string(lo)+"_"+to_string(hi)+"_"+samplename+".root").c_str(),"RECREATE");
  outf->cd();
  for(int i=0; i<ncut; ++i)
    {
      h3_pt_adt_t[i]->Write();
      h3_pt_dt_t[i]->Write();
      h3_pt_dt_lt[i]->Write();
      h3_pt_dtem_tem[i]->Write();
      h3_pt_dtoh_toh[i]->Write();
      h2_pt_t[i]->Write();
      h2_lt_st[i]->Write();
    }

  outf->Write();

  return 0;
      
}
      
