#include <cmath>

int timing_hists(int lo, int hi)
{
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
  TChain* wft = new TChain("wft");
  
  string tempinfilename;
  ifstream chi2list("chi2files.txt");

  for(int i=0; i<lo; ++i)
    {
      std::getline(chi2list,tempinfilename);
    }
  int counter = lo;
  cout << hi << endl;
  while(std::getline(chi2list,tempinfilename))
    {
      cout << "Read: " << tempinfilename << endl;
      jet_tree->Add(tempinfilename.c_str());
      ++counter;
      if(counter >= hi)
	{
	  cout << counter << endl;
	  break;
	}
    }


  
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
  
  //TTree* jet_tree = (TTree*)evtfile->Get("jet_tree");
  jet_tree->SetBranchStatus("*",0);
  jet_tree->SetBranchStatus("jet_n",1);
  jet_tree->SetBranchStatus("failscut",1);
  jet_tree->SetBranchStatus("alljetfrcem",1);
  jet_tree->SetBranchStatus("alljetfrcoh",1);
  jet_tree->SetBranchStatus("jet_et",1);
  jet_tree->SetBranchStatus("jet_pt",1);
  jet_tree->SetBranchStatus("jet_eta",1);
  jet_tree->SetBranchStatus("jet_phi",1);
  jet_tree->SetBranchStatus("zvtx",1);
  jet_tree->SetBranchStatus("jet_t",1);
  
  jet_tree->SetBranchAddress("jet_n",&jet_n);
  jet_tree->SetBranchAddress("failscut",&failscut);
  jet_tree->SetBranchAddress("alljetfrcem",frcem);
  jet_tree->SetBranchAddress("alljetfrcoh",frcoh);
  jet_tree->SetBranchAddress("jet_et",jet_e);
  jet_tree->SetBranchAddress("jet_pt",jet_pt);
  jet_tree->SetBranchAddress("jet_eta",jet_eta);
  jet_tree->SetBranchAddress("jet_phi",jet_phi);
  jet_tree->SetBranchAddress("zvtx",&zvtx);
  jet_tree->SetBranchAddress("jet_t",jet_t);
  
  const int ncut = 4;
  TH3D* h3_pt_dt_t[ncut];
  TH2D* h2_pt_t[ncut];
  string cuttype[ncut] = {"all","dijet","frac","both"};
  for(int i=0; i<ncut; ++i)
    {
      h3_pt_dt_t[i] = new TH3D(("h3_pt_dt_t_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];#Delta t_{l,sl} [samples];t_{jet} [samples]",200,0,200,100,-5,5,100,-5,5);
      h2_pt_t[i] = new TH2D(("h2_pt_t_"+cuttype[i]).c_str(),";p_{T}^{jet} [GeV];t_{jet} [samples]",200,0,200,100,-5,5);
    }
  
  int wfte = 0;
  for(int i=0; i<jet_tree->GetEntries(); ++i)
    {
      if(!(i%100)) cout << i << endl;
      jet_tree->GetEntry(i);
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
      float lt = 0;
      float st = 0;
      for(int j=0; j<jet_n; ++j)
	{
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

      if(sindex > -1)
	{
	  h3_pt_dt_t[0]->Fill(lpt,lt-st,lt);
	  if(failscut==0) h3_pt_dt_t[1]->Fill(lpt,lt-st,lt);
	  if(failscut==1) h3_pt_dt_t[2]->Fill(lpt,lt-st,lt);
	  if(failscut==2)
	    {
	      for(int j=1;j<4;++j) h3_pt_dt_t[j]->Fill(lpt,lt-st,lt);
	    }
	}
      for(int j=0; j<jet_n; ++j)
	{
	  h2_pt_t[0]->Fill(jet_pt[j],jet_t[j]);
	  if(failscut==0) h2_pt_t[1]->Fill(jet_pt[j],jet_t[j]);
	  if(failscut==1) h2_pt_t[2]->Fill(jet_pt[j],jet_t[j]);
	  if(failscut==2)
	    {
	      for(int k=1;k<4;++k) h2_pt_t[k]->Fill(jet_pt[j],jet_t[j]);
	    }
	}
    }

  TFile* outf = TFile::Open(("../timing/timinghists"+to_string(lo)+"_"+to_string(hi)+".root").c_str(),"RECREATE");
  outf->cd();
  for(int i=0; i<ncut; ++i)
    {
      h3_pt_dt_t[i]->Write();
      h2_pt_t[i]->Write();
    }

  outf->Write();

  return 0;
      
}
      
