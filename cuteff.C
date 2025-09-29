float get_dphi(float phi1, float phi2)
{
  float dphi = abs(phi1-phi2);
  if(abs(dphi - 2*M_PI) < dphi) dphi = abs(dphi - 2*M_PI);
  return dphi;
}

float get_emcal_mineta_zcorrected(float zvertex) {
  float minz_EM = -130.23;
  float radius_EM = 93.5;
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float maxz_EM = 130.23;
  float radius_EM = 93.5;
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float minz_IH = -170.299;
  float radius_IH = 127.503;
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float maxz_IH = 170.299;
  float radius_IH = 127.503;
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float minz_OH = -301.683;
  float radius_OH = 225.87;
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float maxz_OH = 301.683;
  float radius_OH = 225.87;
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zvertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zvertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zvertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zvertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zvertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zvertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zvertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

std::vector<vector<float>> match_truth_reco(std::vector<vector<float>> truthjets, std::vector<vector<float>> recojets)
{
  std::vector<vector<float>> matched_pts = {};

  std::sort(truthjets.begin(),truthjets.end(),[] (const std::vector<float> &a, const std::vector<float> &b){return a[0] > b[0];});
  std::sort(recojets.begin(),recojets.end(),[] (const std::vector<float> &a, const std::vector<float> &b){return a[0] > b[0];});

  float lrohf = recojets.at(0).at(4);
  float lremf = recojets.at(0).at(3);

  for(int i=0; i<truthjets.size(); ++i)
    {
      std::vector<float> matched_pt = {};
      for(int j=0; j<recojets.size(); ++j)
	{
	  if(recojets.at(j).at(5) != 0) continue;
	  float deta = abs(truthjets.at(i).at(1) - recojets.at(j).at(1));
	  float dphi = get_dphi(truthjets.at(i).at(2),recojets.at(j).at(2));
	  float dR = sqrt(deta*deta+dphi*dphi);
	  if(dR < 0.3)
	    {
	      matched_pt.push_back(truthjets.at(i).at(0));
	      matched_pt.push_back(recojets.at(j).at(0));
	      matched_pt.push_back(lremf);
	      matched_pt.push_back(lrohf);
	      matched_pt.push_back(recojets.at(j).at(1));
	      recojets.at(j).at(5) = 1;
	      matched_pts.push_back(matched_pt);
	      break;
	    }
	}
    }
  
  return matched_pts;
}

int cuteff(int lo, int hi, int type)
{
  int jet_n, tjet_n, failscut;
  float frcem[100];
  float frcoh[100];
  float jet_pt[100];
  float tjet_pt[100];
  float jet_eta[100];
  float tjet_eta[100];
  float jet_phi[100];
  float tjet_phi[100];
  float zvtx;
  
  TChain* jet_tree = new TChain("jet_tree");
  
  float scalefactor = 7.2695e-9;
  if(type==2) scalefactor = 1.034e-11;
  if(type==3) scalefactor = 2.502e-6;

  string tempinfilename;
  string listname = "chi2files";
  string samplename;
  if(type == 1) samplename = "50";
  if(type == 2) samplename = "70";
  if(type == 3) samplename = "30";
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
      jet_tree->Add(tempinfilename.c_str());
      ++counter;
      if(counter >= hi)
        {
          cout << counter << endl;
          break;
        }
    }

  jet_tree->SetBranchStatus("*",0);
  jet_tree->SetBranchStatus("jet_n",1);
  jet_tree->SetBranchStatus("failscut",1);
  jet_tree->SetBranchStatus("tjet_n",1);
  jet_tree->SetBranchStatus("alljetfrcem",1);
  jet_tree->SetBranchStatus("alljetfrcoh",1);
  jet_tree->SetBranchStatus("jet_pt",1);
  jet_tree->SetBranchStatus("jet_eta",1);
  jet_tree->SetBranchStatus("jet_phi",1);
  jet_tree->SetBranchStatus("tjet_pt",1);
  jet_tree->SetBranchStatus("tjet_eta",1);
  jet_tree->SetBranchStatus("tjet_phi",1);
  jet_tree->SetBranchStatus("zvtx",1);

  jet_tree->SetBranchAddress("failscut",&failscut);
  jet_tree->SetBranchAddress("jet_n",&jet_n);
  jet_tree->SetBranchAddress("tjet_n",&tjet_n);
  jet_tree->SetBranchAddress("alljetfrcem",frcem);
  jet_tree->SetBranchAddress("alljetfrcoh",frcoh);
  jet_tree->SetBranchAddress("jet_pt",jet_pt);
  jet_tree->SetBranchAddress("jet_eta",jet_eta);
  jet_tree->SetBranchAddress("jet_phi",jet_phi);
  jet_tree->SetBranchAddress("tjet_pt",tjet_pt);
  jet_tree->SetBranchAddress("tjet_eta",tjet_eta);
  jet_tree->SetBranchAddress("tjet_phi",tjet_phi);
  jet_tree->SetBranchAddress("zvtx",&zvtx);

  TH3D* h3_pt_lem_loh = new TH3D("h3_pt_lem_loh",";p_{T}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_tpt_lem_loh = new TH3D("h3_tpt_lem_loh",";p_{T}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);

  TH3D* h3_pt_tpt_loh = new TH3D("h3_pt_tpt_loh",";p_{T}^{reco} [GeV];p_{T}^{truth};E_{reco}^{OH}/E_{reco}",100,0,100,100,0,100,120,-0.1,1.1);
  TH3D* h3_pt_tpt_lem = new TH3D("h3_pt_tpt_lem",";p_{T}^{truth} [GeV];p_{T}^{truth};E_{reco}^{EM}/E_{reco} Matched",100,0,100,100,0,100,120,-0.1,1.1);

  TH3D* h3_pt_tpt_loh_etacut = new TH3D("h3_pt_tpt_loh_etacut",";p_{T}^{reco} [GeV];p_{T}^{truth};E_{reco}^{OH}/E_{reco}",100,0,100,100,0,100,120,-0.1,1.1);
  TH3D* h3_pt_tpt_lem_etacut = new TH3D("h3_pt_tpt_lem_etacut",";p_{T}^{truth} [GeV];p_{T}^{truth};E_{reco}^{EM}/E_{reco} Matched",100,0,100,100,0,100,120,-0.1,1.1);

  TH3D* h3_pt_lem_loh_etacut = new TH3D("h3_pt_lem_loh_etacut",";p_{T}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_tpt_lem_loh_etacut = new TH3D("h3_tpt_lem_loh_etacut",";p_{T}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);

  TH3D* h3_pt_lem_loh_z100 = new TH3D("h3_pt_lem_loh_z100",";p_{T}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_tpt_lem_loh_z100 = new TH3D("h3_tpt_lem_loh_z100",";p_{T}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);

  TH3D* h3s_pt_lem_loh = new TH3D("h3s_pt_lem_loh",";p_{T}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3s_tpt_lem_loh = new TH3D("h3s_tpt_lem_loh",";p_{T}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);

  
  TH3D* h3_lpt_lem_loh = new TH3D("h3_lpt_lem_loh",";p_{T,lead}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_ltpt_lem_loh = new TH3D("h3_ltpt_lem_loh",";p_{T,lead}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);

  TH3D* h3_pt_lem_loh_dijet = new TH3D("h3_pt_lem_loh_dijet",";p_{T}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_tpt_lem_loh_dijet = new TH3D("h3_tpt_lem_loh_dijet",";p_{T}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_lpt_lem_loh_dijet = new TH3D("h3_lpt_lem_loh_dijet",";p_{T,lead}^{reco} [GeV];E_{reco}^{EM}/E_{reco};E_{reco}^{OH}/E_{reco}",100,0,100,120,-0.1,1.1,120,-0.1,1.1);
  TH3D* h3_ltpt_lem_loh_dijet = new TH3D("h3_ltpt_lem_loh_dijet",";p_{T,lead}^{truth} [GeV];E_{reco}^{EM}/E_{reco} Matched;E_{reco}^{OH}/E_{reco} Matched",100,0,100,120,-0.1,1.1,120,-0.1,1.1);

  
  
  for(int i=0; i<jet_tree->GetEntries(); ++i)
    {
      if((i+1)%100 == 0) cout << i << endl;
      jet_tree->GetEntry(i);
      std::vector<vector<float>> recojets = {};
      std::vector<vector<float>> truthjets = {};

      for(int j=0; j<jet_n; ++j)
	{
	  std::vector<float> jet = {jet_pt[j], jet_eta[j], jet_phi[j], frcem[j], frcoh[j], 0};
	  recojets.push_back(jet);
	}
      float lpt = 0;
      for(int j=0; j<tjet_n; ++j)
	{
	  std::vector<float> jet = {tjet_pt[j], tjet_eta[j], tjet_phi[j]};
	  if(tjet_pt[j] > lpt) lpt = tjet_pt[j];
	  truthjets.push_back(jet);
	}
      if(type==1 && (lpt < 52 || lpt > 71)) continue;
      if(type==2 && lpt < 71) continue;
      if(type==3 && lpt > 52) continue;
      
      std::vector<vector<float>> matched_jets = match_truth_reco(truthjets,recojets);

      if(matched_jets.size() > 0)
	{
	  h3_lpt_lem_loh->Fill(matched_jets.at(0).at(1),matched_jets.at(0).at(2), matched_jets.at(0).at(3),scalefactor);
	  h3_ltpt_lem_loh->Fill(matched_jets.at(0).at(0),matched_jets.at(0).at(2), matched_jets.at(0).at(3),scalefactor);
	  
	  if(failscut%2==0)
	    {
	      h3_lpt_lem_loh_dijet->Fill(matched_jets.at(0).at(1),matched_jets.at(0).at(2), matched_jets.at(0).at(3),scalefactor);
	      h3_ltpt_lem_loh_dijet->Fill(matched_jets.at(0).at(0),matched_jets.at(0).at(2), matched_jets.at(0).at(3),scalefactor);
	    }
	  
	  for(int j=0; j<matched_jets.size(); ++j)
	    {
	      h3_pt_lem_loh->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
	      h3_tpt_lem_loh->Fill(matched_jets.at(j).at(0),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);


	      h3_pt_tpt_loh->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(0), matched_jets.at(j).at(3),scalefactor);
	      h3_pt_tpt_lem->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(0), matched_jets.at(j).at(2),scalefactor);


	      if(abs(zvtx) < 100)
		{
		  h3_pt_lem_loh_z100->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
		  h3_tpt_lem_loh_z100->Fill(matched_jets.at(j).at(0),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
		}
	      
	      if(check_bad_jet_eta(matched_jets.at(j).at(4),zvtx,0.4)) continue;

	      h3_pt_tpt_loh_etacut->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(0), matched_jets.at(j).at(3),scalefactor);
	      h3_pt_tpt_lem_etacut->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(0), matched_jets.at(j).at(2),scalefactor);
	      
	      h3_pt_lem_loh_etacut->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
	      h3_tpt_lem_loh_etacut->Fill(matched_jets.at(j).at(0),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
	      if(failscut%2==0)
		{
		  h3_pt_lem_loh_dijet->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
		  h3_tpt_lem_loh_dijet->Fill(matched_jets.at(j).at(0),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
		}

	      if(type==1 && matched_jets.at(j).at(1) > 74) continue;
	      if(type==3 && matched_jets.at(j).at(1) > 55) continue;
      
	      h3s_pt_lem_loh->Fill(matched_jets.at(j).at(1),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
	      h3s_tpt_lem_loh->Fill(matched_jets.at(j).at(0),matched_jets.at(j).at(2), matched_jets.at(j).at(3),scalefactor);
	      
	    }
	}
    }

  TFile* outf = TFile::Open(("../cuteff/cuteff_"+to_string(lo)+"_"+to_string(hi)+"_"+samplename+".root").c_str(),"RECREATE");

  outf->cd();

  h3_pt_lem_loh->Write();
  h3_tpt_lem_loh->Write();

  h3_pt_tpt_lem->Write();
  h3_pt_tpt_loh->Write();
  h3_pt_tpt_lem_etacut->Write();
  h3_pt_tpt_loh_etacut->Write();

  h3_pt_lem_loh_z100->Write();
  h3_tpt_lem_loh_z100->Write();

  h3_pt_lem_loh_etacut->Write();
  h3_tpt_lem_loh_etacut->Write();

  h3s_pt_lem_loh->Write();
  h3s_tpt_lem_loh->Write();
  
  h3_lpt_lem_loh->Write();
  h3_ltpt_lem_loh->Write();

  h3_pt_lem_loh_dijet->Write();
  h3_tpt_lem_loh_dijet->Write();
  h3_lpt_lem_loh_dijet->Write();
  h3_ltpt_lem_loh_dijet->Write();
  
  return 0;
}
