#include <iomanip>
#include "Chi2checker.h"
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <g4main/PHG4VtxPoint.h>
#include <mbd/MbdPmtHit.h>
#include <jetbackground/TowerBackgroundv1.h>
#include <cmath>
#include <mbd/MbdOut.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerv1.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos
#include <phool/recoConsts.h>
#include <iostream>
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <TH2D.h>
#include <TColor.h>
#include <TMarker.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TH1D.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffarawobjects/Gl1Packetv3.h>
#include <TLorentzVector.h>
#include <ffaobjects/EventHeader.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <calobase/TowerInfoDefs.h>
#include <ffarawobjects/CaloPacketv1.h>
#include <phparameter/PHParameters.h>
using namespace std;
static const float radius_EM = 93.5;
static const float minz_EM = -130.23;
static const float maxz_EM = 130.23;

static const float radius_IH = 127.503;
static const float minz_IH = -170.299;
static const float maxz_IH = 170.299;

static const float radius_OH = 225.87;
static const float minz_OH = -301.683;
static const float maxz_OH = 301.683;
float get_emcal_mineta_zcorrected(float zvertex) {
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
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

float Chi2checker::bintoeta_hc(int etabin)
{
  return (2.2*etabin)/24 - 1.1;
}
float Chi2checker::bintophi_hc(int phibin)
{
  return (2*M_PI*phibin)/64;
}

float bintoeta_em(int etabin)
{
  return (2.2*etabin)/96-1.1;
}

float bintophi_em(int phibin)
{
  return (2*M_PI*phibin)/256;
}





//____________________________________________________________________________..
Chi2checker::Chi2checker(const std::string &filename, const std::string &name, const int debug, const std::string &wfilename, const int dowf, const bool isdat, const int doall60, const int dotruthpar):
  SubsysReco(name), _cutParams(name)
{
  _dotruthpar = dotruthpar;
  _isdat = isdat;
  _name = name;
  _debug = debug;
  _filename = filename;
  _wfilename = wfilename;
  _nprocessed = 0;
  _dowf = dowf;
  jet_tree = new TTree("jet_tree","a persevering date tree");
  if(_dowf) _wft = new TTree("wft","a stupid waveform tree");
  _doall60 = doall60;
  if(_doall60 && _dowf) _wff = new TFile(_wfilename.c_str(),"RECREATE");
  if(_doall60 && _dowf) _wff->cd();
  if(_doall60 && _dowf) _wft->SetDirectory(_wff);
  //if(_doall60 || !_isdat)
  _f = new TFile(_filename.c_str(), "RECREATE");
  _f->cd();
  jet_tree->SetDirectory(_f);
}

//____________________________________________________________________________..
Chi2checker::~Chi2checker()
{

}

//____________________________________________________________________________..
int Chi2checker::Init(PHCompositeNode *topNode)
{
  
  if(_debug > 1) cout << "Begin init: " << endl;

  jet_tree->Branch("zvtx",&_zvtx,"zvtx/F");
  if(!_isdat) jet_tree->Branch("tzvtx",&_tzvtx,"tzvtx/F");
  jet_tree->Branch("triggervec",&_triggervec,"triggervec/l");
  jet_tree->Branch("runnum",&_runnum,"runnum/I");
  jet_tree->Branch("evtnum",&_evtnum,"evtnum/I");
  jet_tree->Branch("failscut",&_failscut,"failscut/I");
  jet_tree->Branch("mbdavgt",_mbdavgt,"mbdavgt[2]/F");
  jet_tree->Branch("mbdhit",_mbdhit,"mbdhit[2]/i");
  //jet_tree->Branch("bbfqavec",&_bbfqavec,"bbfqavec/i");

  jet_tree->Branch("r02_jet_n",&_r02_jet_n,"r02_jet_n/I");
  jet_tree->Branch("calib_r02_jet_n",&_calib_r02_jet_n,"calib_r02_jet_n/I");
  jet_tree->Branch("allr02_jetfrcoh",_allr02_jetfrcoh,"allr02_jetfrcoh[r02_jet_n]/F");
  jet_tree->Branch("allr02_jetfrcem",_allr02_jetfrcem,"allr02_jetfrcem[r02_jet_n]/F");
  jet_tree->Branch("r02_jet_et",_r02_jet_et,"r02_jet_et[r02_jet_n]/F");
  jet_tree->Branch("r02_jet_etrans",_r02_jet_etrans,"r02_jet_etrans[r02_jet_n]/F");
  jet_tree->Branch("r02_jet_pt",_r02_jet_pt,"r02_jet_pt[r02_jet_n]/F");
  jet_tree->Branch("r02_jet_t",_r02_jet_t,"r02_jet_t[r02_jet_n]/F");
  jet_tree->Branch("r02_jet_pt_calib",_r02_jet_pt_calib,"r02_jet_pt_calib[calib_r02_jet_n]/F");
  jet_tree->Branch("r02_jet_eta",_r02_jet_eta,"r02_jet_eta[r02_jet_n]/F");
  jet_tree->Branch("r02_jet_phi",_r02_jet_phi,"r02_jet_phi[r02_jet_n]/F");

  jet_tree->Branch("r03_jet_n",&_r03_jet_n,"r03_jet_n/I");
  jet_tree->Branch("calib_r03_jet_n",&_calib_r03_jet_n,"calib_r03_jet_n/I");
  jet_tree->Branch("allr03_jetfrcoh",_allr03_jetfrcoh,"allr03_jetfrcoh[r03_jet_n]/F");
  jet_tree->Branch("allr03_jetfrcem",_allr03_jetfrcem,"allr03_jetfrcem[r03_jet_n]/F");
  jet_tree->Branch("r03_jet_et",_r03_jet_et,"r03_jet_et[r03_jet_n]/F");
  jet_tree->Branch("r03_jet_etrans",_r03_jet_etrans,"r03_jet_etrans[r03_jet_n]/F");
  jet_tree->Branch("r03_jet_pt",_r03_jet_pt,"r03_jet_pt[r03_jet_n]/F");
  jet_tree->Branch("r03_jet_t",_r03_jet_t,"r03_jet_t[r03_jet_n]/F");
  jet_tree->Branch("r03_jet_pt_calib",_r03_jet_pt_calib,"r03_jet_pt_calib[calib_r03_jet_n]/F");
  jet_tree->Branch("r03_jet_eta",_r03_jet_eta,"r03_jet_eta[r03_jet_n]/F");
  jet_tree->Branch("r03_jet_phi",_r03_jet_phi,"r03_jet_phi[r03_jet_n]/F");


  jet_tree->Branch("jet_n",&_jet_n,"jet_n/I");
  jet_tree->Branch("calib_jet_n",&_calib_jet_n,"calib_jet_n/I");
  jet_tree->Branch("alljetfrcoh",_alljetfrcoh,"alljetfrcoh[jet_n]/F");
  jet_tree->Branch("alljetfrcem",_alljetfrcem,"alljetfrcem[jet_n]/F");
  jet_tree->Branch("jet_et",_jet_et,"jet_et[jet_n]/F");
  jet_tree->Branch("jet_etrans",_jet_etrans,"jet_etrans[jet_n]/F");
  jet_tree->Branch("jet_pt",_jet_pt,"jet_pt[jet_n]/F");
  jet_tree->Branch("jet_t",_jet_t,"jet_t[jet_n]/F");
  jet_tree->Branch("jet_pt_calib",_jet_pt_calib,"jet_pt_calib[calib_jet_n]/F");
  jet_tree->Branch("jet_eta",_jet_eta,"jet_eta[jet_n]/F");
  jet_tree->Branch("jet_phi",_jet_phi,"jet_phi[jet_n]/F");


  jet_tree->Branch("r06_jet_n",&_r06_jet_n,"r06_jet_n/I");
  jet_tree->Branch("calib_r06_jet_n",&_calib_r06_jet_n,"calib_r06_jet_n/I");
  jet_tree->Branch("allr06_jetfrcoh",_allr06_jetfrcoh,"allr06_jetfrcoh[r06_jet_n]/F");
  jet_tree->Branch("allr06_jetfrcem",_allr06_jetfrcem,"allr06_jetfrcem[r06_jet_n]/F");
  jet_tree->Branch("r06_jet_et",_r06_jet_et,"r06_jet_et[r06_jet_n]/F");
  jet_tree->Branch("r06_jet_etrans",_r06_jet_etrans,"r06_jet_etrans[r06_jet_n]/F");
  jet_tree->Branch("r06_jet_pt",_r06_jet_pt,"r06_jet_pt[r06_jet_n]/F");
  jet_tree->Branch("r06_jet_t",_r06_jet_t,"r06_jet_t[r06_jet_n]/F");
  jet_tree->Branch("r06_jet_pt_calib",_r06_jet_pt_calib,"r06_jet_pt_calib[calib_r06_jet_n]/F");
  jet_tree->Branch("r06_jet_eta",_r06_jet_eta,"r06_jet_eta[r06_jet_n]/F");
  jet_tree->Branch("r06_jet_phi",_r06_jet_phi,"r06_jet_phi[r06_jet_n]/F");


  jet_tree->Branch("r08_jet_n",&_r08_jet_n,"r08_jet_n/I");
  jet_tree->Branch("calib_r08_jet_n",&_calib_r08_jet_n,"calib_r08_jet_n/I");
  jet_tree->Branch("allr08_jetfrcoh",_allr08_jetfrcoh,"allr08_jetfrcoh[r08_jet_n]/F");
  jet_tree->Branch("allr08_jetfrcem",_allr08_jetfrcem,"allr08_jetfrcem[r08_jet_n]/F");
  jet_tree->Branch("r08_jet_et",_r08_jet_et,"r08_jet_et[r08_jet_n]/F");
  jet_tree->Branch("r08_jet_etrans",_r08_jet_etrans,"r08_jet_etrans[r08_jet_n]/F");
  jet_tree->Branch("r08_jet_pt",_r08_jet_pt,"r08_jet_pt[r08_jet_n]/F");
  jet_tree->Branch("r08_jet_t",_r08_jet_t,"r08_jet_t[r08_jet_n]/F");
  jet_tree->Branch("r08_jet_pt_calib",_r08_jet_pt_calib,"r08_jet_pt_calib[calib_r08_jet_n]/F");
  jet_tree->Branch("r08_jet_eta",_r08_jet_eta,"r08_jet_eta[r08_jet_n]/F");
  jet_tree->Branch("r08_jet_phi",_r08_jet_phi,"r08_jet_phi[r08_jet_n]/F");

  if(!_isdat)
    {
      jet_tree->Branch("r02_tjet_n",&_r02_tjet_n,"r02_tjet_n/I");
      jet_tree->Branch("r02_tjet_e",_r02_tjet_e,"r02_tjet_e[r02_tjet_n]/F");
      jet_tree->Branch("r02_tjet_pt",_r02_tjet_pt,"r02_tjet_pt[r02_tjet_n]/F");
      jet_tree->Branch("r02_tjet_eta",_r02_tjet_eta,"r02_tjet_eta[r02_tjet_n]/F");
      jet_tree->Branch("r02_tjet_phi",_r02_tjet_phi,"r02_tjet_phi[r02_tjet_n]/F");

      jet_tree->Branch("r03_tjet_n",&_r03_tjet_n,"r03_tjet_n/I");
      jet_tree->Branch("r03_tjet_e",_r03_tjet_e,"r03_tjet_e[r03_tjet_n]/F");
      jet_tree->Branch("r03_tjet_pt",_r03_tjet_pt,"r03_tjet_pt[r03_tjet_n]/F");
      jet_tree->Branch("r03_tjet_eta",_r03_tjet_eta,"r03_tjet_eta[r03_tjet_n]/F");
      jet_tree->Branch("r03_tjet_phi",_r03_tjet_phi,"r03_tjet_phi[r03_tjet_n]/F");

      jet_tree->Branch("tjet_n",&_tjet_n,"tjet_n/I");
      jet_tree->Branch("tjet_e",_tjet_e,"tjet_e[tjet_n]/F");
      jet_tree->Branch("tjet_pt",_tjet_pt,"tjet_pt[tjet_n]/F");
      jet_tree->Branch("tjet_eta",_tjet_eta,"tjet_eta[tjet_n]/F");
      jet_tree->Branch("tjet_phi",_tjet_phi,"tjet_phi[tjet_n]/F");

      jet_tree->Branch("r06_tjet_n",&_r06_tjet_n,"r06_tjet_n/I");
      jet_tree->Branch("r06_tjet_e",_r06_tjet_e,"r06_tjet_e[r06_tjet_n]/F");
      jet_tree->Branch("r06_tjet_pt",_r06_tjet_pt,"r06_tjet_pt[r06_tjet_n]/F");
      jet_tree->Branch("r06_tjet_eta",_r06_tjet_eta,"r06_tjet_eta[r06_tjet_n]/F");
      jet_tree->Branch("r06_tjet_phi",_r06_tjet_phi,"r06_tjet_phi[r06_tjet_n]/F");

      jet_tree->Branch("r08_tjet_n",&_r08_tjet_n,"r08_tjet_n/I");
      jet_tree->Branch("r08_tjet_e",_r08_tjet_e,"r08_tjet_e[r08_tjet_n]/F");
      jet_tree->Branch("r08_tjet_pt",_r08_tjet_pt,"r08_tjet_pt[r08_tjet_n]/F");
      jet_tree->Branch("r08_tjet_eta",_r08_tjet_eta,"r08_tjet_eta[r08_tjet_n]/F");
      jet_tree->Branch("r08_tjet_phi",_r08_tjet_phi,"r08_tjet_phi[r08_tjet_n]/F");
      if(_dotruthpar) jet_tree->Branch("truthparenergy",&_truthparenergy);
      if(_dotruthpar) jet_tree->Branch("truthpareta",&_truthpareta);
      if(_dotruthpar) jet_tree->Branch("truthparphi",&_truthparphi);
      if(_dotruthpar) jet_tree->Branch("truthparpt",&_truthparpt);
      if(_dotruthpar) jet_tree->Branch("truthparid",&_truthparid,"truthparid/I");
    }

  


  
  _mbevt = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Chi2checker::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

unsigned int get_towerindex(int packet, int pchan)
{
  return (192*(packet%1000-1))+pchan;
}

int get_bindex(int packet, int pchan, string what)
{
  unsigned int key = -1;
  if(packet/1000 == 6)
    {
      key = TowerInfoDefs::encode_emcal(get_towerindex(packet,pchan));
    }
  else if(packet/1000 == 7 || packet / 1000 == 8)
    {
      key = TowerInfoDefs::encode_hcal(get_towerindex(packet,pchan));
    }
  else
    {
      return -1;
    }

  if(what == "eta") return TowerInfoDefs::getCaloTowerEtaBin(key);
  else if(what == "phi") return TowerInfoDefs::getCaloTowerPhiBin(key);
  else return -1;
}

float Chi2checker::getEtaFromBin(int binEta)
{
  return ((2.2*binEta)/nx)-1.1;
}

float Chi2checker::getPhiFromBin(int binPhi)
{
  return ((2*M_PI*binPhi)/ny);}


float getEtaFromBinEM(int binEta)
{
  return ((2.2*binEta)/64)-1.1;
}

float getPhiFromBinEM(int binPhi)
{
  return ((2*M_PI*binPhi)/256);
}
//____________________________________________________________________________..

void print_debug(float jet_eta, float jet_phi, float tower_eta, float tower_phi, float dphi, float deta)
{
  cout << "printing debug info for dphi deta:" << endl;
  cout << "jet eta/phi: " << jet_eta << " " << jet_phi << endl;
  cout << "tower eta/phi: " << tower_eta << " " << tower_phi << endl;
  cout << "deta dphi:" << deta << " " << dphi << endl;
}



int Chi2checker::fill_jet_quantities(PHCompositeNode* topNode, std::string jet_nodename, int& njet, int& calib_njet, float* jet_pt, float* jet_pt_calib, float* jet_eta, float* jet_phi, float* jet_e, float* jet_et, float* jet_t, float* jet_frcem, float* jet_frcoh)
{
  JetContainer* jets = findNode::getClass<JetContainerv1>(topNode,jet_nodename);
  JetContainer* cjets = findNode::getClass<JetContainerv1>(topNode,(jet_nodename+"_calib"));

  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  
  std::vector<int> record_index = {};
  int index_counter = 0;
  calib_njet = 0;
  njet = 0;
  if(cjets)
    {
      for(unsigned int i=0; i<cjets->size(); ++i)
	{
	  Jet* jet = cjets->get_jet(i);
	  if(jet->get_pt() < 1)
	    {
	      ++index_counter;
	      continue;
	    }
	  record_index.push_back(index_counter);
	  ++index_counter;
	  jet_pt_calib[calib_njet] = jet->get_pt();
	  ++calib_njet;
	  if(calib_njet > 98) break;
	}
    }
  else
    {
      cout << "no calib jets" << endl;
      return 1;
    }

  index_counter = 0;

  if(jets)
    {
      for(unsigned int i=0; i<jets->size(); ++i)
	{
	  if (std::find(record_index.begin(), record_index.end(), index_counter) == record_index.end())
            {
              ++index_counter;
              continue;
            }
	  ++index_counter;
	  Jet* jet = jets->get_jet(i);

	  if(jet)
	    {
	      jet_frcem[njet] = 0;
	      jet_frcoh[njet] = 0;
	      jet_pt[njet] = jet->get_pt();
	      jet_eta[njet] = jet->get_eta();
	      jet_phi[njet] = jet->get_phi();
	      jet_e[njet] = jet->get_e();
	      jet_et[njet] = jet_e[njet]/cosh(jet_eta[njet]);
	      float jet_t_esum = 0;

	      for(auto comp: jet->get_comp_vec())
		{
		  unsigned int channel = comp.second;
		  TowerInfo* tower;
		  if(comp.first==5 || comp.first == 26)
		    {
		      tower = towers[1]->get_tower_at_channel(channel);
		      float towerE = tower->get_energy();
		      if(towerE > 0.1)
			{
			  jet_t_esum += towerE;
			  jet_t[njet] += towerE*tower->get_time();
			}
		    }
		  else if(comp.first==7 || comp.first==27)
		    {
		      tower = towers[2]->get_tower_at_channel(channel);
		      float towerE = tower->get_energy();
		      if(towerE > 0.1)
			{
			  jet_t_esum += towerE;
			  jet_t[njet] += towerE*tower->get_time();
			}
		      jet_frcoh[njet] += towerE;
		    }
		  else if(comp.first==13 || comp.first == 28 || comp.first == 25)
		    {
		      tower = towers[0]->get_tower_at_channel(channel);
		      float towerE = tower->get_energy();
		      if(towerE > 0.1)
			{
			  jet_t_esum += towerE;
			  jet_t[njet] += towerE*tower->get_time();
			}
		      jet_frcem[njet] += towerE;
		    }
		}
	      jet_frcem[njet] /= jet_e[njet];
	      jet_frcoh[njet] /= jet_e[njet];
	      jet_t[njet] /= jet_t_esum;
	      ++njet;
	      if(njet > 98) break;
	    }
	}
    }
  else
    {
      return 2;
    }

  return 0;
  
}


int Chi2checker::fill_tjet_quantities(PHCompositeNode* topNode, std::string jet_nodename, int& njet, float* jet_pt, float* jet_eta, float* jet_phi, float* jet_e)
{
  JetContainer* truthjets = findNode::getClass<JetContainerv1>(topNode,jet_nodename);

  njet = 0;
  if(truthjets)
    {
      for(unsigned int i=0; i<truthjets->size(); ++i)
	{
	  Jet* jet = truthjets->get_jet(i);
	  jet_pt[njet] = jet->get_pt();
	  if(jet_pt[njet] < 1) continue;
	  jet_eta[njet] = jet->get_eta();
	  jet_phi[njet] = jet->get_phi();
	  jet_e[njet] = jet->get_e();
	  ++njet;
	  if(njet > 98) break;
	}
    }
  else
    {
      return 1;
    }

  return 0;
}


int Chi2checker::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) cout << endl << endl << endl << "Chi2checker: Beginning event processing" << endl;
  if(_nprocessed % 1000 == 0) cout << "processing event " << _nprocessed << endl;


  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
  if (!_isdat)
    {
      PHG4VtxPoint *gvertex = truthinfo->GetPrimaryVtx(truthinfo->GetPrimaryVertexIndex());
      _tzvtx = gvertex->get_z();
    }

  int isjettrig = 0;
  int isjmbtrig = 0;
  
  if(_isdat)
    {
      
      Gl1Packetv3* gl1 = findNode::getClass<Gl1Packetv3>(topNode, "14001");
      if(!gl1)
	{
	  cout << "No trigger info!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
      if(_debug > 1) cout << "Chi2checker: Getting gl1 trigger vector from: " << gl1 << endl;
      _triggervec = gl1->getScaledVector();
      
      isjettrig = (_triggervec >> 22) & 1;
      isjmbtrig = (_triggervec >> 18) & 1;
      if(_nprocessed == 0)
	{
	  _prevraw18 = gl1->lValue(18,0);
	  _prevraw22 = gl1->lValue(22,0);
	  _prevlive18 = gl1->lValue(18,1);
	  _prevlive22 = gl1->lValue(22,1);
	  _isbadlive = 0;
	}
      long long unsigned int currentlive18, currentlive22, currentraw18, currentraw22;
      currentlive18 = gl1->lValue(18,1);
      currentlive22 = gl1->lValue(22,1);
      currentraw18 = gl1->lValue(18,0);
      currentraw22 = gl1->lValue(22,0);
      long long unsigned int live18diff = currentlive18 - _prevlive18;
      long long unsigned int live22diff = currentlive22 - _prevlive22;
      long long unsigned int raw18diff = currentraw18 - _prevraw18;
      long long unsigned int raw22diff = currentraw22 - _prevraw22;
      if(_nprocessed != 0)
	{
	  if(isjettrig && (((float)live22diff)/raw22diff < 0.1 || ((float)live18diff)/raw18diff < 0.1)) _isbadlive = 1;
	  else _isbadlive = 0;
	}
      else
	{
	  _isbadlive = 0;
	}
      
      _prevraw18 = currentraw18;
      _prevraw22 = currentraw22;
      _prevlive18 = currentlive18;
      _prevlive22 = currentlive22;
    }
  ++_nprocessed;
  if(_isdat)
    {
      int runnumber = 0;
      int evtnum = 0;
      EventHeader *runheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
      if (!runheader)
	{
	  std::cout << "can't find runheader" << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      runnumber = runheader->get_RunNumber();
      evtnum = runheader->get_EvtSequence();
      _runnum = runnumber;
      _evtnum = evtnum;

    }
  /*
  PHNodeIterator itNode(topNode);
  PHCompositeNode* parNode = dynamic_cast<PHCompositeNode*>(itNode.findFirst("PHCompositeNode","PAR"));
  PdbParameterMap* flagNode;
  */

  float zvtx = NAN;
  GlobalVertexMap* vertexmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  if(vertexmap)
    {
      std::vector<GlobalVertex*> vertices = vertexmap->get_gvtxs_with_type({GlobalVertex::MBD});
      if(!vertices.empty())
	{
	  if(vertices.at(0))
	    {
	      zvtx = vertices.at(0)->get_z();
	    }
	}
    }
  else
    {
      cout << "no global vertex map! Abort run!" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  _zvtx = zvtx;
  if(std::isnan(_zvtx) || std::isnan(zvtx))
    {
      _zvtx = 0;
      zvtx = 0;
    }

  int truth_rpass[5] = {0};
  int reco_rpass[5] = {0};

  if(!_isdat)
    {
      truth_rpass[0] = fill_tjet_quantities(topNode,"AntiKt_Truth_r02",_r02_tjet_n,_r02_tjet_pt,_r02_tjet_eta,_r02_tjet_phi,_r02_tjet_e);
      truth_rpass[1] = fill_tjet_quantities(topNode,"AntiKt_Truth_r03",_r03_tjet_n,_r03_tjet_pt,_r03_tjet_eta,_r03_tjet_phi,_r03_tjet_e);
      truth_rpass[2] = fill_tjet_quantities(topNode,"AntiKt_Truth_r04",_tjet_n,_tjet_pt,_tjet_eta,_tjet_phi,_tjet_e);
      truth_rpass[3] = fill_tjet_quantities(topNode,"AntiKt_Truth_r06",_r06_tjet_n,_r06_tjet_pt,_r06_tjet_eta,_r06_tjet_phi,_r06_tjet_e);
      truth_rpass[4] = fill_tjet_quantities(topNode,"AntiKt_Truth_r08",_r08_tjet_n,_r08_tjet_pt,_r08_tjet_eta,_r08_tjet_phi,_r08_tjet_e);

      for(int i=0; i<5; ++i)
	{
	  if(truth_rpass[i])
	    {
	      cout << "some truth jet node not found!" << endl;
	      return Fun4AllReturnCodes::ABORTRUN;
	    }
	}
    }

  reco_rpass[0] = fill_jet_quantities(topNode,"AntiKt_unsubtracted_r02",_r02_jet_n,_calib_r02_jet_n,_r02_jet_pt,_r02_jet_pt_calib,_r02_jet_eta,_r02_jet_phi,_r02_jet_et,_r02_jet_etrans,_r02_jet_t,_allr02_jetfrcem,_allr02_jetfrcoh);
  reco_rpass[1] = fill_jet_quantities(topNode,"AntiKt_unsubtracted_r03",_r03_jet_n,_calib_r03_jet_n,_r03_jet_pt,_r03_jet_pt_calib,_r03_jet_eta,_r03_jet_phi,_r03_jet_et,_r03_jet_etrans,_r03_jet_t,_allr03_jetfrcem,_allr03_jetfrcoh);
  reco_rpass[2] = fill_jet_quantities(topNode,"AntiKt_unsubtracted_r04",_jet_n,_calib_jet_n,_jet_pt,_jet_pt_calib,_jet_eta,_jet_phi,_jet_et,_jet_etrans,_jet_t,_alljetfrcem,_alljetfrcoh);
  reco_rpass[3] = fill_jet_quantities(topNode,"AntiKt_unsubtracted_r06",_r06_jet_n,_calib_r06_jet_n,_r06_jet_pt,_r06_jet_pt_calib,_r06_jet_eta,_r06_jet_phi,_r06_jet_et,_r06_jet_etrans,_r06_jet_t,_allr06_jetfrcem,_allr06_jetfrcoh);
  reco_rpass[4] = fill_jet_quantities(topNode,"AntiKt_unsubtracted_r08",_r08_jet_n,_calib_r08_jet_n,_r08_jet_pt,_r08_jet_pt_calib,_r08_jet_eta,_r08_jet_phi,_r08_jet_et,_r08_jet_etrans,_r08_jet_t,_allr08_jetfrcem,_allr08_jetfrcoh);

  for(int i=0; i<5; ++i)
    {
      if(reco_rpass[i] == 1)
	{
	  cout << "some reco calib jet node not found!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
      if(reco_rpass[i] == 2)
	{
	  cout << "some reco uncalib jet node not found!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
      
    }

  

  MbdPmtContainer* mbdpmt = findNode::getClass<MbdPmtContainer>(topNode,"MbdPmtContainer");


  _mbdavgt[0] = 0;
  _mbdavgt[1] = 0;
  _mbdhit[0] = 0;
  _mbdhit[1] = 0;
  
  if(mbdpmt)
    {
      for(int i=0; i<128; ++i)
	{
	  MbdPmtHit* pmt = mbdpmt->get_pmt(i);
	  if(pmt)
	    {
	      if(pmt->get_q() > 0.4)
		{
		  ++_mbdhit[i/64];
		}
	    }
	}
    }
  else
    {
      if(_isdat)
	{
	  cout << "no MBD PMTs!!!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
    }

  MbdOut * mbdout = static_cast<MbdOut*>(findNode::getClass<MbdOut>(topNode,"MbdOut"));
  if(mbdout){
    _mbdavgt[1] = mbdout->get_time(0); // south side
    _mbdavgt[0] = mbdout->get_time(1); // north side
  }

  if(std::isnan(_mbdavgt[0])) _mbdavgt[0] = -9999;
  if(std::isnan(_mbdavgt[1])) _mbdavgt[1] = -9999;

  float maxjetpt = 0;
  float subjetpt = 0;
  float maxjete = 0;
  float subjete = 0;
  float maxjetphi = 0;
  float subjetphi = 0;
  
  for(int i=0; i<_jet_n; ++i)
    {
      if(_jet_pt[i] > maxjetpt)
	{
	  subjetpt = maxjetpt;
	  subjete = maxjete;
	  subjetphi = maxjetphi;
	  maxjetpt = _jet_pt[i];
	  maxjete = _jet_et[i];
	  maxjetphi = _jet_phi[i];
	}
      else if(_jet_pt[i] > subjetpt)
	{
	  subjetpt = _jet_pt[i];
	  subjete = _jet_et[i];
	  subjetphi = _jet_phi[i];
	}

    }
  
  float dphi = abs(maxjetphi - subjetphi);
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  if(subjete > 0.3*7) _isdijet = 1;
  else _isdijet = 0;
  bool dPhiCut = (dphi < 3*M_PI/4 || !_isdijet || maxjete*0.3 > subjete);
  bool loETCut = _frcem > 0.9 || _frcem < 0.1 || _frcoh < 0.1 || _frcoh > 0.9 || (1.-_frcem-_frcoh) > 0.9;// ((_frcem < 0.1) && (_jet_ET > (50*_frcem+20))) && (dPhiCut || !_isdijet);
  int failsall = -1;//fullCut?1:0;
  if(!loETCut && !dPhiCut)
    {
      failsall = 2;
    }
  else if(!loETCut)
    {
      failsall = 1;
    }
  else if(!dPhiCut)
    {
      failsall = 0;
    }

  _failscut = failsall;

  if(_dotruthpar)
    {
      _truthparenergy.clear();
      _truthpareta.clear();
      _truthparphi.clear();
      _truthparpt.clear();
      _truthparid.clear();
      PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
      if(!truthinfo)
	{
	  std::cout
	    << "MyJetAnalysis::process_event - Error can not find DST Truth Info node "
	    << "G4TruthInfo" << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      // Get the primary particle range
      PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();
      
      // Loop over the G4 truth (stable) particles
      for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	   iter != range.second;
	   ++iter)
	{
	  // Get this truth particle
	  const PHG4Particle *truth = iter->second;
	  
	  /// Get this particles momentum, etc.
	  float m_truthpx = truth->get_px();
	  float m_truthpy = truth->get_py();
	  float m_truthpz = truth->get_pz();
	  float m_truthpt_temp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);
	  if(m_truthpt_temp < 0.5) continue;
	  // float m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
	  float m_tempe = truth->get_e();
	  
	  _truthparenergy.push_back(m_tempe);
	  _truthparphi.push_back(atan(m_truthpy / m_truthpx));
	  _truthparpt.push_back(m_truthpt_temp);
	  float m_tempeta= atanh(m_truthpz / m_tempe);
	  /// Check for nans
	  if (!std::isfinite(m_tempeta))
	    {
	      m_tempeta = -99;
	    }
	  _truthpareta.push_back(m_tempeta);
	  _truthparid.push_back(truth->get_pid());	      
	}
    }
  
  
  if(maxjetpt > _minjetthresh || !_isdat || isjettrig || isjmbtrig)
    {
      jet_tree->Fill();
    }
  if(_calib_jet_n != _jet_n) cout << "calib_jet_n != jet_n!!!" << _calib_jet_n << " != " << _jet_n << endl;
  if(_debug > 3) cout << "end event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int Chi2checker::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "Chi2checker::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Chi2checker::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "Chi2checker::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Chi2checker::End(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "Chi2checker::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  cout << "ending run" << endl;
  if(jet_tree->GetEntries() > 0)
    {
      //if(!_doall60 && _isdat) _f = new TFile(_filename.c_str(), "RECREATE");
      cout << "file created" << endl;
      //if(!_doall60 && _isdat) _f->cd();
      cout << "cded to file" << endl;
      jet_tree->SetDirectory(_f);
      _f->cd();
      cout << "tree set to directory of file" << endl;
      jet_tree->Write();
      cout << "tree written" << endl;
      _f->Write();
      cout << "file written" << endl;
      _f->Close();
      cout << "file closed" << endl;
    }
  if(_dowf)
    {
      if(_wft->GetEntries() > 0)
	{
	  if(!_doall60) _wff = new TFile(_wfilename.c_str(),"RECREATE");
	  if(!_doall60) _wff->cd();
	  if(!_doall60) _wft->SetDirectory(_wff);
	  _wft->Write();
	  _wff->Write();
	  _wff->Close();
	}
    }
  cout << "file saved if necessary" << endl;
  //delete jet_tree;
  //delete _f;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Chi2checker::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "Chi2checker::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void Chi2checker::Print(const std::string &what) const
{
  std::cout << "Chi2checker::Print(const std::string &what) const Printing info for " << what << std::endl;
}
