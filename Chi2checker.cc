#include <iomanip>
#include "Chi2checker.h"
#include <ffaobjects/EventHeaderv1.h>
#include <calobase/RawCluster.h>
#include <calobase/RawClusterContainer.h>
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
#include <g4main/PHG4VtxPoint.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <mbd/MbdPmtContainer.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phhepmc/PHHepMCGenEvent.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <phhepmc/PHHepMCGenEventMap.h>
#include <ffaobjects/EventHeaderv1.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <HepMC/GenEvent.h>
#include <mbd/MbdPmtHit.h>
#include <jetbackground/TowerBackgroundv1.h>
#include <cmath>
#include <mbd/MbdOut.h>
#include <TLorentzVector.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <g4centrality/PHG4CentralityReco.h>
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

#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>

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
#include <ffarawobjects/Gl1Packetv2.h>
using namespace std;

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
}

void sqrt_snn_text(float xp = 0.7, float yp = 0.8, bool isRightAlign=0, double textsize = 0.04)
{
  drawText("#sqrt{S_{NN}} = 200 GeV",xp,yp,isRightAlign,kBlack,textsize);
}

void sphenixtext(float xpos = 0.7, float ypos = 0.96, int ra = 0, float textsize = 0.04)
{
  drawText("#bf{#it{sPHENIX}} Internal", xpos, ypos, ra, kBlack, textsize);
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

void Chi2checker::drawCalo(TowerInfoContainer** towers, float* jet_e, float* jet_et, float* jet_ph, int jet_n, float jet_ecc, float jet_lfrac, RawTowerGeomContainer** geom, float zvtx, int failscut)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int ncircle = 64;
  int ncol = 9;
  double red[ncol] = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[ncol] = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[ncol] = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  TColor::CreateGradientColorTable(ncol, stp, red, grn, blu, ncol);
  
  TCanvas* c = new TCanvas("","",1900,600);
  c->Divide(4,1,0,0);
  
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",48,-2.2,2.2,64,0,2*M_PI);
  TH2D* event_disrt[3];
  for(int i=0; i<3; ++i)
    {
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",48,-2.2,2.2,64,0,2*M_PI);
    }

  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(0.75);
      event_disrt[i]->GetYaxis()->SetTitleOffset(1);
      event_disrt[i]->GetZaxis()->SetTitle("Uncorrected Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(0.1,5);
      event_disrt[i]->GetXaxis()->SetNdivisions(4,kFALSE);
      event_disrt[i]->GetXaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetYaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetZaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetXaxis()->SetTitleOffset(2);
      event_disrt[i]->GetXaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetYaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetZaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetXaxis()->SetLabelOffset(0.02);
    }
  event_sum->GetXaxis()->SetTitle("Calo Sum #eta");
  event_sum->GetYaxis()->SetTitle("Calo Sum #phi");
  event_sum->GetZaxis()->SetTitleOffset(0.75);
  event_sum->GetYaxis()->SetTitleOffset(1);
  event_sum->GetZaxis()->SetTitle("Uncorrected Tower Energy [GeV]");
  event_sum->GetZaxis()->SetRangeUser(0.1,5);
  event_sum->GetXaxis()->SetNdivisions(4,kFALSE);
  event_sum->GetXaxis()->SetTitleSize(0.04);
  event_sum->GetYaxis()->SetTitleSize(0.04);
  event_sum->GetZaxis()->SetTitleSize(0.04);
  event_sum->GetXaxis()->SetTitleOffset(2);
  event_sum->GetXaxis()->SetLabelSize(0.04);
  event_sum->GetYaxis()->SetLabelSize(0.04);
  event_sum->GetZaxis()->SetLabelSize(0.04);
  event_sum->GetXaxis()->SetLabelOffset(0.02);
  
  event_sum->Reset();
  for(int j=0; j<3; ++j)
    {
      event_disrt[j]->Reset();
      for(int k=0; k<1536; ++k)
	{
	  TowerInfo* tower = towers[j]->get_tower_at_channel(k);
	  int key = towers[j]->encode_key(k);
	  float eta = 0;
	  float phi = 0;
	  if(j==1)
	    {
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[j]->getTowerEtaBin(key), towers[j]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[j]->get_tower_geometry(geomkey); //encode tower geometry
	      float newx = tower_geom->get_center_x();
	      float newy = tower_geom->get_center_y();
	      float newz = tower_geom->get_center_z() - zvtx;
	      eta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
			  
	      phi = tower_geom->get_phi()+M_PI;//bintophi_hc(towers[j]->getTowerPhiBin(key))+0.048;
	    }
	  else if(j==2)
	    {
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towers[j]->getTowerEtaBin(key), towers[j]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[j]->get_tower_geometry(geomkey); //encode tower geometry
	      float newx = tower_geom->get_center_x();
	      float newy = tower_geom->get_center_y();
	      float newz = tower_geom->get_center_z() - zvtx;
	      eta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
	      phi = tower_geom->get_phi()+M_PI;//
	    }
	  else
	    {
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[j]->getTowerEtaBin(key), towers[j]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry
	      float newx = tower_geom->get_center_x();
	      float newy = tower_geom->get_center_y();
	      float newz = tower_geom->get_center_z() - zvtx;
	      eta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
	      phi = tower_geom->get_phi()+M_PI;//bintophi_em(towers[j]->getTowerPhiBin(key))+0.012;
	    }
	  event_disrt[j]->Fill(eta,phi,tower->get_energy());
	  event_sum->Fill(eta,phi,tower->get_energy());
	}


      c->cd(j+1);
      gPad->SetLogz();
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.2);                                                               
      event_disrt[j]->Draw("COLZ");

    }
  c->cd(4);
  gPad->SetLogz();                                                                   
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  event_sum->Draw("COLZ");
  c->cd(0);
  sphenixtext(0.96,0.96,1,0.04);

  std::stringstream ecc_stream;
  ecc_stream << std::fixed << std::setprecision(3) << jet_ecc;
  std::string ecc_string = ecc_stream.str();
  drawText(("#epsilon_{max jet}="+ecc_string).c_str(),0.1,0.95);

  std::stringstream lfrac_stream;
  lfrac_stream << std::fixed << std::setprecision(3) << jet_lfrac;
  std::string lfrac_string = lfrac_stream.str();
  drawText(("E_{max layer,jet}/E_{max jet}="+lfrac_string).c_str(),0.4,0.95);

  std::stringstream z_stream;
  z_stream << std::fixed << std::setprecision(1) << zvtx;
  std::string z_string = z_stream.str();
  drawText(("z_{vtx}="+z_string).c_str(),0.25,0.95);
  
  std::string fails = failscut?"Fails OHCal cut":"Passes OHCal cut";
  drawText(fails.c_str(),0.7,0.95);

  c->cd(4);
  float maxJetEta = 0;
  float maxJetPhi = 0;
  float maxJetE = 0;
  for(int k=0; k<jet_n; ++k)
    {
      if(maxJetE < jet_e[k])
	{
	  maxJetE = jet_e[k];
	  maxJetEta = jet_et[k];
	  maxJetPhi = jet_ph[k];
	}
      for(int l=0; l<ncircle; ++l)
	{
	  float eta = jet_et[k]+0.4*cos(2*l*M_PI/ncircle);
	  float phi = jet_ph[k]+0.4*sin(2*l*M_PI/ncircle);
	  phi += M_PI;
	  if(eta > 2.2 || eta < -2.2) continue;
	  if(phi > 2*M_PI) phi -= 2*M_PI;
	  if(phi < 0) phi += 2*M_PI;
	  TMarker* circlemarker = new TMarker(eta,phi,20);
	  circlemarker->SetMarkerSize(0.3);
	  circlemarker->SetMarkerColor(kBlue);
	  circlemarker->Draw();
	}
      std::stringstream e_stream;
      e_stream << std::fixed << std::setprecision(2) << jet_e[k]/cosh(jet_et[k]);
      std::string e_string = e_stream.str();
      drawText((e_string+" GeV").c_str(),jet_et[k],jet_ph[k]+(jet_ph[k]+M_PI>3.84?-0.53:0.43)+M_PI,(jet_et[k]>0?1:0),kBlack,0.04,42,false);

    }
  c->SaveAs(("./output/smg/candidate_"+_name+"_supersuperhighE_eccentricityfinder_"+to_string(cancount)+".png").c_str());
  cout << "Saved" << endl;
  cancount++;

}

//____________________________________________________________________________..
Chi2checker::Chi2checker(const std::string &filename, const std::string &name, const int debug):
  SubsysReco("test")//).c_str())
{
  _name = name;
  _debug = debug;
  _filename = filename;
  _f = new TFile(filename.c_str(), "RECREATE");
  jet_tree = new TTree("jet_tree","a persevering date tree");
  /*
  for(int h=0; h<3; ++h)
    {
      if(h<2) h1_dphi[h] = new TH1D(("h1_dphi"+to_string(h)).c_str(),"",32,0,M_PI);
      for(int i=0; i<6; ++i)
	{
	  h2_ecc_layer[h][i] = new TH2D(("h2_ecc_layer"+to_string(h)+"_"+to_string(i)).c_str(),"",120,0,1.2,100,0,2);
	  h2_ecc_angle[h][i] = new TH2D(("h2_ecc_angle"+to_string(h)+"_"+to_string(i)).c_str(),"",120,0,1.2,32,-M_PI/2,M_PI/2);
	}
      h2_ecc_E[h] = new TH2D(("h2_ecc_E"+to_string(h)).c_str(),"",120,0,1.2,60,0,60);
      h2_g20_ecc_angle[h] = new TH2D(("h2_g20_ecc_angle"+to_string(h)).c_str(),"",120,0,1.2,32,-M_PI/2,M_PI/2);
      h2_g20_ecc_frcoh[h] = new TH2D(("h2_g20_ecc_frcoh"+to_string(h)).c_str(),"",120,0,1.2,100,0,2);
      h2_g20_ecc_frcem[h] = new TH2D(("h2_g20_ecc_frcem"+to_string(h)).c_str(),"",120,0,1.2,100,0,2);
      h1_jet_eta[h] = new TH1D(("h1_jet_eta"+to_string(h)).c_str(),"",48,-2.2,2.2);
      h1_jet_phi[h] = new TH1D(("h1_jet_phi"+to_string(h)).c_str(),"",64,0,2*M_PI);
      }
  */
}

//____________________________________________________________________________..
Chi2checker::~Chi2checker()
{

}

//____________________________________________________________________________..
int Chi2checker::Init(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Begin init: " << endl;
  jet_tree->Branch("nBadChi2",&_nBadChi2,"nBadChi2/I");
  jet_tree->Branch("ecc",&_eccentricity,"ecc/F");
  jet_tree->Branch("frcoh",&_frcoh,"frcoh/F");
  jet_tree->Branch("frcem",&_frcem,"frcem/F");
  jet_tree->Branch("eta",&_eta,"eta/F");
  jet_tree->Branch("phi",&_phi,"phi/F");
  jet_tree->Branch("jet_ET",&_jet_ET,"jet_ET/F");
  jet_tree->Branch("dphi",&_dphi,"dphi/F");
  jet_tree->Branch("isdijet",&_isdijet,"isdijet/I");
  jet_tree->Branch("subjet_ET",&_subjet_ET,"subjet_ET/F");
  jet_tree->Branch("zvtx",&_zvtx,"zvtx/F");
  //jet_tree->Branch("jetcompE",_jetcompE,"jetcompE[3][512]/F");
  //jet_tree->Branch("jetcompEta",_jetcompEta,"jetcompEta[3][512]/F");
  //jet_tree->Branch("jetcompPhi",_jetcompPhi,"jetcompPhi[3][512]/F");
  jet_tree->Branch("maxTowChi2",_maxTowChi2,"maxTowChi2[3]/F");
  jet_tree->Branch("maxTowE",&_maxTowE,"maxTowE/F");
  jet_tree->Branch("subTowE",&_subTowE,"subTowE/F");
  jet_tree->Branch("maxTowDiff",&_maxTowDiff,"maxTowDiff/F");
  jet_tree->Branch("maxETowChi2",&_maxETowChi2,"maxETowChi2/F");
  jet_tree->Branch("ohPhiBinMaxFrac",&_ohPhiBinMaxFrac,"ohPhiBinMaxFrac/F");
  jet_tree->Branch("maxETowIsZS",&_maxETowIsZS,"maxETowIsZS/I");
  jet_tree->Branch("maxETowChi2Det",&_maxETowChi2Det,"maxETowChi2Det/I");
  jet_tree->Branch("triggervec",&_triggervec,"triggervec/g");
  jet_tree->Branch("bbfqavec",&_bbfqavec,"bbfqavec/i");
  jet_tree->Branch("elmbgvec",&_elmbgvec,"elmbgvec/i");
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int Chi2checker::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

float Chi2checker::getEtaFromBin(int binEta)
{
  return ((2.2*binEta)/nx)-1.1;
}

float Chi2checker::getPhiFromBin(int binPhi)
{
  return ((2*M_PI*binPhi)/ny);
}

float getEtaFromBinEM(int binEta)
{
  return ((2.2*binEta)/96)-1.1;
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

int Chi2checker::process_event(PHCompositeNode *topNode)
{

  if(_debug > 1) cout << endl << endl << endl << "Chi2checker: Beginning event processing" << endl;
  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r04");

  if(_debug > 2) cout << towers[0] << " " << towers[1] << " " << towers[2] << endl;

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");

  //TowerInfoContainer* emcrt = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");

  MbdVertexMap* mbdvtxmap = findNode::getClass<MbdVertexMapv1>(topNode, "MbdVertexMap");
  GlobalVertexMap* gvtxmap = findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");
  Gl1Packetv2* gl1 = gl1 = findNode::getClass<Gl1Packetv2>(topNode, "GL1Packet");
  if(!gl1)
    {
      cout << "No trigger info!" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  if(_debug > 1) cout << "Getting gl1 trigger vector from: " << gl1 << endl;
  _triggervec = gl1->getScaledVector();

  float zvtx = NAN;
  if(!gvtxmap)
    {
      for(auto iter = mbdvtxmap->begin(); iter != mbdvtxmap->end(); ++iter)
        {
          MbdVertex* mbdvtx = iter->second;
          zvtx = mbdvtx->get_z();
          break;
        }
    }
  else
    {
      auto iter = gvtxmap->begin();
      while(iter != gvtxmap->end())
        {
          GlobalVertex* gvtx = iter->second;
          zvtx = gvtx->get_z();
          iter++;
          break;
        }
    }
  if(std::isnan(zvtx) || abs(zvtx) > 150)
    {
      /*
      if(std::isnan(zvtx))
	{
	  cout << "no zvtx found! " << n_nozvtx << endl;
	  ++n_nozvtx;
	}
      */
      return Fun4AllReturnCodes::EVENT_OK;
    }
  _zvtx = zvtx;
  //float fracEM = 0;
  //float fracOH = 0;
  _maxTowE = 0;
  _subTowE = 0;
  _maxTowDiff = 0;
  _nBadChi2 = 0;
  float maxJetE = 0;
  float maxJetEta = 0;
  float maxJetPhi = 0;
  float subJetE = 0;
  float subJetPhi = 0;
  float maxLayerE[2] = {0};
  float maxEoverTot = 0;
  float SumEdEtadPhi = 0;
  float SumEdEtadPhiDiff = 0;
  float v1 = 0;
  float v2 = 0;
  float theta = 0;
  int jet_n = 0;
  float jet_e[10];
  float jet_eta[10];
  float jet_phi[10];
  float jet_ecc = -1;
  float jet_lfrac = -1;
  //float Etot = 0;
  //bool isPerimeter[nx][ny] = {0};
  //float calE[3][nx][ny] = {0};
  MinimumBiasInfov1* mb3 = findNode::getClass<MinimumBiasInfov1>(topNode,"mbc_bkgd3");
  //int maxCalN = -1;
  //float maxCaleE = 0;
  //float eTot = 0;

  float eccentricity = 0;
  //int failscut = 0;
  if(jets)
    {
      int tocheck = jets->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *jet = jets->get_jet(i);
          if(jet)
            {
	      //float Etot = 0;
	      if(_debug > 5) cout << "getting jet E/eta" << endl;
	      float testJetE = jet->get_e()/cosh(jet->get_eta());
	      float squareDif = 0;
	      float squareSum = 0;
	      float testJetPhi = jet->get_phi();
	      if(_debug > 5) cout << "jet E/eta: " << testJetE  << " " << jet->get_eta() << endl;
	      if(testJetE < 8) continue;
	      if(_debug > 3) cout << "got a candidate jet" << endl;
	      jet_e[jet_n] = testJetE;
	      jet_eta[jet_n] = jet->get_eta();
	      //if(abs(jet_eta[jet_n]) > 0.9) continue;
	      jet_phi[jet_n] = testJetPhi;//(jet->get_phi()>0?jet->get_phi():jet->get_phi()+2*M_PI);
	      
	      if(testJetE > subJetE && testJetE < maxJetE)
		{
		  subJetE = testJetE;
		  subJetPhi = testJetPhi;
		}
	      if(_debug > 2) cout << "found a good jet!" << endl;
	      if(testJetE > maxJetE)
		{
		  if(maxJetE > subJetE)
		    {
		      subJetE = maxJetE;
		      subJetPhi = maxJetPhi;
		    }
		  int ncomp = 0;
		  _maxTowE = 0;
		  _subTowE = 0;
		  _maxTowDiff = 0;
		  _nBadChi2 = 0;
		  _maxETowChi2 = 0;
		  _maxETowIsZS = -1;
		  _maxETowChi2Det = -1;
		  eccentricity = 0;
		  maxLayerE[0] = 0;
		  maxLayerE[1] = 0;
		  SumEdEtadPhi = 0;
		  SumEdEtadPhiDiff = 0;
		  v1 = 0;
		  v2 = 0;
		  maxJetE = testJetE;
		  maxJetEta = jet->get_eta();
		  maxJetPhi = jet->get_phi();//(jet->get_phi()>0?jet->get_phi():jet->get_phi()+2*M_PI);
		  float sigEtaEta = 0;
		  float sigEtaPhi = 0;
		  float sigPhiPhi = 0;
		  if(_debug > 4) cout << "getting comp vec" << endl;
		  
		  int subcomp[3] = {0};
		  for(int i=0; i<3; ++i)
		    {
		      _maxTowChi2[i] = -1;
		      /*
		      for(int k=0; k<512; ++k)
			{
			  _jetcompE[i][k] = 0;
			  _jetcompEta[i][k] = 0;
			  _jetcompPhi[i][k] = 0;
			}
		      */
		    }
		  
		  int skipflag = 0;
		  for(auto comp: jet->get_comp_vec())
		    {
		      ++ncomp;
		      unsigned int channel = comp.second;
		      TowerInfo* tower;
		      //cout << "type: " << comp.first << endl;
		      if(comp.first == 5 || comp.first == 26)
			{
			  tower = towers[1]->get_tower_at_channel(channel);
			  float towerE = tower->get_energy();
			  float chi2 = tower->get_chi2();
			  //if(abs(chi2 - 0.08) < 0.00001 && towerE > 8) failscut = 2;
			  if(chi2 > _maxTowChi2[1]) _maxTowChi2[1] = chi2;
			  if(tower->get_isBadChi2()) _nBadChi2++;
			  //Etot += towerE;
			  int key = towers[1]->encode_key(channel);
			  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[1]->getTowerEtaBin(key), towers[1]->getTowerPhiBin(key));
			  if(_debug > 6) cout << "encoding tower geom" << endl;
			  RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry
			  if(_debug > 6) cout << "encoded" << endl;
			  float newx = tower_geom->get_center_x();
			  float newy = tower_geom->get_center_y();
			  float newz = tower_geom->get_center_z() - zvtx;
			  if(_debug > 6) cout << "got center values of towers" << endl;
			  float towerEta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
			  
			  float towerPhi = tower_geom->get_phi();//getPhiFromBinEM(towers[1]->getTowerPhiBin(key))+0.048;//tower_geom->get_phicenter();;
			  
			  /*
			  _jetcompE[1][subcomp[1]] = towerE;
			  _jetcompEta[1][subcomp[1]] = towerPhi;
			  _jetcompPhi[1][subcomp[1]] = towerEta;
			  subcomp[1]++;
			  */
			  towerE /= cosh(towerEta);
			  if(towerE > _maxTowE)
			    {
			      _subTowE = _maxTowE;
			      _maxTowE = towerE;
			      _maxETowChi2Det = 1;
			      _maxETowIsZS = tower->get_isZS();
			    }
			  if(towerE < 0) continue;
			  //towerE = sqrt(towerE);
			  float dPhi = towerPhi - maxJetPhi;
			  if(dPhi > M_PI) dPhi -= 2*M_PI;
			  if(dPhi < -M_PI) dPhi += 2*M_PI;
			  float dEta = towerEta - maxJetEta;
			  if(_debug > 4) cout << "IHCal" << endl;
			  if(_debug > 4) print_debug(maxJetEta, maxJetPhi, towerEta, towerPhi, dPhi, dEta);
			  sigEtaEta += towerE*dEta*dEta;
			  sigEtaPhi += towerE*dEta*dPhi;
			  sigPhiPhi += towerE*dPhi*dPhi;
			  SumEdEtadPhi += 2*towerE*dPhi*dEta;
			  SumEdEtadPhiDiff += towerE*(dPhi*dPhi-dEta*dEta);
			}
		      else if(comp.first == 7 || comp.first == 27)
			{
			  tower = towers[2]->get_tower_at_channel(channel);
			  float towerE = tower->get_energy();
			  
			  float chi2 = tower->get_chi2();
			  //if(abs(chi2 - 0.08) < 0.00001 && towerE > 8) failscut = 2;
			  if(chi2 > _maxTowChi2[2]) _maxTowChi2[2] = chi2;
			  if(tower->get_isBadChi2()) _nBadChi2++;
			  //fracOH += towerE;
			  //Etot += towerE;
			  int key = towers[2]->encode_key(channel);
			  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towers[2]->getTowerEtaBin(key), towers[2]->getTowerPhiBin(key));
			  RawTowerGeom *tower_geom = geom[2]->get_tower_geometry(geomkey); //encode tower geometry
			  float newx = tower_geom->get_center_x();
			  float newy = tower_geom->get_center_y();
			  float newz = tower_geom->get_center_z() - zvtx;
			  float towerEta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
			  float towerPhi = tower_geom->get_phi();//getPhiFromBinEM(towers[2]->getTowerPhiBin(key))+0.048;//tower_geom->get_phicenter(towers[2]->getTowerPhiBin(key));//getPhiFromBinEM()+0.012;
			  /*
			  _jetcompE[2][subcomp[2]] = towerE;
			  _jetcompEta[2][subcomp[2]] = towerPhi;
			  _jetcompPhi[2][subcomp[2]] = towerEta;
			  subcomp[2]++;
			  */
			  towerE /= cosh(towerEta);
			  if(towerE > _maxTowE)
			    {
			      _subTowE = _maxTowE;
			      _maxTowE = towerE;
			      _maxETowChi2Det = 2;
			      _maxETowIsZS = tower->get_isZS();
			    }
			  maxLayerE[1] += towerE;
			  if(towerE < 0) continue;
			  //towerE = sqrt(towerE);
			  float dPhi = towerPhi - maxJetPhi;
			  if(dPhi > M_PI) dPhi -= 2*M_PI;
			  if(dPhi < -M_PI) dPhi += 2*M_PI;
			  float dEta = towerEta - maxJetEta;
			  if(_debug > 4) cout << "OHCal" << endl;
			  if(_debug > 4) print_debug(maxJetEta, maxJetPhi, towerEta, towerPhi, dPhi, dEta);
			  sigEtaEta += towerE*dEta*dEta;
			  sigEtaPhi += towerE*dEta*dPhi;
			  sigPhiPhi += towerE*dPhi*dPhi;
			  SumEdEtadPhi += 2*towerE*dPhi*dEta;
			  SumEdEtadPhiDiff += towerE*(dPhi*dPhi-dEta*dEta);

			}
		      else if(comp.first == 13 || comp.first == 28 || comp.first == 25)
			{
			  tower = towers[0]->get_tower_at_channel(channel);
			  float towerE = tower->get_energy();
			  float chi2 = tower->get_chi2();
			  //if(abs(chi2 - 0.08) < 0.00001 && towerE > 8) failscut = 2;
			  if(chi2 > _maxTowChi2[0]) _maxTowChi2[0] = chi2;
      			  if(tower->get_isBadChi2()) _nBadChi2++;
			  //Etot += towerE;
			  //fracEM += towerE;
			  int key = towers[0]->encode_key(channel);
			  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[0]->getTowerEtaBin(key), towers[0]->getTowerPhiBin(key));
			  RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry
			  float newx = tower_geom->get_center_x();
			  float newy = tower_geom->get_center_y();
			  float newz = tower_geom->get_center_z() - zvtx;
			  float towerEta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
			  float towerPhi = tower_geom->get_phi();//getPhiFromBinEM(towers[0]->getTowerPhiBin(key))+0.012;//tower_geom->get_phicenter(towers[0]->getTowerPhiBin(key));//getPhiFromBinEM()+0.012;
			  /*
			  _jetcompE[0][subcomp[0]] = towerE;
			  _jetcompEta[0][subcomp[0]] = towerPhi;
			  _jetcompPhi[0][subcomp[0]] = towerEta;
			  subcomp[0]++;
			  */
			  towerE /= cosh(towerEta);
			  if(towerE > _maxTowE)
			    {
			      _subTowE = _maxTowE;
			      _maxTowE = towerE;
			      _maxETowChi2Det = 0;
			      _maxETowIsZS = tower->get_isZS();
			    }
			  maxLayerE[0] += towerE;
			  if(towerE < 0) continue;
			  //towerE = sqrt(towerE);
			  float dPhi = towerPhi - maxJetPhi;
			  if(dPhi > M_PI) dPhi -= 2*M_PI;
			  if(dPhi < -M_PI) dPhi += 2*M_PI;
			  float dEta = towerEta - maxJetEta;
			  if(_debug > 4) cout << "EMCal" << endl;
			  if(_debug > 4) print_debug(maxJetEta, maxJetPhi, towerEta, towerPhi, dPhi, dEta);
			  sigEtaEta += towerE*dEta*dEta;
			  sigEtaPhi += towerE*dEta*dPhi;
			  sigPhiPhi += towerE*dPhi*dPhi;
			  SumEdEtadPhi += 2*towerE*dPhi*dEta;
			  SumEdEtadPhiDiff += towerE*(dPhi*dPhi-dEta*dEta);
			}
		    }

		  float sqrtMaxJetE = sqrt(maxJetE);
		  sigEtaEta /= maxJetE;
		  sigEtaPhi /= maxJetE;
		  sigPhiPhi /= maxJetE;
		  float a = sigEtaEta;
		  float b = sigEtaPhi;
		  float c = sigEtaPhi;
		  float d = sigPhiPhi;
		  float lam1 = (a+d+sqrt((a+d)*(a+d)-4*(a*d-b*c)))/2;
		  float lam2 = (a+d-sqrt((a+d)*(a+d)-4*(a*d-b*c)))/2;
		  v1 /= ncomp;
		  v2 /= ncomp;
		  maxEoverTot = (maxLayerE[0] > maxLayerE[1]?maxLayerE[0]:maxLayerE[1])/maxJetE;
		  if(_debug > 3 && (lam1 < 0 || lam2 < 0)) cout << "lam1 or lam2 < 0, printing:" << lam1 << " " << lam2 << endl;
		  if(lam1 > lam2) eccentricity = 1-lam2/lam1;
		  else eccentricity = 1-lam1/lam2;
		  if(_debug > 3) cout << "ecc/layer: " << eccentricity << " " << maxEoverTot << endl;
		  //if(_debug > 3) cout << "ecc entries: " << h2_ecc_layer->GetEntries() << endl;
		}
	      ++jet_n;
	    }
          else
            {
              continue;
            }
	}
      _maxETowChi2 = mb3->getMyVal();

      /*
      if(_debug > 4) cout << "ended jets" << endl;
      if (gDirectory->Get("hcal_phi"))
	{
	  gDirectory->Delete("hcal_phi");
	  delete gDirectory->Get("hcal_phi");
	}
      TH1F* hcal_phi = new TH1F("hcal_phi", "", 64, -0.5, 63.5);
      hcal_phi->SetDirectory(nullptr);
      int size = towers[2]->size();
      for (int channel = 0; channel < size;channel++)
	{
	  TowerInfo *tower = towers[2]->get_tower_at_channel(channel);
	  float energy = tower->get_energy();
	  unsigned int towerkey = towers[2]->encode_key(channel);

	  int iphi = towers[2]->getTowerPhiBin(towerkey);
	  short good = (tower->get_isGood() ? 1:0);
	  if (!good) continue;
	  hcal_phi->Fill(iphi,energy);
	}

      int bmax = hcal_phi->GetMaximumBin();
      float binmax = hcal_phi->GetBinContent(bmax);

      //float binsum = binsup+binslo+binmax;                                        
      float binsum = hcal_phi->Integral();
      delete hcal_phi;
      hcal_phi = nullptr;
      if (_f->Get("hcal_phi"))
	{
	  _f->Delete("hcal_phi");
	  delete _f->Get("hcal_phi");
	}
      _ohPhiBinMaxFrac = binmax/binsum;
      */
      MinimumBiasInfov1* mbinfo = findNode::getClass<MinimumBiasInfov1>(topNode, "mbc_bkgd2");
      int failscut = 0;
      if(_debug > 3) cout << "mbinfo: " << mbinfo << endl; 
      if(mbinfo)
	{
	  //cout << "mbinfo exists - setting failscut to" << mbinfo->isAuAuMinimumBias() << endl;
	  _elmbgvec = mbinfo->getBkgdType(); //THIS IS ACTUALLY WHETHER OR NOT IT FAILS HANPU'S CUT//(binmax/binsum > 0.9)?1:0;
	}
      else
	{
	  if(_debug > 3) cout << "NO MBINFO NODE!" << endl;
	}

      MinimumBiasInfov1* mbinfo2 = findNode::getClass<MinimumBiasInfov1>(topNode, "mbc_bkgd");
      if(mbinfo2)
	{
	  _bbfqavec = mbinfo2->getBkgdType();
	}
      _maxTowDiff = _maxTowE - _subTowE;
      //fracEM /= Etot;
      //fracOH /= Etot;
      int fillnum = 0;
      int hfillnum = 2;
      float dphi = abs(maxJetPhi-subJetPhi);
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      //if(subJetE < 8) hfillnum=2;
      if(subJetE > 8) _isdijet = 1;
      else _isdijet = 0;
      /*
      if(failscut && _isdijet) _isdijet = 3;
      else if(failscut) _isdijet = 2;
      */
      //_isdijet += failscut;
      //{
      if(_debug > 4) cout << "assigning vars now" << endl;
      _eccentricity = eccentricity;
      _theta = theta;
      _frcoh = maxLayerE[1]/maxJetE;
      _frcem = maxLayerE[0]/maxJetE;
      _eta = maxJetEta;
      _phi = maxJetPhi;
      _jet_ET = maxJetE;
      _dphi = dphi;
      _subjet_ET = subJetE;
      if(maxJetE > 8) jet_tree->Fill();
      //}
      jet_ecc = eccentricity;
      jet_lfrac = maxEoverTot;
      
      if(maxJetE > 45)
	{
	  drawCalo(towers, jet_e, jet_eta, jet_phi, jet_n, jet_ecc, jet_lfrac, geom, zvtx, failscut);
	  cout << "drew calo" << endl;
	}
      
    }

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
  _f->cd();
  /*
  for(int h=0; h<3; ++h)
    {
      if(h<2) _f->WriteObject(h1_dphi[h],h1_dphi[h]->GetName());
      for(int i=0; i<6; ++i)
	{
	  _f->WriteObject(h2_ecc_layer[h][i],h2_ecc_layer[h][i]->GetName());
	  _f->WriteObject(h2_ecc_angle[h][i],h2_ecc_angle[h][i]->GetName());
	}
      _f->WriteObject(h2_ecc_E[h],h2_ecc_E[h]->GetName());
      _f->WriteObject(h2_g20_ecc_angle[h],h2_g20_ecc_angle[h]->GetName());
      _f->WriteObject(h2_g20_ecc_frcoh[h],h2_g20_ecc_frcoh[h]->GetName());
      _f->WriteObject(h2_g20_ecc_frcem[h],h2_g20_ecc_frcem[h]->GetName());
    }
  */
  jet_tree->Write();
  _f->Write();
  _f->Close();


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
