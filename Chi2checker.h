#ifndef CHI2TREEMAKER_H
#define CHI2TREEMAKER_H
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include <calobase/RawTowerGeomContainer.h>

class PHCompositeNode;
class CentralityInfo;
class Chi2checker : public SubsysReco
{
 public:

  Chi2checker(const std::string &filename = "/sphenix/user/jocl/projects/run2024_earlydata/run/output/temphists/debug.root", const std::string &name = "Chi2checker", const int debug = 0);

  virtual ~Chi2checker();

  float bintoeta_hc(int etabin);

  float bintophi_hc(int phibin);

  void drawCalo(TowerInfoContainer** towers, float* jet_e, float* jet_et, float* jet_ph, int jet_n, float jet_ecc, float jet_lfrac, RawTowerGeomContainer** geom, float zvtx, int failscut);
  float getEtaFromBin(int binEta);

  float getPhiFromBin(int binPhi);

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;


 private:
  int cancount = 0;
  TTree* jet_tree;
  TTree* mbtree;
  TH2D* h2_ecc_layer[3][6];
  TH2D* h2_ecc_angle[3][6];
  TH2D* h2_ecc_E[3];
  TH2D* h2_g20_ecc_angle[3];
  TH2D* h2_g20_ecc_frcoh[3];
  TH2D* h2_g20_ecc_frcem[3];
  TH1D* h1_dphi[2];
  TH1D* h1_jet_eta[3];
  TH1D* h1_jet_phi[3];
  int _debug;
  std::string _filename;
  TFile* _f;
  const int nx = 24;
  const int ny = 64;
  const int nt = 1536;
  std::string _name;
  float _eccentricity;
  float _theta;
  float _frcoh;
  float _frcem;
  float _eta;
  float _phi;
  float _jet_ET;
  float _dphi;
  float _subjet_ET;
  int _isdijet;
  float _jetcompE[3][512];
  float _jetcompPhi[3][512];
  float _jetcompEta[3][512];
  float _maxTowChi2[3];
  float _maxTowE;
  float _subTowE;
  float _maxTowDiff;
  int n_nozvtx = 0;
  int _nBadChi2;
  float _maxETowChi2;
  int _maxETowChi2Det;
  int _maxETowIsZS;
  float _ohPhiBinMaxFrac;
  float _zvtx;
  long long unsigned int _triggervec;
  unsigned int _bbfqavec;
  unsigned int _elmbgvec;
  int _mbevt;
  int _jet_n;
  float _jet_et[10];
  float _jet_eta[10];
  float _jet_phi[10];
};

#endif // CHI2TREEMAKER
