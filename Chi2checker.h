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
#include <phparameter/PHParameters.h>
#include <globalvertex/GlobalVertex.h>
class PHCompositeNode;
class CentralityInfo;
class Chi2checker : public SubsysReco
{
 public:

  Chi2checker(const std::string &filename = "/sphenix/user/jocl/projects/run2024_earlydata/run/output/temphists/debug.root", const std::string &name = "Chi2checker", const int debug = 0, const std::string &wfilename = "", const int dowf = 0, const int doall60 = 1);

  virtual ~Chi2checker();

  float bintoeta_hc(int etabin);

  float bintophi_hc(int phibin);

  void drawCalo(TowerInfoContainer** towers, float* jet_e, float* jet_et, float* jet_ph, int jet_n, float jet_ecc, float jet_lfrac, RawTowerGeomContainer** geom, float zvtx, int failscut, int runnum, int evtnum, float frcoh, float frcem, float emtot, float ohtot, float maxJetE);
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
  int _doall60;
  bool _printedPhi = false;
  int cancount = 0;
  PHParameters _cutParams;
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
  int _failscut;
  int _runnum;
  int _evtnum;
  float _jet_et[100];
  float _jet_pt[100];
  float _jet_t[100];
  float _jet_eta[100];
  float _jet_phi[100];
  float _alljetfrcem[100];
  float _alljetfrcoh[100];
  float _emtow[96][256];
  float _ihtow[24][64];
  float _ohtow[24][64];
  int _isbadem[96][256];
  int _isbadih[24][64];
  int _isbadoh[24][64];
  int _ishotem[96][256];
  int _ishotih[24][64];
  int _ishotoh[24][64];
  int _nocalem[96][256];
  int _nocalih[24][64];
  int _nocaloh[24][64];
  float _jconem[24][64];
  float _jconih[24][64];
  float _jconoh[24][64];
  float _chi2em[96][256];
  float _chi2ih[24][64];
  float _chi2oh[24][64];
  float _dPhi2pc[1000];
  float _dEta2pc[1000];
  float _dPhiLayer[100];
  float _l2pcEta;
  int _nprocessed;
  long long unsigned int _prevraw22;
  long long unsigned int _prevraw18;
  long long unsigned int _prevlive22;
  long long unsigned int _prevlive18;
  int _isbadlive;
  /*
  float _emLayerJetPhi[10];
  float _ohLayerJetPhi[10];
  float _emLayerJetEta[10];
  float _ohLayerJetEta[10];
  float _emLayerJetET[10];
  float _ohLayerJetET[10];
  int _nLayerEm;
  int _nLayerOh;
  */
  int _n2pc;
  GlobalVertex::VTXTYPE _vtxtype = GlobalVertex::MBD;

  int _dowf;
  std::string _wfilename;
  TTree* _wft;
  TFile* _wff;
  unsigned int _emwf[96][256][12];
  unsigned int _ihwf[24][64][12];
  unsigned int _ohwf[24][64][12];
  int _emieta[96][256];
  int _ihieta[24][64];
  int _ohieta[24][64];

  int _emiphi[96][256];
  int _ihiphi[24][64];
  int _ohiphi[24][64];

  float _emt[96][256];
  float _iht[24][64];
  float _oht[24][64];
  
  unsigned int _mbdhit[2];
  float _mbdavgt[2];

  float _emadcfit[96][256];
  float _ihadcfit[24][64];
  float _ohadcfit[24][64];
};

#endif // CHI2TREEMAKER
