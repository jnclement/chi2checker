#ifndef CHI2TREEMAKER_H
#define CHI2TREEMAKER_H
#include <fun4all/SubsysReco.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TH2D.h"
#include <phparameter/PHParameters.h>
#include <globalvertex/GlobalVertex.h>
class PHCompositeNode;
class CentralityInfo;
class Chi2checker : public SubsysReco
{
 public:

  Chi2checker(const std::string &filename = "/sphenix/user/jocl/projects/run2024_earlydata/run/output/temphists/debug.root", const std::string &name = "Chi2checker", const int debug = 0, const std::string &wfilename = "", const int dowf = 0, const bool isdat = 1, const int doall60 = 1, const int dotruthpar = 0);

  virtual ~Chi2checker();

  float bintoeta_hc(int etabin);

  float bintophi_hc(int phibin);

  float getEtaFromBin(int binEta);

  float getPhiFromBin(int binPhi);

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int fill_jet_quantities(PHCompositeNode* topNode, std::string jet_nodename, int& njet, int& calib_njet, float* jet_pt, float* jet_pt_calib, float* jet_eta, float* jet_phi, float* jet_e, float* jet_et, float* jet_t, float* jet_frcem, float* jet_frcoh);

  int fill_tjet_quantities(PHCompositeNode* topNode, std::string jet_nodename, int& njet, float* jet_pt, float* jet_eta, float* jet_phi, float* jet_e);

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void set_minjetthresh(float newthresh)
  {
    _minjetthresh = newthresh;
  }
 
 private:
  float _minjetthresh{7};
  int _dotruthpar;
  bool _isdat;
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
  static const int nx = 24;
  static const int ny = 64;
  static const int nt = 1536;
  static const int nemx = 96;
  static const int nemy = 256;
  static const int nemt = 24576;
  int _nx = 24;
  int _ny = 64;
  int _nemx = 96;
  int _nemy = 256;
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
  float _tzvtx;
  long long unsigned int _triggervec;
  unsigned int _bbfqavec;
  unsigned int _elmbgvec;
  int _mbevt;
  int _failscut;
  int _runnum;
  int _evtnum;


  int _r02_jet_n;
  int _calib_r02_jet_n;
  float _r02_jet_et[100];
  float _r02_jet_etrans[100];
  float _r02_jet_pt[100];
  float _r02_jet_pt_calib[100];
  float _r02_jet_t[100];
  float _r02_jet_t_em[100];
  float _r02_jet_t_ih[100];
  float _r02_jet_t_oh[100];
  float _r02_jet_eta[100];
  float _r02_jet_phi[100];

  int _r02_tjet_n;
  float _r02_tjet_e[100];
  float _r02_tjet_pt[100];
  float _r02_tjet_eta[100];
  float _r02_tjet_phi[100];
  
  float _allr02_jetfrcem[100];
  float _allr02_jetfrcoh[100];



  int _r03_jet_n;
  int _calib_r03_jet_n;
  float _r03_jet_et[100];
  float _r03_jet_etrans[100];
  float _r03_jet_pt[100];
  float _r03_jet_pt_calib[100];
  float _r03_jet_t[100];
  float _r03_jet_t_em[100];
  float _r03_jet_t_ih[100];
  float _r03_jet_t_oh[100];
  float _r03_jet_eta[100];
  float _r03_jet_phi[100];

  int _r03_tjet_n;
  float _r03_tjet_e[100];
  float _r03_tjet_pt[100];
  float _r03_tjet_eta[100];
  float _r03_tjet_phi[100];
  
  float _allr03_jetfrcem[100];
  float _allr03_jetfrcoh[100];



  int _jet_n;
  int _calib_jet_n;
  float _jet_et[100];
  float _jet_etrans[100];
  float _jet_pt[100];
  float _jet_pt_calib[100];
  float _jet_t[100];
  float _jet_t_em[100];
  float _jet_t_ih[100];
  float _jet_t_oh[100];
  float _jet_eta[100];
  float _jet_phi[100];

  int _tjet_n;
  float _tjet_e[100];
  float _tjet_pt[100];
  float _tjet_eta[100];
  float _tjet_phi[100];
  
  float _alljetfrcem[100];
  float _alljetfrcoh[100];


  int _r06_jet_n;
  int _calib_r06_jet_n;
  float _r06_jet_et[100];
  float _r06_jet_etrans[100];
  float _r06_jet_pt[100];
  float _r06_jet_pt_calib[100];
  float _r06_jet_t[100];
  float _r06_jet_t_em[100];
  float _r06_jet_t_ih[100];
  float _r06_jet_t_oh[100];
  float _r06_jet_eta[100];
  float _r06_jet_phi[100];

  int _r06_tjet_n;
  float _r06_tjet_e[100];
  float _r06_tjet_pt[100];
  float _r06_tjet_eta[100];
  float _r06_tjet_phi[100];
  
  float _allr06_jetfrcem[100];
  float _allr06_jetfrcoh[100];




  int _r08_jet_n;
  int _calib_r08_jet_n;
  float _r08_jet_et[100];
  float _r08_jet_etrans[100];
  float _r08_jet_pt[100];
  float _r08_jet_pt_calib[100];
  float _r08_jet_t[100];
  float _r08_jet_t_em[100];
  float _r08_jet_t_ih[100];
  float _r08_jet_t_oh[100];
  float _r08_jet_eta[100];
  float _r08_jet_phi[100];

  int _r08_tjet_n;
  float _r08_tjet_e[100];
  float _r08_tjet_pt[100];
  float _r08_tjet_eta[100];
  float _r08_tjet_phi[100];
  
  float _allr08_jetfrcem[100];
  float _allr08_jetfrcoh[100];





  float _emtow[nemx][nemy];
  float _ihtow[nx][ny];
  float _ohtow[nx][ny];
  int _isbadem[nemx][nemy];
  int _isbadih[nx][ny];
  int _isbadoh[nx][ny];
  int _ishotem[nemx][nemy];
  int _ishotih[nx][ny];
  int _ishotoh[nx][ny];
  int _nocalem[nemx][nemy];
  int _nocalih[nx][ny];
  int _nocaloh[nx][ny];
  float _jconem[nx][ny];
  float _jconih[nx][ny];
  float _jconoh[nx][ny];
  float _chi2em[nemx][nemy];
  float _chi2ih[nx][ny];
  float _chi2oh[nx][ny];
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
  unsigned int _emwf[nemx][nemy][12];
  unsigned int _ihwf[nx][ny][12];
  unsigned int _ohwf[nx][ny][12];
  int _emieta[nemx][nemy];
  int _ihieta[nx][ny];
  int _ohieta[nx][ny];

  int _emiphi[nemx][nemy];
  int _ihiphi[nx][ny];
  int _ohiphi[nx][ny];

  float _emt[nemx][nemy];
  float _iht[nx][ny];
  float _oht[nx][ny];
  
  unsigned int _mbdhit[2];
  float _mbdavgt[2];

  float _emadcfit[nemx][nemy];
  float _ihadcfit[nx][ny];
  float _ohadcfit[nx][ny];

  std::vector<float> _truthparenergy;
  std::vector<float> _truthpareta;
  std::vector<float> _truthparphi;
  std::vector<float> _truthparpt;
  std::vector<int> _truthparid;
};

#endif // CHI2TREEMAKER
