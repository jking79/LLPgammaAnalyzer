#ifndef ReducedNtuple_h
#define ReducedNtuple_h

#include "NtupleBase.hh"
#include "RestFrames/RestFrames.hh"

template class std::vector<std::vector<int> >;

using namespace RestFrames;

template <class Base>
class ReducedNtuple : public NtupleBase<Base> {

public:
  ReducedNtuple(TTree* tree = 0);
  virtual ~ReducedNtuple();

private:
  TTree* InitOutputTree(const string& sample);
  void FillOutputTree(TTree* tree);

  void ClearVariables();

  // common variables for output tree
  double m_weight;
  
  double m_MET;
  double m_MET_phi;

  double m_genMET;
  double m_genMET_phi;

  double m_HT;

  int m_Nele;
  int m_Nmu;
  
  int m_Nlep;
  vector<double> m_PT_lep;
  vector<double> m_Eta_lep;
  vector<double> m_Phi_lep;
  vector<double> m_M_lep;
  vector<int>    m_Charge_lep;
  vector<int>    m_PDGID_lep;
  vector<int>    m_RelIso_lep;
  vector<int>    m_MiniIso_lep;
  vector<int>    m_ID_lep;
  vector<int>    m_Index_lep;

  int m_Njet;
  int m_Nbjet;
  vector<double> m_PT_jet;
  vector<double> m_Eta_jet;
  vector<double> m_Phi_jet;
  vector<double> m_M_jet;
  vector<double> m_Btag_jet;
  vector<double> m_Flavor_jet;

  int m_genNele;
  int m_genNmu;

  int m_genNlep;
  vector<double> m_genPT_lep;
  vector<double> m_genEta_lep;
  vector<double> m_genPhi_lep;
  vector<double> m_genM_lep;
  vector<int>    m_genCharge_lep;
  vector<int>    m_genPDGID_lep;
  vector<int>    m_genMomPDGID_lep;
  vector<int>    m_genIndex_lep;

  int m_genNnu;
  vector<double> m_genPT_nu;
  vector<double> m_genEta_nu;
  vector<double> m_genPhi_nu;
  vector<int>    m_genPDGID_nu;
  vector<int>    m_genMomPDGID_nu;
  
  int m_genNboson;
  vector<double> m_genPT_boson;
  vector<double> m_genEta_boson;
  vector<double> m_genPhi_boson;
  vector<double> m_genM_boson;
  vector<int>    m_genPDGID_boson;
  vector<int>    m_genMomPDGID_boson;
  
  int m_genNsusy;
  vector<double> m_genPT_susy;
  vector<double> m_genEta_susy;
  vector<double> m_genPhi_susy;
  vector<double> m_genM_susy;
  vector<int>    m_genPDGID_susy;
  vector<int>    m_genMomPDGID_susy;

  //////////////////////
  // derived observables
  //////////////////////

  
  
  // Sparticle pair-production trees analysis
  vector<int> m_Njet_a;
  vector<int> m_Njet_b;
  vector<int> m_Nbjet_a;
  vector<int> m_Nbjet_b;
  vector<int> m_Nlep_a;
  vector<int> m_Nlep_b;
  vector<int> m_Njet_ga;
  vector<int> m_Njet_gb;
  vector<int> m_Nbjet_ga;
  vector<int> m_Nbjet_gb;
  vector<int> m_Nlep_ga;
  vector<int> m_Nlep_gb;
  vector<vector<int> > m_index_jet_a;
  vector<vector<int> > m_index_jet_b;
  vector<vector<int> > m_index_lep_a;
  vector<vector<int> > m_index_lep_b;
  vector<vector<int> > m_index_jet_ga;
  vector<vector<int> > m_index_jet_gb;
  vector<vector<int> > m_index_lep_ga;
  vector<vector<int> > m_index_lep_gb;

  vector<double> m_MSS;
  vector<double> m_PSS;
  vector<double> m_cosSS;
  vector<double> m_dphiSS;
  vector<double> m_PTSS;
  vector<double> m_PzSS;

  vector<double> m_MCa;
  vector<double> m_cosCa;
  vector<double> m_MCb;
  vector<double> m_cosCb;

  vector<double> m_MGCa;
  vector<double> m_cosGCa;
  vector<double> m_MGCb;
  vector<double> m_cosGCb;

  vector<double> m_MVa;
  vector<double> m_PVa;
  vector<double> m_cosVa;
  vector<double> m_MVb;
  vector<double> m_PVb;
  vector<double> m_cosVb;

  vector<double> m_H11SS;
  vector<double> m_H21SS;
  vector<double> m_HT21SS;
  vector<double> m_H22SS;
  vector<double> m_HT22SS;
  vector<double> m_H42SS;
  vector<double> m_HT42SS;
  
  vector<double> m_H11Ca;
  vector<double> m_H11Cb;
  vector<double> m_H21Ca;
  vector<double> m_H21Cb;


  // ISR trees analysis
  vector<int> m_Njet_ISR;
  vector<int> m_Njet_S;
  vector<int> m_Nbjet_ISR;
  vector<int> m_Nbjet_S;
  vector<int> m_Nlep_ISR;
  vector<int> m_Nlep_S;
  vector<vector<int> > m_index_jet_ISR;
  vector<vector<int> > m_index_jet_S;
  vector<vector<int> > m_index_lep_ISR;
  vector<vector<int> > m_index_lep_S;
  vector<double> m_PTISR;
  vector<double> m_PTCM;
  vector<double> m_RISR;
  vector<double> m_cosCM;
  vector<double> m_cosS;
  vector<double> m_MISR;
  vector<double> m_MS;
  vector<double> m_MV;
  vector<double> m_ML;
  vector<double> m_dphiCMI;
  vector<double> m_dphiSI;
  vector<double> m_dphiISRI;

  // which tree are we using for PAIR?
  // vanilla - index 0
  bool m_Is_1L_2J;
  bool m_Is_2L_2J;
  bool m_Is_1L_1L;
  bool m_Is_2L_1L;
  bool m_Is_2L_2L;
  // b-aware - index 0
  bool m_Is_1L_B;
  bool m_Is_2L_B;
  bool m_Is_1LB_1LB;
  bool m_Is_3L_B;
  


  // RestFrames frames and friends

  // ISR-style trees
  LabRecoFrame*       LAB_ISR[3];
  DecayRecoFrame*     CM_ISR[3];
  DecayRecoFrame*     S_ISR[3];
  VisibleRecoFrame*   ISR_ISR[3];
  VisibleRecoFrame*   V_ISR[3];
  VisibleRecoFrame*   L_ISR[3];
  InvisibleRecoFrame* I_ISR[3];

  InvisibleGroup*      INV_ISR[3];
  SetMassInvJigsaw*    InvM_ISR[3];
  CombinatoricGroup*   COMB_ISR[3];
  MinMassesCombJigsaw* CombSplit_ISR[3];
  
  // Sparticle pair-production trees
  LabRecoFrame*            LAB_PAIR[5];
  DecayRecoFrame*          S_PAIR[5];
  DecayRecoFrame*          Ca_PAIR[5];
  DecayRecoFrame*          Cb_PAIR[5];
  DecayRecoFrame*          GCa_PAIR[5];
  DecayRecoFrame*          GCb_PAIR[5];
  SelfAssemblingRecoFrame* VSAa_PAIR[5];
  SelfAssemblingRecoFrame* VSAb_PAIR[5];
  VisibleRecoFrame*        Va_PAIR[5];
  VisibleRecoFrame*        Vb_PAIR[5];
  SelfAssemblingRecoFrame* GVSAa_PAIR[5];
  SelfAssemblingRecoFrame* GVSAb_PAIR[5];
  VisibleRecoFrame*        GVa_PAIR[5];
  VisibleRecoFrame*        GVb_PAIR[5];
  InvisibleRecoFrame*      Ia_PAIR[5];
  InvisibleRecoFrame*      Ib_PAIR[5];

  InvisibleGroup*       INV_PAIR[5];
  SetMassInvJigsaw*     InvM_PAIR[5];
  SetRapidityInvJigsaw* InvEta_PAIR[5];
  MinMassesSqInvJigsaw* InvSplit_PAIR[5];
  CombinatoricGroup*    COMB_PAIR[5];
  MinMassesCombJigsaw*  CombSplit_PAIR[5];
  MinMassesCombJigsaw*  CombSplita_PAIR[5];
  MinMassesCombJigsaw*  CombSplitb_PAIR[5];
  CombinatoricGroup*    COMBa_PAIR[5];
  CombinatoricGroup*    COMBb_PAIR[5];

};

#endif
