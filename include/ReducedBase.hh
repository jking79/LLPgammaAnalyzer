//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 11 12:53:08 2019 by ROOT version 6.14/04
// from TTree KUAnalysis/KUAnalysis
// found on file: TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8All.root
//////////////////////////////////////////////////////////

#ifndef ReducedBase_h
#define ReducedBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
//#include "vector"

using std::vector;

class ReducedBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        weight;
   Double_t        MET;
   Double_t        MET_phi;
   Double_t        genMET;
   Double_t        genMET_phi;
   Double_t        HT;
   Int_t           Nele;
   Int_t           Nmu;
   Int_t           Nlep;
   vector<double>  *PT_lep;
   vector<double>  *Eta_lep;
   vector<double>  *Phi_lep;
   vector<double>  *M_lep;
   vector<int>     *Charge_lep;
   vector<int>     *PDGID_lep;
   vector<int>     *RelIso_lep;
   vector<int>     *MiniIso_lep;
   vector<int>     *ID_lep;
   vector<int>     *Index_lep;
   Int_t           Njet;
   Int_t           Nbjet;
   vector<double>  *PT_jet;
   vector<double>  *Eta_jet;
   vector<double>  *Phi_jet;
   vector<double>  *M_jet;
   vector<double>  *Btag_jet;
   Int_t           genNele;
   Int_t           genNmu;
   Int_t           genNlep;
   vector<double>  *genPT_lep;
   vector<double>  *genEta_lep;
   vector<double>  *genPhi_lep;
   vector<double>  *genM_lep;
   vector<int>     *genCharge_lep;
   vector<int>     *genPDGID_lep;
   vector<int>     *genIndex_lep;
   Int_t           genNnu;
   vector<double>  *genPT_nu;
   vector<double>  *genEta_nu;
   vector<double>  *genPhi_nu;
   vector<int>     *genPDGID_nu;
   Int_t           genNboson;
   vector<double>  *genPT_boson;
   vector<double>  *genEta_boson;
   vector<double>  *genPhi_boson;
   vector<double>  *genM_boson;
   vector<int>     *genPDGID_boson;
   Int_t           genNsusy;
   vector<double>  *genPT_susy;
   vector<double>  *genEta_susy;
   vector<double>  *genPhi_susy;
   vector<double>  *genM_susy;
   vector<int>     *genPDGID_susy;
   vector<int>     *Njet_a;
   vector<int>     *Njet_b;
   vector<int>     *Nbjet_a;
   vector<int>     *Nbjet_b;
   vector<int>     *Nlep_a;
   vector<int>     *Nlep_b;
   vector<int>     *Njet_ga;
   vector<int>     *Njet_gb;
   vector<int>     *Nbjet_ga;
   vector<int>     *Nbjet_gb;
   vector<int>     *Nlep_ga;
   vector<int>     *Nlep_gb;
   vector<vector<int> > *index_jet_a;
   vector<vector<int> > *index_jet_b;
   vector<vector<int> > *index_lep_a;
   vector<vector<int> > *index_lep_b;
   vector<vector<int> > *index_jet_ga;
   vector<vector<int> > *index_jet_gb;
   vector<vector<int> > *index_lep_ga;
   vector<vector<int> > *index_lep_gb;
   vector<double>  *MSS;
   vector<double>  *PSS;
   vector<double>  *cosSS;
   vector<double>  *dphiSS;
   vector<double>  *PTSS;
   vector<double>  *PzSS;
   vector<double>  *MCa;
   vector<double>  *cosCa;
   vector<double>  *MCb;
   vector<double>  *cosCb;
   vector<double>  *MGCa;
   vector<double>  *cosGCa;
   vector<double>  *MGCb;
   vector<double>  *cosGCb;
   vector<double>  *H11SS;
   vector<double>  *H21SS;
   vector<double>  *HT21SS;
   vector<double>  *H22SS;
   vector<double>  *HT22SS;
   vector<double>  *H42SS;
   vector<double>  *HT42SS;
   vector<double>  *H11Ca;
   vector<double>  *H11Cb;
   vector<double>  *H21Ca;
   vector<double>  *H21Cb;
   vector<double>  *MVa;
   vector<double>  *PVa;
   vector<double>  *cosVa;
   vector<double>  *MVb;
   vector<double>  *PVb;
   vector<double>  *cosVb;
   Bool_t          Is_1L_2J;
   Bool_t          Is_2L_2J;
   Bool_t          Is_1L_1L;
   Bool_t          Is_2L_1L;
   Bool_t          Is_2L_2L;
   Bool_t          Is_1L_B;
   Bool_t          Is_2L_B;
   Bool_t          Is_1LB_1LB;
   Bool_t          Is_3L_B;
  
   vector<int>     *Njet_ISR;
   vector<int>     *Njet_S;
   vector<int>     *Nbjet_ISR;
   vector<int>     *Nbjet_S;
   vector<int>     *Nlep_ISR;
   vector<int>     *Nlep_S;
   vector<vector<int> > *index_jet_ISR;
   vector<vector<int> > *index_jet_S;
   vector<vector<int> > *index_lep_ISR;
   vector<vector<int> > *index_lep_S;
   vector<double>  *PTISR;
   vector<double>  *PTCM;
   vector<double>  *RISR;
   vector<double>  *cosCM;
   vector<double>  *cosS;
   vector<double>  *MISR;
   vector<double>  *MS;
   vector<double>  *MV;
   vector<double>  *ML;
   vector<double>  *dphiCMI;
   vector<double>  *dphiSI;
   vector<double>  *dphiISRI;

   // List of branches
   TBranch        *b_weight;   //!
   TBranch        *b_MET;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_genMET;   //!
   TBranch        *b_genMET_phi;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_Nele;   //!
   TBranch        *b_Nmu;   //!
   TBranch        *b_Nlep;   //!
   TBranch        *b_PT_lep;   //!
   TBranch        *b_Eta_lep;   //!
   TBranch        *b_Phi_lep;   //!
   TBranch        *b_M_lep;   //!
   TBranch        *b_Charge_lep;   //!
   TBranch        *b_PDGID_lep;   //!
   TBranch        *b_RelIso_lep;   //!
   TBranch        *b_MiniIso_lep;   //!
   TBranch        *b_ID_lep;   //!
   TBranch        *b_Index_lep;   //!
   TBranch        *b_Njet;   //!
   TBranch        *b_Nbjet;   //!
   TBranch        *b_PT_jet;   //!
   TBranch        *b_Eta_jet;   //!
   TBranch        *b_Phi_jet;   //!
   TBranch        *b_M_jet;   //!
   TBranch        *b_Btag_jet;   //!
   TBranch        *b_genNele;   //!
   TBranch        *b_genNmu;   //!
   TBranch        *b_genNlep;   //!
   TBranch        *b_genPT_lep;   //!
   TBranch        *b_genEta_lep;   //!
   TBranch        *b_genPhi_lep;   //!
   TBranch        *b_genM_lep;   //!
   TBranch        *b_genCharge_lep;   //!
   TBranch        *b_genPDGID_lep;   //!
   TBranch        *b_genIndex_lep;   //!
   TBranch        *b_genNnu;   //!
   TBranch        *b_genPT_nu;   //!
   TBranch        *b_genEta_nu;   //!
   TBranch        *b_genPhi_nu;   //!
   TBranch        *b_genPDGID_nu;   //!
   TBranch        *b_genNboson;   //!
   TBranch        *b_genPT_boson;   //!
   TBranch        *b_genEta_boson;   //!
   TBranch        *b_genPhi_boson;   //!
   TBranch        *b_genM_boson;   //!
   TBranch        *b_genPDGID_boson;   //!
   TBranch        *b_genNsusy;   //!
   TBranch        *b_genPT_susy;   //!
   TBranch        *b_genEta_susy;   //!
   TBranch        *b_genPhi_susy;   //!
   TBranch        *b_genM_susy;   //!
   TBranch        *b_genPDGID_susy;   //!
   TBranch        *b_Njet_a;   //!
   TBranch        *b_Njet_b;   //!
   TBranch        *b_Nbjet_a;   //!
   TBranch        *b_Nbjet_b;   //!
   TBranch        *b_Nlep_a;   //!
   TBranch        *b_Nlep_b;   //!
   TBranch        *b_Njet_ga;   //!
   TBranch        *b_Njet_gb;   //!
   TBranch        *b_Nbjet_ga;   //!
   TBranch        *b_Nbjet_gb;   //!
   TBranch        *b_Nlep_ga;   //!
   TBranch        *b_Nlep_gb;   //!
   TBranch        *b_index_jet_a;   //!
   TBranch        *b_index_jet_b;   //!
   TBranch        *b_index_lep_a;   //!
   TBranch        *b_index_lep_b;   //!
   TBranch        *b_index_jet_ga;   //!
   TBranch        *b_index_jet_gb;   //!
   TBranch        *b_index_lep_ga;   //!
   TBranch        *b_index_lep_gb;   //!
   TBranch        *b_MSS;   //!
   TBranch        *b_PSS;   //!
   TBranch        *b_cosSS;   //!
   TBranch        *b_dphiSS;   //!
   TBranch        *b_PTSS;   //!
   TBranch        *b_PzSS;   //!
   TBranch        *b_MCa;   //!
   TBranch        *b_cosCa;   //!
   TBranch        *b_MCb;   //!
   TBranch        *b_cosCb;   //!
   TBranch        *b_MGCa;   //!
   TBranch        *b_cosGCa;   //!
   TBranch        *b_MGCb;   //!
   TBranch        *b_cosGCb;   //!
   TBranch        *b_H11SS;   //!
   TBranch        *b_H21SS;   //!
   TBranch        *b_HT21SS;   //!
   TBranch        *b_H22SS;   //!
   TBranch        *b_HT22SS;   //!
   TBranch        *b_H42SS;   //!
   TBranch        *b_HT42SS;   //!
   TBranch        *b_H11Ca;   //!
   TBranch        *b_H11Cb;   //!
   TBranch        *b_H21Ca;   //!
   TBranch        *b_H21Cb;   //!
   TBranch        *b_MVa;   //!
   TBranch        *b_PVa;   //!
   TBranch        *b_cosVa;   //!
   TBranch        *b_MVb;   //!
   TBranch        *b_PVb;   //!
   TBranch        *b_cosVb;   //!
   TBranch        *b_Is_1L_2J;   //!
   TBranch        *b_Is_2L_2J;   //!
   TBranch        *b_Is_1L_1L;   //!
   TBranch        *b_Is_2L_1L;   //!
   TBranch        *b_Is_2L_2L;   //!
   TBranch        *b_Is_1L_B;   //!
   TBranch        *b_Is_2L_B;   //!
   TBranch        *b_Is_1LB_1LB;   //!
   TBranch        *b_Is_3L_B;   //!
   TBranch        *b_Njet_ISR;   //!
   TBranch        *b_Njet_S;   //!
   TBranch        *b_Nbjet_ISR;   //!
   TBranch        *b_Nbjet_S;   //!
   TBranch        *b_Nlep_ISR;   //!
   TBranch        *b_Nlep_S;   //!
   TBranch        *b_index_jet_ISR;   //!
   TBranch        *b_index_jet_S;   //!
   TBranch        *b_index_lep_ISR;   //!
   TBranch        *b_index_lep_S;   //!
   TBranch        *b_PTISR;   //!
   TBranch        *b_PTCM;   //!
   TBranch        *b_RISR;   //!
   TBranch        *b_cosCM;   //!
   TBranch        *b_cosS;   //!
   TBranch        *b_MISR;   //!
   TBranch        *b_MS;   //!
   TBranch        *b_MV;   //!
   TBranch        *b_ML;   //!
   TBranch        *b_dphiCMI;   //!
   TBranch        *b_dphiSI;   //!
   TBranch        *b_dphiISRI;   //!

   ReducedBase(TTree *tree=0);
   virtual ~ReducedBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif


inline ReducedBase::ReducedBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      //TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8All.root");
      //if (!f || !f->IsOpen()) {
      //   f = new TFile("TTJets_TuneCUETP8M2T4_13TeV-amcatnloFXFX-pythia8All.root");
      std::cout << "No Input tree Chain,  Crashing now!!!! " << std::endl;
      //}
      //f->GetObject("KUAnalysis",tree);
   }
   Init(tree);
}

inline ReducedBase::~ReducedBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t ReducedBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t ReducedBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

inline void ReducedBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   PT_lep = 0;
   Eta_lep = 0;
   Phi_lep = 0;
   M_lep = 0;
   Charge_lep = 0;
   PDGID_lep = 0;
   RelIso_lep = 0;
   MiniIso_lep = 0;
   ID_lep = 0;
   Index_lep = 0;
   PT_jet = 0;
   Eta_jet = 0;
   Phi_jet = 0;
   M_jet = 0;
   Btag_jet = 0;
   genPT_lep = 0;
   genEta_lep = 0;
   genPhi_lep = 0;
   genM_lep = 0;
   genCharge_lep = 0;
   genPDGID_lep = 0;
   genIndex_lep = 0;
   genPT_nu = 0;
   genEta_nu = 0;
   genPhi_nu = 0;
   genPDGID_nu = 0;
   genPT_boson = 0;
   genEta_boson = 0;
   genPhi_boson = 0;
   genM_boson = 0;
   genPDGID_boson = 0;
   genPT_susy = 0;
   genEta_susy = 0;
   genPhi_susy = 0;
   genM_susy = 0;
   genPDGID_susy = 0;
   Njet_a = 0;
   Njet_b = 0;
   Nbjet_a = 0;
   Nbjet_b = 0;
   Nlep_a = 0;
   Nlep_b = 0;
   Njet_ga = 0;
   Njet_gb = 0;
   Nbjet_ga = 0;
   Nbjet_gb = 0;
   Nlep_ga = 0;
   Nlep_gb = 0;
   index_jet_a = 0;
   index_jet_b = 0;
   index_lep_a = 0;
   index_lep_b = 0;
   index_jet_ga = 0;
   index_jet_gb = 0;
   index_lep_ga = 0;
   index_lep_gb = 0;
   MSS = 0;
   PSS = 0;
   cosSS = 0;
   dphiSS = 0;
   PTSS = 0;
   PzSS = 0;
   MCa = 0;
   cosCa = 0;
   MCb = 0;
   cosCb = 0;
   MGCa = 0;
   cosGCa = 0;
   MGCb = 0;
   cosGCb = 0;
   H11SS = 0;
   H21SS = 0;
   HT21SS = 0;
   H22SS = 0;
   HT22SS = 0;
   H42SS = 0;
   HT42SS = 0;
   H11Ca = 0;
   H11Cb = 0;
   H21Ca = 0;
   H21Cb = 0;
   MVa = 0;
   PVa = 0;
   cosVa = 0;
   MVb = 0;
   PVb = 0;
   cosVb = 0;
   Njet_ISR = 0;
   Njet_S = 0;
   Nbjet_ISR = 0;
   Nbjet_S = 0;
   Nlep_ISR = 0;
   Nlep_S = 0;
   index_jet_ISR = 0;
   index_jet_S = 0;
   index_lep_ISR = 0;
   index_lep_S = 0;
   PTISR = 0;
   PTCM = 0;
   RISR = 0;
   cosCM = 0;
   cosS = 0;
   MISR = 0;
   MS = 0;
   MV = 0;
   ML = 0;
   dphiCMI = 0;
   dphiSI = 0;
   dphiISRI = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("MET", &MET, &b_MET);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("genMET", &genMET, &b_genMET);
   fChain->SetBranchAddress("genMET_phi", &genMET_phi, &b_genMET_phi);
   fChain->SetBranchAddress("HT", &HT, &b_HT);
   fChain->SetBranchAddress("Nele", &Nele, &b_Nele);
   fChain->SetBranchAddress("Nmu", &Nmu, &b_Nmu);
   fChain->SetBranchAddress("Nlep", &Nlep, &b_Nlep);
   fChain->SetBranchAddress("PT_lep", &PT_lep, &b_PT_lep);
   fChain->SetBranchAddress("Eta_lep", &Eta_lep, &b_Eta_lep);
   fChain->SetBranchAddress("Phi_lep", &Phi_lep, &b_Phi_lep);
   fChain->SetBranchAddress("M_lep", &M_lep, &b_M_lep);
   fChain->SetBranchAddress("Charge_lep", &Charge_lep, &b_Charge_lep);
   fChain->SetBranchAddress("PDGID_lep", &PDGID_lep, &b_PDGID_lep);
   fChain->SetBranchAddress("RelIso_lep", &RelIso_lep, &b_RelIso_lep);
   fChain->SetBranchAddress("MiniIso_lep", &MiniIso_lep, &b_MiniIso_lep);
   fChain->SetBranchAddress("ID_lep", &ID_lep, &b_ID_lep);
   fChain->SetBranchAddress("Index_lep", &Index_lep, &b_Index_lep);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
   fChain->SetBranchAddress("Nbjet", &Nbjet, &b_Nbjet);
   fChain->SetBranchAddress("PT_jet", &PT_jet, &b_PT_jet);
   fChain->SetBranchAddress("Eta_jet", &Eta_jet, &b_Eta_jet);
   fChain->SetBranchAddress("Phi_jet", &Phi_jet, &b_Phi_jet);
   fChain->SetBranchAddress("M_jet", &M_jet, &b_M_jet);
   fChain->SetBranchAddress("Btag_jet", &Btag_jet, &b_Btag_jet);
   fChain->SetBranchAddress("genNele", &genNele, &b_genNele);
   fChain->SetBranchAddress("genNmu", &genNmu, &b_genNmu);
   fChain->SetBranchAddress("genNlep", &genNlep, &b_genNlep);
   fChain->SetBranchAddress("genPT_lep", &genPT_lep, &b_genPT_lep);
   fChain->SetBranchAddress("genEta_lep", &genEta_lep, &b_genEta_lep);
   fChain->SetBranchAddress("genPhi_lep", &genPhi_lep, &b_genPhi_lep);
   fChain->SetBranchAddress("genM_lep", &genM_lep, &b_genM_lep);
   fChain->SetBranchAddress("genCharge_lep", &genCharge_lep, &b_genCharge_lep);
   fChain->SetBranchAddress("genPDGID_lep", &genPDGID_lep, &b_genPDGID_lep);
   fChain->SetBranchAddress("genIndex_lep", &genIndex_lep, &b_genIndex_lep);
   fChain->SetBranchAddress("genNnu", &genNnu, &b_genNnu);
   fChain->SetBranchAddress("genPT_nu", &genPT_nu, &b_genPT_nu);
   fChain->SetBranchAddress("genEta_nu", &genEta_nu, &b_genEta_nu);
   fChain->SetBranchAddress("genPhi_nu", &genPhi_nu, &b_genPhi_nu);
   fChain->SetBranchAddress("genPDGID_nu", &genPDGID_nu, &b_genPDGID_nu);
   fChain->SetBranchAddress("genNboson", &genNboson, &b_genNboson);
   fChain->SetBranchAddress("genPT_boson", &genPT_boson, &b_genPT_boson);
   fChain->SetBranchAddress("genEta_boson", &genEta_boson, &b_genEta_boson);
   fChain->SetBranchAddress("genPhi_boson", &genPhi_boson, &b_genPhi_boson);
   fChain->SetBranchAddress("genM_boson", &genM_boson, &b_genM_boson);
   fChain->SetBranchAddress("genPDGID_boson", &genPDGID_boson, &b_genPDGID_boson);
   fChain->SetBranchAddress("genNsusy", &genNsusy, &b_genNsusy);
   fChain->SetBranchAddress("genPT_susy", &genPT_susy, &b_genPT_susy);
   fChain->SetBranchAddress("genEta_susy", &genEta_susy, &b_genEta_susy);
   fChain->SetBranchAddress("genPhi_susy", &genPhi_susy, &b_genPhi_susy);
   fChain->SetBranchAddress("genM_susy", &genM_susy, &b_genM_susy);
   fChain->SetBranchAddress("genPDGID_susy", &genPDGID_susy, &b_genPDGID_susy);
   fChain->SetBranchAddress("Njet_a", &Njet_a, &b_Njet_a);
   fChain->SetBranchAddress("Njet_b", &Njet_b, &b_Njet_b);
   fChain->SetBranchAddress("Nbjet_a", &Nbjet_a, &b_Nbjet_a);
   fChain->SetBranchAddress("Nbjet_b", &Nbjet_b, &b_Nbjet_b);
   fChain->SetBranchAddress("Nlep_a", &Nlep_a, &b_Nlep_a);
   fChain->SetBranchAddress("Nlep_b", &Nlep_b, &b_Nlep_b);
   fChain->SetBranchAddress("Njet_ga", &Njet_ga, &b_Njet_ga);
   fChain->SetBranchAddress("Njet_gb", &Njet_gb, &b_Njet_gb);
   fChain->SetBranchAddress("Nbjet_ga", &Nbjet_ga, &b_Nbjet_ga);
   fChain->SetBranchAddress("Nbjet_gb", &Nbjet_gb, &b_Nbjet_gb);
   fChain->SetBranchAddress("Nlep_ga", &Nlep_ga, &b_Nlep_ga);
   fChain->SetBranchAddress("Nlep_gb", &Nlep_gb, &b_Nlep_gb);
   fChain->SetBranchAddress("index_jet_a", &index_jet_a, &b_index_jet_a);
   fChain->SetBranchAddress("index_jet_b", &index_jet_b, &b_index_jet_b);
   fChain->SetBranchAddress("index_lep_a", &index_lep_a, &b_index_lep_a);
   fChain->SetBranchAddress("index_lep_b", &index_lep_b, &b_index_lep_b);
   fChain->SetBranchAddress("index_jet_ga", &index_jet_ga, &b_index_jet_ga);
   fChain->SetBranchAddress("index_jet_gb", &index_jet_gb, &b_index_jet_gb);
   fChain->SetBranchAddress("index_lep_ga", &index_lep_ga, &b_index_lep_ga);
   fChain->SetBranchAddress("index_lep_gb", &index_lep_gb, &b_index_lep_gb);
   fChain->SetBranchAddress("MSS", &MSS, &b_MSS);
   fChain->SetBranchAddress("PSS", &PSS, &b_PSS);
   fChain->SetBranchAddress("cosSS", &cosSS, &b_cosSS);
   fChain->SetBranchAddress("dphiSS", &dphiSS, &b_dphiSS);
   fChain->SetBranchAddress("PTSS", &PTSS, &b_PTSS);
   fChain->SetBranchAddress("PzSS", &PzSS, &b_PzSS);
   fChain->SetBranchAddress("MCa", &MCa, &b_MCa);
   fChain->SetBranchAddress("cosCa", &cosCa, &b_cosCa);
   fChain->SetBranchAddress("MCb", &MCb, &b_MCb);
   fChain->SetBranchAddress("cosCb", &cosCb, &b_cosCb);
   fChain->SetBranchAddress("MGCa", &MGCa, &b_MGCa);
   fChain->SetBranchAddress("cosGCa", &cosGCa, &b_cosGCa);
   fChain->SetBranchAddress("MGCb", &MGCb, &b_MGCb);
   fChain->SetBranchAddress("cosGCb", &cosGCb, &b_cosGCb);
   fChain->SetBranchAddress("H11SS", &H11SS, &b_H11SS);
   fChain->SetBranchAddress("H21SS", &H21SS, &b_H21SS);
   fChain->SetBranchAddress("HT21SS", &HT21SS, &b_HT21SS);
   fChain->SetBranchAddress("H22SS", &H22SS, &b_H22SS);
   fChain->SetBranchAddress("HT22SS", &HT22SS, &b_HT22SS);
   fChain->SetBranchAddress("H42SS", &H42SS, &b_H42SS);
   fChain->SetBranchAddress("HT42SS", &HT42SS, &b_HT42SS);
   fChain->SetBranchAddress("H11Ca", &H11Ca, &b_H11Ca);
   fChain->SetBranchAddress("H11Cb", &H11Cb, &b_H11Cb);
   fChain->SetBranchAddress("H21Ca", &H21Ca, &b_H21Ca);
   fChain->SetBranchAddress("H21Cb", &H21Cb, &b_H21Cb);
   fChain->SetBranchAddress("MVa", &MVa, &b_MVa);
   fChain->SetBranchAddress("PVa", &PVa, &b_PVa);
   fChain->SetBranchAddress("cosVa", &cosVa, &b_cosVa);
   fChain->SetBranchAddress("MVb", &MVb, &b_MVb);
   fChain->SetBranchAddress("PVb", &PVb, &b_PVb);
   fChain->SetBranchAddress("cosVb", &cosVb, &b_cosVb);
   fChain->SetBranchAddress("Is_1L_2J", &Is_1L_2J, &b_Is_1L_2J);
   fChain->SetBranchAddress("Is_2L_2J", &Is_2L_2J, &b_Is_2L_2J);
   fChain->SetBranchAddress("Is_1L_1L", &Is_1L_1L, &b_Is_1L_1L);
   fChain->SetBranchAddress("Is_2L_1L", &Is_2L_1L, &b_Is_2L_1L);
   fChain->SetBranchAddress("Is_2L_2L", &Is_2L_2L, &b_Is_2L_2L);
   fChain->SetBranchAddress("Is_1L_B", &Is_1L_B, &b_Is_1L_B);
   fChain->SetBranchAddress("Is_2L_B", &Is_2L_B, &b_Is_2L_B);
   fChain->SetBranchAddress("Is_1LB_1LB", &Is_1LB_1LB, &b_Is_1LB_1LB);
   fChain->SetBranchAddress("Is_3L_B", &Is_3L_B, &b_Is_3L_B);
   fChain->SetBranchAddress("Njet_ISR", &Njet_ISR, &b_Njet_ISR);
   fChain->SetBranchAddress("Njet_S", &Njet_S, &b_Njet_S);
   fChain->SetBranchAddress("Nbjet_ISR", &Nbjet_ISR, &b_Nbjet_ISR);
   fChain->SetBranchAddress("Nbjet_S", &Nbjet_S, &b_Nbjet_S);
   fChain->SetBranchAddress("Nlep_ISR", &Nlep_ISR, &b_Nlep_ISR);
   fChain->SetBranchAddress("Nlep_S", &Nlep_S, &b_Nlep_S);
   fChain->SetBranchAddress("index_jet_ISR", &index_jet_ISR, &b_index_jet_ISR);
   fChain->SetBranchAddress("index_jet_S", &index_jet_S, &b_index_jet_S);
   fChain->SetBranchAddress("index_lep_ISR", &index_lep_ISR, &b_index_lep_ISR);
   fChain->SetBranchAddress("index_lep_S", &index_lep_S, &b_index_lep_S);
   fChain->SetBranchAddress("PTISR", &PTISR, &b_PTISR);
   fChain->SetBranchAddress("PTCM", &PTCM, &b_PTCM);
   fChain->SetBranchAddress("RISR", &RISR, &b_RISR);
   fChain->SetBranchAddress("cosCM", &cosCM, &b_cosCM);
   fChain->SetBranchAddress("cosS", &cosS, &b_cosS);
   fChain->SetBranchAddress("MISR", &MISR, &b_MISR);
   fChain->SetBranchAddress("MS", &MS, &b_MS);
   fChain->SetBranchAddress("MV", &MV, &b_MV);
   fChain->SetBranchAddress("ML", &ML, &b_ML);
   fChain->SetBranchAddress("dphiCMI", &dphiCMI, &b_dphiCMI);
   fChain->SetBranchAddress("dphiSI", &dphiSI, &b_dphiSI);
   fChain->SetBranchAddress("dphiISRI", &dphiISRI, &b_dphiISRI);
   Notify();

   // Turn off/on different branches to improve processing speed
   fChain->SetBranchStatus("*",0);
   fChain->SetBranchStatus("weight", 1);
   fChain->SetBranchStatus("MET", 1);
   fChain->SetBranchStatus("N*_ISR", 1);
   fChain->SetBranchStatus("N*_S", 1);
   fChain->SetBranchStatus("PTISR", 1);
   fChain->SetBranchStatus("RISR", 1);
   fChain->SetBranchStatus("MS", 1);
   fChain->SetBranchStatus("MV", 1);
   fChain->SetBranchStatus("dphiISRI", 1);
   fChain->SetBranchStatus("dphiCMI", 1);
   fChain->SetBranchStatus("PTCM", 1);
   fChain->SetBranchStatus("ID_lep",1);
   fChain->SetBranchStatus("PDGID_lep",1);
}

inline Bool_t ReducedBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void ReducedBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t ReducedBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

