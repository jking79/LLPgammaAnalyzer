//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct  6 12:02:13 2022 by ROOT version 6.14/09
// from TTree llpgtree/llpgtree
// found on file: output_21.root
//////////////////////////////////////////////////////////

//#ifndef llpgana_ntuple_base_h
//#define llpgana_ntuple_base_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class llpgana_ntuple_base {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   //Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Int_t           nVtx;
   Float_t         vtxX;
   Float_t         vtxY;
   Float_t         vtxZ;
   Float_t         jetHt;
   UInt_t          nJets;
   UInt_t          nGoodDrJets;
   UInt_t          nGoodScJets;
   UInt_t          nGoodBcJets;
   UInt_t          nUnJets;
   vector<float>   *jetE;
   vector<float>   *jetPt;
   vector<float>   *jetEta;
   vector<float>   *jetPhi;
   vector<float>   *jetNHF;
   vector<float>   *jetNEMF;
   vector<float>   *jetCHF;
   vector<float>   *jetCEMF;
   vector<float>   *jetMUF;
   vector<float>   *jetNHM;
   vector<float>   *jetCHM;
   vector<int>     *jetPHM;
   vector<int>     *jetELM;
   vector<float>   *jetPHE;
   vector<float>   *jetPHEF;
   vector<float>   *jetELE;
   vector<float>   *jetELEF;
   vector<float>   *jetMUE;
   vector<int>     *jetID;
   vector<unsigned int> *cljBcCnt;
   UInt_t          nCaloJets;
   vector<float>   *cljSeedTOFTime;
   vector<float>   *cljCMeanTime;
   vector<float>   *cljBc3dEx;
   vector<float>   *cljBc3dEy;
   vector<float>   *cljBc3dEz;
   vector<float>   *cljBc3dEv;
   vector<float>   *cljBc3dEslope;
   vector<float>   *cljBc3dEchisp;
   vector<float>   *cljBc2dEx;
   vector<float>   *cljBc2dEy;
   vector<float>   *cljBc2dEv;
   vector<float>   *cljBc2dEslope;
   vector<float>   *cljBc2dEchisp;
   vector<float>   *cljCDrMeanTime;
   vector<double>  *cljPt;
   vector<double>  *cljEnergy;
   vector<double>  *cljPhi;
   vector<double>  *cljEta;
   vector<double>  *cljPx;
   vector<double>  *cljPy;
   vector<double>  *cljPz;
   UInt_t          nPhotons;
   vector<float>   *phoSeedTOFTime;
   vector<float>   *phoCMeanTime;
   vector<float>   *phoSc3dEx;
   vector<float>   *phoSc3dEy;
   vector<float>   *phoSc3dEz;
   vector<float>   *phoSc3dEv;
   vector<float>   *phoSc3dEslope;
   vector<float>   *phoSc3dEchisp;
   vector<float>   *phoSc2dEx;
   vector<float>   *phoSc2dEy;
   vector<float>   *phoSc2dEv;
   vector<float>   *phoSc2dEslope;
   vector<float>   *phoSc2dEchisp;
   vector<double>  *phoPt;
   vector<double>  *phoEnergy;
   vector<double>  *phoPhi;
   vector<double>  *phoEta;
   vector<double>  *phoPx;
   vector<double>  *phoPy;
   vector<double>  *phoPz;
   UInt_t          nOotPhotons;
   vector<float>   *ootPhoSeedTOFTime;
   vector<float>   *ootPhoCMeanTime;
   vector<float>   *ootPhoSc3dEx;
   vector<float>   *ootPhoSc3dEy;
   vector<float>   *ootPhoSc3dEz;
   vector<float>   *ootPhoSc3dEv;
   vector<float>   *ootPhoSc3dEslope;
   vector<float>   *ootPhoSc3dEchisp;
   vector<float>   *ootPhoSc2dEx;
   vector<float>   *ootPhoSc2dEy;
   vector<float>   *ootPhoSc2dEv;
   vector<float>   *ootPhoSc2dEslope;
   vector<float>   *ootPhoSc2dEchisp;
   vector<double>  *ootPhoPt;
   vector<double>  *ootPhoEnergy;
   vector<double>  *ootPhoPhi;
   vector<double>  *ootPhoEta;
   vector<double>  *ootPhoPx;
   vector<double>  *ootPhoPy;
   vector<double>  *ootPhoPz;
   UInt_t          nElectrons;
   vector<float>   *eleSeedTOFTime;
   vector<float>   *eleCMeanTime;
   vector<float>   *eleSc3dEx;
   vector<float>   *eleSc3dEy;
   vector<float>   *eleSc3dEz;
   vector<float>   *eleSc3dEv;
   vector<float>   *eleSc3dEslope;
   vector<float>   *eleSc3dEchisp;
   vector<float>   *eleSc2dEx;
   vector<float>   *eleSc2dEy;
   vector<float>   *eleSc2dEv;
   vector<float>   *eleSc2dEslope;
   vector<float>   *eleSc2dEchisp;
   vector<double>  *elePt;
   vector<double>  *eleEnergy;
   vector<double>  *elePhi;
   vector<double>  *eleEta;
   vector<double>  *elePx;
   vector<double>  *elePy;
   vector<double>  *elePz;
   vector<float>   *jetSumEPFrac;
   vector<float>   *jetEPEnergy;
   vector<float>   *jetEMEnergy;
   vector<float>   *jetEMEnrFrac;
   vector<float>   *jetEPEnrFrac;
   vector<float>   *jetDrLeadEta;
   vector<float>   *jetDrLeadPhi;
   vector<float>   *jetDrLeadEnr;
   vector<float>   *sJetDrRHEnergy;
   vector<float>   *jetDrEMF;
   vector<unsigned int> *jetDrRhCnt;
   vector<unsigned int> *jetBcTimesCnt;
   vector<float>   *jetBcSumRHEnr;
   vector<float>   *jetBcEMFr;
   vector<unsigned int> *jetBcRhCnt;
   vector<unsigned int> *jetBcGrpCnt;
   vector<unsigned int> *nJetScMatch;
   vector<unsigned int> *jetScRhCnt;
   vector<float>   *sJetScEnergy;
   vector<float>   *sJetScPhEnergy;
   vector<float>   *sJetScRhEnergy;
   vector<float>   *jetScEMF;
   vector<float>   *jetDRMuTime;
   vector<float>   *jetDRTimeError;
   vector<float>   *jetDRTimeRMS;
   vector<float>   *jetDRMedTime;
   vector<float>   *jetCDRMuTime;
   vector<float>   *jetCDRMedTime;
   vector<float>   *jetSCMuTime;
   vector<float>   *jetSCMedTime;
   vector<float>   *jetCSCMuTime;
   vector<float>   *jetCSCMedTime;
   vector<float>   *jetCBCMuTime;
   vector<float>   *jetCBCMedTime;
   vector<float>   *jetGenImpactAngle;
   vector<float>   *jetGenTime;
   vector<float>   *jetGenPt;
   vector<float>   *jetGenEta;
   vector<float>   *jetGenEnergy;
   vector<float>   *jetGenEMFrac;
   vector<float>   *jetGenDrMatch;
   vector<float>   *jetGenTimeVar;
   vector<float>   *jetGenTimeLLP;
   vector<float>   *jetGenLLPPurity;
   vector<float>   *jetGenNextBX;
   vector<float>   *jetGenNKids;
   vector<float>   *jetGenTOF;
   vector<float>   *jetImpactAngle;
   vector<float>   *jetSc3dEx;
   vector<float>   *jetSc3dEy;
   vector<float>   *jetSc3dEz;
   vector<float>   *jetSc3dEv;
   vector<float>   *jetSc3dEslope;
   vector<float>   *jetSc3dEchisp;
   vector<float>   *jetSc2dEx;
   vector<float>   *jetSc2dEy;
   vector<float>   *jetSc2dEv;
   vector<float>   *jetSc2dEslope;
   vector<float>   *jetSc2dEchisp;
   vector<float>   *jetSc2dEslope2;
   vector<float>   *jetSc2dEchisp2;
   vector<float>   *jetSc2dErangle;
   vector<float>   *jetSc2dEnxsum;
   Int_t           nRecHits;
   vector<float>   *rhPosX;
   vector<float>   *rhPosY;
   vector<float>   *rhPosZ;
   vector<float>   *rhPosEta;
   vector<float>   *rhPosPhi;
   vector<float>   *rhEnergy;
   vector<float>   *rhTime;
   vector<float>   *rhTimeErr;
   vector<float>   *rhTOF;
   vector<unsigned int> *rhID;
   vector<int>     *rhXtalI1;
   vector<int>     *rhXtalI2;
   vector<int>     *rhSubdet;
   vector<bool>    *rhisOOT;
   vector<bool>    *rhisGS6;
   vector<bool>    *rhisGS1;
   vector<float>   *rhadcToGeV;
   vector<float>   *rhped12;
   vector<float>   *rhped6;
   vector<float>   *rhped1;
   vector<float>   *rhpedrms12;
   vector<float>   *rhpedrms6;
   vector<float>   *rhpedrms1;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_jetHt;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_nGoodDrJets;   //!
   TBranch        *b_nGoodScJets;   //!
   TBranch        *b_nGoodBcJets;   //!
   TBranch        *b_nUnJets;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetPt;   //!
   TBranch        *b_jetEta;   //!
   TBranch        *b_jetPhi;   //!
   TBranch        *b_jetNHF;   //!
   TBranch        *b_jetNEMF;   //!
   TBranch        *b_jetCHF;   //!
   TBranch        *b_jetCEMF;   //!
   TBranch        *b_jetMUF;   //!
   TBranch        *b_jetNHM;   //!
   TBranch        *b_jetCHM;   //!
   TBranch        *b_jetPHM;   //!
   TBranch        *b_jetELM;   //!
   TBranch        *b_jetPHE;   //!
   TBranch        *b_jetPHEF;   //!
   TBranch        *b_jetELE;   //!
   TBranch        *b_jetELEF;   //!
   TBranch        *b_jetMUE;   //!
   TBranch        *b_jetID;   //!
   TBranch        *b_cljBcCnt;   //!
   TBranch        *b_nCaloJets;   //!
   TBranch        *b_cljSeedTOFTime;   //!
   TBranch        *b_cljCMeanTime;   //!
   TBranch        *b_cljBc3dEx;   //!
   TBranch        *b_cljBc3dEy;   //!
   TBranch        *b_cljBc3dEz;   //!
   TBranch        *b_cljBc3dEv;   //!
   TBranch        *b_cljBc3dEslope;   //!
   TBranch        *b_cljBc3dEchisp;   //!
   TBranch        *b_cljBc2dEx;   //!
   TBranch        *b_cljBc2dEy;   //!
   TBranch        *b_cljBc2dEv;   //!
   TBranch        *b_cljBc2dEslope;   //!
   TBranch        *b_cljBc2dEchisp;   //!
   TBranch        *b_cljCDrMeanTime;   //!
   TBranch        *b_cljPt;   //!
   TBranch        *b_cljEnergy;   //!
   TBranch        *b_cljPhi;   //!
   TBranch        *b_cljEta;   //!
   TBranch        *b_cljPx;   //!
   TBranch        *b_cljPy;   //!
   TBranch        *b_cljPz;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_phoSeedTOFTime;   //!
   TBranch        *b_phoCMeanTime;   //!
   TBranch        *b_phoSc3dEx;   //!
   TBranch        *b_phoSc3dEy;   //!
   TBranch        *b_phoSc3dEz;   //!
   TBranch        *b_phoSc3dEv;   //!
   TBranch        *b_phoSc3dEslope;   //!
   TBranch        *b_phoSc3dEchisp;   //!
   TBranch        *b_phoSc2dEx;   //!
   TBranch        *b_phoSc2dEy;   //!
   TBranch        *b_phoSc2dEv;   //!
   TBranch        *b_phoSc2dEslope;   //!
   TBranch        *b_phoSc2dEchisp;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEnergy;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_nOotPhotons;   //!
   TBranch        *b_ootPhoSeedTOFTime;   //!
   TBranch        *b_ootPhoCMeanTime;   //!
   TBranch        *b_ootPhoSc3dEx;   //!
   TBranch        *b_ootPhoSc3dEy;   //!
   TBranch        *b_ootPhoSc3dEz;   //!
   TBranch        *b_ootPhoSc3dEv;   //!
   TBranch        *b_ootPhoSc3dEslope;   //!
   TBranch        *b_ootPhoSc3dEchisp;   //!
   TBranch        *b_ootPhoSc2dEx;   //!
   TBranch        *b_ootPhoSc2dEy;   //!
   TBranch        *b_ootPhoSc2dEv;   //!
   TBranch        *b_ootPhoSc2dEslope;   //!
   TBranch        *b_ootPhoSc2dEchisp;   //!
   TBranch        *b_ootPhoPt;   //!
   TBranch        *b_ootPhoEnergy;   //!
   TBranch        *b_ootPhoPhi;   //!
   TBranch        *b_ootPhoEta;   //!
   TBranch        *b_ootPhoPx;   //!
   TBranch        *b_ootPhoPy;   //!
   TBranch        *b_ootPhoPz;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_eleSeedTOFTime;   //!
   TBranch        *b_eleCMeanTime;   //!
   TBranch        *b_eleSc3dEx;   //!
   TBranch        *b_eleSc3dEy;   //!
   TBranch        *b_eleSc3dEz;   //!
   TBranch        *b_eleSc3dEv;   //!
   TBranch        *b_eleSc3dEslope;   //!
   TBranch        *b_eleSc3dEchisp;   //!
   TBranch        *b_eleSc2dEx;   //!
   TBranch        *b_eleSc2dEy;   //!
   TBranch        *b_eleSc2dEv;   //!
   TBranch        *b_eleSc2dEslope;   //!
   TBranch        *b_eleSc2dEchisp;   //!
   TBranch        *b_elePt;   //!
   TBranch        *b_eleEnergy;   //!
   TBranch        *b_elePhi;   //!
   TBranch        *b_eleEta;   //!
   TBranch        *b_elePx;   //!
   TBranch        *b_elePy;   //!
   TBranch        *b_elePz;   //!
   TBranch        *b_jetSumEPFrac;   //!
   TBranch        *b_jetEPEnergy;   //!
   TBranch        *b_jetEMEnergy;   //!
   TBranch        *b_jetEMEnrFrac;   //!
   TBranch        *b_jetEPEnrFrac;   //!
   TBranch        *b_jetDrLeadEta;   //!
   TBranch        *b_jetDrLeadPhi;   //!
   TBranch        *b_jetDrLeadEnr;   //!
   TBranch        *b_sJetDrRHEnergy;   //!
   TBranch        *b_jetDrEMF;   //!
   TBranch        *b_jetDrRhCnt;   //!
   TBranch        *b_jetBcTimesCnt;   //!
   TBranch        *b_jetBcSumRHEnr;   //!
   TBranch        *b_jetBcEMFr;   //!
   TBranch        *b_jetBcRhCnt;   //!
   TBranch        *b_jetBcGrpCnt;   //!
   TBranch        *b_nJetScMatch;   //!
   TBranch        *b_jetScRhCnt;   //!
   TBranch        *b_sJetScEnergy;   //!
   TBranch        *b_sJetScPhEnergy;   //!
   TBranch        *b_sJetScRhEnergy;   //!
   TBranch        *b_jetScEMF;   //!
   TBranch        *b_jetDRMuTime;   //!
   TBranch        *b_jetDRTimeError;   //!
   TBranch        *b_jetDRTimeRMS;   //!
   TBranch        *b_jetDRMedTime;   //!
   TBranch        *b_jetCDRMuTime;   //!
   TBranch        *b_jetCDRMedTime;   //!
   TBranch        *b_jetSCMuTime;   //!
   TBranch        *b_jetSCMedTime;   //!
   TBranch        *b_jetCSCMuTime;   //!
   TBranch        *b_jetCSCMedTime;   //!
   TBranch        *b_jetCBCMuTime;   //!
   TBranch        *b_jetCBCMedTime;   //!
   TBranch        *b_jetGenImpactAngle;   //!
   TBranch        *b_jetGenTime;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenEnergy;   //!
   TBranch        *b_jetGenEMFrac;   //!
   TBranch        *b_jetGenDrMatch;   //!
   TBranch        *b_jetGenTimeVar;   //!
   TBranch        *b_jetGenTimeLLP;   //!
   TBranch        *b_jetGenLLPPurity;   //!
   TBranch        *b_jetGenNextBX;   //!
   TBranch        *b_jetGenNKids;   //!
   TBranch        *b_jetGenTOF;   //!
   TBranch        *b_jetImpactAngle;   //!
   TBranch        *b_jetSc3dEx;   //!
   TBranch        *b_jetSc3dEy;   //!
   TBranch        *b_jetSc3dEz;   //!
   TBranch        *b_jetSc3dEv;   //!
   TBranch        *b_jetSc3dEslope;   //!
   TBranch        *b_jetSc3dEchisp;   //!
   TBranch        *b_jetSc2dEx;   //!
   TBranch        *b_jetSc2dEy;   //!
   TBranch        *b_jetSc2dEv;   //!
   TBranch        *b_jetSc2dEslope;   //!
   TBranch        *b_jetSc2dEchisp;   //!
   TBranch        *b_jetSc2dEslope2;   //!
   TBranch        *b_jetSc2dEchisp2;   //!
   TBranch        *b_jetSc2dErangle;   //!
   TBranch        *b_jetSc2dEnxsum;   //!
   TBranch        *b_nRecHits;   //!
   TBranch        *b_rhPosX;   //!
   TBranch        *b_rhPosY;   //!
   TBranch        *b_rhPosZ;   //!
   TBranch        *b_rhPosEta;   //!
   TBranch        *b_rhPosPhi;   //!
   TBranch        *b_rhEnergy;   //!
   TBranch        *b_rhTime;   //!
   TBranch        *b_rhTimeErr;   //!
   TBranch        *b_rhTOF;   //!
   TBranch        *b_rhID;   //!
   TBranch        *b_rhXtalI1;   //!
   TBranch        *b_rhXtalI2;   //!
   TBranch        *b_rhSubdet;   //!
   TBranch        *b_rhisOOT;   //!
   TBranch        *b_rhisGS6;   //!
   TBranch        *b_rhisGS1;   //!
   TBranch        *b_rhadcToGeV;   //!
   TBranch        *b_rhped12;   //!
   TBranch        *b_rhped6;   //!
   TBranch        *b_rhped1;   //!
   TBranch        *b_rhpedrms12;   //!
   TBranch        *b_rhpedrms6;   //!
   TBranch        *b_rhpedrms1;   //!

   //llpgana_ntuple_base(TTree *tree=0);
   //virtual ~llpgana_ntuple_base();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef llpgana_ntuple_base_cxx
/*
llpgana_ntuple_base::llpgana_ntuple_base(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_21.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_21.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("output_21.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

llpgana_ntuple_base::~llpgana_ntuple_base()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t llpgana_ntuple_base::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t llpgana_ntuple_base::LoadTree(Long64_t entry)
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
*/

void llpgana_ntuple_base::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jetE = 0;
   jetPt = 0;
   jetEta = 0;
   jetPhi = 0;
   jetNHF = 0;
   jetNEMF = 0;
   jetCHF = 0;
   jetCEMF = 0;
   jetMUF = 0;
   jetNHM = 0;
   jetCHM = 0;
   jetPHM = 0;
   jetELM = 0;
   jetPHE = 0;
   jetPHEF = 0;
   jetELE = 0;
   jetELEF = 0;
   jetMUE = 0;
   jetID = 0;
   cljBcCnt = 0;
   cljSeedTOFTime = 0;
   cljCMeanTime = 0;
   cljBc3dEx = 0;
   cljBc3dEy = 0;
   cljBc3dEz = 0;
   cljBc3dEv = 0;
   cljBc3dEslope = 0;
   cljBc3dEchisp = 0;
   cljBc2dEx = 0;
   cljBc2dEy = 0;
   cljBc2dEv = 0;
   cljBc2dEslope = 0;
   cljBc2dEchisp = 0;
   cljCDrMeanTime = 0;
   cljPt = 0;
   cljEnergy = 0;
   cljPhi = 0;
   cljEta = 0;
   cljPx = 0;
   cljPy = 0;
   cljPz = 0;
   phoSeedTOFTime = 0;
   phoCMeanTime = 0;
   phoSc3dEx = 0;
   phoSc3dEy = 0;
   phoSc3dEz = 0;
   phoSc3dEv = 0;
   phoSc3dEslope = 0;
   phoSc3dEchisp = 0;
   phoSc2dEx = 0;
   phoSc2dEy = 0;
   phoSc2dEv = 0;
   phoSc2dEslope = 0;
   phoSc2dEchisp = 0;
   phoPt = 0;
   phoEnergy = 0;
   phoPhi = 0;
   phoEta = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   ootPhoSeedTOFTime = 0;
   ootPhoCMeanTime = 0;
   ootPhoSc3dEx = 0;
   ootPhoSc3dEy = 0;
   ootPhoSc3dEz = 0;
   ootPhoSc3dEv = 0;
   ootPhoSc3dEslope = 0;
   ootPhoSc3dEchisp = 0;
   ootPhoSc2dEx = 0;
   ootPhoSc2dEy = 0;
   ootPhoSc2dEv = 0;
   ootPhoSc2dEslope = 0;
   ootPhoSc2dEchisp = 0;
   ootPhoPt = 0;
   ootPhoEnergy = 0;
   ootPhoPhi = 0;
   ootPhoEta = 0;
   ootPhoPx = 0;
   ootPhoPy = 0;
   ootPhoPz = 0;
   eleSeedTOFTime = 0;
   eleCMeanTime = 0;
   eleSc3dEx = 0;
   eleSc3dEy = 0;
   eleSc3dEz = 0;
   eleSc3dEv = 0;
   eleSc3dEslope = 0;
   eleSc3dEchisp = 0;
   eleSc2dEx = 0;
   eleSc2dEy = 0;
   eleSc2dEv = 0;
   eleSc2dEslope = 0;
   eleSc2dEchisp = 0;
   elePt = 0;
   eleEnergy = 0;
   elePhi = 0;
   eleEta = 0;
   elePx = 0;
   elePy = 0;
   elePz = 0;
   jetSumEPFrac = 0;
   jetEPEnergy = 0;
   jetEMEnergy = 0;
   jetEMEnrFrac = 0;
   jetEPEnrFrac = 0;
   jetDrLeadEta = 0;
   jetDrLeadPhi = 0;
   jetDrLeadEnr = 0;
   sJetDrRHEnergy = 0;
   jetDrEMF = 0;
   jetDrRhCnt = 0;
   jetBcTimesCnt = 0;
   jetBcSumRHEnr = 0;
   jetBcEMFr = 0;
   jetBcRhCnt = 0;
   jetBcGrpCnt = 0;
   nJetScMatch = 0;
   jetScRhCnt = 0;
   sJetScEnergy = 0;
   sJetScPhEnergy = 0;
   sJetScRhEnergy = 0;
   jetScEMF = 0;
   jetDRMuTime = 0;
   jetDRTimeError = 0;
   jetDRTimeRMS = 0;
   jetDRMedTime = 0;
   jetCDRMuTime = 0;
   jetCDRMedTime = 0;
   jetSCMuTime = 0;
   jetSCMedTime = 0;
   jetCSCMuTime = 0;
   jetCSCMedTime = 0;
   jetCBCMuTime = 0;
   jetCBCMedTime = 0;
   jetGenImpactAngle = 0;
   jetGenTime = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenEnergy = 0;
   jetGenEMFrac = 0;
   jetGenDrMatch = 0;
   jetGenTimeVar = 0;
   jetGenTimeLLP = 0;
   jetGenLLPPurity = 0;
   jetGenNextBX = 0;
   jetGenNKids = 0;
   jetGenTOF = 0;
   jetImpactAngle = 0;
   jetSc3dEx = 0;
   jetSc3dEy = 0;
   jetSc3dEz = 0;
   jetSc3dEv = 0;
   jetSc3dEslope = 0;
   jetSc3dEchisp = 0;
   jetSc2dEx = 0;
   jetSc2dEy = 0;
   jetSc2dEv = 0;
   jetSc2dEslope = 0;
   jetSc2dEchisp = 0;
   jetSc2dEslope2 = 0;
   jetSc2dEchisp2 = 0;
   jetSc2dErangle = 0;
   jetSc2dEnxsum = 0;
   rhPosX = 0;
   rhPosY = 0;
   rhPosZ = 0;
   rhPosEta = 0;
   rhPosPhi = 0;
   rhEnergy = 0;
   rhTime = 0;
   rhTimeErr = 0;
   rhTOF = 0;
   rhID = 0;
   rhXtalI1 = 0;
   rhXtalI2 = 0;
   rhSubdet = 0;
   rhisOOT = 0;
   rhisGS6 = 0;
   rhisGS1 = 0;
   rhadcToGeV = 0;
   rhped12 = 0;
   rhped6 = 0;
   rhped1 = 0;
   rhpedrms12 = 0;
   rhpedrms6 = 0;
   rhpedrms1 = 0;

   // Set branch addresses and branch pointers
   //if (!tree) return;
   fChain = tree;
   //fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("jetHt", &jetHt, &b_jetHt);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("nGoodDrJets", &nGoodDrJets, &b_nGoodDrJets);
   fChain->SetBranchAddress("nGoodScJets", &nGoodScJets, &b_nGoodScJets);
   fChain->SetBranchAddress("nGoodBcJets", &nGoodBcJets, &b_nGoodBcJets);
   fChain->SetBranchAddress("nUnJets", &nUnJets, &b_nUnJets);
   fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
   fChain->SetBranchAddress("jetPt", &jetPt, &b_jetPt);
   fChain->SetBranchAddress("jetEta", &jetEta, &b_jetEta);
   fChain->SetBranchAddress("jetPhi", &jetPhi, &b_jetPhi);
   fChain->SetBranchAddress("jetNHF", &jetNHF, &b_jetNHF);
   fChain->SetBranchAddress("jetNEMF", &jetNEMF, &b_jetNEMF);
   fChain->SetBranchAddress("jetCHF", &jetCHF, &b_jetCHF);
   fChain->SetBranchAddress("jetCEMF", &jetCEMF, &b_jetCEMF);
   fChain->SetBranchAddress("jetMUF", &jetMUF, &b_jetMUF);
   fChain->SetBranchAddress("jetNHM", &jetNHM, &b_jetNHM);
   fChain->SetBranchAddress("jetCHM", &jetCHM, &b_jetCHM);
   fChain->SetBranchAddress("jetPHM", &jetPHM, &b_jetPHM);
   fChain->SetBranchAddress("jetELM", &jetELM, &b_jetELM);
   fChain->SetBranchAddress("jetPHE", &jetPHE, &b_jetPHE);
   fChain->SetBranchAddress("jetPHEF", &jetPHEF, &b_jetPHEF);
   fChain->SetBranchAddress("jetELE", &jetELE, &b_jetELE);
   fChain->SetBranchAddress("jetELEF", &jetELEF, &b_jetELEF);
   fChain->SetBranchAddress("jetMUE", &jetMUE, &b_jetMUE);
   fChain->SetBranchAddress("jetID", &jetID, &b_jetID);
   fChain->SetBranchAddress("cljBcCnt", &cljBcCnt, &b_cljBcCnt);
   fChain->SetBranchAddress("nCaloJets", &nCaloJets, &b_nCaloJets);
   fChain->SetBranchAddress("cljSeedTOFTime", &cljSeedTOFTime, &b_cljSeedTOFTime);
   fChain->SetBranchAddress("cljCMeanTime", &cljCMeanTime, &b_cljCMeanTime);
   fChain->SetBranchAddress("cljBc3dEx", &cljBc3dEx, &b_cljBc3dEx);
   fChain->SetBranchAddress("cljBc3dEy", &cljBc3dEy, &b_cljBc3dEy);
   fChain->SetBranchAddress("cljBc3dEz", &cljBc3dEz, &b_cljBc3dEz);
   fChain->SetBranchAddress("cljBc3dEv", &cljBc3dEv, &b_cljBc3dEv);
   fChain->SetBranchAddress("cljBc3dEslope", &cljBc3dEslope, &b_cljBc3dEslope);
   fChain->SetBranchAddress("cljBc3dEchisp", &cljBc3dEchisp, &b_cljBc3dEchisp);
   fChain->SetBranchAddress("cljBc2dEx", &cljBc2dEx, &b_cljBc2dEx);
   fChain->SetBranchAddress("cljBc2dEy", &cljBc2dEy, &b_cljBc2dEy);
   fChain->SetBranchAddress("cljBc2dEv", &cljBc2dEv, &b_cljBc2dEv);
   fChain->SetBranchAddress("cljBc2dEslope", &cljBc2dEslope, &b_cljBc2dEslope);
   fChain->SetBranchAddress("cljBc2dEchisp", &cljBc2dEchisp, &b_cljBc2dEchisp);
   fChain->SetBranchAddress("cljCDrMeanTime", &cljCDrMeanTime, &b_cljCDrMeanTime);
   fChain->SetBranchAddress("cljPt", &cljPt, &b_cljPt);
   fChain->SetBranchAddress("cljEnergy", &cljEnergy, &b_cljEnergy);
   fChain->SetBranchAddress("cljPhi", &cljPhi, &b_cljPhi);
   fChain->SetBranchAddress("cljEta", &cljEta, &b_cljEta);
   fChain->SetBranchAddress("cljPx", &cljPx, &b_cljPx);
   fChain->SetBranchAddress("cljPy", &cljPy, &b_cljPy);
   fChain->SetBranchAddress("cljPz", &cljPz, &b_cljPz);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("phoSeedTOFTime", &phoSeedTOFTime, &b_phoSeedTOFTime);
   fChain->SetBranchAddress("phoCMeanTime", &phoCMeanTime, &b_phoCMeanTime);
   fChain->SetBranchAddress("phoSc3dEx", &phoSc3dEx, &b_phoSc3dEx);
   fChain->SetBranchAddress("phoSc3dEy", &phoSc3dEy, &b_phoSc3dEy);
   fChain->SetBranchAddress("phoSc3dEz", &phoSc3dEz, &b_phoSc3dEz);
   fChain->SetBranchAddress("phoSc3dEv", &phoSc3dEv, &b_phoSc3dEv);
   fChain->SetBranchAddress("phoSc3dEslope", &phoSc3dEslope, &b_phoSc3dEslope);
   fChain->SetBranchAddress("phoSc3dEchisp", &phoSc3dEchisp, &b_phoSc3dEchisp);
   fChain->SetBranchAddress("phoSc2dEx", &phoSc2dEx, &b_phoSc2dEx);
   fChain->SetBranchAddress("phoSc2dEy", &phoSc2dEy, &b_phoSc2dEy);
   fChain->SetBranchAddress("phoSc2dEv", &phoSc2dEv, &b_phoSc2dEv);
   fChain->SetBranchAddress("phoSc2dEslope", &phoSc2dEslope, &b_phoSc2dEslope);
   fChain->SetBranchAddress("phoSc2dEchisp", &phoSc2dEchisp, &b_phoSc2dEchisp);
   fChain->SetBranchAddress("phoPt", &phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEnergy", &phoEnergy, &b_phoEnergy);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("nOotPhotons", &nOotPhotons, &b_nOotPhotons);
   fChain->SetBranchAddress("ootPhoSeedTOFTime", &ootPhoSeedTOFTime, &b_ootPhoSeedTOFTime);
   fChain->SetBranchAddress("ootPhoCMeanTime", &ootPhoCMeanTime, &b_ootPhoCMeanTime);
   fChain->SetBranchAddress("ootPhoSc3dEx", &ootPhoSc3dEx, &b_ootPhoSc3dEx);
   fChain->SetBranchAddress("ootPhoSc3dEy", &ootPhoSc3dEy, &b_ootPhoSc3dEy);
   fChain->SetBranchAddress("ootPhoSc3dEz", &ootPhoSc3dEz, &b_ootPhoSc3dEz);
   fChain->SetBranchAddress("ootPhoSc3dEv", &ootPhoSc3dEv, &b_ootPhoSc3dEv);
   fChain->SetBranchAddress("ootPhoSc3dEslope", &ootPhoSc3dEslope, &b_ootPhoSc3dEslope);
   fChain->SetBranchAddress("ootPhoSc3dEchisp", &ootPhoSc3dEchisp, &b_ootPhoSc3dEchisp);
   fChain->SetBranchAddress("ootPhoSc2dEx", &ootPhoSc2dEx, &b_ootPhoSc2dEx);
   fChain->SetBranchAddress("ootPhoSc2dEy", &ootPhoSc2dEy, &b_ootPhoSc2dEy);
   fChain->SetBranchAddress("ootPhoSc2dEv", &ootPhoSc2dEv, &b_ootPhoSc2dEv);
   fChain->SetBranchAddress("ootPhoSc2dEslope", &ootPhoSc2dEslope, &b_ootPhoSc2dEslope);
   fChain->SetBranchAddress("ootPhoSc2dEchisp", &ootPhoSc2dEchisp, &b_ootPhoSc2dEchisp);
   fChain->SetBranchAddress("ootPhoPt", &ootPhoPt, &b_ootPhoPt);
   fChain->SetBranchAddress("ootPhoEnergy", &ootPhoEnergy, &b_ootPhoEnergy);
   fChain->SetBranchAddress("ootPhoPhi", &ootPhoPhi, &b_ootPhoPhi);
   fChain->SetBranchAddress("ootPhoEta", &ootPhoEta, &b_ootPhoEta);
   fChain->SetBranchAddress("ootPhoPx", &ootPhoPx, &b_ootPhoPx);
   fChain->SetBranchAddress("ootPhoPy", &ootPhoPy, &b_ootPhoPy);
   fChain->SetBranchAddress("ootPhoPz", &ootPhoPz, &b_ootPhoPz);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("eleSeedTOFTime", &eleSeedTOFTime, &b_eleSeedTOFTime);
   fChain->SetBranchAddress("eleCMeanTime", &eleCMeanTime, &b_eleCMeanTime);
   fChain->SetBranchAddress("eleSc3dEx", &eleSc3dEx, &b_eleSc3dEx);
   fChain->SetBranchAddress("eleSc3dEy", &eleSc3dEy, &b_eleSc3dEy);
   fChain->SetBranchAddress("eleSc3dEz", &eleSc3dEz, &b_eleSc3dEz);
   fChain->SetBranchAddress("eleSc3dEv", &eleSc3dEv, &b_eleSc3dEv);
   fChain->SetBranchAddress("eleSc3dEslope", &eleSc3dEslope, &b_eleSc3dEslope);
   fChain->SetBranchAddress("eleSc3dEchisp", &eleSc3dEchisp, &b_eleSc3dEchisp);
   fChain->SetBranchAddress("eleSc2dEx", &eleSc2dEx, &b_eleSc2dEx);
   fChain->SetBranchAddress("eleSc2dEy", &eleSc2dEy, &b_eleSc2dEy);
   fChain->SetBranchAddress("eleSc2dEv", &eleSc2dEv, &b_eleSc2dEv);
   fChain->SetBranchAddress("eleSc2dEslope", &eleSc2dEslope, &b_eleSc2dEslope);
   fChain->SetBranchAddress("eleSc2dEchisp", &eleSc2dEchisp, &b_eleSc2dEchisp);
   fChain->SetBranchAddress("elePt", &elePt, &b_elePt);
   fChain->SetBranchAddress("eleEnergy", &eleEnergy, &b_eleEnergy);
   fChain->SetBranchAddress("elePhi", &elePhi, &b_elePhi);
   fChain->SetBranchAddress("eleEta", &eleEta, &b_eleEta);
   fChain->SetBranchAddress("elePx", &elePx, &b_elePx);
   fChain->SetBranchAddress("elePy", &elePy, &b_elePy);
   fChain->SetBranchAddress("elePz", &elePz, &b_elePz);
   fChain->SetBranchAddress("jetSumEPFrac", &jetSumEPFrac, &b_jetSumEPFrac);
   fChain->SetBranchAddress("jetEPEnergy", &jetEPEnergy, &b_jetEPEnergy);
   fChain->SetBranchAddress("jetEMEnergy", &jetEMEnergy, &b_jetEMEnergy);
   fChain->SetBranchAddress("jetEMEnrFrac", &jetEMEnrFrac, &b_jetEMEnrFrac);
   fChain->SetBranchAddress("jetEPEnrFrac", &jetEPEnrFrac, &b_jetEPEnrFrac);
   fChain->SetBranchAddress("jetDrLeadEta", &jetDrLeadEta, &b_jetDrLeadEta);
   fChain->SetBranchAddress("jetDrLeadPhi", &jetDrLeadPhi, &b_jetDrLeadPhi);
   fChain->SetBranchAddress("jetDrLeadEnr", &jetDrLeadEnr, &b_jetDrLeadEnr);
   fChain->SetBranchAddress("sJetDrRHEnergy", &sJetDrRHEnergy, &b_sJetDrRHEnergy);
   fChain->SetBranchAddress("jetDrEMF", &jetDrEMF, &b_jetDrEMF);
   fChain->SetBranchAddress("jetDrRhCnt", &jetDrRhCnt, &b_jetDrRhCnt);
   fChain->SetBranchAddress("jetBcTimesCnt", &jetBcTimesCnt, &b_jetBcTimesCnt);
   fChain->SetBranchAddress("jetBcSumRHEnr", &jetBcSumRHEnr, &b_jetBcSumRHEnr);
   fChain->SetBranchAddress("jetBcEMFr", &jetBcEMFr, &b_jetBcEMFr);
   fChain->SetBranchAddress("jetBcRhCnt", &jetBcRhCnt, &b_jetBcRhCnt);
   fChain->SetBranchAddress("jetBcGrpCnt", &jetBcGrpCnt, &b_jetBcGrpCnt);
   fChain->SetBranchAddress("nJetScMatch", &nJetScMatch, &b_nJetScMatch);
   fChain->SetBranchAddress("jetScRhCnt", &jetScRhCnt, &b_jetScRhCnt);
   fChain->SetBranchAddress("sJetScEnergy", &sJetScEnergy, &b_sJetScEnergy);
   fChain->SetBranchAddress("sJetScPhEnergy", &sJetScPhEnergy, &b_sJetScPhEnergy);
   fChain->SetBranchAddress("sJetScRhEnergy", &sJetScRhEnergy, &b_sJetScRhEnergy);
   fChain->SetBranchAddress("jetScEMF", &jetScEMF, &b_jetScEMF);
   fChain->SetBranchAddress("jetDRMuTime", &jetDRMuTime, &b_jetDRMuTime);
   fChain->SetBranchAddress("jetDRTimeError", &jetDRTimeError, &b_jetDRTimeError);
   fChain->SetBranchAddress("jetDRTimeRMS", &jetDRTimeRMS, &b_jetDRTimeRMS);
   fChain->SetBranchAddress("jetDRMedTime", &jetDRMedTime, &b_jetDRMedTime);
   fChain->SetBranchAddress("jetCDRMuTime", &jetCDRMuTime, &b_jetCDRMuTime);
   fChain->SetBranchAddress("jetCDRMedTime", &jetCDRMedTime, &b_jetCDRMedTime);
   fChain->SetBranchAddress("jetSCMuTime", &jetSCMuTime, &b_jetSCMuTime);
   fChain->SetBranchAddress("jetSCMedTime", &jetSCMedTime, &b_jetSCMedTime);
   fChain->SetBranchAddress("jetCSCMuTime", &jetCSCMuTime, &b_jetCSCMuTime);
   fChain->SetBranchAddress("jetCSCMedTime", &jetCSCMedTime, &b_jetCSCMedTime);
   fChain->SetBranchAddress("jetCBCMuTime", &jetCBCMuTime, &b_jetCBCMuTime);
   fChain->SetBranchAddress("jetCBCMedTime", &jetCBCMedTime, &b_jetCBCMedTime);
   fChain->SetBranchAddress("jetGenImpactAngle", &jetGenImpactAngle, &b_jetGenImpactAngle);
   fChain->SetBranchAddress("jetGenTime", &jetGenTime, &b_jetGenTime);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenEnergy", &jetGenEnergy, &b_jetGenEnergy);
   fChain->SetBranchAddress("jetGenEMFrac", &jetGenEMFrac, &b_jetGenEMFrac);
   fChain->SetBranchAddress("jetGenDrMatch", &jetGenDrMatch, &b_jetGenDrMatch);
   fChain->SetBranchAddress("jetGenTimeVar", &jetGenTimeVar, &b_jetGenTimeVar);
   fChain->SetBranchAddress("jetGenTimeLLP", &jetGenTimeLLP, &b_jetGenTimeLLP);
   fChain->SetBranchAddress("jetGenLLPPurity", &jetGenLLPPurity, &b_jetGenLLPPurity);
   fChain->SetBranchAddress("jetGenNextBX", &jetGenNextBX, &b_jetGenNextBX);
   fChain->SetBranchAddress("jetGenNKids", &jetGenNKids, &b_jetGenNKids);
   fChain->SetBranchAddress("jetGenTOF", &jetGenTOF, &b_jetGenTOF);
   fChain->SetBranchAddress("jetImpactAngle", &jetImpactAngle, &b_jetImpactAngle);
   fChain->SetBranchAddress("jetSc3dEx", &jetSc3dEx, &b_jetSc3dEx);
   fChain->SetBranchAddress("jetSc3dEy", &jetSc3dEy, &b_jetSc3dEy);
   fChain->SetBranchAddress("jetSc3dEz", &jetSc3dEz, &b_jetSc3dEz);
   fChain->SetBranchAddress("jetSc3dEv", &jetSc3dEv, &b_jetSc3dEv);
   fChain->SetBranchAddress("jetSc3dEslope", &jetSc3dEslope, &b_jetSc3dEslope);
   fChain->SetBranchAddress("jetSc3dEchisp", &jetSc3dEchisp, &b_jetSc3dEchisp);
   fChain->SetBranchAddress("jetSc2dEx", &jetSc2dEx, &b_jetSc2dEx);
   fChain->SetBranchAddress("jetSc2dEy", &jetSc2dEy, &b_jetSc2dEy);
   fChain->SetBranchAddress("jetSc2dEv", &jetSc2dEv, &b_jetSc2dEv);
   fChain->SetBranchAddress("jetSc2dEslope", &jetSc2dEslope, &b_jetSc2dEslope);
   fChain->SetBranchAddress("jetSc2dEchisp", &jetSc2dEchisp, &b_jetSc2dEchisp);
   fChain->SetBranchAddress("jetSc2dEslope2", &jetSc2dEslope2, &b_jetSc2dEslope2);
   fChain->SetBranchAddress("jetSc2dEchisp2", &jetSc2dEchisp2, &b_jetSc2dEchisp2);
   fChain->SetBranchAddress("jetSc2dErangle", &jetSc2dErangle, &b_jetSc2dErangle);
   fChain->SetBranchAddress("jetSc2dEnxsum", &jetSc2dEnxsum, &b_jetSc2dEnxsum);
   fChain->SetBranchAddress("nRecHits", &nRecHits, &b_nRecHits);
   fChain->SetBranchAddress("rhPosX", &rhPosX, &b_rhPosX);
   fChain->SetBranchAddress("rhPosY", &rhPosY, &b_rhPosY);
   fChain->SetBranchAddress("rhPosZ", &rhPosZ, &b_rhPosZ);
   fChain->SetBranchAddress("rhPosEta", &rhPosEta, &b_rhPosEta);
   fChain->SetBranchAddress("rhPosPhi", &rhPosPhi, &b_rhPosPhi);
   fChain->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);
   fChain->SetBranchAddress("rhTime", &rhTime, &b_rhTime);
   fChain->SetBranchAddress("rhTimeErr", &rhTimeErr, &b_rhTimeErr);
   fChain->SetBranchAddress("rhTOF", &rhTOF, &b_rhTOF);
   fChain->SetBranchAddress("rhID", &rhID, &b_rhID);
   fChain->SetBranchAddress("rhXtalI1", &rhXtalI1, &b_rhXtalI1);
   fChain->SetBranchAddress("rhXtalI2", &rhXtalI2, &b_rhXtalI2);
   fChain->SetBranchAddress("rhSubdet", &rhSubdet, &b_rhSubdet);
   fChain->SetBranchAddress("rhisOOT", &rhisOOT, &b_rhisOOT);
   fChain->SetBranchAddress("rhisGS6", &rhisGS6, &b_rhisGS6);
   fChain->SetBranchAddress("rhisGS1", &rhisGS1, &b_rhisGS1);
   fChain->SetBranchAddress("rhadcToGeV", &rhadcToGeV, &b_rhadcToGeV);
   fChain->SetBranchAddress("rhped12", &rhped12, &b_rhped12);
   fChain->SetBranchAddress("rhped6", &rhped6, &b_rhped6);
   fChain->SetBranchAddress("rhped1", &rhped1, &b_rhped1);
   fChain->SetBranchAddress("rhpedrms12", &rhpedrms12, &b_rhpedrms12);
   fChain->SetBranchAddress("rhpedrms6", &rhpedrms6, &b_rhpedrms6);
   fChain->SetBranchAddress("rhpedrms1", &rhpedrms1, &b_rhpedrms1);
   //Notify();
}

/*
Bool_t llpgana_ntuple_base::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void llpgana_ntuple_base::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t llpgana_ntuple_base::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/
//#endif // #ifdef llpgana_ntuple_base_cxx
