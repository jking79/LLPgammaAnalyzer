//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  7 17:06:18 2022 by ROOT version 6.14/09
// from TTree llpgtree/llpgtree
// found on file: llpgana_mc_AODSIM_ntuplizer_test_v16.root
//////////////////////////////////////////////////////////

//#ifndef llpgana_hist_rebase_h
//#define llpgana_hist_rebase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

using namespace std;

class llpgana_hist_rebase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
//   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Int_t           nVtx;
   Float_t         vtxX;
   Float_t         vtxY;
   Float_t         vtxZ;
   Float_t         metSumEt;
   Float_t         metPt;
   Float_t         metPx;
   Float_t         metPy;
   Float_t         metPhi;
   Float_t         metCSumEt;
   Float_t         metCPx;
   Float_t         metCPy;
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
   vector<float>   *cljCDrMeanTime;
   vector<float>   *cljPt;
   vector<float>   *cljEnergy;
   vector<float>   *cljPhi;
   vector<float>   *cljEta;
   vector<float>   *cljPx;
   vector<float>   *cljPy;
   vector<float>   *cljPz;
   UInt_t          nPhotons;
   vector<bool>    *phoIsOotPho;
   vector<bool>    *phoExcluded;
   vector<float>   *phoSeedTOFTime;
   vector<float>   *phoCMeanTime;
   vector<float>   *phoSc2dEv;
   vector<float>   *phoPt;
   vector<float>   *phoEnergy;
   vector<float>   *phoPhi;
   vector<float>   *phoEta;
   vector<float>   *phoPx;
   vector<float>   *phoPy;
   vector<float>   *phoPz;
   vector<vector<unsigned int> > *phoRhIds;
   vector<bool>    *phoIsPFPhoton;
   vector<bool>    *phoIsStdPhoton;
   vector<bool>    *phoHasConTracks;
   vector<bool>    *phoIsPixelSeed;
   vector<bool>    *phoIsPhoton;
   vector<bool>    *phoIsEB;
   vector<bool>    *phoIsEE;
   vector<float>   *phoHadOverEM;
   vector<float>   *phoHadD1OverEM;
   vector<float>   *phoHadD2OverEM;
   vector<float>   *phoHadOverEMVaid;
   vector<float>   *phohadTowOverEM;
   vector<float>   *phohadTowD10OverEM;
   vector<float>   *phohadTowD20OverEM;
   vector<float>   *phohadTowOverEMValid;
   vector<float>   *phoE1x5;
   vector<float>   *phoE2x5;
   vector<float>   *phoE3x3;
   vector<float>   *phoE5x5;
   vector<float>   *phoMaxEnergyXtal;
   vector<float>   *phoSigmaEtaEta;
   vector<float>   *phoSigmaIEtaIEta;
   vector<float>   *phoR1x5;
   vector<float>   *phoR2x5;
   vector<float>   *phoR9;
   vector<float>   *phoFull5x5_e1x5;
   vector<float>   *phoFull5x5_e2x5;
   vector<float>   *phoFull5x5_e3x3;
   vector<float>   *phoFull5x5_e5x5;
   vector<float>   *phoFull5x5_maxEnergyXtal;
   vector<float>   *phoFull5x5_sigmaEtaEta;
   vector<float>   *phoFull5x5_sigmaIEtaIEta;
   vector<float>   *phoFull5x5_r9;
   vector<float>   *phoEcalRHSumEtConeDR04;
   vector<float>   *phoHcalTwrSumEtConeDR04;
   vector<float>   *phoHcalDepth1TowerSumEtConeDR04;
   vector<float>   *phoCalDepth2TowerSumEtConeDR04;
   vector<float>   *phoHcalTowerSumEtBcConeDR04;
   vector<float>   *phoHcalDepth1TowerSumEtBcConeDR04;
   vector<float>   *phoHcalDepth2TowerSumEtBcConeDR04;
   vector<float>   *phoTrkSumPtSolidConeDR04;
   vector<float>   *phoTrkSumPtHollowConeDR04;
   vector<int>     *phoNTrkSolidConeDR04;
   vector<int>     *phoNTrkHollowConeDR04;
   vector<int>     *genPhoIdx;
   vector<float>   *genPhoDr;
   UInt_t          nElectrons;
   vector<float>   *eleSeedTOFTime;
   vector<float>   *eleCMeanTime;
   vector<float>   *elePt;
   vector<float>   *eleEnergy;
   vector<float>   *elePhi;
   vector<float>   *eleEta;
   vector<float>   *elePx;
   vector<float>   *elePy;
   vector<float>   *elePz;
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
   vector<vector<unsigned int> > *jetScRhIds;
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
   vector<float>   *jetSc2dEv;
   UInt_t          nGenParts;
   vector<float>   *genPt;
   vector<float>   *genEnergy;
   vector<float>   *genPhi;
   vector<float>   *genEta;
   vector<float>   *genPx;
   vector<float>   *genPy;
   vector<float>   *genPz;
   vector<int>     *genPdgId;
   vector<int>     *genLLP;
   Int_t           nRecHits;
   vector<float>   *rhPosX;
   vector<float>   *rhPosY;
   vector<float>   *rhPosZ;
   vector<float>   *rhPosEta;
   vector<float>   *rhPosPhi;
   vector<float>   *rhEnergy;
   vector<float>   *rhTime;
   vector<float>   *rhTOF;
   vector<unsigned int> *rhID;
   vector<int>     *rhXtalI1;
   vector<int>     *rhXtalI2;
   vector<int>     *rhSubdet;
   vector<bool>    *rhisOOT;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_nVtx;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_metSumEt;   //!
   TBranch        *b_metPt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_metPhi;   //!
   TBranch        *b_metCSumEt;   //!
   TBranch        *b_metCPx;   //!
   TBranch        *b_metCPy;   //!
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
   TBranch        *b_cljCDrMeanTime;   //!
   TBranch        *b_cljPt;   //!
   TBranch        *b_cljEnergy;   //!
   TBranch        *b_cljPhi;   //!
   TBranch        *b_cljEta;   //!
   TBranch        *b_cljPx;   //!
   TBranch        *b_cljPy;   //!
   TBranch        *b_cljPz;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_phoIsOotPho;   //!
   TBranch        *b_phoExcluded;   //!
   TBranch        *b_phoSeedTOFTime;   //!
   TBranch        *b_phoCMeanTime;   //!
   TBranch        *b_phoSc2dEv;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEnergy;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoRhIds;   //!
   TBranch        *b_phoIsPFPhoton;   //!
   TBranch        *b_phoIsStdPhoton;   //!
   TBranch        *b_phoHasConTracks;   //!
   TBranch        *b_phoIsPixelSeed;   //!
   TBranch        *b_phoIsPhoton;   //!
   TBranch        *b_phoIsEB;   //!
   TBranch        *b_phoIsEE;   //!
   TBranch        *b_phoHadOverEM;   //!
   TBranch        *b_phoHadD1OverEM;   //!
   TBranch        *b_phoHadD2OverEM;   //!
   TBranch        *b_phoHadOverEMVaid;   //!
   TBranch        *b_phohadTowOverEM;   //!
   TBranch        *b_phohadTowD10OverEM;   //!
   TBranch        *b_phohadTowD20OverEM;   //!
   TBranch        *b_phohadTowOverEMValid;   //!
   TBranch        *b_phoE1x5;   //!
   TBranch        *b_phoE2x5;   //!
   TBranch        *b_phoE3x3;   //!
   TBranch        *b_phoE5x5;   //!
   TBranch        *b_phoMaxEnergyXtal;   //!
   TBranch        *b_phoSigmaEtaEta;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoR1x5;   //!
   TBranch        *b_phoR2x5;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoFull5x5_e1x5;   //!
   TBranch        *b_phoFull5x5_e2x5;   //!
   TBranch        *b_phoFull5x5_e3x3;   //!
   TBranch        *b_phoFull5x5_e5x5;   //!
   TBranch        *b_phoFull5x5_maxEnergyXtal;   //!
   TBranch        *b_phoFull5x5_sigmaEtaEta;   //!
   TBranch        *b_phoFull5x5_sigmaIEtaIEta;   //!
   TBranch        *b_phoFull5x5_r9;   //!
   TBranch        *b_phoEcalRHSumEtConeDR04;   //!
   TBranch        *b_phoHcalTwrSumEtConeDR04;   //!
   TBranch        *b_phoHcalDepth1TowerSumEtConeDR04;   //!
   TBranch        *b_phoCalDepth2TowerSumEtConeDR04;   //!
   TBranch        *b_phoHcalTowerSumEtBcConeDR04;   //!
   TBranch        *b_phoHcalDepth1TowerSumEtBcConeDR04;   //!
   TBranch        *b_phoHcalDepth2TowerSumEtBcConeDR04;   //!
   TBranch        *b_phoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_phoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_phoNTrkSolidConeDR04;   //!
   TBranch        *b_phoNTrkHollowConeDR04;   //!
   TBranch        *b_genPhoIdx;   //!
   TBranch        *b_genPhoDr;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_eleSeedTOFTime;   //!
   TBranch        *b_eleCMeanTime;   //!
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
   TBranch        *b_jetScRhIds;   //!
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
   TBranch        *b_jetSc2dEv;   //!
   TBranch        *b_nGenParts;   //!
   TBranch        *b_genPt;   //!
   TBranch        *b_genEnergy;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPx;   //!
   TBranch        *b_genPy;   //!
   TBranch        *b_genPz;   //!
   TBranch        *b_genPdgId;   //!
   TBranch        *b_genLLP;   //!
   TBranch        *b_nRecHits;   //!
   TBranch        *b_rhPosX;   //!
   TBranch        *b_rhPosY;   //!
   TBranch        *b_rhPosZ;   //!
   TBranch        *b_rhPosEta;   //!
   TBranch        *b_rhPosPhi;   //!
   TBranch        *b_rhEnergy;   //!
   TBranch        *b_rhTime;   //!
   TBranch        *b_rhTOF;   //!
   TBranch        *b_rhID;   //!
   TBranch        *b_rhXtalI1;   //!
   TBranch        *b_rhXtalI2;   //!
   TBranch        *b_rhSubdet;   //!
   TBranch        *b_rhisOOT;   //!

//  llpgana_hist_rebase(TTree *tree=0);
//  virtual ~llpgana_hist_rebase();
//   virtual Int_t    Cut(Long64_t entry);
//   virtual Int_t    GetEntry(Long64_t entry);
//   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
//   virtual void     Loop();
//   virtual Bool_t   Notify();
//   virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef llpgana_hist_rebase_cxx
/*
llpgana_hist_rebase::llpgana_hist_rebase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("llpgana_mc_AODSIM_ntuplizer_test_v16.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("llpgana_mc_AODSIM_ntuplizer_test_v16.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("llpgana_mc_AODSIM_ntuplizer_test_v16.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

llpgana_hist_rebase::~llpgana_hist_rebase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t llpgana_hist_rebase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t llpgana_hist_rebase::LoadTree(Long64_t entry)
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

void llpgana_hist_rebase::Init(TTree *tree)
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
   cljCDrMeanTime = 0;
   cljPt = 0;
   cljEnergy = 0;
   cljPhi = 0;
   cljEta = 0;
   cljPx = 0;
   cljPy = 0;
   cljPz = 0;
   phoIsOotPho = 0;
   phoExcluded = 0;
   phoSeedTOFTime = 0;
   phoCMeanTime = 0;
   phoSc2dEv = 0;
   phoPt = 0;
   phoEnergy = 0;
   phoPhi = 0;
   phoEta = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   phoRhIds = 0;
   phoIsPFPhoton = 0;
   phoIsStdPhoton = 0;
   phoHasConTracks = 0;
   phoIsPixelSeed = 0;
   phoIsPhoton = 0;
   phoIsEB = 0;
   phoIsEE = 0;
   phoHadOverEM = 0;
   phoHadD1OverEM = 0;
   phoHadD2OverEM = 0;
   phoHadOverEMVaid = 0;
   phohadTowOverEM = 0;
   phohadTowD10OverEM = 0;
   phohadTowD20OverEM = 0;
   phohadTowOverEMValid = 0;
   phoE1x5 = 0;
   phoE2x5 = 0;
   phoE3x3 = 0;
   phoE5x5 = 0;
   phoMaxEnergyXtal = 0;
   phoSigmaEtaEta = 0;
   phoSigmaIEtaIEta = 0;
   phoR1x5 = 0;
   phoR2x5 = 0;
   phoR9 = 0;
   phoFull5x5_e1x5 = 0;
   phoFull5x5_e2x5 = 0;
   phoFull5x5_e3x3 = 0;
   phoFull5x5_e5x5 = 0;
   phoFull5x5_maxEnergyXtal = 0;
   phoFull5x5_sigmaEtaEta = 0;
   phoFull5x5_sigmaIEtaIEta = 0;
   phoFull5x5_r9 = 0;
   phoEcalRHSumEtConeDR04 = 0;
   phoHcalTwrSumEtConeDR04 = 0;
   phoHcalDepth1TowerSumEtConeDR04 = 0;
   phoCalDepth2TowerSumEtConeDR04 = 0;
   phoHcalTowerSumEtBcConeDR04 = 0;
   phoHcalDepth1TowerSumEtBcConeDR04 = 0;
   phoHcalDepth2TowerSumEtBcConeDR04 = 0;
   phoTrkSumPtSolidConeDR04 = 0;
   phoTrkSumPtHollowConeDR04 = 0;
   phoNTrkSolidConeDR04 = 0;
   phoNTrkHollowConeDR04 = 0;
   genPhoIdx = 0;
   genPhoDr = 0;
   eleSeedTOFTime = 0;
   eleCMeanTime = 0;
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
   jetScRhIds = 0;
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
   jetSc2dEv = 0;
   genPt = 0;
   genEnergy = 0;
   genPhi = 0;
   genEta = 0;
   genPx = 0;
   genPy = 0;
   genPz = 0;
   genPdgId = 0;
   genLLP = 0;
   rhPosX = 0;
   rhPosY = 0;
   rhPosZ = 0;
   rhPosEta = 0;
   rhPosPhi = 0;
   rhEnergy = 0;
   rhTime = 0;
   rhTOF = 0;
   rhID = 0;
   rhXtalI1 = 0;
   rhXtalI2 = 0;
   rhSubdet = 0;
   rhisOOT = 0;
   // Set branch addresses and branch pointers
//  if (!tree) return;
   fChain = tree;
//  fCurrent = -1;
//  fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("nVtx", &nVtx, &b_nVtx);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metCSumEt", &metCSumEt, &b_metCSumEt);
   fChain->SetBranchAddress("metCPx", &metCPx, &b_metCPx);
   fChain->SetBranchAddress("metCPy", &metCPy, &b_metCPy);
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
   fChain->SetBranchAddress("cljCDrMeanTime", &cljCDrMeanTime, &b_cljCDrMeanTime);
   fChain->SetBranchAddress("cljPt", &cljPt, &b_cljPt);
   fChain->SetBranchAddress("cljEnergy", &cljEnergy, &b_cljEnergy);
   fChain->SetBranchAddress("cljPhi", &cljPhi, &b_cljPhi);
   fChain->SetBranchAddress("cljEta", &cljEta, &b_cljEta);
   fChain->SetBranchAddress("cljPx", &cljPx, &b_cljPx);
   fChain->SetBranchAddress("cljPy", &cljPy, &b_cljPy);
   fChain->SetBranchAddress("cljPz", &cljPz, &b_cljPz);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("phoIsOotPho", &phoIsOotPho, &b_phoIsOotPho);
   fChain->SetBranchAddress("phoExcluded", &phoExcluded, &b_phoExcluded);
   fChain->SetBranchAddress("phoSeedTOFTime", &phoSeedTOFTime, &b_phoSeedTOFTime);
   fChain->SetBranchAddress("phoCMeanTime", &phoCMeanTime, &b_phoCMeanTime);
   fChain->SetBranchAddress("phoSc2dEv", &phoSc2dEv, &b_phoSc2dEv);
   fChain->SetBranchAddress("phoPt", &phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEnergy", &phoEnergy, &b_phoEnergy);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoRhIds", &phoRhIds, &b_phoRhIds);
   fChain->SetBranchAddress("phoIsPFPhoton", &phoIsPFPhoton, &b_phoIsPFPhoton);
   fChain->SetBranchAddress("phoIsStdPhoton", &phoIsStdPhoton, &b_phoIsStdPhoton);
   fChain->SetBranchAddress("phoHasConTracks", &phoHasConTracks, &b_phoHasConTracks);
   fChain->SetBranchAddress("phoIsPixelSeed", &phoIsPixelSeed, &b_phoIsPixelSeed);
   fChain->SetBranchAddress("phoIsPhoton", &phoIsPhoton, &b_phoIsPhoton);
   fChain->SetBranchAddress("phoIsEB", &phoIsEB, &b_phoIsEB);
   fChain->SetBranchAddress("phoIsEE", &phoIsEE, &b_phoIsEE);
   fChain->SetBranchAddress("phoHadOverEM", &phoHadOverEM, &b_phoHadOverEM);
   fChain->SetBranchAddress("phoHadD1OverEM", &phoHadD1OverEM, &b_phoHadD1OverEM);
   fChain->SetBranchAddress("phoHadD2OverEM", &phoHadD2OverEM, &b_phoHadD2OverEM);
   fChain->SetBranchAddress("phoHadOverEMVaid", &phoHadOverEMVaid, &b_phoHadOverEMVaid);
   fChain->SetBranchAddress("phohadTowOverEM", &phohadTowOverEM, &b_phohadTowOverEM);
   fChain->SetBranchAddress("phohadTowD10OverEM", &phohadTowD10OverEM, &b_phohadTowD10OverEM);
   fChain->SetBranchAddress("phohadTowD20OverEM", &phohadTowD20OverEM, &b_phohadTowD20OverEM);
   fChain->SetBranchAddress("phohadTowOverEMValid", &phohadTowOverEMValid, &b_phohadTowOverEMValid);
   fChain->SetBranchAddress("phoE1x5", &phoE1x5, &b_phoE1x5);
   fChain->SetBranchAddress("phoE2x5", &phoE2x5, &b_phoE2x5);
   fChain->SetBranchAddress("phoE3x3", &phoE3x3, &b_phoE3x3);
   fChain->SetBranchAddress("phoE5x5", &phoE5x5, &b_phoE5x5);
   fChain->SetBranchAddress("phoMaxEnergyXtal", &phoMaxEnergyXtal, &b_phoMaxEnergyXtal);
   fChain->SetBranchAddress("phoSigmaEtaEta", &phoSigmaEtaEta, &b_phoSigmaEtaEta);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoR1x5", &phoR1x5, &b_phoR1x5);
   fChain->SetBranchAddress("phoR2x5", &phoR2x5, &b_phoR2x5);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoFull5x5_e1x5", &phoFull5x5_e1x5, &b_phoFull5x5_e1x5);
   fChain->SetBranchAddress("phoFull5x5_e2x5", &phoFull5x5_e2x5, &b_phoFull5x5_e2x5);
   fChain->SetBranchAddress("phoFull5x5_e3x3", &phoFull5x5_e3x3, &b_phoFull5x5_e3x3);
   fChain->SetBranchAddress("phoFull5x5_e5x5", &phoFull5x5_e5x5, &b_phoFull5x5_e5x5);
   fChain->SetBranchAddress("phoFull5x5_maxEnergyXtal", &phoFull5x5_maxEnergyXtal, &b_phoFull5x5_maxEnergyXtal);
   fChain->SetBranchAddress("phoFull5x5_sigmaEtaEta", &phoFull5x5_sigmaEtaEta, &b_phoFull5x5_sigmaEtaEta);
   fChain->SetBranchAddress("phoFull5x5_sigmaIEtaIEta", &phoFull5x5_sigmaIEtaIEta, &b_phoFull5x5_sigmaIEtaIEta);
   fChain->SetBranchAddress("phoFull5x5_r9", &phoFull5x5_r9, &b_phoFull5x5_r9);
   fChain->SetBranchAddress("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04, &b_phoEcalRHSumEtConeDR04);
   fChain->SetBranchAddress("phoHcalTwrSumEtConeDR04", &phoHcalTwrSumEtConeDR04, &b_phoHcalTwrSumEtConeDR04);
   fChain->SetBranchAddress("phoHcalDepth1TowerSumEtConeDR04", &phoHcalDepth1TowerSumEtConeDR04, &b_phoHcalDepth1TowerSumEtConeDR04);
   fChain->SetBranchAddress("phoCalDepth2TowerSumEtConeDR04", &phoCalDepth2TowerSumEtConeDR04, &b_phoCalDepth2TowerSumEtConeDR04);
   fChain->SetBranchAddress("phoHcalTowerSumEtBcConeDR04", &phoHcalTowerSumEtBcConeDR04, &b_phoHcalTowerSumEtBcConeDR04);
   fChain->SetBranchAddress("phoHcalDepth1TowerSumEtBcConeDR04", &phoHcalDepth1TowerSumEtBcConeDR04, &b_phoHcalDepth1TowerSumEtBcConeDR04);
   fChain->SetBranchAddress("phoHcalDepth2TowerSumEtBcConeDR04", &phoHcalDepth2TowerSumEtBcConeDR04, &b_phoHcalDepth2TowerSumEtBcConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04, &b_phoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04, &b_phoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("phoNTrkSolidConeDR04", &phoNTrkSolidConeDR04, &b_phoNTrkSolidConeDR04);
   fChain->SetBranchAddress("phoNTrkHollowConeDR04", &phoNTrkHollowConeDR04, &b_phoNTrkHollowConeDR04);
   fChain->SetBranchAddress("genPhoIdx", &genPhoIdx, &b_genPhoIdx);
   fChain->SetBranchAddress("genPhoDr", &genPhoDr, &b_genPhoDr);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("eleSeedTOFTime", &eleSeedTOFTime, &b_eleSeedTOFTime);
   fChain->SetBranchAddress("eleCMeanTime", &eleCMeanTime, &b_eleCMeanTime);
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
   fChain->SetBranchAddress("jetScRhIds", &jetScRhIds, &b_jetScRhIds);
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
   fChain->SetBranchAddress("jetSc2dEv", &jetSc2dEv, &b_jetSc2dEv);
   fChain->SetBranchAddress("nGenParts", &nGenParts, &b_nGenParts);
   fChain->SetBranchAddress("genPt", &genPt, &b_genPt);
   fChain->SetBranchAddress("genEnergy", &genEnergy, &b_genEnergy);
   fChain->SetBranchAddress("genPhi", &genPhi, &b_genPhi);
   fChain->SetBranchAddress("genEta", &genEta, &b_genEta);
   fChain->SetBranchAddress("genPx", &genPx, &b_genPx);
   fChain->SetBranchAddress("genPy", &genPy, &b_genPy);
   fChain->SetBranchAddress("genPz", &genPz, &b_genPz);
   fChain->SetBranchAddress("genPdgId", &genPdgId, &b_genPdgId);
   fChain->SetBranchAddress("genLLP", &genLLP, &b_genLLP);
   fChain->SetBranchAddress("nRecHits", &nRecHits, &b_nRecHits);
   fChain->SetBranchAddress("rhPosX", &rhPosX, &b_rhPosX);
   fChain->SetBranchAddress("rhPosY", &rhPosY, &b_rhPosY);
   fChain->SetBranchAddress("rhPosZ", &rhPosZ, &b_rhPosZ);
   fChain->SetBranchAddress("rhPosEta", &rhPosEta, &b_rhPosEta);
   fChain->SetBranchAddress("rhPosPhi", &rhPosPhi, &b_rhPosPhi);
   fChain->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);
   fChain->SetBranchAddress("rhTime", &rhTime, &b_rhTime);
   fChain->SetBranchAddress("rhTOF", &rhTOF, &b_rhTOF);
   fChain->SetBranchAddress("rhID", &rhID, &b_rhID);
   fChain->SetBranchAddress("rhXtalI1", &rhXtalI1, &b_rhXtalI1);
   fChain->SetBranchAddress("rhXtalI2", &rhXtalI2, &b_rhXtalI2);
   fChain->SetBranchAddress("rhSubdet", &rhSubdet, &b_rhSubdet);
   fChain->SetBranchAddress("rhisOOT", &rhisOOT, &b_rhisOOT);
//   Notify();
}

/*
Bool_t llpgana_hist_rebase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void llpgana_hist_rebase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t llpgana_hist_rebase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/
//#endif // #ifdef llpgana_hist_rebase_cxx
