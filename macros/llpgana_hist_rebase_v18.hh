//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 22 13:42:18 2023 by ROOT version 6.14/09
// from TTree llpgtree/llpgtree
// found on file: llpgana_mc_AODSIM_ntuplizer_EB_cl2mom_v18.root
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
   Float_t         metSumEt;
   Float_t         metPt;
   Float_t         metPx;
   Float_t         metPy;
   Float_t         metPhi;
   Float_t         metCSumEt;
   Float_t         metCPx;
   Float_t         metCPy;
   UInt_t          nPhotons;
   vector<bool>    *phoIsOotPho;
   vector<bool>    *phoExcluded;
   vector<float>   *phoSeedTOFTime;
   vector<float>   *phoCMeanTime;
   vector<float>   *phoPt;
   vector<float>   *phoEnergy;
   vector<float>   *phoPhi;
   vector<float>   *phoEta;
   vector<float>   *phoPx;
   vector<float>   *phoPy;
   vector<float>   *phoPz;
   vector<vector<unsigned int> > *phoRhIds;
   vector<bool>    *phoHasConTracks;
   vector<bool>    *phoIsPixelSeed;
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
   vector<float>   *phoSMaj;
   vector<float>   *phoSMin;
   vector<float>   *phoSAlp;
   vector<float>   *phoCovEtaEta;
   vector<float>   *phoCovEtaPhi;
   vector<float>   *phoCovPhiPhi;
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
   TBranch        *b_nPhotons;   //!
   TBranch        *b_phoIsOotPho;   //!
   TBranch        *b_phoExcluded;   //!
   TBranch        *b_phoSeedTOFTime;   //!
   TBranch        *b_phoCMeanTime;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEnergy;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoRhIds;   //!
   TBranch        *b_phoHasConTracks;   //!
   TBranch        *b_phoIsPixelSeed;   //!
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
   TBranch        *b_phoSMaj;   //!
   TBranch        *b_phoSMin;   //!
   TBranch        *b_phoSAlp;   //!
   TBranch        *b_phoCovEtaEta;   //!
   TBranch        *b_phoCovEtaPhi;   //!
   TBranch        *b_phoCovPhiPhi;   //!
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
   TBranch        *b_rhisOOT;   //!

   //llpgana_hist_rebase(TTree *tree=0);
   //virtual ~llpgana_hist_rebase();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef llpgana_hist_rebase_cxx
/*
llpgana_hist_rebase::llpgana_hist_rebase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("llpgana_mc_AODSIM_ntuplizer_EB_cl2mom_v18.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("llpgana_mc_AODSIM_ntuplizer_EB_cl2mom_v18.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("llpgana_mc_AODSIM_ntuplizer_EB_cl2mom_v18.root:/tree");
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
   phoIsOotPho = 0;
   phoExcluded = 0;
   phoSeedTOFTime = 0;
   phoCMeanTime = 0;
   phoPt = 0;
   phoEnergy = 0;
   phoPhi = 0;
   phoEta = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   phoRhIds = 0;
   phoHasConTracks = 0;
   phoIsPixelSeed = 0;
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
   phoSMaj = 0;
   phoSMin = 0;
   phoSAlp = 0;
   phoCovEtaEta = 0;
   phoCovEtaPhi = 0;
   phoCovPhiPhi = 0;
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
   rhisOOT = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
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
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt);
   fChain->SetBranchAddress("metPt", &metPt, &b_metPt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("metPhi", &metPhi, &b_metPhi);
   fChain->SetBranchAddress("metCSumEt", &metCSumEt, &b_metCSumEt);
   fChain->SetBranchAddress("metCPx", &metCPx, &b_metCPx);
   fChain->SetBranchAddress("metCPy", &metCPy, &b_metCPy);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("phoIsOotPho", &phoIsOotPho, &b_phoIsOotPho);
   fChain->SetBranchAddress("phoExcluded", &phoExcluded, &b_phoExcluded);
   fChain->SetBranchAddress("phoSeedTOFTime", &phoSeedTOFTime, &b_phoSeedTOFTime);
   fChain->SetBranchAddress("phoCMeanTime", &phoCMeanTime, &b_phoCMeanTime);
   fChain->SetBranchAddress("phoPt", &phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEnergy", &phoEnergy, &b_phoEnergy);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoRhIds", &phoRhIds, &b_phoRhIds);
   fChain->SetBranchAddress("phoHasConTracks", &phoHasConTracks, &b_phoHasConTracks);
   fChain->SetBranchAddress("phoIsPixelSeed", &phoIsPixelSeed, &b_phoIsPixelSeed);
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
   fChain->SetBranchAddress("phoSMaj", &phoSMaj, &b_phoSMaj);
   fChain->SetBranchAddress("phoSMin", &phoSMin, &b_phoSMin);
   fChain->SetBranchAddress("phoSAlp", &phoSAlp, &b_phoSAlp);
   fChain->SetBranchAddress("phoCovEtaEta", &phoCovEtaEta, &b_phoCovEtaEta);
   fChain->SetBranchAddress("phoCovEtaPhi", &phoCovEtaPhi, &b_phoCovEtaPhi);
   fChain->SetBranchAddress("phoCovPhiPhi", &phoCovPhiPhi, &b_phoCovPhiPhi);
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
   fChain->SetBranchAddress("rhisOOT", &rhisOOT, &b_rhisOOT);
   //Notify();
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
