//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 17 17:28:00 2023 by ROOT version 6.26/07
// from TTree llpgtree/llpgtree
// found on file: gmsb_AODSIM_KUCMSNtuplizer_v1.root
//////////////////////////////////////////////////////////

#ifndef root_rebase_h
#define root_rebase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class root_rebase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Float_t         vtxX;
   Float_t         vtxY;
   Float_t         vtxZ;
   Float_t         metSumEt;
   Float_t         metPx;
   Float_t         metPy;
   Float_t         metCSumEt;
   Float_t         metCPx;
   Float_t         metCPy;
   Float_t         jetHt;
   std::vector<float>   *jetE;
   std::vector<float>   *jetM;
   std::vector<float>   *jetPt;
   std::vector<float>   *jetEta;
   std::vector<float>   *jetPhi;
   std::vector<float>   *jetNHF;
   std::vector<float>   *jetNEMF;
   std::vector<float>   *jetCHF;
   std::vector<float>   *jetCEMF;
   std::vector<float>   *jetMUF;
   std::vector<float>   *jetNHM;
   std::vector<float>   *jetCHM;
   std::vector<int>     *jetPHM;
   std::vector<int>     *jetParts;
   std::vector<std::vector<unsigned int> > *jetDrRhIds;
   std::vector<float>   *jetGenImpactAngle;
   std::vector<float>   *jetGenTime;
   std::vector<float>   *jetGenPt;
   std::vector<float>   *jetGenEta;
   std::vector<float>   *jetGenEnergy;
   std::vector<float>   *jetGenDrMatch;
   std::vector<float>   *jetGenTimeLLP;
   std::vector<float>   *jetGenTOF;
   std::vector<bool>    *phoIsOotPho;
   std::vector<bool>    *phoExcluded;
   std::vector<float>   *phoSeedTOFTime;
   std::vector<float>   *phoPt;
   std::vector<float>   *phoEnergy;
   std::vector<float>   *phoPhi;
   std::vector<float>   *phoEta;
   std::vector<float>   *phoPx;
   std::vector<float>   *phoPy;
   std::vector<float>   *phoPz;
   std::vector<std::vector<unsigned int> > *phoRhIds;
   std::vector<bool>    *phoIsPixelSeed;
   std::vector<bool>    *phoIsEB;
   std::vector<float>   *phohadTowOverEM;
   std::vector<float>   *phoSigmaIEtaIEta;
   std::vector<float>   *phoR9;
   std::vector<float>   *phoEcalRHSumEtConeDR04;
   std::vector<float>   *phoHcalTowerSumEtBcConeDR04;
   std::vector<float>   *phoTrkSumPtSolidConeDR04;
   std::vector<float>   *phoTrkSumPtHollowConeDR04;
   std::vector<int>     *phoGenIdx;
   std::vector<float>   *phoGenDr;
   std::vector<float>   *phoSMaj;
   std::vector<float>   *phoSMin;
   std::vector<float>   *phoSAlp;
   std::vector<float>   *phoCovEtaEta;
   std::vector<float>   *phoCovEtaPhi;
   std::vector<float>   *phoCovPhiPhi;
   std::vector<float>   *genPt;
   std::vector<float>   *genEnergy;
   std::vector<float>   *genPhi;
   std::vector<float>   *genEta;
   std::vector<float>   *genPx;
   std::vector<float>   *genPy;
   std::vector<float>   *genPz;
   std::vector<int>     *genPdgId;
   std::vector<float>   *rhEnergy;
   std::vector<float>   *rhTime;
   std::vector<float>   *rhTOF;
   std::vector<unsigned int> *rhID;
   std::vector<bool>    *rhisOOT;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_vtxX;   //!
   TBranch        *b_vtxY;   //!
   TBranch        *b_vtxZ;   //!
   TBranch        *b_metSumEt;   //!
   TBranch        *b_metPx;   //!
   TBranch        *b_metPy;   //!
   TBranch        *b_metCSumEt;   //!
   TBranch        *b_metCPx;   //!
   TBranch        *b_metCPy;   //!
   TBranch        *b_jetHt;   //!
   TBranch        *b_jetE;   //!
   TBranch        *b_jetM;   //!
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
   TBranch        *b_jetParts;   //!
   TBranch        *b_jetDrRhIds;   //!
   TBranch        *b_jetGenImpactAngle;   //!
   TBranch        *b_jetGenTime;   //!
   TBranch        *b_jetGenPt;   //!
   TBranch        *b_jetGenEta;   //!
   TBranch        *b_jetGenEnergy;   //!
   TBranch        *b_jetGenDrMatch;   //!
   TBranch        *b_jetGenTimeLLP;   //!
   TBranch        *b_jetGenTOF;   //!
   TBranch        *b_phoIsOotPho;   //!
   TBranch        *b_phoExcluded;   //!
   TBranch        *b_phoSeedTOFTime;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEnergy;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPx;   //!
   TBranch        *b_phoPy;   //!
   TBranch        *b_phoPz;   //!
   TBranch        *b_phoRhIds;   //!
   TBranch        *b_phoIsPixelSeed;   //!
   TBranch        *b_phoIsEB;   //!
   TBranch        *b_phohadTowOverEM;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoEcalRHSumEtConeDR04;   //!
   TBranch        *b_phoHcalTowerSumEtBcConeDR04;   //!
   TBranch        *b_phoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_phoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_phoGenIdx;   //!
   TBranch        *b_phoGenDr;   //!
   TBranch        *b_phoSMaj;   //!
   TBranch        *b_phoSMin;   //!
   TBranch        *b_phoSAlp;   //!
   TBranch        *b_phoCovEtaEta;   //!
   TBranch        *b_phoCovEtaPhi;   //!
   TBranch        *b_phoCovPhiPhi;   //!
   TBranch        *b_genPt;   //!
   TBranch        *b_genEnergy;   //!
   TBranch        *b_genPhi;   //!
   TBranch        *b_genEta;   //!
   TBranch        *b_genPx;   //!
   TBranch        *b_genPy;   //!
   TBranch        *b_genPz;   //!
   TBranch        *b_genPdgId;   //!
   TBranch        *b_rhEnergy;   //!
   TBranch        *b_rhTime;   //!
   TBranch        *b_rhTOF;   //!
   TBranch        *b_rhID;   //!
   TBranch        *b_rhisOOT;   //!

   //root_rebase(TTree *tree=0);
   //virtual ~root_rebase();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef root_rebase_cxx
/*
root_rebase::root_rebase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gmsb_AODSIM_KUCMSNtuplizer_v1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gmsb_AODSIM_KUCMSNtuplizer_v1.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("gmsb_AODSIM_KUCMSNtuplizer_v1.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

root_rebase::~root_rebase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t root_rebase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t root_rebase::LoadTree(Long64_t entry)
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

void root_rebase::Init(TTree *tree)
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
   jetM = 0;
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
   jetParts = 0;
   jetDrRhIds = 0;
   jetGenImpactAngle = 0;
   jetGenTime = 0;
   jetGenPt = 0;
   jetGenEta = 0;
   jetGenEnergy = 0;
   jetGenDrMatch = 0;
   jetGenTimeLLP = 0;
   jetGenTOF = 0;
   phoIsOotPho = 0;
   phoExcluded = 0;
   phoSeedTOFTime = 0;
   phoPt = 0;
   phoEnergy = 0;
   phoPhi = 0;
   phoEta = 0;
   phoPx = 0;
   phoPy = 0;
   phoPz = 0;
   phoRhIds = 0;
   phoIsPixelSeed = 0;
   phoIsEB = 0;
   phohadTowOverEM = 0;
   phoSigmaIEtaIEta = 0;
   phoR9 = 0;
   phoEcalRHSumEtConeDR04 = 0;
   phoHcalTowerSumEtBcConeDR04 = 0;
   phoTrkSumPtSolidConeDR04 = 0;
   phoTrkSumPtHollowConeDR04 = 0;
   phoGenIdx = 0;
   phoGenDr = 0;
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
   rhEnergy = 0;
   rhTime = 0;
   rhTOF = 0;
   rhID = 0;
   rhisOOT = 0;
   // Set branch addresses and branch pointers
   //if (!tree) return;
   fChain = tree;
   //fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("vtxX", &vtxX, &b_vtxX);
   fChain->SetBranchAddress("vtxY", &vtxY, &b_vtxY);
   fChain->SetBranchAddress("vtxZ", &vtxZ, &b_vtxZ);
   fChain->SetBranchAddress("metSumEt", &metSumEt, &b_metSumEt);
   fChain->SetBranchAddress("metPx", &metPx, &b_metPx);
   fChain->SetBranchAddress("metPy", &metPy, &b_metPy);
   fChain->SetBranchAddress("metCSumEt", &metCSumEt, &b_metCSumEt);
   fChain->SetBranchAddress("metCPx", &metCPx, &b_metCPx);
   fChain->SetBranchAddress("metCPy", &metCPy, &b_metCPy);
   fChain->SetBranchAddress("jetHt", &jetHt, &b_jetHt);
   fChain->SetBranchAddress("jetE", &jetE, &b_jetE);
   fChain->SetBranchAddress("jetM", &jetM, &b_jetM);
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
   fChain->SetBranchAddress("jetParts", &jetParts, &b_jetParts);
   fChain->SetBranchAddress("jetDrRhIds", &jetDrRhIds, &b_jetDrRhIds);
   fChain->SetBranchAddress("jetGenImpactAngle", &jetGenImpactAngle, &b_jetGenImpactAngle);
   fChain->SetBranchAddress("jetGenTime", &jetGenTime, &b_jetGenTime);
   fChain->SetBranchAddress("jetGenPt", &jetGenPt, &b_jetGenPt);
   fChain->SetBranchAddress("jetGenEta", &jetGenEta, &b_jetGenEta);
   fChain->SetBranchAddress("jetGenEnergy", &jetGenEnergy, &b_jetGenEnergy);
   fChain->SetBranchAddress("jetGenDrMatch", &jetGenDrMatch, &b_jetGenDrMatch);
   fChain->SetBranchAddress("jetGenTimeLLP", &jetGenTimeLLP, &b_jetGenTimeLLP);
   fChain->SetBranchAddress("jetGenTOF", &jetGenTOF, &b_jetGenTOF);
   fChain->SetBranchAddress("phoIsOotPho", &phoIsOotPho, &b_phoIsOotPho);
   fChain->SetBranchAddress("phoExcluded", &phoExcluded, &b_phoExcluded);
   fChain->SetBranchAddress("phoSeedTOFTime", &phoSeedTOFTime, &b_phoSeedTOFTime);
   fChain->SetBranchAddress("phoPt", &phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEnergy", &phoEnergy, &b_phoEnergy);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPx", &phoPx, &b_phoPx);
   fChain->SetBranchAddress("phoPy", &phoPy, &b_phoPy);
   fChain->SetBranchAddress("phoPz", &phoPz, &b_phoPz);
   fChain->SetBranchAddress("phoRhIds", &phoRhIds, &b_phoRhIds);
   fChain->SetBranchAddress("phoIsPixelSeed", &phoIsPixelSeed, &b_phoIsPixelSeed);
   fChain->SetBranchAddress("phoIsEB", &phoIsEB, &b_phoIsEB);
   fChain->SetBranchAddress("phohadTowOverEM", &phohadTowOverEM, &b_phohadTowOverEM);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04, &b_phoEcalRHSumEtConeDR04);
   fChain->SetBranchAddress("phoHcalTowerSumEtBcConeDR04", &phoHcalTowerSumEtBcConeDR04, &b_phoHcalTowerSumEtBcConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04, &b_phoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04, &b_phoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("phoGenIdx", &phoGenIdx, &b_phoGenIdx);
   fChain->SetBranchAddress("phoGenDr", &phoGenDr, &b_phoGenDr);
   fChain->SetBranchAddress("phoSMaj", &phoSMaj, &b_phoSMaj);
   fChain->SetBranchAddress("phoSMin", &phoSMin, &b_phoSMin);
   fChain->SetBranchAddress("phoSAlp", &phoSAlp, &b_phoSAlp);
   fChain->SetBranchAddress("phoCovEtaEta", &phoCovEtaEta, &b_phoCovEtaEta);
   fChain->SetBranchAddress("phoCovEtaPhi", &phoCovEtaPhi, &b_phoCovEtaPhi);
   fChain->SetBranchAddress("phoCovPhiPhi", &phoCovPhiPhi, &b_phoCovPhiPhi);
   fChain->SetBranchAddress("genPt", &genPt, &b_genPt);
   fChain->SetBranchAddress("genEnergy", &genEnergy, &b_genEnergy);
   fChain->SetBranchAddress("genPhi", &genPhi, &b_genPhi);
   fChain->SetBranchAddress("genEta", &genEta, &b_genEta);
   fChain->SetBranchAddress("genPx", &genPx, &b_genPx);
   fChain->SetBranchAddress("genPy", &genPy, &b_genPy);
   fChain->SetBranchAddress("genPz", &genPz, &b_genPz);
   fChain->SetBranchAddress("genPdgId", &genPdgId, &b_genPdgId);
   fChain->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);
   fChain->SetBranchAddress("rhTime", &rhTime, &b_rhTime);
   fChain->SetBranchAddress("rhTOF", &rhTOF, &b_rhTOF);
   fChain->SetBranchAddress("rhID", &rhID, &b_rhID);
   fChain->SetBranchAddress("rhisOOT", &rhisOOT, &b_rhisOOT);
   //Notify();
}

/*
Bool_t root_rebase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void root_rebase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t root_rebase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/
#endif // #ifdef root_rebase_cxx
