//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar  6 12:37:57 2023 by ROOT version 6.24/07
// from TTree llpgtree/llpgtree
// found on file: ku_KUCC_tt_R2018_126_gammares_v5.root
//////////////////////////////////////////////////////////

//#ifndef egammares_hist_base_h
//#define egammares_hist_base_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"

class egammares_hist_base {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   //Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   std::vector<unsigned int> *rhCaliID;
   std::vector<float>   *rhCaliEnergy;
   std::vector<float>   *rhCaliRtTime;
   std::vector<float>   *rhCaliCCTime;

   std::vector<unsigned int> *resRhID;
   std::vector<float>   *resAmp;
   std::vector<float>   *resE;
   std::vector<float>   *resRtTime;
   std::vector<float>   *resCCTime;
   std::vector<float>   *resTOF;

   std::vector<unsigned int> *rhID;
   std::vector<float>   *rhRtTime;
   std::vector<float>   *rhCCTime;
   std::vector<float>   *rhTOF;
   std::vector<float>   *rhEnergy;
   std::vector<float>   *rhAmp;
   std::vector<bool>    *rhRtisOOT;
   std::vector<bool>    *rhCCisOOT;
   std::vector<bool>    *rhisWeird;
   std::vector<bool>    *rhisDiWeird;
   std::vector<float>   *rhSwCross;
   std::vector<bool>    *rhisGS6;
   std::vector<bool>    *rhisGS1;
   std::vector<float>   *rhadcToGeV;
   std::vector<float>   *rhpedrms12;

   std::vector<float>   *phoEnergy;
   std::vector<std::vector<unsigned int> > *phoRhIds;
   std::vector<float>   *phoPt;
   std::vector<float>   *phoEta;
   std::vector<float>   *phoPhi;
   std::vector<float>   *phoHadOverEM;
   std::vector<float>   *phoSigmaIEtaIEta;
   std::vector<float>   *phoCov2IEtaIEta;
   std::vector<float>   *phoCov2IEtaIPhi;
   std::vector<float>   *phoCov2IPhiIPhi;
   std::vector<float>   *phoEcalRHSumEtConeDR04;
   std::vector<float>   *phoHcalTwrSumEtConeDR04;
   std::vector<float>   *phoTrkSumPtSolidConeDR04;
   std::vector<float>   *phoTrkSumPtHollowConeDR04;
   std::vector<float>   *phoR9;
   std::vector<int>     *phoSelType;

   Float_t         phoDiMass;
   Float_t         phoDiAngle;
   Float_t         phoDiDr;
   Float_t         phoDiPhi;
   Float_t         phoDiEta;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_rhCaliID;   //!
   TBranch        *b_rhCaliEnergy;   //!
   TBranch        *b_rhCaliRtTime;   //!
   TBranch        *b_rhCaliCCTime;   //!
   TBranch        *b_resRhID;   //!
   TBranch        *b_resAmp;   //!
   TBranch        *b_resE;   //!
   TBranch        *b_resRtTime;   //!
   TBranch        *b_resCCTime;   //!
   TBranch        *b_resTOF;   //!
   TBranch        *b_rhID;   //!
   TBranch        *b_rhRtTime;   //!
   TBranch        *b_rhCCTime;   //!
   TBranch        *b_rhTOF;   //!
   TBranch        *b_rhEnergy;   //!
   TBranch        *b_rhAmp;   //!
   TBranch        *b_rhRtisOOT;   //!
   TBranch        *b_rhCCisOOT;   //!
   TBranch        *b_rhisWeird;   //!
   TBranch        *b_rhisDiWeird;   //!
   TBranch        *b_rhSwCross;   //!
   TBranch        *b_rhisGS6;   //!
   TBranch        *b_rhisGS1;   //!
   TBranch        *b_rhadcToGeV;   //!
   TBranch        *b_rhpedrms12;   //!
   TBranch        *b_phoEnergy;   //!
   TBranch        *b_phoRhIds;   //!
   TBranch        *b_phoPt;   //!
   TBranch        *b_phoEta;   //!
   TBranch        *b_phoPhi;   //!
   TBranch        *b_phoHadOverEM;   //!
   TBranch        *b_phoSigmaIEtaIEta;   //!
   TBranch        *b_phoCov2IEtaIEta;   //!
   TBranch        *b_phoCov2IEtaIPhi;   //!
   TBranch        *b_phoCov2IPhiIPhi;   //!
   TBranch        *b_phoEcalRHSumEtConeDR04;   //!
   TBranch        *b_phoHcalTwrSumEtConeDR04;   //!
   TBranch        *b_phoTrkSumPtSolidConeDR04;   //!
   TBranch        *b_phoTrkSumPtHollowConeDR04;   //!
   TBranch        *b_phoR9;   //!
   TBranch        *b_phoSelType;   //!
   TBranch        *b_phoDiMass;   //!
   TBranch        *b_phoDiAngle;   //!
   TBranch        *b_phoDiDr;   //!
   TBranch        *b_phoDiPhi;   //!
   TBranch        *b_phoDiEta;   //!

   //egammares_hist_base(TTree *tree=0);
   //virtual ~egammares_hist_base();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef egammares_hist_base_cxx
/*
egammares_hist_base::egammares_hist_base(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ku_KUCC_tt_R2018_126_gammares_v5.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ku_KUCC_tt_R2018_126_gammares_v5.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("ku_KUCC_tt_R2018_126_gammares_v5.root:/tree");
      dir->GetObject("llpgtree",tree);

   }
   Init(tree);
}

egammares_hist_base::~egammares_hist_base()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t egammares_hist_base::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t egammares_hist_base::LoadTree(Long64_t entry)
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
void egammares_hist_base::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   rhCaliID = 0;
   rhCaliEnergy = 0;
   rhCaliRtTime = 0;
   rhCaliCCTime = 0;
   resRhID = 0;
   resAmp = 0;
   resE = 0;
   resRtTime = 0;
   resCCTime = 0;
   resTOF = 0;
   rhID = 0;
   rhRtTime = 0;
   rhCCTime = 0;
   rhTOF = 0;
   rhEnergy = 0;
   rhAmp = 0;
   rhRtisOOT = 0;
   rhCCisOOT = 0;
   rhisWeird = 0;
   rhisDiWeird = 0;
   rhSwCross = 0;
   rhisGS6 = 0;
   rhisGS1 = 0;
   rhadcToGeV = 0;
   rhpedrms12 = 0;
   phoEnergy = 0;
   phoRhIds = 0;
   phoPt = 0;
   phoEta = 0;
   phoPhi = 0;
   phoHadOverEM = 0;
   phoSigmaIEtaIEta = 0;
   phoCov2IEtaIEta = 0;
   phoCov2IEtaIPhi = 0;
   phoCov2IPhiIPhi = 0;
   phoEcalRHSumEtConeDR04 = 0;
   phoHcalTwrSumEtConeDR04 = 0;
   phoTrkSumPtSolidConeDR04 = 0;
   phoTrkSumPtHollowConeDR04 = 0;
   phoR9 = 0;
   phoSelType = 0;
   // Set branch addresses and branch pointers
   //if (!tree) return;
   fChain = tree;
   //fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("rhCaliID", &rhCaliID, &b_rhCaliID);
   fChain->SetBranchAddress("rhCaliEnergy", &rhCaliEnergy, &b_rhCaliEnergy);
   fChain->SetBranchAddress("rhCaliRtTime", &rhCaliRtTime, &b_rhCaliRtTime);
   fChain->SetBranchAddress("rhCaliCCTime", &rhCaliCCTime, &b_rhCaliCCTime);
   fChain->SetBranchAddress("resRhID", &resRhID, &b_resRhID);
   fChain->SetBranchAddress("resAmp", &resAmp, &b_resAmp);
   fChain->SetBranchAddress("resE", &resE, &b_resE);
   fChain->SetBranchAddress("resRtTime", &resRtTime, &b_resRtTime);
   fChain->SetBranchAddress("resCCTime", &resCCTime, &b_resCCTime);
   fChain->SetBranchAddress("resTOF", &resTOF, &b_resTOF);
   fChain->SetBranchAddress("rhID", &rhID, &b_rhID);
   fChain->SetBranchAddress("rhRtTime", &rhRtTime, &b_rhRtTime);
   fChain->SetBranchAddress("rhCCTime", &rhCCTime, &b_rhCCTime);
   fChain->SetBranchAddress("rhTOF", &rhTOF, &b_rhTOF);
   fChain->SetBranchAddress("rhEnergy", &rhEnergy, &b_rhEnergy);
   fChain->SetBranchAddress("rhAmp", &rhAmp, &b_rhAmp);
   fChain->SetBranchAddress("rhRtisOOT", &rhRtisOOT, &b_rhRtisOOT);
   fChain->SetBranchAddress("rhCCisOOT", &rhCCisOOT, &b_rhCCisOOT);
   fChain->SetBranchAddress("rhisWeird", &rhisWeird, &b_rhisWeird);
   fChain->SetBranchAddress("rhisDiWeird", &rhisDiWeird, &b_rhisDiWeird);
   fChain->SetBranchAddress("rhSwCross", &rhSwCross, &b_rhSwCross);
   fChain->SetBranchAddress("rhisGS6", &rhisGS6, &b_rhisGS6);
   fChain->SetBranchAddress("rhisGS1", &rhisGS1, &b_rhisGS1);
   fChain->SetBranchAddress("rhadcToGeV", &rhadcToGeV, &b_rhadcToGeV);
   fChain->SetBranchAddress("rhpedrms12", &rhpedrms12, &b_rhpedrms12);
   fChain->SetBranchAddress("phoEnergy", &phoEnergy, &b_phoEnergy);
   fChain->SetBranchAddress("phoRhIds", &phoRhIds, &b_phoRhIds);
   fChain->SetBranchAddress("phoPt", &phoPt, &b_phoPt);
   fChain->SetBranchAddress("phoEta", &phoEta, &b_phoEta);
   fChain->SetBranchAddress("phoPhi", &phoPhi, &b_phoPhi);
   fChain->SetBranchAddress("phoHadOverEM", &phoHadOverEM, &b_phoHadOverEM);
   fChain->SetBranchAddress("phoSigmaIEtaIEta", &phoSigmaIEtaIEta, &b_phoSigmaIEtaIEta);
   fChain->SetBranchAddress("phoCov2IEtaIEta", &phoCov2IEtaIEta, &b_phoCov2IEtaIEta);
   fChain->SetBranchAddress("phoCov2IEtaIPhi", &phoCov2IEtaIPhi, &b_phoCov2IEtaIPhi);
   fChain->SetBranchAddress("phoCov2IPhiIPhi", &phoCov2IPhiIPhi, &b_phoCov2IPhiIPhi);
   fChain->SetBranchAddress("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04, &b_phoEcalRHSumEtConeDR04);
   fChain->SetBranchAddress("phoHcalTwrSumEtConeDR04", &phoHcalTwrSumEtConeDR04, &b_phoHcalTwrSumEtConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04, &b_phoTrkSumPtSolidConeDR04);
   fChain->SetBranchAddress("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04, &b_phoTrkSumPtHollowConeDR04);
   fChain->SetBranchAddress("phoR9", &phoR9, &b_phoR9);
   fChain->SetBranchAddress("phoSelType", &phoSelType, &b_phoSelType);
   fChain->SetBranchAddress("phoDiMass", &phoDiMass, &b_phoDiMass);
   fChain->SetBranchAddress("phoDiAngle", &phoDiAngle, &b_phoDiAngle);
   fChain->SetBranchAddress("phoDiDr", &phoDiDr, &b_phoDiDr);
   fChain->SetBranchAddress("phoDiPhi", &phoDiPhi, &b_phoDiPhi);
   fChain->SetBranchAddress("phoDiEta", &phoDiEta, &b_phoDiEta);
   //Notify();
}

/*
Bool_t egammares_hist_base::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void egammares_hist_base::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t egammares_hist_base::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/
//#endif // #ifdef egammares_hist_base_cxx
