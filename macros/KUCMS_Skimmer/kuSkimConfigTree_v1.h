//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  1 17:38:55 2023 by ROOT version 6.14/09
// from TTree kuSkimConfigTree/config root file for kUCMSSkimmer
// found on file: kuntuple_gmsb_L150_AODSIM_Ntuple_v2_LLPgama_Skim_v6.root
//////////////////////////////////////////////////////////

#ifndef kuSkimConfigTree_v1_h
#define kuSkimConfigTree_v1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"

class kuSkimConfigTree_v1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          nEvents;
   UInt_t          nSelectedEvents;
   string          *sKey;
   Float_t         sCrossSection;
   Float_t         sGMSBGravMass;
   Float_t         sGMSBChi1Mass;
   Float_t         sMCWgt;
   Int_t           sMCType;

   // List of branches
   TBranch        *b_nEvents;   //!
   TBranch        *b_nSelectedEvents;   //!
   TBranch        *b_sKey;   //!
   TBranch        *b_sCrossSection;   //!
   TBranch        *b_sGMSBGravMass;   //!
   TBranch        *b_sGMSBChi1Mass;   //!
   TBranch        *b_sMCWgt;   //!
   TBranch        *b_sMCType;   //!

   kuSkimConfigTree_v1(TTree *tree=0);
   virtual ~kuSkimConfigTree_v1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef kuSkimConfigTree_v1_cxx
kuSkimConfigTree_v1::kuSkimConfigTree_v1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("kuntuple_gmsb_L150_AODSIM_Ntuple_v2_LLPgama_Skim_v6.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("kuntuple_gmsb_L150_AODSIM_Ntuple_v2_LLPgama_Skim_v6.root");
      }
      f->GetObject("kuSkimConfigTree",tree);

   }
   Init(tree);
}

kuSkimConfigTree_v1::~kuSkimConfigTree_v1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kuSkimConfigTree_v1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kuSkimConfigTree_v1::LoadTree(Long64_t entry)
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

void kuSkimConfigTree_v1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   sKey = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
   fChain->SetBranchAddress("nSelectedEvents", &nSelectedEvents, &b_nSelectedEvents);
   fChain->SetBranchAddress("sKey", &sKey, &b_sKey);
   fChain->SetBranchAddress("sCrossSection", &sCrossSection, &b_sCrossSection);
   fChain->SetBranchAddress("sGMSBGravMass", &sGMSBGravMass, &b_sGMSBGravMass);
   fChain->SetBranchAddress("sGMSBChi1Mass", &sGMSBChi1Mass, &b_sGMSBChi1Mass);
   fChain->SetBranchAddress("sMCWgt", &sMCWgt, &b_sMCWgt);
   fChain->SetBranchAddress("sMCType", &sMCType, &b_sMCType);
   Notify();
}

Bool_t kuSkimConfigTree_v1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kuSkimConfigTree_v1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kuSkimConfigTree_v1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef kuSkimConfigTree_v1_cxx