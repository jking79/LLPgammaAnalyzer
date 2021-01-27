//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 10 20:44:28 2019 by ROOT version 6.14/04
// from TTree Runs/Runs
// found on file: prod2016MC_NANO_1-10.root
//////////////////////////////////////////////////////////

#ifndef SUSYNANORuneBase_h
#define SUSYNANORuneBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class SUSYNANORuneBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          run;
   Long64_t        genEventCount;
   Double_t        genEventSumw;
   Double_t        genEventSumw2;
   UInt_t          nLHEScaleSumw;
   Double_t        LHEScaleSumw[9];   //[nLHEScaleSumw]
   UInt_t          nLHEPdfSumw;
   Double_t        LHEPdfSumw[102];   //[nLHEPdfSumw]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_genEventCount;   //!
   TBranch        *b_genEventSumw;   //!
   TBranch        *b_genEventSumw2;   //!
   TBranch        *b_nLHEScaleSumw;   //!
   TBranch        *b_LHEScaleSumw;   //!
   TBranch        *b_nLHEPdfSumw;   //!
   TBranch        *b_LHEPdfSumw;   //!

   SUSYNANORuneBase(TTree *tree=0);
   virtual ~SUSYNANORuneBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

inline SUSYNANORuneBase::SUSYNANORuneBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("prod2016MC_NANO_1-10.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("prod2016MC_NANO_1-10.root");
      }
      f->GetObject("Runs",tree);

   }
   Init(tree);
}

inline SUSYNANORuneBase::~SUSYNANORuneBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t SUSYNANORuneBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t SUSYNANORuneBase::LoadTree(Long64_t entry)
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

inline void SUSYNANORuneBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("genEventCount", &genEventCount, &b_genEventCount);
   fChain->SetBranchAddress("genEventSumw", &genEventSumw, &b_genEventSumw);
   fChain->SetBranchAddress("genEventSumw2", &genEventSumw2, &b_genEventSumw2);
   fChain->SetBranchAddress("nLHEScaleSumw", &nLHEScaleSumw, &b_nLHEScaleSumw);
   fChain->SetBranchAddress("LHEScaleSumw", LHEScaleSumw, &b_LHEScaleSumw);
   fChain->SetBranchAddress("nLHEPdfSumw", &nLHEPdfSumw, &b_nLHEPdfSumw);
   fChain->SetBranchAddress("LHEPdfSumw", LHEPdfSumw, &b_LHEPdfSumw);
   Notify();
}

inline Bool_t SUSYNANORuneBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void SUSYNANORuneBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t SUSYNANORuneBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
