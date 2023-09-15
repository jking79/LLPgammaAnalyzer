//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Aug  1 17:24:15 2023 by ROOT version 6.14/09
// from TTree kuSkimTree/output root file for kUCMSSkimmer
// found on file: kuntuple_gmsb_L150_AODSIM_Ntuple_v2_LLPgama_Skim_v6.root
//////////////////////////////////////////////////////////

#ifndef kuSkimTree_v1_h
#define kuSkimTree_v1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "string"
#include "vector"

class kuSkimTree_v1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   std::string     *DataSetKey;
   Float_t         selCMet;
   Float_t         selCMetPx;
   Float_t         selCMetPy;
   UInt_t          leadSelPho;
   UInt_t          nPhotons;
   UInt_t          nSelPhotons;
   std::vector<float>   *selPhoClstrRn;
   std::vector<float>   *selPhoEnergy;
   std::vector<float>   *selPhoEta;
   std::vector<float>   *selPhoGeoEgnVal;
   std::vector<float>   *selPhoGeoSMaj;
   std::vector<float>   *selPhoGeoSMin;
   std::vector<unsigned int> *selPhoNrh;
   std::vector<float>   *selPhoPhi;
   std::vector<float>   *selPhoPt;
   std::vector<int>     *selPhoQuality;
   std::vector<float>   *selPhoR9;
   std::vector<float>   *selPhoSMaj;
   std::vector<float>   *selPhoSMin;
   std::vector<float>   *selPhoSieie;
   std::vector<float>   *selPhoTime;
   UInt_t          subLeadSelPho;
   UInt_t          nJets;
   UInt_t          nSelJets;
   std::vector<float>   *selJetEnergy;
   std::vector<float>   *selJetEta;
   std::vector<float>   *selJetMass;
   std::vector<float>   *selJetPhi;
   std::vector<float>   *selJetPt;
   std::vector<int>     *selJetQuality;
   std::vector<float>   *selJetTime;
   std::vector<float>   *SCosA;
   std::vector<float>   *SMass;
   std::vector<float>   *X1aCosA;
   std::vector<float>   *X1aMass;
   std::vector<float>   *X1bCosA;
   std::vector<float>   *X1bMass;
   std::vector<float>   *X2aCosA;
   std::vector<float>   *X2aMass;
   std::vector<float>   *X2bCosA;
   std::vector<float>   *X2bMass;

   // List of branches
   TBranch        *b_DataSetKey;   //!
   TBranch        *b_selCMet;   //!
   TBranch        *b_selCMetPx;   //!
   TBranch        *b_selCMetPy;   //!
   TBranch        *b_leadSelPho;   //!
   TBranch        *b_nPhotons;   //!
   TBranch        *b_nSelPhotons;   //!
   TBranch        *b_selPhoClstrRn;   //!
   TBranch        *b_selPhoEnergy;   //!
   TBranch        *b_selPhoEta;   //!
   TBranch        *b_selPhoGeoEgnVal;   //!
   TBranch        *b_selPhoGeoSMaj;   //!
   TBranch        *b_selPhoGeoSMin;   //!
   TBranch        *b_selPhoNrh;   //!
   TBranch        *b_selPhoPhi;   //!
   TBranch        *b_selPhoPt;   //!
   TBranch        *b_selPhoQuality;   //!
   TBranch        *b_selPhoR9;   //!
   TBranch        *b_selPhoSMaj;   //!
   TBranch        *b_selPhoSMin;   //!
   TBranch        *b_selPhoSieie;   //!
   TBranch        *b_selPhoTime;   //!
   TBranch        *b_subLeadSelPho;   //!
   TBranch        *b_nJets;   //!
   TBranch        *b_nSelJets;   //!
   TBranch        *b_selJetEnergy;   //!
   TBranch        *b_selJetEta;   //!
   TBranch        *b_selJetMass;   //!
   TBranch        *b_selJetPhi;   //!
   TBranch        *b_selJetPt;   //!
   TBranch        *b_selJetQuality;   //!
   TBranch        *b_selJetTime;   //!
   TBranch        *b_SCosA;   //!
   TBranch        *b_SMass;   //!
   TBranch        *b_X1aCosA;   //!
   TBranch        *b_X1aMass;   //!
   TBranch        *b_X1bCosA;   //!
   TBranch        *b_X1bMass;   //!
   TBranch        *b_X2aCosA;   //!
   TBranch        *b_X2aMass;   //!
   TBranch        *b_X2bCosA;   //!
   TBranch        *b_X2bMass;   //!

   //kuSkimTree_v1(TTree *tree=0);
   //virtual ~kuSkimTree_v1();
   //virtual Int_t    Cut(Long64_t entry);
   //virtual Int_t    GetEntry(Long64_t entry);
   //virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     getBranches(Long64_t entry);
   //virtual void     Loop();
   //virtual Bool_t   Notify();
   //virtual void     Show(Long64_t entry = -1);
};

//#endif

//#ifdef kuSkimTree_v1_cxx

/*
kuSkimTree_v1::kuSkimTree_v1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("kuntuple_gmsb_L150_AODSIM_Ntuple_v2_LLPgama_Skim_v6.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("kuntuple_gmsb_L150_AODSIM_Ntuple_v2_LLPgama_Skim_v6.root");
      }
      f->GetObject("kuSkimTree",tree);

   }
   Init(tree);
}

kuSkimTree_v1::~kuSkimTree_v1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t kuSkimTree_v1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t kuSkimTree_v1::LoadTree(Long64_t entry)
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

void kuSkimTree_v1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   DataSetKey = 0;
   selPhoClstrRn = 0;
   selPhoEnergy = 0;
   selPhoEta = 0;
   selPhoGeoEgnVal = 0;
   selPhoGeoSMaj = 0;
   selPhoGeoSMin = 0;
   selPhoNrh = 0;
   selPhoPhi = 0;
   selPhoPt = 0;
   selPhoQuality = 0;
   selPhoR9 = 0;
   selPhoSMaj = 0;
   selPhoSMin = 0;
   selPhoSieie = 0;
   selPhoTime = 0;
   selJetEnergy = 0;
   selJetEta = 0;
   selJetMass = 0;
   selJetPhi = 0;
   selJetPt = 0;
   selJetQuality = 0;
   selJetTime = 0;
   SCosA = 0;
   SMass = 0;
   X1aCosA = 0;
   X1aMass = 0;
   X1bCosA = 0;
   X1bMass = 0;
   X2aCosA = 0;
   X2aMass = 0;
   X2bCosA = 0;
   X2bMass = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   //fCurrent = -1;
   //fChain->SetMakeClass(1);

   fChain->SetBranchAddress("DataSetKey", &DataSetKey, &b_DataSetKey);
   fChain->SetBranchAddress("selCMet", &selCMet, &b_selCMet);
   fChain->SetBranchAddress("selCMetPx", &selCMetPx, &b_selCMetPx);
   fChain->SetBranchAddress("selCMetPy", &selCMetPy, &b_selCMetPy);
   fChain->SetBranchAddress("leadSelPho", &leadSelPho, &b_leadSelPho);
   fChain->SetBranchAddress("nPhotons", &nPhotons, &b_nPhotons);
   fChain->SetBranchAddress("nSelPhotons", &nSelPhotons, &b_nSelPhotons);
   fChain->SetBranchAddress("selPhoClstrRn", &selPhoClstrRn, &b_selPhoClstrRn);
   fChain->SetBranchAddress("selPhoEnergy", &selPhoEnergy, &b_selPhoEnergy);
   fChain->SetBranchAddress("selPhoEta", &selPhoEta, &b_selPhoEta);
   fChain->SetBranchAddress("selPhoGeoEgnVal", &selPhoGeoEgnVal, &b_selPhoGeoEgnVal);
   fChain->SetBranchAddress("selPhoGeoSMaj", &selPhoGeoSMaj, &b_selPhoGeoSMaj);
   fChain->SetBranchAddress("selPhoGeoSMin", &selPhoGeoSMin, &b_selPhoGeoSMin);
   fChain->SetBranchAddress("selPhoNrh", &selPhoNrh, &b_selPhoNrh);
   fChain->SetBranchAddress("selPhoPhi", &selPhoPhi, &b_selPhoPhi);
   fChain->SetBranchAddress("selPhoPt", &selPhoPt, &b_selPhoPt);
   fChain->SetBranchAddress("selPhoQuality", &selPhoQuality, &b_selPhoQuality);
   fChain->SetBranchAddress("selPhoR9", &selPhoR9, &b_selPhoR9);
   fChain->SetBranchAddress("selPhoSMaj", &selPhoSMaj, &b_selPhoSMaj);
   fChain->SetBranchAddress("selPhoSMin", &selPhoSMin, &b_selPhoSMin);
   fChain->SetBranchAddress("selPhoSieie", &selPhoSieie, &b_selPhoSieie);
   fChain->SetBranchAddress("selPhoTime", &selPhoTime, &b_selPhoTime);
   fChain->SetBranchAddress("subLeadSelPho", &subLeadSelPho, &b_subLeadSelPho);
   fChain->SetBranchAddress("nJets", &nJets, &b_nJets);
   fChain->SetBranchAddress("nSelJets", &nSelJets, &b_nSelJets);
   fChain->SetBranchAddress("selJetEnergy", &selJetEnergy, &b_selJetEnergy);
   fChain->SetBranchAddress("selJetEta", &selJetEta, &b_selJetEta);
   fChain->SetBranchAddress("selJetMass", &selJetMass, &b_selJetMass);
   fChain->SetBranchAddress("selJetPhi", &selJetPhi, &b_selJetPhi);
   fChain->SetBranchAddress("selJetPt", &selJetPt, &b_selJetPt);
   fChain->SetBranchAddress("selJetQuality", &selJetQuality, &b_selJetQuality);
   fChain->SetBranchAddress("selJetTime", &selJetTime, &b_selJetTime);
   fChain->SetBranchAddress("SCosA", &SCosA, &b_SCosA);
   fChain->SetBranchAddress("SMass", &SMass, &b_SMass);
   fChain->SetBranchAddress("X1aCosA", &X1aCosA, &b_X1aCosA);
   fChain->SetBranchAddress("X1aMass", &X1aMass, &b_X1aMass);
   fChain->SetBranchAddress("X1bCosA", &X1bCosA, &b_X1bCosA);
   fChain->SetBranchAddress("X1bMass", &X1bMass, &b_X1bMass);
   fChain->SetBranchAddress("X2aCosA", &X2aCosA, &b_X2aCosA);
   fChain->SetBranchAddress("X2aMass", &X2aMass, &b_X2aMass);
   fChain->SetBranchAddress("X2bCosA", &X2bCosA, &b_X2bCosA);
   fChain->SetBranchAddress("X2bMass", &X2bMass, &b_X2bMass);
   //Notify();
}

void kuSkimTree_v1::getBranches(Long64_t entry){

   b_DataSetKey->GetEntry(entry);    //!
   b_selCMet->GetEntry(entry);    //!
   b_selCMetPx->GetEntry(entry);    //!
   b_selCMetPy->GetEntry(entry);    //!
   b_leadSelPho->GetEntry(entry);    //!
   b_nPhotons->GetEntry(entry);    //!
   b_nSelPhotons->GetEntry(entry);    //!
   b_selPhoClstrRn->GetEntry(entry);    //!
   b_selPhoEnergy->GetEntry(entry);    //!
   b_selPhoEta->GetEntry(entry);    //!
   b_selPhoGeoEgnVal->GetEntry(entry);    //!
   b_selPhoGeoSMaj->GetEntry(entry);    //!
   b_selPhoGeoSMin->GetEntry(entry);    //!
   b_selPhoNrh->GetEntry(entry);    //!
   b_selPhoPhi->GetEntry(entry);    //!
   b_selPhoPt->GetEntry(entry);    //!
   b_selPhoQuality->GetEntry(entry);    //!
   b_selPhoR9->GetEntry(entry);    //!
   b_selPhoSMaj->GetEntry(entry);    //!
   b_selPhoSMin->GetEntry(entry);    //!
   b_selPhoSieie->GetEntry(entry);    //!
   b_selPhoTime->GetEntry(entry);    //!
   b_subLeadSelPho->GetEntry(entry);    //!
   b_nJets->GetEntry(entry);    //!
   b_nSelJets->GetEntry(entry);    //!
   b_selJetEnergy->GetEntry(entry);    //!
   b_selJetEta->GetEntry(entry);    //!
   b_selJetMass->GetEntry(entry);    //!
   b_selJetPhi->GetEntry(entry);    //!
   b_selJetPt->GetEntry(entry);    //!
   b_selJetQuality->GetEntry(entry);    //!
   b_selJetTime->GetEntry(entry);    //!
   b_SCosA->GetEntry(entry);    //!
   b_SMass->GetEntry(entry);    //!
   b_X1aCosA->GetEntry(entry);    //!
   b_X1aMass->GetEntry(entry);    //!
   b_X1bCosA->GetEntry(entry);    //!
   b_X1bMass->GetEntry(entry);    //!
   b_X2aCosA->GetEntry(entry);    //!
   b_X2aMass->GetEntry(entry);    //!
   b_X2bCosA->GetEntry(entry);    //!
   b_X2bMass->GetEntry(entry);    //!

}

/*
Bool_t kuSkimTree_v1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void kuSkimTree_v1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t kuSkimTree_v1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
*/

#endif // #ifdef kuSkimTree_v1_cxx
