//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Feb 10 02:08:57 2019 by ROOT version 6.14/04
// from TTree AUX/AUX
// found on file: SIG/TChiWZ.root
//////////////////////////////////////////////////////////

#ifndef StopNtupleTree_h
#define StopNtupleTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TLeafElement.h>
#include <TLorentzVector.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

using std::vector;
using std::string;

class StopNtupleTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

  bool m_USEFLOAT;
  // "special" declarations
  Float_t         met_f;
  Float_t         metphi_f;
  Float_t         calomet_f;
  Float_t         calometphi_f;
  Float_t         stored_weight_f;
  Float_t         evtWeight_f;

  double          met_d;
  double          metphi_d;
  double          calomet_d;
  double          calometphi_d;
  double          stored_weight_d;
  double          evtWeight_d;
  
   // Declaration of leaf types
   UInt_t          run;
   UInt_t          lumi;
   ULong64_t       event;
   Float_t         genmet;
   Float_t         genmetphi;
   Float_t         mht;
   Float_t         mhtphi;
   Float_t         ht;
   Float_t         dPhi0_CUT;
   Float_t         dPhi1_CUT;
   Float_t         dPhi2_CUT;
   Float_t         tru_npv;
   Float_t         avg_npv;
   Int_t           METFilters;
   Int_t           CSCTightHaloFilter;
   Int_t           globalSuperTightHalo2016Filter;
   Int_t           goodVerticesFilter;
   Int_t           ecalBadCalibFilter;
   Int_t           HBHENoiseIsoFilter;
   Int_t           EcalDeadCellTriggerPrimitiveFilter;
   Int_t           noBadMuonsFilter;
   Int_t           badMuonsFilter;
   Int_t           duplicateMuonsFilter;
   Int_t           nMuons_CUT;
   Int_t           nMuons;
   Int_t           nElectrons_CUT;
   Int_t           nElectrons;
   Int_t           NJetsISR;
   Int_t           loose_nIsoTrks;
   Int_t           nIsoTrks_CUT;
   Int_t           nJets_CUT;
   Int_t           vtxSize;
   Int_t           npv;
   Int_t           nm1;
   Int_t           n0;
   Int_t           np1;
   UInt_t          looseJetID;
   UInt_t          tightJetID;
   UInt_t          tightlepvetoJetID;
   UInt_t          looseJetID_NoLep;
   UInt_t          tightJetID_NoLep;
   UInt_t          tightlepvetoJetID_NoLep;
   UInt_t          BadChargedCandidateFilter;
   UInt_t          BadPFMuonFilter;
   UInt_t          HBHENoiseFilter;
   UInt_t          HBHEIsoNoiseFilter;
   vector<float>   *muonsCharge;
   vector<float>   *muonsMtw;
   vector<float>   *muonsRelIso;
   vector<float>   *muonsMiniIso;
   vector<float>   *muonspfActivity;
   vector<float>   *pfGammaIso;
   vector<float>   *isEB;
   vector<float>   *genMatched;
   vector<float>   *hadTowOverEM;
   vector<float>   *sigmaIetaIeta;
   vector<float>   *pfChargedIso;
   vector<float>   *pfNeutralIso;
   vector<float>   *pfChargedIsoRhoCorr;
   vector<float>   *pfNeutralIsoRhoCorr;
   vector<float>   *pfGammaIsoRhoCorr;
   vector<float>   *hasPixelSeed;
   vector<float>   *passElectronVeto;
   vector<float>   *photonPt;
   vector<float>   *photonEta;
   vector<float>   *photonPhi;
   vector<float>   *svPT;
   vector<float>   *svETA;
   vector<float>   *svPhi;
   vector<float>   *svMass;
   vector<float>   *svNTracks;
   vector<float>   *svChi2;
   vector<float>   *svNDF;
   vector<float>   *svDXY;
   vector<float>   *svDXYerr;
   vector<float>   *svD3D;
   vector<float>   *svD3Derr;
   vector<float>   *svCosThetaSVPS;
   vector<float>   *elesCharge;
   vector<float>   *elesMtw;
   vector<float>   *elesRelIso;
   vector<float>   *elesMiniIso;
   vector<float>   *elespfActivity;
   vector<float>   *recoJetsJecUnc;
   vector<float>   *recoJetsJecScaleRawToFull;
   vector<float>   *qgLikelihood;
   vector<float>   *qgPtD;
   vector<float>   *qgAxis2;
   vector<float>   *recoJetschargedHadronEnergyFraction;
   vector<float>   *recoJetschargedEmEnergyFraction;
   vector<float>   *recoJetsneutralEmEnergyFraction;
   vector<float>   *recoJetsmuonEnergyFraction;
   vector<float>   *recoJetsBtag_0;
   vector<float>   *recoJetsCharge_0;
   vector<float>   *DeepCSVb;
   vector<float>   *DeepCSVc;
   vector<float>   *DeepCSVl;
   vector<float>   *DeepCSVbb;
   vector<float>   *DeepCSVcc;
   vector<float>   *DeepCSVbN;
   vector<float>   *DeepCSVcN;
   vector<float>   *DeepCSVlN;
   vector<float>   *DeepCSVbbN;
   vector<float>   *DeepCSVccN;
   vector<float>   *DeepCSVbP;
   vector<float>   *DeepCSVcP;
   vector<float>   *DeepCSVlP;
   vector<float>   *DeepCSVbbP;
   vector<float>   *DeepCSVccP;
   vector<float>   *puppitau1;
   vector<float>   *puppitau2;
   vector<float>   *puppitau3;
   vector<float>   *puppisoftDropMass;
   vector<float>   *puppiSubJetsBdisc;
   vector<float>   *recoJetsJecUncLepCleaned;
   vector<float>   *prodJetsNoLep_qgLikelihood;
   vector<float>   *prodJetsNoLep_qgPtD;
   vector<float>   *prodJetsNoLep_qgAxis2;
   vector<float>   *recoJetschargedHadronEnergyFractionLepCleaned;
   vector<float>   *recoJetsneutralEmEnergyFractionLepCleaned;
   vector<float>   *recoJetschargedEmEnergyFractionLepCleaned;
   vector<float>   *recoJetsmuonEnergyFractionLepCleaned;
   vector<float>   *recoJetsBtag_0_LepCleaned;
   vector<float>   *recoJetsCharge_0_LepCleaned;
   vector<float>   *recoJetsJecScaleRawToFull_LepCleaned;
   vector<float>   *prodJetsNoLep_puppisoftDropMass;
   vector<float>   *prodJetsNoLep_puppitau1;
   vector<float>   *prodJetsNoLep_puppitau2;
   vector<float>   *prodJetsNoLep_puppitau3;
   vector<float>   *prodJetsNoLep_puppiSubJetsBdisc;
   vector<float>   *W_emu_pfActivityVec;
   vector<float>   *W_tau_emu_pfActivityVec;
   vector<float>   *W_tau_prongs_pfActivityVec;
   vector<float>   *trksForIsoVetocharge;
   vector<float>   *trksForIsoVetodz;
   vector<float>   *trksForIsoVetoiso;
   vector<float>   *trksForIsoVetopfActivity;
   vector<float>   *loose_isoTrks_charge;
   vector<float>   *loose_isoTrks_dz;
   vector<float>   *loose_isoTrks_iso;
   vector<float>   *loose_isoTrks_mtw;
   vector<float>   *loose_isoTrks_pfActivity;
   vector<float>   *metMagUp;
   vector<float>   *metMagDown;
   vector<float>   *metPhiUp;
   vector<float>   *metPhiDown;
   vector<int>     *PassTrigger;
   vector<int>     *TriggerPrescales;
   vector<int>     *muonsFlagMedium;
   vector<int>     *muonsFlagTight;
   vector<int>     *elesFlagMedium;
   vector<int>     *elesFlagVeto;
   vector<int>     *recoJetsFlavor;
   vector<int>     *qgMult;
   vector<int>     *muMatchedJetIdx;
   vector<int>     *eleMatchedJetIdx;
   vector<int>     *looseisoTrksMatchedJetIdx;
   vector<int>     *trksForIsoVetoMatchedJetIdx;
   vector<int>     *prodJetsNoLep_qgMult;
   vector<int>     *genDecayIdxVec;
   vector<int>     *genDecayPdgIdVec;
   vector<int>     *genDecayMomIdxVec;
   vector<int>     *genDecayMomRefVec;
   vector<int>     *W_emuVec;
   vector<int>     *W_tauVec;
   vector<int>     *W_tau_emuVec;
   vector<int>     *W_tau_prongsVec;
   vector<int>     *W_tau_nuVec;
   vector<int>     *selPDGid;
   vector<int>     *trksForIsoVetopdgId;
   vector<int>     *trksForIsoVetoidx;
   vector<int>     *loose_isoTrks_pdgId;
   vector<int>     *loose_isoTrks_idx;
   vector<int>     *forVetoIsoTrksidx;
   vector<unsigned int> *loosePhotonID;
   vector<unsigned int> *mediumPhotonID;
   vector<unsigned int> *tightPhotonID;
   vector<unsigned int> *nonPrompt;
   vector<unsigned int> *elesisEB;
   vector<unsigned int> *vetoElectronID;
   vector<unsigned int> *looseElectronID;
   vector<unsigned int> *mediumElectronID;
   vector<unsigned int> *tightElectronID;
   vector<string>  *TriggerNames;
   vector<string>  *genDecayStrVec;
   vector<TLorentzVector> *muonsLVec;
   vector<TLorentzVector> *gammaLVec;
   vector<TLorentzVector> *gammaLVecGen;
   vector<TLorentzVector> *genPartonLVec;
   vector<TLorentzVector> *svSoftLVec;
   vector<TLorentzVector> *svLVec;
   vector<TLorentzVector> *elesLVec;
   vector<TLorentzVector> *jetsLVec;
   vector<TLorentzVector> *puppiJetsLVec;
   vector<TLorentzVector> *puppiSubJetsLVec;
   vector<TLorentzVector> *jetsLVecLepCleaned;
   vector<TLorentzVector> *prodJetsNoLep_puppiJetsLVec;
   vector<TLorentzVector> *prodJetsNoLep_puppiSubJetsLVec;
   vector<TLorentzVector> *genDecayLVec;
   vector<TLorentzVector> *selGenParticle;
   vector<TLorentzVector> *genjetsLVec;
   vector<TLorentzVector> *loose_isoTrksLVec;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_event;   //!
   TBranch        *b_genmet;   //!
   TBranch        *b_genmetphi;   //!
   TBranch        *b_mht;   //!
   TBranch        *b_mhtphi;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_calomet;   //!
   TBranch        *b_calometphi;   //!
   TBranch        *b_dPhi0_CUT;   //!
   TBranch        *b_dPhi1_CUT;   //!
   TBranch        *b_dPhi2_CUT;   //!
   TBranch        *b_tru_npv;   //!
   TBranch        *b_avg_npv;   //!
   TBranch        *b_stored_weight;   //!
   TBranch        *b_evtWeight;   //!
   TBranch        *b_METFilters;   //!
   TBranch        *b_CSCTightHaloFilter;   //!
   TBranch        *b_globalSuperTightHalo2016Filter;   //!
   TBranch        *b_goodVerticesFilter;   //!
   TBranch        *b_ecalBadCalibFilter;   //!
   TBranch        *b_HBHENoiseIsoFilter;   //!
   TBranch        *b_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_noBadMuonsFilter;   //!
   TBranch        *b_badMuonsFilter;   //!
   TBranch        *b_duplicateMuonsFilter;   //!
   TBranch        *b_nMuons_CUT;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_nElectrons_CUT;   //!
   TBranch        *b_nElectrons;   //!
   TBranch        *b_NJetsISR;   //!
   TBranch        *b_loose_nIsoTrks;   //!
   TBranch        *b_nIsoTrks_CUT;   //!
   TBranch        *b_nJets_CUT;   //!
   TBranch        *b_vtxSize;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_nm1;   //!
   TBranch        *b_n0;   //!
   TBranch        *b_np1;   //!
   TBranch        *b_looseJetID;   //!
   TBranch        *b_tightJetID;   //!
   TBranch        *b_tightlepvetoJetID;   //!
   TBranch        *b_looseJetID_NoLep;   //!
   TBranch        *b_tightJetID_NoLep;   //!
   TBranch        *b_tightlepvetoJetID_NoLep;   //!
   TBranch        *b_BadChargedCandidateFilter;   //!
   TBranch        *b_BadPFMuonFilter;   //!
   TBranch        *b_HBHENoiseFilter;   //!
   TBranch        *b_HBHEIsoNoiseFilter;   //!
   TBranch        *b_muonsCharge;   //!
   TBranch        *b_muonsMtw;   //!
   TBranch        *b_muonsRelIso;   //!
   TBranch        *b_muonsMiniIso;   //!
   TBranch        *b_muonspfActivity;   //!
   TBranch        *b_pfGammaIso;   //!
   TBranch        *b_isEB;   //!
   TBranch        *b_genMatched;   //!
   TBranch        *b_hadTowOverEM;   //!
   TBranch        *b_sigmaIetaIeta;   //!
   TBranch        *b_pfChargedIso;   //!
   TBranch        *b_pfNeutralIso;   //!
   TBranch        *b_pfChargedIsoRhoCorr;   //!
   TBranch        *b_pfNeutralIsoRhoCorr;   //!
   TBranch        *b_pfGammaIsoRhoCorr;   //!
   TBranch        *b_hasPixelSeed;   //!
   TBranch        *b_passElectronVeto;   //!
   TBranch        *b_photonPt;   //!
   TBranch        *b_photonEta;   //!
   TBranch        *b_photonPhi;   //!
   TBranch        *b_svPT;   //!
   TBranch        *b_svETA;   //!
   TBranch        *b_svPhi;   //!
   TBranch        *b_svMass;   //!
   TBranch        *b_svNTracks;   //!
   TBranch        *b_svChi2;   //!
   TBranch        *b_svNDF;   //!
   TBranch        *b_svDXY;   //!
   TBranch        *b_svDXYerr;   //!
   TBranch        *b_svD3D;   //!
   TBranch        *b_svD3Derr;   //!
   TBranch        *b_svCosThetaSVPS;   //!
   TBranch        *b_elesCharge;   //!
   TBranch        *b_elesMtw;   //!
   TBranch        *b_elesRelIso;   //!
   TBranch        *b_elesMiniIso;   //!
   TBranch        *b_elespfActivity;   //!
   TBranch        *b_recoJetsJecUnc;   //!
   TBranch        *b_recoJetsJecScaleRawToFull;   //!
   TBranch        *b_qgLikelihood;   //!
   TBranch        *b_qgPtD;   //!
   TBranch        *b_qgAxis2;   //!
   TBranch        *b_recoJetschargedHadronEnergyFraction;   //!
   TBranch        *b_recoJetschargedEmEnergyFraction;   //!
   TBranch        *b_recoJetsneutralEmEnergyFraction;   //!
   TBranch        *b_recoJetsmuonEnergyFraction;   //!
   TBranch        *b_recoJetsBtag_0;   //!
   TBranch        *b_recoJetsCharge_0;   //!
   TBranch        *b_DeepCSVb;   //!
   TBranch        *b_DeepCSVc;   //!
   TBranch        *b_DeepCSVl;   //!
   TBranch        *b_DeepCSVbb;   //!
   TBranch        *b_DeepCSVcc;   //!
   TBranch        *b_DeepCSVbN;   //!
   TBranch        *b_DeepCSVcN;   //!
   TBranch        *b_DeepCSVlN;   //!
   TBranch        *b_DeepCSVbbN;   //!
   TBranch        *b_DeepCSVccN;   //!
   TBranch        *b_DeepCSVbP;   //!
   TBranch        *b_DeepCSVcP;   //!
   TBranch        *b_DeepCSVlP;   //!
   TBranch        *b_DeepCSVbbP;   //!
   TBranch        *b_DeepCSVccP;   //!
   TBranch        *b_puppitau1;   //!
   TBranch        *b_puppitau2;   //!
   TBranch        *b_puppitau3;   //!
   TBranch        *b_puppisoftDropMass;   //!
   TBranch        *b_puppiSubJetsBdisc;   //!
   TBranch        *b_recoJetsJecUncLepCleaned;   //!
   TBranch        *b_prodJetsNoLep_qgLikelihood;   //!
   TBranch        *b_prodJetsNoLep_qgPtD;   //!
   TBranch        *b_prodJetsNoLep_qgAxis2;   //!
   TBranch        *b_recoJetschargedHadronEnergyFractionLepCleaned;   //!
   TBranch        *b_recoJetsneutralEmEnergyFractionLepCleaned;   //!
   TBranch        *b_recoJetschargedEmEnergyFractionLepCleaned;   //!
   TBranch        *b_recoJetsmuonEnergyFractionLepCleaned;   //!
   TBranch        *b_recoJetsBtag_0_LepCleaned;   //!
   TBranch        *b_recoJetsCharge_0_LepCleaned;   //!
   TBranch        *b_recoJetsJecScaleRawToFull_LepCleaned;   //!
   TBranch        *b_prodJetsNoLep_puppisoftDropMass;   //!
   TBranch        *b_prodJetsNoLep_puppitau1;   //!
   TBranch        *b_prodJetsNoLep_puppitau2;   //!
   TBranch        *b_prodJetsNoLep_puppitau3;   //!
   TBranch        *b_prodJetsNoLep_puppiSubJetsBdisc;   //!
   TBranch        *b_W_emu_pfActivityVec;   //!
   TBranch        *b_W_tau_emu_pfActivityVec;   //!
   TBranch        *b_W_tau_prongs_pfActivityVec;   //!
   TBranch        *b_trksForIsoVetocharge;   //!
   TBranch        *b_trksForIsoVetodz;   //!
   TBranch        *b_trksForIsoVetoiso;   //!
   TBranch        *b_trksForIsoVetopfActivity;   //!
   TBranch        *b_loose_isoTrks_charge;   //!
   TBranch        *b_loose_isoTrks_dz;   //!
   TBranch        *b_loose_isoTrks_iso;   //!
   TBranch        *b_loose_isoTrks_mtw;   //!
   TBranch        *b_loose_isoTrks_pfActivity;   //!
   TBranch        *b_metMagUp;   //!
   TBranch        *b_metMagDown;   //!
   TBranch        *b_metPhiUp;   //!
   TBranch        *b_metPhiDown;   //!
   TBranch        *b_PassTrigger;   //!
   TBranch        *b_TriggerPrescales;   //!
   TBranch        *b_muonsFlagMedium;   //!
   TBranch        *b_muonsFlagTight;   //!
   TBranch        *b_elesFlagMedium;   //!
   TBranch        *b_elesFlagVeto;   //!
   TBranch        *b_recoJetsFlavor;   //!
   TBranch        *b_qgMult;   //!
   TBranch        *b_muMatchedJetIdx;   //!
   TBranch        *b_eleMatchedJetIdx;   //!
   TBranch        *b_looseisoTrksMatchedJetIdx;   //!
   TBranch        *b_trksForIsoVetoMatchedJetIdx;   //!
   TBranch        *b_prodJetsNoLep_qgMult;   //!
   TBranch        *b_genDecayIdxVec;   //!
   TBranch        *b_genDecayPdgIdVec;   //!
   TBranch        *b_genDecayMomIdxVec;   //!
   TBranch        *b_genDecayMomRefVec;   //!
   TBranch        *b_W_emuVec;   //!
   TBranch        *b_W_tauVec;   //!
   TBranch        *b_W_tau_emuVec;   //!
   TBranch        *b_W_tau_prongsVec;   //!
   TBranch        *b_W_tau_nuVec;   //!
   TBranch        *b_selPDGid;   //!
   TBranch        *b_trksForIsoVetopdgId;   //!
   TBranch        *b_trksForIsoVetoidx;   //!
   TBranch        *b_loose_isoTrks_pdgId;   //!
   TBranch        *b_loose_isoTrks_idx;   //!
   TBranch        *b_forVetoIsoTrksidx;   //!
   TBranch        *b_loosePhotonID;   //!
   TBranch        *b_mediumPhotonID;   //!
   TBranch        *b_tightPhotonID;   //!
   TBranch        *b_nonPrompt;   //!
   TBranch        *b_elesisEB;   //!
   TBranch        *b_vetoElectronID;   //!
   TBranch        *b_looseElectronID;   //!
   TBranch        *b_mediumElectronID;   //!
   TBranch        *b_tightElectronID;   //!
   TBranch        *b_TriggerNames;   //!
   TBranch        *b_genDecayStrVec;   //!
   TBranch        *b_muonsLVec;   //!
   TBranch        *b_gammaLVec;   //!
   TBranch        *b_gammaLVecGen;   //!
   TBranch        *b_genPartonLVec;   //!
   TBranch        *b_svSoftLVec;   //!
   TBranch        *b_svLVec;   //!
   TBranch        *b_elesLVec;   //!
   TBranch        *b_jetsLVec;   //!
   TBranch        *b_puppiJetsLVec;   //!
   TBranch        *b_puppiSubJetsLVec;   //!
   TBranch        *b_jetsLVecLepCleaned;   //!
   TBranch        *b_prodJetsNoLep_puppiJetsLVec;   //!
   TBranch        *b_prodJetsNoLep_puppiSubJetsLVec;   //!
   TBranch        *b_genDecayLVec;   //!
   TBranch        *b_selGenParticle;   //!
   TBranch        *b_genjetsLVec;   //!
   TBranch        *b_loose_isoTrksLVec;   //!

   StopNtupleTree(TTree *tree=0);
   virtual ~StopNtupleTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

inline StopNtupleTree::StopNtupleTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SIG/TChiWZ.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SIG/TChiWZ.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("SIG/TChiWZ.root:/stopTreeMaker");
      dir->GetObject("AUX",tree);

   }
   Init(tree);
}

inline StopNtupleTree::~StopNtupleTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

inline Int_t StopNtupleTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
inline Long64_t StopNtupleTree::LoadTree(Long64_t entry)
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

inline void StopNtupleTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   muonsCharge = 0;
   muonsMtw = 0;
   muonsRelIso = 0;
   muonsMiniIso = 0;
   muonspfActivity = 0;
   pfGammaIso = 0;
   isEB = 0;
   genMatched = 0;
   hadTowOverEM = 0;
   sigmaIetaIeta = 0;
   pfChargedIso = 0;
   pfNeutralIso = 0;
   pfChargedIsoRhoCorr = 0;
   pfNeutralIsoRhoCorr = 0;
   pfGammaIsoRhoCorr = 0;
   hasPixelSeed = 0;
   passElectronVeto = 0;
   photonPt = 0;
   photonEta = 0;
   photonPhi = 0;
   svPT = 0;
   svETA = 0;
   svPhi = 0;
   svMass = 0;
   svNTracks = 0;
   svChi2 = 0;
   svNDF = 0;
   svDXY = 0;
   svDXYerr = 0;
   svD3D = 0;
   svD3Derr = 0;
   svCosThetaSVPS = 0;
   elesCharge = 0;
   elesMtw = 0;
   elesRelIso = 0;
   elesMiniIso = 0;
   elespfActivity = 0;
   recoJetsJecUnc = 0;
   recoJetsJecScaleRawToFull = 0;
   qgLikelihood = 0;
   qgPtD = 0;
   qgAxis2 = 0;
   recoJetschargedHadronEnergyFraction = 0;
   recoJetschargedEmEnergyFraction = 0;
   recoJetsneutralEmEnergyFraction = 0;
   recoJetsmuonEnergyFraction = 0;
   recoJetsBtag_0 = 0;
   recoJetsCharge_0 = 0;
   DeepCSVb = 0;
   DeepCSVc = 0;
   DeepCSVl = 0;
   DeepCSVbb = 0;
   DeepCSVcc = 0;
   DeepCSVbN = 0;
   DeepCSVcN = 0;
   DeepCSVlN = 0;
   DeepCSVbbN = 0;
   DeepCSVccN = 0;
   DeepCSVbP = 0;
   DeepCSVcP = 0;
   DeepCSVlP = 0;
   DeepCSVbbP = 0;
   DeepCSVccP = 0;
   puppitau1 = 0;
   puppitau2 = 0;
   puppitau3 = 0;
   puppisoftDropMass = 0;
   puppiSubJetsBdisc = 0;
   recoJetsJecUncLepCleaned = 0;
   prodJetsNoLep_qgLikelihood = 0;
   prodJetsNoLep_qgPtD = 0;
   prodJetsNoLep_qgAxis2 = 0;
   recoJetschargedHadronEnergyFractionLepCleaned = 0;
   recoJetsneutralEmEnergyFractionLepCleaned = 0;
   recoJetschargedEmEnergyFractionLepCleaned = 0;
   recoJetsmuonEnergyFractionLepCleaned = 0;
   recoJetsBtag_0_LepCleaned = 0;
   recoJetsCharge_0_LepCleaned = 0;
   recoJetsJecScaleRawToFull_LepCleaned = 0;
   prodJetsNoLep_puppisoftDropMass = 0;
   prodJetsNoLep_puppitau1 = 0;
   prodJetsNoLep_puppitau2 = 0;
   prodJetsNoLep_puppitau3 = 0;
   prodJetsNoLep_puppiSubJetsBdisc = 0;
   W_emu_pfActivityVec = 0;
   W_tau_emu_pfActivityVec = 0;
   W_tau_prongs_pfActivityVec = 0;
   trksForIsoVetocharge = 0;
   trksForIsoVetodz = 0;
   trksForIsoVetoiso = 0;
   trksForIsoVetopfActivity = 0;
   loose_isoTrks_charge = 0;
   loose_isoTrks_dz = 0;
   loose_isoTrks_iso = 0;
   loose_isoTrks_mtw = 0;
   loose_isoTrks_pfActivity = 0;
   metMagUp = 0;
   metMagDown = 0;
   metPhiUp = 0;
   metPhiDown = 0;
   PassTrigger = 0;
   TriggerPrescales = 0;
   muonsFlagMedium = 0;
   muonsFlagTight = 0;
   elesFlagMedium = 0;
   elesFlagVeto = 0;
   recoJetsFlavor = 0;
   qgMult = 0;
   muMatchedJetIdx = 0;
   eleMatchedJetIdx = 0;
   looseisoTrksMatchedJetIdx = 0;
   trksForIsoVetoMatchedJetIdx = 0;
   prodJetsNoLep_qgMult = 0;
   genDecayIdxVec = 0;
   genDecayPdgIdVec = 0;
   genDecayMomIdxVec = 0;
   genDecayMomRefVec = 0;
   W_emuVec = 0;
   W_tauVec = 0;
   W_tau_emuVec = 0;
   W_tau_prongsVec = 0;
   W_tau_nuVec = 0;
   selPDGid = 0;
   trksForIsoVetopdgId = 0;
   trksForIsoVetoidx = 0;
   loose_isoTrks_pdgId = 0;
   loose_isoTrks_idx = 0;
   forVetoIsoTrksidx = 0;
   loosePhotonID = 0;
   mediumPhotonID = 0;
   tightPhotonID = 0;
   nonPrompt = 0;
   elesisEB = 0;
   vetoElectronID = 0;
   looseElectronID = 0;
   mediumElectronID = 0;
   tightElectronID = 0;
   TriggerNames = 0;
   genDecayStrVec = 0;
   muonsLVec = 0;
   gammaLVec = 0;
   gammaLVecGen = 0;
   genPartonLVec = 0;
   svSoftLVec = 0;
   svLVec = 0;
   elesLVec = 0;
   jetsLVec = 0;
   puppiJetsLVec = 0;
   puppiSubJetsLVec = 0;
   jetsLVecLepCleaned = 0;
   prodJetsNoLep_puppiJetsLVec = 0;
   prodJetsNoLep_puppiSubJetsLVec = 0;
   genDecayLVec = 0;
   selGenParticle = 0;
   genjetsLVec = 0;
   loose_isoTrksLVec = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   m_USEFLOAT = string(fChain->GetBranch("evtWeight")->GetLeaf("evtWeight")->GetTypeName()) == "Float_t";
   // "special" declarations
   if(m_USEFLOAT){
     fChain->SetBranchAddress("met", &met_f, &b_met);
     fChain->SetBranchAddress("metphi", &metphi_f, &b_metphi);
     fChain->SetBranchAddress("calomet", &calomet_f, &b_calomet);
     fChain->SetBranchAddress("calometphi", &calometphi_f, &b_calometphi);
     fChain->SetBranchAddress("stored_weight", &stored_weight_f, &b_stored_weight);
     fChain->SetBranchAddress("evtWeight", &evtWeight_f, &b_evtWeight);
   } else {
     fChain->SetBranchAddress("met", &met_d, &b_met);
     fChain->SetBranchAddress("metphi", &metphi_d, &b_metphi);
     fChain->SetBranchAddress("calomet", &calomet_d, &b_calomet);
     fChain->SetBranchAddress("calometphi", &calometphi_d, &b_calometphi);
     fChain->SetBranchAddress("stored_weight", &stored_weight_d, &b_stored_weight);
     fChain->SetBranchAddress("evtWeight", &evtWeight_d, &b_evtWeight);
   }
   
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("genmet", &genmet, &b_genmet);
   fChain->SetBranchAddress("genmetphi", &genmetphi, &b_genmetphi);
   fChain->SetBranchAddress("mht", &mht, &b_mht);
   fChain->SetBranchAddress("mhtphi", &mhtphi, &b_mhtphi);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("dPhi0_CUT", &dPhi0_CUT, &b_dPhi0_CUT);
   fChain->SetBranchAddress("dPhi1_CUT", &dPhi1_CUT, &b_dPhi1_CUT);
   fChain->SetBranchAddress("dPhi2_CUT", &dPhi2_CUT, &b_dPhi2_CUT);
   fChain->SetBranchAddress("tru_npv", &tru_npv, &b_tru_npv);
   fChain->SetBranchAddress("avg_npv", &avg_npv, &b_avg_npv);
   fChain->SetBranchAddress("METFilters", &METFilters, &b_METFilters);
   fChain->SetBranchAddress("CSCTightHaloFilter", &CSCTightHaloFilter, &b_CSCTightHaloFilter);
   fChain->SetBranchAddress("globalSuperTightHalo2016Filter", &globalSuperTightHalo2016Filter, &b_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("goodVerticesFilter", &goodVerticesFilter, &b_goodVerticesFilter);
   fChain->SetBranchAddress("ecalBadCalibFilter", &ecalBadCalibFilter, &b_ecalBadCalibFilter);
   fChain->SetBranchAddress("HBHENoiseIsoFilter", &HBHENoiseIsoFilter, &b_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("EcalDeadCellTriggerPrimitiveFilter", &EcalDeadCellTriggerPrimitiveFilter, &b_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("noBadMuonsFilter", &noBadMuonsFilter, &b_noBadMuonsFilter);
   fChain->SetBranchAddress("badMuonsFilter", &badMuonsFilter, &b_badMuonsFilter);
   fChain->SetBranchAddress("duplicateMuonsFilter", &duplicateMuonsFilter, &b_duplicateMuonsFilter);
   fChain->SetBranchAddress("nMuons_CUT", &nMuons_CUT, &b_nMuons_CUT);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("nElectrons_CUT", &nElectrons_CUT, &b_nElectrons_CUT);
   fChain->SetBranchAddress("nElectrons", &nElectrons, &b_nElectrons);
   fChain->SetBranchAddress("NJetsISR", &NJetsISR, &b_NJetsISR);
   fChain->SetBranchAddress("loose_nIsoTrks", &loose_nIsoTrks, &b_loose_nIsoTrks);
   fChain->SetBranchAddress("nIsoTrks_CUT", &nIsoTrks_CUT, &b_nIsoTrks_CUT);
   fChain->SetBranchAddress("nJets_CUT", &nJets_CUT, &b_nJets_CUT);
   fChain->SetBranchAddress("vtxSize", &vtxSize, &b_vtxSize);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("nm1", &nm1, &b_nm1);
   fChain->SetBranchAddress("n0", &n0, &b_n0);
   fChain->SetBranchAddress("np1", &np1, &b_np1);
   fChain->SetBranchAddress("looseJetID", &looseJetID, &b_looseJetID);
   fChain->SetBranchAddress("tightJetID", &tightJetID, &b_tightJetID);
   fChain->SetBranchAddress("tightlepvetoJetID", &tightlepvetoJetID, &b_tightlepvetoJetID);
   fChain->SetBranchAddress("looseJetID_NoLep", &looseJetID_NoLep, &b_looseJetID_NoLep);
   fChain->SetBranchAddress("tightJetID_NoLep", &tightJetID_NoLep, &b_tightJetID_NoLep);
   fChain->SetBranchAddress("tightlepvetoJetID_NoLep", &tightlepvetoJetID_NoLep, &b_tightlepvetoJetID_NoLep);
   fChain->SetBranchAddress("BadChargedCandidateFilter", &BadChargedCandidateFilter, &b_BadChargedCandidateFilter);
   fChain->SetBranchAddress("BadPFMuonFilter", &BadPFMuonFilter, &b_BadPFMuonFilter);
   fChain->SetBranchAddress("HBHENoiseFilter", &HBHENoiseFilter, &b_HBHENoiseFilter);
   fChain->SetBranchAddress("HBHEIsoNoiseFilter", &HBHEIsoNoiseFilter, &b_HBHEIsoNoiseFilter);
   fChain->SetBranchAddress("muonsCharge", &muonsCharge, &b_muonsCharge);
   fChain->SetBranchAddress("muonsMtw", &muonsMtw, &b_muonsMtw);
   fChain->SetBranchAddress("muonsRelIso", &muonsRelIso, &b_muonsRelIso);
   fChain->SetBranchAddress("muonsMiniIso", &muonsMiniIso, &b_muonsMiniIso);
   fChain->SetBranchAddress("muonspfActivity", &muonspfActivity, &b_muonspfActivity);
   fChain->SetBranchAddress("pfGammaIso", &pfGammaIso, &b_pfGammaIso);
   fChain->SetBranchAddress("isEB", &isEB, &b_isEB);
   fChain->SetBranchAddress("genMatched", &genMatched, &b_genMatched);
   fChain->SetBranchAddress("hadTowOverEM", &hadTowOverEM, &b_hadTowOverEM);
   fChain->SetBranchAddress("sigmaIetaIeta", &sigmaIetaIeta, &b_sigmaIetaIeta);
   fChain->SetBranchAddress("pfChargedIso", &pfChargedIso, &b_pfChargedIso);
   fChain->SetBranchAddress("pfNeutralIso", &pfNeutralIso, &b_pfNeutralIso);
   fChain->SetBranchAddress("pfChargedIsoRhoCorr", &pfChargedIsoRhoCorr, &b_pfChargedIsoRhoCorr);
   fChain->SetBranchAddress("pfNeutralIsoRhoCorr", &pfNeutralIsoRhoCorr, &b_pfNeutralIsoRhoCorr);
   fChain->SetBranchAddress("pfGammaIsoRhoCorr", &pfGammaIsoRhoCorr, &b_pfGammaIsoRhoCorr);
   fChain->SetBranchAddress("hasPixelSeed", &hasPixelSeed, &b_hasPixelSeed);
   fChain->SetBranchAddress("passElectronVeto", &passElectronVeto, &b_passElectronVeto);
   fChain->SetBranchAddress("photonPt", &photonPt, &b_photonPt);
   fChain->SetBranchAddress("photonEta", &photonEta, &b_photonEta);
   fChain->SetBranchAddress("photonPhi", &photonPhi, &b_photonPhi);
   fChain->SetBranchAddress("svPT", &svPT, &b_svPT);
   fChain->SetBranchAddress("svETA", &svETA, &b_svETA);
   fChain->SetBranchAddress("svPhi", &svPhi, &b_svPhi);
   fChain->SetBranchAddress("svMass", &svMass, &b_svMass);
   fChain->SetBranchAddress("svNTracks", &svNTracks, &b_svNTracks);
   fChain->SetBranchAddress("svChi2", &svChi2, &b_svChi2);
   fChain->SetBranchAddress("svNDF", &svNDF, &b_svNDF);
   fChain->SetBranchAddress("svDXY", &svDXY, &b_svDXY);
   fChain->SetBranchAddress("svDXYerr", &svDXYerr, &b_svDXYerr);
   fChain->SetBranchAddress("svD3D", &svD3D, &b_svD3D);
   fChain->SetBranchAddress("svD3Derr", &svD3Derr, &b_svD3Derr);
   fChain->SetBranchAddress("svCosThetaSVPS", &svCosThetaSVPS, &b_svCosThetaSVPS);
   fChain->SetBranchAddress("elesCharge", &elesCharge, &b_elesCharge);
   fChain->SetBranchAddress("elesMtw", &elesMtw, &b_elesMtw);
   fChain->SetBranchAddress("elesRelIso", &elesRelIso, &b_elesRelIso);
   fChain->SetBranchAddress("elesMiniIso", &elesMiniIso, &b_elesMiniIso);
   fChain->SetBranchAddress("elespfActivity", &elespfActivity, &b_elespfActivity);
   fChain->SetBranchAddress("recoJetsJecUnc", &recoJetsJecUnc, &b_recoJetsJecUnc);
   fChain->SetBranchAddress("recoJetsJecScaleRawToFull", &recoJetsJecScaleRawToFull, &b_recoJetsJecScaleRawToFull);
   fChain->SetBranchAddress("qgLikelihood", &qgLikelihood, &b_qgLikelihood);
   fChain->SetBranchAddress("qgPtD", &qgPtD, &b_qgPtD);
   fChain->SetBranchAddress("qgAxis2", &qgAxis2, &b_qgAxis2);
   fChain->SetBranchAddress("recoJetschargedHadronEnergyFraction", &recoJetschargedHadronEnergyFraction, &b_recoJetschargedHadronEnergyFraction);
   fChain->SetBranchAddress("recoJetschargedEmEnergyFraction", &recoJetschargedEmEnergyFraction, &b_recoJetschargedEmEnergyFraction);
   fChain->SetBranchAddress("recoJetsneutralEmEnergyFraction", &recoJetsneutralEmEnergyFraction, &b_recoJetsneutralEmEnergyFraction);
   fChain->SetBranchAddress("recoJetsmuonEnergyFraction", &recoJetsmuonEnergyFraction, &b_recoJetsmuonEnergyFraction);
   fChain->SetBranchAddress("recoJetsBtag_0", &recoJetsBtag_0, &b_recoJetsBtag_0);
   fChain->SetBranchAddress("recoJetsCharge_0", &recoJetsCharge_0, &b_recoJetsCharge_0);
   fChain->SetBranchAddress("DeepCSVb", &DeepCSVb, &b_DeepCSVb);
   fChain->SetBranchAddress("DeepCSVc", &DeepCSVc, &b_DeepCSVc);
   fChain->SetBranchAddress("DeepCSVl", &DeepCSVl, &b_DeepCSVl);
   fChain->SetBranchAddress("DeepCSVbb", &DeepCSVbb, &b_DeepCSVbb);
   fChain->SetBranchAddress("DeepCSVcc", &DeepCSVcc, &b_DeepCSVcc);
   fChain->SetBranchAddress("DeepCSVbN", &DeepCSVbN, &b_DeepCSVbN);
   fChain->SetBranchAddress("DeepCSVcN", &DeepCSVcN, &b_DeepCSVcN);
   fChain->SetBranchAddress("DeepCSVlN", &DeepCSVlN, &b_DeepCSVlN);
   fChain->SetBranchAddress("DeepCSVbbN", &DeepCSVbbN, &b_DeepCSVbbN);
   fChain->SetBranchAddress("DeepCSVccN", &DeepCSVccN, &b_DeepCSVccN);
   fChain->SetBranchAddress("DeepCSVbP", &DeepCSVbP, &b_DeepCSVbP);
   fChain->SetBranchAddress("DeepCSVcP", &DeepCSVcP, &b_DeepCSVcP);
   fChain->SetBranchAddress("DeepCSVlP", &DeepCSVlP, &b_DeepCSVlP);
   fChain->SetBranchAddress("DeepCSVbbP", &DeepCSVbbP, &b_DeepCSVbbP);
   fChain->SetBranchAddress("DeepCSVccP", &DeepCSVccP, &b_DeepCSVccP);
   fChain->SetBranchAddress("puppitau1", &puppitau1, &b_puppitau1);
   fChain->SetBranchAddress("puppitau2", &puppitau2, &b_puppitau2);
   fChain->SetBranchAddress("puppitau3", &puppitau3, &b_puppitau3);
   fChain->SetBranchAddress("puppisoftDropMass", &puppisoftDropMass, &b_puppisoftDropMass);
   fChain->SetBranchAddress("puppiSubJetsBdisc", &puppiSubJetsBdisc, &b_puppiSubJetsBdisc);
   fChain->SetBranchAddress("recoJetsJecUncLepCleaned", &recoJetsJecUncLepCleaned, &b_recoJetsJecUncLepCleaned);
   fChain->SetBranchAddress("prodJetsNoLep_qgLikelihood", &prodJetsNoLep_qgLikelihood, &b_prodJetsNoLep_qgLikelihood);
   fChain->SetBranchAddress("prodJetsNoLep_qgPtD", &prodJetsNoLep_qgPtD, &b_prodJetsNoLep_qgPtD);
   fChain->SetBranchAddress("prodJetsNoLep_qgAxis2", &prodJetsNoLep_qgAxis2, &b_prodJetsNoLep_qgAxis2);
   fChain->SetBranchAddress("recoJetschargedHadronEnergyFractionLepCleaned", &recoJetschargedHadronEnergyFractionLepCleaned, &b_recoJetschargedHadronEnergyFractionLepCleaned);
   fChain->SetBranchAddress("recoJetsneutralEmEnergyFractionLepCleaned", &recoJetsneutralEmEnergyFractionLepCleaned, &b_recoJetsneutralEmEnergyFractionLepCleaned);
   fChain->SetBranchAddress("recoJetschargedEmEnergyFractionLepCleaned", &recoJetschargedEmEnergyFractionLepCleaned, &b_recoJetschargedEmEnergyFractionLepCleaned);
   fChain->SetBranchAddress("recoJetsmuonEnergyFractionLepCleaned", &recoJetsmuonEnergyFractionLepCleaned, &b_recoJetsmuonEnergyFractionLepCleaned);
   fChain->SetBranchAddress("recoJetsBtag_0_LepCleaned", &recoJetsBtag_0_LepCleaned, &b_recoJetsBtag_0_LepCleaned);
   fChain->SetBranchAddress("recoJetsCharge_0_LepCleaned", &recoJetsCharge_0_LepCleaned, &b_recoJetsCharge_0_LepCleaned);
   fChain->SetBranchAddress("recoJetsJecScaleRawToFull_LepCleaned", &recoJetsJecScaleRawToFull_LepCleaned, &b_recoJetsJecScaleRawToFull_LepCleaned);
   fChain->SetBranchAddress("prodJetsNoLep_puppisoftDropMass", &prodJetsNoLep_puppisoftDropMass, &b_prodJetsNoLep_puppisoftDropMass);
   fChain->SetBranchAddress("prodJetsNoLep_puppitau1", &prodJetsNoLep_puppitau1, &b_prodJetsNoLep_puppitau1);
   fChain->SetBranchAddress("prodJetsNoLep_puppitau2", &prodJetsNoLep_puppitau2, &b_prodJetsNoLep_puppitau2);
   fChain->SetBranchAddress("prodJetsNoLep_puppitau3", &prodJetsNoLep_puppitau3, &b_prodJetsNoLep_puppitau3);
   fChain->SetBranchAddress("prodJetsNoLep_puppiSubJetsBdisc", &prodJetsNoLep_puppiSubJetsBdisc, &b_prodJetsNoLep_puppiSubJetsBdisc);
   fChain->SetBranchAddress("W_emu_pfActivityVec", &W_emu_pfActivityVec, &b_W_emu_pfActivityVec);
   fChain->SetBranchAddress("W_tau_emu_pfActivityVec", &W_tau_emu_pfActivityVec, &b_W_tau_emu_pfActivityVec);
   fChain->SetBranchAddress("W_tau_prongs_pfActivityVec", &W_tau_prongs_pfActivityVec, &b_W_tau_prongs_pfActivityVec);
   fChain->SetBranchAddress("trksForIsoVetocharge", &trksForIsoVetocharge, &b_trksForIsoVetocharge);
   fChain->SetBranchAddress("trksForIsoVetodz", &trksForIsoVetodz, &b_trksForIsoVetodz);
   fChain->SetBranchAddress("trksForIsoVetoiso", &trksForIsoVetoiso, &b_trksForIsoVetoiso);
   fChain->SetBranchAddress("trksForIsoVetopfActivity", &trksForIsoVetopfActivity, &b_trksForIsoVetopfActivity);
   fChain->SetBranchAddress("loose_isoTrks_charge", &loose_isoTrks_charge, &b_loose_isoTrks_charge);
   fChain->SetBranchAddress("loose_isoTrks_dz", &loose_isoTrks_dz, &b_loose_isoTrks_dz);
   fChain->SetBranchAddress("loose_isoTrks_iso", &loose_isoTrks_iso, &b_loose_isoTrks_iso);
   fChain->SetBranchAddress("loose_isoTrks_mtw", &loose_isoTrks_mtw, &b_loose_isoTrks_mtw);
   fChain->SetBranchAddress("loose_isoTrks_pfActivity", &loose_isoTrks_pfActivity, &b_loose_isoTrks_pfActivity);
   fChain->SetBranchAddress("metMagUp", &metMagUp, &b_metMagUp);
   fChain->SetBranchAddress("metMagDown", &metMagDown, &b_metMagDown);
   fChain->SetBranchAddress("metPhiUp", &metPhiUp, &b_metPhiUp);
   fChain->SetBranchAddress("metPhiDown", &metPhiDown, &b_metPhiDown);
   fChain->SetBranchAddress("PassTrigger", &PassTrigger, &b_PassTrigger);
   fChain->SetBranchAddress("TriggerPrescales", &TriggerPrescales, &b_TriggerPrescales);
   fChain->SetBranchAddress("muonsFlagMedium", &muonsFlagMedium, &b_muonsFlagMedium);
   fChain->SetBranchAddress("muonsFlagTight", &muonsFlagTight, &b_muonsFlagTight);
   fChain->SetBranchAddress("elesFlagMedium", &elesFlagMedium, &b_elesFlagMedium);
   fChain->SetBranchAddress("elesFlagVeto", &elesFlagVeto, &b_elesFlagVeto);
   fChain->SetBranchAddress("recoJetsFlavor", &recoJetsFlavor, &b_recoJetsFlavor);
   fChain->SetBranchAddress("qgMult", &qgMult, &b_qgMult);
   fChain->SetBranchAddress("muMatchedJetIdx", &muMatchedJetIdx, &b_muMatchedJetIdx);
   fChain->SetBranchAddress("eleMatchedJetIdx", &eleMatchedJetIdx, &b_eleMatchedJetIdx);
   fChain->SetBranchAddress("looseisoTrksMatchedJetIdx", &looseisoTrksMatchedJetIdx, &b_looseisoTrksMatchedJetIdx);
   fChain->SetBranchAddress("trksForIsoVetoMatchedJetIdx", &trksForIsoVetoMatchedJetIdx, &b_trksForIsoVetoMatchedJetIdx);
   fChain->SetBranchAddress("prodJetsNoLep_qgMult", &prodJetsNoLep_qgMult, &b_prodJetsNoLep_qgMult);
   fChain->SetBranchAddress("genDecayIdxVec", &genDecayIdxVec, &b_genDecayIdxVec);
   fChain->SetBranchAddress("genDecayPdgIdVec", &genDecayPdgIdVec, &b_genDecayPdgIdVec);
   fChain->SetBranchAddress("genDecayMomIdxVec", &genDecayMomIdxVec, &b_genDecayMomIdxVec);
   fChain->SetBranchAddress("genDecayMomRefVec", &genDecayMomRefVec, &b_genDecayMomRefVec);
   fChain->SetBranchAddress("W_emuVec", &W_emuVec, &b_W_emuVec);
   fChain->SetBranchAddress("W_tauVec", &W_tauVec, &b_W_tauVec);
   fChain->SetBranchAddress("W_tau_emuVec", &W_tau_emuVec, &b_W_tau_emuVec);
   fChain->SetBranchAddress("W_tau_prongsVec", &W_tau_prongsVec, &b_W_tau_prongsVec);
   fChain->SetBranchAddress("W_tau_nuVec", &W_tau_nuVec, &b_W_tau_nuVec);
   fChain->SetBranchAddress("selPDGid", &selPDGid, &b_selPDGid);
   fChain->SetBranchAddress("trksForIsoVetopdgId", &trksForIsoVetopdgId, &b_trksForIsoVetopdgId);
   fChain->SetBranchAddress("trksForIsoVetoidx", &trksForIsoVetoidx, &b_trksForIsoVetoidx);
   fChain->SetBranchAddress("loose_isoTrks_pdgId", &loose_isoTrks_pdgId, &b_loose_isoTrks_pdgId);
   fChain->SetBranchAddress("loose_isoTrks_idx", &loose_isoTrks_idx, &b_loose_isoTrks_idx);
   fChain->SetBranchAddress("forVetoIsoTrksidx", &forVetoIsoTrksidx, &b_forVetoIsoTrksidx);
   fChain->SetBranchAddress("loosePhotonID", &loosePhotonID, &b_loosePhotonID);
   fChain->SetBranchAddress("mediumPhotonID", &mediumPhotonID, &b_mediumPhotonID);
   fChain->SetBranchAddress("tightPhotonID", &tightPhotonID, &b_tightPhotonID);
   fChain->SetBranchAddress("nonPrompt", &nonPrompt, &b_nonPrompt);
   fChain->SetBranchAddress("elesisEB", &elesisEB, &b_elesisEB);
   fChain->SetBranchAddress("vetoElectronID", &vetoElectronID, &b_vetoElectronID);
   fChain->SetBranchAddress("looseElectronID", &looseElectronID, &b_looseElectronID);
   fChain->SetBranchAddress("mediumElectronID", &mediumElectronID, &b_mediumElectronID);
   fChain->SetBranchAddress("tightElectronID", &tightElectronID, &b_tightElectronID);
   fChain->SetBranchAddress("TriggerNames", &TriggerNames, &b_TriggerNames);
   fChain->SetBranchAddress("genDecayStrVec", &genDecayStrVec, &b_genDecayStrVec);
   fChain->SetBranchAddress("muonsLVec", &muonsLVec, &b_muonsLVec);
   fChain->SetBranchAddress("gammaLVec", &gammaLVec, &b_gammaLVec);
   fChain->SetBranchAddress("gammaLVecGen", &gammaLVecGen, &b_gammaLVecGen);
   fChain->SetBranchAddress("genPartonLVec", &genPartonLVec, &b_genPartonLVec);
   fChain->SetBranchAddress("svSoftLVec", &svSoftLVec, &b_svSoftLVec);
   fChain->SetBranchAddress("svLVec", &svLVec, &b_svLVec);
   fChain->SetBranchAddress("elesLVec", &elesLVec, &b_elesLVec);
   fChain->SetBranchAddress("jetsLVec", &jetsLVec, &b_jetsLVec);
   fChain->SetBranchAddress("puppiJetsLVec", &puppiJetsLVec, &b_puppiJetsLVec);
   fChain->SetBranchAddress("puppiSubJetsLVec", &puppiSubJetsLVec, &b_puppiSubJetsLVec);
   fChain->SetBranchAddress("jetsLVecLepCleaned", &jetsLVecLepCleaned, &b_jetsLVecLepCleaned);
   fChain->SetBranchAddress("prodJetsNoLep_puppiJetsLVec", &prodJetsNoLep_puppiJetsLVec, &b_prodJetsNoLep_puppiJetsLVec);
   fChain->SetBranchAddress("prodJetsNoLep_puppiSubJetsLVec", &prodJetsNoLep_puppiSubJetsLVec, &b_prodJetsNoLep_puppiSubJetsLVec);
   fChain->SetBranchAddress("genDecayLVec", &genDecayLVec, &b_genDecayLVec);
   fChain->SetBranchAddress("selGenParticle", &selGenParticle, &b_selGenParticle);
   fChain->SetBranchAddress("genjetsLVec", &genjetsLVec, &b_genjetsLVec);
   fChain->SetBranchAddress("loose_isoTrksLVec", &loose_isoTrksLVec, &b_loose_isoTrksLVec);
   Notify();
}

inline Bool_t StopNtupleTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

inline void StopNtupleTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
inline Int_t StopNtupleTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

