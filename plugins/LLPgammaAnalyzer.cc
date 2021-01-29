// -*- C++ -*-
//
// Package:    LLPgamma/LLPgammaAnalyzer
// Class:      LLPgammaAnalyzer
//
/**\class LLPgammaAnalyzer LLPgammaAnalyzer.cc LLPgamma/LLPgammaAnalyzer/plugins/LLPgammaAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//----------------------------------------  cc file   --------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------

#include "LLPgamma/LLPgammaAnalyzer/plugins/LLPgammaAnalyzer.hh"
using namespace std;

//
// constructors and destructor
//
LLPgammaAnalyzer::LLPgammaAnalyzer(const edm::ParameterSet& iConfig):

  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

// -- declare tags ----------------------------------------------------------
  // triggers
  triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
  triggerObjectsTag(iConfig.getParameter<edm::InputTag>("triggerObjects")),

  // tracks
  tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),

  // vertices
  verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),

  // rho
  rhoTag(iConfig.getParameter<edm::InputTag>("rho")),

  // mets
  metsTag(iConfig.getParameter<edm::InputTag>("mets")),  

  // jets
  jetsTag(iConfig.getParameter<edm::InputTag>("jets")), 

  // electrons
  electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),  

  // muons
  muonsTag(iConfig.getParameter<edm::InputTag>("muons")),  

  // recHits
  recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),  
  recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

  // gedphotons
  gedPhotonsTag(iConfig.getParameter<edm::InputTag>("gedPhotons")),

  // ootPhotons
  ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons")),

// -- end of tag declarations ---------------------------------------
{

  usesResource();
  usesResource("TFileService");

// -- consume tags ------------------------------------------------------------

  // Triggers
  triggerResultsToken_	=consumes<edm::TriggerResults>(triggerResultsTag);
  triggerObjectsToken_	=consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectsTag);

  // tracks 
  tracksToken_				=consumes<std::vector<reco::Track>>(tracksTag);

  // vertices
  verticesToken_			=consumes<std::vector<reco::Vertex>>(verticesTag);
  
  // rho
  rhoToken_					=consumes<double>(rhoTag);

  // mets
  metsToken_				=consumes<std::vector<pat::MET>>(metsTag);
  
  // jets
  jetsToken_				=consumes<std::vector<pat::Jet>>(jetsTag);
  
  // leptons
  electronsToken_			=consumes<std::vector<pat::Electron>>(electronsTag);
  muonsToken_				=consumes<std::vector<pat::Muon>>(muonsTag);

  // rechits
  recHitsEBToken_			=consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
  recHitsEEToken_			=consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

  // photons
  gedPhotonsToken_		=consumes<std::vector<pat::Photon>>(gedPhotonsTag);
  ootPhotonsToken_		=consumes<std::vector<pat::Photon>>(ootPhotonsTag);

// ---------------------------------------------------------------------------------

}


LLPgammaAnalyzer::~LLPgammaAnalyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

int GetPFJetID(const pat::Jet & jet){ // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017

    const auto eta = std::abs(jet.eta());
    
    const auto NHF  = jet.neutralHadronEnergyFraction();
    const auto NEMF = jet.neutralEmEnergyFraction();
    const auto CHF  = jet.chargedHadronEnergyFraction();
    const auto CEMF = jet.chargedEmEnergyFraction();
    const auto NHM  = jet.neutralMultiplicity();
    const auto CHM  = jet.chargedMultiplicity();
    const auto SHM  = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    const auto MUF  = jet.muonEnergyFraction();
    
    // 2 == TightLepVeto
    // 1 == Tight

    if (eta <= 2.4)
    {
      if      ((NHF < 0.90) && (NEMF < 0.90) && (SHM > 1) && (MUF < 0.80) && (CHF > 0) && (CHM > 0) && (CEMF < 0.80)) return 2;
      else if ((NHF < 0.90) && (NEMF < 0.90) && (SHM > 1) &&                 (CHF > 0) && (CHM > 0))                  return 1;
      else                                                                                                            return 0; 
    }
    else if (eta > 2.4 && eta <= 2.7)
    {
      if      ((NHF < 0.90) && (NEMF < 0.90) && (SHM > 1) && (MUF < 0.80)) return 2;
      else if ((NHF < 0.90) && (NEMF < 0.90) && (SHM > 1))                 return 1;
      else                                                                 return 0; 
    }
    else if (eta > 2.7 && eta <= 3.0)
    {
      if   ((NEMF > 0.02) && (NEMF < 0.99) && (NHM > 2)) return 1;
      else                                               return 0; 
    }
    else 
    {
      if   ((NEMF < 0.90) && (NHF > 0.02) && (NHM > 10)) return 1;
      else                                               return 0; 
    }

    return -1; // should not happen
}



// ------------ method called for each event  ------------
void LLPgammaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

// -- Consume Tokens --------------------------------------------
  // TRIGGER
  iEvent.getByToken(triggerResultsToken,triggerResults_);
  iEvent.getByToken(triggerObjectsToken,triggerObjects_);


  // TRACKS
  iEvent.getByToken(tracksToken,tracks_);

  // VERTICES
  iEvent.getByToken(verticesToken,vertices_);

  // RHO
  iEvent.getByToken(rhoToken,rho_);

  // METS
  iEvent.getByToken(metsToken,mets_);

  // JETS
  iEvent.getByToken(jetsToken,jets_);

  // LEPTONS
  iEvent.getByToken(electronsToken,electrons_);
  iEvent.getByToken(muonsToken,muons_);

  // ECAL RECHITS
  iEvent.getByToken(recHitsEBToken,recHitsEB_);
  iEvent.getByToken(recHitsEEToken,recHitsEE_);

// -- Process Objects ------------------------------------------

// ---- process jets --------------------------

// ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
// ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 

   jets.clear(); 
   jets.reserve(jets_->size());

   for(const auto& jet : *jets_ ){ // Filters jet collection & sorts by pt

      if (jet.pt() < jetpTmin) continue;
      if (std::abs(jet.eta()) > jetEtamax) continue;
      
      const auto jetID = oot::GetPFJetID(jet);
      if (jetID < jetIDmin) continue;

      // save the jets, and then store the ID
      jets.emplace_back(jet);
      jets.back().addUserInt("jetID",jetID);
      
      std::sort(jets.begin(),jets.end(),oot::sortByPt);
   }

   nJets = jets.size();
   auto stJets = nJets; // allows to set the number of jets to skim ( = nJets for all )

   jetE.clear();
   jetpt.clear(); 
   jetphi.clear(); 
   jeteta.clear(); 
   jetID	.clear();
   jetNHF.clear();
   jetNEMF.clear();	
   jetCHF.clear();
   jetCEMF.clear();
   jetMUF.clear();
   jetNHM.clear();
   jetCHM.clear();
 
   for (auto i = 0; i < stJets; i++){

      const auto & jet = jets[i];

      jetE.push_back(	jet.energy());
      jetpt.push_back(	jet.pt());
      jetphi.push_back(	jet.phi());
      jeteta.push_back(	jet.eta());
      jetID.push_back(	jet.userInt("jetID"));
      jetNHF.push_back(	jet.neutralHadronEnergyFraction());
      jetNEMF.push_back(jet.neutralEmEnergyFraction());
      jetCHF.push_back(	jet.chargedHadronEnergyFraction());
      jetCEMF.push_back(jet.chargedEmEnergyFraction());
      jetMUF.push_back(	jet.muonEnergyFraction());
      jetNHM.push_back(	jet.neutralMultiplicity());
      jetCHM.push_back(	jet.chargedMultiplicity());
 
   } 

// ** end of jets  *************************************



// -- Fill output trees ------------------------------------------

  outTree->Fill();

// -- EOFun ------------------------------------------------------
//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//   ESHandle<SetupData> pSetup;
//   iSetup.get<SetupRecord>().get(pSetup);
//#endif
}


// ------------ method called once each job just before starting event loop  ------------
void LLPgammaAnalyzer::beginJob()
{

  edm::Service<TFileService> fs;
  outTree = fs->make<TTree>("llpgtree","llpgtree");

  // Run, Lumi, Event info
  outTree->Branch("run", 		&run);
  outTree->Branch("lumi", 		&lumi);
  outTree->Branch("event", 	&event, "event/l");

  // Jet info
  outTree->Branch("njets", 	&nJets);
  outTree->Branch("jetE", 		&jetE);
  outTree->Branch("jetpt", 	&jetpt);
  outTree->Branch("jeteta", 	&jeteta);
  outTree->Branch("jetphi", 	&jetphi);
  outTree->Branch("jetID", 	&jetID);
  outTree->Branch("jetNHF", 	&jetNHF);
  outTree->Branch("jetNEMF", 	&jetNEMF);  
  outTree->Branch("jetCHF", 	&jetCHF);
  outTree->Branch("jetCEMF", 	&jetCEMF);
  outTree->Branch("jetMUF", 	&jetMUF);
  outTree->Branch("jetNHM", 	&jetNHM);
  outTree->Branch("jetCHM", 	&jetCHM);

}

// ------------ method called once each job just after ending the event loop  ------------
void LLPgammaAnalyzer::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void LLPgammaAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(LLPgammaAnalyzer);
