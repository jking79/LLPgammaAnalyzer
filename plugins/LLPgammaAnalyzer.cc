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
LLPgammaAnalyzer::LLPgammaAnalyzer(const edm::ParameterSet& iConfig) :

  //tracksToken_(consumes<TrackCollection>(iConfig.getUntrackedParameter<edm::InputTag>("tracks")))

// -- declare tags ----------------------------------------------------------
  // triggers
  //triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
  //triggerObjectsTag(iConfig.getParameter<edm::InputTag>("triggerObjects")),

  // tracks
  tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),

  // vertices
  verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),

  // rho
  //rhoTag(iConfig.getParameter<edm::InputTag>("rho")),

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
  ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons"))

// -- end of tag declarations ---------------------------------------
{

  usesResource();
  usesResource("TFileService");

// -- consume tags ------------------------------------------------------------

  // Triggers
  //triggerResultsToken_	=consumes<edm::TriggerResults>(triggerResultsTag);
  //triggerObjectsToken_	=consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectsTag);

  // tracks 
  tracksToken_				=consumes<std::vector<reco::Track>>(tracksTag);

  // vertices
  verticesToken_			=consumes<std::vector<reco::Vertex>>(verticesTag);
  
  // rho
  //rhoToken_					=consumes<double>(rhoTag);

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
    
    int tightLepVeto = 2;
    int tight = 1;
    int loose = 0;

    bool nhf9  = NHF  < 0.90;
    bool nhf2  = NHF  > 0.02;
    bool nemf9 = NEMF < 0.90;
    bool nemf1 = NEMF < 0.99;
    bool nemf2 = NEMF > 0.02;
    bool shm1  = SHM  > 1;
    bool muf8  = MUF  < 0.80;
    bool chf0  = CHF  > 0;
    bool chm0  = CHM  > 0;
    bool cemf8 = CEMF > 0.80;
    bool nhm2  = NHM  > 2;
    bool nhm10 = NHM  > 10;

    bool eta24 = eta <= 2.4;
    bool eta27 = eta <= 2.7;
    bool eta30 = eta <= 3.0;

    if (eta24){

      	if      (nhf9 && nemf9 && shm1 && muf8 && chf0 && chm0 && cemf8) return tightLepVeto;
      	else if (nhf9 && nemf9 && shm1 && chf0 && chm0) return tight;
      	else    return loose;

    } else if (!eta24 && eta27 ){

      	if      (nhf9 && nemf9 && shm1 && muf8) return tightLepVeto;
      	else if (nhf9 && nemf9 && shm1) return tight;
      	else    return loose; 

    } else if (!eta27 && eta30){

      	if      (nemf2 && nemf1 && nhm2) return tight;
      	else    return loose; 

    } else {

      	if      (nemf9 && nhf2 && nhm10) return tight;
      	else    return loose; 

    }

    return -1; // should not happen
}



// ------------ method called for each event  ------------
void LLPgammaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

// -- Consume Tokens --------------------------------------------
  // TRIGGER
  //iEvent.getByToken(triggerResultsToken_,triggerResults_);
  //iEvent.getByToken(triggerObjectsToken_,triggerObjects_);


  // TRACKS
  //iEvent.getByToken(tracksToken_,tracks_);

  // VERTICES
  //iEvent.getByToken(verticesToken_,vertices_);

  // RHO
  //iEvent.getByToken(rhoToken_,rho_);

  // METS
  //iEvent.getByToken(metsToken_,mets_);

  // JETS
  iEvent.getByToken(jetsToken_,jets_);

  // LEPTONS
  //iEvent.getByToken(electronsToken_,electrons_);
  //iEvent.getByToken(muonsToken_,muons_);

  // ECAL RECHITS
  //iEvent.getByToken(recHitsEBToken_,recHitsEB_);
  //iEvent.getByToken(recHitsEEToken_,recHitsEE_);

// -- Process Objects ------------------------------------------

// ---- process jets --------------------------

// ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
// ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 

   jets.clear(); 
   jets.reserve(jets_->size());

   auto jetpTmin  = 100.0;
   auto jetIDmin  = 0;
   auto jetEtamax = 2.4;

   for(const auto& jet : *jets_ ){ // Filters jet collection & sorts by pt

      if (jet.pt() < jetpTmin) continue;
      if (std::abs(jet.eta()) > jetEtamax) continue;
      
      //const auto jetID = GetPFJetID(jet);
      const auto jetID = 1;
      if (jetID < jetIDmin) continue;

      // save the jets, and then store the ID
      jets.emplace_back(jet);
      jets.back().addUserInt("jetID",jetID);
      
      std::sort(jets.begin(),jets.end(),sortByPt);
   }

   nJets = jets.size();
   // set the number of leading jets to skim ( = nJets for all )
   auto stJets = nJets; 

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

      jetE.push_back(jet.energy());
      jetpt.push_back(jet.pt());
      jetphi.push_back(jet.phi());
      jeteta.push_back(jet.eta());
      jetID.push_back(jet.userInt("jetID"));
      jetNHF.push_back(jet.neutralHadronEnergyFraction());
      jetNEMF.push_back(jet.neutralEmEnergyFraction());
      jetCHF.push_back(jet.chargedHadronEnergyFraction());
      jetCEMF.push_back(jet.chargedEmEnergyFraction());
      jetMUF.push_back(jet.muonEnergyFraction());
      jetNHM.push_back(jet.neutralMultiplicity());
      jetCHM.push_back(jet.chargedMultiplicity());
 
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
