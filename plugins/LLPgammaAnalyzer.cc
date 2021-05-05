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

  // pfcands
  pfcandTag(iConfig.getParameter<edm::InputTag>("pfcandidates")),

  // vertices
  verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),

  // rho
  //rhoTag(iConfig.getParameter<edm::InputTag>("rho")),

  // mets
  metsTag(iConfig.getParameter<edm::InputTag>("mets")),  

  // supercluster
  superClusterCollectionEBTag(iConfig.getParameter<edm::InputTag>("superClusters")),
  superClusterCollectionEETag(iConfig.getParameter<edm::InputTag>("ootSuperClusters")),

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

  // pfcandidates
  pfcand_token_         =consumes<CandidateView>(pfcandTag);

  // vertices
  verticesToken_			=consumes<std::vector<reco::Vertex>>(verticesTag);

  // rho
  //rhoToken_					=consumes<double>(rhoTag);

  // mets
  metsToken_				=consumes<std::vector<pat::MET>>(metsTag);

  // supercluster
  ebScToken_            =consumes<reco::SuperClusterCollection>(superClusterCollectionEBTag);
  eeScToken_            =consumes<reco::SuperClusterCollection>(superClusterCollectionEETag); 

  // jets
  jetsToken_				=consumes<std::vector<pat::Jet>>(jetsTag);
  
  // leptons
  electronsToken_			=consumes<std::vector<pat::Electron>>(electronsTag);
  muonsToken_				=consumes<std::vector<pat::Muon>>(muonsTag);

  // rechits
  recHitsEBToken_			=consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
  recHitsEEToken_			=consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

  // photons
  gedPhotonsToken_ 		= consumes<std::vector<pat::Photon>>(gedPhotonsTag);
  ootPhotonsToken_ 		= consumes<std::vector<pat::Photon>>(ootPhotonsTag);

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

int LLPgammaAnalyzer::GetPFJetID(const pat::Jet & jet){ // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017

    const auto eta = std::abs(jet.eta());
    
    const auto NHF  = jet.neutralHadronEnergyFraction();
    const auto NEMF = jet.neutralEmEnergyFraction();
    const auto CHF  = jet.chargedHadronEnergyFraction();
    const auto CEMF = jet.chargedEmEnergyFraction();
    const auto NHM  = jet.neutralMultiplicity();
    const auto CHM  = jet.chargedMultiplicity();
    const auto SHM  = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    const auto MUF  = jet.muonEnergyFraction();

    int tighter = 3;    
    int tightLepVeto = 0;
    int tight = 2;
    int loose = 1;


    bool nhfup  = NHF  < 0.90;
    bool nhflw  = NHF  > 0.2;

    bool nemfup1 = NEMF < 0.90;
    bool nemfup2 = NEMF < 0.99;
    bool nemf80 = NEMF < 0.80;
    bool nemflw = NEMF > 0.01;
    bool nemf10 = NEMF > 0.1;

    bool shm1  = SHM  > 1;
    bool muf8  = MUF  < 0.80;
    bool chf0  = CHF  > 0;
    bool chf10  = CHF  > 0.10;
    bool chm0  = CHM  > 0;
    bool cemf8 = CEMF > 0.80;
    bool nhm2  = NHM  > 1;
    bool nhm10 = NHM  > 10;


    bool eta1 = eta <= 2.6;
    bool eta2 = eta <= 2.7;
    bool eta3 = eta <= 3.0;

    if (eta1){

      	if      (nhfup && nemfup1 && shm1 && muf8 && chf0 && chm0 && cemf8) return tightLepVeto;
         else if (nhfup && nemf80 && shm1 && chf10 && chm0) return tighter;
      	else if (nhfup && nemfup1 && shm1 && chf0 && chm0) return tight;
      	else    return loose;

    } else if (!eta1 && eta2 ){

      	if      (nhfup && nemfup2 && chm0 && muf8 && cemf8) return tightLepVeto;
         else if (nhfup && nemf80 && chm0) return tighter;
      	else if (nhfup && nemfup2 && chm0) return tight;
      	else    return loose; 

    } else if (!eta2 && eta3){

         if      (nemf10 && nemf80 && nhm2) return tighter;
      	else if (nemflw && nemfup2 && nhm2) return tight;
      	else    return loose; 

    } else {

      	if      (nhflw && nemfup1 && nhm10) return tight;
      	else    return loose; 

    }

    return -1; // should not happen
}

vector<float>  LLPgammaAnalyzer::GetRecHitdRMatchedTime( const recHitCol * rheb, const recHitCol * rhee, float eta, float phi, float drmin ){

	 vector<float> result;

            for (const auto recHit : *rheb ){

               const auto recHitId(recHit.detid());
               const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
               const auto rhPosX = recHitPos.x();
               const auto rhPosY = recHitPos.y();
               const auto rhPosZ = recHitPos.z();
               const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
               const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
               const auto tof = (d_rh-d_pv)/sol;
               const auto dr = sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
               if( dr > drmin ) continue;
               auto rht = recHit.time();
               if( rht == 0.0 ) continue;
               result.push_back(rht-tof);
            }

            for (const auto recHit : *rhee ){

               const auto recHitId(recHit.detid());
               const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
               const auto rhPosX = recHitPos.x();
               const auto rhPosY = recHitPos.y();
               const auto rhPosZ = recHitPos.z();
               const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
               const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
               const auto tof = (d_rh-d_pv)/sol;
               const auto dr = sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
               if( dr > drmin ) continue;
               auto rht = recHit.time();
               if( rht == 0.0 ) continue;
               result.push_back(rht-tof);
            }
    
	 return result;
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
  iEvent.getByToken(tracksToken_,tracks_);

  // PFCANDIDATES
  iEvent.getByToken(pfcand_token_, pfcands_);

  // VERTICES
  iEvent.getByToken(verticesToken_,vertices_);

  // RHO
  //iEvent.getByToken(rhoToken_,rho_);

  // METS
  //iEvent.getByToken(metsToken_,mets_);

  // SUPERCLUSTERS
  iEvent.getByToken(ebScToken_, superClusterEB_);  
  iEvent.getByToken(eeScToken_, superClusterEE_);

  // JETS
  iEvent.getByToken(jetsToken_,jets_);

  // LEPTONS & PHOTONS
  iEvent.getByToken(electronsToken_,electrons_);
  iEvent.getByToken(muonsToken_,muons_);

  // PHOTONS
  iEvent.getByToken(gedPhotonsToken_,gedPhotons_);
  iEvent.getByToken(ootPhotonsToken_,ootPhotons_);

  // ECAL RECHITS
  iEvent.getByToken(recHitsEBToken_,recHitsEB_);
  iEvent.getByToken(recHitsEEToken_,recHitsEE_);

  // GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
  iSetup.get<CaloGeometryRecord>().get(caloGeo_); 
  barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
  endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap); 

// -- Process Objects ------------------------------------------

// -- Process Prime Vertix
  const auto & primevtx = vertices_->front();

  auto vtxX = primevtx.position().x();
  auto vtxY = primevtx.position().y();
  auto vtxZ = primevtx.position().z();


// ---- process jets --------------------------

// ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
// ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 

   std::vector<pat::Jet> fjets;
   //jets.clear(); 
   //jets.reserve(jets_->size());

   // deltaRmin from ecalelf SC Associtor : 0.2
   auto deltaRminJet  = 0.4;
   auto deltaRminKid  = 0.2;
   auto jetPTmin  = 200.0;
   auto jetIDmin  = 3; //3;
   auto jetETAmax = 1.5; //1.5;
   unsigned int minRHcnt = 32; //32;
   auto minRHenr = 0;

   //std::cout << "Filter Jets" << std::endl;
   for(const auto& jet : *jets_ ){ // Filters jet collection & sorts by pt

      if (jet.pt() < jetPTmin) continue;
      if (std::abs(jet.eta()) > jetETAmax) continue;
      
      const auto jetID = GetPFJetID(jet);
      if (jetID < jetIDmin) continue;

      // save the jets, and then store the ID
      fjets.emplace_back(jet);
      fjets.back().addUserInt("jetID",jetID);
      
      std::sort(fjets.begin(),fjets.end(),sortByPt);
   }

   nJets = fjets.size();
   // set the number of leading jets to skim ( = nJets for all )
   //auto stJets = nJets; 
   //std::cout << "Init for Jet Loop with " << nJets << " jets"<< std::endl;

   jetE.clear();
   jetPt.clear(); 
   jetPhi.clear(); 
   jetEta.clear(); 
   jetID	.clear();
   jetNHF.clear();
   jetNEMF.clear();	
   jetCHF.clear();
   jetCEMF.clear();
   jetMUF.clear();
   jetNHM.clear();
   jetCHM.clear();
   jetTime.clear();
   jetTimeError.clear();
   jetTimeRMS.clear();
   jetMedTime.clear();
	jetPHM.clear();
	jetELM.clear();
	jetC.clear();
	jetPHE.clear();
	jetPHEF.clear();
	jetELE.clear();
	jetELEF.clear();
	jetMUE.clear();
   jetCharge.clear();

   njetRecHits.clear();
   jetRecHitOfJet.clear();
   jetRecHitId.clear();

   njetKids.clear();
   jetKidOfJet.clear();
   jetKidE.clear();
   jetKidPt.clear();
   jetKidPhi.clear();
   jetKidEta.clear();
	jetKidPdgID.clear();
	jetKidCharge.clear();
	jetKid3Charge.clear();
	jetKidLLP.clear();
	jetKidMass.clear();
	jetKidVx.clear();
	jetKidVy.clear();
	jetKidVz.clear();

   njetSubs.clear();

   auto jetHt(0.0f);

//const reco::TrackRefVector& associatedTracks() const;
//const reco::PixelClusterTagInfo* tagInfoPixelCluster(const std::string& label = "") const; // from b tagging info methods
//reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }

   hist1d18->Fill(nJets);
   //std::cout << "Starting Jet Loop" << std::endl; 
   for (auto i = 0; i < nJets; i++){ // places jet info in output tree

      const auto & jet = fjets[i];

      // jetID in jet.h ?

      jetHt += jet.pt();

      jetE.push_back(jet.energy());
      jetPt.push_back(jet.pt());
      jetPhi.push_back(jet.phi());
      jetEta.push_back(jet.eta());
      jetID.push_back(jet.userInt("jetID"));
      jetNHF.push_back(jet.neutralHadronEnergyFraction());
      jetNEMF.push_back(jet.neutralEmEnergyFraction());
      jetCHF.push_back(jet.chargedHadronEnergyFraction());
      jetCEMF.push_back(jet.chargedEmEnergyFraction());
      jetMUF.push_back(jet.muonEnergyFraction());
      jetNHM.push_back(jet.neutralMultiplicity());
      jetCHM.push_back(jet.chargedMultiplicity());
		jetCharge.push_back(jet.jetCharge());

		// const reco::TrackRefVector& associatedTracks() const;
		jetPHE.push_back(jet.photonEnergy());
		jetPHEF.push_back(jet.photonEnergyFraction()); 
		jetELE.push_back(jet.electronEnergy());
		jetELEF.push_back(jet.electronEnergyFraction());
		jetMUE.push_back(jet.muonEnergy());
		jetPHM.push_back(jet.photonMultiplicity());
		jetELM.push_back(jet.electronMultiplicity());
		// reco::JetID const& jetID() 
		// 0509     CaloTowerFwdPtrVector const& caloTowersFwdPtr() const { return caloTowersFwdPtr_; }
		// 0510     reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }
		// 0511     edm::FwdRef<reco::GenJetCollection> const& genJetFwdRef() const { return genJetFwdRef_; }
		// 0512     TagInfoFwdPtrCollection const& tagInfosFwdPtr() const { return tagInfosFwdPtr_; }

      //std::cout << "Fill jet pt/phi/eta Histograms" << std::endl;
 
      hist1d12->Fill(jet.pt());
      hist1d13->Fill(jet.phi());
      hist1d14->Fill(jet.eta());

      //std::cout << "Initing clusterTimes" << std::endl;
      TH1D clusterTimes( "temp", "temp", 500, -25, 25 );
		float sc_eta(0.0);
      float sc_phi(0.0);
      float sc_enr(0.0);

      //std::cout << "Starting RecHit Loop" << std::endl;
 
      unsigned int rhCount = 0;
	   for (const auto recHit : *recHitsEB_ ){

         const auto recHitId(recHit.detid());
         const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
			const auto rhPosX = recHitPos.x();
         const auto rhPosY = recHitPos.y();
         const auto rhPosZ = recHitPos.z();
         const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
         const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
			const auto tof = (d_rh-d_pv)/sol;
			const auto dr = sqrt(reco::deltaR2(jet.eta(), jet.phi(), recHitPos.eta(), recHitPos.phi()));
		   if( dr > deltaRminJet ) continue;
         auto rht = recHit.time();
         if( rht == 0.0 ) continue;
         auto enr = recHit.energy();
         if( enr < minRHenr ) continue;
         if( enr > sc_enr ){ sc_enr = enr; sc_eta = recHitPos.eta(); sc_phi = recHitPos.phi(); }
			rhCount++;
		   jetRecHitOfJet.push_back(i);
         jetRecHitId.push_back(recHitId);	
         clusterTimes.Fill(rht-tof);
         jetRHTimeHist->Fill(rht-tof);
         hist2d38->Fill(rht-tof,enr);
      } 

      for (const auto recHit : *recHitsEE_ ){

         const auto recHitId(recHit.detid());
         const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
         const auto rhPosX = recHitPos.x();
         const auto rhPosY = recHitPos.y();
         const auto rhPosZ = recHitPos.z();
         const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
         const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
         const auto tof = (d_rh-d_pv)/sol;
         const auto dr = sqrt(reco::deltaR2(jet.eta(), jet.phi(), recHitPos.eta(), recHitPos.phi()));
         if( dr > deltaRminJet ) continue;
         auto rht = recHit.time();
         if( rht == 0.0 ) continue;
         auto enr = recHit.energy();
         if( enr < minRHenr ) continue;
         if( enr > sc_enr ){ sc_enr = enr; sc_eta = recHitPos.eta(); sc_phi = recHitPos.phi(); }
         rhCount++;
         jetRecHitOfJet.push_back(i);
         jetRecHitId.push_back(recHitId);
         clusterTimes.Fill(rht-tof);
         jetRHTimeHist->Fill(rht-tof);
         hist2d38->Fill(rht-tof,enr);
      }

	   //  make jettime varible
	   auto jtime = clusterTimes.GetMean();
      auto jterr = clusterTimes.GetMeanError();
      auto jtrms = clusterTimes.GetRMS();
      //std::cout << "Mean Jet Time : " << jtime << std::endl;
      double medtime(0.0), quant(0.5); // 0.5 for "median"
      if( jtime == 0.0 ){
			if( rhCount == 0 ){
            jtime = -9.5;
            jterr = -9.5;
            jtrms = -9.5;
            medtime = -9.5;
			} else { 
         	jtime = -9.0;
         	jterr = -9.0;
         	jtrms = -9.0;
         	medtime = -9.0;
			}
		} else {
         clusterTimes.ComputeIntegral(); // just a precaution
         clusterTimes.GetQuantiles(1, &medtime, &quant);
		}

      if( rhCount < minRHcnt ) continue;     
      njetRecHits.push_back(rhCount);
      jetTime.push_back(jtime);
      jetTimeError.push_back(jterr);
      jetTimeRMS.push_back(jtrms);
      jetMedTime.push_back(medtime);

      jetTimeHist->Fill(jtime);
      hist1d1->Fill(rhCount);
      hist1d2->Fill(jterr);
      hist1d3->Fill(jtrms);
      hist1d4->Fill(medtime);

      //std::cout << "Filling 2D Histos" << std::endl;

      hist2d1->Fill(jtime,jet.pt());
      hist2d2->Fill(jtime,jet.userInt("jetID"));
      hist2d3->Fill(jtime,jet.neutralHadronEnergyFraction());
      hist2d4->Fill(jtime,jet.chargedHadronEnergyFraction());
      hist2d5->Fill(jtime,jet.neutralEmEnergyFraction());
      hist2d6->Fill(jtime,jet.chargedEmEnergyFraction());
      hist2d7->Fill(jtime,jet.muonEnergyFraction());
      hist2d8->Fill(jtime,jet.neutralMultiplicity());
      hist2d9->Fill(jtime,jet.chargedMultiplicity());

      hist2d10->Fill(jtime,medtime);
      hist2d24->Fill(jtime,rhCount);
      hist2d25->Fill(medtime,rhCount);
      hist2d11->Fill(jtime,jtrms);
      hist2d12->Fill(jtime,jterr);

      hist2d32->Fill(jtime,sc_eta);
      hist2d33->Fill(jtime,sc_phi);
      hist2d34->Fill(jtime,sc_enr);
      hist2d35->Fill(medtime,sc_eta);
      hist2d36->Fill(medtime,sc_phi);
      hist2d37->Fill(medtime,sc_enr);


      hist2d13->Fill(medtime,jet.pt());
      hist2d14->Fill(medtime,jet.userInt("jetID"));
      hist2d15->Fill(medtime,jet.neutralHadronEnergyFraction());
      hist2d16->Fill(medtime,jet.chargedHadronEnergyFraction());
      hist2d17->Fill(medtime,jet.neutralEmEnergyFraction());
      hist2d18->Fill(medtime,jet.chargedEmEnergyFraction());
      hist2d19->Fill(medtime,jet.muonEnergyFraction());
      hist2d20->Fill(medtime,jet.neutralMultiplicity());
      hist2d21->Fill(medtime,jet.chargedMultiplicity());

      //std::cout << "Pulling Kids Info" << std::endl;
      TH1D kidRHmTime( "temp", "temp", 500, -25, 25 );
      TH1D kidSCTime( "temp", "temp", 500, -25, 25 );
	   //int nkidsc = 0;
	   vector<int> phosmatched;
		float sumetime(0.0);
		float sumenergy(0.0);
      auto nKids = jet.numberOfDaughters();     
		njetKids.push_back(nKids);
      if( nKids > 0 ) {
         for( const auto kid : jet.daughterPtrVector() ){
            jetKidOfJet.push_back(i);
				jetKidE.push_back(kid->energy());
            jetKidPt.push_back(kid->pt());
            jetKidPhi.push_back(kid->phi());
            jetKidEta.push_back(kid->eta());
            jetKidPdgID.push_back(kid->pdgId()); //int
            jetKidLLP.push_back(kid->longLived()); // bool
            jetKidCharge.push_back(kid->charge()); // int
            jetKid3Charge.push_back(kid->threeCharge()); // int
            //boostToCM() vector
            jetKidMass.push_back(kid->mass()); // double
            jetKidVx.push_back(kid->vx()); // double
            jetKidVy.push_back(kid->vy()); // double
            jetKidVz.push_back(kid->vz()); // double

				auto mindr = 0.01;
				//auto thismatch = 0;
				auto phocnt = 0;
				//float foundtime = 0;
				//float foundenergy = 0;
				//for( const auto photon : *gedPhotons_ ){
					//phocnt++;
               //const auto dr = sqrt(reco::deltaR2(kid->eta(), kid->phi(), photon.eta(), photon.phi()));
					//if( dr < mindr ){ 
					//	mindr = dr;
					//	thismatch = phocnt;
						//foundtime = photon.superCluster()
						//foundenergy		
					//}
               //std::cout << "dr : " << kid->eta() << " " << kid->phi() << " " << photon.eta() << " " << photon.phi() << " " << dr << std::endl;
					//std::cout << "Photon energy % : " << photon.energy()/kid->energy() << std::endl;
					// if( abs( 1 - photon.energy()/kid->energy() ) < 0.01  ) ph_match_cnt++;
				//} 
				//if( mindr < 0.01 ) std::cout << " Match found : " << mindr << std::endl;
				//else std::cout << " Match Not found " << std::endl;
				// std::cout << " Number of photon matches : " << ph_match_cnt << std::endl;

            // NOPE -> auto kidsupclust = superClusterEB_->ptrAt(kid.key());
				// NOPE -> auto kidPhoton = gedPhotons_->ptrAt(kid.key());
				auto kidcand = pfcands_->ptrAt(kid.key());	
            //std::cout << "Kidcand pt is : " << kidcand->pt() << " kid pt is : " << kid->pt() << std::endl;	
				const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(kidcand.get());
            //const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(&(*kidcand));
            if( packed_cand != NULL ){
					std::cout << "Packed Cand cast worked. " << std::endl;
               //std::cout << " time : " << packed_cand->time()  << std::endl;
               std::cout << " status : " << packed_cand->status()  << std::endl;
					//std::cout << " isGoodEgamma : "	<< packed_cand->isGoodEgamma() << std::endl;
					//std::cout << " MasterClone : " << packed_cand->hasMasterClone() << std::endl;
					//std::cout << " pdgId : " << packed_cand->pdgId()  << std::endl;
					int nass = 0;
					//edm::Ptr<pat::PackedCandidate> pfCandPtr = edm::refToPtr(packed_cand);

					for( const auto photon : *gedPhotons_ ){
    					edm::RefVector<pat::PackedCandidateCollection> associated =  photon.associatedPackedPFCandidates();
    					for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
        					edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
        					if( associatedPtr.get() == packed_cand ) nass++;
    					}
					}
					std::cout << "PackedCand matches to gedPhotons found : " << nass << std::endl;

               for( const auto photon : *ootPhotons_ ){
                  edm::RefVector<pat::PackedCandidateCollection> associated =  photon.associatedPackedPFCandidates();
                  for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
                     edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
                     if( associatedPtr.get() == packed_cand ) nass++;
                  }
               }
               std::cout << "PackedCand matches to ootPhotons found : " << nass << std::endl;

               for( const auto electron : *electrons_ ){
                  edm::RefVector<pat::PackedCandidateCollection> associated =  electron.associatedPackedPFCandidates();
                  for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
                     edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
                     if( associatedPtr.get() == packed_cand ) nass++;
                  }
               }
               std::cout << "PackedCand matches to Electrons found : " << nass << std::endl;

               //for( const auto muon : *muons_ ){
               //   edm::RefVector<pat::PackedCandidateCollection> associated =  muon.embeddedPFCandidate_();
               //   for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
               //      edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
               //      if( associatedPtr.get() == packed_cand ) nass++;
               //   }
               //}
               //std::cout << "PackedCand matches to Muons found : " << nass << std::endl;

				} else std::cout << "PackedCand cast failed." << std::endl;
            //const reco::PFCandidate* pfcand = dynamic_cast<const reco::PFCandidate*>(kid);
				//const auto pfcand = dynamic_cast<const reco::PFCandidate*>(kidcand.get());
            // Nope -> const auto *pfcand = dynamic_cast<const reco::PFCandidate*>(&(*kidcand));
				//if( pfcand != NULL ){
				//	std::cout << "PFCand castworked." << std::endl;
				//	auto scref = pfcand->superClusterRef();
				//	if( scref.isNonnull() ) std::cout << "Found SC Reffrence " << std::endl;		
				//} else std::cout << "PFCand cast failed." << std::endl;
				hist1d5->Fill(kid->pdgId());
            hist1d7->Fill(kid->energy());
            hist1d8->Fill(kid->pt());
            hist1d9->Fill(kid->eta());
            hist1d10->Fill(kid->phi());
            hist1d19->Fill(kid->mass());
            hist1d20->Fill(kid->charge());
            hist1d21->Fill(kid->vx());
            hist1d22->Fill(kid->vy());
            hist1d23->Fill(kid->vz());

				if( kid->hasMasterClone() ) hist1d31->Fill(1); //std::cout << "Has master clone : " << std::endl;
            if( kid->hasMasterClonePtr() ) hist1d31->Fill(2); //std::cout << "Has master clone ptr : " << std::endl;
				if( kid->numberOfMothers() ) hist1d31->Fill(3); //std::cout << "Has : " << kid->numberOfMothers() << " mothers: " <<  std::endl;
            if( kid->numberOfDaughters() ) hist1d31->Fill(4); //std::cout << "Has : " << kid->numberOfDaughters() << " daughters: " <<  std::endl;
            if( kid->numberOfSourceCandidatePtrs() )hist1d31->Fill(5); // std::cout << "Has : " << kid->numberOfSourceCandidatePtrs() i

            if( kid->isElectron() ) hist1d30->Fill(1); //std::cout << "Is an Electron : " << std::endl;
            if( kid->isMuon() ) hist1d30->Fill(2); //std::cout << "Is a Muon : " << std::endl;
            if( kid->isStandAloneMuon() ) hist1d30->Fill(3); //std::cout << "Is a StandAloneMuon : " << std::endl;
            if( kid->isGlobalMuon() ) hist1d30->Fill(4); //std::cout << "Is a GlobalMuon : " << std::endl;
            if( kid->isTrackerMuon() ) hist1d30->Fill(5); //std::cout << "Is an TrackerMuon : " << std::endl;
            if( kid->isCaloMuon() ) hist1d30->Fill(6); //std::cout << "Is an CaloMuon : " << std::endl;
            if( kid->isPhoton() ) hist1d30->Fill(7); //std::cout << "Is an Photon : " << std::endl;
            if( kid->isConvertedPhoton() ) hist1d30->Fill(8); //std::cout << "Is an ConvertedPhoton : " << std::endl;
            if( kid->isJet() ) hist1d30->Fill(9); //std::cout << "Is an Jet : " << std::endl;


            //std::cout << "Initing clusterTimes" << std::endl;
            TH1D kidClTimes( "temp", "temp", 500, -25, 25 );
      
            //std::cout << "finding kid RecHits" << std::endl;
      
            unsigned int kidRhCount = 0;
            for (const auto recHit : *recHitsEB_ ){
      
               const auto recHitId(recHit.detid());
               const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
         		const auto rhPosX = recHitPos.x();
         		const auto rhPosY = recHitPos.y();
         		const auto rhPosZ = recHitPos.z();
         		const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
         		const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
         		const auto tof = (d_rh-d_pv)/sol;
               const auto dr = sqrt(reco::deltaR2(kid->eta(), kid->phi(), recHitPos.eta(), recHitPos.phi()));
               if( dr > deltaRminKid ) continue;
               auto rht = recHit.time();
               if( rht == 0.0 ) continue;
               kidRhCount++;
               //jetRecHitOfJet.push_back(i);
               //jetRecHitId.push_back(recHitId);
               kidClTimes.Fill(rht-tof);
            }
      
            for (const auto recHit : *recHitsEE_ ){
      
               const auto recHitId(recHit.detid());
               const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
               const auto rhPosX = recHitPos.x();
               const auto rhPosY = recHitPos.y();
               const auto rhPosZ = recHitPos.z();
               const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
               const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
               const auto tof = (d_rh-d_pv)/sol;
               const auto dr = sqrt(reco::deltaR2(kid->eta(), kid->phi(), recHitPos.eta(), recHitPos.phi()));
               if( dr > deltaRminKid ) continue;
               auto rht = recHit.time();
               if( rht == 0.0 ) continue;
               kidRhCount++;
               //jetRecHitOfJet.push_back(i);
               //jetRecHitId.push_back(recHitId);
               kidClTimes.Fill(rht-tof);
            }

            //  make kidtime varible
            //std::cout << "Finding kid rh mean time" << std::endl;
            if( kidRhCount != 0 ) { 
            	auto ktime = kidClTimes.GetMean();
            	//auto ktErr = kidClTimes.GetMeanError();
            	//auto ktRms = kidClTimes.GetRMS();
            	//std::cout << "Mean Jet Time : " << jtime << std::endl;
            	double kmedtime(0.0), quant(0.5); // 0.5 for "median"
            	if( ktime == 0.0 ){
						if( kidRhCount == 0 ){
            	   	ktime = -9.5;
            	   	//kterr = -9.5;
            	   	//ktrms = -9.5;
            	   	kmedtime = -9.5;
						} else {
                     ktime = -9.0;
                     //kterr = -9.0;
                     //ktrms = -9.0;
                     kmedtime = -9.0;
						}
            	} else {
            	   kidClTimes.ComputeIntegral(); // just a precaution
            	   kidClTimes.GetQuantiles(1, &kmedtime, &quant);
            	}
				
					//jetKidTime.push_back(ktime);
            	//jetKidMedTime.push_back(kmedtime);
					hist1d27->Fill(ktime);
					hist1d28->Fill(kmedtime);
            	hist1d29->Fill(kidRhCount);
				   kidRHmTime.Fill(ktime);	
				} // eof kid time from rechit collection

				// get time from SC collction
            //std::cout << "Finding kid SC" << std::endl;
            const reco::SuperCluster* kidSC(NULL);
            //auto kidSCDeltaRMin = 0.8; //deltaRminKid;
            //for (const auto sceb : *superClusterEB_ ){ 
            //for(reco::SuperClusterCollection::const_iterator scIt = superClusterEB_->begin(); scIt != superClusterEB_->end(); scIt++ ) {
            //   const auto dr = sqrt(reco::deltaR2(kid->eta(), kid->phi(), scIt->eta(), scIt->phi()));
				//	std::cout << "dr : " << kid->eta() << " " << kid->phi() << " " << scIt->eta() << " " << scIt->phi() << " " << dr << std::endl;
            //   if( dr < kidSCDeltaRMin ){
            //      kidSCDeltaRMin = dr;
	         //         kidSC = &*scIt;
            //   }
            //}
            //for (const auto scee : *superClusterEE_ ){
            //   const auto dr = sqrt(reco::deltaR2(kid->eta(), kid->phi(), scee.eta(), scee.phi()));
            //   if( dr < kidSCDeltaRMin ){
            //      kidSCDeltaRMin = dr;
            //      kidSC = &scee;
            //   }
            //}

				if( kidSC != NULL ){ 	
	            //std::cout << "Finding kid SC seed crystal" << std::endl;
					const auto & seedDetId = kidSC->seed()->seed(); // seed detid
					const auto seedRawId = seedDetId.rawId(); // crystal number
	            //std::cout << "Running Rh for kid SC seed crystal" << std::endl;
					for (const auto recHit : *recHitsEB_ ){
	               	const auto recHitId(recHit.detid());
							if( seedRawId != recHitId ) continue;
	               	const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
	               	const auto rhPosX = recHitPos.x();
	               	const auto rhPosY = recHitPos.y();
	               	const auto rhPosZ = recHitPos.z();
	               	const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
	               	const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
	               	const auto tof = (d_rh-d_pv)/sol;
	               	auto rht = recHit.time();
	               	if( rht == 0.0 ) continue;
					      //nkidsc++;				
	               	kidSCTime.Fill(rht-tof);
	                  hist1d32->Fill(rht-tof);
					}
	            //for (const auto recHit : *recHitsEE_ ){
	            //      const auto recHitId(recHit.detid());
	            //      if( seedRawId != recHitId ) continue;
	            //      const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
	            //      const auto rhPosX = recHitPos.x();
	            //      const auto rhPosY = recHitPos.y();
	            //      const auto rhPosZ = recHitPos.z();
	            //      const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
	            //      const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
	            //      const auto tof = (d_rh-d_pv)/sol;
	            //      auto rht = recHit.time();
	            //      if( rht == 0.0 ) continue;
	            //      nkidsc++;
	            //      kidSCTime.Fill(rht-tof);
					//		hist1d32->Fill(rht-tof);
	            //}
				}
			}// kid loop

			//std::cout << "Finding jet time from kids" << std::endl;
         auto sctime = kidSCTime.GetMean();
			if( sctime == 0.0 ) sctime = -9.0;
         hist1d33->Fill(sctime);
			auto rhmtime = kidRHmTime.GetMean();
         if( rhmtime == 0.0 ) rhmtime = -9.0;
         hist1d34->Fill(rhmtime);
		}// has kids

		hist1d6->Fill(nKids);
      // const reco::Candidate* daughter(size_t i)
      // reco::CandidatePtr daughterPtr(size_t i)
      // const reco::CompositePtrCandidate::daughters& daughterPtrVector()


      //const auto & pfcons = jet.getPFConstituents();
      //auto nSubs = pfcons.size();
      //auto nSubs = jet.nSubjetCollections();
      //njetSubs.push_back(nSubs);
      // pat::JetPtrCollection const& subjets(unsigned int index = 0)
      // pat::JetPtrCollection const& subjets(std::string const& label)
      // bool hasSubjets(std::string const& label)
      // std::vector<std::string> const& subjetCollectionNames()
      // double groomedMass(unsigned int index = 0)
      // double groomedMass(std::string const& label) 
		std::cout << "Next Jet .......................... " << std::endl; 
   }
// ** end of jets  *************************************

   hist1d17->Fill(jetHt);

   // d jetTime for back-to-back high pt jets
   //auto dijetIdCut = 1;
   auto dijetPtMin = 200.0;
   auto difPtLmt = 0.8;
   auto htPctLmt = 0.8;
   auto dPhiLmt = 2.8;

   //std::cout << "Finding jetTimes" << std::endl;
   for (auto q = 0; q < nJets; q++){
      for (auto p = q+1; p < nJets; p++){
         const auto & qjet = fjets[q];
         const auto & pjet = fjets[p];
         if( qjet.pt() < dijetPtMin ) continue;
         auto diffPt = pjet.pt()/qjet.pt();
         hist1d24->Fill(diffPt);
         if( diffPt < difPtLmt ) continue;
         auto htPct= (qjet.pt()+pjet.pt())/jetHt;
			hist1d25->Fill(htPct);
         if( htPct < htPctLmt ) continue;
         auto dPhi = reco::deltaPhi(qjet.phi(),pjet.phi());
         hist1d26->Fill(dPhi);
         if( dPhi < dPhiLmt ) continue;
			auto dTmu = jetTime[q] - jetTime[p];
         auto dTmed = jetMedTime[q] - jetMedTime[p];
         if( dTmu == 0.0 )  dTmu = -5.5;
         if( dTmed == 0.0 )  dTmed = -5.5;
         if( jetTime[q] == 0.0 || jetTime[p] == 0.0 ) dTmu = -5.0;
         if( jetMedTime[q] == 0.0 || jetMedTime[p] == 0.0 ) dTmed = -5.0;
			hist1d15->Fill(dTmu);
			hist2d22->Fill(dTmu,nJets);
         hist2d26->Fill(dTmu,diffPt);
         hist2d27->Fill(dTmu,htPct);
         hist2d28->Fill(dTmu,dPhi);
			hist1d16->Fill(dTmed);
         hist2d23->Fill(dTmed,nJets);
         hist2d29->Fill(dTmu,diffPt);
         hist2d30->Fill(dTmu,htPct);
         hist2d31->Fill(dTmu,dPhi);
		}
	}

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

  // Book output files and trees
  edm::Service<TFileService> fs;
  outTree = fs->make<TTree>("llpgtree","llpgtree");

  // Book histograms
  jetRHTimeHist = fs->make<TH1D>("jetRHTime" , "jetRHTime" , 2000 , -100 , 100 );
  jetTimeHist = fs->make<TH1D>("jetMuTime" , "jetMuTime" , 200 , -10 , 10 );
  hist1d1 = fs->make<TH1D>("jetRHMulti" , "jetRHMulti" , 250 , 0 , 250 );
  hist1d2 = fs->make<TH1D>("jetTimeError" , "jetTimeError" , 300 , 0 , 3 );
  hist1d3 = fs->make<TH1D>("jetTimeRMS" , "jetTimeRMS" , 200 , 0 , 20 );
  hist1d4 = fs->make<TH1D>("jetMedTime" , "jetMedTime" , 200 , -10 , 10 );
  hist1d5 = fs->make<TH1D>("kidPdgID" , "kidPdgID" , 300 , 0 , 300 );
  hist1d6 = fs->make<TH1D>("nKids" , "nKids" , 100 , 0 , 100 );
  hist1d7 = fs->make<TH1D>("kidE" , "kidE" , 500 , 0 , 500 );
  hist1d8 = fs->make<TH1D>("kidPt" , "kidPt" , 500, 0 , 500 );
  hist1d9 = fs->make<TH1D>("kidEta" , "kidEta" , 700 , -3.5 , 3.5 );
  hist1d10 = fs->make<TH1D>("kidPhi" , "kidPhi" , 700 , -3.5 , 3.5 );
  //hist1d11 = fs->make<TH1D>("t0jetRHMulti" , "t0jetRHMulti" , 250 , 0 , 250 );

  hist1d12 = fs->make<TH1D>("jetPt" , "jetPt" , 500, 0 , 500 );
  hist1d13 = fs->make<TH1D>("jetPhi" , "jetPhi" , 700 , -3.5 , 3.5 );
  hist1d14 = fs->make<TH1D>("jetEta" , "jetEta" , 700 , -3.5 , 3.5 );
  hist1d15 = fs->make<TH1D>("jetdtmu" , "jetdtmu" , 120 , -6 , 6 );
  hist1d16 = fs->make<TH1D>("jetdtmed" , "jetdtmed" , 120 , -6 , 6 );
  hist1d17 = fs->make<TH1D>("jetHt" , "jetHt" , 1000 , 0 , 1000 );
  hist1d18 = fs->make<TH1D>("nJet" , "nJets" , 20 , 0 , 20 );
  hist1d19 = fs->make<TH1D>("kidMass" , "kidMass" , 20, -1.0 , 1.0 );
  hist1d20 = fs->make<TH1D>("kidChrg" , "kidChrg" , 5, -2 , 2 );
  hist1d21 = fs->make<TH1D>("kidVx" , "kidVx" , 200, -10 , 10 );
  hist1d22 = fs->make<TH1D>("kidVy" , "kidVy" , 200, -10 , 10 );
  hist1d23 = fs->make<TH1D>("kidVz" , "kidVz" , 200, -10 , 10 );

  hist1d24 = fs->make<TH1D>("diffPt" , "diffPt" , 100, 0 , 1 );
  hist1d25 = fs->make<TH1D>("htPct" , "htPct" , 100, 0 , 1 );
  hist1d26 = fs->make<TH1D>("dPhi" , "dPhi" , 32, 0 , 3.2 );

  hist1d27 = fs->make<TH1D>("kidTime" , "kidTime" , 200 , -10 , 10 );
  hist1d28 = fs->make<TH1D>("kidMedTime" , "kidMedTime" , 200 , -10 , 10 );
  hist1d29 = fs->make<TH1D>("kidRhCount", "kidRhCount", 250 , 0 , 250 );
  hist1d30 = fs->make<TH1D>("kidleafptrs", "kidleafptrs", 11 , 0 , 10 );
  hist1d31 = fs->make<TH1D>("kidleaftype", "kidleaftype", 11 , 0 , 10 );
  hist1d32 = fs->make<TH1D>("kidSCTime" , "kidSCTime" , 500 , -25 , 25 );
  hist1d33 = fs->make<TH1D>("jetSCTime" , "jetSCTime" , 500 , -25 , 25 );
  hist1d34 = fs->make<TH1D>("jetRHmTime" , "jetRHmTime" , 500 , -25 , 25 );

//------------------------------------------------------------------------------------
  hist2d1 = fs->make<TH2D>("jt_pt" , "jt_pt" , 200 , -10 , 10 , 500 , 0 , 500 );
  hist2d2 = fs->make<TH2D>("jt_id" , "jt_id" , 200 , -10 , 10 , 5 , 0 , 5 );
  hist2d3 = fs->make<TH2D>("jt_nhf" , "jt_nhf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d4 = fs->make<TH2D>("jt_chf" , "jt_chf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d5 = fs->make<TH2D>("jt_nemf" , "jt_nemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d6 = fs->make<TH2D>("jt_cemf" , "jt_cemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d7 = fs->make<TH2D>("jt_muf" , "jt_muf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d8 = fs->make<TH2D>("jt_nhm" , "jt_nhm" , 200 , -10 , 10 , 40 , 0 , 40 );
  hist2d9 = fs->make<TH2D>("jt_chm" , "jt_chm" , 200 , -10 , 10 , 40 , 0 , 40 );

  hist2d10 = fs->make<TH2D>("jt_medt" , "jt_medt" , 200 , -10 , 10 , 200 , -10 , 10 );
  hist2d11 = fs->make<TH2D>("jt_rms" , "jt_rms" , 200 , -10 , 10 , 200 , 0 , 20 );
  hist2d12 = fs->make<TH2D>("jt_err" , "jt_err" , 200 , -10 , 10 , 300 , 0 , 3 );

  hist2d13 = fs->make<TH2D>("medt_pt" , "medt_pt" , 200 , -10 , 10 , 500 , 0 , 500 );
  hist2d14 = fs->make<TH2D>("medt_id" , "medt_id" , 200 , -10 , 10 , 5 , 0 , 5 );
  hist2d15 = fs->make<TH2D>("medt_nhf" , "medt_nhf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d16 = fs->make<TH2D>("medt_chf" , "medt_chf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d17 = fs->make<TH2D>("medt_nemf" , "medt_nemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d18 = fs->make<TH2D>("medt_cemf" , "medt_cemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d19 = fs->make<TH2D>("medt_muf" , "medt_muf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d20 = fs->make<TH2D>("medt_nhm" , "medt_nhm" , 200 , -10 , 10 , 40 , 0 , 40 );
  hist2d21 = fs->make<TH2D>("medt_chm" , "medt_chm" , 200 , -10 , 10 , 40 , 0 , 40 );

  hist2d22 = fs->make<TH2D>("jdtmu_nJets" , "jdtmu_nJets" , 120 , -6 , 6 , 6 , 2 , 8 );
  hist2d23 = fs->make<TH2D>("jdtmed_nJets" , "jdtmed_nJets" , 120 , -6 , 6 , 6 , 2 , 8 );

  hist2d24 = fs->make<TH2D>("jt_nrh" , "jt_nrh" , 200 , -10 , 10 , 250 , 0 , 250 );
  hist2d25 = fs->make<TH2D>("medt_nrh" , "medt_nrh" , 200 , -10 , 10 , 250 , 0 , 250 );

  hist2d26 = fs->make<TH2D>("jdtmu_diffPt" , "jdtmu_diffPt" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d27 = fs->make<TH2D>("jdtmu_htPct" , "jdtmu_htPct" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d28 = fs->make<TH2D>("jdtmu_dPhi" , "jdtmu_dPhi" , 120 , -6 , 6 , 400 , 2.8 , 3.2 );
  hist2d29 = fs->make<TH2D>("jdtmed_diffPt" , "jdtmed_diffPt" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d30 = fs->make<TH2D>("jdtmed_htPct" , "jdtmed_htPct" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d31 = fs->make<TH2D>("jdtmed_dPhi" , "jdtmed_dPhi" , 120 , -6 , 6 , 400 , 2.8 , 3.2 );

  hist2d32 = fs->make<TH2D>("jt_sceta" , "jt_sceta", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d33 = fs->make<TH2D>("jt_scphi" ,"jt_scphi", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d34 = fs->make<TH2D>("jt_scenr" , "jt_scenr", 200 , -10 , 10 , 1000 , 0 , 1000 );

  hist2d35 = fs->make<TH2D>("medt_sceta" , "medt_sceta", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d36 = fs->make<TH2D>("medt_scphi" ,"medt_scphi", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d37 = fs->make<TH2D>("medt_scenr" , "medt_scenr", 200 , -10 , 10 , 1000 , 0 , 1000 );

  hist2d38 = fs->make<TH2D>("rht_rhe" , "rht_rhe", 200 , -10 , 10 , 1000 , 0 , 1000 );

  std::cout << "Histograms Booked" << std::endl;

  // Run, Lumi, Event info
  outTree->Branch("run", &run);
  outTree->Branch("lumi", &lumi);
  outTree->Branch("event", &event, "event/l");

  // Jet info
  outTree->Branch("njets", &nJets);
  outTree->Branch("jetE", &jetE);
  outTree->Branch("jetPt", &jetPt);
  outTree->Branch("jetEta", &jetEta);
  outTree->Branch("jetPhi", &jetPhi);
  outTree->Branch("jetID", &jetID);
  outTree->Branch("jetNHF", &jetNHF);
  outTree->Branch("jetNEMF", &jetNEMF);  
  outTree->Branch("jetCHF", &jetCHF);
  outTree->Branch("jetCEMF", &jetCEMF);
  outTree->Branch("jetMUF", &jetMUF);
  outTree->Branch("jetNHM", &jetNHM);
  outTree->Branch("jetCHM", &jetCHM);
  outTree->Branch("jetPHM", &jetPHM);
  outTree->Branch("jetELM", &jetELM);
  outTree->Branch("jetC", &jetC);
  outTree->Branch("jetPHE", &jetPHE);
  outTree->Branch("jetPHEF", &jetPHEF);
  outTree->Branch("jetELE", &jetELE);
  outTree->Branch("jetELEF", &jetELEF);
  outTree->Branch("jetMUE", &jetMUE);
  outTree->Branch("jetTime", &jetTime);
  outTree->Branch("jetTimeError", &jetTimeError);
  outTree->Branch("jetTimeRMS", &jetTimeRMS);
  outTree->Branch("jetMedTime", &jetMedTime);

  outTree->Branch("njetRecHits", &njetRecHits);
  outTree->Branch("jetRecHitOfJet", &jetRecHitOfJet);
  outTree->Branch("jetRecHitId", &jetRecHitId);

  outTree->Branch("njetKids", &njetKids);
  outTree->Branch("jetKidOfJet", &jetKidOfJet);
  outTree->Branch("jetKidE", &jetKidE);
  outTree->Branch("jetKidPt", &jetKidPt);
  outTree->Branch("jetKidPhi", &jetKidPhi);
  outTree->Branch("jetKidEta", &jetKidEta);
  outTree->Branch("jetKidPdgID", &jetKidPdgID);
  outTree->Branch("jetKidCharge", &jetKidCharge);
  outTree->Branch("jetKid3Charge", &jetKid3Charge);
  outTree->Branch("jetKidMass", &jetKidMass);
  outTree->Branch("jetKidVx", &jetKidVx);
  outTree->Branch("jetKidVy", &jetKidVy);
  outTree->Branch("jetKidVz", &jetKidVz);
  outTree->Branch("jetKidTime", &jetKidTime);
  outTree->Branch("jetKidMedTime", &jetKidMedTime);

  outTree->Branch("njetSubs", &njetSubs);

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
