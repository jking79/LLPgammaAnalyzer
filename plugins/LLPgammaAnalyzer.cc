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
  superClusterCollectionTag(iConfig.getParameter<edm::InputTag>("superClusters")),
  ootSuperClusterCollectionTag(iConfig.getParameter<edm::InputTag>("ootSuperClusters")),

  // caloclusters
  caloClusterTag(iConfig.getParameter<edm::InputTag>("caloClusters")),

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
  scToken_              =consumes<reco::SuperClusterCollection>(superClusterCollectionTag);
  ootScToken_           =consumes<reco::SuperClusterCollection>(ootSuperClusterCollectionTag); 

  // caloClusters
  ccToken_			      =consumes<std::vector<reco::CaloCluster>>(caloClusterTag);	

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

int LLPgammaAnalyzer::getPFJetID(const pat::Jet & jet){ // https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVRun2017

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

///home/t3-ku/jaking/llpgamma/CMSSW_10_6_20/src/LLPgamma/LLPgammaAnalyzer/plugins/LLPgammaAnalyzer.cc:211:9: note:   
//no known conversion for argument 1 from 'edm::Handle<edm::SortedCollection<EcalRecHit> >' to 'const recHitCol* {aka const edm::SortedCollection<EcalRecHit>*}'

rhGroup LLPgammaAnalyzer::getRHGroup( const recHitCol rheb, const recHitCol rhee, float eta, float phi, float drmin, float minenr = 0.0 ){

	 rhGroup result;

    for (const auto recHit : rheb ){

       auto enr = recHit.energy();
       if( enr <= minenr ) continue;
       const auto recHitId(recHit.detid());
       const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
       const auto dr = sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
       if( dr > drmin ) continue;
       result.push_back(recHit);
    }       

    for (const auto recHit : rhee ){

       auto enr = recHit.energy();
       if( enr <= minenr ) continue;
       const auto recHitId(recHit.detid());
       const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
       const auto dr = sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
       if( dr > drmin ) continue;
       result.push_back(recHit);
    }

    return result;
}

rhGroup LLPgammaAnalyzer::getRHGroup( const recHitCol rheb, const recHitCol rhee, DetId detid ){

    rhGroup result;

// .rawId() change needed

    for (const auto recHit : rheb ){

       const auto recHitId(recHit.detid());
       if( detid != recHitId ) continue;
       result.push_back(recHit);
    }

    for (const auto recHit : rhee ){

       const auto recHitId(recHit.detid());
       if( detid != recHitId ) continue;
       result.push_back(recHit);
    }

    return result;
}

rhGroup LLPgammaAnalyzer::getRHGroup( const recHitCol rheb, const recHitCol rhee ){

    rhGroup result;

    for (const auto recHit : rheb ){

       const auto recHitId(recHit.detid());
       result.push_back(recHit);
    }        
       
    for (const auto recHit : rhee ){
       
       const auto recHitId(recHit.detid());
       result.push_back(recHit);
    }

    return result;
}
	
vector<float>  LLPgammaAnalyzer::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    vector<float> result;

    for (const auto recHit : recHits ){

       const auto rht = recHit.time();
		 //std::cout << " ----- Get TOF Time rh time: " << rht << std::endl;
       //auto enr = recHit.energy();
		 //std::cout << " ----- Get TOF Time rh energy: " << enr << std::endl;
       const auto recHitId(recHit.detid());
       const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
       const auto rhPosX = recHitPos.x();
       const auto rhPosY = recHitPos.y();
       const auto rhPosZ = recHitPos.z();
		 //std::cout << " ----- Get TOF Time rh POS: " << rhPosX << " " <<  rhPosY << " " << rhPosZ << std::endl;
       const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
       const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
       const auto tof = (d_rh-d_pv)/sol;
       //std::cout << " ----- Get TOF Time rh tof: " << tof << std::endl;
       result.push_back(rht-tof);
    }

    return result;
}

EcalRecHit LLPgammaAnalyzer::getLeadRh( rhGroup recHits ){

    EcalRecHit result;
	 float enr(0.0);

    for (const auto recHit : recHits ){

		auto rhenr = recHit.energy();
		if( rhenr < enr ) continue;
		enr = rhenr;
		result = recHit;
	}

	return result;
}

vector<float>  LLPgammaAnalyzer::getRecHitdRMatchedTime( const recHitCol * rheb, const recHitCol * rhee, float eta, float phi, float drmin ){

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
  iEvent.getByToken(scToken_, superCluster_);  
  iEvent.getByToken(ootScToken_, ootSuperCluster_);

  // CALOCLUSTERS
  iEvent.getByToken(ccToken_, caloCluster_);

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
   auto deltaRminJet  = 0.3;//0.4
   auto deltaRminKid  = 0.15;//0.2
   auto jetPTmin  = 200.0;
   auto jetIDmin  = 3; //3;
   auto jetETAmax = 1.5; //1.5;
   unsigned int minRHcnt = 32; //32;
   auto minRHenr = 0;

   //std::cout << "Filter Jets" << std::endl;
   for(const auto& jet : *jets_ ){ // Filters jet collection & sorts by pt

      if (jet.pt() < jetPTmin) continue;
      if (std::abs(jet.eta()) > jetETAmax) continue;
      
      const auto jetID = getPFJetID(jet);
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

   hist1d[18]->Fill(nJets);
   //std::cout << "Starting Jet Loop for " << nJets << " jets " << std::endl; 
   for (auto ijet = 0; ijet < nJets; ijet++){ // places jet info in output tree

      const auto & jet = fjets[ijet];

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
 
      hist1d[12]->Fill(jet.pt());
      hist1d[13]->Fill(jet.phi());
      hist1d[14]->Fill(jet.eta());

      //std::cout << "Initing clusterTimes" << std::endl;
      TH1D clusterTimes( "temp", "temp", 500, -25, 25 );

      //std::cout << "Getting jet dR rechit group" << std::endl; 
		auto jetDrRhGroup =  getRHGroup( *recHitsEB_, *recHitsEE_, jet.eta(), jet.phi(), deltaRminJet, minRHenr ); 
		auto rhCount( jetDrRhGroup.size() );
      //std::cout << "rhCount is " << rhCount << std::endl;
		if( rhCount > minRHcnt ){
      	//std::cout << "Getting jet dr rh tof time" << std::endl;
      	auto tofTimes = getRhTofTime( jetDrRhGroup, vtxX, vtxY, vtxZ );
      	//std::cout << "Getting lead jet dR rechit" << std::endl;
      	auto leadJetRh = getLeadRh( jetDrRhGroup );

      	//std::cout << "Starting RecHit Loop" << std::endl;
		   for ( unsigned int irhg = 0; irhg < rhCount; irhg++){
	
				//std::cout << " -- irhg: " << irhg << " rhCount: " << rhCount << std::endl;
			   jetRecHitOfJet.push_back(ijet);
				auto detid = (jetDrRhGroup[irhg]).detid();
				//std::cout << " -- (jetDrRhGroup[irhg]).detid(): " << detid.rawId() << std::endl;
	         jetRecHitId.push_back(detid.rawId());	
				auto rhtime = tofTimes[irhg];
				//std::cout << " -- tofTimes[irhg]: " << rhtime << std::endl;
	         clusterTimes.Fill(rhtime);
	         jetRHTimeHist->Fill(rhtime);
				auto rhe = (jetDrRhGroup[irhg]).energy();
	         //std::cout << " -- jetDrRhGroup[irhg]).energy(): " << rhe << std::endl;
	         hist2d[38]->Fill(rhtime, rhe);
	      } 
	
	      const auto leadJetRhId(leadJetRh.detid());
	      const auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
	      auto sc_eta(leadJetRhIdPos.eta());
	      auto sc_phi(leadJetRhIdPos.phi());
	      auto sc_enr(leadJetRh.energy());
	      //std::cout << "Lead Jet dR RH Group E: " << sc_enr << " eta: " << sc_eta << " phi: " << sc_phi << std::endl;

//  **************   need default values for 0 rhCount case ? *****************************************
	
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
			}//<<>>if( jtime == 0.0 )
	
	      njetRecHits.push_back(rhCount);
	      jetTime.push_back(jtime);
	      jetTimeError.push_back(jterr);
	      jetTimeRMS.push_back(jtrms);
	      jetMedTime.push_back(medtime);
	
	      jetTimeHist->Fill(jtime);
	      hist1d[1]->Fill(rhCount);
	      hist1d[2]->Fill(jterr);
	      hist1d[3]->Fill(jtrms);
	      hist1d[4]->Fill(medtime);
	
	      //std::cout << "Filling 2D Histos" << std::endl;
	
	      hist2d[1]->Fill(jtime,jet.pt());
	      hist2d[2]->Fill(jtime,jet.userInt("jetID"));
	      hist2d[3]->Fill(jtime,jet.neutralHadronEnergyFraction());
	      hist2d[4]->Fill(jtime,jet.chargedHadronEnergyFraction());
	      hist2d[5]->Fill(jtime,jet.neutralEmEnergyFraction());
	      hist2d[6]->Fill(jtime,jet.chargedEmEnergyFraction());
	      hist2d[7]->Fill(jtime,jet.muonEnergyFraction());
	      hist2d[8]->Fill(jtime,jet.neutralMultiplicity());
	      hist2d[9]->Fill(jtime,jet.chargedMultiplicity());
	
	      hist2d[10]->Fill(jtime,medtime);
	      hist2d[24]->Fill(jtime,rhCount);
	      hist2d[25]->Fill(medtime,rhCount);
	      hist2d[11]->Fill(jtime,jtrms);
	      hist2d[12]->Fill(jtime,jterr);
	
	      hist2d[32]->Fill(jtime,sc_eta);
	      hist2d[33]->Fill(jtime,sc_phi);
	      hist2d[34]->Fill(jtime,sc_enr);
	      hist2d[35]->Fill(medtime,sc_eta);
	      hist2d[36]->Fill(medtime,sc_phi);
	      hist2d[37]->Fill(medtime,sc_enr);
	
	
	      hist2d[13]->Fill(medtime,jet.pt());
	      hist2d[14]->Fill(medtime,jet.userInt("jetID"));
	      hist2d[15]->Fill(medtime,jet.neutralHadronEnergyFraction());
	      hist2d[16]->Fill(medtime,jet.chargedHadronEnergyFraction());
	      hist2d[17]->Fill(medtime,jet.neutralEmEnergyFraction());
	      hist2d[18]->Fill(medtime,jet.chargedEmEnergyFraction());
	      hist2d[19]->Fill(medtime,jet.muonEnergyFraction());
	      hist2d[20]->Fill(medtime,jet.neutralMultiplicity());
	      hist2d[21]->Fill(medtime,jet.chargedMultiplicity());

      }//<<>>if( rhCount > minRHcnt)

      //std::cout << "Pulling Kids Info" << std::endl;
		//---------------------------------------------------------------------
      TH1D kidRHmTime( "temp", "temp", 500, -25, 25 );
      TH1D kidSCTime( "temp", "temp", 500, -25, 25 );
	   //int nkidsc = 0;
	   //vector<int> phosmatched;
		//float sumtime(0.0);
		int sum_nRechits(0);
		float sum_energy(0.0);
      auto nKids = jet.numberOfDaughters();     
		njetKids.push_back(nKids);
      if( nKids > 0 ) {
         for( const auto kid : jet.daughterPtrVector() ){
            jetKidOfJet.push_back(ijet);
				jetKidE.push_back(kid->energy());
            jetKidPt.push_back(kid->pt());
            jetKidPhi.push_back(kid->phi());
            jetKidEta.push_back(kid->eta());
            jetKidPdgID.push_back(kid->pdgId()); //int
            jetKidLLP.push_back(kid->longLived()); // bool
            jetKidCharge.push_back(kid->charge()); // int
            jetKid3Charge.push_back(kid->threeCharge()); // int
            jetKidMass.push_back(kid->mass()); // double
            jetKidVx.push_back(kid->vx()); // double
            jetKidVy.push_back(kid->vy()); // double
            jetKidVz.push_back(kid->vz()); // double

				//auto mindr = 0.01;
				//auto phocnt = 0;

				//  ----- get unique SClusters ( a group ) that match jet kids & find corosponding rechit groups --------------- 



				// -----------------  end of get kid scluster group & rechits --------------------------------------------------

				//------------------  kid pf cand matching to photon/electron pf cands ------------------------------------------------------
				auto kidcand = pfcands_->ptrAt(kid.key());	
            //std::cout << "Kidcand Energy is : " << kidcand->energy() << " kid energy is : " << kid->energy() << std::endl;	
				const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(kidcand.get());
            if( packed_cand != NULL ){
					//std::cout << "Packed Cand cast worked. " << std::endl;
               //std::cout << " status : " << packed_cand->status()  << std::endl;

					//bool hasPhoton(false);
					int nGedPhotons(0);
					int nMaybeCands(0);
					for( const auto photon : *gedPhotons_ ){
    					edm::RefVector<pat::PackedCandidateCollection> associated =  photon.associatedPackedPFCandidates();
    					for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
        					edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
							//std::cout << "Checking pho asc " << ipc << " with e : " << associatedPtr->energy() << std::endl; 
        					if( associatedPtr.get() == packed_cand ) {  
								nGedPhotons++; 
								//hasPhoton = true;
								//auto phop4 = photon.p4(); 
								//auto op4e = phop4.energy();
								const auto & phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
								auto & hitsAndFractions = phosc->hitsAndFractions();
								auto nrh = hitsAndFractions.size();
								//auto sce = phosc->energy();
								//const auto & seedDetId = phosc->seed()->seed(); 
								//const auto seedRawId = seedDetId.rawId();
								sum_nRechits += nrh;
								auto phe = photon.energy();
                        sum_energy += phe;
								//auto hemf = photon.hadronicOverEm();
								//auto pke = packed_cand->energy();
								//auto ase = associatedPtr->energy();
                        //std::cout << "kidPh e: " << phe << " sce: " << sce << " nRh: " << nrh << " scid: " << seedRawId << " hoverEm: " << hemf;
								//std::cout << " pke: " << pke << " ace: " << ase << std::endl;									
							}//<<>> if( associatedPtr.get() == packed_cand ) 
							//std::cout << " -- Nxt assPkdPFCand -- " << std::endl;
							nMaybeCands++;
    					}//<<>>for( unsigned int ipc = 0; ipc < associated.size(); ipc++ )
						//std::cout << " -- Nxt Photon -- " << std::endl;
					}//<<>>for( const auto photon : *gedPhotons_ )
					//std::cout << " -- Photons done -- " << std::endl;
					if( nGedPhotons > 0 ) hist1d[35]->Fill(1);
					if( nGedPhotons > 1 ) hist1d[35]->Fill(2);
					//if( nGedPhotons > 1 ) {
					//	std::cout << "PackedCand matches to gedPhotons found : " << nGedPhotons << " of " << nMaybeCands << " possible." <<  std::endl;
					//} // if( nGedPhotons > 1 )
/*  ----------------------------  OOT Photons ---------------------------------------------
               int nOotPhotons = 0;
               for( const auto photon : *ootPhotons_ ){
                  edm::RefVector<pat::PackedCandidateCollection> associated =  photon.associatedPackedPFCandidates();
                  for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
                     edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
                     if( associatedPtr.get() == packed_cand ) nOotPhotons++;
                  }//<<>>for( unsigned int ipc = 0; ipc < associated.size(); ipc++ )
               }//<<>>for( const auto photon : *ootPhotons_ )

               if( nOotPhotons > 0 ) hist1d[35]->Fill(2);
               if( nOotPhotons > 1 ) hist1d[35]->Fill(6);
               //std::cout << "PackedCand matches to ootPhotons found : " << nOotPhotons << std::endl;
*/
               int nElectrons = 0;
               for( const auto electron : *electrons_ ){
                  edm::RefVector<pat::PackedCandidateCollection> associated =  electron.associatedPackedPFCandidates();
                  for( unsigned int ipc = 0; ipc < associated.size(); ipc++ ) {
                     edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
							//std::cout << "Checking ele asc " << ipc << " with e : " << associatedPtr->energy() << std::endl;
                     if( associatedPtr.get() == packed_cand ){ 
								nElectrons++;
                        const auto & phosc = electron.superCluster().isNonnull() ? electron.superCluster() : electron.parentSuperCluster();
                        auto & hitsAndFractions = phosc->hitsAndFractions();
                        auto nrh = hitsAndFractions.size();
                        //auto sce = phosc->energy();
                        //const auto & seedDetId = phosc->seed()->seed();
                        //const auto seedRawId = seedDetId.rawId();
                        sum_nRechits += nrh;
								auto ele = electron.energy();
                        sum_energy += ele;
                        //auto hemf = electron.hcalIso();
								//auto pkpt = packed_cand->pt();
                        //std::cout << "kidEle e: " << ele << " sce: " << sce << " nRh: " << nrh << " pkpt: " << pkpt << " id: " << seedRawId << std::endl;
							}
                  }
					//std::cout << " -- Nxt Electron -- " << std::endl;
               }
               if( nElectrons > 0 ) hist1d[35]->Fill(3);
               if( nElectrons > 1 ) hist1d[35]->Fill(4);
               //if( nElectrons > 1 ) std::cout << "PackedCand matches to Electrons found : " << nElectrons << std::endl;

/* ------------- Muons -------------------------------------------------
               int nMuons = 0;
               for( const auto muon : *muons_ ){
                  //edm::RefVector<pat::PackedCandidateCollection> associated =  muon.pfCandidateRef();
						//edm::Ref<std::vector<reco::PFCandidate> >' and 'const pat::PackedCandidate*
						auto candRef = muon.pfCandidateRef();
						auto candPtr = candRef.get();
						auto candSize = candPtr->size(); 
                  for( unsigned int ipc = 0; ipc < candSize; ipc++ ) {
                     //edm::Ptr<pat::PackedCandidate> associatedPtr = edm::refToPtr( associated[ipc] );
							auto assCand = (*associated)[ipc];
                     //if( associatedPtr.get() == packed_cand ) nMuons++;
                     if( assCand == packed_cand ) nMuons++;
                  }
               }
               //std::cout << "PackedCand matches to Muons found : " << nMuons << std::endl;
*/
					//if( kid->isMuon() || kid->isStandAloneMuon() ) std::cout << "Kid particle has Muon flag : " << nMuons << std::endl; 
               if( nGedPhotons > 0 && nElectrons > 0 ) hist1d[35]->Fill(6);
               if( nGedPhotons > 1 && nElectrons > 1 ) hist1d[35]->Fill(7);
               if( nGedPhotons < 1 && nElectrons < 1 ) hist1d[35]->Fill(5);
               //if( nGedPhotons < 1 && nElectrons < 1 ) std::cout << "Kid particle is not ele/photon " << std::endl;
					//if( nGedPhotons < 1 && nOotPhotons < 1 && nElectrons < 1 ) hist1d[35]->Fill(4);

				} //else std::cout << "PackedCand cast failed." << std::endl;
				//-------------------------------------------------------------------------------------------------

				hist1d[5]->Fill(kid->pdgId());
            hist1d[7]->Fill(kid->energy());
            hist1d[8]->Fill(kid->pt());
            hist1d[9]->Fill(kid->eta());
            hist1d[10]->Fill(kid->phi());
            hist1d[19]->Fill(kid->mass());
            hist1d[20]->Fill(kid->charge());
            hist1d[21]->Fill(kid->vx());
            hist1d[22]->Fill(kid->vy());
            hist1d[23]->Fill(kid->vz());

				if( kid->hasMasterClone() ) hist1d[31]->Fill(1); //std::cout << "Has master clone : " << std::endl;
            if( kid->hasMasterClonePtr() ) hist1d[31]->Fill(2); //std::cout << "Has master clone ptr : " << std::endl;
				if( kid->numberOfMothers() ) hist1d[31]->Fill(3); //std::cout << "Has : " << kid->numberOfMothers() << " mothers: " <<  std::endl;
            if( kid->numberOfDaughters() ) hist1d[31]->Fill(4); //std::cout << "Has : " << kid->numberOfDaughters() << " daughters: " <<  std::endl;
            if( kid->numberOfSourceCandidatePtrs() )hist1d[31]->Fill(5); // std::cout << "Has : " << kid->numberOfSourceCandidatePtrs() i

            if( kid->isElectron() ) hist1d[30]->Fill(1); //std::cout << "Is an Electron : " << std::endl;
            if( kid->isMuon() ) hist1d[30]->Fill(2); //std::cout << "Is a Muon : " << std::endl;
            if( kid->isStandAloneMuon() ) hist1d[30]->Fill(3); //std::cout << "Is a StandAloneMuon : " << std::endl;
            if( kid->isGlobalMuon() ) hist1d[30]->Fill(4); //std::cout << "Is a GlobalMuon : " << std::endl;
            if( kid->isTrackerMuon() ) hist1d[30]->Fill(5); //std::cout << "Is an TrackerMuon : " << std::endl;
            if( kid->isCaloMuon() ) hist1d[30]->Fill(6); //std::cout << "Is an CaloMuon : " << std::endl;
            if( kid->isPhoton() ) hist1d[30]->Fill(7); //std::cout << "Is an Photon : " << std::endl;
            if( kid->isConvertedPhoton() ) hist1d[30]->Fill(8); //std::cout << "Is an ConvertedPhoton : " << std::endl;
            if( kid->isJet() ) hist1d[30]->Fill(9); //std::cout << "Is an Jet : " << std::endl;


            //std::cout << "Initing clusterTimes" << std::endl;
            TH1D kidClTimes( "temp", "temp", 500, -25, 25 );
      
            //std::cout << "finding kid RecHits" << std::endl;
            auto jetKidDrRhGroup =  getRHGroup( *recHitsEB_, *recHitsEE_, jet.eta(), jet.phi(), deltaRminKid, minRHenr );
            auto kidRhCount( jetKidDrRhGroup.size() ); 

            if( kidRhCount != 0 ) {
            	auto kidTofTimes = getRhTofTime( jetKidDrRhGroup, vtxX, vtxY, vtxZ );
            	for (const auto kidTime : kidTofTimes ) kidClTimes.Fill(kidTime);
               //  make kidtime varible
               //std::cout << "Finding kid rh mean time" << std::endl;
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
					hist1d[27]->Fill(ktime);
					hist1d[28]->Fill(kmedtime);
            	hist1d[29]->Fill(kidRhCount);
				   kidRHmTime.Fill(ktime);	
				}//<<>> if( kidRhCount != 0 )  eof kid time from rechit collection

/*
				// get time from SC collction
            std::cout << "Finding kid SC" << std::endl;
            //const reco::SuperCluster* kidSC(NULL);
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
	                  hist1d[32]->Fill(rht-tof);
					}
	            for (const auto recHit : *recHitsEE_ ){
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
							hist1d[32]->Fill(rht-tof);
	            }
				}
*/
			//std::cout << "Next kid" << std::endl;
			}//<<>>for( const auto kid : jet.daughterPtrVector() )

			//std::cout << "Finding jet time from kids" << std::endl;
         auto sctime = kidSCTime.GetMean();
			if( sctime == 0.0 ) sctime = -9.0;
         hist1d[33]->Fill(sctime);
			auto rhmtime = kidRHmTime.GetMean();
         if( rhmtime == 0.0 ) rhmtime = -9.0;
         hist1d[34]->Fill(rhmtime);
         hist1d[41]->Fill(sum_energy);
         hist1d[42]->Fill(sum_nRechits);

		}//<<>>if( nKids > 0 )
		//--- end of jet kids -------------------------------------

		hist1d[6]->Fill(nKids);
		
		// SC rh count -----------------------------------------------

      int iph(0);
      bool matched(false);
		int nMatched(0);
		int sum_nrh(0);
      float sum_sce(0.0);
      float sum_phe(0.0);
      for( const auto photon : *gedPhotons_ ){
         //std::cout << "Proccesssing New Photon :" << std::endl;
         edm::RefVector<pat::PackedCandidateCollection> passociated =  photon.associatedPackedPFCandidates();
         for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ ) {
            //std::cout << "Processing asc pfcand # " << ipcp << std::endl;
            edm::Ptr<pat::PackedCandidate> passociatedPtr = edm::refToPtr( passociated[ipcp] );
            const auto *ascpacked_cand = passociatedPtr.get();
               int ijk(0);
               for( const auto kid : jet.daughterPtrVector() ){
                  //std::cout << "Proccesssing New jetKid :" << std::endl;
                  auto kidcand = pfcands_->ptrAt(kid.key());
                  const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(kidcand.get());
                     if( ascpacked_cand == packed_cand ){

                           //const auto pce = ascpacked_cand->energy();
                           //std::cout << "Match found at ph: " << iph << " asc: " << ipcp << " jet: "; 
									//std::cout << ijet << "/" << nJets  << " kid: " << ijk; // << std::endl;
                           //std::cout << " pce: " << pce << std::endl;
                           matched = true;

                     }//<<>>if( ascpacked_cand == packed_cand )
                  ijk++;
               }//<<>>for( const auto kid : jet.daughterPtrVector() )
         }//<<>>for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ )
         iph++;
         if( matched ){
				nMatched++;
            const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
            const auto &hitsAndFractions = phosc->hitsAndFractions();
            const auto nrh = hitsAndFractions.size();
            const auto sce = phosc->energy();
            //const auto &seedDetId = phosc->seed()->seed();
            //const auto seedRawId = seedDetId.rawId();
            //const auto clusters = phosc->clusters();
            //const auto nBClusts = phosc->clustersSize();
            //std::cout << " Match to photon SC with seedid: " << seedRawId << " nrh: " << nrh << " sce: " << sce << " nClust: " << nBClusts << std::endl;
				//if( nrh != 0 ){
					//hist2d[39]->Fill(rhCount,nrh); 
					sum_nrh += nrh;
            	//hist2d[40]->Fill(jet.energy(),sce);
					sum_sce += sce;
            	//hist2d[41]->Fill(jet.energy(),photon.energy());
					sum_phe += photon.energy();
            	//hist2d[42]->Fill(photon.energy(),sce);
           	//}
            matched = false;
				/*
            auto nClust(0);
            if( nrh == 0 ){ std::cout << "  Skipping SC Cluster - bad reffrence issue" << std::endl;  continue; }
            for( const auto &clustptr : clusters ){
               const auto clust = clustptr.get();
               const auto cle = clust->energy();
               const auto &clDetId = clust->seed();
               const auto clRawId = clDetId.rawId();
               const auto &clHitsAndFractions = clust->hitsAndFractions();
               const auto nclrh = clHitsAndFractions.size();
               std::cout << "  SC Clust # " << nClust << " energy: " << cle << " nrh: " << nclrh << " rawID: " << clRawId << std::endl;
               nClust++;
            }//<<>>for( const auto cluster : clusters )
				*/
         }//<<>>if( matched )
      }//<<>>for( const auto photon : *gedPhotons_ )      
      hist2d[39]->Fill(rhCount,sum_nrh);
      hist2d[40]->Fill(jet.photonEnergy(),sum_sce);
      hist2d[41]->Fill(jet.photonEnergy(),sum_phe);
      hist2d[42]->Fill(sum_phe,sum_sce);
		hist1d[43]->Fill(nMatched);
	//****************************  photon/electron to kid pfcand -> SC matcher ********************************

		//std::cout << "Next Jet .......................... " << std::endl; 
   }//<<>>for (auto i = 0; i < nJets; i++)
// ** end of jets  ***************************************************************************************************

   hist1d[17]->Fill(jetHt);

//-------------------------------------------------------------------------------
// matching packedCanadates in photon and electron collections
   int nMatches = 0;
   int nPhoMatches = 0;
   int nEleMatches = 0;

   for( const auto photon : *gedPhotons_ ){
      edm::RefVector<pat::PackedCandidateCollection> passociated =  photon.associatedPackedPFCandidates();
		hist1d[38]->Fill(passociated.size());
      for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ ) {
         edm::Ptr<pat::PackedCandidate> passociatedPtr = edm::refToPtr( passociated[ipcp] );
			for( const auto electron : *electrons_ ){
				edm::RefVector<pat::PackedCandidateCollection> eassociated =  electron.associatedPackedPFCandidates();
				hist1d[39]->Fill(eassociated.size());
				for( unsigned int ipce = 0; ipce < eassociated.size(); ipce++ ) {
					edm::Ptr<pat::PackedCandidate> eassociatedPtr = edm::refToPtr( eassociated[ipce] );
					if( eassociatedPtr.get() == passociatedPtr.get() ){ nMatches++; nPhoMatches++;}
				}//<<>>for( unsigned int ipce = 0; ipce < eassociated.size(); ipce++ )
			}//<<>>for( const auto electron : *electrons_ )
			hist1d[37]->Fill(nPhoMatches);
			nPhoMatches = 0;
      }//<<>>for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ )
   }//<<>>for( const auto photon : *gedPhotons_ )		

   hist1d[36]->Fill(nMatches);

   for( const auto electron : *electrons_ ){
      edm::RefVector<pat::PackedCandidateCollection> eassociated =  electron.associatedPackedPFCandidates();
      for( unsigned int ipce = 0; ipce < eassociated.size(); ipce++ ) {
         edm::Ptr<pat::PackedCandidate> eassociatedPtr = edm::refToPtr( eassociated[ipce] );
         for( const auto photon : *gedPhotons_ ){
            edm::RefVector<pat::PackedCandidateCollection> passociated =  photon.associatedPackedPFCandidates();
            for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ ) {
               edm::Ptr<pat::PackedCandidate> passociatedPtr = edm::refToPtr( passociated[ipcp] );
               if( passociatedPtr.get() == eassociatedPtr.get() ){ nEleMatches++;}
            }//<<>>for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ )
         }//<<>>for( const auto photon : *gedPhotons_ )
         hist1d[40]->Fill(nEleMatches);
         nEleMatches = 0;
      }//<<>>for( unsigned int ipce = 0; ipce < eassociated.size(); ipce++ )
   }//<<>>for( const auto electron : *electrons_ )
//-------------------------------------------------------------------------------------
//****************************  photon/electron to kid pfcand -> SC matcher ********************************
/*
	std::cout << "Processing New Event : " << std::endl;
   int nCalo(0);
   for( const auto calo : *caloCluster_ ){
		
		const auto energy = calo.energy();
		const auto seedDetId = calo.seed();
		const auto seedRawId = seedDetId.rawId();
      const auto &hitsAndFractions = calo.hitsAndFractions();
      const auto nrh = hitsAndFractions.size();
		std::cout << "Calo # " << nCalo << " energy: " << energy << " nrh: " << nrh << " rawID: " << seedRawId;
      std::cout << " Cl Centroid eta: " << calo.eta() << " phi: " << calo.phi() << std::endl;
		nCalo++; 

   }//<<>>for( const auto photon : *caloCluster_ )

	for (auto ijet = 0; ijet < nJets; ijet++){
		const auto &jet = fjets[ijet];
		std::cout << "Jet # " << ijet << " energy: " << jet.energy() << " eta: " << jet.eta() << " phi: " << jet.phi() << std::endl;
	}

	int iph(0);
   bool matched(false);
   for( const auto photon : *gedPhotons_ ){
		//std::cout << "Proccesssing New Photon :" << std::endl;
      edm::RefVector<pat::PackedCandidateCollection> passociated =  photon.associatedPackedPFCandidates();
      for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ ) {
			//std::cout << "Processing asc pfcand # " << ipcp << std::endl;
         edm::Ptr<pat::PackedCandidate> passociatedPtr = edm::refToPtr( passociated[ipcp] );
			const auto *ascpacked_cand = passociatedPtr.get();

   		for (auto ijet = 0; ijet < nJets; ijet++){
				//std::cout << "Processing jet # " << ijet << std::endl;
				const auto &jet = fjets[ijet];
				//std::cout << "Jet # " << ijet << " energy: " << jet.energy() << " eta: " << jet.eta() << " phi: " << jet.phi() << std::endl;
				int ijk(0);
				for( const auto kid : jet.daughterPtrVector() ){
					//std::cout << "Proccesssing New jetKid :" << std::endl;
            	auto kidcand = pfcands_->ptrAt(kid.key());
            	const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(kidcand.get());
            	//if( packed_cand != NULL ){
						if( ascpacked_cand == packed_cand ){

					         const auto pce = ascpacked_cand->energy();
                        const auto pceta = ascpacked_cand->eta();
                        const auto pcphi = ascpacked_cand->phi();
								std::cout << "Match found at ph: " << iph << " asc: " << ipcp << " jet: " << ijet << "/" << nJets  << " kid: " << ijk; // << std::endl;
								std::cout << " pce: " << pce << " eta: " << pceta << " phi: " << pcphi << std::endl;
								//std::cout << " jet e: " << jet.energy() << " eta: " << jet.eta() << " phi: " << jet.phi() << std::endl;
								matched = true;

						}//<<>>if( ascpacked_cand == packed_cand )
					//}//<<>>if( packed_cand != NULL )
					ijk++;
				}//<<>>for( const auto kid : jet.daughterPtrVector() )
			}//<<>>for (auto ijet = 0; ijet < nJets; ijet++)
      }//<<>>for( unsigned int ipcp = 0; ipcp < passociated.size(); ipcp++ )
		iph++;
		if( matched ){
			const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
			const auto &hitsAndFractions = phosc->hitsAndFractions();
			const auto nrh = hitsAndFractions.size();
			const auto sce = phosc->energy();
			const auto &seedDetId = phosc->seed()->seed();
			const auto seedRawId = seedDetId.rawId();
         const auto clusters = phosc->clusters();
         const auto nBClusts = phosc->clustersSize();
			std::cout << " Match to photon SC with seedid: " << seedRawId << " nrh: " << nrh << " sce: " << sce << " nClust: " << nBClusts;
			std::cout << " eta: " << phosc->eta() << " phi: " << phosc->phi() << std::endl;
			matched = false;
			auto nClust(0);
			if( nrh == 0 ){ std::cout << "  Skipping SC Cluster - bad reffrence issue" << std::endl;  continue; }
			for( const auto &clustptr : clusters ){
				const auto clust = clustptr.get();
      		const auto cle = clust->energy();
      		const auto &clDetId = clust->seed();
				const auto clRawId = clDetId.rawId();
      		const auto &clHitsAndFractions = clust->hitsAndFractions();
      		const auto nclrh = clHitsAndFractions.size();
      		std::cout << "  SC Clust # " << nClust << " energy: " << cle << " nrh: " << nclrh << " rawID: " << clRawId;
            std::cout << " eta: " << clust->eta() << " phi: " << clust->phi() << std::endl;				 
				nClust++;
			}//<<>>for( const auto cluster : clusters )
		}//<<>>if( matched )
   }//<<>>for( const auto photon : *gedPhotons_ )      
*/
//****************************  photon/electron to kid pfcand -> SC matcher ********************************

//------------------------------------------------------------------------------------
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
         hist1d[24]->Fill(diffPt);
         if( diffPt < difPtLmt ) continue;
         auto htPct= (qjet.pt()+pjet.pt())/jetHt;
			hist1d[25]->Fill(htPct);
         if( htPct < htPctLmt ) continue;
         auto dPhi = reco::deltaPhi(qjet.phi(),pjet.phi());
         hist1d[26]->Fill(dPhi);
         if( dPhi < dPhiLmt ) continue;
			auto dTmu = jetTime[q] - jetTime[p];
         auto dTmed = jetMedTime[q] - jetMedTime[p];
         if( dTmu == 0.0 )  dTmu = -5.5;
         if( dTmed == 0.0 )  dTmed = -5.5;
         if( jetTime[q] == 0.0 || jetTime[p] == 0.0 ) dTmu = -5.0;
         if( jetMedTime[q] == 0.0 || jetMedTime[p] == 0.0 ) dTmed = -5.0;
			hist1d[15]->Fill(dTmu);
			hist2d[22]->Fill(dTmu,nJets);
         hist2d[26]->Fill(dTmu,diffPt);
         hist2d[27]->Fill(dTmu,htPct);
         hist2d[28]->Fill(dTmu,dPhi);
			hist1d[16]->Fill(dTmed);
         hist2d[23]->Fill(dTmed,nJets);
         hist2d[29]->Fill(dTmu,diffPt);
         hist2d[30]->Fill(dTmu,htPct);
         hist2d[31]->Fill(dTmu,dPhi);
		}//<<>>for (auto p = q+1; p < nJets; p++)
	}//<<>>for (auto q = 0; q < nJets; q++)
//-------------------------------------------------------------------------------

// -- Fill output trees ------------------------------------------
  //std::cout << "---------- Next Event -----" << std::endl;
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
  hist1d[1] = fs->make<TH1D>("jetRHMulti" , "jetRHMulti" , 250 , 0 , 250 );
  hist1d[2] = fs->make<TH1D>("jetTimeError" , "jetTimeError" , 300 , 0 , 3 );
  hist1d[3] = fs->make<TH1D>("jetTimeRMS" , "jetTimeRMS" , 200 , 0 , 20 );
  hist1d[4] = fs->make<TH1D>("jetMedTime" , "jetMedTime" , 200 , -10 , 10 );
  hist1d[5] = fs->make<TH1D>("kidPdgID" , "kidPdgID" , 300 , 0 , 300 );
  hist1d[6] = fs->make<TH1D>("nKids" , "nKids" , 100 , 0 , 100 );
  hist1d[7] = fs->make<TH1D>("kidE" , "kidE" , 500 , 0 , 500 );
  hist1d[8] = fs->make<TH1D>("kidPt" , "kidPt" , 500, 0 , 500 );
  hist1d[9] = fs->make<TH1D>("kidEta" , "kidEta" , 700 , -3.5 , 3.5 );
  hist1d[10] = fs->make<TH1D>("kidPhi" , "kidPhi" , 700 , -3.5 , 3.5 );
  //hist1d[11] = fs->make<TH1D>("t0jetRHMulti" , "t0jetRHMulti" , 250 , 0 , 250 );

  hist1d[12] = fs->make<TH1D>("jetPt" , "jetPt" , 500, 0 , 500 );
  hist1d[13] = fs->make<TH1D>("jetPhi" , "jetPhi" , 700 , -3.5 , 3.5 );
  hist1d[14] = fs->make<TH1D>("jetEta" , "jetEta" , 700 , -3.5 , 3.5 );
  hist1d[15] = fs->make<TH1D>("jetdtmu" , "jetdtmu" , 120 , -6 , 6 );
  hist1d[16] = fs->make<TH1D>("jetdtmed" , "jetdtmed" , 120 , -6 , 6 );
  hist1d[17] = fs->make<TH1D>("jetHt" , "jetHt" , 1000 , 0 , 1000 );
  hist1d[18] = fs->make<TH1D>("nJet" , "nJets" , 20 , 0 , 20 );
  hist1d[19] = fs->make<TH1D>("kidMass" , "kidMass" , 20, -1.0 , 1.0 );
  hist1d[20] = fs->make<TH1D>("kidChrg" , "kidChrg" , 5, -2 , 2 );
  hist1d[21] = fs->make<TH1D>("kidVx" , "kidVx" , 200, -10 , 10 );
  hist1d[22] = fs->make<TH1D>("kidVy" , "kidVy" , 200, -10 , 10 );
  hist1d[23] = fs->make<TH1D>("kidVz" , "kidVz" , 200, -10 , 10 );

  hist1d[24] = fs->make<TH1D>("diffPt" , "diffPt" , 100, 0 , 1 );
  hist1d[25] = fs->make<TH1D>("htPct" , "htPct" , 100, 0 , 1 );
  hist1d[26] = fs->make<TH1D>("dPhi" , "dPhi" , 32, 0 , 3.2 );

  hist1d[27] = fs->make<TH1D>("kidTime" , "kidTime" , 200 , -10 , 10 );
  hist1d[28] = fs->make<TH1D>("kidMedTime" , "kidMedTime" , 200 , -10 , 10 );
  hist1d[29] = fs->make<TH1D>("kidRhCount", "kidRhCount", 250 , 0 , 250 );
  hist1d[30] = fs->make<TH1D>("kidleafptrs", "kidleafptrs", 10 , 0 , 10 );
  hist1d[31] = fs->make<TH1D>("kidleaftype", "kidleaftype", 10 , 0 , 10 );
  hist1d[32] = fs->make<TH1D>("kidSCTime" , "kidSCTime" , 500 , -25 , 25 );
  hist1d[33] = fs->make<TH1D>("jetSCTime" , "jetSCTime" , 500 , -25 , 25 );
  hist1d[34] = fs->make<TH1D>("jetRHmTime" , "jetRHmTime" , 500 , -25 , 25 );
  hist1d[35] = fs->make<TH1D>("kidpfCand", "kidpfCand", 10 , 0 , 10 );
  hist1d[36] = fs->make<TH1D>("nEGMatches", "nEGMatches", 100 , 0 , 100 ); //]->Fill(nMatches);
  hist1d[37] = fs->make<TH1D>("nEpfPhoton", "nEpfPhoton", 100 , 0 , 100 ); //]->Fill(nPhoMatches);
  hist1d[38] = fs->make<TH1D>("nPPfcands", "nPPfcands", 100 , 0 , 100 ); //]->Fill(passociated.size());
  hist1d[39] = fs->make<TH1D>("nEPfcands", "nEPfcands", 100 , 0 , 100 );//]->Fill(eassociated.size());
  hist1d[40] = fs->make<TH1D>("nPpfElectron", "nPpfElectron", 100 , 0 , 100 ); //]->Fill(nEleMatches);
  hist1d[41] = fs->make<TH1D>("sumKidPhEnergy", "sumKidPhEnergy", 1000 , 0 , 1000 );
  hist1d[42] = fs->make<TH1D>("sumKidNPhRechits", "sumKidNPhRechits", 100 , 0 , 100 );
  hist1d[43] = fs->make<TH1D>("nPhotonsPerJet","nPhotonsPerJet", 20 , 0 , 20 );

//------------------------------------------------------------------------------------
  hist2d[1] = fs->make<TH2D>("jt_pt" , "jt_pt" , 200 , -10 , 10 , 500 , 0 , 500 );
  hist2d[2] = fs->make<TH2D>("jt_id" , "jt_id" , 200 , -10 , 10 , 5 , 0 , 5 );
  hist2d[3] = fs->make<TH2D>("jt_nhf" , "jt_nhf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[4] = fs->make<TH2D>("jt_chf" , "jt_chf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[5] = fs->make<TH2D>("jt_nemf" , "jt_nemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[6] = fs->make<TH2D>("jt_cemf" , "jt_cemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[7] = fs->make<TH2D>("jt_muf" , "jt_muf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[8] = fs->make<TH2D>("jt_nhm" , "jt_nhm" , 200 , -10 , 10 , 40 , 0 , 40 );
  hist2d[9] = fs->make<TH2D>("jt_chm" , "jt_chm" , 200 , -10 , 10 , 40 , 0 , 40 );

  hist2d[10] = fs->make<TH2D>("jt_medt" , "jt_medt" , 200 , -10 , 10 , 200 , -10 , 10 );
  hist2d[11] = fs->make<TH2D>("jt_rms" , "jt_rms" , 200 , -10 , 10 , 200 , 0 , 20 );
  hist2d[12] = fs->make<TH2D>("jt_err" , "jt_err" , 200 , -10 , 10 , 300 , 0 , 3 );

  hist2d[13] = fs->make<TH2D>("medt_pt" , "medt_pt" , 200 , -10 , 10 , 500 , 0 , 500 );
  hist2d[14] = fs->make<TH2D>("medt_id" , "medt_id" , 200 , -10 , 10 , 5 , 0 , 5 );
  hist2d[15] = fs->make<TH2D>("medt_nhf" , "medt_nhf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[16] = fs->make<TH2D>("medt_chf" , "medt_chf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[17] = fs->make<TH2D>("medt_nemf" , "medt_nemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[18] = fs->make<TH2D>("medt_cemf" , "medt_cemf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[19] = fs->make<TH2D>("medt_muf" , "medt_muf" , 200 , -10 , 10 , 100 , 0 , 1 );
  hist2d[20] = fs->make<TH2D>("medt_nhm" , "medt_nhm" , 200 , -10 , 10 , 40 , 0 , 40 );
  hist2d[21] = fs->make<TH2D>("medt_chm" , "medt_chm" , 200 , -10 , 10 , 40 , 0 , 40 );

  hist2d[22] = fs->make<TH2D>("jdtmu_nJets" , "jdtmu_nJets" , 120 , -6 , 6 , 6 , 2 , 8 );
  hist2d[23] = fs->make<TH2D>("jdtmed_nJets" , "jdtmed_nJets" , 120 , -6 , 6 , 6 , 2 , 8 );

  hist2d[24] = fs->make<TH2D>("jt_nrh" , "jt_nrh" , 200 , -10 , 10 , 250 , 0 , 250 );
  hist2d[25] = fs->make<TH2D>("medt_nrh" , "medt_nrh" , 200 , -10 , 10 , 250 , 0 , 250 );

  hist2d[26] = fs->make<TH2D>("jdtmu_diffPt" , "jdtmu_diffPt" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d[27] = fs->make<TH2D>("jdtmu_htPct" , "jdtmu_htPct" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d[28] = fs->make<TH2D>("jdtmu_dPhi" , "jdtmu_dPhi" , 120 , -6 , 6 , 400 , 2.8 , 3.2 );
  hist2d[29] = fs->make<TH2D>("jdtmed_diffPt" , "jdtmed_diffPt" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d[30] = fs->make<TH2D>("jdtmed_htPct" , "jdtmed_htPct" , 120 , -6 , 6 , 200 , 0.8 , 1 );
  hist2d[31] = fs->make<TH2D>("jdtmed_dPhi" , "jdtmed_dPhi" , 120 , -6 , 6 , 400 , 2.8 , 3.2 );

  hist2d[32] = fs->make<TH2D>("jt_sceta" , "jt_sceta", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d[33] = fs->make<TH2D>("jt_scphi" ,"jt_scphi", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d[34] = fs->make<TH2D>("jt_scenr" , "jt_scenr", 200 , -10 , 10 , 1000 , 0 , 1000 );

  hist2d[35] = fs->make<TH2D>("medt_sceta" , "medt_sceta", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d[36] = fs->make<TH2D>("medt_scphi" ,"medt_scphi", 200 , -10 , 10 , 700 , -3.5 , 3.5 );
  hist2d[37] = fs->make<TH2D>("medt_scenr" , "medt_scenr", 200 , -10 , 10 , 1000 , 0 , 1000 );

  hist2d[38] = fs->make<TH2D>("rht_rhe" , "rht_rhe", 200 , -10 , 10 , 1000 , 0 , 1000 );

  hist2d[39] = fs->make<TH2D>("njrh_nscrh" , "njrh_nscrh" , 200 , 0 , 200 , 200 , 0 , 200 );
  hist2d[40] = fs->make<TH2D>("jphe_sce" , "jphe_sce" , 500 , 0 , 500 , 500 , 0 , 500 );
  hist2d[41] = fs->make<TH2D>("jphe_phe" , "jphe_phe" , 500 , 0 , 500 , 500 , 0 , 500 );
  hist2d[42] = fs->make<TH2D>("phe_sce" , "phe_sce" , 500 , 0 , 500 , 500 , 0 , 500 );

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
