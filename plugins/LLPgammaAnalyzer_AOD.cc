// -*- C++ -*-
//
// Package:    LLPgammaAnalyzer_AOD
// Class:      LLPgammaAnalyzer_AOD
//
/**\class LLPgammaAnalyzer_AOD LLPgammaAnalyzer_AOD.cc LLPgammaAnalyzer/plugins/LLPgammaAnalyzer_AOD.cc

 		Description: [one line class summary]

 		Implementation: [Notes on implementation]

*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//
//
//----------------------------------------  cc file   --------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "LLPGamma/LLPgammaAnalyzer/plugins/LLPgammaAnalyzer_AOD.hh"
using namespace std;

//#define DEBUG true
#define DEBUG false

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// constructors and destructor
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

LLPgammaAnalyzer_AOD::LLPgammaAnalyzer_AOD(const edm::ParameterSet& iConfig) :

// -- declare tags ----------------------------------------------------------

	// flags
	hasGenInfo(iConfig.existsAs<bool>("hasGenInfo")  ? iConfig.getParameter<bool>("hasGenInfo")  : true),

	// triggers
	//triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
	//triggerObjectsTag(iConfig.getParameter<edm::InputTag>("triggerObjects")),

	// tracks
	tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),
	
	// pfcands
	pfcandTag(iConfig.getParameter<edm::InputTag>("pfcandidates")),
    pfCanTag(iConfig.getParameter<edm::InputTag>("particleflow")),
	pfCanPhoMapTag(iConfig.getParameter<edm::InputTag>("pfcanphomap")),
    pfCanOOTPhoMapTag(iConfig.getParameter<edm::InputTag>("pfcanootphomap")),
    pfCanEleMapTag(iConfig.getParameter<edm::InputTag>("pfcanelemap")),
	
	// vertices
	verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),

	// rho
	//rhoTag(iConfig.getParameter<edm::InputTag>("rho")),

	// mets
	pfmetsTag(iConfig.getParameter<edm::InputTag>("mets")),  

	// supercluster
	superClusterCollectionTag(iConfig.getParameter<edm::InputTag>("superClusters")),
	ootSuperClusterCollectionTag(iConfig.getParameter<edm::InputTag>("ootSuperClusters")),

	// caloclusters
	caloClusterTag(iConfig.getParameter<edm::InputTag>("caloClusters")),

	// jets
	jetsTag(iConfig.getParameter<edm::InputTag>("jets")), 

    // calojets
    caloJetsTag(iConfig.getParameter<edm::InputTag>("calojets")),

    // genJets
    genJetsTag(iConfig.getParameter<edm::InputTag>("genjets")),

	// electrons
	electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),  

	// muons
	muonsTag(iConfig.getParameter<edm::InputTag>("muons")),  

	// recHits
	recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),  
	recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

	// gedphotons
	gedPhotonsTag(iConfig.getParameter<edm::InputTag>("gedPhotons")),
    phoCBIDLooseMapTag(iConfig.getParameter<edm::InputTag>("phoCBIDLooseMap")),

	// ootPhotons
	ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons")),

	// pfcand ref
	//reco2pfTag(iConfig.getParameter<edm::InputTag>("recoToPFMap")),
	//reco2pf_(mayConsume<edm::ValueMap<std::vector<reco::PFCandidateRef> > >(iConfig.getParameter<edm::InputTag>("recoToPFMap"))),

	// genParticles
	genEvtInfoTag(iConfig.getParameter<edm::InputTag>("genEvt")),
  	gent0Tag(iConfig.getParameter<edm::InputTag>("gent0")),
  	genxyz0Tag(iConfig.getParameter<edm::InputTag>("genxyz0")),
  	pileupInfosTag(iConfig.getParameter<edm::InputTag>("pileups")),
  	genParticlesTag(iConfig.getParameter<edm::InputTag>("genParticles"))

// -- end of tag declarations ---------------------------------------
{ //<<<< LLPgammaAnalyzer_AOD::LLPgammaAnalyzer_AOD(const edm::ParameterSet& iConfig) :

	usesResource();
	usesResource("TFileService");

// -- consume tags ------------------------------------------------------------
	if( DEBUG ) std::cout << "In constructor for LLPgammaAnalyzer_AOD - tag and tokens" << std::endl;

	// Triggers
	//triggerResultsToken_ 	= consumes<edm::TriggerResults>(triggerResultsTag);
	//triggerObjectsToken_ 	= consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectsTag);

	// tracks 
	tracksToken_				= consumes<std::vector<reco::Track>>(tracksTag);

	// genparticles
	// genpart_token_              = consumes<std::vector<pat::PackedGenParticle>>(genpartTag);

	// pfcandidates
	pfCan_token_        		= consumes<std::vector<reco::PFCandidate>>(pfCanTag);
	pfCanPhotonMap_token_		= consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate>>>(pfCanPhoMapTag);
    pfCanOOTPhotonMap_token_    = consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate>>>(pfCanOOTPhoMapTag);
    pfCanElectronMap_token_     = consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate>>>(pfCanEleMapTag);

	// pfcandidates view
    pfcand_token_               = consumes<CandidateView>(pfcandTag);

	// vertices
	verticesToken_				= consumes<std::vector<reco::Vertex>>(verticesTag);

	// rho
	//rhoToken_					= consumes<double>(rhoTag);

	// mets
	pfmetsToken_					= consumes<std::vector<reco::PFMET>>(pfmetsTag);

	// supercluster
	scToken_             		= consumes<reco::SuperClusterCollection>(superClusterCollectionTag);
	ootScToken_          		= consumes<reco::SuperClusterCollection>(ootSuperClusterCollectionTag); 

	// caloClusters
	ccToken_			     	= consumes<std::vector<reco::CaloCluster>>(caloClusterTag);	

	// jets
	jetsToken_					= consumes<std::vector<reco::PFJet>>(jetsTag);
    caloJetsToken_              = consumes<std::vector<reco::CaloJet>>(caloJetsTag);

	// leptons
	electronsToken_				= consumes<std::vector<reco::GsfElectron>>(electronsTag);

	// muons
	muonsToken_					= consumes<std::vector<reco::Muon>>(muonsTag);

	// rechits
	recHitsEBToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
	recHitsEEToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

	// photons
	gedPhotonsToken_ 			= consumes<std::vector<reco::Photon>>(gedPhotonsTag);
	phoCBIDLooseMapToken_		= consumes<edm::ValueMap<bool>>(phoCBIDLooseMapTag); 
	ootPhotonsToken_ 			= consumes<std::vector<reco::Photon>>(ootPhotonsTag);

	// pfref
	//reco2pf_					= consumes<edm::ValueMap<std::vector<reco::PFCandidateRef>>>(reco2pfTag);

	//genParticles & genJets
	if( hasGenInfo ){

		genJetsToken_ 			= consumes<std::vector<reco::GenJet>>(genJetsTag);
		genEvtInfoToken_   		= consumes<GenEventInfoProduct>(genEvtInfoTag);
    	gent0Token_        		= consumes<float>(gent0Tag);
    	genxyz0Token_      		= consumes<Point3D>(genxyz0Tag);
    	pileupInfosToken_  		= consumes<std::vector<PileupSummaryInfo>>(pileupInfosTag);
    	genParticlesToken_ 		= consumes<std::vector<reco::GenParticle>>(genParticlesTag);

	}//<<>>if( hasGenInfo )

// ---------------------------------------------------------------------------------
}//>>>>LLPgammaAnalyzer_AOD::LLPgammaAnalyzer_AOD(const edm::ParameterSet& iConfig)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

LLPgammaAnalyzer_AOD::~LLPgammaAnalyzer_AOD(){
	///////////////////////////////////////////////////////////////////
	// do anything here that needs to be done at desctruction time   //
	// (e.g. close files, deallocate resources etc.)                 //
	///////////////////////////////////////////////////////////////////
}//>>>>LLPgammaAnalyzer_AOD::~LLPgammaAnalyzer_AOD()


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// member functions
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
detIdMap LLPgammaAnalyzer_AOD::SetupDetIDs(){

	detIdMap DetIDMap;
	std::string detIDConfig(ecal_config_path);
	uInt  cmsswId, dbID;
	int hashedId, iphi, ieta, absieta, Fed, SM, TT25, iTT, strip, Xtal, phiSM, etaSM;
	int side, ix, iy, SC, iSC, TTCCU, quadrant;
	TString pos;
	auto detIDConfigEB = detIDConfig + "fullinfo_detids_EB.txt";
	std::ifstream infileEB( detIDConfigEB, std::ios::in);
	if( DEBUG ) std::cout << "Setting up EB DetIDs with " << &infileEB << std::endl;
	while( infileEB >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> Fed >> SM >> TT25 >> iTT >> strip >> Xtal >> phiSM >> etaSM )
    	{ DetIDMap[cmsswId] = {iphi,ieta,TT25,ECAL::EB};}
	auto detIDConfigEE = detIDConfig + "fullinfo_detids_EE.txt";
	std::ifstream infileEE( detIDConfigEE, std::ios::in);
	if( DEBUG ) std::cout << "Setting up EE DetIDs with " << &infileEE << std::endl;
	while( infileEE >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> pos >> TTCCU >> strip >> Xtal >> quadrant )
    	{ DetIDMap[cmsswId] = {ix,iy,TTCCU,((side>0) ? ECAL::EP : ECAL::EM )};}

	return DetIDMap;

}//>>>>detIdMap LLPgammaAnalyzer_AOD::SetupDetIDs()
*/

int LLPgammaAnalyzer_AOD::getPFJetID(const reco::PFJet & jet){ 
	
	const auto eta  = std::abs(jet.eta());
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
	} else if (!eta1 && eta2 ){ //<<>>if (eta1)
	 	if      (nhfup && nemfup2 && chm0 && muf8 && cemf8) return tightLepVeto;
	    else if (nhfup && nemf80 && chm0) return tighter;
	 	else if (nhfup && nemfup2 && chm0) return tight;
		else    return loose; 
	} else if (!eta2 && eta3){ //<<>>if (eta1) : else if
	   	if      (nemf10 && nemf80 && nhm2) return tighter;
		else if (nemflw && nemfup2 && nhm2) return tight;
	 	else    return loose; 
	} else { //<<>>if (eta1) : else if : else if
	 	if      (nhflw && nemfup1 && nhm10) return tight;
	 	else    return loose; 
	}//<<>>if (eta1) : else if : else if : else

	return -1; // should not happen

}//>>>>int LLPgammaAnalyzer_AOD::getPFJetID(const reco::PFJet & jet)

rhGroup LLPgammaAnalyzer_AOD::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 ){

	rhGroup result;
	for (const auto recHit : *recHitsEB_ ){
		if( recHit.checkFlag(EcalRecHit::kHasSwitchToGain6) ) continue;
		auto enr = recHit.energy();
		if( enr <= minenr ) continue;
	   	const auto recHitId(recHit.detid());
	   	const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
	   	const auto dr = std::sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
	   	if( dr > drmin ) continue;
	   	result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEB_ )
	for (const auto recHit : *recHitsEE_ ){
	   	if( recHit.checkFlag(EcalRecHit::kHasSwitchToGain6) ) continue;
	   	auto enr = recHit.energy();
	   	if( enr <= minenr ) continue;
	   	const auto recHitId(recHit.detid());
	   	const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
	   	const auto dr = std::sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
	   	if( dr > drmin ) continue;
	   	result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEE_ )

	return result;

}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 )

rhGroup LLPgammaAnalyzer_AOD::getRHGroup( uInt detid ){

	rhGroup result;
	for (const auto recHit : *recHitsEB_ ){
	   	const auto recHitId = recHit.detid();
		const auto rawId = recHitId.rawId();
	   	if( detid != rawId ) continue;
	   	result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEB_ )
	for (const auto recHit : *recHitsEE_ ){
	   	const auto recHitId = recHit.detid();
	   	const auto rawId = recHitId.rawId();
	   	if( detid != rawId ) continue;
	   	result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEE_ )

	return result;
}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup( uInt detid )

rhGroup LLPgammaAnalyzer_AOD::getRHGroup(){

	rhGroup result;
	for (const auto recHit : *recHitsEB_ ) result.push_back(recHit);
	for (const auto recHit : *recHitsEE_ ) result.push_back(recHit);

	return result;

}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup()


rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr ){

	rhGroup result;
	vector<uInt> rawIds;
	for ( const auto superCluster : superClusterGroup ){
	 	auto & hitsAndFractions = superCluster.hitsAndFractions();
	 	const auto nHAF = hitsAndFractions.size();
	 	for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
	    	const auto detId = hitsAndFractions[iHAF].first;
	    	const auto rawId = detId.rawId();
	    	if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
	 	}//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
	}//<<>>for ( const auto superCluster : superClusterGroup )	
	for (const auto recHit : *recHitsEB_ ){
	  	auto enr = recHit.energy();
	  	if( enr <= minenr ) continue;
	  	const auto recHitId = recHit.detid();
	  	const auto rawId = recHitId.rawId();
	  	if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEB_ )
	for (const auto recHit : *recHitsEE_ ){
	  	auto enr = recHit.energy();
	  	if( enr <= minenr ) continue;
	  	const auto recHitId = recHit.detid();
	  	const auto rawId = recHitId.rawId();
	  	if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEE_ )
			
	return result;

}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr = 0.0 )

rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax ){

    rhGroup result;
    vector<uInt> rawIds;

	float dRcut = 0.05;

    int iter = -1;
    for ( const auto superCluster : superClusterGroup ){
		iter++;
		if( phDr[iter] > dRcut ) continue;
		if( phEnergy[iter] != phEnMax ) continue; 
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
            if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClusterGroup )   
    for (const auto recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;

}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax )

rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax, TH1D* rhTimeHist, TH1D* ootHist, TH1D* matchHist ){

    rhGroup result;
    vector<uInt> rawIds;

    float dRcut = 0.05;

    int iter = -1;
	int nRecHits(0);
	int matched(0);
	int nOOT(0);
    for ( const auto superCluster : superClusterGroup ){
        iter++;
        if( phDr[iter] > dRcut ) continue;
        if( phEnergy[iter] != phEnMax ) continue;
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
            if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ){ rawIds.push_back(rawId); nRecHits++; }
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClusterGroup )   
    for (const auto recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() )
                //{ result.push_back(recHit); if(recHit.checkFlag(EcalRecHit::kOutOfTime)) nOOT++; matched++; }
				{ result.push_back(recHit); if(recHit.checkFlag(EcalRecHit::kOutOfTime)) nOOT++; rhTimeHist->Fill(recHit.time()); matched++; }
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() )
				//{ result.push_back(recHit); if(recHit.checkFlag(EcalRecHit::kOutOfTime)) nOOT++; matched++; }
                { result.push_back(recHit); if(recHit.checkFlag(EcalRecHit::kOutOfTime)) nOOT++; rhTimeHist->Fill(recHit.time()); matched++; }
    }//<<>>for (const auto recHit : *recHitsEE_ )
	if( nRecHits > 0 && matched > 0 ){
		matchHist->Fill(matched/nRecHits);
		ootHist->Fill(nOOT/matched);
	} else {
		matchHist->Fill(-0.05);
        ootHist->Fill(-0.05);
	}//<<>>if( nRecHits > 0 )

    return result;

}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax, TH1D ootHist, ... )

rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr, TH1D* rhTimeHist, TH1D* ootHist, TH1D* matchHist ){

    rhGroup result;
    vector<uInt> rawIds;

    int iter = -1;
    float nRecHits(0);
    float matched(0);
    float nOOT(0);
    for ( const auto superCluster : superClusterGroup ){
        iter++;
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
            if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ){ rawIds.push_back(rawId); nRecHits++; }
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClusterGroup )   

//    reverse loop order ? ( loop rawids and find rechit instead of loooping rechits and finding rawids )  : think current method quicker
//    for (const auto rawid : rawIds ){
//		auto recHit = (*recHitsEB_)->find(recHitId);

    for (const auto recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        //if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ){ 
			matched++;
			if(recHit.checkFlag(EcalRecHit::kOutOfTime)) nOOT++;
			if( enr > minenr ){ result.push_back(recHit); rhTimeHist->Fill(recHit.time()); }//<<>>if( enr <= minenr ) 
		}//<<>>if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() )
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        //if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ){ 
            matched++;
            if(recHit.checkFlag(EcalRecHit::kOutOfTime)) nOOT++;
            if( enr > minenr ){ result.push_back(recHit); rhTimeHist->Fill(recHit.time()); }//<<>>if( enr <= minenr ) 
        }//<<>>if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() )
    }//<<>>for (const auto recHit : *recHitsEE_ )
	//std::cout << " --- rh cl stats : nRH " << nRecHits << " matched : " << matched << " nOOT " << nOOT << std::endl;

    if( nRecHits > 0 && matched > 0 ){
        matchHist->Fill(matched/nRecHits);
        ootHist->Fill(nOOT/matched);
	} else if( nRecHits > 0 ){
		matchHist->Fill(matched/nRecHits);
    } else {
        matchHist->Fill(-0.05);
        ootHist->Fill(-0.05);
    }//<<>>if( nRecHits > 0 )

    return result;

}//>>>>rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const scGroup superClusterGroup, float minenr, TH1D ootHist, ... )

rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 ){

	rhGroup result;
	vector<uInt> rawIds;
	auto & hitsAndFractions = basicCluster.hitsAndFractions();
	const auto nHAF = hitsAndFractions.size();
	//std::cout << " --- getRHGroup for basic cluster with " << nHAF << " rechits " << std::endl;
	for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
		const auto detId = hitsAndFractions[iHAF].first;
		const auto rawId = detId.rawId();
		if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
	}//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
	for (const auto recHit : *recHitsEB_ ){
	  	auto enr = recHit.energy();
	  	if( enr <= minenr ) continue;
	  	const auto recHitId = recHit.detid();
	  	const auto rawId = recHitId.rawId();
	  	if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEB_ )
	for (const auto recHit : *recHitsEE_ ){
	  	auto enr = recHit.energy();
	  	if( enr <= minenr ) continue;
	  	const auto recHitId = recHit.detid();
	  	const auto rawId = recHitId.rawId();
	  	if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
	}//<<>>for (const auto recHit : *recHitsEE_ )		

	return result;

}////rhGroup LLPgammaAnalyzer_AOD::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 )

	
vector<float> LLPgammaAnalyzer_AOD::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

	vector<float> result;
	for (const auto recHit : recHits ){		
		const auto rht = recHit.time();
		//std::cout << " ----- Get TOF Time rh time: " << rht << std::endl;
	   	const auto recHitId(recHit.detid());
	   	const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
	   	const auto rhPosX = recHitPos.x();
	   	const auto rhPosY = recHitPos.y();
	   	const auto rhPosZ = recHitPos.z();
		//std::cout << " ----- Get TOF Time rh POS: " << rhPosX << " " <<	rhPosY << " " << rhPosZ << std::endl;
	   	const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
	   	const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
	   	const auto tof = (d_rh-d_pv)/SOL;
	   	//std::cout << " ----- Get TOF Time rh tof: " << tof << std::endl;
	   	result.push_back(rht-tof);
	}//<<>>for (const auto recHit : recHits )	

	return result;

}//>>>>vector<float>  LLPgammaAnalyzer_AOD::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float LLPgammaAnalyzer_AOD::getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ ){

    const auto recHitId(rechit.detid());
    const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
	const auto rhPosX = recHitPos.x();
	const auto rhPosY = recHitPos.y();
	const auto rhPosZ = recHitPos.z();
	const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
	const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
	const auto tof = (d_rh-d_pv)/SOL;
	return tof;	

}//>>>>>float getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ )

EcalRecHit LLPgammaAnalyzer_AOD::getLeadRh( rhGroup recHits ){
	
	EcalRecHit result;
	float enr(0.0);
	for (const auto recHit : recHits ){
		auto rhenr = recHit.energy();
		if( rhenr < enr ) continue;
		enr = rhenr;
		result = recHit;
	}//<<>>for (const auto recHit : recHits )

	return result;

}//>>>>EcalRecHit LLPgammaAnalyzer_AOD::getLeadRh( rhGroup recHits )


vector<float> LLPgammaAnalyzer_AOD::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

	vector<float> result;
	auto lrh = getLeadRh(recHits);
	const auto recHitId(lrh.detid());
	const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
	const auto rhPosX = recHitPos.x();
	const auto rhPosY = recHitPos.y();
	const auto rhPosZ = recHitPos.z();
	const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
	const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
	const auto tof = (d_rh-d_pv)/SOL;
	for (const auto recHit : recHits ){result.push_back(recHit.time()-tof);}
	return result;

}//>>>>vector<float>  LLPgammaAnalyzer_AOD::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float LLPgammaAnalyzer_AOD::getSeedTofTime( reco::SuperCluster sprclstr, double vtxX, double vtxY, double vtxZ ){

	float rhTime(-100.0);
	const auto seedDetId = sprclstr.seed().get()->seed(); // seed detid
	const auto seedRawId = seedDetId.rawId(); // crystal number
    const auto isEB = (seedDetId.subdetId() == EcalSubdetector::EcalBarrel); // which subdet
    const auto recHits = (isEB ? *recHitsEB_ : *recHitsEE_ );
	for( const auto rechit : recHits ){
        const auto recHitId = rechit.detid();
        const auto rawId = recHitId.rawId();
        if( rawId == seedRawId ){ rhTime = rechit.time(); continue; }
    }//<<>>for (const auto recHit : *recHitsEE_ ) 

    const auto recHitPos = barrelGeometry->getGeometry(seedDetId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
	const auto seedTofTime = rhTime - tof;
    return seedTofTime;

}//>>>>float  LLPgammaAnalyzer_AOD::getSeedTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float LLPgammaAnalyzer_AOD::getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    auto lrh = getLeadRh(recHits);
    const auto recHitId(lrh.detid());
    const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    return lrh.time()-tof;

}//>>>>float  LLPgammaAnalyzer_AOD::getSeedTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

vector<float> LLPgammaAnalyzer_AOD::getTimeDistStats( vector<float> times ){

	vector<float> results;
	int size = times.size();
	
	if( size == 0 ){		
	   	std::cout << "Time group is empty" << std::endl; 
	   	results.push_back(-29.25);//0
	   	results.push_back(99.9);//1
	   	results.push_back(-29.25);//2
	   	results.push_back(99.9);//3
	   	results.push_back(99.9);//4
	   	results.push_back(9.9);//5
	   	results.push_back(-29.25);//6
	   	results.push_back(99.9);//7
	   	results.push_back(99.9);//8
	   	results.push_back(99.9);//9
	   	results.push_back(-29.25);//10
	   	results.push_back(99.9);//11
	   	return results;
	}//<<>>if( size == 0 )

	TH1D timeDist( "temp", "temp", 500, -25, 25 );
	for( const auto time : times ) timeDist.Fill(time);
	auto nEntries = timeDist.GetEntries();
	double medtime(0.0);
	if( nEntries > 0 ){	
	   	double quant(0.5);
	   	timeDist.ComputeIntegral();
	   	timeDist.GetQuantiles(1, &medtime, &quant);
	} else medtime = -9.9;

	auto error = timeDist.GetMeanError();
	results.push_back(timeDist.GetMean());//0
	results.push_back(error);//1
	results.push_back(medtime);//2
	results.push_back(1.2533*error);//3 valid for noraml distibutions
	results.push_back(timeDist.GetRMS());//4
	results.push_back(timeDist.GetSkewness());//5

	auto rmu = mean(times);
	results.push_back(rmu);//6 mean
	auto rsd = stdev(times, rmu);
	results.push_back(rsd);//7 stdev
	results.push_back(rms(times));//8 rms
	auto err = rsd/std::sqrt(size);
	results.push_back(err);//9 error of mean

	std::sort(times.begin(),times.end());
	if( size%2 == 0 ) results.push_back((times[(size/2)-1] + times[(size/2)])/2);//10
	else results.push_back(times[size/2]);//10 median
	results.push_back(1.2533*err);//11 error of median
	
	return results;
}//>>>> vector<float> LLPgammaAnalyzer_AOD::getTimeDistStats( vector<float> times )


vector<float> LLPgammaAnalyzer_AOD::getTimeDistStats( vector<float> times, vector<float> wts ){

	vector<float> results;
	int size = times.size();	
	int wtsize = wts.size();
	if( size != wtsize || size == 0 ){ 
		std::cout << "Time & rechit groups not same or empty" << std::endl; 
		results.push_back(-29.25);//0 mu
	   	results.push_back(99.9);//1
	   	results.push_back(-29.25);//2
	   	results.push_back(99.9);//3
	   	results.push_back(99.9);//4
	   	results.push_back(9.9);//5
	   	results.push_back(-29.25);//6
	   	results.push_back(99.9);//7
	   	results.push_back(99.9);//8
	   	results.push_back(99.9);//9
		results.push_back(-29.25);//10
	   	results.push_back(9.9);//11
	   	results.push_back(-10.0);//12
		return results; 
	}//<<>>if( size != rechits.size() || size == 0 )
			
	TH1D timeDist( "temp", "temp", 500, -10, 10 );
	vector<float> wtimes;
	float wtot(0.0);
	for( int it(0); it < size; it++ ){
		timeDist.Fill(times[it], wts[it]);
		wtimes.push_back(times[it]*wts[it]);
		wtot += wts[it];
	}//<<>>for( uInt it(0); it < size; it++ )
   	//auto nEntries = timeDist.GetEntries();
	double hmedtime(9.9);
	//if( nEntries > 0 ){
	//	double quant(0.5);
	// 	timeDist.ComputeIntegral();
	// 	timeDist.GetQuantiles(1, &hmedtime, &quant);
	//} else hmedtime = -9.9; 
	
	auto herror = timeDist.GetMeanError();
	auto hmutime = timeDist.GetMean();
	auto hrms = timeDist.GetRMS();
	results.push_back(hmutime);//0
	results.push_back(herror);//1
	results.push_back(hmedtime);//2
	results.push_back(1.2533*herror);//3 - valid for noraml distibutions
	results.push_back(hrms);//4
	results.push_back(timeDist.GetSkewness());//5

	auto rmu = mean(wtimes,wtot);
	results.push_back(rmu);//6 mean
	auto rsd = stdev(times, rmu, wts, wtot);
	results.push_back(rsd);//7 stdev
	results.push_back(rms(wtimes));//8 rms
	auto err = rsd/std::sqrt(size);
	results.push_back(err);//9 error of mean

	std::sort(times.begin(),times.end());
	if( size%2 == 0 ) results.push_back((times[(size/2)-1]+times[(size/2)])/2);//10
	else results.push_back(times[size/2]);//10 median
	results.push_back(1.2533*err);//11 error of median	
	results.push_back(wtot);//12 tot e of all rechits
	results[2] = results[6] - results[10];

	return results;
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getTimeDistStats( vector<float> times, vector<float> weights )

vector<float> LLPgammaAnalyzer_AOD::getTimeDistStats( vector<float> times, rhGroup rechits ){

   	// N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

   	float N(3.64);
   	float C(0.3000);

   	vector<float> wts;
   	for( uInt it(0); it < rechits.size(); it++ ){
	//auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
      	auto wt = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
		wts.push_back(wt);
   	}//<<>>for( uInt it(0); it < rechits.size(); it++ )

   	return getTimeDistStats( times, wts );
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getTimeDistStats( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen( vector<float> xs, vector<float> wts ){
//spherical

    vector<float> results;

    auto ts2 = wsin2( xs, wts );
    auto tc2 = wcos2( xs, wts );
    auto tsc = wsincos( xs, wts );
    double array[] = { ts2, tsc, tsc, tc2 };
    TMatrixDSym mat(2,array);
    TMatrixDSymEigen eigen(mat);
    const TVectorD& eginVal = eigen.GetEigenValues();
    const TMatrixD& eginVec = eigen.GetEigenVectors();
    TVectorD evaules(eginVal);

    int index(1);
    if( (eginVal(0) >= eginVal(1)) ){
        index = 0;
        //if( (eginVal(0) == eginVal(1)) ) std::cout << " -- rhGrp is Spherical" << std::endl;
    }//<<>>if( (eginVal(0) >= eginVal(1)) )
    //std::cout << "Largest Eigin Value Found" << std::endl;

    auto ex = eginVec(index,0);
    auto ey = eginVec(index,1);
    auto ev = eginVal(index);
    //std::cout << "Eigin Angles Found" << std::endl;

    results.push_back(ex);
    results.push_back(ey);
    results.push_back(ev);

    return results;
}//<<>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen2D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen( vector<float> xs, vector<float> ys, vector<float> wts ){

    vector<float> results;

    auto mean_x = mean( xs, wts );
    auto mean_y = mean( ys, wts );
    auto swts = vfsum( wts );
    auto var_x = var( xs, mean_x, wts, swts );
    auto var_y = var( ys, mean_y, wts, swts );
    auto var_xy = cvar( xs, mean_x, ys, mean_y, wts, swts );
    //std::cout << "Varencies and Means Calculated" << std::endl;

    //TMatrixDSym mat(3,3);
    double array[] = { var_x, var_xy, var_xy, var_y };
    TMatrixDSym mat(2,array);
    //std::cout << "Input matrix created" << std::endl;
    //mat.Print();

    TMatrixDSymEigen eigen(mat);
    const TVectorD& eginVal = eigen.GetEigenValues();
    const TMatrixD& eginVec = eigen.GetEigenVectors();
    TVectorD evaules(eginVal);

    int index(1);
    if( (eginVal(0) >= eginVal(1)) ){
        index = 0;
        //if( (eginVal(0) == eginVal(1)) ) std::cout << " -- rhGrp is Spherical" << std::endl;
    }//<<>>if( (eginVal(0) >= eginVal(1)) )
    //std::cout << "Largest Eigin Value Found" << std::endl;

    auto ex = eginVec(index,0);
    auto ey = eginVec(index,1);
    auto ev = eginVal(index);
    //std::cout << "Eigin Angles Found" << std::endl;

    results.push_back(ex);
    results.push_back(ey);
    results.push_back(ev);

    return results;
}//<<>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen2D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts ){
// ieipt

	vector<float> results;

	auto mean_x = mean( xs, wts );
    auto mean_y = mean( ys, wts );
    auto mean_z = mean( zs, wts );
    auto swts = vfsum( wts );
    auto var_x = var( xs, mean_x, wts, swts );
    auto var_y = var( ys, mean_y, wts, swts );
    auto var_z = var( zs, mean_z, wts, swts );
    auto var_xy = cvar( xs, mean_x, ys, mean_y, wts, swts );
    auto var_xz = cvar( xs, mean_x, zs, mean_z, wts, swts );
    auto var_yz = cvar( ys, mean_y, zs, mean_z, wts, swts );

	//TMatrixDSym mat(3,3);
	double array[] = { var_x, var_xy, var_xz, var_xy, var_y, var_yz, var_xz, var_yz, var_z };
    TMatrixDSym mat(3,array);
    //std::cout << "Input matrix created" << std::endl;
	//mat.Print();

	TMatrixDSymEigen eigen(mat);
	const TVectorD& eginVal = eigen.GetEigenValues();
	const TMatrixD& eginVec = eigen.GetEigenVectors();
	TVectorD evaules(eginVal);	

	int index(2);
	int index2(2);
	if( (eginVal(0) >= eginVal(1)) && (eginVal(0) >= eginVal(2)) ){ 
		index = 0;
		if( eginVal(1) >= eginVal(2) ) index2 = 1; 
		if( (eginVal(0) == eginVal(1)) && (eginVal(0) == eginVal(2)) ) std::cout << " -- rhGrp is Spherical" << std::endl;
    } else if( eginVal(1) >= eginVal(2) ) {
		index = 1;
		if( eginVal(0) >= eginVal(2) ) index2 = 0;
		if( eginVal(1) == eginVal(2) ) std::cout << " -- rhGrp is a Flatend Sphere" << std::endl;
	} else { if( eginVal(0) >= eginVal(1) ) index2 = 0; else index2 = 1; } //<<>>if( (eginVal(0) >= eginVal(1)) && (eginVal(0) >= eginVal(2)) )
    //std::cout << "Largest Eigin Value Found" << std::endl;

	auto index3 = 3 - index - index2;
	auto ev23hypo = hypo(eginVal[index2], eginVal[index3]);
	auto ev1Shypo = hypo(eginVal[index], ev23hypo);
	auto speginval = sq2(eginVal[index]/ev1Shypo);

	auto ex = eginVec(index,0);
    auto ey = eginVec(index,1);
    auto ez = eginVec(index,2);
    auto ev = speginval; //eginVal(index);
    //std::cout << "Eigin Angles Found" << std::endl;

	results.push_back(ex);//0
    results.push_back(ey);//1
    results.push_back(ez);//2
    results.push_back(ev);//3

	return results;
}//<<>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen3D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_ieipt( vector<float> times, rhGroup rechits ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);
    float tdis(2.2/SOL);

    vector<float> wts;
    vector<float> etas;
    vector<float> phis;
    vector<float> ebtimes;
    vector<float> emepReturn(8,-9);

	float totWts(0);
    for( uInt it(0); it < rechits.size(); it++ ){

		const auto recHitId = getRawID(rechits[it]);
		//const auto & idinfo = DetIDMap[recHitID];
		//auto isEB = getIsEB(rechits[it]);
		//if( idinfo.ecal == ECAL::EB ){
        if( getIsEB(rechits[it]) ){
    		//auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
        	auto resolution = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
        	wts.push_back(resolution);
			totWts += resolution;
        	const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        	//const auto rhEtaPos = tdis*(idinfo.i2);//recHitPos.eta();
        	//etas.push_back(rhEtaPos);
            etas.push_back(recHitPos.z());
        	//const auto rhPhiPos = tdis*idinfo.i1;//recHitPos.phi();
        	//phis.push_back(tdis*(idinfo.i1);// tdis * below with wrap correction
			float rawPhi = PI+recHitPos.phi();
            phis.push_back(tdis*rawPhi/0.01745329);
			ebtimes.push_back(times[it]);
		} else { return emepReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
		}//<<>>if( idinfo.ecal == ECAL::EB )

    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

	vector<float> adjphis;
	auto mphi = ceil(meanIPhi(phis));
	for( uInt it(0); it < phis.size(); it++ ){ 
		auto dstphi = dIPhi(phis[it],mphi);
		adjphis.push_back(tdis*(mphi+dstphi));
	}//<<>>for( auto it : phis )

    auto eigens = getRhGrpEigen( etas, adjphis, ebtimes, wts );//0 x, 1 y, 2 z, 3 values

    vector<float> angles;
    auto meta = mean(etas,wts);
    auto adjmphi = mean(adjphis,wts);
    auto mtime = mean(ebtimes,wts);

    auto eprotangle = getAngle(eigens[0], eigens[1]);
    float epeignsin = std::sin(eprotangle);
    float epeigncos = std::cos(eprotangle);
    vector<float> xs;
    float xsum(0.0);
    for( uInt it(0); it < ebtimes.size(); it++ ){

        float leta = etas[it] - meta;
        float lphi = adjphis[it] - adjmphi;
        float ltime = ebtimes[it] - mtime;
        auto epxcor = epeigncos*(leta) - epeignsin*(lphi);
		auto sl = (ltime)/(epxcor);
		xsum += sl*wts[it];
        xs.push_back(sl);

    }//<<>>for( uInt it(0); it < wts.size(); it++ )

	auto slope = xsum/totWts;
	auto nWts = wts.size();

    auto varsl = var(xs,slope,wts,totWts);
    auto chi2 = chisq(xs,slope,varsl);
    auto chi2pf = 1 - TMath::Prob(chi2, nWts);

    eigens.push_back(slope);//4
    eigens.push_back(chi2pf);//5

    return eigens;
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_xyz( vector<float> times, rhGroup rechits ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);

    vector<float> wts;
    vector<float> etas;
    vector<float> phis;
    for( uInt it(0); it < rechits.size(); it++ ){
    //auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
        auto wt = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
        wts.push_back(wt);
        const auto recHitId(rechits[it].detid());
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        const auto rhEtaPos = recHitPos.eta();
        etas.push_back(rhEtaPos);
        const auto rhPhiPos = recHitPos.phi();
        phis.push_back(rhPhiPos);
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    return getRhGrpEigen( etas, phis, times, wts );
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_sph( vector<float> times, rhGroup rechits ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);
    //float twoPI(2*PI);

    vector<float> wts;
    vector<float> etas;
    vector<float> phis;
    vector<float> ebtimes;
    vector<float> angles;
    //vector<float> energies;
    vector<float> resolutions;
    //float sumenergy = 0;
    vector<float> emepReturn(8,-9);

	//bool ebp = true;
    //bool ebn = true;
	auto nRecHits = rechits.size();
    for( uInt it(0); it < nRecHits; it++ ){

        const auto recHitId = getRawID(rechits[it]);
		const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        //const auto & idinfo = DetIDMap[recHitID];
        //if( idinfo.ecal == ECAL::EB ){ 
        if( getIsEB(rechits[it]) ){
        	//const auto rhEtaPos = idinfo.i2;//recHitPos.eta();
        	//etas.push_back(rhEtaPos);
            etas.push_back(recHitPos.z());
        	//const auto rhPhiPos = idinfo.i1;//recHitPos.phi();
        	//phis.push_back(rhPhiPos);
            phis.push_back(2.2*(PI+recHitPos.phi())/0.01745329);
        	ebtimes.push_back(times[it]);	
			auto rhenergy = rechits[it].energy();
			auto resolution = 1/(sq2(N/rhenergy)+2*C*C);
            resolutions.push_back(resolution);
			//std::cout << "In getRhGrpEigen_sph w/ rheta " << rhEtaPos << " : rhphi " << rhPhiPos << " : rht " << times[it] << std::endl;		
		} else { return emepReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
		}//<<>>if( idinfo.ecal != ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    auto meta = mean(etas,resolutions);
    auto mphi = meanIPhi(phis,resolutions);
	auto mtime = mean(ebtimes,resolutions);

	//std::cout << "In getRhGrpEigen_sph w/ meta " << meta << " : mphi " << mphi << " : mt " << mtime << std::endl;

	for( uInt it(0); it < ebtimes.size(); it++ ){

		float leta = etas[it]-meta;
        float lphi = dIPhi(phis[it],mphi);
		if( leta == 0 && lphi == 0 ) continue;
        float ltim = ebtimes[it]-mtime;
		float angle = getAngle( leta, lphi );
		angles.push_back(angle);
	//std::cout << "In getRhGrpEigen_sph w/1 leta " << leta << " : lphi " 
	//			<< lphi << " : ltim " << ltim << " : mag " << resolutions[it] << " : angle " << angle << std::endl;
        wts.push_back(ltim*ltim*resolutions[it]);
	}//<<>>for( uInt it(0); it < ebtimes.size(); it++ )

    auto eigens =  getRhGrpEigen( angles, wts );//0 x, 1 y, 2 values
    auto rotangle = getAngle(eigens[0], eigens[1]);
	float eignsin = std::sin(rotangle);
	float eigncos = std::cos(rotangle);
	auto orgsin = eignsin;
	auto orgcos = eigncos;

	// ----------------------
	// aligning eigen to point same direction
	// ----------------------
    float totRes(0.0);
	float ltsum(0.0);
	int slopeCorr(1);
	for( uInt it(0); it < wts.size(); it++ ){

        float leta = etas[it] - meta;
        float lphi = dIPhi(phis[it],mphi);
		float ltime = ebtimes[it]-mtime;
		auto xcor = orgcos*(leta) - orgsin*(lphi);
		//if( xcor > 0 ) ltsum += ltime*resolutions[it];
		ltsum += ltime*resolutions[it]/xcor;	
        totRes += resolutions[it];	

	}//<<>>for( uInt it(0); it < wts.size(); it++ )

	if( ltsum/totRes < 0 ){ 

    	rotangle = getAngle(-1*eigens[0], -1*eigens[1]);
		eignsin = std::sin(rotangle);
		eigncos = std::cos(rotangle);
		slopeCorr = -1;

	}//if( ebp && ebn && ltsum < 0 )


	// -----------------------------------------
	// finding nemo ( slope )
	// -----------------------------------------
    auto nWts = resolutions.size();

	vector<float> xs;
    vector<float> xs2;
	float xsum(0.0);
    float xsum2(0.0);
    for( uInt it(0); it < nWts; it++ ){

        float leta = etas[it] - meta;
        float lphi = dIPhi(phis[it],mphi);
        float ltime = ebtimes[it]-mtime;
        auto xcor = eigncos*(leta*2.2) - eignsin*(lphi*2.2);
        auto ycor = eignsin*(leta*2.2) + eigncos*(lphi*2.2);
        auto sxcor = (orgcos*(leta*2.2) - orgsin*(lphi*2.2))*slopeCorr;
        auto sycor = orgsin*(leta*2.2) + orgcos*(lphi*2.2);
        if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << leta << " : lphi " << lphi 
								<< " : xcor " << xcor << " : ycor " << ycor << " : dt " << wts[it] << std::endl;
        if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << leta << " : lphi " << lphi 
								<< " : sxcor " << sxcor << " : sycor " << sycor << " : dt " << wts[it] << std::endl;
        //auto sl = (ltime)/(sxcor);
        auto sl = (ltime)/(sxcor)*slopeCorr;
		//auto sl2 = (ltime)/(leta*2.2);// changed to eta to look at slope in eta direction
        auto sl2 = (ltime)/(leta*2.2)*slopeCorr;
		xs.push_back(sl);
        xs2.push_back(sl2);
        xsum += sl*resolutions[it];
		xsum2 += sl2*resolutions[it]; 
    }//<<>>for( uInt it(0); it < wts.size(); it++ )

	//float xsum(0); for( uInt it(0); it < nWts; it++ ){ xsum += xs[it]*resolutions[it]; //hist1d[100]->Fill(xs[it]);}
	auto nXSum = xs.size();
	auto slope = xsum/totRes;
    auto slope2 = xsum2/totRes;
	auto varsl = var(xs,slope,resolutions,totRes);
    //auto varsl2 = var(xs2,slope2,resolutions,totRes);
	auto chi2 = chisq(xs,slope,varsl);
    auto chi22 = chisq(xs,slope,varsl);
    auto chi2pf = 1 - TMath::Prob(chi2, nWts);
    auto chi2pf2 = 1 - TMath::Prob(chi22, nWts);
	eigens.push_back(slope);//3
    eigens.push_back(chi2pf);//4
    eigens.push_back(slope2);//5
    eigens.push_back(chi2pf2);//6
    eigens.push_back(rotangle);//7
    eigens.push_back(nXSum);//8

	return eigens;
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_sph( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_ep( vector<float> times, rhGroup rechits ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);
    float tdis(2.2/SOL);

    vector<float> wts;
    vector<float> etas;
    vector<float> phis;
    vector<float> ebtimes;

    for( uInt it(0); it < rechits.size(); it++ ){

        const auto recHitId = getRawID(rechits[it]);
		const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        //const auto & idinfo = DetIDMap[recHitID];
        //if( idinfo.ecal == ECAL::EB ){
        if( getIsEB(rechits[it]) ){
            //auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
            auto wt = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
            wts.push_back(wt);
            //const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
            //const auto rhEtaPos = tdis*idinfo.i2;//recHitPos.eta();
            //etas.push_back(rhEtaPos);
            etas.push_back(recHitPos.z());
			////hist1d[69]->Fill(rhEtaPos);
            //const auto rhPhiPos = tdis*idinfo.i1;//recHitPos.phi();
            ////hist1d[70]->Fill(rhPhiPos);
            //phis.push_back(rhPhiPos);
			float rawPhi = PI+recHitPos.phi();
            phis.push_back(tdis*(rawPhi)/0.01745329);
            ebtimes.push_back(times[it]);
        }//<<>>if( idinfo.ecal == ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    return getRhGrpEigen( etas, phis, wts );
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_ep( vector<float> times, rhGroup rechits ){

float LLPgammaAnalyzer_AOD::getdt( float t1, float t2 ){

	auto dt = t1 - t2;
	if( dt == 0.0 )	dt = -5.5;// to be commented out in final version
	if( t1 == 0.0 || t2 == 0.0 ) dt = -5.0;// to be copmmented out in final version
	if( t1 < -28.9 || t2 < -28.9 ) dt = -3.75;
	
	return dt;
}//>>>>float LLPgammaAnalyzer_AOD::getdt( float t1, float t2 )


void LLPgammaAnalyzer_AOD::mrgRhGrp( rhGroup & x, rhGroup & y ){

	//std::cout << " --- In mrgRhGrp " << std::endl;
	bool matched(false);
	for(const auto rhy : y ){ 
		for(const auto rhx : x ){ if( rhMatch(rhx,rhy) ) matched = true; }
		if( not matched ) x.push_back(rhy);
	}//<<>>for(const auto rhx : y )
   //std::cout << " ---- Finished mrgRhGrp " << std::endl;

	return;
}//>>>>rhGroup LLPgammaAnalyzer_AOD::mrgRhGrp(const rhGroup x, const rhGroup y)


bool LLPgammaAnalyzer_AOD::reduceRhGrps( vector<rhGroup> & x ){

	//std::cout << " --- In reduceRhGrps " << std::endl;
	uInt s = x.size();
	uInt a(0), b(0);
	bool match(false);
   	for( uInt ita(0); ita < s; ita++ ){ 
		for( auto itb = ita+1; itb < s; itb++ ){ 
			if( dupRhFnd(x[ita], x[itb]) ){ a = ita; b = itb; match = true; } 
			if(match) break;
   	}//<<>>for( int itb = ita+1; itb < s; itb++ ) 
		if(match) break;
	}//<<>>for( int ita = 0; ita < nBcRhGroups; ita++ )
	if(match){
		mrgRhGrp( x[a], x[b] );
		x.erase(x.begin()+b);
	} //<<>>if(match) : else
	//std::cout << " ---- Finished reduceRhGrps " << std::endl;

	return match;
}//>>>>vector<rhGroup> LLPgammaAnalyzer_AOD::reduceRhGrps(const vector<rhGroup> x)

vector<float> LLPgammaAnalyzer_AOD::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  ){

    vector<float> result;
    vector<float> kidtime;
    vector<float> kide;
	vector<float> impact;
	bool llpcut(false);
	float nextBX(0.0);
	float llpkid(0.0);
    for( auto kid : kids ){

        bool done(false);
        bool first(true);
        bool llp(false);
        bool bquark(false);
        auto mom = kid->mother();
        vector<double> stepp{kid->p()};
        vector<double> stepe{kid->energy()};
        vector<double> orignvx;
        vector<double> orignvy;
        vector<double> orignvz;
        bool top(true);
		int steps(0);
        while( not done ){
			auto pdgid = std::abs(mom->pdgId());
			if( pdgid == 6000113 ) llp = true;
			if( pdgid == 5 ) bquark = true;
			if( std::floor(pdgid/100) == 5 ) bquark = true;
            if( mom->pt() == 0 ) done = true;
            else {
                auto gmom = mom->mother();
                if( top ){
                    if( first ) first = false;
                    else {
                        stepp.push_back(mom->p());
                        stepe.push_back(mom->energy());
                    }//<<>>if( first )
					steps++;
					//if( llp || bquark ) llpstep.push_back(1.0); else llpstep.push_back(0.0);
                    orignvx.push_back(mom->vx());
                    orignvy.push_back(mom->vy());
                    orignvz.push_back(mom->vz());
                }//<<>>if( top )
                if( gmom->pt() == 0 ) done = true;
                else {
                    if( mom->vx() == gmom->vx() ) top = false;
                    else top = true;
                    mom = gmom;
                }//<<>>if( gmom->pt() == 0 )
            }//<<>>if( mom->pt() == 0 )                 
        }//<<>>while( not done )
        first = true;

        if( DEBUG && false ) std::cout << " ---- jetGenTime Calc Steps : " << steps << std::endl;
		bool stepcut( steps == 0 );
		if( DEBUG && false ) std::cout << " ---- jetGenTime Calc llpcut : " << llpcut << std::endl;
		if( not stepcut ){
			float maxe(0.0);
            float totaltime(0.0);
            auto destvx(cx);
            auto destvy(cy);
            auto destvz(cz);
            for( int it(0); it < steps; it++ ){
				if( stepe[it] > maxe ) maxe = stepe[it];
                auto beta = stepp[it]/stepe[it];
                auto legtime = hypo( (destvx-orignvx[it]), (destvy-orignvy[it]), (destvz-orignvz[it]) )/(SOL*beta);
				totaltime += legtime;
				//std::cout << "Sum path legs beta: " << beta << " x: " << (destvx-orignvx[it]) << " y: " << (destvy-orignvy[it]);
				//std::cout << " z: " << (destvz-orignvz[it]) << " legtime: " << legtime << " totaltime: " << totaltime << std::endl;
                //std::cout << " -- dest : " << destvx << ", " << destvy << ", " << destvz; 
                //std::cout << " orig : " << orignvx[it] << ", " << orignvy[it] << ", " << orignvz[it] << std::endl;
                destvx = orignvx[it];
                destvy = orignvy[it];
                destvz = orignvz[it];
            }//<<>>for( int it(0); it < steps; it++ )
			auto r = hypo((cx-orignvx[0]), (cy-orignvy[0]));
			auto iangle = getAngle((cz-orignvz[0]),r);	
			//std::cout << "getAngle: " << iangle << " r: " << r << " z: " << cz-orignvz[0] << " x: " << cx-orignvx[0];
			//std::cout << " y: " << cy-orignvy[0] << std::endl; 
			if( totaltime > 25.0 ) {
				nextBX = 1.0;
			} else {
				impact.push_back(iangle);
            	kidtime.push_back(totaltime);
				kide.push_back(maxe);
				llpcut = llpcut || llp || bquark;
				if( llp || bquark ) llpkid++;
			}//<<>>if( totaltime > 25.0 )
        }//<<>>if( steps == 0 )

    }//<<>>for( auto kid : kids )

	auto mtime(-30.0);
	auto mangle(7.0);
	auto tvar(-2.0);
	auto ktsize = kidtime.size();
	//if( llpcut && ( ktsize > 0 ) ){
    if( ktsize > 0 ){
		if( ktsize > 1 ){
			mtime = mean(kidtime,kide);
			mangle = mean(impact,kide);
			tvar = var(kidtime,mtime,kide);
		} else {
			mtime = kidtime[0];
			mangle = impact[0];
			tvar = -1.0;
		}//<<>>if( kidtime.size() > 1 )
	}//<<>>if( llpcut )
	result.push_back(mtime);//0
	result.push_back(mangle);//1
	result.push_back(tvar);//2
	result.push_back(nextBX);//3
	if( llpcut ) result.push_back(1.0); else result.push_back(0.0);//4
	result.push_back(llpkid/ktsize);//5
    result.push_back(ktsize);//6
    return result;

}//>>>>vector<float> LLPgammaAnalyzer_AOD::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  )

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ method called for each event	------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void LLPgammaAnalyzer_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	using namespace edm;

// -- Consume Tokens --------------------------------------------
	if( DEBUG ) std::cout << "Consume Tokens -------------------------------------------- " << std::endl;

	// TRIGGER
	//iEvent.getByToken(triggerResultsToken_,triggerResults_);
	//iEvent.getByToken(triggerObjectsToken_,triggerObjects_);

	// TRACKS
	iEvent.getByToken(tracksToken_, tracks_);

	// PFCANDIDATES
	iEvent.getByToken(pfcand_token_, pfcands_);
	iEvent.getByToken(pfCan_token_, pfCans_);
	iEvent.getByToken(pfCanPhotonMap_token_, pfCanPhotonMap_);
    iEvent.getByToken(pfCanOOTPhotonMap_token_, pfCanOOTPhotonMap_);
    iEvent.getByToken(pfCanElectronMap_token_, pfCanElectronMap_);

	// VERTICES
	iEvent.getByToken(verticesToken_, vertices_);

	// RHO
	//iEvent.getByToken(rhoToken_, rho_);

	// METS
	iEvent.getByToken(pfmetsToken_, pfmets_);

	// SUPERCLUSTERS
	iEvent.getByToken(scToken_, superCluster_);  
	iEvent.getByToken(ootScToken_, ootSuperCluster_);

	// CALOCLUSTERS
	iEvent.getByToken(ccToken_, caloCluster_);

	// JETS
	iEvent.getByToken(jetsToken_, jets_);

    // CALOJETS
    iEvent.getByToken(caloJetsToken_, caloJets_);

	// GENPARTICLES & GENJETS
	if( hasGenInfo ){

		iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);
    	iEvent.getByToken(gent0Token_, gent0_);
    	iEvent.getByToken(genxyz0Token_, genxyz0_);
    	iEvent.getByToken(pileupInfosToken_, pileupInfos_);
    	iEvent.getByToken(genParticlesToken_, genParticles_);
    	iEvent.getByToken(genJetsToken_, genJets_);

	}//<<>>if( hasGenInfo )

	// ELECTRONS
	iEvent.getByToken(electronsToken_, electrons_);

	// PHOTONS
	iEvent.getByToken(gedPhotonsToken_, gedPhotons_);
    iEvent.getByToken(phoCBIDLooseMapToken_, phoCBIDLooseMap_);
	iEvent.getByToken(ootPhotonsToken_, ootPhotons_);

	// MUONS
    iEvent.getByToken(muonsToken_, muons_);

	// ECAL RECHITS
	iEvent.getByToken(recHitsEBToken_, recHitsEB_);
	iEvent.getByToken(recHitsEEToken_, recHitsEE_);

	// GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
	iSetup.get<CaloGeometryRecord>().get(caloGeo_); 
	barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
	endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap); 

  	// Laser constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	iSetup.get<EcalLaserDbRecord>().get(laserH);
	evTime = iEvent.time();

  	// Intercalibration constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	iSetup.get<EcalIntercalibConstantsRcd>().get(interCalibH);
	interCalibMap = &interCalibH->getMap();  

  	// ADCToGeV : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	iSetup.get<EcalADCToGeVConstantRcd>().get(adcToGeVH);
	adcToGeVEB = adcToGeVH->getEBValue();
	adcToGeVEE = adcToGeVH->getEEValue();

  	// Pedestals : https://github.com/ferriff/usercode/blob/master/DBDump/plugins/DBDump.cc
  	iSetup.get<EcalPedestalsRcd>().get(pedestalsH);

	// PFCand refs
	//iEvent.getByToken(reco2pf_, reco2pf);

// -- event information ------------------------------------------------------------------------------------------------------------------------

   	run   = iEvent.id().run();
   	lumi  = iEvent.luminosityBlock();
   	event = iEvent.id().event();
    if( DEBUG ) std::cout << "******************************************************************************************************" << std::endl;
	if( DEBUG ) std::cout << "Processing event: " << event << " in run: " << run << " and lumiblock: " << lumi << std::endl;

// -- Process Prime Vertix -----------------------------------
	const auto & primevtx = vertices_->front();
	
	vtxX = primevtx.position().x();
	vtxY = primevtx.position().y();
	vtxZ = primevtx.position().z();
	
// -- Process Objects ---------------------------------------    

// ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
// ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 

	std::vector<reco::PFJet> 		fjets;
  	std::vector<int> 				fjetsID;
  	std::vector<reco::CaloCluster> 	fbclusts;
  	std::vector<EcalRecHit> 		frechits;
  	std::vector<reco::Photon> 		fphotons;
  	std::vector<reco::Photon> 		footphotons;
  	std::vector<reco::GsfElectron> 	felectrons;
  	std::vector<reco::CaloJet> 		fcalojets;
	std::vector<reco::Muon> 		fmuons;
	std::vector<reco::PFMET> 		fpfmet;
	std::vector<reco::GenParticle>  fgenparts;
  	std::vector<int>				fgenpartllp;

	//jets.clear(); 
	//jets.reserve(jets_->size());

	// Selection criteria declartation
	// deltaRmin from ecalelf SC Associtor : 0.2
	float	deltaRminJet	= 0.4;//0.4
	//float deltaRminKid	= 0.15;//0.2
	//float 	jetPTmin		= 200.0;
    float   jetPTmin        = 0.0;// for energy/time comp
	int 	jetIDmin		= 2; //3;
    //int     jetIDmin        = 1; //3; for energy/time comp
	float 	jetETAmax		= 1.5; //1.5;
	uInt 	minRHcnt		= 5; //32;
    uInt    minObjRHcnt     = 3; //32;
	float 	minRHenr		= 2.0;
    //float   minRHenr        = 0.0;
   	float 	bcMinEnergy		= 0.667;
    //float   bcMinEnergy     = 0.0;
   	uInt 	bcMinRHGrpSize	= 3;
	float 	minEmf			= 0.0;//0.2

	// -----------------------------------------------------------------------------------------------
	// ---  Collecting objests ( preprocessing object pruning ) ---------------------------------------
	// -------------------------------------------------------------------------------------------------

    if( DEBUG ) std::cout << "Collecting ECAL RecHits" << std::endl;
    for( const auto recHit : *recHitsEB_ ){ frechits.push_back(recHit); }
    for( const auto recHit : *recHitsEE_ ){ frechits.push_back(recHit); }

    if( DEBUG ) std::cout << "Collecting Muons" << std::endl;
    for( const auto muon : *muons_ ){ fmuons.push_back(muon); }

    if( DEBUG ) std::cout << "Collecting Met" << std::endl;
    for( const auto met : *pfmets_ ){ fpfmet.push_back(met); }

    if( DEBUG ) std::cout << "Collecting caloClusters" << std::endl;
    for( const auto &bclust : *caloCluster_ ){ fbclusts.push_back(bclust); }

	if( DEBUG ) std::cout << "Filter Jets" << std::endl;
	for(const auto &jet : *jets_ ){ // Filters jet collection & sorts by pt

		if (jet.pt() < jetPTmin) continue;
	  	if (std::abs(jet.eta()) > jetETAmax) continue;
	  
	  	const auto jetID = getPFJetID(jet);
	  	if (jetID < jetIDmin) continue;

	  	// save the jets, and then store the ID
	  	fjets.emplace_back(jet);
	  	//fjets.back().addUserInt("jetID",jetID);
		fjetsID.push_back(jetID);	  

	  	std::sort(fjets.begin(),fjets.end(),sortByPt);
	}//<<>>for(const auto& jet : *jets_ )

	//genjets	
	//for(const auto& genjet : *genJets_ ){}

	//calojets
	if( DEBUG ) std::cout << "Collecting CaloJets" << std::endl;
	for(const auto &calojet : *caloJets_ ){ fcalojets.push_back(calojet);  }

	if( DEBUG ) std::cout << "Collecting Calo Clusters" << std::endl;
	for( const auto &bclust : *caloCluster_ ){ fbclusts.push_back(bclust); }

    if( DEBUG ) std::cout << "Collecting Photons/Electrons" << std::endl;


    // if( DEBUG ) std::cout << "Finding oot/ged photon matches " << std::endl;
    //----------------- oot/ged photon matching ------------------

    phoExcluded.clear();
    phoIsOotPho.clear();

    for( const auto &ootPho : *ootPhotons_ ){
        double minDr(0.5);
        double dRmatch(10.0);
        float matchpt(0);
        auto oEta = ootPho.eta();
        auto oPhi = ootPho.phi();
		auto oPt = ootPho.pt();
        for( const auto gedPho : *gedPhotons_ ){
            auto pEta = gedPho.eta();
            auto pPhi = gedPho.phi();
            auto pPt = gedPho.pt();
            dRmatch = deltaR( pEta, oEta, pPhi, oPhi );
            if( dRmatch < minDr ){ minDr = dRmatch; matchpt = pPt; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
		fphotons.push_back(ootPho);
		phoIsOotPho.push_back(true);
        if( dRmatch < 0.1 && oPt < matchpt ) phoExcluded.push_back(true); 
		else phoExcluded.push_back(false);
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )

    for( const auto gedPho : *gedPhotons_ ){
        double minDr(0.5);
        double dRmatch(10.0);
        float matchpt(0);
		auto pEta = gedPho.eta();
		auto pPhi = gedPho.phi();
		auto pPt = gedPho.pt();
        for( const auto &ootPho : *ootPhotons_ ){
			auto oEta = ootPho.eta();
			auto oPhi = ootPho.phi();
			auto oPt = ootPho.pt(); 
            dRmatch = deltaR( pEta, oEta, pPhi, oPhi );
            if( dRmatch < minDr ){ minDr = dRmatch; matchpt = oPt; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        fphotons.push_back(gedPho);
        phoIsOotPho.push_back(false);
        if( dRmatch < 0.1 && pPt < matchpt ) phoExcluded.push_back(true);
        else phoExcluded.push_back(false); 
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )


    for( const auto electron : *electrons_ ){ felectrons.push_back(electron); }

	// !!!!!!!  collect gen particles
    //if( DEBUG ) 
	std::cout << "Collecting Gen Particles" << std::endl;
    //auto fgenXYZ0 = *genxyz0_;
    //auto fgenT0 = *gent0_;
    int nGenPart(0);
    std::cout << "----------------------------------------------------------------------------------------------" << std::endl;
/*
    for( const auto genPart : *genParticles_ ){

        nGenPart++;
		int nMoms = genPart.numberOfMothers();
        std::cout << " GenPart " << nGenPart << " : " << genPart.isLastCopy() << " : " << genPart.pdgId();
        std::cout << " M " << nMoms << " A " << &genPart << std::endl;

    }//<<>>for(const auto& genPart : *genParticles_ )
*/
    std::cout << "----------------------------------------------------------------------------------------------" << std::endl;  
    nGenPart = 0;
	for( const auto genPart : *genParticles_ ){ 

		nGenPart++;
		if( genPart.status() == 1 ){
		//if( genPart.isLastCopy() ){
			fgenparts.push_back(genPart);
			auto genPartID = llpGenChaseP(genPart,0);
/*
			auto genPartPdgID = genPart.pdgId();
			bool match(false);
			auto nGenMom = genPart.numberOfMothers();
			if( nGenMom ){
				auto genmom = genPart.mother(0);
				auto momID = genmom->pdgId(); 
				match = ( momID > 1000020 ) && ( std::abs(momID) < 1000040 );
			}//<<>>if( nGenMom )
			if( genPartPdgID == 22 && match ){ 
				std::cout << " ------- Matched Ino Photon: " << genPartPdgID << " ID: " << genPartID << std::endl;
			// }//<<>>if( genPart.pdgId() == 22 ){
			} else if( genPartID < 400 ) std::cout << " --------- Other GenPart: " << genPartPdgID << " ID: " <<genPartID << std::endl;
*/
			fgenpartllp.push_back(genPartID);
		}//<<>>if( genPart.isLastCopy() )

	}//<<>>for(const auto& genPart : *genParticles_ )
    std::cout << "Finished Collecting Gen Particles" << std::endl;

    //------------------------------------------------------------------------------------
	// ----   Object processing ----------------------------------------------------------
    //------------------------------------------------------------------------------------

    if( DEBUG ) std::cout << "Processing GenParticles" << std::endl;

	nGenParts = 0;
	genPt.clear();
	genEnergy.clear();
	genPhi.clear();
	genEta.clear();
	genPx.clear();
	genPy.clear();
	genPz.clear();
	genPdgId.clear();
	genLLP.clear();

    if( DEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
	auto ngens(0);
    for (const auto genpart : fgenparts ){
		
		ngens++;
		genPt.push_back(genpart.pt());
		genEnergy.push_back(genpart.energy());
        genPhi.push_back(genpart.phi());
        genEta.push_back(genpart.eta());
        genPdgId.push_back(genpart.pdgId());
        genLLP.push_back(fgenpartllp[ngens-1]);
        genPx.push_back(genpart.px());
        genPy.push_back(genpart.py());
        genPz.push_back(genpart.pz());

	}//<<>> for (const auto genpart : fgenparts )
	nGenParts = ngens;

    if( DEBUG ) std::cout << "Processing RecHits" << std::endl;

    auto nRecHitCnt(0);
    nRecHits = 0;
    rhPosX.clear();
    rhPosY.clear();
    rhPosZ.clear();
    rhPosEta.clear();
    rhPosPhi.clear();
    rhID.clear();
    rhXtalI1.clear();
    rhXtalI2.clear();
    rhSubdet.clear();
    rhEnergy.clear();
    rhTime.clear();
    rhTimeErr.clear();
    rhTOF.clear();
    rhisOOT.clear();
    rhisGS6.clear();
    rhisGS1.clear();
    rhisWeird.clear();
    rhisDiWeird.clear();
    rhadcToGeV.clear();
    rhSwCross.clear();
    rhped12.clear();
    rhped6.clear();
    rhped1.clear();
    rhpedrms12.clear();
    rhpedrms6.clear();
    rhpedrms1.clear();


    if( DEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
    for (const auto recHit : frechits ){

        //if( DEBUG ) std::cout << " -- proccesing ID info" << std::endl;
        // something in this section is seg faluting after several rechits for crab jobs
        const auto recHitID = getRawID(recHit);
        auto isEB = getIsEB(recHit); // which subdet
        //const auto & idinfo = DetIDMap[recHitID];
	    //if( DEBUG ) std::cout << " -- proccesing EBEE info" << std::endl;
        //if( DEBUG ) std::cout << " -- proccesing GEO info" << std::endl;
        const auto geometry( isEB ? barrelGeometry : endcapGeometry );
        //const auto geometry( ( idinfo.ecal == ECAL::EB ) ? barrelGeometry : endcapGeometry );
        auto recHitPos = geometry->getGeometry(recHit.detid())->getPosition();
        //if( DEBUG ) std::cout << " -- proccesing POSITION info" << std::endl;
        const auto rhX = recHitPos.x();
        const auto rhY = recHitPos.y();
        const auto rhZ = recHitPos.z();
        //if( DEBUG ) std::cout << " -- proccesing TOF info" << std::endl;
        const auto d_rh = hypo(rhX,rhY,rhZ);
        const auto d_pv = hypo(rhX-vtxX,rhY-vtxY,rhZ-vtxZ);
        const auto tof = (d_rh-d_pv)/SOL;
        //if( DEBUG ) std::cout << " -- proccesing SWISSCROSS info" << std::endl;
        float swisscross(0.0);
        if( isEB ) swisscross = EcalTools::swissCross(recHitID, *recHitsEB_, 0.0, true);
        else swisscross = EcalTools::swissCross(recHitID, *recHitsEE_, 0.0, true);

        //if( DEBUG ) std::cout << " -- proccesing LASER info" << std::endl;
        // adcToGeVInfo : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#0204
        const auto laser = laserH->getLaserCorrection(recHitID,evTime);
        const auto interCalibIter = interCalibMap->find(recHitID);
        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        //if( DEBUG ) std::cout << " -- proccesing ADC info" << std::endl;
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) rhadcToGeV[pos] = (laser*interCalib*adcToGeV);
        const float adcToGeV( isEB ? adcToGeVEB : adcToGeVEE );
        //if( DEBUG ) std::cout << " -- proccesing PED info" << std::endl;
        // pedestal info
        const auto & pediter = pedestalsH->find(recHitID);

        //if( DEBUG ) std::cout << " -- storing values BASE" << std::endl;
        rhID.push_back(recHitID);
        rhPosX.push_back(rhX);
        rhPosY.push_back(rhY);
        rhPosZ.push_back(rhZ);
        rhTOF.push_back(tof);
        rhPosEta.push_back(recHitPos.eta());
        rhPosPhi.push_back(recHitPos.phi());
        rhTime.push_back(recHit.time());
        rhTimeErr.push_back(recHit.timeError());
        //if( DEBUG ) std::cout << " -- storing values FLAGS" << std::endl;
        //rhSubdet.push_back( isEB ? 0 : ( ( idinfo.ecal == ECAL::EP ) ? 1 : 2 ) );
        //rhXtalI1.push_back(idinfo.i1);
        //rhXtalI2.push_back(idinfo.i2);
        rhisOOT.push_back(recHit.checkFlag(EcalRecHit::kOutOfTime));
        rhEnergy.push_back(recHit.energy());
        //energyError()
        rhSwCross.push_back(swisscross);
        rhisWeird.push_back(recHit.checkFlag(EcalRecHit::kWeird));
        rhisDiWeird.push_back(recHit.checkFlag(EcalRecHit::kDiWeird));
        rhisGS6.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain6));
        rhisGS1.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain1));
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) 
        //if( DEBUG ) std::cout << " -- storing values PED" << std::endl;
        rhadcToGeV.push_back(laser*interCalib*adcToGeV);
        //else rhadcToGeV.push_back(0.0);
        if (pediter != pedestalsH->end()){
            const auto & ped = (*pediter);
            rhped12.push_back(ped.mean(1));
            rhped6.push_back(ped.mean(2));
            rhped1.push_back(ped.mean(3));
            rhpedrms12.push_back(ped.rms(1));
            rhpedrms6.push_back(ped.rms(2));
            rhpedrms1.push_back(ped.rms(3));
        } else {
            rhped12.push_back(0.0);
            rhped6.push_back(0.0);
            rhped1.push_back(0.0);
            rhpedrms12.push_back(0.0);
            rhpedrms6.push_back(0.0);
            rhpedrms1.push_back(0.0);
        }//<<>>if (pediter != pedestalsH->end())
        nRecHitCnt++;
        //if( DEBUG ) std::cout << " -- next rechit" << std::endl;

    }//<<>>for (const auto recHit : *recHitsEB_ )   
    nRecHits = nRecHitCnt;

    //-------------------------------------------------------------------------------------
    if( DEBUG ) std::cout << "Processing MET" << std::endl;

	auto t1pfMET = (*pfmets_)[0];
	metSumEt = t1pfMET.sumEt(); 
	metPt = t1pfMET.pt();
	metPx = t1pfMET.px();
	metPy = t1pfMET.py();
	metPhi = t1pfMET.phi();

	metCSumEt = metSumEt;
	metCPx = metPx;
	metCPy = metPy;

	// ???  need corrected values for photon Et and Pt
	int nphidx = fphotons.size();
	for( auto phidx = 0; phidx < nphidx; phidx++ ){

		auto photon = fphotons[phidx];
		if( phoIsOotPho[phidx] && not phoExcluded[phidx] ){ 
			metCPx -= ((photon.pt())*std::cos(photon.phi())); 
			metCPy -= ((photon.pt())*std::sin(photon.phi())); 
			metCSumEt -= photon.et();
		}//<<>>if( phoIsOotPho[phidx] && not phoExcluded[phidx] )
		if( not phoIsOotPho[phidx] && phoExcluded[phidx] ){
            metCPx += ((photon.pt())*std::cos(photon.phi()));
            metCPy += ((photon.pt())*std::sin(photon.phi()));
            metCSumEt += photon.et();
		}//<<>>if( phoIsOotPho[phidx] && not phoExcluded[phidx] )	

	}//<<>>for( auto phidx = 0; phidx < fphotons.size(); phidx++ )

	//-------------------------------------------------------------------------------------
	if( DEBUG ) std::cout << "Processing CaloJets" << std::endl;

    nCaloJets = 0;
    //cljSeedTOFTime.clear();
    cljBcCnt.clear();
    cljCMeanTime.clear();
    cljBc3dEx.clear();
    cljBc3dEy.clear();
    cljBc3dEz.clear();
    cljBc3dEv.clear();
    cljBc3dEslope.clear();
    cljBc3dEchisp.clear();
    cljBc2dEx.clear();
    cljBc2dEy.clear();
    cljBc2dEv.clear();
    cljBc2dEslope.clear();
    cljBc2dEchisp.clear();
	cljCDrMeanTime.clear();
    cljPt.clear();
    cljEnergy.clear();
    cljPhi.clear();
    cljEta.clear();
    cljPx.clear();
    cljPy.clear();
    cljPz.clear();
	//cljRhIds.clear();

	cljEMFrac.clear();

	int iCaloJet(0);
	for (const auto caloJet : fcalojets ){

        cljPt.push_back(caloJet.pt());
        cljEnergy.push_back(caloJet.energy());
        cljPhi.push_back(caloJet.phi());
        cljEta.push_back(caloJet.eta());
        cljPx.push_back(caloJet.px());
        cljPy.push_back(caloJet.py());
        cljPz.push_back(caloJet.pz());
		cljEMFrac.push_back(caloJet.emEnergyFraction());
        iCaloJet++;

        rhGroup cljRhGroup;
		int nCljBcs(0);
		bool first = true;
		auto cjeta = caloJet.eta();
    	auto cjphi = caloJet.phi();	
		for( const auto cluster : fbclusts ){
			auto cleta = cluster.eta();
			auto clphi = cluster.phi();
			auto cjcldr = std::sqrt(reco::deltaR2( cjeta, cjphi, cleta, clphi ));		
			if( cjcldr < 0.4 ){ //std::cout << "CJ - BC Match : " << cjcldr << std::endl;
				nCljBcs++;
				if( first ){ cljRhGroup = getRHGroup( cluster, bcMinEnergy/2 ); first = false; }
				else { auto rhgrp = getRHGroup( cluster, bcMinEnergy/2 ); mrgRhGrp( cljRhGroup, rhgrp ); }
			}//<<>>if( cjcldr < 0.1 )
		}//<<>>for( const auto cluster : fbclusts )
        auto cljRhIdsGroup = getRhGrpIDs( cljRhGroup );
        cljRhIds.push_back(cljRhIdsGroup);
        if( cljRhGroup.size() < 3 ){ 
			cljCMeanTime.push_back( -29.75 ); 
			cljSeedTOFTime.push_back( -29.75 );
			cljBc3dEx.push_back(0.0);
			cljBc3dEy.push_back(0.0);
			cljBc3dEz.push_back(0.0);
			cljBc3dEv.push_back(0.0);
			cljBc3dEslope.push_back(0.0);
			cljBc3dEchisp.push_back(0.0);
			cljBc2dEx.push_back(0.0);
			cljBc2dEy.push_back(0.0);
			cljBc2dEv.push_back(0.0);
			cljBc2dEslope.push_back(0.0);
			cljBc2dEchisp.push_back(0.0);
			//if( DEBUG ) std::cout << " - caloJet 2d/3d: below min RH cnt : " << nCljBcs << " w/ " << cljRhGroup.size() << std::endl;
			//continue; 
		} else { //<<>>if( cljRhGroup.size() < minRHcnt ) ***** IF cljRhGroup.size() > minRHcnt BELOW THIS POINT IN LOOP *****
        	auto tofTimes = getLeadTofRhTime( cljRhGroup, vtxX, vtxY, vtxZ );
        	auto timeStats = getTimeDistStats( tofTimes, cljRhGroup );
			auto leadJetRh = getLeadRh( cljRhGroup );
			auto seedTOFTime = leadJetRh.time();
        	//auto seedTOFTime = getSeedTofTime( *scptr, vtxX, vtxY, vtxZ );
        	auto cljBCEigen3D = getRhGrpEigen_ieipt( tofTimes, cljRhGroup );
        	auto cljBCEigen2D = getRhGrpEigen_sph( tofTimes, cljRhGroup );

        	cljBc3dEx.push_back(cljBCEigen3D[0]);
        	cljBc3dEy.push_back(cljBCEigen3D[1]);
        	cljBc3dEz.push_back(cljBCEigen3D[2]);
        	cljBc3dEv.push_back(cljBCEigen3D[3]);
        	cljBc3dEslope.push_back(cljBCEigen3D[4]);
        	cljBc3dEchisp.push_back(cljBCEigen3D[5]);

        	cljBc2dEx.push_back(cljBCEigen2D[0]);
        	cljBc2dEy.push_back(cljBCEigen2D[1]);
        	cljBc2dEv.push_back(cljBCEigen2D[2]);
        	cljBc2dEslope.push_back(cljBCEigen2D[3]);
        	cljBc2dEchisp.push_back(cljBCEigen2D[4]);
	
        	cljSeedTOFTime.push_back(seedTOFTime);
        	cljCMeanTime.push_back(timeStats[6]);
			cljBcCnt.push_back(nCljBcs);
//        	if( DEBUG ) std::cout << " - caloJet 2d/3d: " << cljBCEigen2D[0] << " " << cljBCEigen2D[1] << " " << cljBCEigen2D[2] << " " 
//						<< cljBCEigen2D[3] << " " << cljBCEigen2D[4] << " / " << cljBCEigen3D[0] << " " << cljBCEigen3D[1] << " " << cljBCEigen3D[2] << " "
//                      << cljBCEigen3D[3] << " " << cljBCEigen3D[4] << " " << cljBCEigen3D[5] << " time: " << timeStats[6] << " #" 
//						<< nCljBcs << " w/ " << cljRhGroup.size() << std::endl;
		}//<<>>if( cljRhGroup.size() < minRHcnt )

        const auto cljDrGroup = getRHGroup( cjeta, cjphi, 0.4, 2.0 );
		if( cljRhGroup.size() < 13 ){ cljCDrMeanTime.push_back( -29.75 ); }
		else { // <<>> if( cljRhGroup.size() < 3 )
			auto tofTimes = getLeadTofRhTime( cljRhGroup, vtxX, vtxY, vtxZ );
			auto timeStats = getTimeDistStats( tofTimes, cljRhGroup );
			cljCDrMeanTime.push_back(timeStats[6]);
		}//<<>>if( cljRhGroup.size() < 3 )


	}//<<>>for (const auto caloJet : fcalojets )
    nCaloJets = iCaloJet;

	//-------------------------------------------------------------------------------------------
    if( DEBUG ) std::cout << "Processing Photons" << std::endl;

    int iGedPhos(0);
	nPhotons = 0;
	phoSeedTOFTime.clear();
	phoCMeanTime.clear();
    phoSc3dEx.clear();
	phoSc3dEy.clear();
    phoSc3dEz.clear();
    phoSc3dEv.clear();
    phoSc3dEslope.clear();
    phoSc3dEchisp.clear();
    phoSc2dEx.clear();
    phoSc2dEy.clear();
    phoSc2dEv.clear();
    phoSc2dEslope.clear();
    phoSc2dEchisp.clear();
	phoPt.clear(); 
	phoEnergy.clear();
	phoPhi.clear();
	phoEta.clear();
	phoPx.clear();
	phoPy.clear();
	phoPz.clear();
	phoRhIds.clear();

    phoIsPFPhoton.clear();
    phoIsStdPhoton.clear();
    phoHasConTracks.clear();
    phoIsPixelSeed.clear();
    phoIsPhoton.clear();
    phoIsEB.clear();
    phoIsEE.clear();
    phoIsEBGap.clear();
    phoIsEBEtaGap.clear();
    phoIsEBPhiGap.clear();
    phoIsEEGap.clear();
    phoIsEERingGap.clear();
    phoIsEEDeeGap.clear();
    phoIsEBEEGap.clear();

    phoHadOverEM.clear();
    phoHadD1OverEM.clear();
    phoHadD2OverEM.clear();
    phoHadOverEMVaid.clear();
    phohadTowOverEM.clear();
    phohadTowD10OverEM.clear();
    phohadTowD20OverEM.clear();
    phohadTowOverEMValid.clear();
    phoE1x5.clear();
    phoE2x5.clear();
    phoE3x3.clear();
    phoE5x5.clear();
    phoMaxEnergyXtal.clear();
    phoSigmaEtaEta.clear();
    phoSigmaIEtaIEta.clear();

    phoR1x5.clear();
    phoR2x5.clear();
    phoR9.clear();
    phoFull5x5_e1x5.clear();
    phoFull5x5_e2x5.clear();
    phoFull5x5_e3x3.clear();
    phoFull5x5_e5x5.clear();
    phoFull5x5_maxEnergyXtal.clear();
    phoFull5x5_sigmaEtaEta.clear();
    phoFull5x5_sigmaIEtaIEta.clear();
    phoFull5x5_r1x5.clear();
    phoFull5x5_r2x5.clear();
    phoFull5x5_r9.clear();
    phoNSatXtals.clear();
    phoIsSeedSat.clear();

    phoMipChi2.clear();
    phoMipTotEnergy.clear();
    phoMipSlope.clear();
    phoMipInter.clear();
    phoMipNHitCone.clear();
    phoMipIsHalo.clear();

    phoEcalRHSumEtConeDR04.clear();
    phoHcalTwrSumEtConeDR04.clear();
    phoHcalDepth1TowerSumEtConeDR04.clear();
    phoCalDepth2TowerSumEtConeDR04.clear();
    phoHcalTowerSumEtBcConeDR04.clear();
    phoHcalDepth1TowerSumEtBcConeDR04.clear();
    phoHcalDepth2TowerSumEtBcConeDR04.clear();
    phoTrkSumPtSolidConeDR04.clear();
    phoTrkSumPtHollowConeDR04.clear();
    phoNTrkSolidConeDR04.clear();
    phoNTrkHollowConeDR04.clear();

	genPhoIdx.clear();
    genPhoDr.clear();

    for( const auto photon : fphotons ){

        phoPt.push_back(photon.pt());
        phoEnergy.push_back(photon.energy());
        phoPhi.push_back(photon.phi());
        phoEta.push_back(photon.eta());
        phoPx.push_back(photon.px());
        phoPy.push_back(photon.py());
        phoPz.push_back(photon.pz());

        phoIsPFPhoton.push_back(photon.isPFlowPhoton()) ;
        phoIsStdPhoton.push_back(photon.isStandardPhoton());
        phoHasConTracks.push_back(photon.hasConversionTracks());
        phoIsPixelSeed.push_back(photon.hasPixelSeed());
        phoIsPhoton.push_back(photon.isPhoton());
        phoIsEB.push_back(photon.isEB());
        phoIsEE.push_back(photon.isEE());
        phoIsEBGap.push_back(photon.isEBGap());
        phoIsEBEtaGap.push_back(photon.isEBEtaGap());
        phoIsEBPhiGap.push_back(photon.isEBPhiGap());
        phoIsEEGap.push_back(photon.isEEGap());
        phoIsEERingGap.push_back(photon.isEERingGap());
        phoIsEEDeeGap.push_back(photon.isEEDeeGap());
        phoIsEBEEGap.push_back(photon.isEBEEGap());

        phoHadOverEM.push_back(photon.hadronicOverEm());
        phoHadD1OverEM.push_back(photon.hadronicDepth1OverEm() );
        phoHadD2OverEM.push_back(photon.hadronicDepth2OverEm() );
        phoHadOverEMVaid.push_back(photon.hadronicOverEmValid());
        phohadTowOverEM.push_back(photon.hadTowOverEm());
        phohadTowD10OverEM.push_back(photon.hadTowDepth1OverEm());
        phohadTowD20OverEM.push_back(photon.hadTowDepth2OverEm() );
        phohadTowOverEMValid.push_back(photon.hadTowOverEmValid());
        phoE1x5.push_back(photon.e1x5());
        phoE2x5.push_back(photon.e2x5());
        phoE3x3.push_back(photon.e3x3());
        phoE5x5.push_back(photon.e5x5());
        phoMaxEnergyXtal.push_back(photon.maxEnergyXtal());
        phoSigmaEtaEta.push_back(photon.sigmaEtaEta() );
        phoSigmaIEtaIEta.push_back(photon.sigmaIetaIeta() );

        phoR1x5.push_back(photon.r1x5());
        phoR2x5.push_back(photon.r2x5());
        phoR9.push_back(photon.r9());
        phoFull5x5_e1x5.push_back(photon.full5x5_e1x5());
        phoFull5x5_e2x5.push_back(photon.full5x5_e2x5());
        phoFull5x5_e3x3.push_back(photon.full5x5_e3x3());
        phoFull5x5_e5x5.push_back(photon.full5x5_e5x5());
        phoFull5x5_maxEnergyXtal.push_back(photon.full5x5_maxEnergyXtal());
        phoFull5x5_sigmaEtaEta.push_back(photon.full5x5_sigmaEtaEta());
        phoFull5x5_sigmaIEtaIEta.push_back(photon.full5x5_sigmaIetaIeta());
        //phoFull5x5_r1x5.push_back(photon.full5x5_r1x5());
        //phoFull5x5_r2x5.push_back(photon.full5x5_r2x5());
        phoFull5x5_r9.push_back(photon.full5x5_r9());

        phoNSatXtals.push_back(photon.nSaturatedXtals());
        phoIsSeedSat.push_back(photon.isSeedSaturated());

        phoMipChi2.push_back(photon.mipChi2());
        phoMipTotEnergy.push_back(photon.mipTotEnergy());
        phoMipSlope.push_back(photon.mipSlope());
        phoMipInter.push_back(photon.mipIntercept());

        phoMipNHitCone.push_back(photon.mipNhitCone());
        phoMipIsHalo.push_back(photon.mipIsHalo());

        phoEcalRHSumEtConeDR04.push_back(photon.ecalRecHitSumEtConeDR04());
        phoHcalTwrSumEtConeDR04.push_back(photon.hcalTowerSumEtConeDR04());
        phoHcalDepth1TowerSumEtConeDR04.push_back(photon.hcalDepth1TowerSumEtConeDR04());
        phoCalDepth2TowerSumEtConeDR04.push_back(photon.hcalDepth2TowerSumEtConeDR04());
        phoHcalTowerSumEtBcConeDR04.push_back(photon.hcalTowerSumEtBcConeDR04() );
        phoHcalDepth1TowerSumEtBcConeDR04.push_back(photon.hcalDepth1TowerSumEtBcConeDR04());
        phoHcalDepth2TowerSumEtBcConeDR04.push_back(photon.hcalDepth2TowerSumEtBcConeDR04());
        phoTrkSumPtSolidConeDR04.push_back(photon.trkSumPtSolidConeDR04());
		phoTrkSumPtHollowConeDR04.push_back(photon.trkSumPtHollowConeDR04());
        phoNTrkSolidConeDR04.push_back(photon.nTrkSolidConeDR04());
        phoNTrkHollowConeDR04.push_back(photon.nTrkHollowConeDR04());

		if( DEBUG ) std::cout << " --- Proccesssing : " << photon << std::endl;
        const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
        const auto scptr = phosc.get();
        scGroup phoSCGroup{*scptr};
        //auto phoRhGroup = getRHGroup( phoSCGroup, 2.0, //hist1d[137], hist1d[138], hist1d[139] );
        auto phoRhGroup = getRHGroup( phoSCGroup, 0.0 );//bcMinEnergy );
        auto phoRhIdsGroup = getRhGrpIDs( phoRhGroup );
		phoRhIds.push_back(phoRhIdsGroup);
        if( DEBUG ) std::cout << " -- gedPhotons : " << scptr << " #: " << phoRhGroup.size() << std::endl;
        if( phoRhGroup.size() < minObjRHcnt ){ 
            phoCMeanTime.push_back( -29.75 ); 
			phoSeedTOFTime.push_back( -29.75 );
            phoSc3dEx.push_back(0.0);
            phoSc3dEy.push_back(0.0);
            phoSc3dEz.push_back(0.0);
            phoSc3dEv.push_back(0.0);
            phoSc3dEslope.push_back(0.0);
            phoSc3dEchisp.push_back(0.0);
            phoSc2dEx.push_back(0.0);
            phoSc2dEy.push_back(0.0);
            phoSc2dEv.push_back(0.0);
            phoSc2dEslope.push_back(0.0);
            phoSc2dEchisp.push_back(0.0);
			if( DEBUG ) std::cout << " - Photon 2d/3d: below min RH cnt: " << std::endl;
            continue;
        }//<<>>if( phoRhGroup.size() < minRHcnt ) ***** IF phoRhGroup.size() > minRHcnt BELOW THIS POINT IN LOOP *****
        auto tofTimes = getLeadTofRhTime( phoRhGroup, vtxX, vtxY, vtxZ );
        auto timeStats = getTimeDistStats( tofTimes, phoRhGroup );
		auto seedTOFTime = getSeedTofTime( *scptr, vtxX, vtxY, vtxZ );
        //auto phoLeadTOFTime =  getLeadTofTime( phoRhGroup, vtxX, vtxY, vtxZ );
		auto phoSCEigen3D = getRhGrpEigen_ieipt( tofTimes, phoRhGroup );
		auto phoSCEigen2D = getRhGrpEigen_sph( tofTimes, phoRhGroup );
		
        phoSc3dEx.push_back(phoSCEigen3D[0]);
    	phoSc3dEy.push_back(phoSCEigen3D[1]);
    	phoSc3dEz.push_back(phoSCEigen3D[2]);
        phoSc3dEv.push_back(phoSCEigen3D[3]);
        phoSc3dEslope.push_back(phoSCEigen3D[4]);
        phoSc3dEchisp.push_back(phoSCEigen3D[5]);

        phoSc2dEx.push_back(phoSCEigen2D[0]);
        phoSc2dEy.push_back(phoSCEigen2D[1]);
        phoSc2dEv.push_back(phoSCEigen2D[2]);
        phoSc2dEslope.push_back(phoSCEigen2D[3]);
        phoSc2dEchisp.push_back(phoSCEigen2D[4]);

		phoSeedTOFTime.push_back(seedTOFTime);
        phoCMeanTime.push_back(timeStats[6]);
		if( DEBUG ) std::cout << " - Photon 2d/3d: " << phoSCEigen2D[0] << " " << phoSCEigen2D[1] << " " << phoSCEigen2D[2] << " " << phoSCEigen2D[3] << " "
							  << phoSCEigen2D[4] << " / " << phoSCEigen3D[0] << " " << phoSCEigen3D[1] << " " << phoSCEigen3D[2] << " " 
                              << phoSCEigen3D[3] << " " << phoSCEigen3D[4] << " " << phoSCEigen3D[5] << " time: " << timeStats[6] << std::endl;


        // GenParticle Info for photon  -------------------------------------------------------------------
        if( DEBUG ) std::cout << "Getting phoGenParton Information" << std::endl;

		// set defaults for no match here
        auto mgenDr = -1.0;
	
        if( hasGenInfo ){
            if( DEBUG ) std::cout << " -- Pulling pho gen info " << std::endl;

            bool matchfound(false);
            float goodDr(0.4);
            int matchedIdx(-1);
            int index(0);
			float meta(9);
            float mphi(9);
            float geta(9);
            float gphi(9);

            for(const auto& genPart : fgenparts ){

                const auto rhX = scptr->x();
                const auto rhY = scptr->y();
                const auto rhZ = scptr->z();
                const auto gnX = genPart.vx();
                const auto gnY = genPart.vy();
                const auto gnZ = genPart.vz();
				auto cphoEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
				auto cphoPhi = std::atan2(rhY-gnY,rhX-gnX);
                auto dr = std::sqrt(reco::deltaR2(genPart.eta(), genPart.phi(), cphoEta, cphoPhi ));
                auto isPhoton = std::abs(genPart.pdgId()) == 22;

                if( ( dr < goodDr ) && isPhoton ){
                    matchfound = true;
                    goodDr = dr;
                    matchedIdx = index;
					meta = cphoEta;
					mphi = cphoPhi;
					geta = genPart.eta();
					gphi = genPart.phi();
                }//<<>>if( jtgjdr <= goodDr )
                index++;

            }//<<>>for(const auto& genJet : *genJets_ ) 

            if( matchedIdx >= 0 ){ mgenDr = goodDr; }//<<>>if( matchedIdx >= 0 ) 
            if( not matchfound ){ if( DEBUG ) std::cout << " - genpho GenTime : genPho == 0 " << std::endl; } //genPhoTime = -50.0; }
			if( DEBUG) std::cout << "------------------------- " << std::endl;
			if( DEBUG) std::cout << " Matched Eta - ph: " <<  photon.eta() << " gen: " << geta << " match: " << meta; 
			if( DEBUG) std::cout << " oot: " << phoIsOotPho[iGedPhos]  << std::endl;
            if( DEBUG) std::cout << " Matched Eta - ph: " <<  photon.phi() << " gen: " << gphi << " match: " << mphi << std::endl;
    		genPhoIdx.push_back(matchedIdx);
			genPhoDr.push_back(mgenDr);

            if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;

        }//<<>>if( hasGenInfo )

        iGedPhos++;

    }//<<>>for( const auto photon : *gedPhotons_ )
	nPhotons = iGedPhos;

	//------------------------------------------------------------------------------------------

    nElectrons = 0;
    eleSeedTOFTime.clear();
    eleCMeanTime.clear();
    eleSc3dEx.clear();
    eleSc3dEy.clear();
    eleSc3dEz.clear();
    eleSc3dEv.clear();
    eleSc3dEslope.clear();
    eleSc3dEchisp.clear();
    eleSc2dEx.clear();
    eleSc2dEy.clear();
    eleSc2dEv.clear();
    eleSc2dEslope.clear();
    eleSc2dEchisp.clear();
    elePt.clear();
    eleEnergy.clear();
    elePhi.clear();
    eleEta.clear();
    elePx.clear();
    elePy.clear();
    elePz.clear();
	//eleRhIdsGroup.clear();

    int iElectros(0);
	scGroup jetEleSCGroup;
    if( DEBUG ) std::cout << "Processing Electrons" << std::endl;
	for( const auto electron : felectrons ){

        elePt.push_back(electron.pt());
        eleEnergy.push_back(electron.energy());
        elePhi.push_back(electron.phi());
        eleEta.push_back(electron.eta());
        elePx.push_back(electron.px());
        elePy.push_back(electron.py());
        elePz.push_back(electron.pz());

        iElectros++;

		if( DEBUG ) std::cout << " --- Proccesssing : " << electron << std::endl;
    	const auto &elesc = electron.superCluster().isNonnull() ? electron.superCluster() : electron.parentSuperCluster();
        const auto scptr = elesc.get();
        scGroup eleSCGroup{*scptr};
    	//auto eleRhGroup = getRHGroup( eleSCGroup, 2.0, //hist1d[143], hist1d[144], hist1d[145] );
        auto eleRhGroup = getRHGroup( eleSCGroup, bcMinEnergy );
        //auto eleRhIdsGroup = getRhGrpIDs( eleRhGroup );
        //eleRhIds.push_back(eleRhIdsGroup);
        if( DEBUG ) std::cout << " -- Electrons : " << scptr << " #: " << eleRhGroup.size() << std::endl;
        if( eleRhGroup.size() < minObjRHcnt ){
            eleCMeanTime.push_back( -29.75 );
            eleSeedTOFTime.push_back( -29.75 );
            eleSc3dEx.push_back(0.0);
            eleSc3dEy.push_back(0.0);
            eleSc3dEz.push_back(0.0);
            eleSc3dEv.push_back(0.0);
            eleSc3dEslope.push_back(0.0);
            eleSc3dEchisp.push_back(0.0);
            eleSc2dEx.push_back(0.0);
            eleSc2dEy.push_back(0.0);
            eleSc2dEv.push_back(0.0);
            eleSc2dEslope.push_back(0.0);
            eleSc2dEchisp.push_back(0.0);
			if( DEBUG ) std::cout << " - Photon 2d/3d: below min RH cnt" << std::endl;
            continue;
        }//<<>>if( eleRhGroup.size() < minRHcnt ) ***** IF eleRhGroup.size() > minRHcnt BELOW THIS POINT IN LOOP *****
        auto tofTimes = getLeadTofRhTime( eleRhGroup, vtxX, vtxY, vtxZ );
        auto timeStats = getTimeDistStats( tofTimes, eleRhGroup );
        auto seedTOFTime = getSeedTofTime( *scptr, vtxX, vtxY, vtxZ );
        auto eleSCEigen3D = getRhGrpEigen_ieipt( tofTimes, eleRhGroup );
        auto eleSCEigen2D = getRhGrpEigen_sph( tofTimes, eleRhGroup );

        eleSc3dEx.push_back(eleSCEigen3D[0]);
        eleSc3dEy.push_back(eleSCEigen3D[1]);
        eleSc3dEz.push_back(eleSCEigen3D[2]);
        eleSc3dEv.push_back(eleSCEigen3D[3]);
        eleSc3dEslope.push_back(eleSCEigen3D[4]);
        eleSc3dEchisp.push_back(eleSCEigen3D[5]);

        eleSc2dEx.push_back(eleSCEigen2D[0]);
        eleSc2dEy.push_back(eleSCEigen2D[1]);
        eleSc2dEv.push_back(eleSCEigen2D[2]);
        eleSc2dEslope.push_back(eleSCEigen2D[3]);
        eleSc2dEchisp.push_back(eleSCEigen2D[4]);

        eleSeedTOFTime.push_back(seedTOFTime);
        eleCMeanTime.push_back(timeStats[6]);
        if( DEBUG ) std::cout << " - Electron 2d/3d: " << eleSCEigen2D[0] << " " << eleSCEigen2D[1] << " " << eleSCEigen2D[2] << " " << eleSCEigen2D[3] << " "
                              << eleSCEigen2D[4] << " / " << eleSCEigen3D[0] << " " << eleSCEigen3D[1] << " " << eleSCEigen3D[2] << " " 
                              << eleSCEigen3D[3] << " " << eleSCEigen3D[4] << " " << eleSCEigen3D[5] << " time: " << timeStats[6] << std::endl;

	}//<<>>for( const auto electron : *electrons_ )
	nElectrons = iElectros;

	//------------------------------------------------------------------------------------------
    nUnJets = (*jets_).size();
    nJets = fjets.size();
    if( DEBUG ) std::cout << "Init for Jet Loop with " << nJets << " jets"<< std::endl;

	nGoodDrJets = 0;
   	nGoodScJets = 0;
   	nGoodBcJets = 0;

	vector<uInt> noIdGrp(1,0);

	//set the number of leading jets to skim ( = nJets for all )
	//auto stJets = nJets; 

    jetHt = 0.0;

	jetE.clear();
	jetPt.clear(); 
	jetPhi.clear(); 
	jetEta.clear(); 
	jetEtaetaMmt.clear();
	jetPhiphiMnt.clear();
	jetEtaphiMnt.clear();
	jetMaxD.clear();
	jetConPtDis.clear();
	jetConEtaPhiSprd.clear();
	jetArea.clear();
	jetNCarry.clear();
	jetNConst.clear();

	jetID.clear();
	jetNHF.clear();
	jetNEMF.clear();	
	jetCHF.clear();
	jetCEMF.clear();
	jetMUF.clear();
	jetNHM.clear();
	jetCHM.clear();

	jetMuTime.clear();
	jetTimeError.clear();
	jetTimeRMS.clear();
	jetMedTime.clear();
	jetSCMuTime.clear();
	jetSCMedTime.clear();
	jetCMuTime.clear();
	jetCMedTime.clear();
	jetCSCMuTime.clear();
	jetCSCMedTime.clear();
	jetCBCMuTime.clear();
	jetCBCMedTime.clear();
    jetPhMuTime.clear();
    jetOOTPhMuTime.clear();
    jetEleMuTime.clear();

    //jetDrRhIds.clear();
	//jetScRhIds.clear();

	jetPHM.clear();
	jetELM.clear();
	//jetC.clear();
	jetPHE.clear();
	jetPHEF.clear();
	jetELE.clear();
	jetELEF.clear();
	jetMUE.clear();
	//jetCharge.clear();

	njetRecHits.clear();
	jetRecHitOfJet.clear();
	jetRecHitId.clear();

	jetSumEPFrac.clear();
    jetEPEnergy.clear();
    jetEMEnergy.clear();
    jetEMEnrFrac.clear();
    jetEPEnrFrac.clear();
    
    jetDrLeadEta.clear();
    jetDrLeadPhi.clear();
    jetDrLeadEnr.clear();
    
    sJetDrRHEnergy.clear();
    jetDrEMF.clear();
    jetDrRhCnt.clear();
    
    jetGenImpactAngle.clear();
    jetGenTime.clear();
    jetGenPt.clear();
    jetGenEta.clear();
    jetGenEnergy.clear();
    jetGenEMFrac.clear();
    jetGenDrMatch.clear();
    jetGenTimeVar.clear();
    jetGenTimeLLP.clear();
    jetGenLLPPurity.clear();
    jetGenNextBX.clear();
    jetGenNKids.clear();
    jetGenTOF.clear();
    
    nJetScMatch.clear();
    sJetScEnergy.clear();
    sJetScPhEnergy.clear();

	jetScRhCnt.clear();    
    sJetScRhEnergy.clear();
    jetScEMF.clear();
    
    jetImpactAngle.clear();
    jetBcTimesCnt.clear();
    jetBcSumRHEnr.clear();
    jetBcEMFr.clear();
    jetBcRhCnt.clear();
	jetBcGrpCnt.clear();

    jetSc3dEx.clear();
    jetSc3dEy.clear();
    jetSc3dEz.clear();
    jetSc3dEv.clear();
    jetSc3dEslope.clear();
    jetSc3dEchisp.clear();

    jetSc2dEx.clear();
    jetSc2dEy.clear();
    jetSc2dEv.clear();
    jetSc2dEslope.clear();
    jetSc2dEchisp.clear();
    jetSc2dEslope2.clear();
    jetSc2dEchisp2.clear();
    jetSc2dErangle.clear();
    jetSc2dEnxsum.clear();

//////////  reffrence information
// **** const reco::TrackRefVector& associatedTracks() const;
// **** const reco::PixelClusterTagInfo* tagInfoPixelCluster(const std::string& label = "") const; // from b tagging info methods
// **** reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }

	////hist1d[18]->Fill(nJets);
	bool goodJetEvent(false);
	if( DEBUG ) std::cout << "-----------------------------------------------------------" << std::endl;
	if( DEBUG ) std::cout << "Starting Jet Loop for " << nJets << " jets " << std::endl; 
	for ( uInt ijet(0); ijet < nJets; ijet++ ){ 
		// places jet info in output tree
	
		std::vector<reco::CaloCluster> unusedfbclusts(fbclusts);	
		//std::cout << "Processing with jet " << ijet << std::endl;
	   	const auto & jet = fjets[ijet];
		const auto jetid = fjetsID[ijet];
		if( DEBUG ) std::cout << " --- Proccesssing : " << jet  << std::endl;

	   	// jetID in jet.h ?

	   	jetHt += jet.pt();

	   	jetE.push_back(jet.energy());
	   	jetPt.push_back(jet.pt());
	   	jetPhi.push_back(jet.phi());
	   	jetEta.push_back(jet.eta());

		jetEtaetaMmt.push_back(jet.etaetaMoment());
		jetPhiphiMnt.push_back(jet.phiphiMoment());
		jetEtaphiMnt.push_back(jet.etaphiMoment());
		jetMaxD.push_back(jet.maxDistance());
		jetConPtDis.push_back(jet.constituentPtDistribution());
		jetConEtaPhiSprd.push_back(jet.constituentEtaPhiSpread());
		jetArea.push_back(jet.jetArea());
		jetNCarry.push_back(jet.nCarrying(0.1));
		jetNConst.push_back(jet.nConstituents());

	   	//jetID.push_back(jet.userInt("jetID"));
		jetID.push_back(jetid);
	   	jetNHF.push_back(jet.neutralHadronEnergyFraction());
	   	jetNEMF.push_back(jet.neutralEmEnergyFraction());
	   	jetCHF.push_back(jet.chargedHadronEnergyFraction());
	   	jetCEMF.push_back(jet.chargedEmEnergyFraction());
	   	jetMUF.push_back(jet.muonEnergyFraction());
	   	jetNHM.push_back(jet.neutralMultiplicity());
		jetCHM.push_back(jet.chargedMultiplicity());
		//jetCharge.push_back(jet.jetCharge());

		jetPHE.push_back(jet.photonEnergy());
		jetPHEF.push_back(jet.photonEnergyFraction()); 
		jetELE.push_back(jet.electronEnergy());
		jetELEF.push_back(jet.electronEnergyFraction());
		jetMUE.push_back(jet.muonEnergy());
		jetPHM.push_back(jet.photonMultiplicity());
		jetELM.push_back(jet.electronMultiplicity());

   //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
   
	   	//if( DEBUG ) std::cout << "Fill jet pt/phi/eta Histograms" << std::endl;
 		
      	const auto jetepafrac 	= jet.photonEnergyFraction() + jet.electronEnergyFraction();
      	const auto jetepe 		= jet.photonEnergy() + jet.electronEnergy();
		const auto jeteme 		= jet.chargedEmEnergy() + jet.neutralEmEnergy();
      	const auto jetemfrac 	= jeteme/jet.energy();
      	const auto jetepfrac 	= jetepe/jet.energy();

		jetSumEPFrac.push_back(jetepafrac);
		jetEPEnergy.push_back(jetepe);
		jetEMEnergy.push_back(jeteme);
		jetEMEnrFrac.push_back(jetemfrac);
		jetEPEnrFrac.push_back(jetepfrac);

   //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

        // dR matched RH group  -----------------------------------------------
		// --------------------------------------------------------------------

	   	if( DEBUG ) std::cout << "Getting jet dR rechit group" << std::endl; 
		auto jetDrRhGroup = getRHGroup( jet.eta(), jet.phi(), deltaRminJet, minRHenr ); 
        auto jetDrRhIdsGroup = getRhGrpIDs( jetDrRhGroup );
        jetDrRhIds.push_back( jetDrRhIdsGroup );
		auto rhCount = jetDrRhGroup.size();

	   	//std::cout << "rhCount is " << rhCount << std::endl;
	   	auto sumdrrhe = getRhGrpEnr( jetDrRhGroup );
		auto dremf = sumdrrhe/jet.energy();
		sJetDrRHEnergy.push_back(sumdrrhe);
		jetDrEMF.push_back(dremf);
        jetDrRhCnt.push_back(rhCount);

		//float jetDrTime(-99.9);
		if( rhCount >= minRHcnt && dremf > minEmf ){

			if( DEBUG ) std::cout << " - get jetDRtofTimes " << std::endl;
	   		auto tofTimes = getLeadTofRhTime( jetDrRhGroup, vtxX, vtxY, vtxZ );
	   		auto leadJetRh = getLeadRh( jetDrRhGroup );


	   		if( DEBUG ) std::cout << "Starting RecHit Loop" << std::endl;
			for ( uInt irhg = 0; irhg < rhCount; irhg++){
			
				//std::cout << " -- irhg: " << irhg << " rhCount: " << rhCount << std::endl;
				jetRecHitOfJet.push_back(ijet);
				auto detid = (jetDrRhGroup[irhg]).detid();
				//std::cout << " -- (jetDrRhGroup[irhg]).detid(): " << detid.rawId() << std::endl;
		      	jetRecHitId.push_back(detid.rawId());	
				//std::cout << " -- tofTimes[irhg]: " << rhtime << std::endl;
		      	//std::cout << " -- jetDrRhGroup[irhg]).energy(): " << rhe << std::endl;
		   	}//<<>>for ( uInt irhg = 0; irhg < rhCount; irhg++)
			
		   	const auto leadJetRhId = leadJetRh.detid();
		   	const auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
		   	auto sc_eta = leadJetRhIdPos.eta();
		   	auto sc_phi = leadJetRhIdPos.phi();
		   	auto sc_enr = leadJetRh.energy();
		   	//std::cout << "Lead Jet dR RH Group E: " << sc_enr << " eta: " << sc_eta << " phi: " << sc_phi << std::endl;

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<if( rhCount >= minRHcnt ){

			//  make jettime varible
	      	//for( auto t : tofTimes ) //hist1d[5]->Fill(t);
			auto jetTimeStats = getTimeDistStats( tofTimes, jetDrRhGroup ); 
			auto jmutime = jetTimeStats[0];
		   	auto jterr = jetTimeStats[1];
		   	auto jtrms = jetTimeStats[4];
			auto jmedtime = jetTimeStats[2];
	      	//auto mederr = jetTimeStats[3];
			auto jcmutime = jetTimeStats[6];
	      	auto jcmedtime = jetTimeStats[10];

		   	jetMuTime.push_back(jmutime);
		   	jetTimeError.push_back(jterr);
		   	jetTimeRMS.push_back(jtrms);
		   	jetMedTime.push_back(jmedtime);
	      	jetCMuTime.push_back(jcmutime);	
	      	jetCMedTime.push_back(jcmedtime);

			jetDrLeadEta.push_back(sc_eta);
            jetDrLeadPhi.push_back(sc_phi);
            jetDrLeadEnr.push_back(sc_enr);

			if( jcmutime > -28.9 ) nGoodDrJets++;			


   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<if( rhCount >= minRHcnt ){


		} else { //<<>>if( rhCount > minRHcnt && dremf > minEmf )

			njetRecHits.push_back(0);
	      	jetMuTime.push_back(-29.25);
	      	jetTimeError.push_back(99.9);
	      	jetTimeRMS.push_back(99.9);
	      	jetMedTime.push_back(-29.25);
	      	jetCMuTime.push_back(-29.25);
	      	jetCMedTime.push_back(-29.25);
			
		}//<<>>if( rhCount > minRHcnt && dremf > minEmf ) : else
		
	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

		// GenJet Info for MC  -------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

        if( DEBUG ) std::cout << "Getting jetGenParton Information" << std::endl;

        //const reco::GenParticle * jetGenParton(0);
        float genImpactAngle(7.0);
        float genTime(-99.9);
        float genPt(-1.0);
        float genEta(10.0);
        float genEnergy(-1.0);
        float genEMFrac(-1.0);
        float genDrMatch(-1.0);
		float genTimeVar(-1.0);
        float genTimeLLP(0.0);
        float genLLPPurity(-1.0);
        float genNextBX(-1.0);
        float genNKids(-1.0);
        float genTOF(-99.9);

        if( hasGenInfo ){

            //if( DEBUG ) std::cout << " -- Pulling jet gen info " << std::endl;

			bool matchfound(false);
			float goodDr(0.1);
			int matchedIdx(-1);
			int index(0);
			for(const auto& genJet : *genJets_ ){

            	auto gjeta = genJet.eta();
            	auto gjphi = genJet.phi();
                auto jtgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jet.eta(), jet.phi() ));
				if( jtgjdr <= goodDr ){
					matchfound = true;
					goodDr = jtgjdr;	
					matchedIdx = index;
				}//<<>>if( jtgjdr <= goodDr )
				index++;

            }//<<>>for(const auto& genJet : *genJets_ ) 

            if( matchedIdx >= 0 && rhCount > 0 ){

				auto genJet = (*genJets_)[matchedIdx];
				if( DEBUG ) std::cout << " --- Jet-GenJet dR match : " << goodDr << std::endl;

            	//if( DEBUG ) std::cout << " ---- genParton : " << genParton << " genJet : " << genJet << std::endl;
            	// size_t numberOfDaughters() const override;
            	// size_t numberOfMothers() const override;
            	// const Candidate * daughter(size_type) const override;
            	// const Candidate * mother(size_type = 0) const override;
            	// Candidates :
            	//  virtual double vx() const  = 0;
            	//  virtual double vy() const  = 0;
            	//  virtual double vz() const  = 0;
            	//  virtual int pdgId() const  = 0;
            	//  + mass/energy/momentum info ....

            	//auto nMother = genJet.numberOfMothers();
            	//auto nDaughter = genJet.numberOfDaughters();
            	auto nSources = genJet.numberOfSourceCandidatePtrs();
            	if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;
            	//std::cout << " - genJet mothers : " << nMother << " daughters : " << nDaughter << " sources : " << nSources << std::endl;
            	if( DEBUG ) std::cout << " - genJet srcs : " << nSources << " PV (" << vtxX << "," << vtxY << "," << vtxZ << ")" << std::endl;
            	auto kids = genJet.daughterPtrVector();
            	//std::cout << bigKidChase( kids, vtxX ) << std::endl;
            	//kidChase( kids, vtxX, vtxY, vtxZ );
				auto leadJetRh = getLeadRh( jetDrRhGroup );
            	auto leadJetRhId = leadJetRh.detid();
            	auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
				auto cx = leadJetRhIdPos.x();
				auto cy = leadJetRhIdPos.y();
				auto cz = leadJetRhIdPos.z();
				auto tofcor = hypo( cx, cy, cz )/SOL;
				//if( DEBUG ) kidChase( kids, vtxX, vtxY, vtxZ );
            	auto genKidInfo = kidTOFChain( kids, cx, cy, cz );
				if( DEBUG ) std::cout << " - genJet GenTime noTOF : " << genKidInfo[0] << " rhPos: " << cx << "," << cy << "," << cz << std::endl;
				genEta = genJet.eta();
				if( genKidInfo[0] > 25.0 ) genTime = -28.0;
				else if( genKidInfo[0] > -25.0 ) genTime = genKidInfo[0]-tofcor;
				else genTime = -27.0;
				genImpactAngle = genKidInfo[1];
				if( DEBUG ) std::cout << " - genJet GenTime : " << genTime << " Angle: " << genImpactAngle << std::endl;
				genPt = genJet.pt();
            	genEnergy = genJet.energy();
				genEMFrac = (genJet.chargedEmEnergy() + genJet.neutralEmEnergy())/genEnergy;
				genDrMatch = goodDr; //std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), genJet.eta(), genJet.phi()));
				genTimeVar = genKidInfo[2];
            	genNextBX = genKidInfo[3];
            	genTimeLLP = genKidInfo[4];
            	genLLPPurity = genKidInfo[5];
				genNKids = genKidInfo[6];
				genTOF = tofcor;
				if( DEBUG ) std::cout << " -- Energy : " << genEnergy << " Pt : " << genPt << " EMfrac : " << genEMFrac << std::endl;
					
            }//<<>>if( matchedIdx >= 0 && rhCount > 0 ) 

			if( not matchfound ){ if( DEBUG ) std::cout << " - genJet GenTime : genJet == 0 " << std::endl; genTime = -50.0; }

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
		//<<<<if( hasGenInfo )
    
			// load event level vectors for this jet with gen info

            jetGenImpactAngle.push_back(genImpactAngle);
            jetGenTime.push_back(genTime);
            jetGenPt.push_back(genPt);
            jetGenEta.push_back(genEta);
            jetGenEnergy.push_back(genEnergy);
            jetGenEMFrac.push_back(genEMFrac);
            jetGenDrMatch.push_back(genDrMatch);
            jetGenTimeVar.push_back(genTimeVar);
            jetGenTimeLLP.push_back(genTimeLLP);
            jetGenLLPPurity.push_back(genLLPPurity);
            jetGenNextBX.push_back(genNextBX);
            jetGenNKids.push_back(genNKids);
            jetGenTOF.push_back(genTOF);

			if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;

        }//<<>>if( hasGenInfo )


		// Super Cluster group	-----------------------------------------------
		// --------------------------------------------------------------------		

		if( DEBUG ) std::cout << " - Get jet SC Group " << std::endl;
		int nMatched(0);
		int sum_nrh(0);
	   	float sum_sce(0.0);
	   	float sum_phe(0.0);
        //auto jetSCTime(-29.25);
		scGroup jetSCGroup;
        scGroup jetPhSCGroup;
        scGroup jetOOTPhSCGroup;
        scGroup jetEleSCGroup;
		bcGroup jetBCGroup;
		vector<float> phEnergy;
        vector<float> phDr;
        vector<float> ootPhEnergy;
        vector<float> ootPhDr;
        vector<float> eleEnergy;
        vector<float> eleDr;

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

	if( true ) { //---------------------------------- gedPhotons lock ------------------------------------------------
        int iph(0);
        bool pmatched(false);
		if( DEBUG ) std::cout << "Proccesssing Photon :" << std::endl;
	   	for( const auto photon : *gedPhotons_ ){

			if( DEBUG ) std::cout << " --- Proccesssing : " << photon  << std::endl;

			//----------    refrence to cluster from supercluster
			//if( DEBUG ) std::cout << "----------------------------------" << std::endl;

        	auto pheta = photon.eta();
        	auto phphi = photon.phi();
			for( const auto kid : jet.daughterPtrVector() ){
				auto kidcand = (pfcands_->ptrAt(kid.key())).get();
				//std::cout << " - Matching : " << pfcand.eta() << " = " << kidcand->eta() << std::endl;
				auto kdeta = kidcand->eta();
				auto kdphi = kidcand->phi();
				auto kdpfdr = std::sqrt(reco::deltaR2(pheta, phphi, kdeta, kdphi ));
				if( kdpfdr <= 0.001 ){ pmatched = true; if( DEBUG ) std::cout << " --- Matched : " << kdpfdr << std::endl; }
			}//<<>>for( const auto kid : jet.daughterPtrVector() )


   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<for( const auto photon : *gedPhotons_ ){ 
  
	      	if( pmatched ){
				if( DEBUG ) std:: cout << " ----- Photon Match !!!! " << std::endl;
				iph++;
				nMatched++;
	         	const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
				const auto scptr = phosc.get();
				jetSCGroup.push_back(*scptr);

				const auto clusters = phosc->clusters();
	      		//const auto nBClusts = phosc->clustersSize();
	         	const auto &hitsAndFractions = phosc->hitsAndFractions();
	         	const auto nrh = hitsAndFractions.size();
				//std::cout << " -- SC match found with nBClusts: " << nBClusts << " nClusters: "<< clusters.size() << std::endl;
				//std::cout << " -- SC has nRecHits: " << nrh << std::endl;
				if( nrh < minRHcnt ) continue;
				if( nrh != 0 ){
					for( const auto &clustptr : clusters ){
						if( not clustptr.isAvailable() ) continue;
	         			const auto clust = clustptr.get();
						jetBCGroup.push_back(*clust);
						//std::cout << " --- Adding cluster " << std::endl;
					}//<<>>for( const auto &clustptr : clusters ){
				}//<<>>if( nrh != 0 ){
	         	const auto sce = phosc->energy();
				sum_nrh += nrh;
				sum_sce += sce;
				sum_phe += photon.energy();

                pmatched = false;
	    	}//<<>>if( matched )

		}//<<>>for( const auto photon : *gedPhotons_ ) 
	} // ------------ gedPhotons lock ------------------------------------------------------------------

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

    if( true ) { //---------------------------------- ootPhotons lock ------------------------------------------------
        int iootph(0);
        auto ootmatched = false;
		if( DEBUG ) std::cout << "Proccesssing OOTPhoton :" << std::endl;
  		for( const auto ootphoton : *ootPhotons_ ){

            if( DEBUG ) std::cout << " --- Proccesssing : " << ootphoton  << std::endl;
            auto pheta = ootphoton.eta();
            auto phphi = ootphoton.phi();
            for( const auto kid : jet.daughterPtrVector() ){
                auto kidcand = (pfcands_->ptrAt(kid.key())).get();
                //std::cout << " - Matching : " << pfcand.eta() << " = " << kidcand->eta() << std::endl;
                auto kdeta = kidcand->eta();
                auto kdphi = kidcand->phi();
                auto kdpfdr = std::sqrt(reco::deltaR2(pheta, phphi, kdeta, kdphi ));
                if( kdpfdr <= 0.001 ){ ootmatched = true; if( DEBUG ) std::cout << " --- Matched : " << kdpfdr << std::endl; }
            }//<<>>for( const auto kid : jet.daughterPtrVector() )

            if( ootmatched ){
				if( DEBUG ) std:: cout << " ----- OOT Photon Match !!!! " << std::endl;
                iootph++;
                nMatched++;
                const auto &ootphosc = ootphoton.superCluster().isNonnull() ? ootphoton.superCluster() : ootphoton.parentSuperCluster();
                const auto scptr = ootphosc.get();
                jetSCGroup.push_back(*scptr);
                const auto clusters = ootphosc->clusters();
                const auto &hitsAndFractions = ootphosc->hitsAndFractions();
                const auto nrh = hitsAndFractions.size();
                //std::cout << " -- SC match found with nBClusts: " << nBClusts << " nClusters: "<< clusters.size() << std::endl;
                //std::cout << " -- SC has nRecHits: " << nrh << std::endl;
                if( nrh < minRHcnt ) continue;
                if( nrh != 0 ){
                    for( const auto &clustptr : clusters ){
                        if( not clustptr.isAvailable() ) continue;
                        const auto clust = clustptr.get();
                        jetBCGroup.push_back(*clust);
                        //std::cout << " --- Adding cluster " << std::endl;
                    }//<<>>for( const auto &clustptr : clusters ){
                }//<<>>if( nrh != 0 ){
                const auto sce = ootphosc->energy();
                sum_nrh += nrh;
                sum_sce += sce;
                sum_phe += ootphoton.energy();

                ootmatched = false;
            }//<<>>if( matched )

        }//<<>>for( const auto photon : *ootPhotons_ ) 

	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
    } // ------------ ootPhotons lock ------------------------------------------------------------------		

    if( true ) { //---------------------------------- electrons lock ------------------------------------------------
        int iel(0);
        bool ematched = false;
		if( DEBUG ) std::cout << "Proccesssing Electron :" << std::endl;
	   	for( const auto electron : *electrons_ ){

			if( DEBUG ) std::cout << " --- Proccesssing : " << electron  << std::endl;

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<for( const auto electron : *electrons_ ){

            auto eleeta = electron.eta();
            auto elephi = electron.phi();
            for( const auto kid : jet.daughterPtrVector() ){
                auto kidcand = (pfcands_->ptrAt(kid.key())).get();
                //std::cout << " - Matching : " << pfcand.eta() << " = " << kidcand->eta() << std::endl;
                auto kdeta = kidcand->eta();
                auto kdphi = kidcand->phi();
                auto kdpfdr = std::sqrt(reco::deltaR2(eleeta, elephi, kdeta, kdphi ));
                if( kdpfdr <= 0.001 ){ ematched = true; if( DEBUG ) std::cout << " --- Matched : " << kdpfdr << std::endl; }
            }//<<>>for( const auto kid : jet.daughterPtrVector() )

	      	if( ematched ){  //  makes electron SC group and adds to SC Group if not already present -  leaves BC groups unchanged
				if( DEBUG ) std:: cout << " ----- Electron Match !!!! " << std::endl;
				bool found(false);
	         	const auto &elesc = electron.superCluster().isNonnull() ? electron.superCluster() : electron.parentSuperCluster();
	         	const auto scptr = elesc.get();

				for( const auto scit : jetSCGroup ){ 
					if( scit == *scptr ){ 
						found = true;
						break;
					}//<<>>if( scit == *scptr ) 
				}//<<>>for( const auto scit : jetSCGroup )

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<for( const auto electron : *electrons_ ){

				if( not found ){ // in list of SC already found
					iel++;
	         		jetSCGroup.push_back(*scptr);
					const auto clusters = elesc->clusters();
	            	const auto &hitsAndFractions = elesc->hitsAndFractions();
					const auto nrh = hitsAndFractions.size();
					if( nrh < minRHcnt ) continue;
 					//const auto nBClusts = elesc->clustersSize();
	            	//std::cout << " -- SC match found with nBClusts: " << nBClusts << " nClusters: "<< clusters.size() << std::endl;
	            	//std::cout << " -- SC has nRecHits: " << nrh << std::endl;
					if( nrh != 0 ){ 
	         			for( const auto &clustptr : clusters ){
							if( not clustptr.isAvailable() ) continue;
	            			const auto clust = clustptr.get();
	            			jetBCGroup.push_back(*clust);
							//std::cout << " --- Adding cluster " << std::endl;
	            		}//<<>>for( const auto &clustptr : clusters )
					}//<<>>if( nrh != 0 )
	         		const auto sce = elesc->energy();
	         		nMatched++;
	         		sum_nrh += nrh;
	         		sum_sce += sce;
	         		sum_phe += electron.energy();

				}//<<>>if( not found )
	         	ematched = false;
	    	}//<<>>if( matched )
	   	}//<<>>for( const auto electron : *electrons_ )
    } // ------------ electrons lock ------------------------------------------------------------------
        // SC group creation finished <<<<< -----------------------------------------------

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

		nJetScMatch.push_back(nMatched);
		sJetScEnergy.push_back(sum_sce);
		sJetScPhEnergy.push_back(sum_phe);

		if( DEBUG ) std::cout << " -- get jetScRhGroup " << std::endl;
		auto jetScRhGroup = getRHGroup( jetSCGroup, minRHenr );
        auto jetScRhIdsGroup = getRhGrpIDs( jetScRhGroup );
		jetScRhIds.push_back( jetScRhIdsGroup );
        //auto jetScRhIdsGroup = getRhGrpIDs( jetScRhGroup );
        //jetScRhIds.push_back(jetScRhIdsGroup);
		if( not isRhGrpEx(jetScRhGroup) ) std::cout << " --- !!!!! jetScRhGroup is not exclusive !!! " << std::endl;  
		//std::cout << " Num SC rechits : " << jetScRhGroup.size() << std::endl;
            
		// get times for SC and BC jet times
		auto sumscrhe = getRhGrpEnr(jetScRhGroup);
		auto scemf = sumscrhe/jet.energy();

		sJetScRhEnergy.push_back(sumscrhe);
        jetScEMF.push_back(scemf);

		if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf ){

			// Get and fill SC times ----------------------------------------------
			if( DEBUG ) std::cout << " --- get jetSCtofTimes for " << jetScRhGroup.size() << " rechits " << std::endl;
	 		auto jetSCtofTimes = getLeadTofRhTime( jetScRhGroup, vtxX, vtxY, vtxZ );
	 		auto jetSCTimeStats = getTimeDistStats( jetSCtofTimes, jetScRhGroup );

			//if( useJetSC ){  //  use for jetSC eigen info
				//std::cout << " --- get jetSC eigen vectors " << std::endl;
            auto jetSCEigen3D = getRhGrpEigen_ieipt( jetSCtofTimes, jetScRhGroup );
    	    auto jetSCEigen2D = getRhGrpEigen_sph( jetSCtofTimes, jetScRhGroup );
            auto impangle = getATan2( hypo( jet.px(), jet.py()), jet.pz());
		
			jetImpactAngle.push_back(impangle);
			jetSc3dEx.push_back(jetSCEigen3D[0]);
            jetSc3dEy.push_back(jetSCEigen3D[1]);
            jetSc3dEz.push_back(jetSCEigen3D[2]);
            jetSc3dEv.push_back(jetSCEigen3D[3]);
            jetSc3dEslope.push_back(jetSCEigen3D[4]);
            jetSc3dEchisp.push_back(jetSCEigen3D[5]);	

            jetSc2dEx.push_back(jetSCEigen2D[0]);
            jetSc2dEy.push_back(jetSCEigen2D[1]);
            jetSc2dEv.push_back(jetSCEigen2D[2]);
            jetSc2dEslope.push_back(jetSCEigen2D[3]);
            jetSc2dEchisp.push_back(jetSCEigen2D[4]);
            jetSc2dEslope2.push_back(jetSCEigen2D[5]);
            jetSc2dEchisp2.push_back(jetSCEigen2D[6]);
            jetSc2dErangle.push_back(jetSCEigen2D[7]);
            jetSc2dEnxsum.push_back(jetSCEigen2D[8]);

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
    	//<<<<if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf )

			jetScRhCnt.push_back(jetScRhGroup.size());

			if( jetSCTimeStats[6] > -28.9 ) nGoodScJets++;

            if( DEBUG ) std::cout << "JetDist Vars : " << jetSCTimeStats[7] << " " << jet.eta() << " " << jet.etaetaMoment() << " " << jet.phiphiMoment();
        	if( DEBUG ) std::cout << " " << jet.etaphiMoment() << " " << jet.maxDistance(); 
			if( DEBUG ) std::cout << " " << jet.constituentPtDistribution() << " " << jet.constituentEtaPhiSpread();
        	if( DEBUG ) std::cout << " " << jet.jetArea() << " " << jet.nCarrying(0.1) << " " << jet.nConstituents() << std::endl; 	
        	if( DEBUG ) std::cout << " - fill vars " << std::endl;

            jetSCMedTime.push_back(jetSCTimeStats[2]);
            jetSCMuTime.push_back(jetSCTimeStats[0]);
            //std::cout << "fill phCSCMuTimeTemp : " << jetSCTimeStats[6] << std::endl;
            jetCSCMuTime.push_back(jetSCTimeStats[6]);
            jetCSCMedTime.push_back(jetSCTimeStats[10]);

			// Get and fill BC times ----------------------------------------------
			// jetBCGroup  -- get BC rechit groups
			if( DEBUG ) std::cout << " --- get jetBCtofTimes for " << jetBCGroup.size() << " clusters " << std::endl;				
			vector<float> bcTimes;
	        vector<rhGroup> bcRhGroups; 
			vector<float> bcEnergies;
			for( auto bc : jetBCGroup ){
	            auto bcRhGroup = getRHGroup( bc, bcMinEnergy );
				if( not isRhGrpEx( bcRhGroup ) ) std::cout << " --- !!!!! bcRhGroup is not exclusive !!! " << std::endl;
				//int bcRhGroupSize = bcRhGroup.size();
				//if( bcRhGroupSize < bcMinRHGrpSize ) continue; //std::cout << " ---- bcRhGroup empty : skip " << std::endl; continue; }
				if( bcRhGroup.size() < bcMinRHGrpSize ) continue; //std::cout << " ---- bcRhGroup empty : skip " << std::endl; continue; }
				bcRhGroups.push_back(bcRhGroup);
				bcEnergies.push_back(bc.energy());
				uInt it(0); for(auto fbc : unusedfbclusts){if(bc.seed() == fbc.seed()){unusedfbclusts.erase(unusedfbclusts.begin()+it); break;}it++;}
			}//<<>>for( auto bc : jetBCGroup )

        	// search for duplicate rhits in cluster rh collections
        	if( DEBUG ) std::cout << " --- premerge bcRhGroup dupCnt: " << getDupCnt(bcRhGroups) << " in " << bcRhGroups.size() <<std::endl;
        	// condense clusters with duplicate rhits in rh collections 
			while(reduceRhGrps(bcRhGroups)){}
        	// search for duplicate rhits in cluster rh collections
        	if( DEBUG ) std::cout << " --- postmerge bcRhGroup dupCnt: " << getDupCnt(bcRhGroups) << " in " << bcRhGroups.size() << std::endl;

			if( DEBUG ) std::cout << " -- Get energy and rh count for matched basics clusters" << std::endl;
			vector<float> bcRhGrpEnergy;
        	vector<uInt> bcRhGrpCnt;
			uInt bcRhCnt(0);
        	for( auto bcRhGroup : bcRhGroups ){
				auto bcRhGrpSize = bcRhGroup.size();
				auto bcRhGrpEnr = getRhGrpEnr(bcRhGroup);
				bcRhCnt += bcRhGrpSize;
               	//std::cout << " ---- get times for bcRhGroup w/ nRecHits: " << bcRhGroup.size() << std::endl;
               	auto bcRhTimes = getLeadTofRhTime( bcRhGroup, vtxX, vtxY, vtxZ );
               	auto bcRhTimeStats = getTimeDistStats( bcRhTimes, bcRhGroup );
               	bcTimes.push_back(bcRhTimeStats[6]);//c mean
				bcRhGrpEnergy.push_back(bcRhGrpEnr);
				bcRhGrpCnt.push_back(bcRhGrpSize);
            }//<<>>for( auto bc : jetBCGroup )

  	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){

			if( DEBUG ) std::cout << " -- Get energy and rh count for unused basics clusters" << std::endl;
			vector<float> fbcRhGrpEnergy;
           	vector<uInt> fbcRhGrpCnt;
			for( auto obc : unusedfbclusts ){
           		float bcMinEnergy(0.0);
            	auto bcRhGroup = getRHGroup( obc, bcMinEnergy );				
				auto bcRhGrpEnr = getRhGrpEnr(bcRhGroup);
				auto bcRhGrpSize = bcRhGroup.size();
            	fbcRhGrpEnergy.push_back(bcRhGrpEnr);
            	fbcRhGrpCnt.push_back(bcRhGrpSize);
			}//<<>>for( auto obc : fbclusts )

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){
				
			// Fill BC based jet times
			if( DEBUG ) std::cout << "Fill BC based jet times ------------------" << std::endl;
			auto nBCTimes = bcTimes.size();
			auto sumbcrhe = vfsum(bcRhGrpEnergy);
			auto bcemf = sumbcrhe/jet.energy();

			jetBcTimesCnt.push_back(nBCTimes);
			jetBcSumRHEnr.push_back(sumbcrhe);
			jetBcEMFr.push_back(bcemf);
			jetBcRhCnt.push_back(bcRhCnt);
			jetBcGrpCnt.push_back(bcRhGroups.size());


			if( nBCTimes != 0 && bcemf > minEmf ){
				auto jetBCTimeStats = getTimeDistStats( bcTimes, bcRhGrpEnergy );
            	//auto jetBCTimeStats = getTimeDistStats( bcTimes, bcEnergies );
				//auto jetBCTimeStats = getTimeDistStats( bcTimes );
				if( jetBCTimeStats[6] > -28.9 ) nGoodBcJets++;
				jetCBCMedTime.push_back(jetBCTimeStats[10]);//c med
	        	jetCBCMuTime.push_back(jetBCTimeStats[6]);//c mu
        	} else { //<<>>if( nBCTimes == 0 )
               		jetCBCMedTime.push_back(-29.75);
               		jetCBCMuTime.push_back(-29.75);
			}//<<>>if( nBCTimes == 0 ) : else

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){

	   	} else { //<<>>if( jetSCGroup.size() >= minRHcnt)

			if( DEBUG ) std::cout << " - jetScRhGroup is < minRHcnt " << std::endl;
	        jetSCMuTime.push_back(-29.25);
	   	   	jetCSCMuTime.push_back(-29.25);
	   	   	jetSCMedTime.push_back(-29.25);
	        jetCSCMedTime.push_back(-29.25);
            jetPhMuTime.push_back(-29.25);
            jetEleMuTime.push_back(-29.25);
	    	jetCBCMuTime.push_back(-29.25);
	    	jetCBCMedTime.push_back(-29.25);
			//jetCSBCMuTime.push_back(-29.25);

            jetImpactAngle.push_back(-99);
            jetSc3dEx.push_back(-99);
            jetSc3dEy.push_back(-99);
            jetSc3dEz.push_back(-99);
            jetSc3dEv.push_back(-99);
            jetSc3dEslope.push_back(-99);
            jetSc3dEchisp.push_back(-99);

            jetSc2dEx.push_back(-99);
            jetSc2dEy.push_back(-99);
            jetSc2dEv.push_back(-99);
            jetSc2dEslope.push_back(-99);
            jetSc2dEchisp.push_back(-99);
            jetSc2dEslope2.push_back(-99);
            jetSc2dEchisp2.push_back(-99);
            jetSc2dErangle.push_back(-99);
            jetSc2dEnxsum.push_back(-99);
            jetScRhCnt.push_back(-99);

            jetBcTimesCnt.push_back(-99);
            jetBcSumRHEnr.push_back(-99);
            jetBcEMFr.push_back(-99);
            jetBcRhCnt.push_back(-99);
            jetBcGrpCnt.push_back(-99);
				
	   	}//<<>>if( jetSCGroup.size() >= minRHcnt) : else			

		if( DEBUG ) std::cout << "Finished with get SC Med times --------------------------- " << std::endl;
			
		if( DEBUG ) std::cout << "Next Jet .......................... " << std::endl; 	
	}//<<>>for ( uInt ijet = 0; ijet < nJets; ijet++ )
	// ** end of jets	***************************************************************************************************

	if( goodJetEvent ) nGoodJetEvents++;

	// -- Fill output trees ------------------------------------------
	if( DEBUG ) std::cout << "---------- Next Event -----" << std::endl;
	outTree->Fill();

	// -- EOFun ------------------------------------------------------
	//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//	 ESHandle<SetupData> pSetup;
	//	 iSetup.get<SetupRecord>().get(pSetup);
	//#endif
}//>>>>void LLPgammaAnalyzer_AOD::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ method called once each job just before starting event loop	------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void LLPgammaAnalyzer_AOD::beginJob(){

	// Global Varibles
	nGoodJetEvents = 0;

   	// Set up DetIdMap
   	//DetIDMap = SetupDetIDs();

	// Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	// Book //histograms
	
	//std::cout << "Histograms Booked" << std::endl;

	// Create output Tree branches -----------------------------
	//     outTree->Branch("",&);

	// Run, Lumi, Event info
    //    ulInt               event; // technically unsigned long long in Event.h...
    //    uInt                run, lumi;

	outTree->Branch("run", &run);
	outTree->Branch("lumi", &lumi);
	outTree->Branch("event", &event, "event/l");

    //    int                 nVtx;
    //    float               vtxX, vtxY, vtxZ;

    outTree->Branch("nVtx", &nVtx);
    outTree->Branch("vtxX", &vtxX);
    outTree->Branch("vtxY", &vtxY);
    outTree->Branch("vtxZ", &vtxZ);

    //    float  metSumEt, metPt, metPx, metPy, metPhi;

    outTree->Branch("metSumEt", &metSumEt);
    outTree->Branch("metPt", &metPt);
    outTree->Branch("metPx", &metPx);
    outTree->Branch("metPy", &metPy);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("metCSumEt", &metCSumEt);
    outTree->Branch("metCPx", &metCPx);
    outTree->Branch("metCPy", &metCPy);

    //    uInt                nJets, nGoodDrJets, nGoodScJets, nGoodBcJets, nUnJets;
    //    float               jetHt;
    // xxx   uInt                nGoodJetEvents;

    outTree->Branch("jetHt", &jetHt);
    outTree->Branch("nJets", &nJets);
//    outTree->Branch("nGoodDrJets", &nGoodDrJets);
//    outTree->Branch("nGoodScJets", &nGoodScJets);
//    outTree->Branch("nGoodBcJets", &nGoodBcJets);
//    outTree->Branch("nUnJets", &nUnJets);

	// Jet info
    //    std::vector<float>  jetE, jetPt, jetPhi, jetEta;
	// xxx   std::vector<float>  jetEtaetaMmt, jetPhiphiMnt, jetEtaphiMnt, jetMaxD, jetConPtDis, jetConEtaPhiSprd, jetArea;
    // xxx   std::vector<uInt>   jetNCarry, jetNConst;
    //    std::vector<float>  jetNHF, jetNEMF, jetCHF, jetCEMF, jetMUF, jetNHM, jetCHM, jetC, jetPHE, jetPHEF;
    //    std::vector<float>  jetELE, jetELEF, jetMUE, jetCharge;
	
	outTree->Branch("jetE", &jetE);
	outTree->Branch("jetPt", &jetPt);
	outTree->Branch("jetEta", &jetEta);
	outTree->Branch("jetPhi", &jetPhi);
//	outTree->Branch("jetNHF", &jetNHF);
//	outTree->Branch("jetNEMF", &jetNEMF);  
//	outTree->Branch("jetCHF", &jetCHF);
//	outTree->Branch("jetCEMF", &jetCEMF);
//	outTree->Branch("jetMUF", &jetMUF);
//	outTree->Branch("jetNHM", &jetNHM);
//	outTree->Branch("jetCHM", &jetCHM);
//	outTree->Branch("jetPHM", &jetPHM);
//	outTree->Branch("jetELM", &jetELM);
	//outTree->Branch("jetC", &jetC);
//	outTree->Branch("jetPHE", &jetPHE);
//	outTree->Branch("jetPHEF", &jetPHEF);
//	outTree->Branch("jetELE", &jetELE);
//	outTree->Branch("jetELEF", &jetELEF);
//  outTree->Branch("jetMUE", &jetMUE);
    //outTree->Branch("jetCharge", &jetCharge);
    outTree->Branch("jetID", &jetID);

    //    std::vector<int>    jetID, njetKids, jetKidOfJet, njetSubs, njetRecHits, jetRecHitOfJet;
    //    std::vector<int>    jetKidPdgID, jetKidCharge, jetKid3Charge, jetPHM, jetELM;
    //    std::vector<uInt>   jetRecHitId;
    //    std::vector<bool>   jetKidLLP;
    //    std::vector<double> jetKidMass, jetKidVx, jetKidVy, jetKidVz;
    //    std::vector<float>  jetKidE, jetKidPt, jetKidPhi, jetKidEta, jetKidTime, jetKidMedTime;

    // xxx    std::vector<float>  jetPhMuTime, jetOOTPhMuTime, jetEleMuTime;

    //    uInt                nCaloJets;
    //    std::vector<float>  cljSeedTOFTimes, cljCMeanTimes;
    //    std::vector<float>  cljBc3dEx, cljBc3dEy, cljBc3dEz, cljBc3dEv, cljBc3dEslope, cljBc3dEchisp;
    //    std::vector<float>  cljBc2dEx, cljBc2dEy, cljBc2dEv, cljBc2dEslope, cljBc2dEchisp;
    //    std::vector<double> cljPt, cljEnergy, cljPhi, cljEta, cljPx, cljPy, cljPz;

/*
    outTree->Branch("cljBcCnt", &cljBcCnt);
    outTree->Branch("nCaloJets", &nCaloJets);
    //outTree->Branch("cljRhIds", &cljRhIds);
    outTree->Branch("cljSeedTOFTime", &cljSeedTOFTime);
    outTree->Branch("cljCMeanTime", &cljCMeanTime);
    //outTree->Branch("cljBc3dEx", &cljBc3dEx);
    //outTree->Branch("cljBc3dEy", &cljBc3dEy);
    //outTree->Branch("cljBc3dEz", &cljBc3dEz);
    //outTree->Branch("cljBc3dEv", &cljBc3dEv);
    //outTree->Branch("cljBc3dEslope", &cljBc3dEslope);
    //outTree->Branch("cljBc3dEchisp", &cljBc3dEchisp);
    //outTree->Branch("cljBc2dEx", &cljBc2dEx);
    //outTree->Branch("cljBc2dEy", &cljBc2dEy);
    //outTree->Branch("cljBc2dEv", &cljBc2dEv);
    //outTree->Branch("cljBc2dEslope", &cljBc2dEslope);
    //outTree->Branch("cljBc2dEchisp", &cljBc2dEchisp);
    outTree->Branch("cljCDrMeanTime", &cljCDrMeanTime);
    outTree->Branch("cljPt", &cljPt);
    outTree->Branch("cljEnergy", &cljEnergy);
    outTree->Branch("cljPhi", &cljPhi);
    outTree->Branch("cljEta", &cljEta);
    outTree->Branch("cljPx", &cljPx);
    outTree->Branch("cljPy", &cljPy);
    outTree->Branch("cljPz", &cljPz);
*/

    //    uInt                nPhotons;
    //    std::vector<float>  phoSeedTOFTimes, phoCMeanTimes;
    //    std::vector<float>  phoSc3dEx, phoSc3dEy, phoSc3dEz, phoSc3dEv, phoSc3dEslope, phoSc3dEchisp;
    //    std::vector<float>  phoSc2dEx, phoSc2dEy, phoSc2dEv, phoSc2dEslope, phoSc2dEchisp;
    //    std::vector<double> phoPt, phoEnergy, phoPhi, phoEta, phoPx, phoPy, phoPz;

    outTree->Branch("nPhotons", &nPhotons);
    outTree->Branch("phoIsOotPho", &phoIsOotPho);
    outTree->Branch("phoExcluded", &phoExcluded);
    outTree->Branch("phoSeedTOFTime", &phoSeedTOFTime);
    outTree->Branch("phoCMeanTime", &phoCMeanTime);
    //outTree->Branch("phoSc3dEx", &phoSc3dEx);
    //outTree->Branch("phoSc3dEy", &phoSc3dEy);
    //outTree->Branch("phoSc3dEz", &phoSc3dEz);
    //outTree->Branch("phoSc3dEv", &phoSc3dEv);
    //outTree->Branch("phoSc3dEslope", &phoSc3dEslope);
    //outTree->Branch("phoSc3dEchisp", &phoSc3dEchisp);
    //outTree->Branch("phoSc2dEx", &phoSc2dEx);
    //outTree->Branch("phoSc2dEy", &phoSc2dEy);
    outTree->Branch("phoSc2dEv", &phoSc2dEv);
    //outTree->Branch("phoSc2dEslope", &phoSc2dEslope);
    //outTree->Branch("phoSc2dEchisp", &phoSc2dEchisp);
    outTree->Branch("phoPt", &phoPt);
    outTree->Branch("phoEnergy", &phoEnergy);
    outTree->Branch("phoPhi", &phoPhi);
    outTree->Branch("phoEta", &phoEta);
    outTree->Branch("phoPx", &phoPx);
    outTree->Branch("phoPy", &phoPy);
    outTree->Branch("phoPz", &phoPz);
    outTree->Branch("phoRhIds", &phoRhIds);

    outTree->Branch("phoIsPFPhoton", &phoIsPFPhoton);
    outTree->Branch("phoIsStdPhoton", &phoIsStdPhoton);
    outTree->Branch("phoHasConTracks", &phoHasConTracks);
    outTree->Branch("phoIsPixelSeed", &phoIsPixelSeed);
    outTree->Branch("phoIsPhoton", &phoIsPhoton);
    outTree->Branch("phoIsEB", &phoIsEB);
    outTree->Branch("phoIsEE", &phoIsEE);
    //outTree->Branch("phoIsEBGap", &phoIsEBGap);
    //outTree->Branch("phoIsEBEtaGap", &phoIsEBEtaGap);
    //outTree->Branch("phoIsEBPhiGap", &phoIsEBPhiGap);
    //outTree->Branch("phoIsEEGap", &phoIsEEGap);
    //outTree->Branch("phoIsEERingGap", &phoIsEERingGap);
    //outTree->Branch("phoIsEEDeeGap", &phoIsEEDeeGap);
    //outTree->Branch("phoIsEBEEGap", &phoIsEBEEGap);

    outTree->Branch("phoHadOverEM", &phoHadOverEM);
    outTree->Branch("phoHadD1OverEM", &phoHadD1OverEM);
    outTree->Branch("phoHadD2OverEM", &phoHadD2OverEM);
    outTree->Branch("phoHadOverEMVaid", &phoHadOverEMVaid);
    outTree->Branch("phohadTowOverEM", &phohadTowOverEM);
    outTree->Branch("phohadTowD10OverEM", &phohadTowD10OverEM);
    outTree->Branch("phohadTowD20OverEM", &phohadTowD20OverEM);
    outTree->Branch("phohadTowOverEMValid", &phohadTowOverEMValid);
    outTree->Branch("phoE1x5", &phoE1x5);
    outTree->Branch("phoE2x5", &phoE2x5);
    outTree->Branch("phoE3x3", &phoE3x3);
    outTree->Branch("phoE5x5", &phoE5x5);
    outTree->Branch("phoMaxEnergyXtal", &phoMaxEnergyXtal);
    outTree->Branch("phoSigmaEtaEta", &phoSigmaEtaEta);
    outTree->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);

    outTree->Branch("phoR1x5", &phoR1x5);
    outTree->Branch("phoR2x5", &phoR2x5);
    outTree->Branch("phoR9", &phoR9);
    outTree->Branch("phoFull5x5_e1x5", &phoFull5x5_e1x5);
    outTree->Branch("phoFull5x5_e2x5", &phoFull5x5_e2x5);
    outTree->Branch("phoFull5x5_e3x3", &phoFull5x5_e3x3);
    outTree->Branch("phoFull5x5_e5x5", &phoFull5x5_e5x5);
    outTree->Branch("phoFull5x5_maxEnergyXtal", &phoFull5x5_maxEnergyXtal);
    outTree->Branch("phoFull5x5_sigmaEtaEta", &phoFull5x5_sigmaEtaEta);
    outTree->Branch("phoFull5x5_sigmaIEtaIEta", &phoFull5x5_sigmaIEtaIEta);
    //outTree->Branch("phoFull5x5_r1x5 ", &phoFull5x5_r1x5);
    //outTree->Branch("phoFull5x5_r2x5 ", &phoFull5x5_r2x5);
    outTree->Branch("phoFull5x5_r9", &phoFull5x5_r9);

    //outTree->Branch("phoNSatXtals", &phoNSatXtals);
    //outTree->Branch("phoIsSeedSat", &phoIsSeedSat);

    //outTree->Branch("phoMipChi2", &phoMipChi2);
    //outTree->Branch("phoMipTotEnergy", &phoMipTotEnergy);
    //outTree->Branch("phoMipSlope", &phoMipSlope);
    //outTree->Branch("phoMipInter", &phoMipInter);

    //outTree->Branch("phoMipNHitCone", &phoMipNHitCone);
    //outTree->Branch("phoMipIsHalo", &phoMipIsHalo);

    outTree->Branch("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04);
    outTree->Branch("phoHcalTwrSumEtConeDR04", &phoHcalTwrSumEtConeDR04);
    outTree->Branch("phoHcalDepth1TowerSumEtConeDR04", &phoHcalDepth1TowerSumEtConeDR04);
    outTree->Branch("phoCalDepth2TowerSumEtConeDR04", &phoCalDepth2TowerSumEtConeDR04);
    outTree->Branch("phoHcalTowerSumEtBcConeDR04", &phoHcalTowerSumEtBcConeDR04);
    outTree->Branch("phoHcalDepth1TowerSumEtBcConeDR04", &phoHcalDepth1TowerSumEtBcConeDR04);
    outTree->Branch("phoHcalDepth2TowerSumEtBcConeDR04", &phoHcalDepth2TowerSumEtBcConeDR04);
    outTree->Branch("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04);
    outTree->Branch("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04);
    outTree->Branch("phoNTrkSolidConeDR04", &phoNTrkSolidConeDR04);
    outTree->Branch("phoNTrkHollowConeDR04", &phoNTrkHollowConeDR04);

    //outTree->Branch("genPhoPt", &genPhoPt);
    //outTree->Branch("genPhoEnergy", &genPhoEnergy);
    //outTree->Branch("genPhoPhi", &genPhoPhi);
    //outTree->Branch("genPhoEta", &genPhoEta);
    //outTree->Branch("genPhoPx", &genPhoPx);
    //outTree->Branch("genPhoPy", &genPhoPy);
    //outTree->Branch("genPhoPz", &genPhoPz);
    //outTree->Branch("genPhoPdgId", &genPhoPdgId);
    //outTree->Branch("genPhoLLP", &genPhoLlp);
    outTree->Branch("genPhoIdx", &genPhoIdx);
    outTree->Branch("genPhoDr", &genPhoDr);

    //    uInt                nElectrons;
    //    std::vector<float>  eleSeedTOFTimes, eleCMeanTimes;
    //    std::vector<float>  eleSc3dEx, eleSc3dEy, eleSc3dEz, eleSc3dEv, eleSc3dEslope, eleSc3dEchisp;
    //    std::vector<float>  eleSc2dEx, eleSc2dEy, eleSc2dEv, eleSc2dEslope, eleSc2dEchisp;
    //    std::vector<double> elePt, eleEnergy, elePhi, eleEta, elePx, elePy, elePz;

/*
    outTree->Branch("nElectrons", &nElectrons);
    //outTree->Branch("eleRhIds", &eleRhIds); 
    outTree->Branch("eleSeedTOFTime", &eleSeedTOFTime);
    outTree->Branch("eleCMeanTime", &eleCMeanTime);
    //outTree->Branch("eleSc3dEx", &eleSc3dEx);
    //outTree->Branch("eleSc3dEy", &eleSc3dEy);
    //outTree->Branch("eleSc3dEz", &eleSc3dEz);
    //outTree->Branch("eleSc3dEv", &eleSc3dEv);
    //outTree->Branch("eleSc3dEslope", &eleSc3dEslope);
    //outTree->Branch("eleSc3dEchisp", &eleSc3dEchisp);
    //outTree->Branch("eleSc2dEx", &eleSc2dEx);
    //outTree->Branch("eleSc2dEy", &eleSc2dEy);
    //outTree->Branch("eleSc2dEv", &eleSc2dEv);
    //outTree->Branch("eleSc2dEslope", &eleSc2dEslope);
    //outTree->Branch("eleSc2dEchisp", &eleSc2dEchisp);
    outTree->Branch("elePt", &elePt);
    outTree->Branch("eleEnergy", &eleEnergy);
    outTree->Branch("elePhi", &elePhi);
    outTree->Branch("eleEta", &eleEta);
    outTree->Branch("elePx", &elePx);
    outTree->Branch("elePy", &elePy);
    outTree->Branch("elePz", &elePz);

    //    std::vector<float>  jetSumEPFrac, jetEPEnergy, jetEMEnergy, jetEMEnrFrac, jetEPEnrFrac;

    outTree->Branch("jetSumEPFrac", &jetSumEPFrac);
    outTree->Branch("jetEPEnergy", &jetEPEnergy);
    outTree->Branch("jetEMEnergy", &jetEMEnergy);
    outTree->Branch("jetEMEnrFrac", &jetEMEnrFrac);
    outTree->Branch("jetEPEnrFrac", &jetEPEnrFrac);

    //    std::vector<float>  jetDrLeadEta, jetDrLeadPhi, jetDrLeadEnr;
    //    std::vector<uInt>   jetDrRhCnt;
    //    std::vector<float>  sJetDrRHEnergy, jetDrEMF;

    //    std::vector<uInt>   nJetScMatch, jetScRhCnt;
    //    std::vector<float>  sJetScEnergy, sJetScPhEnergy, sJetScRhEnergy, jetScEMF;

    //    std::vector<uInt>   jetBcTimesCnt, jetBcRhCnt, jetBcGrpCnt;
    //    std::vector<float>  jetBcSumRHEnr, jetBcEMFr;

    //    std::vector<float>  jetMuTime, jetTimeError, jetTimeRMS, jetMedTime, jetCMuTime, jetCMedTime;
    //    std::vector<float>  jetSCMuTime, jetSCMedTime, jetCSCMuTime, jetCSCMedTime, jetCBCMuTime, jetCBCMedTime;

    outTree->Branch("jetDrLeadEta", &jetDrLeadEta);
    outTree->Branch("jetDrLeadPhi", &jetDrLeadPhi);
    outTree->Branch("jetDrLeadEnr", &jetDrLeadEnr);

    outTree->Branch("sJetDrRHEnergy", &sJetDrRHEnergy);
    outTree->Branch("jetDrEMF", &jetDrEMF);
    outTree->Branch("jetDrRhCnt", &jetDrRhCnt);
    //outTree->Branch("jetDrRhIds", &jetDrRhIds);

    outTree->Branch("jetBcTimesCnt", &jetBcTimesCnt);
    outTree->Branch("jetBcSumRHEnr", &jetBcSumRHEnr);
    outTree->Branch("jetBcEMFr", &jetBcEMFr);
    outTree->Branch("jetBcRhCnt", &jetBcRhCnt);
    outTree->Branch("jetBcGrpCnt", &jetBcGrpCnt);

    outTree->Branch("nJetScMatch", &nJetScMatch);
    outTree->Branch("jetScRhCnt", &jetScRhCnt);
    outTree->Branch("jetScRhIds", &jetScRhIds);
    outTree->Branch("sJetScEnergy", &sJetScEnergy);
    outTree->Branch("sJetScPhEnergy", &sJetScPhEnergy);
    outTree->Branch("sJetScRhEnergy", &sJetScRhEnergy);
    outTree->Branch("jetScEMF", &jetScEMF);

    //outTree->Branch("jetDRMuTime", &jetMuTime);
    //outTree->Branch("jetDRTimeError", &jetTimeError);
    //outTree->Branch("jetDRTimeRMS", &jetTimeRMS);
    //outTree->Branch("jetDRMedTime", &jetMedTime);
    outTree->Branch("jetCDRMuTime", &jetCMuTime);
    outTree->Branch("jetCDRMedTime", &jetCMedTime);
    //outTree->Branch("jetSCMuTime", &jetSCMuTime);
    //outTree->Branch("jetSCMedTime", &jetSCMedTime);
    outTree->Branch("jetCSCMuTime", &jetCSCMuTime);
    outTree->Branch("jetCSCMedTime", &jetCSCMedTime);
    outTree->Branch("jetCBCMuTime", &jetCBCMuTime);
    outTree->Branch("jetCBCMedTime", &jetCBCMedTime);

    //    std::vector<float>  jetGenImpactAngle, jetGenTime, jetGenPt, jetGenEta, jetGenEnergy, jetGenEMFrac, jetGenDrMatch;
    //    std::vector<float>  jetGenTimeVar, jetGenTimeLLP, jetGenLLPPurity, jetGenNextBX, jetGenNKids, jetGenTOF;

    outTree->Branch("jetGenImpactAngle", &jetGenImpactAngle);
    outTree->Branch("jetGenTime", &jetGenTime);
    outTree->Branch("jetGenPt", &jetGenPt);
    outTree->Branch("jetGenEta", &jetGenEnergy);
    outTree->Branch("jetGenEnergy", &jetGenEnergy);
    outTree->Branch("jetGenEMFrac", &jetGenEMFrac);
    outTree->Branch("jetGenDrMatch", &jetGenDrMatch);
    outTree->Branch("jetGenTimeVar", &jetGenTimeVar);
    outTree->Branch("jetGenTimeLLP", &jetGenTimeLLP);
    outTree->Branch("jetGenLLPPurity", &jetGenLLPPurity);
    outTree->Branch("jetGenNextBX", &jetGenNextBX);
    outTree->Branch("jetGenNKids", &jetGenNKids);
    outTree->Branch("jetGenTOF", &jetGenTOF);

    //    std::vector<float>  jetImpactAngle;
    //    std::vector<float>  jetSc3dEx, jetSc3dEy, jetSc3dEz, jetSc3dEv, jetSc3dEslope, jetSc3dEchisp;

    //outTree->Branch("jetImpactAngle", &jetImpactAngle);
    //outTree->Branch("jetSc3dEx", &jetSc3dEx);
    //outTree->Branch("jetSc3dEy", &jetSc3dEy);
    //outTree->Branch("jetSc3dEz", &jetSc3dEz);
    //outTree->Branch("jetSc3dEv", &jetSc3dEv);
    //outTree->Branch("jetSc3dEslope", &jetSc3dEslope);
    //outTree->Branch("jetSc3dEchisp", &jetSc3dEchisp);

    //    std::vector<float>  jetSc2dEx, jetSc2dEy, jetSc2dEv, jetSc2dEslope, jetSc2dEchisp;
    //    std::vector<float>  jetSc2dEslope2, jetSc2dEchisp2, jetSc2dErangle, jetSc2dEnxsum;

    //outTree->Branch("jetSc2dEx", &jetSc2dEx);
    //outTree->Branch("jetSc2dEy", &jetSc2dEy);
    //outTree->Branch("jetSc2dEv", &jetSc2dEv);
    //outTree->Branch("jetSc2dEslope", &jetSc2dEslope);
    //outTree->Branch("jetSc2dEchisp", &jetSc2dEchisp);
    //outTree->Branch("jetSc2dEslope2", &jetSc2dEslope2);
    //outTree->Branch("jetSc2dEchisp2", &jetSc2dEchisp2);
    //outTree->Branch("jetSc2dErangle", &jetSc2dErangle);
    //outTree->Branch("jetSc2dEnxsum", &jetSc2dEnxsum);
*/

    //    uInt nGenParts;
    //    std::vector<double> genPt, genEnergy, genPhi, genEta, genPx, genPy, genPz;
    //    std::vector<int>    genPdgId, genLLP;

    outTree->Branch("nGenParts", &nGenParts);
    outTree->Branch("genPt", &genPt);
    outTree->Branch("genEnergy", &genEnergy);
    outTree->Branch("genPhi", &genPhi);
    outTree->Branch("genEta", &genEta);
    outTree->Branch("genPx", &genPx);
    outTree->Branch("genPy", &genPy);
    outTree->Branch("genPz", &genPz);
    outTree->Branch("genPdgId", &genPdgId);
    outTree->Branch("genLLP", &genLLP);


    //    int nRecHits;
    //    std::vector<float>  rhPosX, rhPosY, rhPosZ, rhPosEta, rhPosPhi;
    //    std::vector<uInt>   rhID, rhXtalI1, rhXtalI2, rhSubdet;
    //    std::vector<float>  rhEnergy, rhTime, rhTimeErr, rhTOF;
    //    std::vector<bool>   rhisOOT, rhisGS6, rhisGS1;
    //    std::vector<float>  rhadcToGeV;
	//    std::vector<float>  rhped12, rhped6, rhped1;
    //    std::vector<float>  rhpedrms12, rhpedrms6, rhpedrms1;

    outTree->Branch("nRecHits", &nRecHits);
    outTree->Branch("rhPosX", &rhPosX);
    outTree->Branch("rhPosY", &rhPosY);
    outTree->Branch("rhPosZ", &rhPosZ);
    outTree->Branch("rhPosEta", &rhPosEta);
    outTree->Branch("rhPosPhi", &rhPosPhi);
    outTree->Branch("rhEnergy", &rhEnergy);
    outTree->Branch("rhTime", &rhTime);
    //outTree->Branch("rhTimeErr", &rhTimeErr);
    outTree->Branch("rhTOF", &rhTOF);
    outTree->Branch("rhID", &rhID);
    //outTree->Branch("rhXtalI1", &rhXtalI1);
    //outTree->Branch("rhXtalI2", &rhXtalI2);
    //outTree->Branch("rhSubdet", &rhSubdet);
    outTree->Branch("rhisOOT", &rhisOOT);
	//outTree->Branch("rhAmplitude", &rhAmplitude); // <<<< ---------   calc and save this todo !!!!1
    //outTree->Branch("rhisGS6", &rhisGS6);
    //outTree->Branch("rhisGS1", &rhisGS1);
    //outTree->Branch("rhadcToGeV", &rhadcToGeV);
    //outTree->Branch("rhped12", &rhped12);
    //outTree->Branch("rhped6", &rhped6);
    //outTree->Branch("rhped1", &rhped1);
    //outTree->Branch("rhpedrms12", &rhpedrms12);
    //outTree->Branch("rhpedrms6", &rhpedrms6);
    //outTree->Branch("rhpedrms1", &rhpedrms1);

}//>>>>void LLPgammaAnalyzer_AOD::beginJob()


// ------------ method called once each job just after ending the event loop	------------
void LLPgammaAnalyzer_AOD::endJob(){

}//>>>>void LLPgammaAnalyzer_AOD::endJob()


// ------------ method fills 'descriptions' with the allowed parameters for the module	------------
void LLPgammaAnalyzer_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
}//>>>>void LLPgammaAnalyzer_AOD::fillDescriptions(edm::ConfigurationDescriptions& descriptions)


//define this as a plug-in
DEFINE_FWK_MODULE(LLPgammaAnalyzer_AOD);
