// -*- C++ -*-
//
// Package:    LLPgammaAnalyzer
// Class:      LLPgammaAnalyzer
//
/**\class LLPgammaAnalyzer LLPgammaAnalyzer.cc LLPgammaAnalyzer/plugins/LLPgammaAnalyzer.cc

 		Description: [one line class summary]

 		Implementation: [Notes on implementation]

*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//----------------------------------------  cc file   --------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------

#include "LLPGamma/LLPgammaAnalyzer/plugins/LLPgammaAnalyzer.hh"
using namespace std;

//#define DEBUG true
#define DEBUG false

//
// constructors and destructor
//
LLPgammaAnalyzer::LLPgammaAnalyzer(const edm::ParameterSet& iConfig) :

// -- declare tags ----------------------------------------------------------

	// flags
	hasGenInfo (iConfig.existsAs<bool>("hasGenInfo")  ? iConfig.getParameter<bool>("hasGenInfo")  : false),

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
{ //<<<< LLPgammaAnalyzer::LLPgammaAnalyzer(const edm::ParameterSet& iConfig) :

	usesResource();
	usesResource("TFileService");

// -- consume tags ------------------------------------------------------------
	if( DEBUG ) std::cout << "In constructor for LLPgammaAnalyzer - tag and tokens" << std::endl;

	// Triggers
	//triggerResultsToken_ 	= consumes<edm::TriggerResults>(triggerResultsTag);
	//triggerObjectsToken_ 	= consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectsTag);

	// tracks 
	tracksToken_				= consumes<std::vector<reco::Track>>(tracksTag);

	// genparticles
	// genpart_token_              = consumes<std::vector<pat::PackedGenParticle>>(genpartTag);

	// pfcandidates
	pfcand_token_        		= consumes<CandidateView>(pfcandTag);

	// vertices
	verticesToken_				= consumes<std::vector<reco::Vertex>>(verticesTag);

	// rho
	//rhoToken_					= consumes<double>(rhoTag);

	// mets
	metsToken_					= consumes<std::vector<pat::MET>>(metsTag);

	// supercluster
	scToken_             		= consumes<reco::SuperClusterCollection>(superClusterCollectionTag);
	ootScToken_          		= consumes<reco::SuperClusterCollection>(ootSuperClusterCollectionTag); 

	// caloClusters
	ccToken_			     	= consumes<std::vector<reco::CaloCluster>>(caloClusterTag);	

	// jets
	jetsToken_					= consumes<std::vector<pat::Jet>>(jetsTag);
	
	// leptons
	electronsToken_				= consumes<std::vector<pat::Electron>>(electronsTag);
	muonsToken_					= consumes<std::vector<pat::Muon>>(muonsTag);

	// rechits
	recHitsEBToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
	recHitsEEToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

	// photons
	gedPhotonsToken_ 			= consumes<std::vector<pat::Photon>>(gedPhotonsTag);
	ootPhotonsToken_ 			= consumes<std::vector<pat::Photon>>(ootPhotonsTag);

// ---------------------------------------------------------------------------------
}//>>>>LLPgammaAnalyzer::LLPgammaAnalyzer(const edm::ParameterSet& iConfig)


LLPgammaAnalyzer::~LLPgammaAnalyzer(){
	///////////////////////////////////////////////////////////////////
	// do anything here that needs to be done at desctruction time   //
	// (e.g. close files, deallocate resources etc.)                 //
	///////////////////////////////////////////////////////////////////
}//>>>>LLPgammaAnalyzer::~LLPgammaAnalyzer()


//
// member functions
//

detIdMap LLPgammaAnalyzer::SetupDetIDs(){

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

}//>>>>detIdMap LLPgammaAnalyzer::SetupDetIDs()


int LLPgammaAnalyzer::getPFJetID(const pat::Jet & jet){ 
	
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

}//>>>>int LLPgammaAnalyzer::getPFJetID(const pat::Jet & jet)


rhGroup LLPgammaAnalyzer::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 ){

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

}//>>>>rhGroup LLPgammaAnalyzer::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 )

rhGroup LLPgammaAnalyzer::getRHGroup( uInt detid ){

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
}//>>>>rhGroup LLPgammaAnalyzer::getRHGroup( uInt detid )

rhGroup LLPgammaAnalyzer::getRHGroup(){

	rhGroup result;
	for (const auto recHit : *recHitsEB_ ) result.push_back(recHit);
	for (const auto recHit : *recHitsEE_ ) result.push_back(recHit);

	return result;

}//>>>>rhGroup LLPgammaAnalyzer::getRHGroup()


rhGroup LLPgammaAnalyzer::getRHGroup( const scGroup superClusterGroup, float minenr ){

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

}//>>>>rhGroup LLPgammaAnalyzer::getRHGroup( const scGroup superClusterGroup, float minenr = 0.0 )

rhGroup LLPgammaAnalyzer::getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax ){

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

}//>>>>rhGroup LLPgammaAnalyzer::getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax )

rhGroup LLPgammaAnalyzer::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 ){

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

}////rhGroup LLPgammaAnalyzer::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 )

	
vector<float>	LLPgammaAnalyzer::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

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

}//>>>>vector<float>  LLPgammaAnalyzer::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )


EcalRecHit LLPgammaAnalyzer::getLeadRh( rhGroup recHits ){
	
	EcalRecHit result;
	float enr(0.0);
	for (const auto recHit : recHits ){
		auto rhenr = recHit.energy();
		if( rhenr < enr ) continue;
		enr = rhenr;
		result = recHit;
	}//<<>>for (const auto recHit : recHits )

	return result;

}//>>>>EcalRecHit LLPgammaAnalyzer::getLeadRh( rhGroup recHits )


vector<float>	LLPgammaAnalyzer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

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

}//>>>>vector<float>  LLPgammaAnalyzer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

vector<float> LLPgammaAnalyzer::getTimeDistStats( vector<float> times ){

	vector<float> results;
	int size = times.size();
	
	if( size == 0 ){		
	   	std::cout << "Time group is empty" << std::endl; 
	   	results.push_back(-9.25);//0
	   	results.push_back(99.9);//1
	   	results.push_back(-9.25);//2
	   	results.push_back(99.9);//3
	   	results.push_back(99.9);//4
	   	results.push_back(9.9);//5
	   	results.push_back(-9.25);//6
	   	results.push_back(99.9);//7
	   	results.push_back(99.9);//8
	   	results.push_back(99.9);//9
	   	results.push_back(-9.25);//10
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
}//>>>> vector<float> LLPgammaAnalyzer::getTimeDistStats( vector<float> times )


vector<float> LLPgammaAnalyzer::getTimeDistStats( vector<float> times, vector<float> wts ){

	vector<float> results;
	int size = times.size();	
	int wtsize = wts.size();
	if( size != wtsize || size == 0 ){ 
		std::cout << "Time & rechit groups not same or empty" << std::endl; 
		results.push_back(-9.25);//0 mu
	   	results.push_back(99.9);//1
	   	results.push_back(-9.25);//2
	   	results.push_back(99.9);//3
	   	results.push_back(99.9);//4
	   	results.push_back(9.9);//5
	   	results.push_back(-9.25);//6
	   	results.push_back(99.9);//7
	   	results.push_back(99.9);//8
	   	results.push_back(99.9);//9
		results.push_back(-9.25);//10
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
}//>>>>vector<float> LLPgammaAnalyzer::getTimeDistStats( vector<float> times, vector<float> weights )

vector<float> LLPgammaAnalyzer::getTimeDistStats( vector<float> times, rhGroup rechits ){

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
}//>>>>vector<float> LLPgammaAnalyzer::getTimeDistStats( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer::getRhGrpEigen( vector<float> xs, vector<float> wts ){
//spherical

    vector<float> results;

    auto ts2 = wsin2( xs, wts );
    auto tc2 = wcos2( xs, wts );
    auto tsc = wsincos( xs, wts );
	//std::cout << "Sin2 : " << ts2 << " Cos2 : " << tc2 << " SinCos : " << tsc << std::endl;
    //std::cout << "Varencies and Means Calculated" << std::endl;

    //TMatrixDSym mat(3,3);
    double array[] = { ts2, tsc, tsc, tc2 };
    //for( auto thing : array ) std::cout << thing;
    //cout << std::endl;
    //TArrayD data(9,array);
    //mat.SetMatrixArray(data.GetArray());
    //mat.SetMatrixArray(array);
    TMatrixDSym mat(2,array);
    //std::cout << "Input matrix created" << std::endl;
    //mat.Print();

    TMatrixDSymEigen eigen(mat);
    const TVectorD& eginVal = eigen.GetEigenValues();
    const TMatrixD& eginVec = eigen.GetEigenVectors();
    TVectorD evaules(eginVal);
    //evaules.ResizeTo(3);
    //std::cout << "Eigin Values and Vectors Found " << std::endl;
    //eginVal.Print();
    //evaules.Print();
    //eginVec.Print();

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
}//<<>>vector<float> LLPgammaAnalyzer::getRhGrpEigen2D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

vector<float> LLPgammaAnalyzer::getRhGrpEigen( vector<float> xs, vector<float> ys, vector<float> wts ){

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
    //for( auto thing : array ) std::cout << thing;
    //cout << std::endl;
    //TArrayD data(9,array);
    //mat.SetMatrixArray(data.GetArray());
    //mat.SetMatrixArray(array);
    TMatrixDSym mat(2,array);
    //std::cout << "Input matrix created" << std::endl;
    //mat.Print();

    TMatrixDSymEigen eigen(mat);
    const TVectorD& eginVal = eigen.GetEigenValues();
    const TMatrixD& eginVec = eigen.GetEigenVectors();
    TVectorD evaules(eginVal);
    //evaules.ResizeTo(3);
    //std::cout << "Eigin Values and Vectors Found " << std::endl;
    //eginVal.Print();
    //evaules.Print();
    //eginVec.Print();

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
}//<<>>vector<float> LLPgammaAnalyzer::getRhGrpEigen2D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

vector<float> LLPgammaAnalyzer::getRhGrpEigen( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts ){
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
	//std::cout << "Varencies and Means Calculated" << std::endl;
    //std::cout << " xx- " << var_x  << " xy- " << var_xy << " xz- " << var_xz << std::endl; 
    //std::cout << " xy- " << var_xy << " yy- " << var_y  << " yz- " << var_yz << std::endl;
    //std::cout << " xz- " << var_xz << " yz- " << var_yz << " zz- " << var_z  << std::endl;

	//TMatrixDSym mat(3,3);
	double array[] = { var_x, var_xy, var_xz, var_xy, var_y, var_yz, var_xz, var_yz, var_z };
    //for( auto thing : array ) std::cout << thing;
	//cout << std::endl;
   	//TArrayD data(9,array);
	//mat.SetMatrixArray(data.GetArray());
    //mat.SetMatrixArray(array);
    TMatrixDSym mat(3,array);
    //std::cout << "Input matrix created" << std::endl;
	//mat.Print();

	TMatrixDSymEigen eigen(mat);
	const TVectorD& eginVal = eigen.GetEigenValues();
	const TMatrixD& eginVec = eigen.GetEigenVectors();
	TVectorD evaules(eginVal);	
    //evaules.ResizeTo(3);
    //std::cout << "Eigin Values and Vectors Found " << std::endl;
    //eginVal.Print();
    //evaules.Print();
    //eginVec.Print();

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

	results.push_back(ex);
    results.push_back(ey);
    results.push_back(ez);
    results.push_back(ev);

	return results;
}//<<>>vector<float> LLPgammaAnalyzer::getRhGrpEigen3D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

vector<float> LLPgammaAnalyzer::getRhGrpEigen_ieipt( vector<float> times, rhGroup rechits ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);
    float tdis(2.2/SOL);

    vector<float> wts;
    vector<float> etas;
    vector<float> phis;
    vector<float> ebtimes;
    vector<float> emepReturn(4,-999);

    //bool ebp = false;
    for( uInt it(0); it < rechits.size(); it++ ){

		const auto recHitID = getRawID(rechits[it]);
		const auto & idinfo = DetIDMap[recHitID];
		if( idinfo.ecal == ECAL::EB ){
    		//auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
        	auto resolution = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
        	wts.push_back(resolution);
        	//const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        	const auto rhEtaPos = tdis*(idinfo.i2);//recHitPos.eta();
			//if( rhEtaPos > 0 ) ebp = true;
        	etas.push_back(rhEtaPos);
            hist1d[66]->Fill(idinfo.i2);
        	//const auto rhPhiPos = tdis*idinfo.i1;//recHitPos.phi();
        	phis.push_back(idinfo.i1);
            hist1d[67]->Fill(idinfo.i1);
			ebtimes.push_back(times[it]);
            hist1d[68]->Fill(times[it]);
		} else { return emepReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
		}//<<>>if( idinfo.ecal == ECAL::EB )

    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

	vector<float> adjphis;
	auto mphi = ceil(meanPhi(phis));
	for( uInt it(0); it < phis.size(); it++ ){ 
		auto dstphi = dPhi(phis[it],mphi);
		adjphis.push_back(tdis*(mphi+dstphi));
	}//<<>>for( auto it : phis )

    auto eigens = getRhGrpEigen( etas, adjphis, ebtimes, wts );

    vector<float> angles;
    auto meta = mean(etas,wts);
    auto adjmphi = mean(adjphis,wts);
    auto mtime = mean(ebtimes,wts);
    //for( uInt it(0); it < ebtimes.size(); it++ ){
    //    float leta = etas[it] - meta;
    //    float lphi = adjphis[it] - amphi;
    //    if( leta == 0 && lphi == 0 ) continue;
    //    //float ltim = ebtimes[it]-mtime;
    //    float angle = getAngle( leta, lphi );
    //    angles.push_back(angle);
    //}//<<>>for( uInt it(0); it < ebtimes.size(); it++ )

    auto eprotangle = getAngle(eigens[0], eigens[1]);
	auto ephypo = hypo(eigens[0], eigens[1]);
    auto etrotangle = getAngle(ephypo, eigens[2]);
    //if( not ebp ){ rotangle += PI; if( rotangle > TWOPI ) rotangle -= TWOPI; }
    float epeignsin = std::sin(eprotangle);
    float epeigncos = std::cos(eprotangle);
    float eteignsin = std::sin(etrotangle);
    float eteigncos = std::cos(etrotangle);
    for( uInt it(0); it < ebtimes.size(); it++ ){

        float leta = etas[it] - meta;
        float lphi = adjphis[it] - adjmphi;
        float ltime = ebtimes[it] - mtime;
		//if( not ebp ) leta *= -1;
        hist1d[75]->Fill(leta);
        hist1d[76]->Fill(lphi);
        hist1d[77]->Fill(ltime);
        auto epxcor = epeigncos*(leta) - epeignsin*(lphi);
        auto epycor = epeignsin*(leta) + epeigncos*(lphi);
        auto etxcor = eteigncos*(epxcor) - eteignsin*(ltime);
        auto etycor = eteignsin*(epxcor) + eteigncos*(ltime);
        //if( not ebp ) xcor *= -1;
        //std::cout << "In getRhGrpEigen_sph w/2 leta " << leta << " : lphi " << lphi << " : xcor " << xcor << " : ycor " << ycor << " : dt " << wts[it] << std::endl;
        //if( abs(wts[it]) < 8 )
        auto fill = ltime*wts[it];
		//if( std::abs(fill) < 1 ){
        hist2d[77]->Fill(epxcor,epycor,fill);
        hist2d[78]->Fill(epxcor,epycor,wts[it]);
        //}//<<>>if( std::abs(fill) < 20 )
        hist2d[79]->Fill(etxcor,etycor,fill);
        hist2d[84]->Fill(etxcor,etycor,wts[it]);
        hist2d[80]->Fill(epxcor,ltime,fill);
        hist2d[85]->Fill(epxcor,ltime,wts[it]);

    }//<<>>for( uInt it(0); it < wts.size(); it++ )

    return eigens;
}//>>>>vector<float> LLPgammaAnalyzer::getRhGrpEigen( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer::getRhGrpEigen_xyz( vector<float> times, rhGroup rechits ){

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
}//>>>>vector<float> LLPgammaAnalyzer::getRhGrpEigen( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer::getRhGrpEigen_sph( vector<float> times, rhGroup rechits ){

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
    vector<float> emepReturn(4,-999);

	//bool ebp = true;
    //bool ebn = true;
	auto nRecHits = rechits.size();
    for( uInt it(0); it < nRecHits; it++ ){

        const auto recHitID = getRawID(rechits[it]);
        const auto & idinfo = DetIDMap[recHitID];
        if( idinfo.ecal == ECAL::EB ){ 
        	//auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
        	//auto wt = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
        	//wts.push_back(wt);
        	//const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        	const auto rhEtaPos = idinfo.i2;//recHitPos.eta();
			//if( rhEtaPos < 10 ) ebp = false;
            //if( rhEtaPos > 10 ) ebn = false;
        	etas.push_back(rhEtaPos);
        	//hist1d[69]->Fill(rhEtaPos);
        	const auto rhPhiPos = idinfo.i1;//recHitPos.phi();
        	//hist1d[70]->Fill(rhPhiPos);
        	phis.push_back(rhPhiPos);
        	ebtimes.push_back(times[it]);	
			auto rhenergy = rechits[it].energy();
			auto resolution = 1/(sq2(N/rhenergy)+2*C*C);
			//energies.push_back(sq2(rhenergy));
            resolutions.push_back(resolution);
			//sumenergy += rhenergy;
			//std::cout << "In getRhGrpEigen_sph w/ rheta " << rhEtaPos << " : rhphi " << rhPhiPos << " : rht " << times[it] << std::endl;		
		} else { return emepReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
		}//<<>>if( idinfo.ecal != ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

	//auto meta = std::ceil(mean(etas,resolution));
    auto meta = mean(etas,resolutions);
    hist1d[84]->Fill(meta);
    //auto mphi = std::ceil(meanPhi(phis));
	//auto mphi = std::ceil(meanPhi(phis,resolution));
    auto mphi = meanPhi(phis,resolutions);
	//vector<float> phiTest{25,50,255,240};
	hist1d[83]->Fill(mphi);      //meanPhi(phiTest));
	auto mtime = mean(ebtimes,resolutions);
    hist1d[85]->Fill(mtime);

	//std::cout << "In getRhGrpEigen_sph w/ meta " << meta << " : mphi " << mphi << " : mt " << mtime << std::endl;

	for( uInt it(0); it < ebtimes.size(); it++ ){

		float leta = etas[it]-meta;
        float lphi = dPhi(phis[it],mphi);
		if( leta == 0 && lphi == 0 ) continue;
        float ltim = ebtimes[it]-mtime;
		float angle = getAngle( leta, lphi );
        //if( meta < 0 ) angle += PI;
		//if( ltim < 0 ) angle += PI;
		//if( angle > twoPI ) angle -= twoPI;
		angles.push_back(angle);
	//std::cout << "In getRhGrpEigen_sph w/1 leta " << leta << " : lphi " << lphi << " : ltim " << ltim << " : mag " << resolutions[it] << " : angle " << angle << std::endl;
		//wts.push_back(ltim*energies[it]/sumenergy);
        wts.push_back(ltim*ltim*resolutions[it]);
        hist1d[69]->Fill(leta);
        hist1d[70]->Fill(lphi);
        hist1d[71]->Fill(ltim);
        hist1d[72]->Fill(angle);
	}//<<>>for( uInt it(0); it < ebtimes.size(); it++ )

	//for( int ia = -50; ia < 51; ia++ ){ for( int ja = -50; ja < 51; ja++ ){ auto a=getAngle(ia,ja); hist1d[74]->Fill(a); hist1d[80]->Fill(std::sin(a));}}
    //for( int ia = 1; ia < 360; ia++ ){ for( int ja = 1; ja < 360; ja++ ){ auto a=dPhi(ia,ja); hist1d[82]->Fill(a);}}

    auto eigens =  getRhGrpEigen( angles, wts );

	//float eigncot = eigens[0]/eigens[1];
    //float eignsin = 1/(std::sqrt(1+eigncot*eigncot));
	//std::cout << "In getRhGrpEigen_sph eigen0 " << eigens[0] << " : eigen1 " << eigens[1] << " sin : " << eignsin << std::endl;
    //float eigncos = std::qrt(1.0 - eignsin*eignsin);
    auto rotangle = getAngle(eigens[0], eigens[1]);
    //hist1d[73]->Fill(rotangle);
	//hist1d[78]->Fill(eigens[0]);
	//hist1d[79]->Fill(eigens[1]);
    //if( not ebp ){ rotangle += PI; if( rotangle > TWOPI ) rotangle -= TWOPI; }
	//hist1d[74]->Fill(rotangle);
	float eignsin = std::sin(rotangle);
	float eigncos = std::cos(rotangle);
	float ltsum(0.0);
	for( uInt it(0); it < wts.size(); it++ ){

        float leta = etas[it] - meta;
        float lphi = dPhi(phis[it],mphi);
		float ltime = ebtimes[it]-mtime;
		auto xcor = eigncos*(leta) - eignsin*(lphi);
		if( xcor > 1 ) ltsum += ltime*resolutions[it];
		
	}//<<>>for( uInt it(0); it < wts.size(); it++ )

	if( ltsum < 0 ){ 
		eigens[0] *= -1; 
		eigens[1] *= -1;
    	rotangle = getAngle(eigens[0], eigens[1]);
		eignsin = std::sin(rotangle);
		eigncos = std::cos(rotangle);
	}//if( ebp && ebn && ltsum < 0 )
	//if( not ebn && ebp ){
    //    eigens[0] *= -1;
    //    eigens[1] *= -1;
    //    rotangle = getAngle(eigens[0], eigens[1]);
    //    eignsin = std::sin(rotangle);
    //    eigncos = std::cos(rotangle);
    //}//if( not ebp && ebn ) 

	hist1d[73]->Fill(rotangle);
	hist1d[78]->Fill(eigens[0]);
	hist1d[79]->Fill(eigens[1]);

	//bool pfcut = eigens[2] < 0.95 && eigens[2] > 0.6;
    bool pfcut = eigens[2] > 0.6 && wts.size() > 3;

    for( uInt it(0); it < wts.size(); it++ ){

        float leta = etas[it] - meta;
        float lphi = dPhi(phis[it],mphi);
        float ltime = ebtimes[it]-mtime;
        auto xcor = eigncos*(leta) - eignsin*(lphi);
        auto ycor = eignsin*(leta) + eigncos*(lphi);
        //if( not ebp ) xcor *= -1;
        //std::cout << "In getRhGrpEigen_sph w/2 leta " << leta << " : lphi " << lphi << " : xcor " << xcor << " : ycor " << ycor << " : dt " << wts[it] << std::endl;
        //if( abs(wts[it]) < 8 )
        auto fill = ltime*resolutions[it];
        //if( std::abs(fill) < 1 ){
        hist2d[73]->Fill(xcor,ycor,fill);
        hist2d[74]->Fill(xcor,ycor,resolutions[it]);
        hist2d[75]->Fill(leta,lphi,fill);
        hist2d[76]->Fill(leta,lphi,resolutions[it]);
		if( pfcut ){
			hist2d[86]->Fill(ycor,ltime,resolutions[it]);
			hist2d[87]->Fill(xcor,ltime,resolutions[it]);
        }//<<>>if( pfcut )
    }//<<>>for( uInt it(0); it < wts.size(); it++ )


	return eigens;
}//>>>>vector<float> LLPgammaAnalyzer::getRhGrpEigen_sph( vector<float> times, rhGroup rechits ){

vector<float> LLPgammaAnalyzer::getRhGrpEigen_ep( vector<float> times, rhGroup rechits ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);
    float tdis(2.2/SOL);

    vector<float> wts;
    vector<float> etas;
    vector<float> phis;
    vector<float> ebtimes;

    for( uInt it(0); it < rechits.size(); it++ ){

        const auto recHitID = getRawID(rechits[it]);
        const auto & idinfo = DetIDMap[recHitID];
        if( idinfo.ecal == ECAL::EB ){
            //auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
            auto wt = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
            wts.push_back(wt);
            //const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
            const auto rhEtaPos = tdis*idinfo.i2;//recHitPos.eta();
            etas.push_back(rhEtaPos);
			hist1d[69]->Fill(rhEtaPos);
            const auto rhPhiPos = tdis*idinfo.i1;//recHitPos.phi();
            hist1d[70]->Fill(rhPhiPos);
            phis.push_back(rhPhiPos);
            ebtimes.push_back(times[it]);
        }//<<>>if( idinfo.ecal == ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    return getRhGrpEigen( etas, phis, wts );
}//>>>>vector<float> LLPgammaAnalyzer::getRhGrpEigen_ep( vector<float> times, rhGroup rechits ){

float LLPgammaAnalyzer::getdt( float t1, float t2 ){

	auto dt = t1 - t2;
	if( dt == 0.0 )	dt = -5.5;// to be commented out in final version
	if( t1 == 0.0 || t2 == 0.0 ) dt = -5.0;// to be copmmented out in final version
	if( t1 < -8.9 || t2 < -8.9 ) dt = -3.75;
	
	return dt;
}//>>>>float LLPgammaAnalyzer::getdt( float t1, float t2 )


void LLPgammaAnalyzer::mrgRhGrp( rhGroup & x, rhGroup & y ){

	//std::cout << " --- In mrgRhGrp " << std::endl;
	bool matched(false);
	for(const auto rhy : y ){ 
		for(const auto rhx : x ){ if( rhMatch(rhx,rhy) ) matched = true; }
		if( not matched ) x.push_back(rhy);
	}//<<>>for(const auto rhx : y )
   //std::cout << " ---- Finished mrgRhGrp " << std::endl;

	return;
}//>>>>rhGroup LLPgammaAnalyzer::mrgRhGrp(const rhGroup x, const rhGroup y)


bool LLPgammaAnalyzer::reduceRhGrps( vector<rhGroup> & x ){

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
}//>>>>vector<rhGroup> LLPgammaAnalyzer::reduceRhGrps(const vector<rhGroup> x)

vector<float> LLPgammaAnalyzer::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  ){

    vector<float> result;
    vector<float> kidtime;
    vector<float> kide;
	vector<float> impact;
    for( auto kid : kids ){

        bool done(false);
        bool first(true);
        auto mom = kid->mother();
        vector<double> stepp{kid->p()};
        vector<double> stepe{kid->energy()};
        vector<double> orignvx;
        vector<double> orignvy;
        vector<double> orignvz;
        bool top(true);
		int steps(0);
        while( not done ){
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

        kide.push_back(stepe[0]);
        //std::cout << " ---- jetGenTime Calc Steps : " << steps << std::endl;
        if( steps == 0 ){ 
			kidtime.push_back(999); 
			impact.push_back(999);
			//std::cout << "Steps are Zero" << std::endl; 
		} else {
            float totaltime(0.0);
            auto destvx(cx);
            auto destvy(cy);
            auto destvz(cz);
            for( int it(0); it < steps; it++ ){
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
			impact.push_back(iangle);
            kidtime.push_back(totaltime);
        }//<<>>if( steps == 0 )

    }//<<>>for( auto kid : kids )

    double ste(0.0);
	double se(0.0);
    double sae(0.0);
	int nKids = kidtime.size();
	//std::cout << " ---- jetGenTime Calc : " << std::endl;
    for( int it(0); it < nKids; it++ ){ 
		if( kidtime[it] != 999 ){
			ste += kidtime[it]*kide[it]; 
			sae += impact[it]*kide[it];
			se += kide[it]; 
			//std::cout << " t: " << kidtime[it] << " e: " << kide[it] << " a: " << impact[it] << std::endl;
		}//<<>>if( kidtime[it] != 999 )
	}//<<>>for( int it(0); it < nKids; it++ )
	//std::cout << " ste: " << ste << " se: " << se << " sae: " << sae << std::endl;
    result.push_back(ste/se);
    result.push_back(sae/se);

    return result;

}//>>>>vector<float> LLPgammaAnalyzer::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  )

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ------------ method called for each event	------------
void LLPgammaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

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

    //GEN PARTICLES
	//if( hasGenInfo ) iEvent.getByToken(genpart_token_, genpart_);

	// VERTICES
	iEvent.getByToken(verticesToken_, vertices_);

	// RHO
	//iEvent.getByToken(rhoToken_, rho_);

	// METS
	//iEvent.getByToken(metsToken_, mets_);

	// SUPERCLUSTERS
	iEvent.getByToken(scToken_, superCluster_);  
	iEvent.getByToken(ootScToken_, ootSuperCluster_);

	// CALOCLUSTERS
	iEvent.getByToken(ccToken_, caloCluster_);

	// JETS
	iEvent.getByToken(jetsToken_, jets_);

	// LEPTONS & PHOTONS
	iEvent.getByToken(electronsToken_, electrons_);
	iEvent.getByToken(muonsToken_, muons_);

	// PHOTONS
	iEvent.getByToken(gedPhotonsToken_, gedPhotons_);
	iEvent.getByToken(ootPhotonsToken_, ootPhotons_);

	// ECAL RECHITS
	iEvent.getByToken(recHitsEBToken_, recHitsEB_);
	iEvent.getByToken(recHitsEEToken_, recHitsEE_);

	// GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
	iSetup.get<CaloGeometryRecord>().get(caloGeo_); 
	barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
	endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap); 

// -- Process Objects ------------------------------------------

// -- event information 

   	const auto run   = iEvent.id().run();
   	const auto lumi  = iEvent.luminosityBlock();
   	const auto event = iEvent.id().event();
    if( DEBUG ) std::cout << "******************************************************************************************************" << std::endl;
	if( DEBUG ) std::cout << "Processing event: " << event << " in run: " << run << " and lumiblock: " << lumi << std::endl;

// -- Process Prime Vertix
	const auto & primevtx = vertices_->front();
	
	auto vtxX = primevtx.position().x();
	auto vtxY = primevtx.position().y();
	auto vtxZ = primevtx.position().z();
	
// ---- process jets --------------------------

// ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
// ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 

	std::vector<pat::Jet> fjets;
   	std::vector<reco::CaloCluster> fbclusts;

	//jets.clear(); 
	//jets.reserve(jets_->size());

	// Selection criteria declartation
	// deltaRmin from ecalelf SC Associtor : 0.2
	float	deltaRminJet	= 0.3;//0.4
	//float deltaRminKid	= 0.15;//0.2
	float 	jetPTmin		= 200.0;
	int 	jetIDmin		= 2; //3;
	float 	jetETAmax		= 1.5; //1.5;
	uInt 	minRHcnt		= 5; //32;
	float 	minRHenr		= 2.0;
   	float 	bcMinEnergy		= 0.667;
   	uInt 	bcMinRHGrpSize	= 3;
	float 	minEmf			= 0.0;//0.2
    //bool    useJetSC        = true; // for cluster eigen calcs
    bool    useJetSC		= false; // for cluster eigen calcs use photon instead

	if( DEBUG ) std::cout << "Filter Jets" << std::endl;
	for(const auto& jet : *jets_ ){ // Filters jet collection & sorts by pt
		
		if (jet.pt() < jetPTmin) continue;
	  	if (std::abs(jet.eta()) > jetETAmax) continue;
	  
	  	const auto jetID = getPFJetID(jet);
	  	if (jetID < jetIDmin) continue;

	  	// save the jets, and then store the ID
	  	fjets.emplace_back(jet);
	  	fjets.back().addUserInt("jetID",jetID);
	  
	  	std::sort(fjets.begin(),fjets.end(),sortByPt);
	}//<<>>for(const auto& jet : *jets_ )
	
	for(const auto& bclust : *caloCluster_ ) fbclusts.push_back(bclust);

	auto nUnJets = (*jets_).size();
	nJets = fjets.size();
	uInt nGoodDrJets(0);
   	uInt nGoodScJets(0);
   	uInt nGoodBcJets(0);
	//set the number of leading jets to skim ( = nJets for all )
	//auto stJets = nJets; 
	if( DEBUG ) std::cout << "Init for Jet Loop with " << nJets << " jets"<< std::endl;

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
    jetEleMuTime.clear();

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
	
	float jetHt(0.0);

//////////  reffrence information
// **** const reco::TrackRefVector& associatedTracks() const;
// **** const reco::PixelClusterTagInfo* tagInfoPixelCluster(const std::string& label = "") const; // from b tagging info methods
// **** reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }

	hist1d[18]->Fill(nJets);
	bool goodJetEvent(false);
	if( DEBUG ) std::cout << "-----------------------------------------------------------" << std::endl;
	if( DEBUG ) std::cout << "Starting Jet Loop for " << nJets << " jets " << std::endl; 
	for ( uInt ijet(0); ijet < nJets; ijet++ ){ 
		// places jet info in output tree
		
		//std::cout << "Processing with jet " << ijet << std::endl;
	   	const auto & jet = fjets[ijet];

		float jetGenTime(-9.9);
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

	   	jetID.push_back(jet.userInt("jetID"));
	   	jetNHF.push_back(jet.neutralHadronEnergyFraction());
	   	jetNEMF.push_back(jet.neutralEmEnergyFraction());
	   	jetCHF.push_back(jet.chargedHadronEnergyFraction());
	   	jetCEMF.push_back(jet.chargedEmEnergyFraction());
	   	jetMUF.push_back(jet.muonEnergyFraction());
	   	jetNHM.push_back(jet.neutralMultiplicity());
		jetCHM.push_back(jet.chargedMultiplicity());
		jetCharge.push_back(jet.jetCharge());

		jetPHE.push_back(jet.photonEnergy());
		jetPHEF.push_back(jet.photonEnergyFraction()); 
		jetELE.push_back(jet.electronEnergy());
		jetELEF.push_back(jet.electronEnergyFraction());
		jetMUE.push_back(jet.muonEnergy());
		jetPHM.push_back(jet.photonMultiplicity());
		jetELM.push_back(jet.electronMultiplicity());

   //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
   
	   	if( DEBUG ) std::cout << "Fill jet pt/phi/eta Histograms" << std::endl;
 		
      	const auto jetepafrac 	= jet.photonEnergyFraction() + jet.electronEnergyFraction();
      	const auto jetepe 		= jet.photonEnergy() + jet.electronEnergy();
		const auto jeteme 		= jet.chargedEmEnergy() + jet.neutralEmEnergy();
      	const auto jetemfrac 	= jeteme/jet.energy();
      	const auto jetepfrac 	= jetepe/jet.energy();

		hist2d[61]->Fill(jetepafrac,jetepfrac);
      	hist2d[62]->Fill(jetepfrac,jetemfrac);

	   	hist1d[12]->Fill(jet.pt());
	   	hist1d[13]->Fill(jet.phi());
	   	hist1d[14]->Fill(jet.eta());

   //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

	   	if( DEBUG ) std::cout << "Getting jet dR rechit group" << std::endl; 
		auto jetDrRhGroup = getRHGroup( jet.eta(), jet.phi(), deltaRminJet, minRHenr ); 
		auto rhCount = jetDrRhGroup.size();

	   	//std::cout << "rhCount is " << rhCount << std::endl;
	   	auto sumdrrhe = getRhGrpEnr( jetDrRhGroup );
		auto dremf = sumdrrhe/jet.energy();
		if( rhCount >= minRHcnt && dremf > minEmf ){

			if( DEBUG ) std::cout << " - get jetDRtofTimes " << std::endl;
	   		auto tofTimes = getLeadTofRhTime( jetDrRhGroup, vtxX, vtxY, vtxZ );
	   		auto leadJetRh = getLeadRh( jetDrRhGroup );
			hist2d[52]->Fill( dremf, jetemfrac );			
			hist2d[55]->Fill( sumdrrhe, jeteme );

	   		if( DEBUG ) std::cout << "Starting RecHit Loop" << std::endl;
			for ( uInt irhg = 0; irhg < rhCount; irhg++){
			
				//std::cout << " -- irhg: " << irhg << " rhCount: " << rhCount << std::endl;
				jetRecHitOfJet.push_back(ijet);
				auto detid = (jetDrRhGroup[irhg]).detid();
				//std::cout << " -- (jetDrRhGroup[irhg]).detid(): " << detid.rawId() << std::endl;
		      	jetRecHitId.push_back(detid.rawId());	
				auto rhtime = tofTimes[irhg];
				//std::cout << " -- tofTimes[irhg]: " << rhtime << std::endl;
		      	hist1d[0]->Fill(rhtime);
				auto rhe = (jetDrRhGroup[irhg]).energy();
		      	//std::cout << " -- jetDrRhGroup[irhg]).energy(): " << rhe << std::endl;
		      	hist2d[38]->Fill(rhtime, rhe);
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
	      	for( auto t : tofTimes ) hist1d[5]->Fill(t);
			auto jetTimeStats = getTimeDistStats( tofTimes, jetDrRhGroup ); 
			auto jmutime = jetTimeStats[0];
		   	auto jterr = jetTimeStats[1];
		   	auto jtrms = jetTimeStats[4];
			auto jmedtime = jetTimeStats[2];
	      	//auto mederr = jetTimeStats[3];
			auto jcmutime = jetTimeStats[6];
	      	auto jcmedtime = jetTimeStats[10];

		   	njetRecHits.push_back(rhCount);
		   	jetMuTime.push_back(jmutime);
		   	jetTimeError.push_back(jterr);
		   	jetTimeRMS.push_back(jtrms);
		   	jetMedTime.push_back(jmedtime);
	      	jetCMuTime.push_back(jcmutime);	
	      	jetCMedTime.push_back(jcmedtime);

			if( jcmutime > -8.9 ) nGoodDrJets++;			

		   	hist1d[29]->Fill(jmutime);
		   	hist1d[1]->Fill(rhCount);
		   	hist1d[2]->Fill(jterr);
		   	hist1d[3]->Fill(jtrms);
		   	hist1d[4]->Fill(jmedtime);
	      	hist1d[6]->Fill(jcmutime);
	      	hist1d[7]->Fill(jcmedtime);

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<if( rhCount >= minRHcnt ){

		   	if( DEBUG ) std::cout << "Filling 2D Histos" << std::endl;
	
		   	hist2d[1]->Fill(jmutime,jet.pt());
		   	hist2d[2]->Fill(jmutime,jet.userInt("jetID"));
		   	hist2d[3]->Fill(jmutime,jet.neutralHadronEnergyFraction());
		   	hist2d[4]->Fill(jmutime,jet.chargedHadronEnergyFraction());
		   	hist2d[5]->Fill(jmutime,jet.neutralEmEnergyFraction());
		   	hist2d[6]->Fill(jmutime,jet.chargedEmEnergyFraction());
		   	hist2d[7]->Fill(jmutime,jet.muonEnergyFraction());
		   	hist2d[8]->Fill(jmutime,jet.neutralMultiplicity());
		   	hist2d[9]->Fill(jmutime,jet.chargedMultiplicity());
	
		   	hist2d[10]->Fill(jmutime,jmedtime);
		   	hist2d[24]->Fill(jmutime,rhCount);
		   	hist2d[25]->Fill(jmedtime,rhCount);
		   	hist2d[11]->Fill(jmutime,jtrms);
		   	hist2d[12]->Fill(jmutime,jterr);
		
		   	hist2d[32]->Fill(jmutime,sc_eta);
		   	hist2d[33]->Fill(jmutime,sc_phi);
		   	hist2d[34]->Fill(jmutime,sc_enr);
		   	hist2d[35]->Fill(jmedtime,sc_eta);
		   	hist2d[36]->Fill(jmedtime,sc_phi);
		   	hist2d[37]->Fill(jmedtime,sc_enr);
	
	
		   	hist2d[13]->Fill(jmedtime,jet.pt());
		   	hist2d[14]->Fill(jmedtime,jet.userInt("jetID"));
		   	hist2d[15]->Fill(jmedtime,jet.neutralHadronEnergyFraction());
		   	hist2d[16]->Fill(jmedtime,jet.chargedHadronEnergyFraction());
		   	hist2d[17]->Fill(jmedtime,jet.neutralEmEnergyFraction());
		   	hist2d[18]->Fill(jmedtime,jet.chargedEmEnergyFraction());
		   	hist2d[19]->Fill(jmedtime,jet.muonEnergyFraction());
		   	hist2d[20]->Fill(jmedtime,jet.neutralMultiplicity());
		   	hist2d[21]->Fill(jmedtime,jet.chargedMultiplicity());

		} else { //<<>>if( rhCount > minRHcnt && dremf > minEmf )

			njetRecHits.push_back(0);
	      	jetMuTime.push_back(-9.25);
	      	jetTimeError.push_back(99.9);
	      	jetTimeRMS.push_back(99.9);
	      	jetMedTime.push_back(-9.25);
	      	jetCMuTime.push_back(-9.25);
	      	jetCMedTime.push_back(-9.25);
			
		}//<<>>if( rhCount > minRHcnt && dremf > minEmf ) : else
		
	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

		// GenJet Info for MC  -------------------------------------------------------------------

        if( DEBUG ) std::cout << "Getting jetGenParton Information" << std::endl;
        //const reco::GenParticle * jetGenParton(0);
        if( hasGenInfo ){
            if( DEBUG ) std::cout << " -- Pulling jet gen info " << std::endl;
            auto jetGenParton = jet.genParton();
            auto jetGenJet = jet.genJet();
			if( rhCount > 0 ){
			if( jetGenJet ){
            	if( DEBUG ) std::cout << " ---- jetGenParton : " << jetGenParton << " genJet : " << jetGenJet << std::endl;
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

            	//auto nMother = jetGenJet->numberOfMothers();
            	//auto nDaughter = jetGenJet->numberOfDaughters();
            	auto nSources = jetGenJet->numberOfSourceCandidatePtrs();
            	if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;
            	//std::cout << " - jetGenJet mothers : " << nMother << " daughters : " << nDaughter << " sources : " << nSources << std::endl;
            	if( DEBUG ) std::cout << " - jetGenJet srcs : " << nSources << " PV (" << vtxX << "," << vtxY << "," << vtxZ << ")" << std::endl;
            	auto kids = jetGenJet->daughterPtrVector();
            	//std::cout << bigKidChase( kids, vtxX ) << std::endl;
            	//kidChase( kids, vtxX, vtxY, vtxZ );
				auto leadJetRh = getLeadRh( jetDrRhGroup );
            	auto leadJetRhId = leadJetRh.detid();
            	auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
				auto cx = leadJetRhIdPos.x();
				auto cy = leadJetRhIdPos.y();
				auto cz = leadJetRhIdPos.z();
				auto tofcor = hypo( cx, cy, cz )/SOL;
				if( DEBUG ) kidChase( kids, vtxX, vtxY, vtxZ );
            	auto genTime = kidTOFChain( kids, cx, cy, cz );
				jetGenTime = genTime[0]-tofcor;
				auto impactAngle = genTime[1];
				if( DEBUG ) std::cout << " - jetGenJet GenTime : " << jetGenTime << " rhPos: " << cx << "," << cy << "," << cz << " Angle: " << impactAngle << std::endl;
				hist1d[88]->Fill(jetGenTime);
                hist1d[90]->Fill(impactAngle);
			}//<<>>if( jetGenJet )
			else { if( DEBUG ) std::cout << " - jetGenJet GenTime : jetGenJet == 0 " << std::endl; jetGenTime = -50.0; }
			}//<<>>if( rhCount >= minRHcnt )

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
		//<<<<if( hasGenInfo )
    
			else if( DEBUG ) std::cout << " - jetGenJet GenTime : rhCount == 0 " << std::endl;
			if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;
            //for( auto kid : kids ){
            //  std::string depth(" --");
            //    std::cout << " -- kid > pdgID : " << kid->pdgId() << " pt : " << kid->pt() << " vtx (" << kid->vx() << "," << kid->vy() << "," << kid->vz() << ")";
            //  std::cout << " nMothers : " << kid->numberOfMothers() << std::endl;
            //  motherChase( kid.get(), depth );
            //}//<<>>for( auto kid : kids )

        }//<<>>if( hasGenInfo )


		// Super Cluster group	-----------------------------------------------
		
		if( DEBUG ) std::cout << " - Get jet SC Group " << std::endl;
	   	int iph(0);
	   	bool matched(false);
		int nMatched(0);
		int sum_nrh(0);
	   	float sum_sce(0.0);
	   	float sum_phe(0.0);
		scGroup jetSCGroup;
        scGroup jetPhSCGroup;
        scGroup jetEleSCGroup;
		bcGroup jetBCGroup;
		vector<float> phEnergy;
        vector<float> phDr;
        vector<float> eleEnergy;
        vector<float> eleDr;

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

	   	for( const auto photon : *gedPhotons_ ){

	      	//std::cout << "Proccesssing New Photon :" << std::endl;
	      	edm::RefVector<pat::PackedCandidateCollection> passociated =  photon.associatedPackedPFCandidates();
	      	for( uInt ipcp = 0; ipcp < passociated.size(); ipcp++ ) {
	         	//std::cout << "Processing asc pfcand # " << ipcp << std::endl;
	         	edm::Ptr<pat::PackedCandidate> passociatedPtr = edm::refToPtr( passociated[ipcp] );
	         	const auto *ascpacked_cand = passociatedPtr.get();
	            //int ijk(0);
	            for( const auto kid : jet.daughterPtrVector() ){
	            	//std::cout << "Proccesssing New jetKid :" << std::endl;
	               	auto kidcand = pfcands_->ptrAt(kid.key());
	               	const auto *packed_cand = dynamic_cast<const pat::PackedCandidate *>(kidcand.get());
	            	if( ascpacked_cand == packed_cand ){

	        			//const auto pce = ascpacked_cand->energy();
	                    //std::cout << "Match found at ph: " << iph << " asc: " << ipcp << " jet: "; 
						//std::cout << ijet << "/" << nJets	<< " kid: " << ijk; // << std::endl;
	                    //std::cout << " pce: " << pce << std::endl;
	                    matched = true;

	                }//<<>>if( ascpacked_cand == packed_cand )
	               	//ijk++;
	        	}//<<>>for( const auto kid : jet.daughterPtrVector() )
			}//<<>>for( uInt ipcp = 0; ipcp < passociated.size(); ipcp++ )

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<for( const auto photon : *gedPhotons_ ){ 
  
	      	if( matched ){
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
	         	matched = false;

                jetPhSCGroup.push_back(*scptr);
				phEnergy.push_back(photon.energy());
				phDr.push_back(std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), photon.eta(), photon.phi())));
	    	}//<<>>if( matched )

		}//<<>>for( const auto photon : *gedPhotons_ ) 

	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
		
		int iel(0);
		matched = false;
	   	for( const auto electron : *electrons_ ){
	    	edm::RefVector<pat::PackedCandidateCollection> eassociated = electron.associatedPackedPFCandidates();
	      	for( uInt ipce = 0; ipce < eassociated.size(); ipce++ ) {
	         	edm::Ptr<pat::PackedCandidate> eassociatedPtr = edm::refToPtr( eassociated[ipce] );
	         	const auto *ascpacked_cand = eassociatedPtr.get();
	            //int ijk(0);
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
	               	//ijk++;
	        	}//<<>>for( const auto kid : jet.daughterPtrVector() )
	    	}//<<>>for( uInt ipce = 0; ipce < eassociated.size(); ipce++ )

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
      	//<<<<for( const auto electron : *electrons_ ){

	      	if( matched ){  //  modified to remove electrons from SC group instead of add them to it -  leaves BC groups unchanged
				bool found(false);
	         	const auto &elesc = electron.superCluster().isNonnull() ? electron.superCluster() : electron.parentSuperCluster();
	         	const auto scptr = elesc.get();

                jetEleSCGroup.push_back(*scptr);
				eleEnergy.push_back(electron.energy());
				eleDr.push_back(std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), electron.eta(), electron.phi())));

				for( const auto scit : jetSCGroup ){ 
					if( scit == *scptr ){ 
						found = true;
						break;
					}//<<>>if( scit == *scptr ) 
				}//<<>>for( const auto scit : jetSCGroup )
				//int nscit = 0;
                //for( const auto scit : jetPhSCGroup ){
                //    if( scit == *scptr ){
                //        jetPhSCGroup.erase(jetPhSCGroup.begin()+nscit);
				//		phEnergy.erase(phEnergy.begin()+nscit);
				//		phDr.erase(phDr.begin()+nscit);
				//		break;
                //  }//<<>>for( const auto scit : jetPhSCGroup )
                //    nscit++; 
                //}//<<>>for( const auto scit : jetSCGroup )

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<for( const auto electron : *electrons_ ){

				if( not found ){ // in list of SC already found
					iel++;
	         		jetSCGroup.push_back(*scptr);
					const auto clusters = elesc->clusters();
	            	const auto &hitsAndFractions = elesc->hitsAndFractions();
					const auto nrh = hitsAndFractions.size();
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
	         	matched = false;
	    	}//<<>>if( matched )
	   	}//<<>>for( const auto electron : *electrons_ )
        // SC group creation finished <<<<< -----------------------------------------------

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )

	   	hist1d[43]->Fill( nMatched ); // # of SC matched to a jet
	    hist2d[40]->Fill( jet.energy(), sum_sce );
	    hist2d[41]->Fill( jet.energy(), sum_phe );
	    hist2d[42]->Fill( sum_phe, sum_sce );

		if( DEBUG ) std::cout << " -- get jetScRhGroup " << std::endl;
		auto jetScRhGroup = getRHGroup( jetSCGroup, minRHenr );
		if( not isRhGrpEx(jetScRhGroup) ) std::cout << " --- !!!!! jetScRhGroup is not exclusive !!! " << std::endl;  
		//std::cout << " Num SC rechits : " << jetScRhGroup.size() << std::endl;
            
		// get times for SC and BC jet times
		auto sumscrhe = getRhGrpEnr(jetScRhGroup);
		auto scemf = sumscrhe/jet.energy();
		if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf ){

			// Get and fill SC times ----------------------------------------------
			//std::cout << " --- get jetSCtofTimes " << std::endl;
	 		auto jetSCtofTimes = getLeadTofRhTime( jetScRhGroup, vtxX, vtxY, vtxZ );
	 		auto jetSCTimeStats = getTimeDistStats( jetSCtofTimes, jetScRhGroup );

			if( useJetSC ){  //  use for jetSC eigen info
			//std::cout << " --- get jetSC eigen vectors " << std::endl;
            auto jetSCEigen3D = getRhGrpEigen_ieipt( jetSCtofTimes, jetScRhGroup );
    	    auto jetSCEigen2D = getRhGrpEigen_sph( jetSCtofTimes, jetScRhGroup );
			if( jetSCEigen3D[0] != -999 ){
				auto epanlge = getAngle( jetSCEigen3D[0], jetSCEigen3D[1] );
            	hist1d[63]->Fill(epanlge);//etaphi angle
				auto ephypo3D = hypo( jetSCEigen3D[0], jetSCEigen3D[1] );
            	auto etanlge = getAngle( ephypo3D, jetSCEigen3D[2] );
            	hist1d[64]->Fill(etanlge);//etatim angle
            	hist2d[82]->Fill( jetSCEigen3D[0], jetSCEigen3D[1] );
            	hist2d[83]->Fill( jetSCEigen3D[0], jetSCEigen3D[2] );
				hist1d[86]->Fill(jetSCEigen3D[3]);
			}//<<>>if( jetSCEigen3D[0] != -999 )
			if( jetSCEigen2D[0] != -999 ){
				auto sphanlge = getAngle( jetSCEigen2D[0], jetSCEigen2D[1] );
		    	hist1d[65]->Fill(sphanlge);//eliptical angle
            	hist2d[81]->Fill( jetSCEigen2D[0], jetSCEigen2D[1] );
            	hist1d[81]->Fill(jetSCEigen2D[2]);
				if( jetScRhGroup.size() > 20 ) hist1d[87]->Fill(jetSCEigen2D[2]);
				hist2d[88]->Fill(sphanlge, jetSCEigen2D[2]);
            }//<<>>if( jetSCEigen2D[0] != -999 )
			}//if( true/false ) on/off switch for jetSC or PhotonSC

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
    	//<<<<if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf )

			if( jetSCTimeStats[6] > -8.9 ) nGoodScJets++;
            //std::cout << " - fill hists " << std::endl;
            for( auto t : jetSCtofTimes ) hist1d[47]->Fill(t);
            //if( rhCount >= minRHcnt ) 
            hist2d[39]->Fill( rhCount, jetScRhGroup.size() );
            hist1d[44]->Fill(jetSCTimeStats[2]);//median
            hist1d[45]->Fill(jetSCTimeStats[0]);//mean
            hist1d[46]->Fill(jetSCTimeStats[4]);//rms
            hist1d[50]->Fill(jetSCTimeStats[5]);//skew
            hist1d[8]->Fill(jetSCTimeStats[6]);//c mean
            hist1d[9]->Fill(jetSCTimeStats[10]);//c med   

			hist1d[89]->Fill(jetSCTimeStats[6]-jetGenTime);

        	hist2d[53]->Fill( scemf, jetemfrac );
			hist2d[56]->Fill( sumscrhe, jeteme );
			hist2d[59]->Fill( dremf, scemf );

			hist2d[63]->Fill( jet.eta(), jetSCTimeStats[7] );  
            hist2d[64]->Fill( jet.etaetaMoment(), jetSCTimeStats[7] );
            hist2d[65]->Fill( jet.phiphiMoment(), jetSCTimeStats[7] );
            hist2d[66]->Fill( jet.etaphiMoment(), jetSCTimeStats[7] );
            hist2d[67]->Fill( jet.maxDistance(), jetSCTimeStats[7] );
            hist2d[68]->Fill( jet.constituentPtDistribution(), jetSCTimeStats[7] );
            hist2d[69]->Fill( jet.constituentEtaPhiSpread(), jetSCTimeStats[7] );
            hist2d[70]->Fill( jet.jetArea(), jetSCTimeStats[7] );
            hist2d[71]->Fill( jet.nCarrying(0.1), jetSCTimeStats[7] );
            hist2d[72]->Fill( jet.nConstituents(), jetSCTimeStats[7] );
            if( DEBUG ) std::cout << "JetDist Vars : " << jetSCTimeStats[7] << " " << jet.eta() << " " << jet.etaetaMoment() << " " << jet.phiphiMoment();
        	if( DEBUG ) std::cout << " " << jet.etaphiMoment() << " " << jet.maxDistance(); 
			if( DEBUG ) std::cout << " " << jet.constituentPtDistribution() << " " << jet.constituentEtaPhiSpread();
        	if( DEBUG ) std::cout << " " << jet.jetArea() << " " << jet.nCarrying(0.1) << " " << jet.nConstituents() << std::endl; 	

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf )

        	if( DEBUG ) std::cout << " - fill vars " << std::endl;
            jetSCMedTime.push_back(jetSCTimeStats[2]);
            jetSCMuTime.push_back(jetSCTimeStats[0]);
            //std::cout << "fill phCSCMuTimeTemp : " << jetSCTimeStats[6] << std::endl;
            jetCSCMuTime.push_back(jetSCTimeStats[6]);
            jetCSCMedTime.push_back(jetSCTimeStats[10]);

            if( jetPhSCGroup.size() != 0 ){
                auto jetPhRhGroup = getRHGroup( jetPhSCGroup, minRHenr, phEnergy, phDr, max(phEnergy) );
                if( not isRhGrpEx( jetPhRhGroup ) ) std::cout << " --- !!!!! jetPhRhGroup is not exclusive !!! " << std::endl;
                //std::cout << " Num Ph rechits : " << jetPhRhGroup.size() << std::endl;
               	if( jetPhRhGroup.size() >= minRHcnt ){

              		//std::cout << "getting jetPhtofTimes with jetPhRhGroup size : " << jetPhRhGroup.size() << std::endl;
                    auto jetPhtofTimes = getLeadTofRhTime( jetPhRhGroup, vtxX, vtxY, vtxZ );
                    auto jetPhTimeStats = getTimeDistStats( jetPhtofTimes, jetPhRhGroup );

            		if( not useJetSC ){ // use Photon SC eigen info
            		//std::cout << " --- get jetSC eigen vectors " << std::endl;
            		auto jetSCEigen3D = getRhGrpEigen_ieipt( jetPhtofTimes, jetPhRhGroup );
            		auto jetSCEigen2D = getRhGrpEigen_sph( jetPhtofTimes, jetPhRhGroup );
            		if( jetSCEigen3D[0] != -999 ){
                		auto epanlge = getAngle( jetSCEigen3D[0], jetSCEigen3D[1] );
                		hist1d[63]->Fill(epanlge);//etaphi angle
                		auto ephypo3D = hypo( jetSCEigen3D[0], jetSCEigen3D[1] );
                		auto etanlge = getAngle( ephypo3D, jetSCEigen3D[2] );
                		hist1d[64]->Fill(etanlge);//etatim angle
                		hist2d[82]->Fill( jetSCEigen3D[0], jetSCEigen3D[1] );
                		hist2d[83]->Fill( jetSCEigen3D[0], jetSCEigen3D[2] );
                		hist1d[86]->Fill(jetSCEigen3D[3]);
            		}//<<>>if( jetSCEigen3D[0] != -999 )
            		if( jetSCEigen2D[0] != -999 ){
                		auto sphanlge = getAngle( jetSCEigen2D[0], jetSCEigen2D[1] );
                		hist1d[65]->Fill(sphanlge);//eliptical angle
                		hist2d[81]->Fill( jetSCEigen2D[0], jetSCEigen2D[1] );
                		hist1d[81]->Fill(jetSCEigen2D[2]);
                		if( jetScRhGroup.size() > 20 ) hist1d[87]->Fill(jetSCEigen2D[2]);
                		hist2d[88]->Fill(sphanlge, jetSCEigen2D[2]);
           		 	}//<<>>if( jetSCEigen2D[0] != -999 )
          			}//if( ture/false ) on/off switch for jetSC or PhotonSC
  
                    hist1d[57]->Fill(jetPhTimeStats[6]);//c mean
                    hist1d[58]->Fill(jetPhTimeStats[10]);//c med  
                
                    jetPhMuTime.push_back(jetPhTimeStats[6]);
                    //std::cout << "fill phMuTimeTemp : " << jetPhTimeStats[6] << std::endl;

                } else {
                    jetPhMuTime.push_back( -9.75 );
                    if( DEBUG ) std::cout << " - jetPhRhGroup is < minRHcnt " << std::endl;
                }//<<>>if( jetPhRhGroup.size() >= minRHcnt )
			} else {
				jetPhMuTime.push_back( -9.25 );
            }//<<>>if( jetPhSCGroup.size() != 0 )

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt && scemf > minEmf )
        
            if( jetEleSCGroup.size() != 0 ){
            	auto jetEleRhGroup = getRHGroup( jetEleSCGroup, minRHenr, eleEnergy, eleDr, max(eleEnergy) );
            	if( not isRhGrpEx( jetEleRhGroup ) ) std::cout << " --- !!!!! jetEleRhGroup is not exclusive !!! " << std::endl;
            	//std::cout << " Num Ph rechits : " << jetPhRhGroup.size() << std::endl;
                if( jetEleRhGroup.size() >= minRHcnt ){

                    //std::cout << "getting jetPhtofTimes with jetPhRhGroup size : " << jetPhRhGroup.size() << std::endl;
                    auto jetEletofTimes = getLeadTofRhTime( jetEleRhGroup, vtxX, vtxY, vtxZ );
                    auto jetEleTimeStats = getTimeDistStats( jetEletofTimes, jetEleRhGroup );

                    hist1d[61]->Fill(jetEleTimeStats[6]);//c mean
                    hist1d[62]->Fill(jetEleTimeStats[10]);//c med  

                    jetEleMuTime.push_back(jetEleTimeStats[6]);
                    //std::cout << "fill phMuTimeTemp : " << jetPhTimeStats[6] << std::endl;

                } else {
                    jetEleMuTime.push_back( -9.75 );
                    if( DEBUG ) std::cout << " - jetEleRhGroup is < minRHcnt " << std::endl;
                }//<<>>if( jetPhRhGroup.size() >= minRHcnt )
			} else {
				jetEleMuTime.push_back( -9.25 );
            }//<<>>if( jetEleSCGroup.size() != 0 )

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){

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
				uInt it(0); for(auto fbc : fbclusts){if(bc.seed() == fbc.seed()){fbclusts.erase(fbclusts.begin()+it); break;}it++;}
			}//<<>>for( auto bc : jetBCGroup )
			hist1d[21]->Fill(bcRhGroups.size());

        	// search for duplicate rhits in cluster rh collections
        	if( DEBUG ) std::cout << " --- premerge bcRhGroup dupCnt: " << getDupCnt(bcRhGroups) << " in " << bcRhGroups.size() <<std::endl;
        	hist1d[20]->Fill(getDupCnt(bcRhGroups));
        	// condense clusters with duplicate rhits in rh collections 
			while(reduceRhGrps(bcRhGroups)){}
        	// search for duplicate rhits in cluster rh collections
        	if( DEBUG ) std::cout << " --- postmerge bcRhGroup dupCnt: " << getDupCnt(bcRhGroups) << " in " << bcRhGroups.size() << std::endl;
			hist1d[22]->Fill(bcRhGroups.size());

			if( DEBUG ) std::cout << " -- Get energy and rh count for matched basics clusters" << std::endl;
			vector<float> bcRhGrpEnergy;
        	vector<uInt> bcRhGrpCnt;
			uInt bcRhCnt(0);

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){

        	for( auto bcRhGroup : bcRhGroups ){
				auto bcRhGrpSize = bcRhGroup.size();
				auto bcRhGrpEnr = getRhGrpEnr(bcRhGroup);
               	hist2d[47]->Fill(bcRhGrpSize, bcRhGrpEnr);
				bcRhCnt += bcRhGrpSize;
               	//std::cout << " ---- get times for bcRhGroup w/ nRecHits: " << bcRhGroup.size() << std::endl;
               	auto bcRhTimes = getLeadTofRhTime( bcRhGroup, vtxX, vtxY, vtxZ );
               	auto bcRhTimeStats = getTimeDistStats( bcRhTimes, bcRhGroup );
               	bcTimes.push_back(bcRhTimeStats[6]);//c mean
				bcRhGrpEnergy.push_back(bcRhGrpEnr);
				bcRhGrpCnt.push_back(bcRhGrpSize);
				if( bcRhGrpSize == 1 ) hist1d[54]->Fill(bcRhGrpEnr/jet.energy());
				else if( bcRhGrpSize > 1 ) hist1d[56]->Fill(bcRhGrpEnr/jet.energy());
				for(const auto rh : bcRhGroup ) hist2d[49]->Fill(rh.time(), rh.energy());
            }//<<>>for( auto bc : jetBCGroup )

  	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){

			if( DEBUG ) std::cout << " -- Get energy and rh count for unused basics clusters" << std::endl;
			vector<float> fbcRhGrpEnergy;
           	vector<uInt> fbcRhGrpCnt;
			for( auto obc : fbclusts ){
           		float bcMinEnergy(0.0);
            	auto bcRhGroup = getRHGroup( obc, bcMinEnergy );				
				auto bcRhGrpEnr = getRhGrpEnr(bcRhGroup);
				auto bcRhGrpSize = bcRhGroup.size();
				hist2d[48]->Fill(bcRhGrpSize, bcRhGrpEnr);
            	fbcRhGrpEnergy.push_back(bcRhGrpEnr);
            	fbcRhGrpCnt.push_back(bcRhGrpSize);
			}//<<>>for( auto obc : fbclusts )

			// make ecal maps of jet/cluster rh collections
			if( DEBUG ) std::cout << " -- make ecal maps of jet/cluster rh collections" << std::endl;
			//auto nBcRhGroups = bcRhGroups.size();
			if( nGoodJetEvents < nEBEEMaps && nJets > 2 ){
				goodJetEvent = true;
				if( DEBUG ) std::cout << " --- fill SC ecal map" << std::endl;
				for (const auto recHit : jetScRhGroup ){
					const auto & idinfo = DetIDMap[getRawID(recHit)];
					if( idinfo.ecal == ECAL::EB ) ebeeMapSc[nGoodJetEvents]->Fill( idinfo.i2, idinfo.i1, (ijet+1)*10 );
				}//<<>>for (const auto recHit : jetScRhGroup )
            	for (const auto recHit : jetDrRhGroup ){
              		const auto & idinfo = DetIDMap[getRawID(recHit)];
              		if( idinfo.ecal == ECAL::EB ) ebeeMapDr[nGoodJetEvents]->Fill( idinfo.i2, idinfo.i1, (ijet+1)*10 );
            	}//<<>>for (const auto recHit : jetDrRhGroup )
				if( DEBUG ) std::cout << " --- fill BC/T/E ecal map" << std::endl;

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){
			//<<<<if( nGoodJetEvents < nEBEEMaps && nJets > 2 )

				for( uInt it = 0; it < bcRhGroups.size(); it++ ){
					for (const auto recHit : bcRhGroups[it] ){
						const auto & idinfo = DetIDMap[getRawID(recHit)];
						if( idinfo.ecal == ECAL::EB ){ 
							//std::cout << "  --- Filling ebeb map with : " << ijet << std::endl;
                    		ebeeMapBc[nGoodJetEvents]->Fill( idinfo.i2, idinfo.i1, (it+1)*10 );
                    		//std::cout << "  --- Filling ebeb map with : " << it+1 << std::endl;
                       		ebeeMapT[nGoodJetEvents]->Fill( idinfo.i2, idinfo.i1, recHit.time()+100 );
                        	//std::cout << "  --- Filling ebeb map with : " << recHit.time() << std::endl;
                        	ebeeMapE[nGoodJetEvents]->Fill( idinfo.i2, idinfo.i1, recHit.energy() );
                        	//std::cout << "  --- Filling ebeb map with : " << recHit.energy() << std::endl;
				}	}	}//<<>>for( uInt it = 0; it < bcRhGroups.size(); it++ )

            	if( DEBUG ) std::cout << " --- fill max values for ecal maps" << std::endl;
            	ebeeMapSc[nGoodJetEvents]->Fill( -85.0, 1.0, 100.0 ); 
            	ebeeMapBc[nGoodJetEvents]->Fill( -85.0, 1.0, 100.0 ); 
            	ebeeMapDr[nGoodJetEvents]->Fill( -85.0, 1.0, 100.0 ); 
            	ebeeMapE[nGoodJetEvents]->Fill( -85.0, 1.0, 1000.0 ); 
            	ebeeMapT[nGoodJetEvents]->Fill( -85.0, 1.0, 200.0 );

			}//<<>>if( nGoodJetEvents < nEBEEMaps ? && ... )

   	//<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){
				
			// Fill BC based jet times
			if( DEBUG ) std::cout << "Fill BC based jet times ------------------" << std::endl;
			auto nBCTimes = bcTimes.size();
			auto sumbcrhe = vfsum(bcRhGrpEnergy);
			auto bcemf = sumbcrhe/jet.energy();
			if( nBCTimes != 0 && bcemf > minEmf ){
				//if( rhCount >= minRHcnt ) 
				hist2d[51]->Fill(rhCount,bcRhCnt);
				hist1d[55]->Fill(bcRhGroups.size());
				auto jetBCTimeStats = getTimeDistStats( bcTimes, bcRhGrpEnergy );
            	//auto jetBCTimeStats = getTimeDistStats( bcTimes, bcEnergies );
				//auto jetBCTimeStats = getTimeDistStats( bcTimes );
				if( jetBCTimeStats[6] > -8.9 ) nGoodBcJets++;
				jetCBCMedTime.push_back(jetBCTimeStats[10]);//c med
	        	jetCBCMuTime.push_back(jetBCTimeStats[6]);//c mu
				hist1d[10]->Fill(jetBCTimeStats[10]);//c med
	        	hist1d[19]->Fill(jetBCTimeStats[6]);//c mu
				//std::cout << " - fill dbct hist " << std::endl;
				if( nBCTimes == 1 ){ hist1d[11]->Fill(-3.5); }
				else {//<<>>if( nBCTimes == 1 )
					for( uInt ita = 0; ita < nBCTimes; ita++ ){
                 		hist2d[50]->Fill(bcTimes[ita],bcRhGrpEnergy[ita]);
						//hist2d[50]->Fill(bcTimes[ita],bcEnergies[ita]);
						for( uInt itb = ita+1; itb < nBCTimes; itb++ ){
							auto dt = getdt(bcTimes[ita],bcTimes[itb]);
							hist1d[11]->Fill(dt);
							hist1d[23]->Fill(dt);
							auto effe = effMean(bcRhGrpEnergy[ita],bcRhGrpEnergy[itb]);
							hist2d[45]->Fill(dt, effe);
							hist2d[46]->Fill(dt, effe);
				}	}	}//<<>>if( nBCTimes == 1 ) : else	
				hist2d[54]->Fill( bcemf, jetemfrac );
          		hist2d[57]->Fill( sumbcrhe, jeteme );
           		hist2d[58]->Fill( sumbcrhe, sumscrhe );
            	hist2d[60]->Fill( scemf, bcemf );

        	} else { //<<>>if( nBCTimes == 0 )
               		jetCBCMedTime.push_back(-9.75);
               		jetCBCMuTime.push_back(-9.75);
               		hist1d[10]->Fill(-9.75);
			}//<<>>if( nBCTimes == 0 ) : else

    //<<<<for ( uInt ijet(0); ijet < nJets; ijet++ )
        //<<<<if( jetScRhGroup.size() >= minRHcnt ){

	   	} else { //<<>>if( jetSCGroup.size() >= minRHcnt)

			if( DEBUG ) std::cout << " - jetScRhGroup is < minRHcnt " << std::endl;
	        jetSCMuTime.push_back(-9.25);
	   	   	jetCSCMuTime.push_back(-9.25);
	   	   	jetSCMedTime.push_back(-9.25);
	        jetCSCMedTime.push_back(-9.25);
            jetPhMuTime.push_back(-9.25);
            jetEleMuTime.push_back(-9.25);
	    	jetCBCMuTime.push_back(-9.25);
	    	jetCBCMedTime.push_back(-9.25);
			//jetCSBCMuTime.push_back(-9.25);
				
	   	}//<<>>if( jetSCGroup.size() >= minRHcnt) : else			

		if( DEBUG ) std::cout << "Finished with get SC Med times --------------------------- " << std::endl;
			
		//****************************	photon/electron to kid pfcand -> SC matcher ********************************

		if( DEBUG ) std::cout << "Next Jet .......................... " << std::endl; 	
	}//<<>>for ( uInt ijet = 0; ijet < nJets; ijet++ )
	// ** end of jets	***************************************************************************************************

	hist1d[17]->Fill(jetHt);
   	hist1d[30]->Fill(nGoodDrJets);
   	hist1d[31]->Fill(nGoodScJets);
   	hist1d[32]->Fill(nGoodBcJets);
   	hist1d[33]->Fill(nUnJets);
   	if( nUnJets != 0 ) hist1d[34]->Fill(float(nJets)/nUnJets);
   	if( nGoodScJets != 0 ) hist1d[38]->Fill(float(nGoodBcJets)/nGoodScJets);
	if( nJets != 0 ){
   		hist1d[35]->Fill(float(nGoodDrJets)/nJets);
   		hist1d[36]->Fill(float(nGoodScJets)/nJets);
   		hist1d[37]->Fill(float(nGoodBcJets)/nJets);
	}//<<>>if( nUnJets != 0 )
	
	//-----------------------------------------------------------------------------------------------------
	// ***************************** d jetTime for back-to-back high pt jets	*****************************
	//auto dijetIdCut = 1;
	auto dijetPtMin = 200.0;
	auto difPtLmt = 0.8;
	auto htPctLmt = 0.8;
	auto dPhiLmt = 2.8;
	
/*
	if( jetSCTime.size() ){
		std::cout << "Jet times sizes : jetTime: " << jetTime.size() << " jetMedTime: "; 
		std::cout << jetMedTime.size() << " jetSCTime: " << jetSCTime.size() << std::endl;
	}//<<>>if( jetSCTime.size() )
*/

   	if( jetCSCMuTime.size() != jetCBCMuTime.size() ) hist1d[53]->Fill(-3.25);
	else if( jetCSCMuTime.size() != nJets ) hist1d[53]->Fill(-3.0);
	else for( uInt q = 0; q < nJets; q++ ) hist1d[53]->Fill( getdt(jetCSCMuTime[q],jetCBCMuTime[q]) ); 

	if( DEBUG ) std::cout << "Finding jet dt pairs" << std::endl;
	for ( uInt q = 0; q < nJets; q++ ){
		for ( uInt p = q+1; p < nJets; p++ ){
	
			if( DEBUG ) std::cout << " - filter jet pairs" << std::endl;
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

			if( DEBUG ) std::cout << " - get jet pair dt" << std::endl;
			auto dTmu = getdt( jetMuTime[q], jetMuTime[p] );
            auto dTmed = getdt( jetMedTime[q], jetMedTime[p] );
            auto dTcmu = getdt( jetCMuTime[q], jetCMuTime[p] );
            auto dTcmed = getdt( jetCMedTime[q], jetCMedTime[p] );
            if( DEBUG ) std::cout << "dT dR      : " << dTmu <<  " " << jetMuTime[q] << " " << jetMuTime[p] << std::endl;
	      	auto dTmusc = getdt( jetSCMuTime[q], jetSCMuTime[p] );
	      	auto dTmedsc = getdt( jetSCMedTime[q], jetSCMedTime[p] );
	      	auto dTcmusc = getdt( jetCSCMuTime[q], jetCSCMuTime[p] );
	      	auto dTcmedsc = getdt( jetCSCMedTime[q], jetCSCMedTime[p] );
            if( DEBUG ) std::cout << "dT SC      : " << dTmusc <<  " " << jetSCMuTime[q] << " " << jetSCMuTime[p] << std::endl;
	      	auto dTcmubc = getdt( jetCBCMuTime[q], jetCBCMuTime[p] );
            if( DEBUG ) std::cout << "dT cMu BC  : " << dTcmubc <<  " " << jetCBCMuTime[q] << " " << jetCBCMuTime[p] << std::endl;
	      	auto dTcmedbc = getdt( jetCBCMedTime[q], jetCBCMedTime[p] );
            if( DEBUG ) std::cout << "dT cMed BC : " << dTcmedbc <<  " " << jetCBCMedTime[q] << " " << jetCBCMedTime[p] << std::endl;
            auto dTmuph = getdt( jetPhMuTime[q], jetPhMuTime[p] );
            if( DEBUG ) std::cout << "dT Ph      : " << dTmuph <<  " " << jetPhMuTime[q] << " " << jetPhMuTime[p] << std::endl;
            auto dTmuel = getdt( jetEleMuTime[q], jetEleMuTime[p] );
            if( DEBUG ) std::cout << "dT Ele     : " << dTmuel <<  " " << jetEleMuTime[q] << " " << jetEleMuTime[p] << std::endl;

   	//<<<<for ( uInt q = 0; q < nJets; q++ ){
      	//<<<<for ( uInt p = q+1; p < nJets; p++ ){

			if( DEBUG ) std::cout << " - fill hists" << std::endl;

			auto dtThrs = -2.5;//removes default dt values from getdt
			if( dTmu > dtThrs ) hist1d[15]->Fill(dTmu);
	      	if( dTmed > dtThrs ) hist1d[16]->Fill(dTmed);
            if( dTcmu > dtThrs ) hist1d[39]->Fill(dTcmu);
            if( dTcmed > dtThrs ) hist1d[40]->Fill(dTcmed);

	      	if( dTmedsc > dtThrs ) hist1d[46]->Fill(dTmedsc);
	      	if( dTmusc > dtThrs ) hist1d[49]->Fill(dTmusc);
	      	if( dTcmusc > dtThrs ) hist1d[41]->Fill(dTcmusc);
	      	if( dTcmedsc > dtThrs ) hist1d[42]->Fill(dTcmedsc);

	      	if( dTcmubc > dtThrs ) hist1d[27]->Fill(dTcmubc);
	      	if( dTcmedbc > dtThrs ) hist1d[28]->Fill(dTcmedbc);

            if( dTmuph > dtThrs ) hist1d[59]->Fill(dTmuph);
            if( dTmuel > dtThrs ) hist1d[60]->Fill(dTmuel);

			hist2d[22]->Fill(dTmu,nJets);
	      	hist2d[26]->Fill(dTmu,diffPt);
	      	hist2d[27]->Fill(dTmu,htPct);
	      	hist2d[28]->Fill(dTmu,dPhi);
	     	hist2d[23]->Fill(dTmed,nJets);
	      	hist2d[29]->Fill(dTmu,diffPt);
	      	hist2d[30]->Fill(dTmu,htPct);
	      	hist2d[31]->Fill(dTmu,dPhi);

			if( DEBUG ) std::cout << " - fill dt vs eff e hists" << std::endl;
	      	auto effje = effMean(jetPHE[p],jetPHE[q]);
	      	hist2d[43]->Fill(dTmusc,effje);
	      	hist2d[44]->Fill(dTmu,effje);

		}//<<>>for ( uInt p = q+1; p < nJets; p++ )
	}//<<>>for ( uInt q = 0; q < nJets; q++ )
	//-------------------------------------------------------------------------------

	if( goodJetEvent ) nGoodJetEvents++;

	// -- Fill output trees ------------------------------------------
	if( DEBUG ) std::cout << "---------- Next Event -----" << std::endl;
	outTree->Fill();

	// -- EOFun ------------------------------------------------------
	//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//	 ESHandle<SetupData> pSetup;
	//	 iSetup.get<SetupRecord>().get(pSetup);
	//#endif
}//>>>>void LLPgammaAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)


// ------------ method called once each job just before starting event loop	------------
void LLPgammaAnalyzer::beginJob(){

	// Global Varibles
	nGoodJetEvents = 0;

   	// Set up DetIdMap
   	DetIDMap = SetupDetIDs();

	// Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	// Book histograms
	
	int jtdiv(400);
	float jtran(8);
	int jdtdiv(200);
	float jdtran(4);
   	int jztdiv(100);
   	float jztran(2);
	int rhcnt(80);

   	//------ 1D Hists --------------------------------------------------------------------------

	hist1d[0] = fs->make<TH1D>("jetRHTime", "jetRHTime", 2000, -100, 100);
	hist1d[29] = fs->make<TH1D>("jetMuTime", "jetMuTime", jtdiv, -1*jtran, jtran);
	hist1d[1] = fs->make<TH1D>("jetRHMulti", "jetRHMulti", rhcnt, 0, rhcnt);
	hist1d[2] = fs->make<TH1D>("jetTimeError", "jetTimeError", 300, 0, 3);
	hist1d[3] = fs->make<TH1D>("jetTimeRMS", "jetTimeRMS", 200, 0, 20);
	hist1d[4] = fs->make<TH1D>("jetMedTime", "jetMedTime", jtdiv, -1*jtran, jtran);
	hist1d[5] = fs->make<TH1D>("jetRawTime", "jetRawTime", jtdiv, -1*jtran, jtran);
	hist1d[6] = fs->make<TH1D>("jetCMuTime", "jetCMuTime", jtdiv, -1*jtran, jtran);
	hist1d[7] = fs->make<TH1D>("jetCMedTime", "jetCMedTime", jtdiv, -1*jtran, jtran);
	hist1d[8] = fs->make<TH1D>("jetCSCMuTime", "jetCSCMuTime", jtdiv, -1*jtran, jtran);
	hist1d[9] = fs->make<TH1D>("jetCSCMedTime", "jetCSCMedTime", jtdiv, -1*jtran, jtran);
	hist1d[10] = fs->make<TH1D>("jetCBCMedTime", "jetCBCMedTime", jtdiv, -1*jtran, jtran);
	hist1d[11] = fs->make<TH1D>("jetBCTimeDiff", "jetBCTimeDiff", jdtdiv, -1*jdtran, jdtran);

	hist1d[12] = fs->make<TH1D>("jetPt", "jetPt", 500, 0, 500);
	hist1d[13] = fs->make<TH1D>("jetPhi", "jetPhi", 700, -3.5, 3.5);
	hist1d[14] = fs->make<TH1D>("jetEta", "jetEta", 700, -3.5, 3.5);
	hist1d[15] = fs->make<TH1D>("jetdtmu", "jetdtmu", jdtdiv, -1*jdtran, jdtran);
	hist1d[16] = fs->make<TH1D>("jetdtmed", "jetdtmed", jdtdiv, -1*jdtran, jdtran);
	hist1d[17] = fs->make<TH1D>("jetHt", "jetHt", 1000, 0, 1000);
	hist1d[18] = fs->make<TH1D>("nJet", "nJets", 21, -0.5, 20.5);
	hist1d[19] = fs->make<TH1D>("jetCBCMuTime", "jetCBCMuTime", jtdiv, -1*jtran, jtran);
	hist1d[20] = fs->make<TH1D>("nBCDupRHs", "nBCDupRHs", 51, -0.5, 50.5);
	hist1d[21] = fs->make<TH1D>("nOrgBCRecHits", "nOrgBCRecHits", 51, -0.5, 50.5);
	hist1d[22] = fs->make<TH1D>("nRedBCRecHits", "nRedBCRecHits", 51, -0.5, 50.5);
   	hist1d[23] = fs->make<TH1D>("jetBCTimeDiffZoom", "jetBCTimeDiffZoom", jztdiv, -1*jztran, jztran);

	hist1d[24] = fs->make<TH1D>("diffPt", "diffPt", 100, 0, 1);
	hist1d[25] = fs->make<TH1D>("htPct", "htPct", 100, 0, 1);
	hist1d[26] = fs->make<TH1D>("dPhi", "dPhi", 32, 0, 3.2);

	hist1d[27] = fs->make<TH1D>("jetcmudtbc", "jetcmudtbc", jdtdiv, -1*jdtran, jdtran);
	hist1d[28] = fs->make<TH1D>("jetcmeddtbc", "jetcmeddtbc", jdtdiv, -1*jdtran, jdtran);
	// ----  moved to 2nd from top of list : hist1d[29] ( after hist1d[0] )
	hist1d[30] = fs->make<TH1D>("nGoodDrJets", "nGoodDrJets", 21, -0.5, 20.5);
	hist1d[31] = fs->make<TH1D>("nGoodScJets", "nGoodScJets", 21, -0.5, 20.5);
	hist1d[32] = fs->make<TH1D>("nGoodBcJets", "nGoodBcJets", 21, -0.5, 20.5);
	hist1d[33] = fs->make<TH1D>("nUnJets", "nUnJets", 51, -0.5, 50.5);
	hist1d[34] = fs->make<TH1D>("pJets", "pJets", 110, 0, 1.1);
	hist1d[35] = fs->make<TH1D>("pGoodDrJets", "pGoodDrJets", 110, 0, 1.1);
	hist1d[36] = fs->make<TH1D>("pGoodScJets", "pGoodScJets", 110, 0, 1.1);
	hist1d[37] = fs->make<TH1D>("pGoodBcJets", "pGoodBcJets", 110, 0, 1.1);
	hist1d[38] = fs->make<TH1D>("pGoodBcToScJets", "pGoodBcToScJets", 110, 0, 1.1);

	hist1d[39] = fs->make<TH1D>("jetcmudt", "jetcmudt", jdtdiv, -1*jdtran, jdtran);
	hist1d[40] = fs->make<TH1D>("jetcmeddt", "jetcmeddt", jdtdiv, -1*jdtran, jdtran);
	hist1d[41] = fs->make<TH1D>("jetcmudtsc", "jetcmudtsc", jdtdiv, -1*jdtran, jdtran);
	hist1d[42] = fs->make<TH1D>("jetcmeddtsc", "jetcmeddtsc", jdtdiv, -1*jdtran, jdtran);

	hist1d[43] = fs->make<TH1D>("nPhotonsPerJet","nPhotonsPerJet", 21, -0.5, 20.5);

	hist1d[44] = fs->make<TH1D>("jetSCmedTime", "jetSCmedTime", jtdiv, -1*jtran, jtran);
	hist1d[45] = fs->make<TH1D>("jetSCmuTime", "jetSCmuTime", jtdiv, -1*jtran, jtran);
	hist1d[46] = fs->make<TH1D>("jetSCTimeRms", "jetSCTimeRms", 200, 0, 20);
	hist1d[47] = fs->make<TH1D>("jetSCrawTime", "jetSCrawTime", jtdiv, -1*jtran, jtran);

	hist1d[48] = fs->make<TH1D>("jetmeddtsc", "jetmeddtsc", jdtdiv, -1*jdtran, jdtran);
	hist1d[49] = fs->make<TH1D>("jetmudtsc", "jetmudtsc", jdtdiv, -1*jdtran, jdtran);

	hist1d[50] = fs->make<TH1D>("jetSCTimeSkew", "jetSCTimeSkew", 40, -2.0, 2.0);
	hist1d[51] = fs->make<TH1D>("jetPhotons", "jetPhotons", 21, -0.5, 20.5);
	hist1d[52] = fs->make<TH1D>("jetElectrons", "jetElectrons", 21, -0.5, 20.5);

   	hist1d[53] = fs->make<TH1D>("scbcdt", "scbcdt", jdtdiv, -1*jdtran, jdtran);
   	hist1d[54] = fs->make<TH1D>("bc1rhef", "bc1rhef", 110, 0, 1.1);
   	hist1d[55] = fs->make<TH1D>("nBCinJet", "nBCinJet", 11, -0.5, 10.5);
   	hist1d[56] = fs->make<TH1D>("bcMrhef", "bcMrhef", 110, 0, 1.1);

    hist1d[57] = fs->make<TH1D>("jetCPhMuTime", "jetCPhMuTime", jtdiv, -1*jtran, jtran);
    hist1d[58] = fs->make<TH1D>("jetCPhMedTime", "jetCPhMedTime", jtdiv, -1*jtran, jtran);
    hist1d[59] = fs->make<TH1D>("jetmudtph", "jetmudtph", jdtdiv, -1*jdtran, jdtran);
    hist1d[60] = fs->make<TH1D>("jetmudtel", "jetmudtel", jdtdiv, -1*jdtran, jdtran);
    hist1d[61] = fs->make<TH1D>("jetCEleMuTime", "jetCEleMuTime", jtdiv, -1*jtran, jtran);
    hist1d[62] = fs->make<TH1D>("jetCEleMedTime", "jetCEleMedTime", jtdiv, -1*jtran, jtran);

    hist1d[63] = fs->make<TH1D>("scEtaPhiAngle3D", "scEtaPhiAngle3D", 660, -0.2, 6.4);
    hist1d[64] = fs->make<TH1D>("scEtaTimAngle3D", "scEtaTimAngle3D", 660, -0.2, 6.4);
    hist1d[65] = fs->make<TH1D>("scEtaPhiAngle2D", "scEtaTimAngle2D", 660, -0.2, 6.4);

    hist1d[66] = fs->make<TH1D>("sciEta3D", "sciEta", 171, -85, 85);
    hist1d[67] = fs->make<TH1D>("sciPhi3D", "sciPhi", 361, 0, 360);
    hist1d[68] = fs->make<TH1D>("scTime3D", "scTime", 5000, -25, 25);
    hist1d[69] = fs->make<TH1D>("sciEta2D", "sciEtaDiff", 201, -100, 100 );
    hist1d[70] = fs->make<TH1D>("sciPhi2D", "sciPhiDiff", 201, -100, 100);
    hist1d[71] = fs->make<TH1D>("sciTim2D", "sciTimDiff", 400, -10, 10);
    hist1d[72] = fs->make<TH1D>("sciAngle2D", "sciAngleSph", 660, -0.2, 6.4);
    hist1d[73] = fs->make<TH1D>("scRotAngle", "scRotAngle", 660, -0.2, 6.4);
    hist1d[74] = fs->make<TH1D>("scAngleTest", "scAngleTest", 660, -0.2, 6.4);
    hist1d[75] = fs->make<TH1D>("sciEta3Diff", "sciEta3Diff", 400, -10, 10 );
    hist1d[76] = fs->make<TH1D>("sciPhi3Diff", "sciPhi3Diff", 400, -10, 10);
    hist1d[77] = fs->make<TH1D>("sciTim3Diff", "sciTim3Diff", 400, -10, 10);
    hist1d[78] = fs->make<TH1D>("scSphEgn0", "scSphEgn0", 400, -10, 10);
    hist1d[79] = fs->make<TH1D>("scSphEgn1", "scSphEgn1", 400, -10, 10);
    hist1d[80] = fs->make<TH1D>("scSinTest", "scSinTest", 200, -1, 1);
    hist1d[81] = fs->make<TH1D>("eginValueSph", "eginValueSph", 1500, -2.25, 12.75);
    hist1d[82] = fs->make<TH1D>("dIPhiTest", "dIPhiTest", 1440, -720, 720);
    hist1d[83] = fs->make<TH1D>("meanIPhiTest", "meanIPhiTest", 1440, -720, 720);
    hist1d[84] = fs->make<TH1D>("meanIEtaTest", "meanIEtaTest", 200, -100, 100);
    hist1d[85] = fs->make<TH1D>("meanTimeTest", "meanTimeTest", 2000, -25, 25);
    hist1d[86] = fs->make<TH1D>("eginValue3D", "eginValue3D", 1500, -2.25, 12.75);
    hist1d[87] = fs->make<TH1D>("eginValueSphCut", "eginValueSphCut", 1500, -2.25, 12.75);

    hist1d[88] = fs->make<TH1D>("jetGenTime", "jetGenTime", jdtdiv, -1*jdtran, jdtran);
    hist1d[89] = fs->make<TH1D>("jetmudtgen", "jetmudtgen", jdtdiv, -1*jdtran, jdtran);
    hist1d[90] = fs->make<TH1D>("jetImpactAngle", "jetImpactAngle", 660, -0.2, 6.4);

    // used below : hist1d[87] = fs->make<TH1D>("cluster_etprofile", "Cluster Eta Time Profile Sph", 80, -40, 40);

	//------ 2D Hists --------------------------------------------------------------------------

   //hist2d[0] = fs->make<TH2D>("jt_pt", "jt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[1] = fs->make<TH2D>("jt_pt", "jt_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[2] = fs->make<TH2D>("jt_id", "jt_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
	hist2d[3] = fs->make<TH2D>("jt_nhf", "jt_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[4] = fs->make<TH2D>("jt_chf", "jt_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[5] = fs->make<TH2D>("jt_nemf", "jt_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[6] = fs->make<TH2D>("jt_cemf", "jt_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[7] = fs->make<TH2D>("jt_muf", "jt_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[8] = fs->make<TH2D>("jt_nhm", "jt_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
	hist2d[9] = fs->make<TH2D>("jt_chm", "jt_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

	hist2d[10] = fs->make<TH2D>("jt_medt", "jt_medt", jtdiv, -1*jtran, jtran, 200, -10, 10);
	hist2d[11] = fs->make<TH2D>("jt_rms", "jt_rms", jtdiv, -1*jtran, jtran, 200, 0, 20);
	hist2d[12] = fs->make<TH2D>("jt_err", "jt_err", jtdiv, -1*jtran, jtran, 300, 0, 3);

	hist2d[13] = fs->make<TH2D>("medt_pt", "medt_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
	hist2d[14] = fs->make<TH2D>("medt_id", "medt_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
	hist2d[15] = fs->make<TH2D>("medt_nhf", "medt_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[16] = fs->make<TH2D>("medt_chf", "medt_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[17] = fs->make<TH2D>("medt_nemf", "medt_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[18] = fs->make<TH2D>("medt_cemf", "medt_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[19] = fs->make<TH2D>("medt_muf", "medt_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
	hist2d[20] = fs->make<TH2D>("medt_nhm", "medt_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
	hist2d[21] = fs->make<TH2D>("medt_chm", "medt_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

	hist2d[22] = fs->make<TH2D>("jdtmu_nJets", "jdtmu_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);
	hist2d[23] = fs->make<TH2D>("jdtmed_nJets", "jdtmed_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);

	hist2d[24] = fs->make<TH2D>("jt_nrh", "jt_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);
	hist2d[25] = fs->make<TH2D>("medt_nrh", "medt_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);

	hist2d[26] = fs->make<TH2D>("jdtmu_diffPt", "jdtmu_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[27] = fs->make<TH2D>("jdtmu_htPct", "jdtmu_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[28] = fs->make<TH2D>("jdtmu_dPhi", "jdtmu_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);
	hist2d[29] = fs->make<TH2D>("jdtmed_diffPt", "jdtmed_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[30] = fs->make<TH2D>("jdtmed_htPct", "jdtmed_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
	hist2d[31] = fs->make<TH2D>("jdtmed_dPhi", "jdtmed_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);

	hist2d[32] = fs->make<TH2D>("jt_sceta", "jt_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[33] = fs->make<TH2D>("jt_scphi","jt_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[34] = fs->make<TH2D>("jt_scenr", "jt_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

	hist2d[35] = fs->make<TH2D>("medt_sceta", "medt_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[36] = fs->make<TH2D>("medt_scphi","medt_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
	hist2d[37] = fs->make<TH2D>("medt_scenr", "medt_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

	hist2d[38] = fs->make<TH2D>("rht_rhe", "rht_rhe", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

	hist2d[39] = fs->make<TH2D>("njrh_nscrh", "njrh_nscrh", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);
	hist2d[40] = fs->make<TH2D>("je_sce", "je_sce", 500, 0, 500, 500, 0, 500);
	hist2d[41] = fs->make<TH2D>("je_ege", "je_ege", 500, 0, 500, 500, 0, 500);
	hist2d[42] = fs->make<TH2D>("ege_sce", "ege_sce", 500, 0, 500, 500, 0, 500);
	hist2d[43] = fs->make<TH2D>("scdt_effje", "scdt_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);
	hist2d[44] = fs->make<TH2D>("jdt_effje", "jdt_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);

   	hist2d[45] = fs->make<TH2D>("bcdt_effbce", "bcdt_effbce", jdtdiv, -1*jdtran, jdtran, 250, 0, 250);
   	hist2d[46] = fs->make<TH2D>("bcdt_effbce_zoom", "bcdt_effbce_zoom", jztdiv, -1*jztran, jdtran, 250, 0, 250);
   	hist2d[47] = fs->make<TH2D>("nbcrh_srhe", "nbcrh_srhe", rhcnt, 0, rhcnt, 250, 0, 250);
   	hist2d[48] = fs->make<TH2D>("nbcrh_srhe_other", "nrh_srhe_other", rhcnt, 0, rhcnt, 250, 0, 250);
   	hist2d[49] = fs->make<TH2D>("bcrht_bcrhe", "bcrht_bcrhe", jtdiv, -1*jtran, jtran, 250, 0, 250);
   	hist2d[50] = fs->make<TH2D>("bct_bce", "bct_bce", jtdiv, -1*jtran, jtran, 250, 0, 250);
   	hist2d[51] = fs->make<TH2D>("njrh_nsbcrh", "njrh_nsbcrh", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);

   	hist2d[52] = fs->make<TH2D>("dremf_emf", "dremf_emf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[53] = fs->make<TH2D>("scemf_emf", "scemf_emf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[54] = fs->make<TH2D>("bcemf_emf", "bcemf_emf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[55] = fs->make<TH2D>("dreme_eme", "dreme_eme", 500, 0, 500, 500, 0, 500);
   	hist2d[56] = fs->make<TH2D>("sceme_eme", "sceme_eme", 500, 0, 500, 500, 0, 500);
   	hist2d[57] = fs->make<TH2D>("bceme_eme", "bceme_eme", 500, 0, 500, 500, 0, 500);
   	hist2d[58] = fs->make<TH2D>("sce_bce", "sce_bce", 500, 0, 500, 500, 0, 500);
   	hist2d[59] = fs->make<TH2D>("dremf_scemf", "dremf_scemf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[60] = fs->make<TH2D>("scemf_bcemf", "scemf_bcemf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[61] = fs->make<TH2D>("epaf_epf", "epaf_epf", 110, 0, 1.1, 110, 0, 1.1);
   	hist2d[62] = fs->make<TH2D>("epf_emf", "epf_emf", 110, 0, 1.1, 110, 0, 1.1);

	auto stddiv = 120;
    auto stdtran = 3;
    hist2d[63] = fs->make<TH2D>("jetEta_stdSCt", "jetEta_stdSCt", 200, -2.5, 2.5, stddiv, 0, stdtran);
    hist2d[64] = fs->make<TH2D>("jetEtaetaMmt_stdSCt", "jetEtaetaMmt_stdSCt", 800, 0, 0.02, stddiv, 0, stdtran);
    hist2d[65] = fs->make<TH2D>("jetPhiphiMnt_stdSCt", "jetPhiphiMnt_stdSCt", 800, 0, 0.02, stddiv, 0, stdtran);
    hist2d[66] = fs->make<TH2D>("jetEtaphiMnt_stdSCt", "jetEtaphiMnt_stdSCt", 1600, -0.02, 0.02, stddiv, 0, stdtran);
    hist2d[67] = fs->make<TH2D>("jetMaxD_stdSCt", "jetMaxD_stdSCt", 400, 0, 1, stddiv, 0, stdtran);
    hist2d[68] = fs->make<TH2D>("jetConPtDis_stdSCt", "jetConPtDis_stdSCt", 400, 0, 1, stddiv, 0, stdtran);
    hist2d[69] = fs->make<TH2D>("jetConEtaPhiSprd_stdSCt", "jetConEtaPhiSprd_stdSCt", 600, -0.005, 0.01, stddiv, 0, stdtran);
    hist2d[70] = fs->make<TH2D>("jetArea_stdSCt", "jetArea_stdSCt", 200, 0, 1, stddiv, 0, stdtran);
    hist2d[71] = fs->make<TH2D>("jetNCarry_stdSCt", "jetNCarry_stdSCt", 40, 0, 40, stddiv, 0, stdtran);
    hist2d[72] = fs->make<TH2D>("jetNConst_stdSCt", "jetNConst_stdSCt", 100, 0, 100, stddiv, 0, stdtran);

	auto cldiv = 80;
	auto cltrn = 40;
	hist2d[73] = fs->make<TH2D>("cluster_tmap_rot", "Cluster Time Map Rotated", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[74] = fs->make<TH2D>("cluster_occmap_rot", "Cluster Occ Map Rotated", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[75] = fs->make<TH2D>("cluster_tmap", "Cluster Time Map", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[76] = fs->make<TH2D>("cluster_occmap", "Cluster Occ Map", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);

    hist2d[88] = fs->make<TH2D>("scEtaPhiAngle2D_eginValueSph","scEtaTimAngle(x) v eginValueSph(y)",660, -0.2, 6.4, 160, 0.4, 1.2);

    auto cl3ddiv = 200;
    auto cl3dtrn = 4;
    auto cl3ddiv1 = 200;
    auto cl3dtrn1 = 4;
    hist2d[77] = fs->make<TH2D>("cluster_3D_tmap", "Cluster Time Map 3D", cl3ddiv, -1*cl3dtrn, cl3dtrn, cl3ddiv, -1*cl3dtrn, cl3dtrn);
    hist2d[78] = fs->make<TH2D>("cluster_3D_occmap", "Cluster Occ Map 3D", cl3ddiv, -1*cl3dtrn, cl3dtrn, cl3ddiv, -1*cl3dtrn, cl3dtrn);
    hist2d[79] = fs->make<TH2D>("cluster_3D_etmap", "Cluster Eta(x)Time(y) Map 3D", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1, cl3ddiv, -1*cl3dtrn, cl3dtrn);
    hist2d[80] = fs->make<TH2D>("cluster_3D_etprofl", "Clstr Eta(x)Time(y) Map 3DProfile", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1, cl3ddiv, -1*cl3dtrn, cl3dtrn);

    hist2d[81] = fs->make<TH2D>("scSphEgn01", "scSphEgn01", 200, -5, 5, 200, -5, 5);
    hist2d[82] = fs->make<TH2D>("sc3DEgn01", "sc3DEgn01", 200, -5, 5, 200, -5, 5);
    hist2d[83] = fs->make<TH2D>("sc3DEgn02", "sc3DEgn02", 200, -5, 5, 200, -5, 5);
    hist2d[84] = fs->make<TH2D>("cluster_3D_etwtmap", "Cluster Eta(x)Time(y) Wt Map 3D", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1, cl3ddiv, -1*cl3dtrn, cl3dtrn);
    hist2d[85] = fs->make<TH2D>("cluster_3D_etwtprofl", "Clstr Eta(x)Time(y) WtMap 3DProfile", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1, cl3ddiv, -1*cl3dtrn, cl3dtrn);

    auto clsphdiv = 400;
    auto clsphtrn = 4;
    auto cwdiv = 80;
    auto cwtrn = 40;
    hist2d[86] = fs->make<TH2D>("cluster_ptwtmap", "Cluster Phi(x)Time(y) Wt Map Sph", cwdiv, -1*cwtrn, cwtrn, clsphdiv, -1*clsphtrn, clsphtrn);
    hist2d[87] = fs->make<TH2D>("cluster_etwtmap", "Cluster Eta(x)Time(y) Wt Map Sph", cwdiv, -1*cwtrn, cwtrn, clsphdiv, -1*clsphtrn, clsphtrn);
	hist1d[87] = fs->make<TH1D>("cluster_etprofile", "Cluster Eta Time Profile Sph", cwdiv, -1*cwtrn, cwtrn);
    hist1d[88] = fs->make<TH1D>("cluster_et3Dprofl", "Cluster Eta Time Profile 3D", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1);


   //------ ECAL Map Hists --------------------------------------------------------------------------

	for(int it=0; it<nEBEEMaps; it++){
		string stsc("ebeeMapSc_"+std::to_string(it)); 
		ebeeMapSc[it] = fs->make<TH2D>( stsc.c_str(), stsc.c_str(), 361, -90, 90, 721, 0, 360); 
      	string stbc("ebeeMapBc_"+std::to_string(it));
      	ebeeMapBc[it] = fs->make<TH2D>( stbc.c_str(), stbc.c_str(), 361, -90, 90, 721, 0, 360);
      	string stdr("ebeeMapDr_"+std::to_string(it));
      	ebeeMapDr[it] = fs->make<TH2D>( stdr.c_str(), stdr.c_str(), 361, -90, 90, 721, 0, 360);
      	string stt("ebeeMapTime_"+std::to_string(it));
      	ebeeMapT[it] = fs->make<TH2D>( stt.c_str(), stt.c_str(), 361, -90, 90, 721, 0, 360);
      	string ste("ebeeMapEnr_"+std::to_string(it));
      	ebeeMapE[it] = fs->make<TH2D>( ste.c_str(), ste.c_str(), 361, -90, 90, 721, 0, 360);
	}//<<>>for(int it=0; it<nEBEEMaps; it++)

	std::cout << "Histograms Booked" << std::endl;

	// Create output Tree branches -----------------------------

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

}//>>>>void LLPgammaAnalyzer::beginJob()


// ------------ method called once each job just after ending the event loop	------------
void LLPgammaAnalyzer::endJob(){

    normTH1D(hist1d[15]);
    normTH1D(hist1d[16]);
    normTH1D(hist1d[39]);
    normTH1D(hist1d[40]);

	normTH1D(hist1d[46]);
    normTH1D(hist1d[49]);
    normTH1D(hist1d[41]);
    normTH1D(hist1d[42]);

    normTH1D(hist1d[27]);
    normTH1D(hist1d[28]);

    normTH1D(hist1d[59]);
    normTH1D(hist1d[60]);

	//hist2d[73]->Divide(hist2d[74]);
	thresDivTH2D( hist2d[73], hist2d[74], 100 );
    //hist2d[75]->Divide(hist2d[76]);
    thresDivTH2D( hist2d[75], hist2d[76], 100 );
    //hist2d[77]->Divide(hist2d[78]);
    thresDivTH2D( hist2d[77], hist2d[78], 100 );
    thresDivTH2D( hist2d[79], hist2d[84], 100 );
    thresDivTH2D( hist2d[80], hist2d[85], 100 );
	profileTH2D( hist2d[87], hist1d[87]);
    profileTH2D( hist2d[85], hist1d[88]);

}//>>>>void LLPgammaAnalyzer::endJob()


// ------------ method fills 'descriptions' with the allowed parameters for the module	------------
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
}//>>>>void LLPgammaAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)


//define this as a plug-in
DEFINE_FWK_MODULE(LLPgammaAnalyzer);
