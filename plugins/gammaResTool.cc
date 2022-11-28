// -*- C++ -*-
//
// Package:    GammaResTool
// Class:      GammaResTool
//
/**\class GammaResTool GammaResTool.cc LLPgammaAnalyzer/plugins/GammaResTool.cc

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

#include "LLPGamma/LLPgammaAnalyzer/plugins/GammaResTool.hh"
using namespace std;

//#define DEBUG true
#define DEBUG false

//
// constructors and destructor
//
GammaResTool::GammaResTool(const edm::ParameterSet& iConfig) :

// -- declare tags ----------------------------------------------------------

	// flags
	hasGenInfo (iConfig.existsAs<bool>("hasGenInfo")  ? iConfig.getParameter<bool>("hasGenInfo")  : false),

	// tracks
	tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),

    // vertices
    verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),
	
	// electrons
	electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),  

	// recHits
	recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),  
	recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

	// gedphotons
	gedPhotonsTag(iConfig.getParameter<edm::InputTag>("gedPhotons")),

	// ootPhotons
	ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons"))

// -- end of tag declarations ---------------------------------------
{ //<<<< GammaResTool::GammaResTool(const edm::ParameterSet& iConfig) :

	usesResource();
	usesResource("TFileService");

// -- consume tags ------------------------------------------------------------
	if( DEBUG ) std::cout << "In constructor for GammaResTool - tag and tokens" << std::endl;

	// tracks 
	tracksToken_				= consumes<std::vector<reco::Track>>(tracksTag);

	// vertices
	verticesToken_				= consumes<std::vector<reco::Vertex>>(verticesTag);

	// leptons
	electronsToken_				= consumes<std::vector<pat::Electron>>(electronsTag);

	// rechits
	recHitsEBToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
	recHitsEEToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

	// photons
	gedPhotonsToken_ 			= consumes<std::vector<pat::Photon>>(gedPhotonsTag);
	ootPhotonsToken_ 			= consumes<std::vector<pat::Photon>>(ootPhotonsTag);

// ---------------------------------------------------------------------------------
}//>>>>GammaResTool::GammaResTool(const edm::ParameterSet& iConfig)


GammaResTool::~GammaResTool(){
	///////////////////////////////////////////////////////////////////
	// do anything here that needs to be done at desctruction time   //
	// (e.g. close files, deallocate resources etc.)                 //
	///////////////////////////////////////////////////////////////////
}//>>>>GammaResTool::~GammaResTool()


//
// member functions
//

detIdMap GammaResTool::SetupDetIDs(){

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

}//>>>>detIdMap GammaResTool::SetupDetIDs()

/*
float GammaResTool::getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

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

}//>>>>float  GammaResTool::getSeedTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )
*/

float GammaResTool::getPhotonSeedTime( pat::Photon ){

	const auto & phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
	const auto & seedDetId = phosc->seed()->seed(); // get seed detid
	const auto recHits = ((seedDetId.subdetId() == EcalSubdetector::EcalBarrel) ? recHitsEB_ : recHitsEE_); // which recHits to use
	const auto seedHit = recHits->find(seedDetId); // get the underlying rechit
	const auto seedTime = ((seedHit != recHits->end()) ? seedHit->time() : -9999.f);
	return seedTime;

}//<<>>float GammaResTool::getPhotonSeedTime( pat::Photon )

int GammaResTool::getRhIdx( uInt rhDetID ){

    //b_rhID->GetEntry(entry);
    for( int idx = 0; idx < rhID->size(); idx++ ){ if( rhDetID == (*rhID)[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int makehists::getRhIdx( int rhDetID )


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// ------------ method called for each event	------------
void GammaResTool::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

	using namespace edm;

// -- Consume Tokens --------------------------------------------
	if( DEBUG ) std::cout << "Consume Tokens -------------------------------------------- " << std::endl;

	// TRACKS
	iEvent.getByToken(tracksToken_, tracks_);

	// VERTICES
	iEvent.getByToken(verticesToken_, vertices_);

	// LEPTONS & PHOTONS
	iEvent.getByToken(electronsToken_, electrons_);

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

    std::vector<EcalRecHit>         frechits;
    std::vector<pat::Photon>        fphotons;
	std::vector<bool>				phoOOT;
    std::vector<pat::Electron>  	felectrons;
	
	if( DEBUG ) std::cout << "Processing RecHits" << std::endl;
	for (const auto recHit : *recHitsEB_ ){ if( recHit.energy > 1.0 ) frechits.push_back(recHit); }
    for (const auto recHit : *recHitsEE_ ){ if( recHit.energy > 1.0 ) frechits.push_back(recHit); }

    if( DEBUG ) std::cout << "Processing gedPhotons" << std::endl;
    for( const auto photon : *gedPhotons_ ){
		if( getPhotonSeedTime(photon) > -25.0 ){ fphotons.push_back(photon); phoOOT.push_back(false); }
    }//<<>>for( const auto photon : *gedPhotons_ )

    if( DEBUG ) std::cout << "Processing ootPhotons" << std::endl;
    for( const auto ootphoton : *ootPhotons_ ){
        if( getPhotonSeedTime(ootphoton) > -25.0 ){ fphotons.push_back(photon); phoOOT.push_back(true); }
    }//<<>>for( const auto photon : *gedPhotons_ )

	const auto nPhotons = fphotons.size();
    vector<bool> phoExcluded(nPhotons+1,false);
    for( int io = 0; io < nPhotons; io++ ){
        double minDr(0.1);
		double lowestDr(10.0);
        int match(-1);
		auto ofPhoton = (*fphotons)[io];
        auto oEta = ofPhoton.eta();
        auto oPhi = ofPhoton.phi();
        for( int ip = io; ip < nPhotons; ip++ ){
			auto pfPhoton = (*fphotons)[ip];
            auto pEta = pfPhoton.eta();
            auto pPhi = pfPhoton.phi();
            auto dRmatch = deltaR( pEta, oEta, pPhi, oPhi );
            if( dRmatch < lowestDr ){ lowestDr = dRmatch; match = ip; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( lowestDr < minDr ){
            auto pPt = ((*fphotons)[match]).pt();
            auto oPt = ((*fphotons)[io]).pt();
            if( oPt > pPt ) phoExcluded[match] = true;
            else phoExcluded[io] = true;
        }//<<>>if( dRmatch < 0.3 )
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )

    if( DEBUG ) std::cout << "Processing Electrons" << std::endl;
	for( const auto electron : *electrons_ ){
		felectrons.push_back(electron);
	}//<<>>for( const auto electron : *electrons_ )

    //------------------------------------------------------------------------------------
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

        if( DEBUG ) std::cout << " -- proccesing ID info" << std::endl;
        // something in this section is seg faluting after several rechits for crab jobs
        const auto recHitID = getRawID(recHit);
        auto isEB = getIsEB(recHit); // which subdet
        const auto & idinfo = DetIDMap[recHitID];
        if( DEBUG ) std::cout << " -- proccesing EBEE info" << std::endl;
        if( DEBUG ) std::cout << " -- proccesing GEO info" << std::endl;
        //const auto geometry( isEB ? barrelGeometry : endcapGeometry );
        const auto geometry( ( idinfo.ecal == ECAL::EB ) ? barrelGeometry : endcapGeometry );
        auto recHitPos = geometry->getGeometry(recHit.detid())->getPosition();
        if( DEBUG ) std::cout << " -- proccesing POSITION info" << std::endl;
        const auto rhX = recHitPos.x();
        const auto rhY = recHitPos.y();
        const auto rhZ = recHitPos.z();
		const auto rhTrigTow = idinfo.TT;
        if( DEBUG ) std::cout << " -- proccesing TOF info" << std::endl;
        const auto d_rh = hypo(rhX,rhY,rhZ);
        const auto d_pv = hypo(rhX-vtxX,rhY-vtxY,rhZ-vtxZ);
        const auto tof = (d_rh-d_pv)/SOL;
        if( DEBUG ) std::cout << " -- proccesing SWISSCROSS info" << std::endl;
        float swisscross(0.0);
        if( isEB ) swisscross = EcalTools::swissCross(recHitID, *recHitsEB_, 0.0, true);
        else swisscross = EcalTools::swissCross(recHitID, *recHitsEE_, 0.0, true);

        if( DEBUG ) std::cout << " -- proccesing LASER info" << std::endl;
        // adcToGeVInfo : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#0204
        const auto laser = laserH->getLaserCorrection(recHitID,evTime);
        const auto interCalibIter = interCalibMap->find(recHitID);
        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        if( DEBUG ) std::cout << " -- proccesing ADC info" << std::endl;
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) rhadcToGeV[pos] = (laser*interCalib*adcToGeV);
        const float adcToGeV( isEB ? adcToGeVEB : adcToGeVEE );
        if( DEBUG ) std::cout << " -- proccesing PED info" << std::endl;
        // pedestal info
        const auto & pediter = pedestalsH->find(recHitID);

        if( DEBUG ) std::cout << " -- storing values BASE" << std::endl;
        rhID.push_back(recHitID);
        rhPosX.push_back(rhX);
        rhPosY.push_back(rhY);
        rhPosZ.push_back(rhZ);
        rhTOF.push_back(tof);
		rhTT.push_back(rhTrigTow);
        rhPosEta.push_back(recHitPos.eta());
        rhPosPhi.push_back(recHitPos.phi());
        rhTime.push_back(recHit.time());
        rhTimeErr.push_back(recHit.timeError());
        if( DEBUG ) std::cout << " -- storing values FLAGS" << std::endl;
        rhSubdet.push_back( isEB ? 0 : ( ( idinfo.ecal == ECAL::EP ) ? 1 : 2 ) );
        rhXtalI1.push_back(idinfo.i1);
        rhXtalI2.push_back(idinfo.i2);
        rhisOOT.push_back(recHit.checkFlag(EcalRecHit::kOutOfTime));
        rhEnergy.push_back(recHit.energy());
        //energyError()
        rhSwCross.push_back(swisscross);
        rhisWeird.push_back(recHit.checkFlag(EcalRecHit::kWeird));
        rhisDiWeird.push_back(recHit.checkFlag(EcalRecHit::kDiWeird));
        rhisGS6.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain6));
        rhisGS1.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain1));
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) 
        if( DEBUG ) std::cout << " -- storing values PED" << std::endl;
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
        if( DEBUG ) std::cout << " -- next rechit" << std::endl;

    }//<<>>for (const auto recHit : *recHitsEB_ )   
    nRecHits = nRecHitCnt;

    rhResVecI1.clear();
    rhResVecI2.clear();
    rhResVecEcal.clear();
    rhResVecTT.clear();
    rhResVecE.clear();
    rhResVecadcToGeV.clear();
    rhResVecpedrms12.clear();
    rhResVecTOF.clear();
    rhResVectime.clear();




	// -- Fill output trees ------------------------------------------
	if( DEBUG ) std::cout << "---------- Next Event -----" << std::endl;
	outTree->Fill();

	// -- EOFun ------------------------------------------------------
	//#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//	 ESHandle<SetupData> pSetup;
	//	 iSetup.get<SetupRecord>().get(pSetup);
	//#endif
}//>>>>void GammaResTool::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)


// ------------ method called once each job just before starting event loop	------------
void GammaResTool::beginJob(){

   	// Set up DetIdMap
   	DetIDMap = SetupDetIDs();

	// Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	// Book histograms
	
    //------ 1D Hists --------------------------------------------------------------------------

	//------ 2D Hists --------------------------------------------------------------------------

	std::cout << "Histograms Booked" << std::endl;

	// Create output Tree branches -----------------------------

	// Run, Lumi, Event info
	outTree->Branch("run", &run);
	outTree->Branch("lumi", &lumi);
	outTree->Branch("event", &event, "event/l");


}//>>>>void GammaResTool::beginJob()


// ------------ method called once each job just after ending the event loop	------------
void GammaResTool::endJob(){


}//>>>>void GammaResTool::endJob()


// ------------ method fills 'descriptions' with the allowed parameters for the module	------------
void GammaResTool::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
}//>>>>void GammaResTool::fillDescriptions(edm::ConfigurationDescriptions& descriptions)


//define this as a plug-in
DEFINE_FWK_MODULE(GammaResTool);
