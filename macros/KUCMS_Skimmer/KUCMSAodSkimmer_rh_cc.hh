//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSkimmer.hh"

//-----------------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
//// KUCMSAodSkimmer class ----------------------------------------------------------------------------------------------------
////---------------------------------------------------------------------------------------------------------------------------

//#define DEBUG true
#define DEBUG false

//#define CLSTRMAPS true
#define CLSTRMAPS false

KUCMSAodSkimmer::KUCMSAodSkimmer(){

	LAB = new LabRecoFrame("LAB","LAB");
	S   = new DecayRecoFrame("S","#tilde{S}");
	X2a = new DecayRecoFrame("X2a","#tilde{#chi}_{2a}");
	X2b = new DecayRecoFrame("X2b","#tilde{#chi}_{2b}");
	
	Ja = new VisibleRecoFrame("Ja","jets_{a}");
	Jb = new VisibleRecoFrame("Jb","jets_{b}");
	X1a = new InvisibleRecoFrame("X1a","#tilde{#chi}_{1a}");
	X1b = new InvisibleRecoFrame("X1b","#tilde{#chi}_{1b}");
	
	LAB->SetChildFrame(*S);
	S->AddChildFrame(*X2a);
	S->AddChildFrame(*X2b);
	X2a->AddChildFrame(*X1a);
	X2b->AddChildFrame(*X1b);
	X2a->AddChildFrame(*Ja);
	X2b->AddChildFrame(*Jb);
	
	if(!LAB->InitializeTree()){ std::cout << "Problem initializing tree" << std::endl; }
	
	INV = new InvisibleGroup("INV","Invisible System");
	INV->AddFrame(*X1a);
	INV->AddFrame(*X1b);
	
	InvM = new SetMassInvJigsaw("InvM", "Set inv. system mass");
	INV->AddJigsaw(*InvM);
	
	InvEta = new SetRapidityInvJigsaw("InvEta", "Set inv. system rapidity");
	INV->AddJigsaw(*InvEta);
	InvEta->AddVisibleFrames(S->GetListVisibleFrames());
	
	InvSplit = new MinMassesSqInvJigsaw("InvSplit", "INV -> #tilde{#chi_{1a}}+ #tilde{#chi_{1b}}", 2);
	INV->AddJigsaw(*InvSplit);
	InvSplit->AddVisibleFrame(*Ja, 0);
	InvSplit->AddVisibleFrame(*Jb, 1);
	InvSplit->AddInvisibleFrame(*X1a, 0);
	InvSplit->AddInvisibleFrame(*X1b, 1);
	
	COMB_J =  new CombinatoricGroup("COMB_J", "Combinatoric System of Jets");
	CombSplit_J = new MinMassesSqCombJigsaw("CombSplit_J", "Minimize M_{Va}^{2} + M_{Vb}^{2} ",2,2);
	
	COMB_J->AddFrame(*Ja);
	COMB_J->SetNElementsForFrame(*Ja, 1);
	COMB_J->AddFrame(*Jb);
	COMB_J->SetNElementsForFrame(*Jb, 1);
	
	COMB_J->AddJigsaw(*CombSplit_J);
	CombSplit_J->AddCombFrame(*Ja, 0);
	CombSplit_J->AddCombFrame(*Jb, 1);
	CombSplit_J->AddObjectFrames(X2a->GetListVisibleFrames(), 0);
	CombSplit_J->AddObjectFrames(X2b->GetListVisibleFrames(), 1);// check syntax from example
	 
	if(!LAB->InitializeAnalysis()) std::cout << "Problem initializing analysis tree" << std::endl;

	/*
	  TreePlot tree_plot("TreePlot","TreePlot");

	  for(int t = 0; t < 2; t++){
	  tree_plot.SetTree(*LAB);
	  tree_plot.Draw("ANA_tree", "Reconstruction Tree");

	  tree_plot.SetTree(*COMB_J);
	  tree_plot.Draw("ANA_comb", "Combinatoric Jigsaws for jets");

	  tree_plot.SetTree(*COMB_L);
	  tree_plot.Draw("ANA_comb_L", "Combinatoric Jigsaws for leps");

	  tree_plot.SetTree(*INV);
	  tree_plot.Draw("ANA_inv", "Invisible Jigsaws");
 
	  }
	  tree_plot.WriteOutput("trees.root");
	*/

}//<<>>KUCMSAodSkimmer::KUCMSAodSkimmer()

KUCMSAodSkimmer::~KUCMSAodSkimmer(){

    delete LAB;
    delete S;
    delete X2a;
    delete X2b;
    delete Ja;
    delete Jb;
    delete X1a;
    delete X1b;
      
    delete INV;
    delete InvM;
    delete InvEta;
    delete InvSplit;
      
    delete COMB_J;
    delete CombSplit_J;
      
}//<<>>KUCMSAodSkimmer::~KUCMSAodSkimmer()

void KUCMSAodSkimmer::kucmsAodSkimmer( std::string listdir, std::string eosdir, 
													std::string infilelist, std::string outfilename, bool hasGenInfo, bool genSigPerfect ){

    doGenInfo = hasGenInfo;
    const std::string disphotreename = "tree/llpgtree";
    //std::string inpath, infiles, key, 
	std::string masterstr; 
    //int mct;
	//float crossSection, gmsblam, gmsbct, mcw;
	std::cout << "Processing Input Lists for : " << infilelist << std::endl;
	std::ifstream masterInfile(listdir+infilelist);
	//while( masterInfile >> inpath >> infiles >> key >> crossSection >> gmsblam >> gmsbct >> mcwgt >> mctype ){
    sumEvtGenWgt = 0.0;
    while( std::getline( masterInfile, masterstr ) ){

		if( DEBUG ) std:: cout << masterstr << std::endl;
        if( masterstr[0] == '#' ) continue;
		if( masterstr == " " ) continue;
		auto instrs = splitString( masterstr, " " );
		if( DEBUG ) std:: cout << instrs.size() << std::endl;
        if( instrs.size() < 8 ) continue;

        auto inpath = instrs[0];
		auto infiles = instrs[1];
		auto key = instrs[2];
		auto crossSection = std::stof( instrs[3] );
		auto gmsbgm = std::stof( instrs[4] );
        auto gmsbxm = std::stof( instrs[5] );
        auto mcw = std::stof( instrs[6] );
        auto mct = std::stoi( instrs[7] );		

		std::cout << "Processing Events for : " << infiles << std::endl;
	    std::ifstream infile(listdir+infiles);
	    auto fInTree = new TChain(disphotreename.c_str());
	    std::cout << "Adding files to TChain." << std::endl;
	    std::cout << " - With : " << infiles << " >> " << fInTree << std::endl;
	    std::string str;
		if( not DEBUG ) std::cout << "--  adding files";
	    while( std::getline( infile, str ) ){
	        auto tfilename = eosdir + inpath + str;
			fInTree->Add(tfilename.c_str());
	        if(DEBUG) std::cout << "--  adding file: " << tfilename << std::endl; else std::cout << ".";
			//if(DEBUG) break;
	    }//<<>>while (std::getline(infile,str))
		if( not DEBUG ) std::cout << std::endl;
	
		auto fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
	    auto fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");

        dataSetKey = key;
        xsctn = crossSection;
        gmass = gmsbgm; // = 0 if not gmsb
        xmass = gmsbxm; // = 0 if not gmsb
        mcwgt = mcw; // default 1
        mctype = mct; // 0 = fullsim
	
		KUCMSAodSkimmer::Init( fInTree, doGenInfo );
		initHists();
	    setOutputBranches(fOutTree);	
	
		SetupDetIDsEB(DetIDMap);
		SetupDetIDsEE(DetIDMap);
	
		startJobs();
	
	    std::cout << "Setting up For Main Loop." << std::endl;
		int loopCounter(100000);
	    auto nEntries = fInTree->GetEntries();
	    if(DEBUG){ nEntries = 1000; loopCounter = 100; }
	    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
        nEvents = nEntries;
        int nSelEvts(0);
	    for (Long64_t centry = 0; centry < nEntries; centry++){

	        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
	        auto entry = fInTree->LoadTree(centry);
			if(DEBUG) std::cout << " -- Getting Branches " << std::endl;
			getBranches( entry, doGenInfo );
			geCnts.clear();
            geVars.clear();
			if( genSigPerfect ) geVars.set( "genSigPerfect", 1 ); else geVars.set( "genSigPerfect", 0 );
			if(DEBUG) std::cout << " -- Event Loop " << std::endl;
			auto saveToTree = eventLoop(entry);
			if( saveToTree ){ nSelEvts++; fOutTree->Fill(); }

	    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
        nSelectedEvents = nSelEvts;
		fillConfigTree( fConfigTree ); 

		endJobs();
	
	    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

		auto ext = splitString( infiles, "." );
		std::string extOutFileName( ext[0] + outfilename );
	    TFile* fOutFile = new TFile( extOutFileName.c_str(), "RECREATE" );
	    fOutFile->cd();

		fOutTree->Write();
		fConfigTree->Write();
	
		for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
	    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
	    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }
	
		if( CLSTRMAPS ){
			nMaps = 0;
			for( int it = 0; it < nEBEEMaps; it++ ){ 
				ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
				ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
				ebeeMapR[it]->Write(); delete ebeeMapR[it];
			}//<<>>for( int it = 0; it < nEBEEMaps; it++ )
		}//<<>>f( clsttrMaps )
	
		std::cout << "Finished processing events for : " << infiles << std::endl;
	
	    fOutFile->Close();
	
		delete fInTree;
		delete fOutTree;
		delete fConfigTree;
		delete fOutFile;

	}//<<>>while (std::getline(infile,str))

    std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;

}//<<>>void kucmsSkimmer

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// event loop and startup jobs 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


void KUCMSAodSkimmer::startJobs(){

    std::cout << "Runing StartJobs." << std::endl;

	nEvents = 0;
	nSelectedEvents = 0;

};//<<>>void KUCMSAodSkimmer::startJobs()

bool KUCMSAodSkimmer::eventLoop( Long64_t entry ){

	// counts events and saves event varibles
	// --------------------------------------
    //if( doGenInfo ){ processGenParticles(); }
	//processEvntVars();	
	processRechits();
	//processMet();
	processPhotons();
	processElectrons();
	//processMuons();
	//processJets();
    processRechitsHist();

	// select events to process and store
	//--------------------------------------
	auto saveToTree = eventSelection();	
	//if( saveToTree ) processRJR();
	return saveToTree;

}//<<>>void KUCMSAodSkimmer::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//  do any processing and calulations for objects and save values to output varibles 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

// void KUCMSAodSkimmer::processTemplate(){
// 	
// 		Clear out branch vector varibles &/or initilize other out branch vars
//	------------------------------------------------
// 		Do any calculations / cuts
//  ------------------------------------------------
//		Fill out branch varibles
//
//}//<<>>void KUCMSAodSkimmer::processTemplate()

void KUCMSAodSkimmer::processEvntVars(){

    selEvtVars.clearBranches(); // <<<<<<<   must do

	// calc
    nEvents++;
    selEvtVars.fillBranch( "PVx", PV_x );
    selEvtVars.fillBranch( "PVy", PV_y );
    selEvtVars.fillBranch( "PVz", PV_z ); 
   
	//fill

    selEvtVars.fillBranch( "dsKey", dataSetKey );
    selEvtVars.fillBranch( "evtGenWgt", Evt_genWgt );
    sumEvtGenWgt += Evt_genWgt;

}//<<>>void KUCMSAodSkimmer::processEvntVars()

void KUCMSAodSkimmer::processMet(){

	//intilize
	selMet.clearBranches(); // <<<<<<<   must do

	//calc
	//auto met = std::sqrt(sq2(Met_Cpx)+sq2(Met_Cpy));

    geVars.set("metPx", Met_Cpx );
    geVars.set("metPy", Met_Cpy );

	//fill
	selMet.fillBranch( "met", Met_CPt );
    selMet.fillBranch( "metPx", Met_Cpx );
    selMet.fillBranch( "metPy", Met_Cpy );

}//<<>>void KUCMSAodSkimmer::processMet()

void KUCMSAodSkimmer::processRechits(){

	// initilize
    //selECALRecHit.clearBranches(); // <<<<<<<   must do

	// calc
    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

	rhused.clear();
	rhispho.clear();
	rhisele.clear();
    auto nRecHits = ECALRecHit_ID->size();
    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;
    for( int it = 0; it < nRecHits; it++ ){

		rhispho.push_back(0);
		rhisele.push_back(0);
		rhused.push_back(false);

	}//<<>>for( int it = 0; it < nRecHits; it++ )

	// fill


}//<<>>void KUCMSAodSkimmer::processRechits()

void KUCMSAodSkimmer::processRechitsHist(){

    // initilize
    //selECALRecHit.clearBranches(); // <<<<<<<   must do

    // calc
    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
    //------------ rechits -------------------------

    auto nRecHits = ECALRecHit_ID->size();
    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;
    for( int it = 0; it < nRecHits; it++ ){

        auto id = (*ECALRecHit_ID)[it];
        auto idinfo = DetIDMap[id];
        if( idinfo.ecal == ECAL::EB ){
        //if( true ){

            hist1d[10]->Fill( (*ECALRecHit_energy)[it], 1 );
            hist1d[0]->Fill( (*ECALRecHit_energy)[it], 1 );
            hist1d[1]->Fill( (*ECALRecHit_energy)[it], 1 );
            hist1d[2]->Fill( (*ECALRecHit_energy)[it], 1 );
            //auto radius = hypo( (*rhPosX)[it], (*rhPosY)[it] )

			if( rhispho[it] > 0 ){

				hist2d[0]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_isOOT)[it] );
                hist2d[1]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_swCross)[it] );
                hist2d[2]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_eta)[it] );
                hist2d[3]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_amplitude)[it] );
                hist2d[4]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_time)[it] );
				hist1d[13]->Fill( (*ECALRecHit_energy)[it], rhispho[it] );
				hist1d[12]->Fill( 1 );

			}//<<>>if( rhispho[it] )

            //if( rhisele[it] && not rhispho[it] ){
            if( rhisele[it] > 0 ){

                hist2d[10]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_isOOT)[it] );
                hist2d[11]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_swCross)[it] );
                hist2d[12]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_eta)[it] );
                hist2d[13]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_amplitude)[it] );
                hist2d[14]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_time)[it] );
                hist1d[14]->Fill( (*ECALRecHit_energy)[it], 1 );
                hist1d[15]->Fill( rhisele[it] );

            }//<<>>if( rhispho[it] )

            if( rhispho[it] == 0 && rhisele[it] == 0 ){

                hist2d[20]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_isOOT)[it] );
                hist2d[21]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_swCross)[it] );
                hist2d[22]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_eta)[it] );
                hist2d[23]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_amplitude)[it] );
                hist2d[24]->Fill( (*ECALRecHit_energy)[it], (*ECALRecHit_time)[it] );
				hist1d[17]->Fill( (*ECALRecHit_energy)[it], 1 );

            }//<<>>if( rhispho[it] )

			//if( not rhispho[it] && rhisele[it] ) std::cout << " -- Electron not in Photons !!!!!!!! " << std::endl;

			if( not rhused[it] ){ hist1d[7]->Fill( (*ECALRecHit_energy)[it], 1 ); }


        }//<<>>if( (*rhSubdet)[it] == 0 )

    }//<<>>for( int it = 0; it < nRecHits; it++ )

    // fill

}//<<>>void KUCMSAodSkimmer::processRechits()



void KUCMSAodSkimmer::processGenParticles(){

	// initilize
	selGenPart.clearBranches(); // <<<<<<<   must do

	// calc
    if( DEBUG ) std::cout << "Finding genParticles" << std::endl;
    //------------  genparts ------------------------

	int nSGlue = 0;
	int nSQuark = 0;
	int nGenParts = Gen_pdgId->size();
    for( int it = 0; it < nGenParts; it++ ){

        float energy = (*Gen_energy)[it];
        float eta = (*Gen_eta)[it];
        float phi = (*Gen_phi)[it];
        float pt = (*Gen_pt)[it];
        uInt  pdgId = (*Gen_pdgId)[it];
        int   susId = (*Gen_susId)[it];
        float charge = (*Gen_charge)[it];
        float mass = (*Gen_mass)[it];
        float status = (*Gen_status)[it];
        float vx = (*Gen_vx)[it];
        float vy = (*Gen_vy)[it];
        float vz = (*Gen_vz)[it];
        float px = (*Gen_px)[it];
        float py = (*Gen_py)[it];
        float pz = (*Gen_pz)[it];

		//if( susId == 21 || susId == 20 ) std::cout << " -- Sus pdgId " << pdgId  << " susId " << susId << std::endl;
        //if( pdgId == 1000021 ) std::cout << " -- Sus pdgId " << pdgId  << " susId " << susId << std::endl;
		//if( pdgId > 1000000 && pdgId < 1000007 ) std::cout << " -- Sus pdgId " << pdgId  << " susId " << susId << std::endl;
		if( susId == 97 && ( pdgId > 1000000 && pdgId < 1000007 ) ) nSQuark++;
        if( susId == 97 &&pdgId == 1000021 ) nSGlue++;

		//if( pdgId < 7 ) continue;

		selGenPart.fillBranch( "genPartEnergy", energy );
        selGenPart.fillBranch( "genPartEta", eta );
        selGenPart.fillBranch( "genPartPhi", phi );
        selGenPart.fillBranch( "genPartPt", pt );
        selGenPart.fillBranch( "genPartPdgId", pdgId );
        selGenPart.fillBranch( "genPartSusId", susId );
        selGenPart.fillBranch( "genCharge", charge );
        selGenPart.fillBranch( "genMass", mass );
        selGenPart.fillBranch( "genStatus", status );
        selGenPart.fillBranch( "genVx", vx );
        selGenPart.fillBranch( "genVy", vy );
        selGenPart.fillBranch( "genVz", vz );
        selGenPart.fillBranch( "genPx", px );
        selGenPart.fillBranch( "genPy", py );
        selGenPart.fillBranch( "genPz", pz );

    }//<<>>for( int it = 0; it < nGenParts; it++ )

	// fill
	//std::cout << " ---- Counts nSQuark " << nSQuark << " nSGlue " << nSGlue << std::endl;
	geVars.set( "nSQuark", nSQuark );
	geVars.set( "nSGlue", nSGlue );

}//<<>>void KUCMSAodSkimmer::processGenParticles()

void KUCMSAodSkimmer::processCalojets(){

	// initilize
    // selCaloJets.clearBranches(); // <<<<<<<   must do

	/* calc
	int nCaloJets = 0;
	for( int it = 0; it < nCaloJets; it++ ){


	}//<<>>for( int it = 0; it < nCaloJets; it++ )
	*/
	//fill


}//<<>>void KUCMSAodSkimmer::processCalojets()

void KUCMSAodSkimmer::processPhotons(){

    //bool verbose = true;
    bool verbose = false;

	// intilize
	selPhotons.clearBranches(); // <<<<<<<   must do

	// calc
    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

	std::vector<int> phoOrderIndx;
    uInt nSelPhotons = 0;
	int hfdiff = 0;
    uInt nPhotons = Photon_excluded->size();	
    if( DEBUG || verbose ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( uInt it = 0; it < nPhotons; it++ ){

		if( (*Photon_excluded)[it] ) continue;
        int scIndx = (*Photon_scIndex)[it];
        auto rhlist = (*SuperCluster_rhIds)[scIndx]; // std::vector<std::vector<uInt>>
		//float scetawidth = (*SuperCluster_etaWidth)[scIndx];
		//float scphiwidth = (*SuperCluster_phiWidth)[scIndx];
		//float scdrsize = hypo(scetawidth/2,scphiwidth/2);
		float scdrsize = 0.15;
		//std::cout << " -- pho scdrsize: " << scdrsize << std::endl;
        float eta = (*Photon_eta)[it];
        float phi = (*Photon_phi)[it];
		auto isOOT = (*Photon_isOot)[it];

        int cntRhList = rhlist.size();
        int nRecHits = ECALRecHit_ID->size();
        for( int scrhit = 0; scrhit < cntRhList; scrhit++ ){

            uInt scrhid = rhlist[scrhit];
            for( int rhit = 0; rhit < nRecHits; rhit++ ){

                uInt rhid = (*ECALRecHit_ID)[rhit];
                if( scrhid == rhid ){

                    float rhe = (*ECALRecHit_energy)[rhit];
                    rhused[rhit] = true;
					//rhispho[rhit] = true;
                    rhispho[rhit] += 1;
                    hist1d[3]->Fill( rhe, 1 );
					if( isOOT ){ 

						hist1d[8]->Fill( rhe, 1 );
                   		hist2d[15]->Fill( (*ECALRecHit_energy)[rhit], (*ECALRecHit_isOOT)[rhit] );
                    	hist2d[16]->Fill( (*ECALRecHit_energy)[rhit], (*ECALRecHit_swCross)[rhit] );
                    	hist2d[17]->Fill( (*ECALRecHit_energy)[rhit], (*ECALRecHit_eta)[rhit] );
                    	hist2d[18]->Fill( (*ECALRecHit_energy)[rhit], (*ECALRecHit_amplitude)[rhit] );
                    	hist2d[19]->Fill( (*ECALRecHit_energy)[rhit], (*ECALRecHit_time)[rhit] );

					}//<<>>if( isOOT )


					hfdiff++;
					//break;

                }//<<>> if( rh == id )
				else {

                    float rheta = (*ECALRecHit_eta)[rhit];
                    float rhphi = (*ECALRecHit_phi)[rhit];
                    float dphi12 = dPhi( rhphi, phi );
                    float pho2dr = hypo( rheta-eta, dphi12 );
                    if( pho2dr < scdrsize ){

                    	float rhe = (*ECALRecHit_energy)[rhit];
                        hist1d[5]->Fill( rhe, 1 );

                    }//<<>>if( pho2dr < scdrsize )

				}//<<>> else
            }//<<>>for( int it = 0; it < nRecHits; it+;

        }//<<>>for( auto rh : rhlist )
		hfdiff -= cntRhList;
		//hist1d[100]->Fill( hfdiff, 1 );

		//---------------------------------------------------
    ///////////  pho selection ////////////////////////////////////////////////////////////////////
        if( DEBUG ) std::cout << " -- looping photons : getting pho isEB, has min pt, not excluded, electron veto " << std::endl;
		auto isExcluded = (*Photon_excluded)[it];
		//auto scIndx = (*Photon_scIndex)[it];
        auto isEB = (*SuperCluster_seedIsEB)[scIndx];
		auto hasEleVeto = (*Photon_electronVeto)[it];
        bool hasPixSeed = (*Photon_pixelSeed)[it];

		//if( geVars("genSigPerfect") == 1 ) std::cout << " -- pho sel susid " << (*Gen_susId)[(*Photon_genIdx)[it]] << std::endl;
		bool isGenSig = doGenInfo ? ( (*Gen_susId)[(*Photon_genIdx)[it]] == 22 )  : 0;

        auto pt = (*Photon_pt)[it];
		bool underMinPt = pt < 30;
        //auto eta = (*Photon_eta)[it];
		auto overMaxEta = std::abs(eta) > 1.479;
        //auto phi = (*Photon_phi)[it];

        if( DEBUG ) std::cout << " -- looping photons : getting phojet iso " << std::endl;
        bool isJetPhoton = false;
        for( int jit = 0; jit < Jet_energy->size(); jit++ ){

            bool underMinJPt = (*Jet_pt)[jit] < 75.0;
            bool underMinJQual = getJetQuality(it)  < 2;
			auto jeta = (*Jet_eta)[jit];
            auto jphi = (*Jet_phi)[jit];
			auto overMaxJEta = std::abs(jeta) > 2.4;
            if( underMinJPt || underMinJQual || overMaxJEta ) continue;
			
            float dpjeta = jeta - eta;
            float dpjphi = dPhi( jphi, phi );
            float dr = hypo( dpjeta, dpjphi );
            if( dr < 0.4 ) isJetPhoton = true;

        } // for( int jit = 0; jit < nSelJets; jit++ )

        if( DEBUG ) std::cout << " -- looping photons : getting pho iso " << std::endl;
        auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
        bool passHcalSum = true;
        auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
        bool passTrkSum = tspscdr4 < 6.0; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
		auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
        bool passsEcalRhSum = erhsecdr4 < 10.0;
		auto htoem = (*Photon_hadTowOverEM)[it];
        bool passHOE = htoem < 0.02;
        bool failPhoIso = not ( passHOE && passsEcalRhSum && passTrkSum && passHcalSum );

		bool phoskip = isExcluded || hasPixSeed || overMaxEta || underMinPt || isJetPhoton || failPhoIso;
		//if( geVars("genSigPerfect") == 1 && isGenSig ) std::cout << " -- pho sel: phoskip " << phoskip <<  " isGenSig " << isGenSig << std::endl;  
		//if( geVars("genSigPerfect") == 1 &&  phoskip && isGenSig  ){ 
		//		std::cout << "   -- xsepji: " << isExcluded  << hasPixSeed << overMaxEta << underMinPt << isJetPhoton << failPhoIso << std::endl; }
		if( phoskip ) continue;		
		if( ( geVars("genSigPerfect") == 1 ) &&  ( not isGenSig ) ) continue;

    ///////////  pho selection ////////////////////////////////////////////////////////////////////

        if( DEBUG ) std::cout << " -- pho pull SC info" << std::endl;
        auto scSize = SuperCluster_seedIsEB->size();
        auto rhids = (*SuperCluster_rhIds)[scIndx];
        uInt nrh = rhids.size();

        if( DEBUG ) std::cout << " -- pho pull info" << std::endl;
		//auto isOOT = (*Photon_isOot)[it];
		auto time = (*Photon_seedTOFTime)[it];
		auto smaj = (*SuperCluster_smaj)[scIndx];
        auto smin = (*SuperCluster_smin)[scIndx];
        auto r9 = (*Photon_r9)[it];
        auto sieie = 0.f;//(*Photon_sieie)[it];
        auto energy = (*Photon_energy)[it];
        auto cee = (*SuperCluster_covEtaEta)[scIndx];
        auto cep = (*SuperCluster_covEtaPhi)[scIndx];
        auto cpp = (*SuperCluster_covPhiPhi)[scIndx];
        auto s4 = (*Photon_s4)[it];
        auto salp = (*SuperCluster_salp)[scIndx];
        auto etawide = (*SuperCluster_etaWidth)[scIndx];
        auto phiwide = (*SuperCluster_phiWidth)[scIndx];
        auto sigmaieie = (*Photon_SigmaIEtaIEta)[it];   //!
        auto scx = (*SuperCluster_x_calo)[scIndx];
        auto scy = (*SuperCluster_y_calo)[scIndx];
        auto scz = (*SuperCluster_z_calo)[scIndx];

        int genIdx = -1;
        int momIdx = -1;
        float susId = -1;
		float genpt = -1;
        float momEnergy = -1;   //!
        float momEta = -1;   //!
        float momMass = -1;   //!
        float momPhi = -1;   //!
        float momPt = -1;   //!
        float momPx = -1;   //!
        float momPy = -1;   //!
        float momPz = -1;   //!
        float momVx = -1;   //!
        float momVy = -1;   //!
        float momVz = -1;   //!

		if( doGenInfo ){
			genIdx = (*Photon_genIdx)[it];
			momIdx = (*Photon_genSigMomId)[it];
        	susId = (*Gen_susId)[genIdx];

			if( genIdx > -1.0 ){ 
				genpt = (*Gen_pt)[genIdx]; 
			}//if( genIdx > -1.0 )
			
			if( momIdx > -1.0 ){
                momEnergy = (*Gen_energy)[momIdx];   //!
                momEta = (*Gen_eta)[momIdx];   //!
                momMass = (*Gen_mass)[momIdx];   //!
                momPhi = (*Gen_phi)[momIdx];   //!
                momPt = (*Gen_pt)[momIdx];   //!
                momPx = (*Gen_px)[momIdx];   //!
                momPy = (*Gen_py)[momIdx];   //!
                momPz = (*Gen_pz)[momIdx];   //!
                momVx = (*Gen_vx)[momIdx];   //!
                momVy = (*Gen_vy)[momIdx];   //!
                momVz = (*Gen_vz)[momIdx];   //!
			}//<<>>if( momIdx > -1.0 )

		}//if( doGenInfo )

        auto htsebcdr4 = (*Photon_hcalTowerSumEtBcConeDR04)[it];
        auto tsphcdr3 = (*Photon_trkSumPtHollowConeDR03)[it];
        auto tsphcdr4 = (*Photon_trkSumPtHollowConeDR04)[it];
        //auto tspscdr4 = (*Photon_trkSumPtSolidConeDR04)[it];
        //auto erhsecdr4 = (*Photon_ecalRHSumEtConeDR04)[it];
        //auto htoem = (*Photon_hadTowOverEM)[it];
        auto hoem = (*Photon_hadOverEM)[it];
        auto sieip = 0.f;//(*Photon_sieip)[it];
        auto sipip = 0.f;//(*Photon_sipip)[it];
        auto phoOOT = (*Photon_isOot)[it];
        //auto htsecdr4 = (*Photon_hcalTowerSumEtConeDR04)[it];   //!
        auto nthcdr4 = (*Photon_nTrkHollowConeDR04)[it];   //!
        auto ntscdr4 = (*Photon_nTrkSolidConeDR04)[it];   //!
        auto epfciso = (*Photon_ecalPFClusterIso)[it];   //!
        auto hct = (*hasConversionTracks)[it];   //!
        auto hpfciso = (*Photon_hcalPFClusterIso)[it];   //!
        auto hoepuc = (*Photon_Hoe_PUcorr)[it];   //!
        auto pfciso = (*Photon_pfChargedIso)[it];   //!
        auto pscisopfpv = (*Photon_pfChargedIsoPFPV)[it];   //!
        auto pfcisowv = (*Photon_pfChargedIsoWorstVtx)[it];   //!
        auto pfphoisso = (*Photon_pfPhoIso03)[it];   //!
        auto pfriso3aq = (*pfRelIso03_all_quadratic)[it];   //!
        auto pfrtso3cq = (*pfRelIso03_chg_quadratic)[it];   //!

        //--------------------------------------------------------------
        if( DEBUG ) std::cout << " -- pho get calclated values" << std::endl;
        int phoQuality = getPhoQuality(it);
        auto phoClstrR9 = clstrR9( rhids );
        //auto phoEigens2D = getRhGrpEigenFromAngles( rhids );
        //auto evaluegeo = phoEigens2D[2];
		//auto geosmaj = phoEigens2D[3];
        //auto geosmin = phoEigens2D[4];

        float phoPhoIsoDr = 10.0;
        for( uInt it2 = 0; it2 < nPhotons; it2++ ){
            if( it == it2 ) continue;
        	float eta2 = (*Photon_eta)[it2];
        	float phi2 = (*Photon_phi)[it2];
			auto dphi12 = dPhi( phi, phi2 );
			auto pho2dr = hypo( eta-eta2, dphi12 ); 
            //std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << pho2dr << std::endl;
			if( pho2dr < phoPhoIsoDr ) phoPhoIsoDr = pho2dr;
		}//for( uInt it2 = it+1; it2 < nPhotons; it2++ )

/*
        // pho object selection ------------------------------------------
        if( DEBUG ) std::cout << " -- pho obj selection" << std::endl;
		auto isMinMedQuality = phoQuality > 1;
		auto underMaxSMaj = smaj == 1.3;
        auto underMaxSMin = smin <= 0.4;
		auto overMinR9 = r9 >= 0.9;
		auto underMaxSieie = sieie <= 0.014;
        auto overMinRhCnt = nrh >= 20;

		//auto phoSelected = isMinMedQuality;
        auto phoSelected = true;
		if( not phoSelected ) continue;
*/
		//float ordpt = pt;
		//while( phoOrderIndx.count(ordpt) > 0 ){ ordpt -= 0.000001; }
		//if( DEBUG ) std::cout << " -- setting pho index : " << it << " for pt : " << ordpt << std::endl;
		phoOrderIndx.push_back(it);
        //if( DEBUG ) std::cout << " -- pho index set : " << phoOrderIndx[ordpt] << std::endl;

		// fill ( vectors )
		if( DEBUG ) std::cout << " -- pho fill out branches" << std::endl;
        selPhotons.fillBranch( "selPhoOOT", phoOOT );
        selPhotons.fillBranch( "selPhoEleVeto", hasEleVeto );
		selPhotons.fillBranch( "selPhoQuality", phoQuality );
        selPhotons.fillBranch( "selPhoTime", time );
        selPhotons.fillBranch( "selPhoSCx", scx );
        selPhotons.fillBranch( "selPhoSCy", scy );
        selPhotons.fillBranch( "selPhoSCz", scz );
		//selPhotons.fillBranch( "selPhoGeoEgnVal", evaluegeo );
        selPhotons.fillBranch( "selPhoEta", eta );
        selPhotons.fillBranch( "selPhoPhi", phi );
        selPhotons.fillBranch( "selPhoPt", pt );
        selPhotons.fillBranch( "selPhoGenPt", genpt );
        selPhotons.fillBranch( "selPhoSMaj", smaj );
        selPhotons.fillBranch( "selPhoSMin", smin );
        selPhotons.fillBranch( "selPhoClstrRn", phoClstrR9 );
        selPhotons.fillBranch( "selPhoR9", r9 );
        selPhotons.fillBranch( "selPhoSieie", sieie );
        selPhotons.fillBranch( "selPhoSieip", sieip );
        selPhotons.fillBranch( "selPhoSipip", sipip );
        selPhotons.fillBranch( "selPhoCovEtaEta", cee );
        selPhotons.fillBranch( "selPhoCovEtaPhi", cep );
        selPhotons.fillBranch( "selPhoCovPhiPhi", cpp );
        selPhotons.fillBranch( "selPhoS4", s4 );
        selPhotons.fillBranch( "selPhoSAlp", salp );
        selPhotons.fillBranch( "selPhoEtaWidth", etawide );
        selPhotons.fillBranch( "selPhoPhiWidth", phiwide );
        selPhotons.fillBranch( "selPhoNrh", nrh );
        selPhotons.fillBranch( "selPhoEnergy" , energy );
        selPhotons.fillBranch( "selPhoSusyId" , susId );
        selPhotons.fillBranch( "selPhoSigmaIEtaIEta" , sigmaieie );
        //selPhotons.fillBranch( "selPhoGeoSMaj", geosmaj );
        //selPhotons.fillBranch( "selPhoGeoSMin", geosmin );
        selPhotons.fillBranch( "selPhoHcalTowerSumEtBcConeDR04", htsebcdr4 );
        selPhotons.fillBranch( "selPhoTrkSumPtHollowConeDR03", tsphcdr3 );
        selPhotons.fillBranch( "selPhoTrkSumPtHollowConeDR04", tsphcdr4 );
        selPhotons.fillBranch( "selPhoTrkSumPtSolidConeDR04", tspscdr4 );
        selPhotons.fillBranch( "selPhoPixelSeed", hasPixSeed );
        selPhotons.fillBranch( "selPhoEcalRHSumEtConeDR04", erhsecdr4 );
        selPhotons.fillBranch( "selPhoHadTowOverEM", htoem );
        selPhotons.fillBranch( "selPhoHadOverEM", hoem );
        selPhotons.fillBranch( "selPhoHcalTowerSumEtConeDR04", htsecdr4  );   //!
        selPhotons.fillBranch( "selPhoNTrkHollowConeDR04", nthcdr4  );   //!
        selPhotons.fillBranch( "selPhoNTrkSolidConeDR04", ntscdr4  );   //!
        selPhotons.fillBranch( "selPhoEcalPFClusterIso", epfciso  );   //!
        selPhotons.fillBranch( "selPhoHasConversionTracks", hct  );   //!
        selPhotons.fillBranch( "selPhoHcalPFClusterIso", hpfciso  );   //!
        selPhotons.fillBranch( "selPhoHoe_PUcorr", hoepuc  );   //!
        selPhotons.fillBranch( "selPhoPfChargedIso", pfciso  );   //!
        selPhotons.fillBranch( "selPhoPfChargedIsoPFPV", pscisopfpv  );   //!
        selPhotons.fillBranch( "selPhoPfChargedIsoWorstVtx", pfcisowv  );   //!
        selPhotons.fillBranch( "selPhoPfPhoIso03", pfphoisso  );   //!
        selPhotons.fillBranch( "selPhoPfRelIso03_all_quadratic", pfriso3aq  );   //!
        selPhotons.fillBranch( "selPhoPfRelIso03_chg_quadratic", pfrtso3cq  );   //!

        selPhotons.fillBranch( "selPhoPhoIsoDr", phoPhoIsoDr );
        selPhotons.fillBranch( "selPhoGenIdx", genIdx );
        selPhotons.fillBranch( "selPhoGenSigMomEnergy", momEnergy  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomEta", momEta  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomMass", momMass  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPhi", momPhi  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPt", momPt  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPx", momPx  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPy", momPy  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomPz", momPz  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomVx", momVx  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomVy", momVy  );   //!
        selPhotons.fillBranch( "selPhoGenSigMomVz", momVz  );   //!
		//if( verbose ) std::cout << " -- selPho Pt: " << pt << " phi: " << phi << " geo: " << evaluegeo << " clrn: " << phoClstrR9;
		if( verbose ) std::cout << " nrh: " << nrh << " quality: " << phoQuality << std::endl;

        nSelPhotons++;
    }//<<>>for( int it = 0; it < nPhotons; it++ )
    if( DEBUG ) std::cout << " -- pho loop finished" << std::endl;
	// fill ( other )

/*
	std::vector<unsigned int> phoPtOrder;
	if( phoOrderIndx.size() > 0 ){ 
		for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ ){ 
			phoPtOrder.push_back(phoptit->second); 
			if( DEBUG || verbose ) std::cout << " phoptit pushback : " << phoptit->second << std::endl;
		}//for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ )
	}//if( phoOrderIndx.size() > 0 )
*/

    if( DEBUG || verbose ) std::cout << " - Selected " << nSelPhotons << " photons" << std::endl;
    geCnts.set( "nSelPhotons", nSelPhotons );	
	selPhotons.fillBranch( "nSelPhotons",  nSelPhotons );
    selPhotons.fillBranch( "nPhotons", nPhotons );	
    //selPhotons.fillBranch( "nEleVeto", nEleVeto );
    //selPhotons.fillBranch( "selPhoPtOrder", phoPtOrder );

    //if( DEBUG || verbose ) std::cout << " - Setting Leading/Subleading Photon with " << phoPtOrder.size() << std::endl;
	uInt leadPhoIdx = ( nSelPhotons >= 1 ) ? phoOrderIndx[0]  : 9999;
	uInt subLeadPhoIdx = ( nSelPhotons >= 2 ) ? phoOrderIndx[1] : 9999;
    geCnts.set("leadPho",leadPhoIdx);
    geCnts.set("subLeadPho",subLeadPhoIdx);
	//selPhotons.fillBranch( "leadSelPho", leadPho );
	//selPhotons.fillBranch( "subLeadSelPho", subLeadPho );
	//if( verbose ) std::cout << " -- pho lead & sublead selected : " << leadPho << " - " << subLeadPho << std::endl;
	if( DEBUG ) std::cout << " -- pho lead & sublead idx selected : " << leadPhoIdx << " - " << subLeadPhoIdx << std::endl;	
    float lPhoPt = ( nSelPhotons > 0 ) ? (*Photon_pt)[leadPhoIdx] : 0.f;
    geVars.set( "leadPhoPt", lPhoPt );
    float slPhoPt = ( nSelPhotons > 1 ) ? (*Photon_pt)[subLeadPhoIdx] : 0.f;
    geVars.set( "subLeadPhoPt", slPhoPt);
    float lPhoPhi = ( nSelPhotons > 0 ) ? (*Photon_phi)[leadPhoIdx] : 0.f;
    geVars.set( "leadPhoPhi", lPhoPhi );
    float slPhoPhi = ( nSelPhotons > 1 ) ? (*Photon_phi)[subLeadPhoIdx] : 0.f;
    geVars.set( "subLeadPhoPhi", slPhoPhi );

}//<<>>void KUCMSAodSkimmer::processPhoton(){

void KUCMSAodSkimmer::processElectrons(){

    if( DEBUG ) std::cout << "Finding electrons" << std::endl;
	//-------- electrons --------------------------------------

    // intilize
    selElectrons.clearBranches(); // <<<<<<<   must do

	uInt nElectrons = Electron_energy->size();
	float nEpDrMatch = 0;
    float nEpSeedIdMatch = 0;
    float nEpsidsolo = 0;
    float elePhoIsoMinDr(10.0); 
    if( DEBUG ) std::cout << " -- Looping electrons: " << std::endl;
    for( uInt itr = 0; itr < nElectrons; itr++ ){

        int scIndx = (*Electron_scIndex)[itr];
        auto rhlist = (*SuperCluster_rhIds)[scIndx]; // std::vector<std::vector<uInt>>
        //float scetawidth = (*SuperCluster_etaWidth)[scIndx];
        //float scphiwidth = (*SuperCluster_phiWidth)[scIndx];
        //float scdrsize = hypo(scetawidth/2,scphiwidth/2);
        float scdrsize = 0.15;
        //std::cout << " -- ele scdrsize: " << scdrsize << std::endl;
        float eta = (*Electron_eta)[itr];
        float phi = (*Electron_phi)[itr];

		int cntRhList = rhlist.size();
		int nRecHits = ECALRecHit_ID->size();
        for( int scrhit = 0; scrhit < cntRhList; scrhit++ ){

			uInt scrhid = rhlist[scrhit];
            for( int rhit = 0; rhit < nRecHits; rhit++ ){
				
				float rhe = (*ECALRecHit_energy)[rhit];
                uInt rhid = (*ECALRecHit_ID)[rhit];
                if( scrhid == rhid ){

                    rhused[rhit] = true;
					//rhisele[rhit] = true;
                    rhisele[rhit] += 1;
                    hist1d[4]->Fill( rhe, 1 );

					//break;

                }//<<>> if( rh == id )
                else {

                    float rheta = (*ECALRecHit_eta)[rhit];
                    float rhphi = (*ECALRecHit_phi)[rhit];
                    float dphi12 = dPhi( rhphi, phi );
                    float pho2dr = hypo( rheta-eta, dphi12 );
                    if( pho2dr < scdrsize ){

                        hist1d[6]->Fill( rhe, 1 );

                    }//<<>>if( pho2dr < scdrsize )

                }//<<>> else

            }//<<>>for( int it = 0; it < nRecHits; it+;

        }//<<>>for( auto rh : rhlist )

    	if( DEBUG ) std::cout << " ---- Processing  electron: " << itr << " of: " << nElectrons << std::endl;
		//float eta = (*Electron_eta)[itr];
    	//float phi = (*Electron_phi)[itr];
		auto scIdx = (*Electron_scIndex)[itr];
		if( DEBUG ) std::cout << " -- doing sc sstuff " << std::endl;
		auto esid = getLeadRhID( (*SuperCluster_rhIds)[scIdx] );
		//std::cout << " -- esid : " << esid << std::endl;
		auto eletime = (*Electron_seedTOFTime)[itr];
		if( std::abs(eletime) > 30 ) continue;
    	uInt nPhotons = Photon_excluded->size();
		float elePhoIsoDr(0.5);
		bool epSeedIdMatch = false;
		bool epsidsolo = true;
		if( DEBUG ) std::cout << " -- doing ele pho iso sstuff " << std::endl;
    	for( uInt pit2 = 0; pit2 < nPhotons; pit2++ ){
        	float eta2 = (*Photon_eta)[pit2];
        	float phi2 = (*Photon_phi)[pit2];
			auto scIdx2 = (*Photon_scIndex)[pit2];
			auto dphi12 = dPhi( phi, phi2 );
        	float elephodr = hypo( eta-eta2, dphi12 );
        	//std::cout << " -- Eta : " << eta-eta2 << " Phi : " << dPhi(phi,phi2) << " dr : " << elephodr << std::endl;
        	if( elephodr < elePhoIsoDr ) elePhoIsoDr = elephodr;
            //std::cout << " -- scIdx : " << scIdx2 << std::endl; 
			auto rhIds = (*SuperCluster_rhIds)[scIdx2];
			//std::cout << " -- scRh1 : " << rhIds[0] << std::endl;
			auto psid = getLeadRhID(rhIds);
			//std::cout << " -- psid : " << psid << std::endl;
			if( esid == psid ){ if( epSeedIdMatch == true ) epsidsolo = false ; epSeedIdMatch = true; }
    	}//for( uInt pit2 = 0; pit2 < nPhotons; pit2++ )
		elePhoIsoMinDr = elePhoIsoDr;
		bool epDrMatch = ( elePhoIsoDr < 0.2  ) ? true : false;

		if( epDrMatch ) nEpDrMatch++;
		if( epSeedIdMatch ) nEpSeedIdMatch++;
		if( epsidsolo ) nEpsidsolo++;
		if( DEBUG ) std::cout << " ---- next electrons: " << std::endl;

	}//<<>>for( int itr = 0; itr < nElectrons; itr++ )
    if( DEBUG ) std::cout << " -- Finishd looping electrons " << std::endl;

    selElectrons.fillBranch( "elePhoMinDr", elePhoIsoMinDr );
    selElectrons.fillBranch( "nElectrons", nElectrons );
	if( nElectrons == 0 ) nElectrons = 1;
    selElectrons.fillBranch( "epDrMatch", nEpDrMatch/nElectrons );
    selElectrons.fillBranch( "epSeedIdMatch", nEpSeedIdMatch/nElectrons );
    selElectrons.fillBranch( "epSeedIdUnique", nEpsidsolo/nElectrons );

}//<<>>void KUCMSAodSkimmer::processElectrons

void KUCMSAodSkimmer::processJets(){

	// intilize
	selJets.clearBranches(); // <<<<<<<   must do

	// calc
    //if( DEBUG ) std::cout << "Finding jets" << std::endl;
	//--------- jets --------------------------

	uInt nSelJets = 0;
	uInt nJets = Jet_energy->size();
    std::vector<float> seljetpt;
    std::vector<float> seljeteta;
    std::vector<float> seljetphi;
    std::vector<float> seljetmass;
    if( DEBUG ) std::cout << " - Looping over for " << nJets << " jets" << std::endl;
    for( uInt it = 0; it < nJets; it++ ){

		// pull values ---------------------------------------------------
		auto energy = (*Jet_energy)[it];
		auto mass = (*Jet_mass)[it];
        //if( DEBUG ) std::cout << " - Finding Jet Quality" << std::endl;
		auto quality = getJetQuality(it);
		auto pt = (*Jet_pt)[it];
		auto eta = (*Jet_eta)[it];
        auto phi = (*Jet_phi)[it];
		auto rhids = (*Jet_drRhIds)[it];
		if( DEBUG ) std::cout << " -- jettDrId size: " << rhids.size() << std::endl;		

        auto area = (*Jet_area)[it];
        auto chEmEF = (*Jet_chEmEF)[it];
        auto chHEF = (*Jet_chHEF)[it];
        auto chHM = (*Jet_chHM)[it];
        auto muEF = (*Jet_muEF)[it];
        auto neEmEF = (*Jet_neEmEF)[it];
        auto neHEF = (*Jet_neHEF)[it];
        auto neHM = (*Jet_neHM)[it];

        float susId = 0; //*Jet_genLlpId)[it];
        float jgdpt = 0; //*Jet_genDptMatch)[it];
        float jgdr = 0; //*Jet_genDrMatch)[it];
        float jge = 0; //*Jet_genEnergy)[it];
        float jgeta = 0; //*Jet_genEta)[it];
        float jgimpang = 0; //*Jet_genImpactAngle)[it];
        float llpdp = 0; //*Jet_genLlpDp)[it];
        float llpdr = 0; //*Jet_genLlpDr)[it];
        float jgpt = 0; //*Jet_genPt)[it];
        float jgtof = 0; //*Jet_genTOF)[it];
        float jgt = 0; //(*Jet_genTime)[it];
        float jgllpt = 0; //*Jet_genTimeLLP)[it];
        if( doGenInfo ){
        	susId = (*Jet_genLlpId)[it];
        	jgdpt = (*Jet_genDptMatch)[it];
        	jgdr = (*Jet_genDrMatch)[it];
        	jge = (*Jet_genEnergy)[it];
        	jgeta = (*Jet_genEta)[it];
        	jgimpang = (*Jet_genImpactAngle)[it];
        	llpdp = (*Jet_genLlpDp)[it];
        	llpdr = (*Jet_genLlpDr)[it];
        	jgpt = (*Jet_genPt)[it];
        	jgtof = (*Jet_genTOF)[it];
        	jgt = (*Jet_genTime)[it];
        	jgllpt = (*Jet_genTimeLLP)[it];
		}//if( doGenInfo )

		// get cacluated values-------------------------------------------
		if( DEBUG ) std::cout << " - Finding Jet Time." << std::endl;
		auto rhenergies = getRhGrpEnergies( rhids );
		auto rhtimes = getRhGrpTimes( rhids );
		auto timedist = getDistStats( rhtimes, rhenergies );
		auto time = timedist[6];

		if( DEBUG ) std::cout << " - Jet Obj selection." << std::endl;
		// jet object selection ------------------------------------------
		auto overMinPt = pt > 30; 
		auto underMaxEta = std::abs(eta) <= 2.4;
		auto isMinQuality = quality > 1; // 2 = "tight" 3 = "tighter"

		auto jetSelected = underMaxEta && isMinQuality && overMinPt;
		if( not jetSelected ) continue;
		nSelJets++;

        seljetpt.push_back(pt);
        seljeteta.push_back(eta);
        seljetphi.push_back(phi);
        seljetmass.push_back(mass);

		// fill vectors
		selJets.fillBranch( "selJetSusyId", susId );
		selJets.fillBranch( "selJetQuality", quality );
        selJets.fillBranch( "selJetPt", pt);
        selJets.fillBranch( "selJetMass", mass);
        selJets.fillBranch( "selJetEnergy", energy);
        selJets.fillBranch( "selJetEta", eta);
        selJets.fillBranch( "selJetPhi", phi);
        selJets.fillBranch( "selJetTime", time);

        selJets.fillBranch( "selJetArea", area ); //*Jet_area)[it]; 
        selJets.fillBranch( "selJetChEmEF", chEmEF ); //*Jet_chEmEF)[it]; 
        selJets.fillBranch( "selJetchHEF", chHEF ); //*Jet_chHEF)[it]; 
        selJets.fillBranch( "selJetChHM", chHM ); //*Jet_chHM)[it]; 
        selJets.fillBranch( "selJetMuEF", muEF ); //*Jet_muEF)[it]; 
        selJets.fillBranch( "selJetNeEmEF", neEmEF ); //*Jet_neEmEF)[it]; 
        selJets.fillBranch( "selJetNeHEF", neHEF ); //*Jet_neHEF)[it]; 
        selJets.fillBranch( "selJetNeHM", neHM ); //*Jet_neHM)[it]; 

        selJets.fillBranch( "selGenJetDpt", jgdpt ); //*Jet_genDptMatch)[it]; 
        selJets.fillBranch( "selGenJetdr", jgdr ); //*Jet_genDrMatch)[it]; 
        selJets.fillBranch( "selGenJetEnergy", jge ); //*Jet_genEnergy)[it]; 
        selJets.fillBranch( "selGenJeteta", jgeta ); //*Jet_genEta)[it]; 
        selJets.fillBranch( "selGenJetImpAng", jgimpang ); //*Jet_genImpactAngle)[it]; 
        selJets.fillBranch( "selJetLlpDp", llpdp ); //*Jet_genLlpDp)[it]; 
        selJets.fillBranch( "selJetLlpDr", llpdr ); //*Jet_genLlpDr)[it]; 
        selJets.fillBranch( "selGenJetPt", jgpt ); //*Jet_genPt)[it]; 
        selJets.fillBranch( "selGenJetTof", jgtof ); //*Jet_genTOF)[it]; 
        selJets.fillBranch( "selGenJetTime", jgt ); // (*Jet_genTime)[it]; 
        selJets.fillBranch( "selGenJetLlpTime", jgllpt ); //*Jet_genTimeLLP)[it]; 

	}//<<>>for( int it = 0; it < nJets; it++ )
    if( DEBUG ) std::cout << " - Finished Jet loop." << std::endl;

	// fill other
	geVects.set( "selJetPt", seljetpt );
    geVects.set( "selJetEta", seljeteta );
    geVects.set( "selJetPhi", seljetphi );
    geVects.set( "selJetMass", seljetmass );
	geCnts.set( "nSelJets", nSelJets );
    selJets.fillBranch( "nJets", nJets);
	selJets.fillBranch( "nSelJets", nSelJets );

}//<<>>void KUCMSAodSkimmer::processJets()

//------------------------------------------------------------------------------------------------------------
// process RJR for event
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::processRJR(){

	//bool verbose = true;
    bool verbose = false;

	if( DEBUG || verbose ) std::cout << "Processing RJR event varibles." << std::endl;

	// intilize out branches
	selRjrVars.clearBranches(); // <<<<<<<   must do

	// process event

	LAB->ClearEvent();

	auto leadSelPho = geCnts("leadPho"); //selPhotons.getUIBranchValue("leadSelPho");
    auto subLeadSelPho = geCnts("subLeadPho"); //selPhotons.getUIBranchValue("subLeadSelPho");
    auto leadPhoPt = geVars("leadPhoPt"); //selPhotons.getFLBranchValue( "selPhoPt", leadSelPho );
    auto subLeadPhoPt = geVars("subLeadPhoPt"); //selPhotons.getFLBranchValue( "selPhoPt", subLeadSelPho );
    auto leadPhoPhi = geVars("leadPhoPhi"); //selPhotons.getFLBranchValue( "selPhoPhi", leadSelPho );
    auto subLeadPhoPhi = geVars("subLeadPhoPhi"); //selPhotons.getFLBranchValue( "selPhoPhi", subLeadSelPho );
    auto phoRMetCPx = geVars("metPx"); //selMet.getFLBranchValue("metPx");
    auto phoRMetCPy = geVars("metPy"); //selMet.getFLBranchValue("metPy");
    auto nSelPhotons = geCnts("nSelPhotons");

	int n_selphos = 1;
	if( DEBUG ) std::cout << " - Loading MET." << std::endl;
	phoRMetCPx += leadPhoPt*std::cos(leadPhoPhi); 
    phoRMetCPy += leadPhoPt*std::sin(leadPhoPhi);
    if( nSelPhotons > 1 ){ 
		n_selphos = 2;
		phoRMetCPx += subLeadPhoPt*std::cos(subLeadPhoPhi);
    	phoRMetCPy += subLeadPhoPt*std::sin(subLeadPhoPhi);
	}//<<>>if( nSelPhotons > 1 )
	TVector3 ETMiss(phoRMetCPx,phoRMetCPy,0);
	if( verbose ){
		std::cout << " - Loading MET lPt: " << leadPhoPt << " lPhi: " << leadPhoPhi << std::endl;
    	std::cout << " - Loading MET slPt: " << subLeadPhoPt << " slPhi: " << subLeadPhoPhi << std::endl; 
		std::cout << " - Loading MET x: " << geVars("metPx") << " -> " << phoRMetCPx;
        std::cout << " y: " << geVars("metPy") << " -> " << phoRMetCPy << std::endl;
	}//<<>>if( verbose )
	INV->SetLabFrameThreeVector(ETMiss);

    auto nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
    auto selJetPt = geVects( "selJetPt");
    auto selJetEta = geVects( "selJetEta");
    auto selJetPhi = geVects( "selJetPhi");
    auto selJetMass = geVects( "selJetMass");
	if( DEBUG ) std::cout << " - Loading Jets." << std::endl;
	std::vector<RFKey> jetID;
  	for( uInt it = 0; it < nSelJets; it++ ){
		auto sjetPt = selJetPt[it]; //selJets.getFLBranchValue( "selJetPt", it );
		auto sjetEta = selJetEta[it]; //selJets.getFLBranchValue( "selJetEta", it );
		auto sjetPhi = selJetPhi[it]; //selJets.getFLBranchValue( "selJetPhi", it ); 
        auto sjetMass = ( selJetMass[it] > 0 ) ? selJetMass[it] : 0; //selJets.getFLBranchValue( "selJetMass", it );   
		TLorentzVector jet;
		jet.SetPtEtaPhiM( sjetPt, sjetEta, sjetPhi, sjetMass );
		if( verbose ) std::cout << " - Loading Jet Pt: " << sjetPt << " Eta: " << sjetEta;
		if( verbose ) std::cout << " Phi: " << sjetPhi << " M: " << sjetMass << std::endl;
		jetID.push_back(COMB_J->AddLabFrameFourVector(jet)); 
	}//<<>>for( int i = 0; i < nSelJets; i++ )

  	if( !LAB->AnalyzeEvent() ) std::cout << "Something went wrong with tree event analysis" << std::endl;
	
	if( DEBUG ) std::cout << " - Getting RJR varibles." << std::endl;

  	float m_MS = S->GetMass();
  	//float m_PS = S->GetMomentum(*CM);
  	float m_cosS  = S->GetCosDecayAngle();
  	float m_dphiS = S->GetDeltaPhiDecayAngle();
  	float m_dphiSI  = S->GetDeltaPhiBoostVisible();
  	//float m_PTS = S->GetFourVector().Pt();
  	//float m_PzS = S->GetFourVector().Pz();
  
  	float m_MX2a = X2a->GetMass();
  	float m_cosX2a = X2a->GetCosDecayAngle();
  	float m_MX2b = X2b->GetMass();
  	float m_cosX2b = X2b->GetCosDecayAngle();

  	float m_EVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).E();
  	float m_EVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).E();
  	float m_PVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).P();
  	float m_PVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).P();

  	float m_MX1a = X1a->GetMass();
  	float m_cosX1a = X1a->GetCosDecayAngle();
  	float m_MX1b = X1b->GetMass();
  	float m_cosX1b = X1b->GetCosDecayAngle();  

  	float m_MV = S->GetListVisibleFrames().GetMass();
  	float m_PV = S->GetListVisibleFrames().GetFourVector(*S).P();
  	float m_MVa = X2a->GetListVisibleFrames().GetMass();
  	float m_MVb = X2b->GetListVisibleFrames().GetMass();

  	float m_PV_lab    = S->GetListVisibleFrames().GetFourVector().P();
  	float m_dphiMET_V = S->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

	// fill branches
 	selRjrVars.fillBranch( "X1aMass", m_MX1a );
    selRjrVars.fillBranch( "X1aCosA", m_cosX1a );
    selRjrVars.fillBranch( "X1bMass", m_MX1b );
    selRjrVars.fillBranch( "X1bCosA", m_cosX1b );

    selRjrVars.fillBranch( "X2aMass", m_MX2a );
    selRjrVars.fillBranch( "X2aCosA", m_cosX2a );
    selRjrVars.fillBranch( "X2bMass", m_MX2b );
    selRjrVars.fillBranch( "X2bCosA", m_cosX2b );

    selRjrVars.fillBranch( "SMass", m_MS );
    selRjrVars.fillBranch( "SCosA", m_cosS );

   selRjrVars.fillBranch( "RjrMET", hypo(phoRMetCPx,phoRMetCPy) ); 
   selRjrVars.fillBranch( "nRjrPhotons", n_selphos );

}

//------------------------------------------------------------------------------------------------------------
// decide which events to save to tree
//------------------------------------------------------------------------------------------------------------

bool KUCMSAodSkimmer::eventSelection(){
// select which events to save and fill output branches

	if( DEBUG ) std::cout << "Event selection." << std::endl;
	// determine if we want to save event

	//float selmet; selMet.getBranch( "Met", selmet );
	//std::cout << " Event Met : " << selmet << std::endl;
/*
	int nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
    float nSelPhotons = geCnts("nSelPhotons"); //selPhotons.getUIBranchValue("nSelPhotons");
	float leadSelPho = geCnts("leadPho"); //selPhotons.getUIBranchValue("leadSelPho");
    float subLeadSelPho = geCnts("subLeadPho"); //selPhotons.getUIBranchValue("subLeadSelPho");
	float leadPhoPt = ( nSelPhotons > 0 ) ? geVars("leadPhoPt") : 0;
    float subLeadPhoPt = ( nSelPhotons > 1 ) ? geVars("subLeadPhoPt") : 0;
	if( DEBUG ) std::cout << " - Lead/Sublead Photons: " << leadSelPho << " - " << subLeadSelPho << std::endl;

    bool gt1phos = nSelPhotons >= 1;
    bool gt2jets = nSelJets >= 2;
    bool gt2phos = nSelPhotons >= 2;
	bool leadPhoPt70 = leadPhoPt >= 70;
    bool leadPhoPt30 = leadPhoPt >= 30;
	bool subLeadPhoPt40 = subLeadPhoPt >= 40; 
	
    bool evtSelected = gt2jets && gt1phos && leadPhoPt30;
	//auto evtSelected = leadPhoPt70 && subLeadPhoPt40 && gt2jets && gt2phos;
	//if( geVars( "genSigPerfect" ) == 1 ) std::cout << "  -- Event sel nSGlue : " << geVars( "nSGlue" ) << std::endl;
	if( ( geVars( "genSigPerfect" ) == 1 ) && ( geVars( "nSGlue" ) < 1 ) ) evtSelected = false;
*/

	bool evtSelected = false;
	if( DEBUG ){ 
		if( evtSelected ) std::cout << " ---------------- Event Passed !!!!!!!" << std::endl; 
		else std::cout << " --------------------- Event Failed." << std::endl;}
	//if( evtSelected ) std::cout << " ---------------- Event Passed !!!!!!!" << std::endl;
	if( evtSelected ){ nSelectedEvents++; return true; } 
	else return false;

}//<<>>void KUCMSAodSkimmer::eventSelection()

//------------------------------------------------------------------------------------------------------------
// object quality ids
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getPhoQuality( int it ){

    // determine pog quality class
    // -----------------------------------------------------
    if( DEBUG ) std::cout << " -- pho id" << std::endl;
    //bool rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= ( 0.006*(*Photon_pt)[it] + 4.2 );
    bool rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= 10.0;
    //bool hcalTowIso = (*Photon_hcalTowerSumEtBcConeDR04)[it] <= ( 0.0025*(*Photon_pt)[it] + 2.2 );
    //if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
    //bool hcTrkIsoL = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 3.5 ); //hallow cone track iso
    //bool hcTrkIsoT = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 2 ); //hallow cone track iso
    bool hcTrkIso = (*Photon_trkSumPtSolidConeDR04)[it] <= 6.0;
    //if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
    //bool hadOverE = (*Photon_hadTowOverEM)[it] <= 0.05;
    bool hadOverE = (*Photon_hadTowOverEM)[it] <= 0.02;
    //auto sieie = (*Photon_sieie)[it];
    //bool sigmaIeieEE = sieie <= 0.03; // tight only EE
    //bool sigmaIeieEB = sieie <= 0.013; // tight only EB

    //if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
    //auto baseCut = rhIso && hcalTowIso && hadOverE;
    //auto looseCut = baseCut && hcTrkIsoL;
    //auto tightCut = baseCut && hcTrkIsoT;
    //auto tightEB = tightCut && sigmaIeieEB;
    //auto tightEE = tightCut && sigmaIeieEE;
	bool isocut = rhIso && hcTrkIso && hadOverE;

	auto phoClass = isocut ? 1 : 0;
    //auto phoClass = tightCut ? 3 : looseCut ? 2 : baseCut ? 1 : 0;

	return phoClass;

}//<<>>int KUCMSAodSkimmer::getPhoQuality( int iter )

int KUCMSAodSkimmer::getJetQuality( int it ){

    const auto eta  = std::abs((*Jet_eta)[it]);	     
    const auto NHF  = (*Jet_neHEF)[it];
    const auto NEMF = (*Jet_neEmEF)[it];
    const auto CHF  = (*Jet_chHEF)[it];
    const auto CEMF = (*Jet_chEmEF)[it];
    const auto NHM  = (*Jet_neHM)[it];
    const auto CHM  = (*Jet_chHM)[it];
    const auto SHM  = NHM + CHM;
    const auto MUF  = (*Jet_muEF)[it];

    int tighter = 3;
    int tightLepVeto = 0;
    int tight = 2;
    int loose = 1;

    bool nhfup  = NHF  <= 0.90; 	// delpho
    bool nhflw  = NHF  >= 0.2;

    bool nemfup1 = NEMF <= 0.90; // delpho
    bool nemfup2 = NEMF <= 0.99;
    bool nemf80 = NEMF <= 0.80;
    bool nemflw = NEMF >= 0.01;
    bool nemf10 = NEMF >= 0.1;
	
    bool shm1  = SHM  >= 1;
    bool muf8  = MUF  <= 0.80;
    bool chf0  = CHF  >= 0;		// delpho
    bool chf10  = CHF  >= 0.10;
    bool chm0  = CHM  >= 0;		// delpho
    bool cemf8 = CEMF >= 0.80;	
    bool nhm2  = NHM  >= 1;
    bool nhm10 = NHM  >= 10;

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

}//<<>>int KUCMSAodSkimmer::getJetQuality( int iter )

//------------------------------------------------------------------------------------------------------------
// set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree ){

	//fOutTree->Branch( "RunNumber", &RunNumber );
    selEvtVars.makeBranch( "dsKey", "DataSetKey", STR, "Key for source data set of event" );
    selEvtVars.makeBranch( "evtGenWgt", FLOAT );
    selEvtVars.makeBranch( "PVx", FLOAT );
    selEvtVars.makeBranch( "PVy", FLOAT );
    selEvtVars.makeBranch( "PVz", FLOAT );
    selEvtVars.attachBranches( fOutTree );

	//selMet.makeBranch( "Met", FLOAT );
	selMet.makeBranch( "met", "selCMet", FLOAT, "Magnitude of event Met corrected for OOT photons" );
    selMet.makeBranch( "metPx", "selCMetPx", FLOAT, "Magnitude of event MetPx corrected for OOT photons" );
    selMet.makeBranch( "metPy", "selCMetPy", FLOAT, "Magnitude of event MetPy corrected for OOT photons" );
	selMet.attachBranches( fOutTree );

	selGenPart.makeBranch( "genPartEnergy", VFLOAT );
    selGenPart.makeBranch( "genPartEta", VFLOAT );
    selGenPart.makeBranch( "genPartPhi", VFLOAT );
    selGenPart.makeBranch( "genPartPt", VFLOAT );
    selGenPart.makeBranch( "genPartPdgId", VUINT );
    selGenPart.makeBranch( "genPartSusId", VINT );
    selGenPart.makeBranch( "genCharge", VINT );   //!
    selGenPart.makeBranch( "genMass", VFLOAT );   //!   
    selGenPart.makeBranch( "genStatus", VBOOL );   //!
    selGenPart.makeBranch( "genVx", VFLOAT );   //!
    selGenPart.makeBranch( "genVy", VFLOAT );   //!
    selGenPart.makeBranch( "genVz", VFLOAT );   //!
    selGenPart.makeBranch( "genPx", VFLOAT );   //!
    selGenPart.makeBranch( "genPy", VFLOAT );   //!
    selGenPart.makeBranch( "genPz", VFLOAT );   //!

	selGenPart.attachBranches( fOutTree );

	selElectrons.makeBranch( "nElectrons", UINT );
    selElectrons.makeBranch( "epDrMatch", VFLOAT );
    selElectrons.makeBranch( "epSeedIdMatch", VFLOAT );
    selElectrons.makeBranch( "epSeedIdUnique", VFLOAT );
    selElectrons.makeBranch( "elePhoMinDr", VFLOAT );

    selElectrons.attachBranches( fOutTree );

    selPhotons.makeBranch( "nPhotons", UINT );
    selPhotons.makeBranch( "nSelPhotons", UINT ); 
    selPhotons.makeBranch( "selPhoEleVeto", UINT );
    selPhotons.makeBranch( "selPhoOOT", VBOOL );
    selPhotons.makeBranch( "selPhoSusyId", VFLOAT );
    selPhotons.makeBranch( "selPhoQuality", VINT ); 
    selPhotons.makeBranch( "selPhoTime", VFLOAT ); 

    selPhotons.makeBranch( "selPhoEnergy", VFLOAT );
    selPhotons.makeBranch( "selPhoEta", VFLOAT ); 
    selPhotons.makeBranch( "selPhoPhi", VFLOAT );     
	selPhotons.makeBranch( "selPhoPt", VFLOAT ); 
    selPhotons.makeBranch( "selPhoSCx", VFLOAT );
    selPhotons.makeBranch( "selPhoSCy", VFLOAT );
    selPhotons.makeBranch( "selPhoSCz", VFLOAT );

    selPhotons.makeBranch( "selPhoHcalTowerSumEtBcConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoTrkSumPtHollowConeDR03", VFLOAT );
    selPhotons.makeBranch( "selPhoTrkSumPtHollowConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoTrkSumPtSolidConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoPixelSeed", VBOOL );
    selPhotons.makeBranch( "selPhoEcalRHSumEtConeDR04", VFLOAT );
    selPhotons.makeBranch( "selPhoHadTowOverEM", VFLOAT );
    selPhotons.makeBranch( "selPhoHadOverEM", VFLOAT );

    selPhotons.makeBranch( "selPhoHcalTowerSumEtConeDR04", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoNTrkHollowConeDR04", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoNTrkSolidConeDR04", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoEcalPFClusterIso", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoHasConversionTracks", VBOOL );   //!
    selPhotons.makeBranch( "selPhoHcalPFClusterIso", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoHoe_PUcorr", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfChargedIso", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfChargedIsoPFPV", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfChargedIsoWorstVtx", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfPhoIso03", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfRelIso03_all_quadratic", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoPfRelIso03_chg_quadratic", VFLOAT );   //!

    selPhotons.makeBranch( "selPhoSMaj", VFLOAT );
    selPhotons.makeBranch( "selPhoSMin", VFLOAT );
    selPhotons.makeBranch( "selPhoSAlp", VFLOAT );
    selPhotons.makeBranch( "selPhoClstrRn", VFLOAT );
    selPhotons.makeBranch( "selPhoR9", VFLOAT );
    selPhotons.makeBranch( "selPhoNrh", VUINT );
    selPhotons.makeBranch( "selPhoCovEtaEta", VFLOAT );
    selPhotons.makeBranch( "selPhoCovEtaPhi", VFLOAT );
    selPhotons.makeBranch( "selPhoCovPhiPhi", VFLOAT );
    selPhotons.makeBranch( "selPhoSieie", VFLOAT );
    selPhotons.makeBranch( "selPhoSieip", VFLOAT );
    selPhotons.makeBranch( "selPhoSipip", VFLOAT );
    selPhotons.makeBranch( "selPhoEtaWidth", VFLOAT );
    selPhotons.makeBranch( "selPhoPhiWidth", VFLOAT );
    selPhotons.makeBranch( "selPhoS4", VFLOAT );
    selPhotons.makeBranch( "selPhoSigmaIEtaIEta", VFLOAT );   //!

    selPhotons.makeBranch( "selPhoGenPt", VFLOAT );
    selPhotons.makeBranch( "selPhoPhoIsoDr", VFLOAT );
    selPhotons.makeBranch( "selPhoGenIdx", VINT );

    selPhotons.makeBranch( "selPhoGenSigMomEnergy", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomEta", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomMass", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPhi", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPt", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPx", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPy", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomPz", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomVx", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomVy", VFLOAT );   //!
    selPhotons.makeBranch( "selPhoGenSigMomVz", VFLOAT );   //!

    selPhotons.attachBranches( fOutTree );

    //selJets.makeBranch( "JetHt", &JetHt );
    selJets.makeBranch( "nJets", UINT );
    selJets.makeBranch( "nSelJets", UINT );
    selJets.makeBranch( "selJetSusyId", VFLOAT ); 
    selJets.makeBranch( "selJetQuality", VINT ); 
    selJets.makeBranch( "selJetPt", VFLOAT ); 
    selJets.makeBranch( "selJetEnergy", VFLOAT );
    selJets.makeBranch( "selJetEta", VFLOAT ); 
    selJets.makeBranch( "selJetPhi", VFLOAT );
    selJets.makeBranch( "selJetTime", VFLOAT ); 
    selJets.makeBranch( "selJetMass", VFLOAT );

    selJets.makeBranch( "selJetArea", VFLOAT ); //*Jet_area)[it]; 
    selJets.makeBranch( "selJetChEmEF", VFLOAT ); //*Jet_chEmEF)[it]; 
    selJets.makeBranch( "selJetchHEF", VFLOAT ); //*Jet_chHEF)[it]; 
    selJets.makeBranch( "selJetChHM", VFLOAT ); //*Jet_chHM)[it]; 
    selJets.makeBranch( "selJetMuEF", VFLOAT ); //*Jet_muEF)[it]; 
    selJets.makeBranch( "selJetNeEmEF", VFLOAT ); //*Jet_neEmEF)[it]; 
    selJets.makeBranch( "selJetNeHEF", VFLOAT ); //*Jet_neHEF)[it]; 
    selJets.makeBranch( "selJetNeHM", VFLOAT ); //*Jet_neHM)[it]; 

    selJets.makeBranch( "selGenJetDpt", VFLOAT ); //*Jet_genDptMatch)[it]; 
    selJets.makeBranch( "selGenJetdr", VFLOAT ); //*Jet_genDrMatch)[it]; 
    selJets.makeBranch( "selGenJetEnergy", VFLOAT ); //*Jet_genEnergy)[it]; 
    selJets.makeBranch( "selGenJeteta", VFLOAT ); //*Jet_genEta)[it]; 
    selJets.makeBranch( "selGenJetImpAng", VFLOAT ); //*Jet_genImpactAngle)[it]; 
    selJets.makeBranch( "selJetLlpDp", VFLOAT ); //*Jet_genLlpDp)[it]; 
    selJets.makeBranch( "selJetLlpDr", VFLOAT ); //*Jet_genLlpDr)[it]; 
    selJets.makeBranch( "selGenJetPt", VFLOAT ); //*Jet_genPt)[it]; 
    selJets.makeBranch( "selGenJetTof", VFLOAT ); //*Jet_genTOF)[it]; 
    selJets.makeBranch( "selGenJetTime", VFLOAT ); // (*Jet_genTime)[it]; 
    selJets.makeBranch( "selGenJetLlpTime", VFLOAT ); //*Jet_genTimeLLP)[it]; 

    selJets.attachBranches( fOutTree );

    selRjrVars.makeBranch( "X1aMass", VFLOAT );
    selRjrVars.makeBranch( "X1aCosA", VFLOAT );
    selRjrVars.makeBranch( "X1bMass", VFLOAT );
    selRjrVars.makeBranch( "X1bCosA", VFLOAT );
    selRjrVars.makeBranch( "X2aMass", VFLOAT );
    selRjrVars.makeBranch( "X2aCosA", VFLOAT );
    selRjrVars.makeBranch( "X2bMass", VFLOAT );
    selRjrVars.makeBranch( "X2bCosA", VFLOAT );
    selRjrVars.makeBranch( "SMass", VFLOAT );
    selRjrVars.makeBranch( "SCosA", VFLOAT );
    selRjrVars.makeBranch( "nRjrPhotons", VINT );
    selRjrVars.makeBranch( "RjrMET", VFLOAT );

    selRjrVars.attachBranches( fOutTree );

}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

void KUCMSAodSkimmer::endJobs(){ 

    std::cout << "Running EndJobs." << std::endl;

}//void KUCMSAodSkimmer::endJobs()

void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree ){

    std::cout << "Filling ConfigTree." << std::endl;
	TBranch *nEventBranch = fConfigTree->Branch( "nEvents", &nEvents );
	nEventBranch->Fill();
	TBranch *nSelectedEventsBranch = fConfigTree->Branch( "nSelectedEvents", &nSelectedEvents );
	nSelectedEventsBranch->Fill();
    TBranch *sKeyBranch = fConfigTree->Branch( "sKey", &dataSetKey ); 
	sKeyBranch->Fill();
    TBranch *sCrossSectionBranch = fConfigTree->Branch( "sCrossSection", &xsctn );
	sCrossSectionBranch->Fill();
    TBranch *sLambdaBranch = fConfigTree->Branch( "sGMSBGravMass", &gmass );
    sLambdaBranch->Fill();
    TBranch *sCTauBranch = fConfigTree->Branch( "sGMSBChi1Mass", &xmass );
    sCTauBranch->Fill();
    TBranch *sMCWgtBranch = fConfigTree->Branch( "sMCWgt", &mcwgt );
    sMCWgtBranch->Fill();
    TBranch *sMCTypeBranch = fConfigTree->Branch( "sMCType", &mctype );
    sMCTypeBranch->Fill();
    TBranch *sumEvtGenWgtBranch = fConfigTree->Branch( "sumEvtGenWgt", &sumEvtGenWgt );
    sumEvtGenWgtBranch->Fill();

	fConfigTree->Fill();

}//<<>>void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree, std::string key )

void KUCMSAodSkimmer::initHists(){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

	hist1d[0] = new TH1D("ecalrhenergy0", "RecHit Energy;rechit E [GeV]",200,0,0.2);
    hist1d[1] = new TH1D("ecalrhenergy1", "RecHit Energy;rechit E [GeV]",200,0,20);
    hist1d[2] = new TH1D("ecalrhenergy2", "RecHit Energy;rechit E [GeV]",200,0,2);

    hist1d[3] = new TH1D("inphorhe", "inPhoRecHit;rechit E [GeV]",200,0,2);
    hist1d[4] = new TH1D("inelerhe", "inEleRecHit;rechit E [GeV]",200,0,2);
    hist1d[5] = new TH1D("sbinphorhe", "In dr width of PhoSC;rechit E [GeV]",800,0,200);
    hist1d[6] = new TH1D("sbinelerhe", "In dr width of EleSC;rechit E [GeV]",800,0,200);
    hist1d[7] = new TH1D("notinrhe", "Not in Ele/Pho SC;rechit E [GeV]",200,0,2);
    hist1d[8] = new TH1D("inootphorhe", "inOOTPhoRecHit;rechit E [GeV]",200,0,2);

    hist1d[10] = new TH1D("ecalrhenergy4", "RecHit Energy;rechit E [GeV]",200,0,200);

    hist1d[13] = new TH1D("inphorhe2", "inPhoRecHit;rechit E [GeV]",200,0,2);
    hist1d[12] = new TH1D("inphorhecnt", "inPhoRecHit;rechit cnt",200,0,200);

    hist1d[14] = new TH1D("inelerhe2", "inEleRecHit;rechit E [GeV]",200,0,2);
    hist1d[15] = new TH1D("inelerhecnt", "inEleRecHit;rechit cnt",200,0,200);

    hist1d[17] = new TH1D("notinrhe2", "Not in Ele/Pho SC;rechit E [GeV]",200,0,2);

	hist1d[100] = new TH1D("hfdiff", "Photon H&F Diff;Diff( rh found - rh in H&F )",100,-50,50);
    ////hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    ////hist2d[1] = new TH2D("jetDrMuTime_pt", "jetDrMuTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);

	hist2d[0] = new TH2D("prhevoot", "prhe v isOOT", 200, 0, 2, 3, -1, 2);
    hist2d[1] = new TH2D("prhevswcr", "prhe v swissCross", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[2] = new TH2D("prheveta", "prhe v Eta", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[3] = new TH2D("prhevamp", "prhe v amplitude", 200, 0, 2, 150, 0, 150);
    hist2d[4] = new TH2D("prhevtime", "prhe v time", 200, 0, 1, 150, -15, 15);

    hist2d[10] = new TH2D("erhevoot", "erhe v isOOT", 200, 0, 2, 3, -1, 2);
    hist2d[11] = new TH2D("erhevswcr", "erhe v swissCross", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[12] = new TH2D("erheveta", "erhe v Eta", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[13] = new TH2D("erhevamp", "erhe v amplitude", 200, 0, 2, 150, 0, 150);
    hist2d[14] = new TH2D("erhevtime", "erhe v time", 200, 0, 1, 150, -15, 15);

    hist2d[15] = new TH2D("ootprhevoot", "ootprhe v isOOT", 200, 0, 2, 3, -1, 2);
    hist2d[16] = new TH2D("ootprhevswcr", "ootprhe v swissCross", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[17] = new TH2D("ootprheveta", "ootprhe v Eta", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[18] = new TH2D("ootprhevamp", "ootprhe v amplitude", 200, 0, 2, 150, 0, 150);
    hist2d[19] = new TH2D("ootprhevtime", "ootprhe v time", 200, 0, 1, 150, -15, 15);

    hist2d[20] = new TH2D("urhevoot", "urhe v isOOT", 200, 0, 2, 3, -1, 2);
    hist2d[21] = new TH2D("urhevswcr", "urhe v swissCross", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[22] = new TH2D("urheveta", "urhe v Eta", 200, 0, 2, 150, -1.5, 1.5);
    hist2d[23] = new TH2D("urhevamp", "urhe v amplitude", 200, 0, 2, 150, 0, 150);
    hist2d[24] = new TH2D("urhevtime", "urhe v time", 200, 0, 1, 150, -15, 15);
    
	//------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

    //hist3d[0] = new TH3D("phoNRH_ClR9_phoID","Photon nClRecHits v clR9 v phoId;nRecHits;ClusterR9;PhotonID(fake1,loose2,tight3)",200,0,200,100,0,1,8,0,4);

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------

    if( CLSTRMAPS ){
		nMaps = 0;
		for(int it=0; it<nEBEEMaps; it++){
			fMap[it] = false;
			std::string label(";iEta;iPhi");
    	    std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
    	    ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
    	    std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
    	    ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
			std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
    	    ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
		}//<<>>for(int it=0; it<nEBEEMaps; it++)
	}//<<>>if( CLSTRMAPS )

}//<<>>void KUCMSAodSkimmer::initHists()

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
//	Hlper functions that depended on varibles from GetEntry()
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

int KUCMSAodSkimmer::getRhIdx( uInt rhDetID ){

    for( int idx = 0; idx < ECALRecHit_ID->size(); idx++ ){ if( rhDetID == (*ECALRecHit_ID)[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int KUCMSAodSkimmer::getRhIdx( int rhDetID )

uInt KUCMSAodSkimmer::getLeadRhID( std::vector<uInt> recHitIds ){

    uInt result;
    float enr(0.0);
    for( auto id : recHitIds ){
        auto rhenr = (*ECALRecHit_energy)[getRhIdx(id)];
        if( rhenr > enr ){ enr = rhenr; result = id; }
    }//<<>>for (const auto recHit : recHits )

    return result;

}//>>>>EcalRecHit KUCMSAodSkimmer::getLeadRh( rhGroup recHitsi

float KUCMSAodSkimmer::clstrR9( std::vector<uInt> recHitIds ){

    auto leadRhID = getLeadRhID( recHitIds );
    auto leadRhEn = (*ECALRecHit_energy)[getRhIdx(leadRhID)];
    float sumRhEn(0);
    for ( auto id : recHitIds ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
    return sumRhEn > 0 ? leadRhEn/sumRhEn  : - 1.0;

}//<<>>float KUCMSAodSkimmer::clstrR9( vector<uInt> recHitIds )

std::vector<float> KUCMSAodSkimmer::getRhGrpEnergies( std::vector<uInt> rechitids ){

	std::vector<float> result;
	for ( auto id : rechitids ){ result.push_back((*ECALRecHit_energy)[getRhIdx(id)]); }
	return result;

};//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpEnergies( std::vector<uInt> rechitids )

std::vector<float> KUCMSAodSkimmer::getRhGrpTimes( std::vector<uInt> rechitids ){

    std::vector<float> result;
    for ( auto id : rechitids ){
		auto rhtime = (*ECALRecHit_time)[getRhIdx(id)];
		auto rh0tof = (*ECALRecHit_0TOF)[getRhIdx(id)];
        auto rhpvtof = (*ECALRecHit_pvTOF)[getRhIdx(id)];
		result.push_back(rhtime+rh0tof-rhpvtof); 
	}//<<>>for ( auto id : recHitIds )
    return result;

};//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpTimes( std::vector<uInt> rechitids )

std::vector<float> KUCMSAodSkimmer::getRhGrpEigenFromAngles( std::vector<uInt> rechitids ){

	//bool verbose = true;
    bool verbose = false;

	if( DEBUG || verbose ) std::cout << " Starting getRhGrpEigen_sph  " << std::endl;

    std::vector<float> emptyReturn(9,-9);
    std::vector<float> egwts;
    std::vector<float> rhetas, rhphis;
    std::vector<float> logwtvec, tresvec;
    auto nRecHits = rechitids.size();
    if( nRecHits < 5 ){ if( verbose ) std::cout << " ----  rechit collection has too few rechits" << std::endl; return emptyReturn; }
    float sumRhEn(0);
    for ( auto id : rechitids ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
	if( verbose ) std::cout << " --- EigenAngles sumRhEn : " << sumRhEn << std::endl;
    if( sumRhEn <= 0 ){ if( verbose ) std::cout << " ----  rechit collection has no energy" << std::endl; return emptyReturn; }
    if( DEBUG ) std::cout << "1a, ";
    std::cout << " getRhGrpEigenFromAngles : nRechits: " << nRecHits << " -------------------------- " << std::endl;
    for( uInt it(0); it < nRecHits; it++ ){

        const auto rhIDX = getRhIdx(rechitids[it]);
        auto idinfo = DetIDMap[rechitids[it]];
        auto isEB = idinfo.ecal == ECAL::EB;
        if( rhIDX == -1 ){ if( verbose ) std::cout << " ---- Bad idx !!!!! -- In getRhGrpEigen ---- " << std::endl; return emptyReturn; }
        if( not isEB ){ if( verbose ) std::cout << " ---- rechit group has EE members " << idinfo.ecal << std::endl; return emptyReturn; }

        std::cout << " rhIDX: " << rhIDX << std::endl;
    	const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
        std::cout << " rhEtaPos: " << rhEtaPos << std::endl;
        rhetas.push_back((rhEtaPos>0)?rhEtaPos+84.5:rhEtaPos+85.5);
        const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
        std::cout << " rhPhiPos: " << rhPhiPos << std::endl;
        rhphis.push_back(rhPhiPos-0.5);
        auto rhenergy = (*ECALRecHit_energy)[rhIDX];
        auto logwt = std::max(0.0, 4.2 + log(rhenergy/sumRhEn));// cut at rh energy < 1.5% of cluster
        logwtvec.push_back(logwt);

    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    std::vector<float> detas, dphis, angles, redlogwtvec;
    auto meta = mean( rhetas, logwtvec );
    auto mphi = meanIPhi( rhphis, logwtvec );
    for( uInt it(0); it < rhetas.size(); it++ ){

        float deta = rhetas[it]-meta;
        float dphi = dIPhi( rhphis[it], mphi );
		if( dphi < 0.5 && deta < 0.5 ) continue;
        detas.push_back(deta);
        dphis.push_back(dphi);
        float angle = getAngle( deta, dphi );
        angles.push_back(angle);
		redlogwtvec.push_back(logwtvec[it]);

    }//<<>>for( uInt it(0); it < etas.size(); it++ )

    auto eigens = getRhGrpEigen( angles, redlogwtvec );
	if( verbose ) std::cout << " --- EigenAngles VAlue : " << eigens[2] << std::endl;

    auto phiCorrFactor = 0.8;
    std::cout << " mean eta: " << meta << " mean phi: " << mphi << std::endl; 
    auto sxx = var( detas, meta, redlogwtvec );
    auto syy = var( dphis, mphi, redlogwtvec, accum(redlogwtvec)/phiCorrFactor );
    auto sxy = cvar( detas, meta, dphis, mphi, redlogwtvec, accum(redlogwtvec)/std::sqrt(phiCorrFactor) );
    std::cout << " sxx: " << sxx << " syy: " << syy << " sxy: " << sxy << std::endl;
    auto smaj = (sxx + syy + std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto smin = (sxx + syy - std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto sang = std::atan((sxx-syy+std::sqrt(sq2(syy-sxx)+4.*sq2(sxy)))/(2.*sxy));

    //eigens[0] //0 geoeigan x vec
    //eigens[1] //1 geoeigan y vec
    //eigens[2] //2 geoeigan mag vec
    eigens.push_back(smaj);//3
    eigens.push_back(smin);//4
    eigens.push_back(sang);//5

    if( DEBUG ) std::cout << " Done" << std::endl;;
    return eigens;

}//<<>>std::vector<float> KUCMSAodSkimmer::getRhGrpEigen( std::vector<uInt> rechitids )

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
////!!!!!!!!!!!!!!!!!11  need rhPosX, rhPosY, & rhPosZ to be saved for following functions   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

std::vector<float> KUCMSAodSkimmer::getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ ){

    std::vector<float> result;
    if( recHitIds.size() < 1 ){ result.push_back(-99); return result; }
    auto lrhid = getLeadRhID(recHitIds);
    auto lrhidx = getRhIdx(lrhid);
    auto X = (*ECALRecHit_rhx)[lrhidx];
    auto Y = (*ECALRecHit_rhy)[lrhidx];
    auto Z = (*ECALRecHit_rhz)[lrhidx];
    const auto d_rh = hypo( X, Y, Z);
    const auto d_pv = hypo( X-vtxX, Y-vtxY, Z-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    for( int idx = 0; idx < ECALRecHit_time->size(); idx++ ){result.push_back((*ECALRecHit_time)[idx]-tof);}
    return result;

}//>>>>vector<float> KUCMSAodSkimmer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

std::vector<float> KUCMSAodSkimmer::getRhGrpEigen_sph( std::vector<float> times, std::vector<uInt> rechitids ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2
	if( DEBUG ) std::cout << " Starting getRhGrpEigen_sph  " << std::endl;

    float N(3.64);
    float C(0.3000);

    std::vector<float> eg2wts;
    std::vector<float> rhetas, rhphis;;
    std::vector<float> xs, ys, zs;
    std::vector<float> rhxs, rhys, rhzs;
    std::vector<float> rhtimes;
    std::vector<float> grhtimes;
    std::vector<float> angles;
    std::vector<float> zcAngles;
    std::vector<float> invtresvec, rhinvtresvec;
    std::vector<float> rhlogwtvec, rhtresvec;
    std::vector<float> logwtvec, tresvec;
    std::vector<float> emptyReturn(40,-9);

    // --------- prepare inputs for eigan calcs --------------------------
    if( DEBUG ) std::cout << " getRhGrpEigen_sph 1, ";

    auto nRecHits = rechitids.size();
    if( nRecHits < 16 ) return emptyReturn;
    float sumRhEn(0);
    for ( auto id : rechitids ){ sumRhEn +=  (*ECALRecHit_energy)[getRhIdx(id)]; }
    if( sumRhEn <= 0 ) return emptyReturn;
    if( DEBUG ) std::cout << "1a, ";
    for( uInt it(0); it < nRecHits; it++ ){

        const auto rhIDX = getRhIdx(rechitids[it]);
        auto idinfo = DetIDMap[rechitids[it]];
        auto isEB = idinfo.ecal == ECAL::EB;
        //if( isEB ) hist1d[123]->Fill(times[it]);
        //if( DEBUG ) std::cout << "In getRhGrpEigen_sph w/ idx : " << rhIDX << std::endl;
        if( rhIDX == -1 ){ return emptyReturn; std::cout << " -- Bad idx !!!!! -- In getRhGrpEigen_sph ---- " << std::endl; }
        if( isEB ){
            const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
            rhetas.push_back((rhEtaPos>0)?rhEtaPos+84.5:rhEtaPos+85.5);
            const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
            rhphis.push_back(rhPhiPos-0.5);
            const auto rhXPos = (*ECALRecHit_rhx)[rhIDX];
            rhxs.push_back(rhXPos);
            const auto rhYPos = (*ECALRecHit_rhy)[rhIDX];
            rhys.push_back(rhYPos);
            const auto rhZPos = (*ECALRecHit_rhz)[rhIDX];
            rhzs.push_back(rhZPos);
            rhtimes.push_back(times[it]);
            auto rhenergy = (*ECALRecHit_energy)[rhIDX];
            auto resolution = sq2(N/rhenergy)+2*C*C;
            auto logwt = std::max(0.0, 4.2 + log(rhenergy/sumRhEn));// cut at rh energy < 1.5% of cluster
            rhlogwtvec.push_back(logwt);
            rhinvtresvec.push_back(1/resolution);
            rhtresvec.push_back(resolution);
            //if( DEBUG ) std::cout << "In getRhGrpEigen_sph w/ rheta " << rhEtaPos << " : rhphi " << rhPhiPos; 
            //if( DEBUG ) std::cout << " : rht " << times[it] << std::endl;

        } else {
            return emptyReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
        }//<<>> else - if( idinfo.ecal == ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )
	
    if( rhtimes.size() < 9 ) return emptyReturn;

    bool uselog(true);

    auto rhewgt = (uselog)?rhlogwtvec:rhinvtresvec;
    auto tmtime = mean(rhtimes,rhewgt);
    for( uInt it(0); it < rhtimes.size(); it++ ){
        auto rht = rhtimes[it];
        auto goodDifTime = (rhtimes[it]-tmtime) < 10.0;
        auto isInTime = ( rhtimes[it] > -50.0 ) && ( rhtimes[it] < 50.0 );
        if( isInTime && goodDifTime ){
            grhtimes.push_back(rhtimes[it]);
            xs.push_back(rhxs[it]);
            ys.push_back(rhys[it]);
            zs.push_back(rhzs[it]);
            logwtvec.push_back(rhlogwtvec[it]);
            tresvec.push_back(rhtresvec[it]);
            invtresvec.push_back(rhinvtresvec[it]);
        }//<<>>if( isInTime && goodDifTime )
    }//<<>>for( auto rht : rhtimes )
    if( grhtimes.size() < 9 ) return emptyReturn;

    if( DEBUG ) std::cout << "2, ";
    std::vector<float> letas;
    std::vector<float> lphis;
    auto meta = mean(rhetas,rhewgt);
    auto mphi = meanIPhi(rhphis,rhewgt);
    //auto meta = mean(rhetas);
    //auto mphi = meanIPhi(rhphis);
    for( uInt it(0); it < rhetas.size(); it++ ){
        float leta = rhetas[it]-meta;
        letas.push_back(leta);
        float lphi = dIPhi(rhphis[it],mphi);
        lphis.push_back(lphi);
        float angle = getAngle( leta, lphi );
        angles.push_back(angle);
    }//<<>>for( uInt it(0); it < etas.size(); it++ )

    auto ewgt = (uselog)?logwtvec:invtresvec;
    auto mtime = mean(grhtimes,ewgt);
    auto mx = mean(xs,ewgt);
    auto my = mean(ys,ewgt);
    auto mz = mean(zs,ewgt);
    auto mr = hypo(mx,my);
    auto ma = std::atan2(my,mx);

    //if( DEBUG ) std::cout << "In getRhGrpEigen_sph w/ meta " << meta << " : mphi " << mphi << " : mt " << mtime << std::endl;
    std::vector<float> lzs;
    std::vector<float> lcs;
    std::vector<float> lts;
    std::vector<float> ltds;
    std::vector<float> nolts;
    std::vector<float> invltres;
    float minDr(100.0);
    float vslz(0.0);
    float vslc(0.0);
    float vs3lz(0.0);
    float vs3lc(0.0);
    float vs3lt(0.0);
    for( uInt it(0); it < grhtimes.size(); it++ ){

        float ltim = grhtimes[it]-mtime;
        lts.push_back(ltim);
        auto ltd = ltim*SOL;
        ltds.push_back(ltd);
        invltres.push_back(ewgt[it]);
        auto lt2reswt = ltim*ltim*ewgt[it];
        eg2wts.push_back(lt2reswt);
        nolts.push_back(0.0);

        float lz = zs[it]-mz;
        lzs.push_back(lz);
        float lc = mr*(std::atan2(ys[it],xs[it])-ma);
        lcs.push_back(lc);
        //float zcAngle = getAngle( lz, lc );
        float zcAngle = std::atan2(lc,lz);
        zcAngles.push_back(zcAngle);

        auto dr = hypo(lz,lc);
        auto drt = hypo(ltd,dr);
        if( dr < minDr ) minDr = dr;
        // do dr calc in ieta,iphi cross check

        auto sqrtreswt = (uselog)?ewgt[it]:std::sqrt(ewgt[it]);
        auto ltdsqrtreswt = std::abs(ltd)*sqrtreswt;
        vs3lz += lz*ltdsqrtreswt/drt;
        vs3lc += lc*ltdsqrtreswt/drt;
        vs3lt += ltd*sqrtreswt/drt;
        vslz += lz*ltdsqrtreswt/dr;
        vslc += lc*ltdsqrtreswt/dr;
        //if( (std::abs(lz) < minDl) && (std::abs(lc) < minDl) ) ewgt[it] = 0;

    }//<<>>for( uInt it(0); it < grhtimes.size(); it++ )

    // --------------  get eigan values and vectors ---------------------------------------
    if( DEBUG ) std::cout << "3, ";

    auto eigens =  getRhGrpEigen( zcAngles, eg2wts );//0 x, 1 y, 2 values
    auto d2dot = eigens[0]*vslz + eigens[1]*vslc;
    if( d2dot < 0 ){ eigens[0] *= -1; eigens[1] *= -1; }
    auto eigens3 =  getRhGrpEigen( lzs, lcs, ltds, invltres );//0 x, 1 y, 2 values
    auto d3dot = eigens3[0]*vs3lz + eigens3[1]*vs3lc + eigens3[2]*vs3lt;
    if( d3dot < 0 ){ eigens3[0] *= -1; eigens3[1] *= -1; eigens3[2] *= -1; }
    auto geoeigens =  getRhGrpEigen( zcAngles, invltres );//0 x, 1 y, 2 values
    auto geoddot = geoeigens[0]*vslz + geoeigens[1]*vslc;
    if( geoddot < 0 ){ geoeigens[0] *= -1; geoeigens[1] *= -1; }
    ////auto geoeigens3 =  getRhGrpEigen( lzs, lcs, nolts, invltres );//0 x, 1 y, 2 values

    // --------------  get eigan vector angles ------------------------------------------ 
    if( DEBUG ) std::cout << "4, ";

    float rotangle = getAngle(eigens[0], eigens[1]);
    float e2sin = std::sin(rotangle); //eigens[1];
    float e2cos = std::cos(rotangle); //eigens[0];
    //float rot3angle = getAngle(eigens3[0], eigens3[1]);
    float rot3angle = getAngle(eigens3[0], eigens3[1]);
    float e3sin = std::sin(rot3angle);
    float e3cos = std::cos(rot3angle);
    float rotgangle = getAngle(geoeigens[0], geoeigens[1]);
    float egsin = std::sin(rotgangle);
    float egcos = std::cos(rotgangle);

    // -----------------------------------------
    // finding nemo ( slope )
    // -----------------------------------------
    if( DEBUG ) std::cout << "6, ";

    auto nWts = invltres.size();

    std::vector<float> xs1;
    std::vector<float> xs3;
    std::vector<float> xsgeo;
    std::vector<float> slvars1;
    std::vector<float> slvars3;
    std::vector<float> slvarsgeo;
    float xsum1(0.0);
    float xsum3(0.0);
    float xsumgeo(0.0);

    auto dsxcor = e2cos*(2.2) - e2sin*(2.2);
    auto d3xcor = e3cos*(2.2) - e3sin*(2.2);
    auto dgxcor = (geoeigens[0])*(2.2) - (geoeigens[1])*(2.2);
    auto xscorvar = sq2(dsxcor)/12;
    auto x3corvar = sq2(d3xcor)/12;
    auto xgcorvar = sq2(dgxcor)/12;

    // for pairs method
    std::vector<float> plzs;
    std::vector<float> plcs;
    //vector<float> plws;
    std::vector<float> pl3zs;
    std::vector<float> pl3cs;
    std::vector<float> plgzs;
    std::vector<float> plgcs;

    if( DEBUG ) std::cout << "7, ";
    for( uInt it(0); it < lzs.size(); it++ ){

        auto xscor = e2cos*lzs[it] - e2sin*lcs[it];
        auto yscor = e2sin*lzs[it] + e2cos*lcs[it];
        auto x3cor = e3cos*lzs[it] - e3sin*lcs[it];
        auto y3cor = e3sin*lzs[it] + e3cos*lcs[it];
        auto xgcor = egcos*lzs[it] - egsin*lcs[it];
        auto ygcor = egsin*lzs[it] + egcos*lcs[it];

        // for pairs method && histogram maps
        plzs.push_back(xscor);
        plcs.push_back(yscor);
        pl3zs.push_back(x3cor);
        pl3cs.push_back(y3cor);
        plgzs.push_back(xgcor);
        plgcs.push_back(ygcor);

        //if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << letas[it] << " : lphi " << lphis[it]
        //                        << " : xsor " << xscor << " : ycor " << yscor << " : dt " << eg2wts[it] << std::endl;
        //if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << letas[it] << " : lphi " << lphis[it]
        //                        << " : sxcor " << x3cor << " : sycor " << y3cor << " : dt " << eg2wts[it] << std::endl;

        // calc slope info
        auto sl1 = (lts[it])/(xscor);//*slopeCorr;
        auto sl3 = (lts[it])/(x3cor);//*slopeCorr;
        auto slg = (lts[it])/(xgcor);//*slopeCorr;
        xs1.push_back(sl1);
        xs3.push_back(sl3);
        xsgeo.push_back(slg);
        //slvars1.push_back(1/((notresvec[it]+tottresvec+sq2(sl1)*xscorvar*(1.0+(1.0/nWts)))/sq2(xscor)));
        slvars1.push_back(invltres[it]);
        //slvars3.push_back(1/((notresvec[it]+tottresvec+sq2(sl3)*x3corvar*(1.0+(1.0/nWts)))/sq2(x3cor)));
        slvars3.push_back(invltres[it]);
        //slvarsgeo.push_back(1/((notresvec[it]+tottresvec+sq2(slg)*xgcorvar*(1.0+(1.0/nWts)))/sq2(xgcor)));
        slvarsgeo.push_back(invltres[it]);
        xsum1 += sl1*slvars1[it];
        xsum3 += sl3*slvars3[it];
        xsumgeo += slg*slvarsgeo[it];

    }//<<>>for( uInt it(0); it < wts.size(); it++ )

//===================================================================
    if( DEBUG ) std::cout << "8, ";

    std::vector<float> pszs1;
    std::vector<float> plsrs;
    float plssum(0.0);
    std::vector<float> p3zs1;
    std::vector<float> pl3rs;
    float pl3sum(0.0);
    std::vector<float> pgzs1;
    std::vector<float> plgrs;
    float plgsum(0.0);

    //=============================================================
    // pairs method set up
    //============================================================

    for( uInt it1(0); it1 < plzs.size(); it1++ ){
        for( uInt it2(it1); it2 < plzs.size(); it2++ ){

            auto minDz = 3.3;
            auto plz = plzs[it1]-plzs[it2];
            if( std::abs(plz) < minDz ) plz = 999;
            auto pl3z = pl3zs[it1]-pl3zs[it2];
            if( std::abs(pl3z) < minDz ) pl3z = 999;
            auto plgz = plgzs[it1]-plgzs[it2];
            if( std::abs(plgz) < minDz ) plgz = 999;
            //auto plc = plcs[it1]-plcs[it2]; 
            auto pltime = grhtimes[it1]-grhtimes[it2];

            //auto gaplr = std::sqrt(invltres[it1]*invltres[it2]); 
            auto gaplr = std::sqrt(std::sqrt(invltres[it1])*std::sqrt(invltres[it2]));

            if( plz != 999 ){
                auto psl = pltime/plz;
                pszs1.push_back(psl);
                //auto gaplr = (slvars1[it1] + slvars1[it2]);
                plsrs.push_back( gaplr );
                plssum += psl*gaplr;
            }//<<>>if( plz != 999 )

            if( pl3z != 999 ){
                auto p3sl = pltime/pl3z;
                p3zs1.push_back(p3sl);
                //auto gaplr = (slvars3[it1] + slvars3[it2]);
                pl3rs.push_back( gaplr );
                pl3sum += p3sl*gaplr;
            }//<<>>if( plz != 999 )

            if( plgz != 999 ){
                auto pgsl = pltime/plgz;
                pgzs1.push_back(pgsl);
                //auto gaplr = (slvarsgeo[it1] + slvarsgeo[it2]);
                plgrs.push_back( gaplr );
                plgsum += pgsl*gaplr;
            }//<<>>if( plz != 999 )

        }//<<>>for( uInt it2(it1); it2 < grhtimes.size(); it2++ )
    }//<<>>for( uInt it1(0); it1 < grhtimes.size(); it1++ )

//--------------------------------------------------------------------
    //std::cout << "9, ";


    //find eigan vecgtor aligment
    //float eigensOld2d0(eigens[0]), eigensOld2d1(eigens[1]);
    //if( xsum1 < 0 ){ eigensOld2d0 *= -1; eigensOld2d1 *= -1; xsum1 = std::abs(xsum1);}
    //if( plssum < 0 ){ eigens[0] *= -1; eigens[1] *= -1; plssum = std::abs(plssum);}
    //if( pl3sum < 0 ){ eigens3[0] *= -1; eigens3[1] *= -1; eigens3[2] *= -1; pl3sum = std::abs(pl3sum);}
    //if( plgsum < 0 ){ geoeigens[0] *= -1; geoeigens[1] *= -1; plgsum = std::abs(plgsum);}


    // -----------   compute final outputs --------------------------------

    auto phiCorrFactor = 0.8;
    auto sxx = var( letas, 0., rhlogwtvec );
    auto syy = var( lphis, 0., rhlogwtvec, accum(rhlogwtvec)/phiCorrFactor );
    auto sxy = cvar( letas, 0., lphis, 0., rhlogwtvec, accum(rhlogwtvec)/std::sqrt(phiCorrFactor) );
    auto smaj = (sxx + syy + std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto smin = (sxx + syy - std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto sang = std::atan((sxx-syy+std::sqrt(sq2(syy-sxx)+4.*sq2(sxy)))/(2.*sxy));

    auto ebside = ( mz > 0 ) ? 1 : -1;
    auto taflip = ( ((std::abs(mz) < std::abs(PV_z)) && (mz*PV_z > 0) ) ? -1 : 1 )*ebside;

    //2d egians slope
    auto nXSum = xs1.size();
    auto totSloRes1 = accum(slvars1);
    auto slope1 = xsum1/totSloRes1;
    auto slope1err = std::sqrt(1/totSloRes1);
    auto varsl = var(xs1,slope1,slvars1,totSloRes1);
    auto chi1 = chisq(xs1,slope1,varsl);
    auto slchi2v1 = 0.f; //= chisqv(xs1,slope1,slvars1,varsl);//?????????????????
    auto chi2pf1 = 1 - TMath::Prob(chi1, nWts);

    // 3d eigan slope
    auto totSloRes3 = accum(pl3rs);
    auto slope3 = pl3sum/totSloRes3;
    auto slope3err = std::sqrt(1/totSloRes3);

    //geo eigan slope
    auto totSloResGeo = accum(plgrs);
    auto slopeg = plgsum/totSloResGeo;
    auto slopegerr = std::sqrt(1/totSloResGeo);

    //pairs slope
    auto totSloResPrs = accum(plsrs);
    auto slopeprs = plssum/totSloResPrs;
    auto slopeprserr = std::sqrt(1/totSloResPrs);

    // eigan3d angle slope
    auto angle3d = getAngle(eigens3[1],eigens3[2]);
    auto hypo3d02 = hypo(eigens3[0],eigens3[2]);
    auto hypo3d01 = hypo(eigens3[0],eigens3[1]);
    auto hypo3d12 = hypo(eigens3[1],eigens3[2]);
    auto slope3d = std::atan(eigens3[2]/hypo3d01);
    //auto slope3d = 100*(eigens3[0]/hypo3d12)/SOL;

/*///////////////////////////////////////////////////////////////////////
    // Fill Histograms
    for( uInt it(0); it < lzs.size(); it++ ){

        //plzs.push_back(xscor);
        //plcs.push_back(yscor);
        //pl3zs.push_back(x3cor);
        //pl3cs.push_back(y3cor);
        //plgzs.push_back(xgcor);
        //plgcs.push_back(ygcor);
        if( invltres[it] == 0 ) continue;
        auto sqrtwt = std::sqrt(invltres[it]);
        auto wichsum = slopeprs;
        //auto wichsum = pl3sum;
        //auto wichsum = xsum1;
        //auto wichsum = 1;
        auto xcor = ( wichsum < 0 ) ? -1*plzs[it] : plzs[it];
        auto ycor = ( wichsum < 0 ) ? -1*plcs[it] : plcs[it];
        auto fill = lts[it]*sqrtwt;
        hist2d[206]->Fill(xcor,ycor,fill);
        hist2d[207]->Fill(xcor,ycor,sqrtwt);
        hist2d[208]->Fill(lzs[it],lcs[it],fill);
        hist2d[209]->Fill(lzs[it],lcs[it],sqrtwt);
        hist2d[214]->Fill(ycor,lts[it],sqrtwt);
        hist2d[215]->Fill(xcor,lts[it],sqrtwt);

    }//<<>>for( uInt it(0); it < lzs.size(); it++ )
*/////////////////////////////////////////////////////////////////////

    // Fill results vector
    if( DEBUG ) std::cout << "10, ";
    // eigens 0 = vector x, 1 = vector y, 2 = vector mag
    eigens.push_back(slope1);//3  aligned slope
    eigens.push_back(chi2pf1);//4 aligned slope chi sqr prob
    eigens.push_back(slope3);//5 3d pairs slope
    eigens.push_back(0);//6 
    eigens.push_back(rotangle);//7 aligned rotation angle
    eigens.push_back(nXSum);//8 # of entries ( rechits )
    eigens.push_back(rot3angle);//9 3d rotation angle
    eigens.push_back(std::sqrt(varsl));//10 stdev aligned slope
    eigens.push_back(slope1err);//11 err aligned slope
    eigens.push_back(0);//12 
    eigens.push_back(slope3err);//13 errr unaligned slope
    eigens.push_back(slchi2v1);//14 chisqr like gof aligned slope
    eigens.push_back(minDr);//15 
    eigens.push_back(geoeigens[0]);//16 geoeigan x vec
    eigens.push_back(geoeigens[1]);//17 geoeigan y vec
    eigens.push_back(geoeigens[2]);//18 geoeigan mag vec
    eigens.push_back(ebside);//19 EB side
    eigens.push_back(taflip);//20 towards(+)/away(-) slope sign flip
    eigens.push_back(0);//21
    eigens.push_back(slopeg);//22 geoeigan slope
    eigens.push_back(slopegerr);//23 geoeigan slope error
    eigens.push_back(slopeprs);//24 pairs slope method
    eigens.push_back(slopeprserr);//25 pairs slope method error
    eigens.push_back(0);//26
    eigens.push_back(slope3d);//27 slope from 3d major eiganvector in "time" deminsion
    eigens.push_back(eigens3[3]);//28 "scaled magnitude" of 3d eiganvalue3
    eigens.push_back(angle3d);//29 rotation angle for eigan 3d
    eigens.push_back(eigens3[0]);//30 ev x  "z"
    eigens.push_back(eigens3[1]);//31 ev y  "c"
    eigens.push_back(eigens3[2]);//32 ev z  "time"
    eigens.push_back(eigens3[6]);//33 3dEV 6 = c vs c+t , ? 3d4Value
    eigens.push_back(eigens3[5]);//34 3d Time EV 5 = t vs t+zc oval ? 3dTvalue
//std::cout << "Slope egin : " << slope1 << " " << chi2pf1 << " " << rotangle << " " << std::sqrt(varsl) << " " << slope1err << std::endl;
    eigens.push_back(smaj);//35
    eigens.push_back(smin);//36
    eigens.push_back(sang);//37
    eigens.push_back(eigens3[8]);//38
    eigens.push_back(eigens3[9]);//39
    eigens.push_back(eigens3[10]);//40

    if( DEBUG ) std::cout << " Done" << std::endl;;
    return eigens;
}//>>>>vector<float> KUCMSAodSkimmer::getRhGrpEigen_sph( vector<float> times, rhGroup rechits ){

