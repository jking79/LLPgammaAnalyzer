//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSkimmer.hh"

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// KUCMSAodSkimmer class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

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

void KUCMSAodSkimmer::kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename ){

    const std::string disphotreename = "tree/llpgtree";
    //std::string inpath, infiles, key, 
	std::string masterstr; 
    //int mct;
	//float crossSection, gmsblam, gmsbct, mcw;
	std::cout << "Processing Input Lists for : " << infilelist << std::endl;
	std::ifstream masterInfile(listdir+infilelist);
	//while( masterInfile >> inpath >> infiles >> key >> crossSection >> gmsblam >> gmsbct >> mcwgt >> mctype ){
    while( std::getline( masterInfile, masterstr ) ){

		if( DEBUG ) std:: cout << masterstr << std::endl;
		if( masterstr == " " ) continue;
		auto instrs = splitString( masterstr, " " );
		if( DEBUG ) std:: cout << instrs.size() << std::endl;
        auto inpath = instrs[0];
        if( inpath == "#" || instrs.size() < 8 ) continue;

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
	
		KUCMSAodSkimmer::Init(fInTree);
		initHists();
	    setOutputBranches(fOutTree);	
	
		SetupDetIDsEB(DetIDMap);
		SetupDetIDsEE(DetIDMap);
	
		startJobs();
	
	    std::cout << "Setting up For Main Loop." << std::endl;
		int loopCounter(5000);
	    auto nEntries = fInTree->GetEntries();
	    if(DEBUG){ nEntries = 1000; loopCounter = 100; }
	    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
	    for (Long64_t centry = 0; centry < nEntries; centry++){

	        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
	        auto entry = fInTree->LoadTree(centry);
			getBranches(entry);
			geCnts.clear();
            geVars.clear();
			auto saveToTree = eventLoop(entry);
			if( saveToTree ) fOutTree->Fill();

	    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
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
	processEvntVars();	
	//processRechits();
	processMet();
	processPhotons();
	//processElectrons();
	//processMuons();
	processJets();

	// select events to process and store
	//--------------------------------------
	auto saveToTree = eventSelection();	
	if( saveToTree ) processRJR();
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

    selEvtVars.clearBranches();

	// calc
    nEvents++;
    
	//fill

    selEvtVars.fillBranch( "dsKey", dataSetKey );

}//<<>>void KUCMSAodSkimmer::processEvntVars()

void KUCMSAodSkimmer::processMet(){

	//intilize
	selMet.clearBranches();

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


	// calc
    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

    auto nRecHits = ECALRecHit_ID->size();
    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;
    for( int it = 0; it < nRecHits; it++ ){

		auto id = (*ECALRecHit_ID)[it];
		auto idinfo = DetIDMap[id];
		if( idinfo.ecal == ECAL::EB ){

			//auto radius = hypo( (*rhPosX)[it], (*rhPosY)[it] );

		}//<<>>if( (*rhSubdet)[it] == 0 )

	}//<<>>for( int it = 0; it < nRecHits; it++ )

	// fill


}//<<>>void KUCMSAodSkimmer::processRechits()

void KUCMSAodSkimmer::processGenParticles(){

	// initilize

	// calc
    if( DEBUG ) std::cout << "Finding genParticles" << std::endl;
    //------------  genparts ------------------------

	auto nGenParts = Gen_pdgId->size();
    for( int it = 0; it < nGenParts; it++ ){

        auto genID = (*Gen_pdgId)[it];

    }//<<>>for( int it = 0; it < nGenParts; it++ )

	// fill


}//<<>>void KUCMSAodSkimmer::processGenParticles()

void KUCMSAodSkimmer::processCalojets(){

	// initilize

	// calc
	int nCaloJets = 0;
	for( int it = 0; it < nCaloJets; it++ ){


	}//<<>>for( int it = 0; it < nCaloJets; it++ )

	//fill


}//<<>>void KUCMSAodSkimmer::processCalojets()

void KUCMSAodSkimmer::processPhotons(){

    //bool verbose = true;
    bool verbose = false;

	// intilize
	selPhotons.clearBranches();

	// calc
    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

    int leadPho(-1);
    int subLeadPho(-1);
    float leadPhoPt(-1);
    float subLeadPhoPt(-1);
    uInt nSelPhotons = 0;
    uInt nPhotons = Photon_excluded->size();	
    if( DEBUG || verbose ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( uInt it = 0; it < nPhotons; it++ ){

		//---------------------------------------------------
        if( DEBUG ) std::cout << " -- pho isEB, has min rhs, not excluded" << std::endl;
		auto isExcluded = (*Photon_excluded)[it];
        //auto isCmb = phoExcluded[it];
		if( isExcluded ) continue; 
        auto isEB = (*Photon_seedIsEB)[it];
		if( not isEB ) continue;
        auto rhids = (*Photon_rhIds)[it];
        uInt nrh = rhids.size();
		if( nrh < 5 ) continue;

		//--------------------------------------------------------------
        if( DEBUG ) std::cout << " -- pho pull info" << std::endl;
        auto genIdx = (*Photon_genIdx)[it];
		auto isOOT = (*Photon_isOot)[it];
		auto time = (*Photon_seedTOFTime)[it];
		auto eta = (*Photon_eta)[it];
        auto phi = (*Photon_phi)[it];
        auto pt = (*Photon_pt)[it];
		auto smaj = (*Photon_smaj)[it];
        auto smin = (*Photon_smin)[it];
        auto r9 = (*Photon_r9)[it];
        auto sieie = (*Photon_sieie)[it];
        auto energy = (*Photon_energy)[it];

        //--------------------------------------------------------------
        if( DEBUG ) std::cout << " -- pho get calclated values" << std::endl;
        //auto isSUSY = (genIdx >= 0)?((*genLLP)[genIdx] < 700 ):false;
        //auto isSigPho = (genIdx >= 0)?((*genLLP)[genIdx] == 1):false;
        auto phoQuality = getPhoQuality(it);
        auto phoClstrR9 = clstrR9( rhids );
        auto phoEigens2D = getRhGrpEigenFromAngles( rhids );
        auto evaluegeo = phoEigens2D[2];
		auto geosmaj = phoEigens2D[3];
        auto geosmin = phoEigens2D[4];

        // pho object selection ------------------------------------------
        if( DEBUG ) std::cout << " -- pho obj selection" << std::endl;
		auto isMinMedQuality = phoQuality > 1;
		auto underMaxSMaj = smaj = 1.3;
        auto underMaxSMin = smin <= 0.4;
		auto overMinR9 = r9 >= 0.9;
		auto underMaxSieie = sieie <= 0.014;
        auto overMinRhCnt = nrh >= 20;

		auto phoSelected = isMinMedQuality;
		if( not phoSelected ) continue;
		nSelPhotons++;

        if( pt > leadPhoPt ){ subLeadPho = leadPho; subLeadPhoPt = leadPhoPt; leadPho = it; leadPhoPt = pt; }
        else if( pt > subLeadPhoPt ){ subLeadPho = it; subLeadPhoPt = pt; }

		// fill ( vectors )
		if( DEBUG ) std::cout << " -- pho fill out branches" << std::endl;
		selPhotons.fillBranch( "selPhoQuality", phoQuality );
        selPhotons.fillBranch( "selPhoTime", time );
		selPhotons.fillBranch( "selPhoGeoEgnVal", evaluegeo );
        selPhotons.fillBranch( "selPhoEta", eta );
        selPhotons.fillBranch( "selPhoPhi", phi );
        selPhotons.fillBranch( "selPhoPt", pt );
        selPhotons.fillBranch( "selPhoSMaj", smaj );
        selPhotons.fillBranch( "selPhoSMin", smin );
        selPhotons.fillBranch( "selPhoClstrRn", phoClstrR9 );
        selPhotons.fillBranch( "selPhoR9", r9 );
        selPhotons.fillBranch( "selPhoSieie", sieie );
        selPhotons.fillBranch( "selPhoNrh", nrh );
        selPhotons.fillBranch( "selPhoEnergy" , energy );
        selPhotons.fillBranch( "selPhoGeoSMaj", geosmaj );
        selPhotons.fillBranch( "selPhoGeoSMin", geosmin );
		if( verbose ) std::cout << " -- selPho Pt: " << pt << " phi: " << phi << " geo: " << evaluegeo << " clrn: " << phoClstrR9;
		if( verbose ) std::cout << " nrh: " << nrh << " quality: " << phoQuality << std::endl;

    }//<<>>for( int it = 0; it < nPhotons; it++ )
    if( DEBUG ) std::cout << " -- pho loop finished" << std::endl;
	// fill ( other )

    if( DEBUG || verbose ) std::cout << " - Selected " << nSelPhotons << " photons" << std::endl;
    geCnts.set( "nSelPhotons", nSelPhotons );	
	selPhotons.fillBranch( "nSelPhotons",  nSelPhotons );
	
	uInt leadPhoIdx = ( nSelPhotons > 0 ) ? leadPho : 9999;
	uInt subLeadPhoIdx = ( nSelPhotons > 1 ) ? subLeadPho : 9999;
    geCnts.set("leadPho",leadPhoIdx);
    geCnts.set("subLeadPho",subLeadPhoIdx);
	selPhotons.fillBranch( "leadSelPho", leadPho );
	selPhotons.fillBranch( "subLeadSelPho", subLeadPho );
	if( verbose ) std::cout << " -- pho lead & sublead selected : " << leadPho << " - " << subLeadPho << std::endl;
	if( verbose ) std::cout << " -- pho lead & sublead idx selected : " << leadPhoIdx << " - " << subLeadPhoIdx << std::endl;	
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

	int nElectrons = 0;
    if( DEBUG ) std::cout << " - Looping over for " << nElectrons << " electrons" << std::endl;
    for( int it = 0; it < nElectrons; it++ ){

    }//<<>>for( int it = 0; it < nElectrons; it++ )

}//<<>>void KUCMSAodSkimmer::processElectrons

void KUCMSAodSkimmer::processJets(){

	// intilize
	selJets.clearBranches();

	// calc
    if( DEBUG ) std::cout << "Finding jets" << std::endl;
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
        if( DEBUG ) std::cout << " - Finding Jet Quality" << std::endl;
		auto quality = getJetQuality(it);
		auto pt = (*Jet_pt)[it];
		auto eta = (*Jet_eta)[it];
        auto phi = (*Jet_phi)[it];
		auto rhids = (*Jet_drRhIds)[it];

		// get cacluated values-------------------------------------------
		if( DEBUG ) std::cout << " - Finding Jet Time." << std::endl;
		auto rhenergies = getRhGrpEnergies( rhids );
		auto rhtimes = getRhGrpTimes( rhids );
		auto timedist = getDistStats( rhtimes, rhenergies );
		auto time = timedist[6];

		if( DEBUG ) std::cout << " - Jet Obj selection." << std::endl;
		// jet object selection ------------------------------------------
		auto overMinPt = pt >= 40; 
		auto underMaxEta = std::abs(eta) <= 3.0;
		auto isMinQuality = quality > 1;
	
		auto jetSelected = underMaxEta && isMinQuality && overMinPt;
		if( not jetSelected ) continue;
		nSelJets++;

        seljetpt.push_back(pt);
        seljeteta.push_back(eta);
        seljetphi.push_back(phi);
        seljetmass.push_back(mass);

		// fill vectors
        selJets.fillBranch( "nJets", nJets);
        selJets.fillBranch( "nSelJets", nSelJets);
		selJets.fillBranch( "selJetQuality", quality );
        selJets.fillBranch( "selJetPt", pt);
        selJets.fillBranch( "selJetMass", mass);
        selJets.fillBranch( "selJetEnergy", energy);
        selJets.fillBranch( "selJetEta", eta);
        selJets.fillBranch( "selJetPhi", phi);
        selJets.fillBranch( "selJetTime", time);

	}//<<>>for( int it = 0; it < nJets; it++ )
    if( DEBUG ) std::cout << " - Finished Jet loop." << std::endl;

	// fill other
	geVects.set( "selJetPt", seljetpt );
    geVects.set( "selJetEta", seljeteta );
    geVects.set( "selJetPhi", seljetphi );
    geVects.set( "selJetMass", seljetmass );
	geCnts.set( "nSelJets", nSelJets );
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
	selRjrVars.clearBranches();

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

	if( DEBUG ) std::cout << " - Loading MET." << std::endl;
	phoRMetCPx =+ leadPhoPt*std::cos(leadPhoPhi); 
    if( nSelPhotons > 1 ) phoRMetCPx =+ subLeadPhoPt*std::cos(subLeadPhoPhi);
    phoRMetCPy =+ leadPhoPt*std::sin(leadPhoPhi); 
    if( nSelPhotons > 1 ) phoRMetCPy =+ subLeadPhoPhi*std::sin(subLeadPhoPhi);
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
        auto sjetMass = selJetMass[it]; //selJets.getFLBranchValue( "selJetMass", it );
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

	auto nSelJets = geCnts("nSelJets"); //selJets.getUIBranchValue("nSelJets");
    auto nSelPhotons = geCnts("nSelPhotons"); //selPhotons.getUIBranchValue("nSelPhotons");
	auto leadSelPho = geCnts("leadPho"); //selPhotons.getUIBranchValue("leadSelPho");
    auto subLeadSelPho = geCnts("subLeadPho"); //selPhotons.getUIBranchValue("subLeadSelPho");
	auto leadPhoPt = ( nSelPhotons > 0 ) ? geVars("leadPhoPt") : 0;
    auto subLeadPhoPt = ( nSelPhotons > 1 ) ? geVars("subLeadPhoPt") : 0;
	if( DEBUG ) std::cout << " - Lead/Sublead Photons: " << leadSelPho << " - " << subLeadSelPho << std::endl;

    auto gt1phos = nSelPhotons >= 1;
    auto gt2jets = nSelJets >= 2;
    auto gt2phos = nSelPhotons >= 2;
	auto leadPhoPt70 = leadPhoPt >= 70;
    auto leadPhoPt20 = leadPhoPt >= 20;
	auto subLeadPhoPt40 = subLeadPhoPt >= 40; 
	
    auto evtSelected = gt2jets && gt1phos && leadPhoPt20;
	//auto evtSelected = leadPhoPt70 && subLeadPhoPt40 && gt2jets && gt2phos;

	if( DEBUG ){ if( evtSelected ) std::cout << " - Event Passed." << std::endl; else std::cout << " - Event Failed." << std::endl;}
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
    auto rhIso = (*Photon_ecalRHSumEtConeDR04)[it] <= ( 0.006*(*Photon_pt)[it] + 4.2 );
    auto hcalTowIso = (*Photon_hcalTowerSumEtBcConeDR04)[it] <= ( 0.0025*(*Photon_pt)[it] + 2.2 );
    //if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
    auto hcTrkIsoL = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 3.5 ); //hallow cone track iso
    auto hcTrkIsoT = (*Photon_trkSumPtSolidConeDR04)[it] <= ( 0.001*(*Photon_pt)[it] + 2 ); //hallow cone track iso
    //if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
    auto hadOverE = (*Photon_hadTowOverEM)[it] <= 0.05;
    auto sieie = (*Photon_sieie)[it];
    auto sigmaIeieEE = sieie <= 0.03; // tight only EE
    auto sigmaIeieEB = sieie <= 0.013; // tight only EB

    //if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
    auto baseCut = rhIso && hcalTowIso && hadOverE;
    auto looseCut = baseCut && hcTrkIsoL;
    auto tightCut = baseCut && hcTrkIsoT;
    auto tightEB = tightCut && sigmaIeieEB;
    auto tightEE = tightCut && sigmaIeieEE;

    auto phoClass = tightCut ? 3 : looseCut ? 2 : baseCut ? 1 : 0;

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
    selEvtVars.attachBranches( fOutTree );

	//selMet.makeBranch( "Met", FLOAT );
	selMet.makeBranch( "met", "selCMet", FLOAT, "Magnitude of event Met corrected for OOT photons" );
    selMet.makeBranch( "metPx", "selCMetPx", FLOAT, "Magnitude of event MetPx corrected for OOT photons" );
    selMet.makeBranch( "metPy", "selCMetPy", FLOAT, "Magnitude of event MetPy corrected for OOT photons" );
	selMet.attachBranches( fOutTree );

    selPhotons.makeBranch( "nPhotons", UINT );
    selPhotons.makeBranch( "nSelPhotons", UINT ); 
    selPhotons.makeBranch( "leadSelPho", UINT );
    selPhotons.makeBranch( "subLeadSelPho", UINT );
    selPhotons.makeBranch( "selPhoQuality", VINT ); 
    selPhotons.makeBranch( "selPhoTime", VFLOAT ); 
    selPhotons.makeBranch( "selPhoGeoEgnVal", VFLOAT ); 
    selPhotons.makeBranch( "selPhoEta", VFLOAT ); 
    selPhotons.makeBranch( "selPhoPhi", VFLOAT );     
	selPhotons.makeBranch( "selPhoPt", VFLOAT ); 
    selPhotons.makeBranch( "selPhoSMaj", VFLOAT ); 
    selPhotons.makeBranch( "selPhoSMin", VFLOAT ); 
    selPhotons.makeBranch( "selPhoGeoSMaj", VFLOAT );
    selPhotons.makeBranch( "selPhoGeoSMin", VFLOAT );
    selPhotons.makeBranch( "selPhoClstrRn", VFLOAT );
    selPhotons.makeBranch( "selPhoR9", VFLOAT );
    selPhotons.makeBranch( "selPhoSieie", VFLOAT );
    selPhotons.makeBranch( "selPhoNrh", VUINT );
    selPhotons.makeBranch( "selPhoEnergy", VFLOAT );
    selPhotons.attachBranches( fOutTree );

    //selJets.makeBranch( "JetHt", &JetHt );
    selJets.makeBranch( "nJets", UINT );
    selJets.makeBranch( "nSelJets", UINT ); 
    selJets.makeBranch( "selJetQuality", VINT ); 
    selJets.makeBranch( "selJetPt", VFLOAT ); 
    selJets.makeBranch( "selJetEnergy", VFLOAT );
    selJets.makeBranch( "selJetEta", VFLOAT ); 
    selJets.makeBranch( "selJetPhi", VFLOAT );
    selJets.makeBranch( "selJetTime", VFLOAT ); 
    selJets.makeBranch( "selJetMass", VFLOAT );
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

	fConfigTree->Fill();

}//<<>>void KUCMSAodSkimmer::fillConfigTree( TTree* fConfigTree, std::string key )

void KUCMSAodSkimmer::initHists(){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

    ////hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

    ////hist2d[1] = new TH2D("jetDrMuTime_pt", "jetDrMuTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    
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
    return sumRhEn > 0 ? leadRhEn/sumRhEn  : 1.2;

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
		auto rhtof = (*ECALRecHit_TOF)[getRhIdx(id)];
		result.push_back(rhtime+rhtof); 
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
    for( uInt it(0); it < nRecHits; it++ ){

        const auto rhIDX = getRhIdx(rechitids[it]);
        auto idinfo = DetIDMap[rechitids[it]];
        auto isEB = idinfo.ecal == ECAL::EB;
        if( rhIDX == -1 ){ if( verbose ) std::cout << " ---- Bad idx !!!!! -- In getRhGrpEigen ---- " << std::endl; return emptyReturn; }
        if( not isEB ){ if( verbose ) std::cout << " ---- rechit group has EE members " << idinfo.ecal << std::endl; return emptyReturn; }

    	const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
        rhetas.push_back((rhEtaPos>0)?rhEtaPos+84.5:rhEtaPos+85.5);
        const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
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
    auto sxx = var( detas, 0., redlogwtvec );
    auto syy = var( dphis, 0., redlogwtvec, accum(redlogwtvec)/phiCorrFactor );
    auto sxy = cvar( detas, 0., dphis, 0., redlogwtvec, accum(redlogwtvec)/std::sqrt(phiCorrFactor) );
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

