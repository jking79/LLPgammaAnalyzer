//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "KUCMSAodSkimmer_class.hh"

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// KUCMSAodSkimmer class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

//#define DEBUG true
#define DEBUG false

//#define CLSTRMAPS true
#define CLSTRMAPS false

KUCMSAodSkimmer::KUCMSAodSkimmer(){

    disphotreename = "tree/llpgtree";
    //const string KUCMSAodSkimmer::eosdir = "root://cmseos.fnal.gov//store/user/jaking/";
    eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/";
    listdir = "../llpgana_list_files/";

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

void KUCMSAodSkimmer::kucmsAodSkimmer( std::string indir, std::string infilelist, std::string outfilename, int pct ){

	std::cout << "Processing Events for : " << outfilename << std::endl;
    std::ifstream infile(listdir+infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
		//std::cout << "--  for Fine #" << cnt << " moduls " << cnt%pct << " ";
		if( cnt%pct == 0 ){ 
        	auto tfilename = eosdir + indir + str;
        	std::cout << "--  adding file: " << tfilename << std::endl;
        	fInTree->Add(tfilename.c_str());
		}//<<>>if( cnt%4 == 0 ){
		//else std::cout << " do not add file" << std::endl;
		cnt++;
    }//<<>>while (std::getline(infile,str))

	auto fOutTree = new TTree("kuSkimTree","output root file for kUCMSSkimmer");
    auto fConfigTree = new TTree("kuSkimConfigTree","config root file for kUCMSSkimmer");

	KUCMSAodSkimmer::Init(fInTree);
	initHists();
    setOutputBranches(fOutTree);	

	SetupDetIDsEB(DetIDMap);
	SetupDetIDsEE(DetIDMap);

	startJobs();

    std::cout << "Setting up For Main Loop." << std::endl;
	int loopCounter(10000);
    auto nEntries = fInTree->GetEntries();
    if(DEBUG){ nEntries = 100000; loopCounter = 1; }
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
        auto entry = fInTree->LoadTree(centry);
		getBranches(entry);
		auto saveToTree = eventLoop(entry);
		if( saveToTree ) fOutTree->Fill();
    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

	endJobs(fConfigTree);
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

    std::cout << "KUCMSAodSkimmer : Thats all Folks!!" << std::endl;

    fOutFile->Close();

	delete fInTree;
	delete fOutTree;
	delete fConfigTree;
	delete fOutFile;

}//<<>>void kucmsSkimmer

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------
// event loop and startup jobs 
//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------


void KUCMSAodSkimmer::startJobs(){

	nEvents = 0;
	nSelectedEvents = 0;

};

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

	// initilize
	RunNumber = 0;
	
	// calc
    nEvents++;
    
	//fill
	RunNumber = run;

}//<<>>void KUCMSAodSkimmer::processEvntVars()

void KUCMSAodSkimmer::processMet(){

	//intilize
	Met =0;

	//calc
	auto met = std::sqrt(sq2(metCPx)+sq2(metCPy));

	//fill
	Met = met;

}//<<>>void KUCMSAodSkimmer::processMet()

void KUCMSAodSkimmer::processRechits(){

	// initilize


	// calc
    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

    auto nRecHits = rhID->size();
    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;
    for( int it = 0; it < nRecHits; it++ ){

		auto id = (*rhID)[it];
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

	auto nGenParts = genPdgId->size();
    for( int it = 0; it < nGenParts; it++ ){

        auto genID = (*genPdgId)[it];

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

	// intilize
	nSelPhotons = 0;
	leadSelPho = 0;
	subLeadSelPho = 0;
	selPhoQuality.clear();
	selPhoTime.clear();
	selPhoGEgnVal.clear();
	selPhoEta.clear();
    selPhoPhi.clear();
	selPhoPt.clear();
	selPhoSMaj.clear();
	selPhoSMin.clear();
	selPhoClstrRn.clear();

	// calc
    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

    auto selPhotons = 0;
    auto nPhotons = phoExcluded->size();	
    if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( int it = 0; it < nPhotons; it++ ){


		// detrimine photon classification ( sig/susy/ect... )
		//---------------------------------------------------
		auto isExcluded = (*phoExcluded)[it];
        //auto isCmb = phoExcluded[it];
        auto genIdx = (*phoGenIdx)[it];
        //auto isSUSY = (genIdx >= 0)?((*genLLP)[genIdx] < 700 ):false;
		auto isOOT = (*phoIsOotPho)[it];
		//auto isSigPho = (genIdx >= 0)?((*genLLP)[genIdx] == 1):false;
		auto phoQuality = getPhoQuality(it);

		//----------------------------------------------------------------
		if( DEBUG ) std::cout << " -- pho rechits" << std::endl;
        auto rhids = (*phoRhIds)[it];
		auto nrh = rhids.size();
		auto time = (*phoSeedTOFTime)[it];
		auto eta = (*phoEta)[it];
        auto phi = (*phoPhi)[it];
        auto isEB = (*phoIsEB)[it];
        auto pt = (*phoPt)[it];
		auto smaj = (*phoSMaj)[it];
        auto smin = (*phoSMin)[it];

		// pho object selection ------------------------------------------
		if( DEBUG ) std::cout << " -- pho obj selection" << std::endl;
		auto isEE = not isEB;
		auto isNotMinLoose = phoQuality < 2;
		auto underMinRecHits  = nrh < 15;
		auto overSMajMin = smaj > 1.3 || smin > 0.4;

		auto failSelection = isEE || isNotMinLoose || underMinRecHits || overSMajMin || isExcluded;
		if( failSelection ) continue;
		selPhotons++;

		if( DEBUG ) std::cout << " -- pho get rhgroup calcs clr9" << std::endl;
        auto phoClstrR9 = clstrR9( rhids );
        if( DEBUG ) std::cout << " -- pho get rhgroup calcs tof times" << std::endl;
		//auto tofTimes = getLeadTofRhTime( rhids, vtxX, vtxY, vtxZ );
        if( DEBUG ) std::cout << " -- pho get rhgroup calcs eigans" << std::endl;
		//auto phoEigens2D = getRhGrpEigen_sph( tofTimes, rhids );
		auto evaluegeo = 0;//phoEigens2D[18];
        //auto isEiganGood = phoEigens2D[0] != -9;

		// fill ( vectors )
		if( DEBUG ) std::cout << " -- pho fill out branches" << std::endl;
    	selPhoQuality.push_back(phoQuality);
    	selPhoTime.push_back(time);
		selPhoGEgnVal.push_back(evaluegeo);
		selPhoEta.push_back(eta);
        selPhoPhi.push_back(phi);
		selPhoPt.push_back(pt);
		selPhoSMaj.push_back(smaj);
		selPhoSMin.push_back(smin);
		selPhoClstrRn.push_back(phoClstrR9);

    }//<<>>for( int it = 0; it < nPhotons; it++ )

	// fill ( other )
	
	nSelPhotons = selPhotons;
	leadSelPho = leadIdx( selPhoPt );
	subLeadSelPho = subldIdx( selPhoPt, leadSelPho );

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
	nSelJets = 0;
    JetHt = 0;
    selJetQuality.clear();
    selJetPt.clear();
    selJetMass.clear();
    selJetEnergy.clear();
	selJetEta.clear();
    selJetPhi.clear();
	selJetTime.clear();

	// calc
    if( DEBUG ) std::cout << "Finding jets" << std::endl;
	//--------- jets --------------------------

	int selJets = 0;
	int nJets = jetE->size(); // temp till new ntupes
    if( DEBUG ) std::cout << " - Looping over for " << nJets << " jets" << std::endl;
    for( int it = 0; it < nJets; it++ ){

		auto energy = (*jetE)[it];
		auto mass = (*jetM)[it];	
		auto quality = getJetQuality(it);
		auto pt = (*jetPt)[it];
		auto eta = (*jetEta)[it];
        auto phi = (*jetPhi)[it];
		auto time = -25.0; //(*jetCMeanTime)[it];

		// jet object selection ------------------------------------------
		auto underMinEnergy = energy < 30; // should be pt
		auto overMaxEta = eta > 3.0;
		
		auto failJetSelection = underMinEnergy || overMaxEta;
		if( failJetSelection ) continue;
		selJets++;

		// fill vectors
		selJetQuality.push_back(quality);
		selJetPt.push_back(pt);
        selJetMass.push_back(mass);
		selJetEnergy.push_back(energy);
		selJetEta.push_back(eta);
        selJetPhi.push_back(phi);
		selJetTime.push_back(time);

	}//<<>>for( int it = 0; it < nJets; it++ )

	// fill other
	JetHt = 0; //jetHt;
	nSelJets = selJets;

}//<<>>void KUCMSAodSkimmer::processJets()

//------------------------------------------------------------------------------------------------------------
// process RJR for event
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::processRJR(){

	LAB->ClearEvent();

	auto phoRMetCPx = metCPx;
	phoRMetCPx =+ selPhoPt[leadSelPho]*std::cos(selPhoPhi[leadSelPho]); 
    phoRMetCPx =+ selPhoPt[subLeadSelPho]*std::cos(selPhoPhi[subLeadSelPho]);
    auto phoRMetCPy = metCPy;
    phoRMetCPy =+ selPhoPt[leadSelPho]*std::sin(selPhoPhi[leadSelPho]); 
    phoRMetCPy =+ selPhoPt[subLeadSelPho]*std::sin(selPhoPhi[subLeadSelPho]);
	TVector3 ETMiss(phoRMetCPx,phoRMetCPy,0);
	INV->SetLabFrameThreeVector(ETMiss);

	std::vector<RFKey> jetID;
  	for( int it = 0; it < nSelJets; it++ ){ 
		TLorentzVector jet;
		jet.SetPtEtaPhiM(selJetPt[it],selJetEta[it],selJetPhi[it],selJetMass[it]);
		jetID.push_back(COMB_J->AddLabFrameFourVector(jet)); 
	}//<<>>for( int i = 0; i < nSelJets; i++ )

  	if( !LAB->AnalyzeEvent() ) std::cout << "Something went wrong with tree event analysis" << std::endl;
	
  	auto m_MS = S->GetMass();
  	//auto m_PS = S->GetMomentum(*CM);
  	auto m_cosS  = S->GetCosDecayAngle();
  	auto m_dphiS = S->GetDeltaPhiDecayAngle();
  	auto m_dphiSI  = S->GetDeltaPhiBoostVisible();
  	auto m_PTS = S->GetFourVector().Pt();
  	auto m_PzS = S->GetFourVector().Pz();
  
  	auto m_MX2a = X2a->GetMass();
  	//auto m_cosX3a = X3a->GetCosDecayAngle();
  	auto m_MX2b = X2b->GetMass();
  	//auto m_cosX3b = X3b->GetCosDecayAngle();
  	auto m_EVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).E();
  	auto m_EVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).E();
  	auto m_PVa = X2a->GetListVisibleFrames().GetFourVector(*X2a).P();
  	auto m_PVb = X2b->GetListVisibleFrames().GetFourVector(*X2b).P();

  	auto m_MX1a = X1a->GetMass();
  	//auto m_cosX1a = X1a->GetCosDecayAngle();
  	auto m_MX1b = X1b->GetMass();
  	//auto m_cosX1b = X1b->GetCosDecayAngle();  

  	auto m_MV = S->GetListVisibleFrames().GetMass();
  	auto m_PV = S->GetListVisibleFrames().GetFourVector(*S).P();
  	auto m_MVa = X2a->GetListVisibleFrames().GetMass();
  	auto m_MVb = X2b->GetListVisibleFrames().GetMass();

  	auto m_PV_lab    = S->GetListVisibleFrames().GetFourVector().P();
  	auto m_dphiMET_V = S->GetListVisibleFrames().GetFourVector().Vect().DeltaPhi(ETMiss);

}

//------------------------------------------------------------------------------------------------------------
// decide which events to save to tree
//------------------------------------------------------------------------------------------------------------

bool KUCMSAodSkimmer::eventSelection(){
// select which events to save and fill output branches

	// determine if we want to save event
	
	auto gt2jets = nSelJets > 2;
	auto leadPhoPt70 = ( nSelPhotons > 0 ) ? selPhoPt[leadSelPho] > 70 : false;
	auto subLeadPhoPt40 = ( nSelPhotons > 1 ) ? selPhoPt[subLeadSelPho] > 40 : false; 
	
	auto evtSelected = leadPhoPt70 && subLeadPhoPt40 && gt2jets;
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
    auto rhIso = (*phoEcalRHSumEtConeDR04)[it] < ( 0.006*(*phoPt)[it] + 4.2 );
    auto hcalTowIso = (*phoHcalTowerSumEtBcConeDR04)[it] < ( 0.0025*(*phoPt)[it] + 2.2 );
    if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
    auto hcTrkIsoL = (*phoTrkSumPtSolidConeDR04)[it] < ( 0.001*(*phoPt)[it] + 3.5 ); //hallow cone track iso
    auto hcTrkIsoT = (*phoTrkSumPtSolidConeDR04)[it] < ( 0.001*(*phoPt)[it] + 2 ); //hallow cone track iso
    if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
    auto hadOverE = (*phohadTowOverEM)[it] < 0.05;
    auto sieie = (*phoSigmaIEtaIEta)[it];
    auto sigmaIeieEE = sieie < 0.03; // tight only EE
    auto sigmaIeieEB = sieie < 0.013; // tight only EB

    if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
    auto baseCut = rhIso && hcalTowIso && hadOverE;
    auto looseCut = baseCut && hcTrkIsoL;
    auto tightCut = baseCut && hcTrkIsoT;
    auto tightEB = tightCut && sigmaIeieEB;
    auto tightEE = tightCut && sigmaIeieEE;

    auto phoClass = tightCut?3:looseCut?2:1;

		return phoClass;

}//<<>>int KUCMSAodSkimmer::getPhoQuality( int iter )

int KUCMSAodSkimmer::getJetQuality( int it ){

    const auto eta  = std::abs((*jetEta)[it]);	     
    const auto NHF  = (*jetNHF)[it];
    const auto NEMF = (*jetNEMF)[it];
    const auto CHF  = (*jetCHF)[it];
    const auto CEMF = (*jetCEMF)[it];
    const auto NHM  = (*jetNHM)[it];
    const auto CHM  = (*jetCHM)[it];
    const auto SHM  = NHM + CHM;
    const auto MUF  = (*jetMUF)[it];

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

}//<<>>int KUCMSAodSkimmer::getJetQuality( int iter )

//------------------------------------------------------------------------------------------------------------
// get branches, set output branches, initialize histograms, and endjobs
//------------------------------------------------------------------------------------------------------------

void KUCMSAodSkimmer::getBranches( Long64_t entry ){

	//fChain->GetEntry(entry);

   if( DEBUG ) std::cout << " -- get event" << std::endl;
   b_run->GetEntry(entry);   //!
   b_lumi->GetEntry(entry);   //!
   b_event->GetEntry(entry);   //!

   if( DEBUG ) std::cout << " -- get vtx" << std::endl;
   b_vtxX->GetEntry(entry);   //!
   b_vtxY->GetEntry(entry);   //!
   b_vtxZ->GetEntry(entry);   //!

   if( DEBUG ) std::cout << " -- get met" << std::endl;
   //b_metSumEt->GetEntry(entry);   //!
   //b_metPt->GetEntry(entry);   //!
   //b_metPx->GetEntry(entry);   //!
   //b_metPy->GetEntry(entry);   //!
   //b_metPhi->GetEntry(entry);   //!
   b_metCSumEt->GetEntry(entry);   //!
   b_metCPx->GetEntry(entry);   //!
   b_metCPy->GetEntry(entry);   //!

   b_jetHt->GetEntry(entry);   //!
   b_jetE->GetEntry(entry);   //!
   b_jetM->GetEntry(entry);   //!
   b_jetPt->GetEntry(entry);   //!
   b_jetEta->GetEntry(entry);   //!
   b_jetPhi->GetEntry(entry);   //!
   b_jetNHF->GetEntry(entry);   //!
   b_jetNEMF->GetEntry(entry);   //!
   b_jetCHF->GetEntry(entry);   //!
   b_jetCEMF->GetEntry(entry);   //!
   b_jetMUF->GetEntry(entry);   //!
   b_jetNHM->GetEntry(entry);   //!
   b_jetCHM->GetEntry(entry);   //!
   b_jetPHM->GetEntry(entry);   //!

   if( DEBUG ) std::cout << " -- get pho" << std::endl;
   b_phoIsOotPho->GetEntry(entry);   //!
   b_phoExcluded->GetEntry(entry);   //!
   b_phoSeedTOFTime->GetEntry(entry);   //!
   //b_phoSc2dEv->GetEntry(entry);   //!
   b_phoPt->GetEntry(entry);   //!
   b_phoEnergy->GetEntry(entry);   //!
   b_phoPhi->GetEntry(entry);   //!
   b_phoEta->GetEntry(entry);   //!
   b_phoPx->GetEntry(entry);   //!
   b_phoPy->GetEntry(entry);   //!
   b_phoPz->GetEntry(entry);   //!
   b_phoRhIds->GetEntry(entry);   //!
   //b_phoIsPFPhoton->GetEntry(entry);   //!
   //b_phoIsStdPhoton->GetEntry(entry);   //!
   //b_phoHasConTracks->GetEntry(entry);   //!
   //b_phoIsPixelSeed->GetEntry(entry);   //!
   //b_phoIsPhoton->GetEntry(entry);   //!
   b_phoIsEB->GetEntry(entry);   //!
   //b_phoIsEE->GetEntry(entry);   //!
   //b_phoHadOverEM->GetEntry(entry);   //!
   //b_phoHadD1OverEM->GetEntry(entry);   //!
   //b_phoHadD2OverEM->GetEntry(entry);   //!
   //b_phoHadOverEMVaid->GetEntry(entry);   //!
   b_phohadTowOverEM->GetEntry(entry);   //!
   //b_phohadTowD10OverEM->GetEntry(entry);   //!
   //b_phohadTowD20OverEM->GetEntry(entry);   //!
   //b_phohadTowOverEMValid->GetEntry(entry);   //!
   //b_phoE1x5->GetEntry(entry);   //!
   //b_phoE2x5->GetEntry(entry);   //!
   //b_phoE3x3->GetEntry(entry);   //!
   //b_phoE5x5->GetEntry(entry);   //!
   //b_phoMaxEnergyXtal->GetEntry(entry);   //!
   //b_phoSigmaEtaEta->GetEntry(entry);   //!
   b_phoSigmaIEtaIEta->GetEntry(entry);   //!
   //b_phoR1x5->GetEntry(entry);   //!
   //b_phoR2x5->GetEntry(entry);   //!
   b_phoR9->GetEntry(entry);   //!
   //b_phoFull5x5_e1x5->GetEntry(entry);   //!
   //b_phoFull5x5_e2x5->GetEntry(entry);   //!
   //b_phoFull5x5_e3x3->GetEntry(entry);   //!
   //b_phoFull5x5_e5x5->GetEntry(entry);   //!
   //b_phoFull5x5_maxEnergyXtal->GetEntry(entry);   //!
   //b_phoFull5x5_sigmaEtaEta->GetEntry(entry);   //!
   //b_phoFull5x5_sigmaIEtaIEta->GetEntry(entry);   //!
   //b_phoFull5x5_r9->GetEntry(entry);   //!
   b_phoEcalRHSumEtConeDR04->GetEntry(entry);   //!
   //b_phoHcalTwrSumEtConeDR04->GetEntry(entry);   //!
   //b_phoHcalDepth1TowerSumEtConeDR04->GetEntry(entry);   //!
   //b_phoCalDepth2TowerSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalTowerSumEtBcConeDR04->GetEntry(entry);   //!
   //b_phoHcalDepth1TowerSumEtBcConeDR04->GetEntry(entry);   //!
   //b_phoHcalDepth2TowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoTrkSumPtSolidConeDR04->GetEntry(entry);   //!
   //b_phoTrkSumPtHollowConeDR04->GetEntry(entry);   //!
   //b_phoNTrkSolidConeDR04->GetEntry(entry);   //!
   //b_phoNTrkHollowConeDR04->GetEntry(entry);   //!
   b_phoGenIdx->GetEntry(entry);   //!
   b_phoGenDr->GetEntry(entry);   //!

   b_phoSMaj->GetEntry(entry);   //!
   b_phoSMin->GetEntry(entry);   //!
   //b_phoSAlp->GetEntry(entry);   //!
   //b_phoCovEtaEta->GetEntry(entry);   //!
   //b_phoCovEtaPhi->GetEntry(entry);   //!
   //b_phoCovPhiPhi->GetEntry(entry);   //!

   if( DEBUG ) std::cout << " -- eget gen" << std::endl;
   b_genPt->GetEntry(entry);   //!
   //b_genEnergy->GetEntry(entry);   //!
   //b_genPhi->GetEntry(entry);   //!
   //b_genEta->GetEntry(entry);   //!
   //b_genPx->GetEntry(entry);   //!
   //b_genPy->GetEntry(entry);   //!
   //b_genPz->GetEntry(entry);   //!
   b_genPdgId->GetEntry(entry);   //!
   //b_genLLP->GetEntry(entry);   //!

   if( DEBUG ) std::cout << " -- eget rh" << std::endl;
   b_rhEnergy->GetEntry(entry);   //!
   b_rhTime->GetEntry(entry);   //!
   b_rhTOF->GetEntry(entry);   //!
   b_rhID->GetEntry(entry);   //!
   b_rhisOOT->GetEntry(entry);   //!

}//<<>>void KUCMSAodSkimmer::getBranches( Long64_t entry )

void KUCMSAodSkimmer::setOutputBranches( TTree* fOutTree ){


	//fOutTree->Branch( "RunNumber", &RunNumber );

    fOutTree->Branch( "Met", &Met ); 

    fOutTree->Branch( "nSelPhotons", &nSelPhotons ); 
    fOutTree->Branch( "selPhoQualtiy", &selPhoQuality ); 
    fOutTree->Branch( "selPhoTime", &selPhoTime ); 
    fOutTree->Branch( "selPhoGEgnVal", &selPhoGEgnVal ); 
    fOutTree->Branch( "selPhoEta", &selPhoEta ); 
    fOutTree->Branch( "selPhoPt", &selPhoPt ); 
    fOutTree->Branch( "selPhoSMaj", &selPhoSMaj ); 
    fOutTree->Branch( "selPhoSMin", &selPhoSMin ); 
    fOutTree->Branch( "selPhoClstrRn", &selPhoClstrRn );

    fOutTree->Branch( "JetHt", &JetHt );
    fOutTree->Branch( "nSelJets", &nSelJets ); 
    fOutTree->Branch( "selJetQuality", &selJetQuality ); 
    fOutTree->Branch( "selJetPt", &selJetPt ); 
    fOutTree->Branch( "selJetEnergy", &selJetEnergy );
    fOutTree->Branch( "selJetEta", &selJetEta ); 
    fOutTree->Branch( "selJetTime", &selJetTime ); 


}//<<>>void KUCMSAodSkimmer::setBranches( TTree& fOutTree )

void KUCMSAodSkimmer::endJobs( TTree* fConfigTree ){

	TBranch *nEventBranch = fConfigTree->Branch( "nEvents", &nEvents );
	nEventBranch->Fill();
	TBranch *nSelectedEventsBranch = fConfigTree->Branch( "nSelectedEvents", &nSelectedEvents );
	nSelectedEventsBranch->Fill();
	fConfigTree->Fill();

}//<<>>void KUCMSAodSkimmer::endJobs()

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

    for( int idx = 0; idx < rhID->size(); idx++ ){ if( rhDetID == (*rhID)[idx] ) return idx; }
    //std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
    return -1;

}//<<>>int KUCMSAodSkimmer::getRhIdx( int rhDetID )

uInt KUCMSAodSkimmer::getLeadRhID( std::vector<uInt> recHitIds ){

    uInt result;
    float enr(0.0);
    for( auto id : recHitIds ){
        auto rhenr = (*rhEnergy)[getRhIdx(id)];
        if( rhenr > enr ){ enr = rhenr; result = id; }
    }//<<>>for (const auto recHit : recHits )

    return result;

}//>>>>EcalRecHit KUCMSAodSkimmer::getLeadRh( rhGroup recHitsi

float KUCMSAodSkimmer::clstrR9( std::vector<uInt> recHitIds ){

    auto leadRhID = getLeadRhID( recHitIds );
    auto leadRhEn = (*rhEnergy)[getRhIdx(leadRhID)];
    float sumRhEn(0);
    for ( auto id : recHitIds ){ sumRhEn +=  (*rhEnergy)[getRhIdx(id)]; }
    return sumRhEn > 0 ? leadRhEn/sumRhEn  : 1.2;

}//<<>>float KUCMSAodSkimmer::clstrR9( vector<uInt> recHitIds )

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!11  need rhPosX, rhPosY, & rhPosZ to be saved for following functions   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!

std::vector<float> KUCMSAodSkimmer::getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ ){

    std::vector<float> result;
    if( recHitIds.size() < 1 ){ result.push_back(-99); return result; }
    auto lrhid = getLeadRhID(recHitIds);
    auto lrhidx = getRhIdx(lrhid);
    auto X = 0;//(*rhPosX)[lrhidx];
    auto Y = 0;//(*rhPosY)[lrhidx];
    auto Z = 0;//(*rhPosZ)[lrhidx];
    const auto d_rh = hypo( X, Y, Z);
    const auto d_pv = hypo( X-vtxX, Y-vtxY, Z-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    for( int idx = 0; idx < rhTime->size(); idx++ ){result.push_back((*rhTime)[idx]-tof);}
    return result;

}//>>>>vector<float> KUCMSAodSkimmer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

std::vector<float> KUCMSAodSkimmer::getRhGrpEigen_sph( std::vector<float> times, std::vector<uInt> rechitids ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2
	if( DEBUG ) std::cout << " Starting getRhGrpEigen_sph  " << std::endl;;

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
    for ( auto id : rechitids ){ sumRhEn +=  (*rhEnergy)[getRhIdx(id)]; }
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
            const auto rhXPos = 0;//(*rhPosX)[rhIDX];
            rhxs.push_back(rhXPos);
            const auto rhYPos = 0;//(*rhPosY)[rhIDX];
            rhys.push_back(rhYPos);
            const auto rhZPos = 0;//(*rhPosZ)[rhIDX];
            rhzs.push_back(rhZPos);
            rhtimes.push_back(times[it]);
            auto rhenergy = (*rhEnergy)[rhIDX];
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
        float lphi = dltIPhi(rhphis[it],mphi);
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
    auto taflip = ( ((std::abs(mz) < std::abs(vtxZ)) && (mz*vtxZ > 0) ) ? -1 : 1 )*ebside;

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

