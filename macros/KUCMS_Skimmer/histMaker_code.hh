//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "histMaker_class.hh"

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// HistMaker class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

//HistMaker::HistMaker(){}

//HistMaker::~HistMaker(){}

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, int pct ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("tree/llpgtree");
    //const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
    const std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");
    const std::string listdir("llpgana_list_files/");

	std::cout << "Producing Histograms for : " << outfilename << std::endl;
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

	Init(fInTree);
	initHists();

    std::cout << "Setting up For Main Loop." << std::endl;
	int loopCounter(10000);
    auto nEntries = fInTree->GetEntries();
    if(debug){ nEntries = 10; loopCounter = 1; }
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
        if(debug) std::cout << "*****************************************************************************" << std::endl;
        auto entry = fInTree->LoadTree(centry);
		if(debug) std::cout << " - getBranches " << std::endl;
		getBranches(entry);
		if(debug) std::cout << " - eventLoop " << std::endl;
		eventLoop(entry);
    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop

    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

	endJobs();
	for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }

	nMaps = 0;
	for( int it = 0; it < nEBEEMaps; it++ ){ 
		ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
		ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
		ebeeMapR[it]->Write(); delete ebeeMapR[it];
	}//<<>>for( int it = 0; it < nEBEEMaps; it++ )

    fOutFile->Close();
    std::cout << "histMaker : Thats all Folks!!" << std::endl;
}//<<>>void kucmsSkimmer
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------- event loop -------------------------------------------------------------------------
void HistMaker::eventLoop( Long64_t entry ){

    float   jetPTmin        = 100.0;// for energy/time comp
    int     jetIDmin        = 2; //3;
    float   jetETAmax       = 1.5; //1.5;
    int     minRHcnt        = 15; //32;
    float   minRHenr        = 2.0;
    float   bcMinEnergy     = 0.667;
    int     bcMinRHGrpSize  = 3;
    float   minEmf          = 0.0;//0.2

    //if( DEBUG ) std::cout << "Finding Events" << std::endl;
	//------------ event varibles ------------------

    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;

    for( int it = 0; it < nRecHits; it++ ){

		auto id = (*rhID)[it];
		auto idinfo = DetIDMap[id];
		if( idinfo.ecal == ECAL::EB ){
			auto radius = hypo( (*rhPosX)[it], (*rhPosY)[it] );
			fillTH1( radius, hist1d[300]);
			fillTH1( (*rhTime)[it], hist1d[301]); 
			fillTH1( (*rhEnergy)[it], hist1d[302]); 
			hist2d[350]->Fill( (*rhTime)[it], (*rhEnergy)[it]);
			fillTH1( (*rhisOOT)[it], hist1d[305]);
		} else {
			fillTH1( (*rhTime)[it], hist1d[303]); 
			fillTH1( (*rhEnergy)[it], hist1d[304]); 
			hist2d[351]->Fill( (*rhTime)[it], (*rhEnergy)[it]);
			fillTH1( (*rhisOOT)[it], hist1d[306]);
		}//<<>>if( (*rhSubdet)[it] == 0 )

	}//<<>>for( int it = 0; it < nRecHits; it++ )

    if( true ) { // genparts lock
    if( DEBUG ) std::cout << "Finding genParticles" << std::endl;
    //------------  genparts ------------------------

	int nLlpPho(0);
    for( int it = 0; it < nGenParts; it++ ){

		if( (*genPdgId)[it] == 22 && (*genLLP)[it] == 1 ) nLlpPho++;
        auto genID = (*genPdgId)[it];
        //if( abs((*genPdgId)[it]) > 2000000 ) genID = abs((*genPdgId)[it]) - 2000000 + 200;
        //else if( abs((*genPdgId)[it]) > 1000000 ) genID = abs((*genPdgId)[it]) - 1000000 + 100;
        hist1d[200]->Fill( genID );	

    }//<<>>for( int it = 0; it < nGenParts; it++ )
	hist1d[201]->Fill( nLlpPho );

    }//<<>>if( false ) { // genparts lock


    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------
	
    if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( int it = 0; it < nPhotons; it++ ){


		// detrimine photon classification ( sig/susy/ect... )
		//---------------------------------------------------
		auto isCmb = not (*phoExcluded)[it];
        //auto isCmb = phoExcluded[it];
        auto phoGenIdx = (*genPhoIdx)[it];
        auto isSUSY = (phoGenIdx >= 0)?((*genLLP)[phoGenIdx] < 700 ):false;
		auto isOOT = (*phoIsOotPho)[it];
		auto isSigPho = (phoGenIdx >= 0)?((*genLLP)[phoGenIdx] == 1):false;

        if( DEBUG ) std::cout << " -- looping photons : filling gen" << std::endl;
    
        if(phoGenIdx >= 0){
    
            hist1d[100]->Fill( (*genPt)[phoGenIdx] );
            fillTH1( (*genPhoDr)[it], hist1d[130]); //->Fill( (*genPhoDr)[it] );
            fillTH1( abs((*genPdgId)[phoGenIdx]), hist1d[131]);
            //auto genParentExoID = 98;
            //if( abs((*genLLP)[phoGenIdx]) > 2000000 ) genParentExoID = abs((*genLLP)[phoGenIdx]) - 2000000 + 200;
            //else if( abs((*genLLP)[phoGenIdx]) > 1000000 ) genParentExoID = abs((*genLLP)[phoGenIdx]) - 1000000 + 100;
            //if( abs((*genPdgId)[phoGenIdx]) == 0 ) genParentExoID = 95;
            hist1d[101]->Fill( (*genLLP)[phoGenIdx] );
    
        }//<<>>if(phoGenIdx >= 0)


		// determine pog id class
		// -----------------------------------------------------
		if( DEBUG ) std::cout << " -- pho id" << std::endl;
		auto rhIso = (*phoEcalRHSumEtConeDR04)[it] < ( 0.006*(*phoPt)[it] + 4.2 );
		auto hcalTowIso = (*phoHcalTowerSumEtBcConeDR04)[it] < ( 0.0025*(*phoPt)[it] + 2.2 );
        if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
		auto hcTrkIsoL = (*phoTrkSumPtSolidConeDR04)[it] < ( 0.001*(*phoPt)[it] + 3.5 ); //hallow cone track iso
        auto hcTrkIsoT = (*phoTrkSumPtSolidConeDR04)[it] < ( 0.001*(*phoPt)[it] + 2 ); //hallow cone track iso
        if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
        auto hadOverE = (*phohadTowOverEM)[it] < 0.05;
        auto sigmaIeieEE = (*phoSigmaIEtaIEta)[it] < 0.03; // tight only EE
        auto sigmaIeieEB = (*phoSigmaIEtaIEta)[it] < 0.013; // tight only EB

        if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
		auto baseCut = rhIso && hcalTowIso && hadOverE;
		auto looseCut = baseCut && hcTrkIsoL;
		auto tightCut = baseCut && hcTrkIsoT;
		auto tightEB = tightCut && sigmaIeieEB;
		auto tightEE = tightCut && sigmaIeieEE;

		auto phoClass = tightCut?3:looseCut?2:1;
		hist1d[122]->Fill( phoClass );//122

        if( DEBUG ) std::cout << " -- looping photons : filling ids" << std::endl;
        hist1d[107]->Fill( (*phoIsPixelSeed)[it] );
        hist1d[108]->Fill( (*phoHadOverEM)[it] );
        hist1d[109]->Fill( (*phoR9)[it] );
        hist1d[128]->Fill( (*phoR9)[it] );//
        //hist1d[110]->Fill( (*phoMipTotEnergy)[it] );
        hist1d[111]->Fill( (*phoEcalRHSumEtConeDR04)[it] );
        hist1d[112]->Fill( (*phoHcalTwrSumEtConeDR04)[it] );

		//find intial cut values 
		//----------------------------------------------------------------
		if( DEBUG ) std::cout << " -- pho rechits" << std::endl;
		auto nrh = ((*phoRhIds)[it]).size();
        auto phoClstrR9 = clstrR9( (*phoRhIds)[it] );

        auto isEB = (*phoIsEB)[it];
		auto isTightEB = std::abs((*phoEta)[it]) < 1.45;

		auto isTight = phoClass == 3;
        auto isLoose = phoClass > 1;
        auto isFake = phoClass == 1;

        auto goodRhCnt = nrh > 14;
        auto minPhoPt = (*phoPt)[it] > 30;

        //auto usePho = true;
		//auto usePho = not isSUSY && isFake;
        //auto usePho = not isSUSY && isLoose;
        //auto usePho = not isSUSY && isLoose && isEB && isTightEB;
        //auto usePho = isSUSY && isFake;
        //auto usePho = isSUSY && isLoose && isEB && isTightEB;
        //auto usePho = not isSUSY && isLoose && phoExcluded[it];
        //auto usePho = isLoose && isEB && isTightEB;
        
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isClR9r26;
        //auto usePho = isLoose && isEB && isTightEB && isClR9r68;

        auto usePho = isLoose && isEB && isTightEB && isCmb;
        //auto usePho = isTight && isEB && isTightEB && isCmb && isSigPho;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isSigPho;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && not isSigPho;


		if( DEBUG ) std::cout << " -- photons # " << it << " isSUSY " << isSUSY << " isOOT " << isOOT << " isCmb " << isCmb << std::endl;
		if( usePho && goodRhCnt && minPhoPt ) { //----------------------------------- 
	
	        hist1d[132]->Fill(1);
	        if( isSUSY ) hist1d[132]->Fill(2);
	        if( isOOT ) hist1d[132]->Fill(3);
	        if( isCmb ) hist1d[132]->Fill(4);
	
			if( DEBUG ) std::cout << " -- pho time" << std::endl;
            hist1d[104]->Fill( (*phoPt)[it] );
            hist1d[105]->Fill( (*phoEnergy)[it] );
			hist2d[205]->Fill( (*phoEnergy)[it], (*phoPt)[it] );

	        hist1d[119]->Fill( (*phoCMeanTime)[it] );//c mean
	        hist1d[120]->Fill( (*phoSeedTOFTime)[it] );//lead time 
	        hist1d[121]->Fill( (*phoCMeanTime)[it] - (*phoSeedTOFTime)[it] );//diff
			for( auto rhid : (*phoRhIds)[it] ){ hist1d[115]->Fill((*rhEnergy)[getRhIdx(rhid)]);}

            auto isPho = (phoGenIdx >= 0)?(((*genPdgId)[phoGenIdx] == 22)?1:0):-1;
            hist2d[204]->Fill(isPho,isSUSY);

			if( DEBUG ) std::cout << " Finding Eigans for Photon W/ " << nrh << " rechits." << std::endl;

        }//<<>>if( usePho ) { //----------------------------------- 

    }//<<>>for( int it = 0; it < nPhotons; it++ )

	//-------- electrons --------------------------------------

	
    for( int it = 0; it < nElectrons; it++ ){

        hist1d[350]->Fill( (*eleCMeanTime)[it] );//c mean
        hist1d[351]->Fill( (*eleSeedTOFTime)[it] );//lead time 
        hist1d[352]->Fill( (*eleCMeanTime)[it] - (*eleSeedTOFTime)[it] );//diff

    }//<<>>for( int it = 0; it < nElectrons; it++ )


	//--------- jets --------------------------


	//// ---  base info  ----

	hist1d[9]->Fill(nJets);

    for( int it = 0; it < nJets; it++ ){

        const auto jetepafrac   = (*jetPHEF)[it] + (*jetELEF)[it];;
        const auto jetepe       = (*jetPHE)[it] + (*jetELE)[it];
        const auto jeteme       = (*jetCEMF)[it] + (*jetNEMF)[it];
        const auto jetemfrac    = jeteme/(*jetE)[it];
        const auto jetepfrac    = jetepe/(*jetE)[it];
		
        hist2d[61]->Fill(jetepafrac,jetepfrac);
        hist2d[62]->Fill(jetepfrac,jetemfrac);

        fillTH1((*jetPt)[it],hist1d[0]);//hist1d[12]->Fill(jet.pt());
        fillTH1((*jetPhi)[it],hist1d[1]);//hist1d[13]->Fill(jet.phi());
        fillTH1((*jetEta)[it],hist1d[2]);//hist1d[14]->Fill(jet.eta());

	}//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding genjet" << std::endl;
	//// --- genjet info -----------------------------------

    for( int it = 0; it < nJets; it++ ){

        hist2d[161]->Fill((*jetGenTime)[it],(*jetGenTOF)[it]);
        hist1d[50]->Fill((*jetGenTime)[it]);
    	hist1d[51]->Fill((*jetGenTOF)[it]);

    }//<<>>for( int it = 0; it < nJets; it++ )

}//<<>>void HistMaker::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::getBranches( Long64_t entry ){

	//fChain->GetEntry(entry);

   //b_run->GetEntry(entry);   //!
   //b_lumi->GetEntry(entry);   //!
   //b_event->GetEntry(entry);   //!

   b_nVtx->GetEntry(entry);   //!
   b_vtxX->GetEntry(entry);   //!
   b_vtxY->GetEntry(entry);   //!
   b_vtxZ->GetEntry(entry);   //!

   //b_metSumEt->GetEntry(entry);   //!
   //b_metPt->GetEntry(entry);   //!
   //b_metPx->GetEntry(entry);   //!
   //b_metPy->GetEntry(entry);   //!
   //b_metPhi->GetEntry(entry);   //!
   //b_metCSumEt->GetEntry(entry);   //!
   //b_metCPx->GetEntry(entry);   //!
   //b_metCPy->GetEntry(entry);   //!

   //b_jetHt->GetEntry(entry);   //!
   //b_nJets->GetEntry(entry);   //!
   //b_nGoodDrJets->GetEntry(entry);   //!
   //b_nGoodScJets->GetEntry(entry);   //!
   //b_nGoodBcJets->GetEntry(entry);   //!
   //b_nUnJets->GetEntry(entry);   //!
   //b_jetE->GetEntry(entry);   //!

   b_nPhotons->GetEntry(entry);   //!
   b_phoIsOotPho->GetEntry(entry);   //!
   b_phoExcluded->GetEntry(entry);   //!
   b_phoSeedTOFTime->GetEntry(entry);   //!
   b_phoCMeanTime->GetEntry(entry);   //!
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
   b_phoHasConTracks->GetEntry(entry);   //!
   b_phoIsPixelSeed->GetEntry(entry);   //!
   //b_phoIsPhoton->GetEntry(entry);   //!
   b_phoIsEB->GetEntry(entry);   //!
   b_phoIsEE->GetEntry(entry);   //!
   b_phoHadOverEM->GetEntry(entry);   //!
   b_phoHadD1OverEM->GetEntry(entry);   //!
   b_phoHadD2OverEM->GetEntry(entry);   //!
   b_phoHadOverEMVaid->GetEntry(entry);   //!
   b_phohadTowOverEM->GetEntry(entry);   //!
   b_phohadTowD10OverEM->GetEntry(entry);   //!
   b_phohadTowD20OverEM->GetEntry(entry);   //!
   b_phohadTowOverEMValid->GetEntry(entry);   //!
   b_phoE1x5->GetEntry(entry);   //!
   b_phoE2x5->GetEntry(entry);   //!
   b_phoE3x3->GetEntry(entry);   //!
   b_phoE5x5->GetEntry(entry);   //!
   b_phoMaxEnergyXtal->GetEntry(entry);   //!
   b_phoSigmaEtaEta->GetEntry(entry);   //!
   b_phoSigmaIEtaIEta->GetEntry(entry);   //!
   b_phoR1x5->GetEntry(entry);   //!
   b_phoR2x5->GetEntry(entry);   //!
   b_phoR9->GetEntry(entry);   //!
   b_phoFull5x5_e1x5->GetEntry(entry);   //!
   b_phoFull5x5_e2x5->GetEntry(entry);   //!
   b_phoFull5x5_e3x3->GetEntry(entry);   //!
   b_phoFull5x5_e5x5->GetEntry(entry);   //!
   b_phoFull5x5_maxEnergyXtal->GetEntry(entry);   //!
   b_phoFull5x5_sigmaEtaEta->GetEntry(entry);   //!
   b_phoFull5x5_sigmaIEtaIEta->GetEntry(entry);   //!
   b_phoFull5x5_r9->GetEntry(entry);   //!
   b_phoEcalRHSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalTwrSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalDepth1TowerSumEtConeDR04->GetEntry(entry);   //!
   b_phoCalDepth2TowerSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalTowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoHcalDepth1TowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoHcalDepth2TowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoTrkSumPtSolidConeDR04->GetEntry(entry);   //!
   b_phoTrkSumPtHollowConeDR04->GetEntry(entry);   //!
   b_phoNTrkSolidConeDR04->GetEntry(entry);   //!
   b_phoNTrkHollowConeDR04->GetEntry(entry);   //!
   b_genPhoIdx->GetEntry(entry);   //!
   b_genPhoDr->GetEntry(entry);   //!

   b_phoSMaj->GetEntry(entry);   //!
   b_phoSMin->GetEntry(entry);   //!
   b_phoSAlp->GetEntry(entry);   //!
   b_phoCovEtaEta->GetEntry(entry);   //!
   b_phoCovEtaPhi->GetEntry(entry);   //!
   b_phoCovPhiPhi->GetEntry(entry);   //!

   b_nGenParts->GetEntry(entry);   //!
   b_genPt->GetEntry(entry);   //!
   //b_genEnergy->GetEntry(entry);   //!
   //b_genPhi->GetEntry(entry);   //!
   //b_genEta->GetEntry(entry);   //!
   //b_genPx->GetEntry(entry);   //!
   //b_genPy->GetEntry(entry);   //!
   //b_genPz->GetEntry(entry);   //!
   b_genPdgId->GetEntry(entry);   //!
   b_genLLP->GetEntry(entry);   //!

   b_nRecHits->GetEntry(entry);   //!
   b_rhPosX->GetEntry(entry);   //!
   b_rhPosY->GetEntry(entry);   //!
   b_rhPosZ->GetEntry(entry);   //!
   b_rhPosEta->GetEntry(entry);   //!
   b_rhPosPhi->GetEntry(entry);   //!
   b_rhEnergy->GetEntry(entry);   //!
   b_rhTime->GetEntry(entry);   //!
   b_rhTOF->GetEntry(entry);   //!
   b_rhID->GetEntry(entry);   //!
   b_rhisOOT->GetEntry(entry);   //!

}//<<>>void HistMaker::getBranches( Long64_t entry )

void HistMaker::endJobs(){

/* 1D jets 0-99

    normTH1D(hist1d[3]);
    normTH1D(hist1d[4]);
    normTH1D(hist1d[34]);
    normTH1D(hist1d[35]);
    normTH1D(hist1d[38]);
    normTH1D(hist1d[39]);
    normTH1D(hist1d[36]);
    normTH1D(hist1d[37]);
    normTH1D(hist1d[32]);
    normTH1D(hist1d[33]);
*/

/* 1D cljets 250 - 299

    normTH1D(hist1d[256]);
*/

    thresDivTH2D( hist2d[206], hist2d[207], 0 );
    thresDivTH2D( hist2d[208], hist2d[209], 0 );

    profileTH2D( hist2d[215], hist1d[116], hist1d[118] );


/* jets 0 - 199

    normTH2D(hist2d[100]);
    normTH2D(hist2d[121]);
    normTH2D(hist2d[126]);
    normTH2D(hist2d[127]);
*/

    //1normTH2D(hist2d[202]);
    //normTH2D(hist2d[201]);
    //normTH2D(hist2d[200]);
    //normTH2D(hist2d[226]);
    //normTH2D(hist2d[227]);

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists(){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	// jet time
    //int jtdiv(400);
    //float jtran(8);
    int jtdiv(300);
    float jtran(15);
    int jdtdiv(160);
    float jdtran(4);
    int rhcnt(160);

	// jet id stuff
    //auto stddiv = 120;
    //auto stdtran = 3;

	// eigan cluster 2d maps
    //auto cldiv = 1200;
    //auto cltrn = 30;
    auto cldiv = 36;
    auto cltrn = 19.8;

	// eigan cluster 3d maps
    //auto cl3ddiv = 160;
    //auto cl3dtrn = 4;
    //auto cl3ddiv1 = 160;
    //auto cl3dtrn1 = 4;

	// photon time phi - eta map
    auto clsphdiv = 72/2;
    auto clsphtrn = 19.8/2;

	// time profile maps
    auto cwdiv = 100;
    auto cwtrn = 5/2;

	// time cl slope plots
    auto slmax = 2.0;
    auto slmin = -2.0;
    auto sldiv = 200;

    auto sl3max = 2.0;
    auto sl3min = -2.0;
    auto sl3div = 200;

    //auto chimax = 1.01;
    //auto chimin = 0.91;
    //auto chidiv = 2880;


	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

/*
	//------- jets 0 - 99

    hist1d[0] = new TH1D("jetPt", "jetPt", 500, 0, 5000);
    hist1d[1] = new TH1D("jetPhi", "jetPhi", 700, -3.5, 3.5);
    hist1d[2] = new TH1D("jetEta", "jetEta", 700, -3.5, 3.5);
    hist1d[3] = new TH1D("jetdtmu", "jetdtmu", jdtdiv, -1*jdtran, jdtran);
    hist1d[4] = new TH1D("jetdtmed", "jetdtmed", jdtdiv, -1*jdtran, jdtran);
    hist1d[5] = new TH1D("jetHt", "jetHt", 1000, 0, 5000);
    hist1d[6] = new TH1D("diffPt", "diffPt", 1000, 0, 10);
    hist1d[7] = new TH1D("htPct", "htPct", 100, 0, 1);
    hist1d[8] = new TH1D("dPhi", "dPhi", 70, -3.5, 3.5);

    hist1d[9] = new TH1D("nJet", "nJets", 21, -0.5, 20.5);
    hist1d[10] = new TH1D("nGoodDrJets", "nGoodDrJets", 21, -0.5, 20.5);
    hist1d[11] = new TH1D("nGoodScJets", "nGoodScJets", 21, -0.5, 20.5);
    hist1d[12] = new TH1D("nGoodBcJets", "nGoodBcJets", 21, -0.5, 20.5);
    hist1d[13] = new TH1D("nUnJets", "nUnJets", 101, -0.5, 100.5);
    hist1d[14] = new TH1D("pJets", "pJets", 110, 0, 1.1);
    hist1d[15] = new TH1D("pGoodDrJets", "pGoodDrJets", 110, 0, 1.1);
    hist1d[16] = new TH1D("pGoodScJets", "pGoodScJets", 110, 0, 1.1);
    hist1d[17] = new TH1D("pGoodBcJets", "pGoodBcJets", 110, 0, 1.1);
    hist1d[18] = new TH1D("pGoodBcToScJets", "pGoodBcToScJets", 110, 0, 1.1);

    hist1d[19] = new TH1D("jetDRRHMulti", "jetDRRHMulti", rhcnt, 0, rhcnt);
    hist1d[20] = new TH1D("jetDRTimeError", "jetDRTimeError", 300, 0, 3);
    hist1d[21] = new TH1D("jetDRTimeRMS", "jetDRTimeRMS", 200, 0, 20);
    hist1d[22] = new TH1D("jetDrMuTime", "jetDrMuTime", jtdiv, -1*jtran, jtran);
    hist1d[23] = new TH1D("jetDRMedTime", "jetMedDRTime", jtdiv, -1*jtran, jtran);
    hist1d[24] = new TH1D("jetCDRMuTime", "jetCDRMuTime", jtdiv, -1*jtran, jtran);
    hist1d[25] = new TH1D("jetCDRMedTime", "jetCDRMedTime", jtdiv, -1*jtran, jtran);

    hist1d[26] = new TH1D("jetSCmuTime", "jetSCmuTime", jtdiv, -1*jtran, jtran);
    hist1d[27] = new TH1D("jetSCmedTime", "jetSCmedTime", jtdiv, -1*jtran, jtran);
    hist1d[28] = new TH1D("jetCSCMuTime", "jetCSCMuTime", jtdiv, -1*jtran, jtran);
    hist1d[29] = new TH1D("jetCSCMedTime", "jetCSCMedTime", jtdiv, -1*jtran, jtran);

    hist1d[30] = new TH1D("jetCBCMedTime", "jetCBCMedTime", jtdiv, -1*jtran, jtran);
    hist1d[31] = new TH1D("jetCBCMuTime", "jetCBCMuTime", jtdiv, -1*jtran, jtran);


    hist1d[32] = new TH1D("jetcmudtbc", "jetcmudtbc", jdtdiv, -1*jdtran, jdtran);
    hist1d[33] = new TH1D("jetcmeddtbc", "jetcmeddtbc", jdtdiv, -1*jdtran, jdtran);
    hist1d[34] = new TH1D("jetcmudtdr", "jetcmudtdr", jdtdiv, -1*jdtran, jdtran);
    hist1d[35] = new TH1D("jetcmeddtdr", "jetcmeddtdr", jdtdiv, -1*jdtran, jdtran);
    hist1d[36] = new TH1D("jetcmudtsc", "jetcmudtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[37] = new TH1D("jetcmeddtsc", "jetcmeddtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[38] = new TH1D("jetmeddtsc", "jetmeddtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[39] = new TH1D("jetmudtsc", "jetmudtsc", jdtdiv, -1*jdtran, jdtran);

    hist1d[40] = new TH1D("scbcdt", "scbcdt", jdtdiv, -1*jdtran, jdtran);
    hist1d[41] = new TH1D("nBCinJet", "nBCinJet", 11, -0.5, 10.5);

    hist1d[42] = new TH1D("jetmudtgen", "jetmudtgen", jdtdiv, -1*jdtran, jdtran);
    hist1d[43] = new TH1D("genJetImpactAngle", "genJetImpactAngle", 6600, -0.2, 6.4);

    hist1d[44] = new TH1D("genJetDrMatchJet", "genJetDrMatchJet", 100, 0, 0.1);
    hist1d[45] = new TH1D("genJetSCTimeDiff", "genJetSCTimeDiff", 300, 0, 30);
    hist1d[46] = new TH1D("genJetDrTimeDiff", "genJetSCTimeDiff", 300, 0, 30);

    hist1d[47] = new TH1D("jetGenTimeCut", "jetGenTimeCut", jtdiv, -1*jtran, jtran);
    hist1d[48] = new TH1D("jetGenTimeVar", "jetGenTimeVar", 408, -2, 100);
    hist1d[49] = new TH1D("jetGenTimeNextBX", "jetGenTimeNextBX", 3, -1, 2);
    hist1d[50] = new TH1D("jetGenTime", "jetGenTime", jtdiv, -1*jtran, jtran);
    hist1d[51] = new TH1D("jetGenTOF", "jetGenTOF", 300, 0, 30);
    hist1d[52] = new TH1D("jetGenTimeIsLLP", "jetGenTimeIsLLP", 3, -1, 2);
    hist1d[53] = new TH1D("jetGenTimeLLPPurity", "jetGenTimeLLPPurity", 100, 0, 1);
    hist1d[54] = new TH1D("jetGenNKids", "jetGenNKids", 100, 0, 100);

    hist1d[55] = new TH1D("jetClEtaTimeSlopeRangeA", "Jet Cluster Eta Time Slope Eta > 1.0", 480, -240, 240);
    hist1d[56] = new TH1D("jerClEtaTimeSlopeRangeB", "Jet Cluster Eta Time Slope Eta < 0.5", 480, -240, 240);

    hist1d[57] = new TH1D("jetEginValue3D", "jetEginValue3D", 110, 0, 1.1);
    hist1d[58] = new TH1D("jetEtaPhiAngle3D", "jetEtaPhiAngle3D", 660, -0.2, 6.4);
    hist1d[59] = new TH1D("jetEtaTimAngle3D", "jetEtaTimAngle3D", 660, -0.2, 6.4);
    hist1d[60] = new TH1D("jetEtaPhiAngle2D", "jetEtaTimAngle2D", 660, -0.2, 6.4);
    hist1d[61] = new TH1D("jetEginValueSph", "jetEginValueSph", 150, 0.4, 1.1);
*/

	//----- photons 100 - 249

    hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);
    hist1d[101] = new TH1D("genPhoOriginCode", "genPhoOriginCode;Code",1000,0,1000);

    hist1d[102] = new TH1D("phoClstrR9", "phoClstrR9", 100, 0, 1);

    hist1d[103] = new TH1D("phoNumRecHits", "phoNumRecHits", 80, 0, 80);
    hist1d[104] = new TH1D("phoPtPreCut", "phoPt PreCut", 1000, 0, 1000);
    hist1d[105] = new TH1D("phoEnergyPreCut", "phoEnergy PreCut", 2000, 0, 2000);

    //hist1d[106] = new TH1D("pho3DRotAngle", "pho3DRotAngle", 140, -0.5, 6.5);

    hist1d[107] = new TH1D("phoIsPixelSeed", "phoIsPixelSeed", 3, 0, 2);
    hist1d[108] = new TH1D("phoHadOverEM", "phoHadOverEM", 250, 0, 10);
    hist1d[109] = new TH1D("phoR9", "phoR9", 100, 0, 1);
    //hist1d[110] = new TH1D("phoMipTotEnergy", "phoMipTotEnergy", 250, 0, 250 );
    hist1d[111] = new TH1D("phoEcalRHSumEtConeDR04", "phoEcalRHSumEtConeDR04", 750, 0, 750);
    hist1d[112] = new TH1D("phoHcalTwrSumEtConeDR04", "phoHcalTwrSumEtConeDR04", 750, 0, 750);

    //hist1d[113] = new TH1D("phoEginValueSph", "phoEginValueSph", 150, 0.4, 1.1);
    //hist1d[114] = new TH1D("phoEginSlopeAP", "phoEginSlopeSph 1Way", sldiv, slmin, slmax);
    hist1d[115] = new TH1D("clRhEnergy", "clRhEnergy", 3000, 0, 1000);

    hist1d[116] = new TH1D("pho_etprofile", "Photon rEta Time Profile Sph;MajAxis [cm]", clsphdiv, -1*clsphtrn, clsphtrn);

    //hist1d[117] = new TH1D("phoEginSlopeSphTA", "phoEginSlopeSph TA", sldiv, slmin, slmax);
	hist1d[118] = new TH1D("pho_proFitEtavChi", "Profile Fit Eta v Chi2Prob Sph", cwdiv, -1*cwtrn, cwtrn );

    hist1d[119] = new TH1D("phoClTimePreCut", "phoTime PreCut", jtdiv, -1*jtran, jtran);
    hist1d[120] = new TH1D("phoSeedRhTime", "phoLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[121] = new TH1D("phoSeedTimeDiff", "phoLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[122] = new TH1D("phoClass", "phoClass", 4, 0, 4);

    hist1d[123] = new TH1D("phoClRhTime","phoClRhTime",1000, -50, 50);
    hist1d[124] = new TH1D("phoClSMaj","phoClSMaj",50,0,5);
    hist1d[125] = new TH1D("phoClSMin","phoClSMin",50,0,0.5);
    hist1d[126] = new TH1D("phoSigIEtaIEta","phoSigIEtaIEta",300,0,0.03);

    //hist1d[127] = new TH1D("phoRotAngle","phoRotAngle",140,-0.5,6.5);

    hist1d[128] = new TH1D("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
    //hist1d[129] = new TH1D("phoGeoEginValueSph", "phoGeoEginValueSph", 150, 0.4, 1.1);
    hist1d[130] = new TH1D("phoGenDr", "phoGenDr", 800, 0, 2.0 );
    hist1d[131] = new TH1D("phoGenPdgId", "phoGenPdgId", 50, 0, 50 );
    hist1d[132] = new TH1D("phoType", "phoType(1-pho 2-susy 3-oot 4-cmb)", 6,0,5);

    hist1d[133] = new TH1D("phoGeoValue", "phoGeoValue", 150, 0.4, 1.1);
    //hist1d[134] = new TH1D("phoGeoSlope", "phoGeoSlope", sldiv, slmin, slmax);
    //hist1d[135] = new TH1D("phoGeoAngle","phoGeoAngle",70,0,7);

    //hist1d[136] = new TH1D("pho2DEiganV0","pho2DEiganV0",200,-1,1);
    //hist1d[137] = new TH1D("pho2DEiganV1","pho2DEiganV1",200,-1,1);
    //hist1d[138] = new TH1D("phoG2DEiganV0","phoG2DEiganV0",200,-1,1);
    //hist1d[139] = new TH1D("phoG2DEiganV1","phoG2DEiganV1",200,-1,1);

    hist1d[140] = new TH1D("phoMinDr","phoMinDr",500,0,5);

    //hist1d[141] = new TH1D("pho3DEiganV0","pho3DEiganV0",200,-1,1);
    //hist1d[142] = new TH1D("pho3DEiganV1","pho3DEiganV1",200,-1,1);
    //hist1d[143] = new TH1D("pho3DEiganV2","pho3DEiganV2",200,-1,1);
    //hist1d[144] = new TH1D("pho3d4Value","pho3d4Value",100,0,1);
    //hist1d[145] = new TH1D("pho3dTValue","pho3dTValue",400,0,0.4);

    hist1d[146] = new TH1D("phoPtPostCut", "phoPt PostCut", 1000, 0, 1000);
    hist1d[147] = new TH1D("phoEnergyPostCut", "phoEnergy PostCut", 2000, 0, 2000);
    hist1d[148] = new TH1D("phoClTimePostCut", "phoTime PostCut", jtdiv, -1*jtran, jtran);

    //hist1d[149] = new TH1D("pho3dEV1","pho3dEV1",20,0,20);
    //hist1d[150] = new TH1D("pho3dEV2","pho3dEV2",20,0,20);
    //hist1d[151] = new TH1D("pho3dEV3","pho3dEV3",20,0,20);

    hist1d[152] = new TH1D("phoSMaj", "phoSMaj", 300,0,3);
    hist1d[153] = new TH1D("phoSMin", "phoSMin", 60,0,.6);
    hist1d[154] = new TH1D("phoSAlp", "phoSAlp", 350,-1.75,1.75);
    hist1d[155] = new TH1D("phoCovEtaEta", "phoCovEtaEtaSqrt", 100,0,0.001);
    hist1d[156] = new TH1D("phoCovEtaPhi", "phoCovEtaPhiSqrt", 800,-0.0004,0.0004);
    hist1d[157] = new TH1D("phoCovPhiPhi", "phoCovPhiPhiSqrt", 100,0,0.001);

	//------  genparticles 200 - 249

    hist1d[200] = new TH1D("genPdgId", "genPdgId;PdgId",220,0,220);
    hist1d[201] = new TH1D("nGenLLPPho", "nGenLLPPho", 5, 0, 5 );

	//------  cluster jet 250 - 299
/*
    hist1d[250] = new TH1D("cljTime", "cljTime", jtdiv, -1*jtran, jtran);
    hist1d[251] = new TH1D("cljLeadRhTime", "cljLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[252] = new TH1D("cljClstLeadTimeDiff", "cljClstLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[253] = new TH1D("cljDrTime", "cljDrTime", jtdiv, -1*jtran, jtran);
    hist1d[254] = new TH1D("cljDiffPt", "cljDiffPt", 1000, 0, 10);
    hist1d[255] = new TH1D("cljdPhi", "cljdPhi", 70, -3.5, 3.5);
    hist1d[256] = new TH1D("cljdtmu", "cljdtmu", jdtdiv, -1*jdtran, jdtran);
*/
    //------ ecal rechits 300 - 349
    hist1d[300] = new TH1D("rhRadius", "rhRadius", 300, 128, 131);

    hist1d[301] = new TH1D("ebRhTime", "ebRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[302] = new TH1D("ebRhEnergy", "ebRhEnergy", 1000, 0, 1000);
    hist1d[303] = new TH1D("eeRhTime", "eeRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[304] = new TH1D("eeRhEnergy", "eeRhEnergy", 1000, 0, 1000);
    hist1d[305] = new TH1D("ebRhkOOT", "ebRhkOOT", 3, 0, 1);
    hist1d[306] = new TH1D("eeRhkOOT", "eeRhkOOT", 3, 0, 1);

    //------  electrons 350 - 400

    //hist1d[350] = new TH1D("eleClTime", "eleClTime", jtdiv, -1*jtran, jtran);
    //hist1d[351] = new TH1D("eleSeedRhTime", "eleLeadRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[352] = new TH1D("eleClSeedTimeDiff", "eleClLeadTimeDiff", jtdiv, -1*jtran, jtran);

	//-------- event vars 400 - 450
	
	//hist1d[400] = new TH1D("nphoexclusions", " Number Pho (0) Excl (1) OOT (2) Excl (3) ", 4, 0, 4);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	//------- jets ( time ) 0-49 ------------------------------
/*
	//hist2d[0]
    hist2d[1] = new TH2D("jetDrMuTime_pt", "jetDrMuTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[2] = new TH2D("jetDrMuTime_id", "jetDrMuTime_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
    hist2d[3] = new TH2D("jetDrMuTime_nhf", "jetDrMuTime_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[4] = new TH2D("jetDrMuTime_chf", "jetDrMuTime_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[5] = new TH2D("jetDrMuTime_nemf", "jetDrMuTime_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[6] = new TH2D("jetDrMuTime_cemf", "jetDrMuTime_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[7] = new TH2D("jetDrMuTime_muf", "jetDrMuTime_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[8] = new TH2D("jetDrMuTime_nhm", "jetDrMuTime_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
    hist2d[9] = new TH2D("jetDrMuTime_chm", "jetDrMuTime_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

    hist2d[10] = new TH2D("jetDrMuTime_jetDrMedTime", "jetDrMuTime_jetDrMedTime", jtdiv, -1*jtran, jtran, 200, -10, 10);
    hist2d[11] = new TH2D("jetDrMuTime_rms", "jetDrMuTime_rms", jtdiv, -1*jtran, jtran, 200, 0, 20);
    hist2d[12] = new TH2D("jetDrMuTime_err", "jetDrMuTime_err", jtdiv, -1*jtran, jtran, 300, 0, 3);

    hist2d[13] = new TH2D("jetDrMedTime_pt", "jetDrMedTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[14] = new TH2D("jetDrMedTime_id", "jetDrMedTime_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
    hist2d[15] = new TH2D("jetDrMedTime_nhf", "jetDrMedTime_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[16] = new TH2D("jetDrMedTime_chf", "jetDrMedTime_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[17] = new TH2D("jetDrMedTime_nemf", "jetDrMedTime_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[18] = new TH2D("jetDrMedTime_cemf", "jetDrMedTime_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[19] = new TH2D("jetDrMedTime_muf", "jetDrMedTime_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[20] = new TH2D("jetDrMedTime_nhm", "jetDrMedTime_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
    hist2d[21] = new TH2D("jetDrMedTime_chm", "jetDrMedTime_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

    hist2d[22] = new TH2D("jetDrDtMu_nJets", "jetDrDtMu_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);
    hist2d[23] = new TH2D("jetDrDtMed_nJets", "jetDrDtMed_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);

    hist2d[24] = new TH2D("jetDrMu_nrh", "jetDrMu_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);
    hist2d[25] = new TH2D("jetDrMedTime_nrh", "jetDrMedTime_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);

    hist2d[26] = new TH2D("jetDrDtMu_diffPt", "jetDrDtMu_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[27] = new TH2D("jetDrDtMu_htPct", "jetDrDtMu_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[28] = new TH2D("jetDrDtMu_dPhi", "jetDrDtMu_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);
    hist2d[29] = new TH2D("jetDrDtMed_diffPt", "jetDrDtMed_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[30] = new TH2D("jetDrDtMed_htPct", "jetDrDtMed_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[31] = new TH2D("jetDrDtMed_dPhi", "jetDrDtMed_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);

    //hist2d[32] = new TH2D("jetDrMuTime_sceta", "jetDrMuTime_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[33] = new TH2D("jetDrMuTime_scphi","jetDrMuTime_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[34] = new TH2D("jetDrMuTime_scenr", "jetDrMuTime_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    //hist2d[35] = new TH2D("jetDrMedTime_sceta", "jetDrMedTime_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[36] = new TH2D("jetDrMedTime_scphi","jetDrMedTime_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[37] = new TH2D("jetDrMedTime_scenr", "jetDrMedTime_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    hist2d[38] = new TH2D("nRhDrJet_nRhBcJet", "nRhDrJet_nRhBcJet", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);
    hist2d[39] = new TH2D("nRhDrJet_nRhScJet", "nRhDrJet_nRhScJet", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);

    hist2d[30] = new TH2D("jetScDtMu_effje", "jetScDtMu_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);
    hist2d[31] = new TH2D("jetDrDtMu_effje", "jetDrDtMu_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);


	//---jet id stuff 50 - 99 ---------------------------------------------------
	//hist2d[50]
	//hist2d[51]
	//hist2d[52]
    hist2d[53] = new TH2D("scemf_emf", "scemf_emf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[54] = new TH2D("bcemf_emf", "bcemf_emf", 110, 0, 1.1, 110, 0, 1.1);
    //hist2d[55] = new TH2D("dreme_eme", "dreme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[56] = new TH2D("sceme_eme", "sceme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[57] = new TH2D("bceme_eme", "bceme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[58] = new TH2D("sce_bce", "sce_bce", 500, 0, 500, 500, 0, 500);
    hist2d[59] = new TH2D("dremf_scemf", "dremf_scemf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[60] = new TH2D("scemf_bcemf", "scemf_bcemf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[61] = new TH2D("epaf_epf", "epaf_epf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[62] = new TH2D("epf_emf", "epf_emf", 110, 0, 1.1, 110, 0, 1.1);

	//--- jet cluster slope 100 - 149 --------------------------------------------------
    hist2d[100] = new TH2D("jetEta_Slope", "Jet Eta v Slope;Eta;Slope ps/cm", 120, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[101] = new TH2D("jetImpAngle_Slope", "Jet ImpactAngle v Slope;ImpactAngle;Slope ps/cm", 120, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[102] = new TH2D("jetImpAngle_Slope3D", "Jet ImpactAngle v Slope 3D;ImpactAngle;Slope", 150, 0, 1.5, sl3div, sl3min, sl3max);

    //hist2d[119] = new TH2D("jetSlope_RotSlope", "Jet Slope v rotated Slope;Slope;rotated Slope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    //hist2d[120] = new TH2D("jetSlope_DifSlope", "Jet Slope v dif w/ rotSlope;Slope;difSlope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    hist2d[121] = new TH2D("jetPhi_Slope", "Jet Phi v Slope;Phi;Slope ps/cm", 70, -3.5, 3.5, sldiv, slmin, slmax);

    hist2d[126] = new TH2D("jetEta_Slope3D", "Jet Eta v Slope 3D;Eta;Slope ps/cm", 60, -1.5, 1.5, sl3div, sl3min, sl3max);
    hist2d[127] = new TH2D("jetPhi_Slope3D", "Jet Phi v Slope 3D;Phi;Slope ps/cm", 140, -3.5, 3.5, sl3div, sl3min, sl3max);

	hist2d[128] = new TH2D("jetSphEgnxy", "jetSphEgn x,y;x;y", 200, -5, 5, 200, -5, 5);

	// jet gen stuff 150 - 199 --------------------------------------------------------------
    hist2d[150] = new TH2D("jetE_GenE", "Jet Energy v GenEnergy;JetEnergy;GenEnergy", 100, 0, 1000, 100, 0, 1000 );
    hist2d[151] = new TH2D("jetEGenERatio_GenTime", "Jet E/GenE v GenTime;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[152] = new TH2D("jetEMFrac_SCTime", "Jet EMFrac v SCTime;EMFrac;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[153] = new TH2D("jetGenRatio_GenTimePre", "Jet E/GenE v GenTime Pre;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[154] = new TH2D("jetGenTime_DrJetTime", "GenTime v DrJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[155] = new TH2D("jetEGenERatio_SCTimeDiff", "Jet E/GenE v JetSC GenJet TimeDif;E/GenE;SCTimeDif", 80, 0, 2, 300, 0, 30.0 );
    hist2d[156] = new TH2D("jetSCTime_DrTime", "JetSCTime v JetDrTime;JetSCTime;JetDrTime", 280, -15, 25, 280, -15, 25 );
    hist2d[157] = new TH2D("jetGenTime_SCJetTime", "GenTime v SCJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[158] = new TH2D("jetGenTime_GenEnergy", "GenTime v GenEnergy;GenTime;GenEnergy", 280, -15, 25, 100, 0, 1000 );
    hist2d[159] = new TH2D("jetGenjetDr_SCTimeDiff", "Jet GenJet Dr v SCJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[160] = new TH2D("jetGenjetDr_DRTimeDiff", "Jet GenJet Dr v DRJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[161] = new TH2D("jetGenTime_TOFcorr", "GenTime v TOFcorr;GenTime;TOFcorr", 250, 0, 25, 250, 0, 25 );
    hist2d[162] = new TH2D("jetEGenERatio_SCTime", "Jet E/GenE v SCTime;E/GenE;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[163] = new TH2D("jetEGenERatio_DRTime", "Jet E/GenE v DRTime;E/GenE;DRTime", 80, 0, 2, 40, -15, 25 );
    hist2d[164] = new TH2D("jetGenVar_SCTimeDiff", "Jet GenJet Var v SCJet GenJet TimeDif;Var;TimeDif", 270, -2, 25, 200, 0, 20.0 );
    hist2d[165] = new TH2D("jetGenPurity_SCTimeDiff", "Jet GenJet Purity v SCJet GenJet TimeDif;Purity;TimeDif", 100, 0, 1, 200, 0, 20.0 );
    hist2d[166] = new TH2D("jetGenPurity_GenJet_ar", "Jet GenJet Purity v GenJet Var;Purity;Var", 100, 0, 1, 270, -2, 25 );
    hist2d[167] = new TH2D("jetGenVar_GenJetNKids", "Jet GenJet Var v GenJet nKids;Var;nKids", 270, -2, 25, 100, 0, 100 );
    hist2d[168] = new TH2D("jetGenPurity_GenJetNKids", "Jet GenJet Purity v GenJet nKids;Purity;nKids", 100, 0, 1, 100, 0, 100 );
    hist2d[169] = new TH2D("jetGenSCTimeDiff_DrMatchJet", "genJet SCTimeDiff v DrMatchJet;SCTimeDiff;DrMatchJet", 300, 0, 30, 320, 0, 3.2 );
    hist2d[170] = new TH2D("jetGenTime_JetEMFrac", "GenTime v JetEMFrac;GenTime;JetEMFrac", 280, -15, 25, 150, 0, 1.5 );

    hist2d[171] = new TH2D("jetEmFrac_GenjetDr", "Jet emFrac v GenJet Dr;emFrac;jetGenDr", 40, 0, 1, 200, 0, 0.5 );
    hist2d[172] = new TH2D("jetEGenERatio_GenjetDr", "Jet E/GenE v GenJet Dr;E/GenE;jetGenDr", 80, 0, 2, 200, 0, 0.5 );

    hist2d[173] = new TH2D("jet3DEgnxy", "jet3DEgn x,y;x;y", 200, -5, 5, 200, -5, 5);
    hist2d[174] = new TH2D("jet3DEgnxz", "jet3DEgn x,z;x;z", 200, -5, 5, 200, -5, 5);
*/

	//--- Photons 200 - 349 -------------------------------------------
    //hist2d[200] = new TH2D( "phoEta_TASpSlope", "Photon Eta Vs TA Slope;Eta;Slope [ns/cm]", 60, -1.5, 1.5, sldiv, slmin, slmax);
    //hist2d[201] = new TH2D( "phoEta_2dtValue", "Photon Eta Vs EiganValue2D;Eta;EiganValue",60, -1.5, 1.5, 25, 0.5, 1.0);
    //hist2d[202] = new TH2D( "phoEta_3DTAngle", "Photon Eta Vs 3D TAngle;Eta;TAngle", 60, -1.5, 1.5, 640, -3.2, 3.2);
	//hist2d[203] = new TH2D( "phoSAngle_pho3Dangle","phoRotAngle Vs pho3DAngle",140,-0.5,6.5,140,-0.5,6.5);

    hist2d[204] = new TH2D("phoID_SUSY", "phoIDvSUSY (-1 not matched);isPho;isSUSY", 3, -1, 2, 2, 0, 2);
    hist2d[205] = new TH2D("phoPreE_PrePt", "phoPreEvPrePt;Energy [GeV];Pt [GeV]", 2000, 0, 2000, 1000, 0, 1000);

    hist2d[206] = new TH2D("pho_tmap_rot", "Photon Time Map Rotated;MajAxis [cm];MinAxis [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[207] = new TH2D("pho_occmap_rot", "Photon Occ Map Rotated;MajAxis [cm];MinAxis [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[208] = new TH2D("pho_tmap", "Photon Time Map;Z [cm];C [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[209] = new TH2D("pho_occmap", "Photon Occ Map;Z [cm];C [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);

    //hist2d[210] = new TH2D("phoSphEgnxy", "phoSphEgn z,c;z;c", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[211] = new TH2D("pho3DEgnxy", "pho3DEgn x,y;x;y", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[212] = new TH2D("pho3DEgnxz", "pho3DEgn x,z;x;z", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[213] = new TH2D("pho3DEgnyz", "pho3DEgn y,z;y;z", 220, -1.1, 1.1, 220, -1.1, 1.1);

    hist2d[214] = new TH2D("pho_ptwtmap", "Photon Phi(x) Time(y) WtMap Sph;MinAxis [cm];Time [ns]", clsphdiv, -1*clsphtrn, clsphtrn, cwdiv, -1*cwtrn, cwtrn );
    hist2d[215] = new TH2D("pho_etwtmap", "Photon Eta(x) Time(y) WtMap Sph;MajAxis [cm];Time [ns]", clsphdiv, -1*clsphtrn, clsphtrn, cwdiv, -1*cwtrn, cwtrn );

    //hist2d[216] = new TH2D("pho2dSlope_3dTValue","2D Slope v 3dTValue;Slope;3dTValue",200,-2,2,400,0,0.4);
    hist2d[217] = new TH2D("phoPostE_PostPt", "phoPostEvPostPt;Energy [GeV];Pt [GeV]", 2000, 0, 2000, 1000, 0, 1000);

    hist2d[218] = new TH2D("phoNRH_ClR9","Photon nClRecHits v clR9;nRecHits;ClusterR9",80,0,80,100,0,1);
    hist2d[222] = new TH2D("phoNRH_phoClSMaj","phoNRH_phoClSMaj;#rh;SMaj",80,0,80,200,0,2);
    hist2d[223] = new TH2D("phoNRH_phoClSMin","phoNRH_phoClSMin;#rh;SMin",80,0,80,80,0,0.8);
    hist2d[224] = new TH2D("phoNRH_phoClSigIEtaIEta","phoNRH_phoClSigIEtaIEta;#rh;sieie",80,0,80,250,0,0.025);
    hist2d[219] = new TH2D("phoNRH_phoE","Photon nClRecHits v phoE;nRecHits;Energy",80,0,80,2000,0,2000);
    hist2d[220] = new TH2D("phoNRH_phoClTime","Photon nClRecHits v ClTime;nRecHits;ClTime [ns]",80,0,80,jtdiv, -1*jtran, jtran);
    hist2d[229] = new TH2D("phoNRH_geoValue","Photon nClRecHits v geoValue;nRecHits;geoValue",80,0,80,50,0.5,1);
    hist2d[228] = new TH2D("phoNRH_phoSMaj","Photon nClRecHits v phoSMaj;nRecHits;phoSMaj",80,0,80,250,0,2.5);
    hist2d[230] = new TH2D("phoNRH_phoSMin","Photon nClRecHits v phoSMin;nRecHits;phoSMin",80,0,80,750,0,0.75);
    hist2d[231] = new TH2D("phoNRH_phoSAlp","Photon nClRecHits v phoSAlp;nRecHits;phoSAlp",80,0,80,350,-1.75,1.75);
    hist2d[232] = new TH2D("phoNRH_phoCovEtaEta","Photon nClRecHits v phoCovEtaEta;nRecHits;phoCovEtaEta",80,0,80,100,0,0.001);
    hist2d[233] = new TH2D("phoNRH_phoCovEtaPhi","Photon nClRecHits v phoCovEtaPhi;nRecHits;phoCovEtaPhi",80,0,80,100,-0.0005,0.0005);
    hist2d[234] = new TH2D("phoNRH_phoCovPhiPhi","Photon nClRecHits v phoCovPhiPhi;nRecHits;phoCovPhiPhi",80,0,80,100,0,0.001);

    hist2d[247] = new TH2D("phoClR9_geoValue","Photon ClR9 v geoValue;ClR9;geoValue",100,0,1,50,0.5,1);
    hist2d[252] = new TH2D("phoClSMaj_geoValue","phoClSMaj v geoValue;ClSMaj;geoValue",200,0,2,50,0.5,1);
    hist2d[255] = new TH2D("phoClSMin_geoValue","phoClSMin v geoValue;ClSMin;geoValue",80,0,0.8,50,0.5,1);
    hist2d[258] = new TH2D("phoClSigIEtaIEta_geoValue","phoClSigIEtaIEta v geoValue;ClSigIEtaIEta;geoValue",250,0,0.025,50,0.5,1);
    hist2d[235] = new TH2D("phoE_geoValue","Photon phoE v geoValue;Energy;geoValue",2000,0,2000,50,0.5,1);
    hist2d[236] = new TH2D("phoClTime_geoValue","Photon ClTime v geoValue;ClTime [ns];geoValue",jtdiv,-1*jtran,jtran,50,0.5,1);
    hist2d[237] = new TH2D("phoSMaj_geoValue","Photon SMaj v geoValue;phoSMaj;geoValue",250,0,2.5,50,0.5,1);
    hist2d[238] = new TH2D("phoSMin_geoValue","Photon SMin v geoValue;phoSMin;geoValue",750,0,0.75,50,0.5,1);
    hist2d[239] = new TH2D("phoSAlp_geoValue","Photon SAlp v geoValue;phoSAlp;geoValue",350,-1.75,1.75,50,0.5,1);
    hist2d[240] = new TH2D("phoCovEtaEta_geoValue","Photon CovEtaEta v geoValue;phoCovEtaEta;geoValue",100,0,0.001,50,0.5,1);
    hist2d[241] = new TH2D("phoCovEtaPhi_geoValue","Photon CovEtaPhi v geoValue;phoCovEtaPhi;geoValue",100,-0.0005,0.0005,50,0.5,1);
    hist2d[242] = new TH2D("phoCovPhiPhi_geoValue","Photon CovPhiPhi v geoValue;phoCovPhiPhi;geoValue",100,0,0.001,50,0.5,1);

    hist2d[225] = new TH2D("phoClR9_phoClSMaj","phoClR9_phoClSMaj;clr9;SMaj",100,0,1,200,0,2);
    hist2d[226] = new TH2D("phoClR9_phoClSMin","phoClR9_phoClSMin;clr9;SMin",100,0,1,80,0,0.8);
    hist2d[227] = new TH2D("phoClR9_phoClSigIEtaIEta","phoClR9_phoClSigIEtaIEta;clr9;sieie",100,0,1,250,0,0.025);
    hist2d[245] = new TH2D("phoClR9_phoE","Photon phoClR9 v phoE;clr9;Energy",100,0,1,2000,0,2000);
    hist2d[246] = new TH2D("phoClR9_phoClTime","Photon phoClR9 v ClTime;clr9;ClTime [ns]",100,0,1,jtdiv, -1*jtran, jtran);
    hist2d[248] = new TH2D("phoClR9_phoSMaj","Photon phoClR9 v phoSMaj;clr9;phoSMaj",100,0,1,250,0,2.5);
    hist2d[249] = new TH2D("phoClR9_phoSMin","Photon phoClR9 v phoSMin;clr9;phoSMin",100,0,1,750,0,0.75);
    hist2d[250] = new TH2D("phoClR9_phoSAlp","Photon phoClR9 v phoSAlp;clr9;phoSAlp",100,0,1,350,-1.75,1.75);
    hist2d[251] = new TH2D("phoClR9_phoCovEtaEta","Photon phoClR9 v phoCovEtaEta;clr9;phoCovEtaEta",100,0,1,100,0,0.001);
    hist2d[252] = new TH2D("phoClR9_phoCovEtaPhi","Photon phoClR9 v phoCovEtaPhi;clr9;phoCovEtaPhi",100,0,1,100,-0.0005,0.0005);
    hist2d[253] = new TH2D("phoClR9_phoCovPhiPhi","Photon phoClR9 v phoCovPhiPhi;clr9;phoCovPhiPhi",100,0,1,100,0,0.001);

    hist2d[254] = new TH2D("phoClSMaj_phoClSMin","phoClSMaj v phoClSMin;SMaj;SMin",200,0,2,80,0,0.8);
    hist2d[255] = new TH2D("phoClSMaj_phoClSigIEtaIEta","phoClSMaj v phoClSigIEtaIEta;SMaj;sieie",200,0,2,250,0,0.025);
    hist2d[256] = new TH2D("phoClSMaj_phoE","Photon phoClSMaj v phoE;SMaj;Energy",200,0,2,2000,0,2000);
    hist2d[257] = new TH2D("phoClSMaj_phoClTime","Photon phoClSMaj v ClTime;SMaj;ClTime [ns]",200,0,2,jtdiv, -1*jtran, jtran);
    hist2d[259] = new TH2D("phoClSMaj_phoSMaj","Photon phoClSMaj v phoSMaj;SMaj;phoSMaj",200,0,2,250,0,2.5);
    hist2d[260] = new TH2D("phoClSMaj_phoSMin","Photon phoClSMaj v phoSMin;SMaj;phoSMin",200,0,2,750,0,0.75);
    hist2d[261] = new TH2D("phoClSMaj_phoSAlp","Photon phoClSMaj v phoSAlp;SMaj;phoSAlp",200,0,2,350,-1.75,1.75);
    hist2d[262] = new TH2D("phoClSMaj_phoCovEtaEta","Photon phoClSMaj v phoCovEtaEta;SMaj;phoCovEtaEta",200,0,2,100,0,0.001);
    hist2d[263] = new TH2D("phoClSMaj_phoCovEtaPhi","Photon phoClSMaj v phoCovEtaPhi;SMaj;phoCovEtaPhi",200,0,2,100,-0.0005,0.0005);
    hist2d[264] = new TH2D("phoClSMaj_phoCovPhiPhi","Photon phoClSMaj v phoCovPhiPhi;SMaj;phoCovPhiPhi",200,0,2,100,0,0.001);

    hist2d[265] = new TH2D("phoClSMin_phoClSigIEtaIEta","phoClSMin v phoClSigIEtaIEta;SMin;sieie",80,0,0.8,250,0,0.025);
    hist2d[266] = new TH2D("phoClSMin_phoE","Photon phoClSMin v phoE;SMin;Energy",80,0,0.8,2000,0,2000);
    hist2d[267] = new TH2D("phoClSMin_phoClTime","Photon phoClSMin v ClTime;SMin;ClTime [ns]",80,0,0.8,jtdiv, -1*jtran, jtran);
    hist2d[269] = new TH2D("phoClSMin_phoSMaj","Photon phoClSMin v phoSMaj;SMin;phoSMaj",80,0,0.8,250,0,2.5);
    hist2d[270] = new TH2D("phoClSMin_phoSMin","Photon phoClSMin v phoSMin;SMin;phoSMin",80,0,0.8,750,0,0.75);
    hist2d[271] = new TH2D("phoClSMin_phoSAlp","Photon phoClSMin v phoSAlp;SMin;phoSAlp",80,0,0.8,350,-1.75,1.75);
    hist2d[272] = new TH2D("phoClSMin_phoCovEtaEta","Photon phoClSMin v phoCovEtaEta;SMin;phoCovEtaEta",80,0,0.8,100,0,0.001);
    hist2d[273] = new TH2D("phoClSMin_phoCovEtaPhi","Photon phoClSMin v phoCovEtaPhi;SMin;phoCovEtaPhi",80,0,0.8,100,-0.0005,0.0005);
    hist2d[274] = new TH2D("phoClSMin_phoCovPhiPhi","Photon phoClSMin v phoCovPhiPhi;SMin;phoCovPhiPhi",80,0,0.8,100,0,0.001);

    hist2d[275] = new TH2D("phoClSigIEtaIEta_phoE","Photon phoClSigIEtaIEta v phoE;sieie;Energy",250,0,0.025,2000,0,2000);
    hist2d[276] = new TH2D("phoClSigIEtaIEta_phoClTime","Photon phoClSigIEtaIEta v ClTime;sieie;ClTime [ns]",250,0,0.025,jtdiv, -1*jtran, jtran);
    hist2d[278] = new TH2D("phoClSigIEtaIEta_phoSMaj","Photon phoClSigIEtaIEta v phoSMaj;sieie;phoSMaj",250,0,0.025,250,0,2.5);
    hist2d[279] = new TH2D("phoClSigIEtaIEta_phoSMin","Photon phoClSigIEtaIEta v phoSMin;sieie;phoSMin",250,0,0.025,750,0,0.75);
    hist2d[280] = new TH2D("phoClSigIEtaIEta_phoSAlp","Photon phoClSigIEtaIEta v phoSAlp;sieie;phoSAlp",250,0,0.025,350,-1.75,1.75);
    hist2d[281] = new TH2D("phoClSigIEtaIEta_phoCovEtaEta","Photon phoClSigIEtaIEta v phoCovEtaEta;sieie;phoCovEtaEta",250,0,0.025,100,0,0.001);
    hist2d[282] = new TH2D("phoClSigIEtaIEta_phoCovEtaPhi","Photon phoClSigIEtaIEta v phoCovEtaPhi;sieie;phoCovEtaPhi",250,0,0.025,100,-0.0005,0.0005);
    hist2d[283] = new TH2D("phoClSigIEtaIEta_phoCovPhiPhi","Photon phoClSigIEtaIEta v phoCovPhiPhi;sieie;phoCovPhiPhi",250,0,0.025,100,0,0.001);

    hist2d[284] = new TH2D("phoE_phoClTime","Photon phoE v ClTime;Energy;ClTime [ns]",2000,0,2000,jtdiv, -1*jtran, jtran);
    hist2d[286] = new TH2D("phoE_phoSMaj","Photon phoE v phoSMaj;Energy;phoSMaj",2000,0,2000,250,0,2.5);
    hist2d[287] = new TH2D("phoE_phoSMin","Photon phoE v phoSMin;Energy;phoSMin",2000,0,2000,750,0,0.75);
    hist2d[288] = new TH2D("phoE_phoSAlp","Photon phoE v phoSAlp;Energy;phoSAlp",2000,0,2000,350,-1.75,1.75);
    hist2d[289] = new TH2D("phoE_phoCovEtaEta","Photon phoE v phoCovEtaEta;Energy;phoCovEtaEta",2000,0,2000,100,0,0.001);
    hist2d[290] = new TH2D("phoE_phoCovEtaPhi","Photon phoE v phoCovEtaPhi;Energy;phoCovEtaPhi",2000,0,2000,100,-0.0005,0.0005);
    hist2d[291] = new TH2D("phoE_phoCovPhiPhi","Photon phoE v phoCovPhiPhi;Energy;phoCovPhiPhi",2000,0,2000,100,0,0.001);

    hist2d[293] = new TH2D("phoClTime_phoSMaj","Photon phoClTime v phoSMaj;ClTime [ns];phoSMaj",jtdiv, -1*jtran, jtran,250,0,2.5);
    hist2d[294] = new TH2D("phoClTime_phoSMin","Photon phoClTime v phoSMin;ClTime [ns];phoSMin",jtdiv, -1*jtran, jtran,750,0,0.75);
    hist2d[295] = new TH2D("phoClTime_phoSAlp","Photon phoClTime v phoSAlp;ClTime [ns];phoSAlp",jtdiv, -1*jtran, jtran,350,-1.75,1.75);
    hist2d[296] = new TH2D("phoClTime_phoCovEtaEta","Photon phoClTime v phoCovEtaEta;ClTime [ns];phoCovEtaEta",jtdiv, -1*jtran, jtran,100,0,0.001);
    hist2d[297] = new TH2D("phoClTime_phoCovEtaPhi","Photon phoClTime v phoCovEtaPhi;ClTime [ns];phoCovEtaPhi",jtdiv, -1*jtran, jtran,100,-0.0005,0.0005);
    hist2d[298] = new TH2D("phoClTime_phoCovPhiPhi","Photon phoClTime v phoCovPhiPhi;ClTime [ns];phoCovPhiPhi",jtdiv, -1*jtran, jtran,100,0,0.001);

/* avalible numbers -------
	237
	247
	258
	268
	277
	285
	292
------------------------*/

    //hist2d[] = new TH2D("phoClR9_minDr","Photon clR9 Vs minDr;ClusterR9;minDr",100,0,1,500,0,5);
    //hist2d[] = new TH2D("phoR9_ClR9","phoR9_ClR9;R9;clR9",100,0,1,100,0,1);

    //hist2d[] = new TH2D("phoNRH_R9_fake","Photon Fake nClRecHits v R9;nRecHits;R9",200,0,200,100,0,1);
    //hist2d[] = new TH2D("phoNRH_angle","Photon nClRecHits v Angle;nRecHits;Angle",120,0,120,70,0,7);
    //hist2d[] = new TH2D("phoNRH_","Photon nClRecHits v R9;nRecHits;R9",120,0,120,100,0,1);
    //hist2d[] = new TH2D("phoNRH_value","Photon nClRecHits v 2dtValue;nRecHits;2dtValue",120,0,120,50,0.5,1);
    //hist2d[] = new TH2D("phoNRH_slope","Photon nClRecHits v Slope;nRecHits;Slope",120,0,120,sldiv, slmin, slmax);
    //hist2d[ = new TH2D("phoNRH_geoSlope","Photon nClRecHits v geoSlope;nRecHits;geoSlope",120,0,120,sldiv, slmin, slmax);
    //hist2d[] = new TH2D("phoNRH_geoAngle","Photon nClRecHits v geoAngle;nRecHits;geoAngle",120,0,120,70,0,7);

    //hist2d[228] = new TH2D("phoAngle_phoClSMaj","phoAngle_phoClSMaj;2dtAng;Smaj",70,0,7,500,0,5);
    //hist2d[] = new TH2D("phoAngle_phoClSMin","phoAngle_phoClSMin;2dtAng;SMin",70,0,7,200,0,2);
    //hist2d[260] = new TH2D("phoAngle_phoClSigIEtaIEta","phoAngle_phoClSigIEtaIEta;2dtAng;sieie",70,0,7,500,0,0.05);

    //hist2d[261] = new TH2D("phoGeoAngle_phoClSMaj","phoGeoAngle_phoClSMaj;geoAng;SMaj",70,0,7,500,0,5);
    //hist2d[262] = new TH2D("phoGeoAngle_phoClSMin","phoGeoAngle_phoClSMin;geoAng;SMin",70,0,7,200,0,2);
    //hist2d[263] = new TH2D("phoGeoAngle_phoClSigIEtaIEta","phoGeoAngle_phoClSigIEtaIEta;geoAng;sieie",70,0,7,500,0,0.05);

	// eigan plots value/slope/angle/clr9

    //hist2d[230] = new TH2D("phoClR9_value","Photon ClR9 v tValue;ClR9;tValue",100,0,1,50,0.5,1);
    //hist2d[231] = new TH2D("phoClR9_slope","Photon ClR9 v Slope;ClR9;Slope",100,0,1,sldiv, slmin, slmax);
    //hist2d[232] = new TH2D("phoClR9_Angle","Photon ClR9 v Angle;ClR9;Angle",100,0,1,70,0,7);
    //hist2d[233] = new TH2D("phoValue_slope","Photon Value v Slope;Value;Slope",50,0.5,1,sldiv, slmin, slmax);
    //hist2d[235] = new TH2D("phoAngle_value","Photon Angle v tValue;Angle;Value",70,0,7,50,0.5,1);
    //hist2d[234] = new TH2D("phoAngle_slope","Photon Angle v Slope;Angle;Slope",70,0,7,sldiv, slmin, slmax);

    //hist2d[236] = new TH2D("phoGeoValue_angle","Photon geoValue v angle;geoValue;angle",50,0.5,1,70,0,7);
    //hist2d[237] = new TH2D("phoGeoSlope_angle","Photon geoSlope v angle;geoSlope;angle",sldiv, slmin, slmax,70,0,7);
    //hist2d[238] = new TH2D("phoGeoSlope_value","Photon geoSlope v value;geoSlope;value",sldiv, slmin, slmax,50,0.5,1);

    //hist2d[239] = new TH2D("phoGeoAngle_Angle","Photon geoAngle v Angle;geoAngle;rotAngle",70,0,7,70,0,7);
    //hist2d[242] = new TH2D("phoGeoAngle_Value","Photon geoAngle v tValue;geoAngle;tValue",70,0,7,50,0.5,1);
    //hist2d[241] = new TH2D("phoGeoAngle_slope","Photon geoAngle v Slope;geoAngle;slope",70,0,7,sldiv, slmin, slmax);
    //hist2d[243] = new TH2D("phoGeoValue_value","Photon geoValue v tValue;geoValue;tValue",50,0.5,1,50,0.5,1);
    //hist2d[244] = new TH2D("phoGeoValue_slope","Photon geoValue v slope;geoValue;slope",50,0.5,1,sldiv, slmin, slmax);
    //hist2d[240] = new TH2D("phoGeoSlope_slope","Photon geoSlope v slope;geoSlope;slope",sldiv, slmin, slmax,sldiv, slmin, slmax);

    //hist2d[245] = new TH2D("phoGeoAngle_geoValue","Photon geoAngle v geoValue;geoAngle;geoValue",70,0,7,50,0.5,1);
    //hist2d[246] = new TH2D("phoGeoAngle_geoSlope","Photon geoAngle v geoSlope;geoAngle;geoSlope",70,0,7,sldiv, slmin, slmax);
    //hist2d[248] = new TH2D("phoClR9_geoSlope","Photon ClR9 v geoSlope;ClR9;geoSlope",100,0,1,sldiv, slmin, slmax);
    //hist2d[249] = new TH2D("phoClR9_geoAngle","Photon ClR9 v geoAngle;ClR9;geoAngle",100,0,1,70,0,7);
    //hist2d[250] = new TH2D("phoGeoValue_geoSlope","Photon geoValue v geoSlope;geoValue;geoSlope",50,0.5,1,sldiv, slmin, slmax);

    //hist2d[251] = new TH2D("phoClSMaj_2dtValue","phoClSMaj_2dtValue;ClSMaj;2dtValue",500,0,5,50,0.5,1);
    //hist2d[253] = new TH2D("phoClSMaj_3d4Value","phoClSMaj_3d4Value;ClSMaj;3d4Value",500,0,5,50,0.5,1);
    //hist2d[254] = new TH2D("phoClSMin_2dtValue","phoClSMin_2dtValue;ClSMin;2dtValue",200,0,2,50,0.5,1);
    //hist2d[256] = new TH2D("phoClSMin_3d4Value","phoClSMin_3d4Value;ClSMin;3d4Value",200,0,2,50,0.5,1);
    //hist2d[257] = new TH2D("phoClSigIEtaIEta_2dtValue","phoClSigIEtaIEta_2dtValue;ClSigIEtaIEta;2dtValue",500,0,0.05,50,0.5,1);
    //hist2d[259] = new TH2D("phoClSigIEtaIEta_3d4Value","phoClSigIEtaIEta_3d4Value;ClSigIEtaIEta;3d4Value",500,0,0.05,50,0.5,1);
	//60 - 63

	//--- rechit collections 350 - 399 -------------------------------------------------
    hist2d[350] = new TH2D("ebRhTime_Energy", "ebRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );
    hist2d[351] = new TH2D("eeRhTime_Energy", "eeRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );

    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

    //hist3d[0] = new TH3D("phoNRH_ClR9_phoID","Photon nClRecHits v clR9 v phoId;nRecHits;ClusterR9;PhotonID(fake1,loose2,tight3)",200,0,200,100,0,1,8,0,4);

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------
	nMaps = 0;
	for(int it=0; it<nEBEEMaps; it++){
		fMap[it] = false;
		string label(";iEta;iPhi");
        string stt1("ebeeMapPhoCluster_"+std::to_string(it));
        ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
        string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
        ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
		string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
        ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
	}//<<>>for(int it=0; it<nEBEEMaps; it++)

}//<<>>void HistMaker::initHists()

