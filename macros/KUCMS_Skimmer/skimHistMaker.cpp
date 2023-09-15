//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "skimHistMaker.hh"

//#define DEBUG true
#define DEBUG false
#define doEBEEmaps false

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// HistMaker class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

//HistMaker::HistMaker(){}

//HistMaker::~HistMaker(){}

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("kuSkimTree");
    //const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
    const std::string eosdir("");
    const std::string listdir("");

	std::cout << "Producing Histograms for : " << outfilename << std::endl;
    std::ifstream infile(listdir+infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
		//std::cout << "--  for Fine #" << cnt << " moduls " << cnt%pct << " ";
		//if( cnt%pct == 0 ){ 
        auto tfilename = eosdir + indir + str;
        std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
		//}//<<>>if( cnt%4 == 0 ){
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
   
    if(debug) std::cout << " - Creating output file " << std::endl;
    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

	endJobs();
	for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }

	nMaps = 0;
	if( doEBEEmaps ){ for( int it = 0; it < nEBEEMaps; it++ ){ 
		ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
		ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
		ebeeMapR[it]->Write(); delete ebeeMapR[it];
	}}//<<>>for( int it = 0; it < nEBEEMaps; it++ )

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

/*
    if( false ) { // rechits lock
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
    }//if( false ) { // rechits lock
*/
/*
    if( false ) { // genparts lock
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
*/

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------
	
    if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( int it = 0; it < nSelPhotons; it++ ){

/*
		// detrimine photon classification ( sig/susy/ect... )
		//---------------------------------------------------
		auto isCmb = not (*phoExcluded)[it];
        auto isCmb = phoExcluded[it];
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
*/

		// determine pog id class
		// -----------------------------------------------------

		auto phoClass = (*selPhoQuality)[it];
		hist1d[122]->Fill( phoClass );//122

        if( DEBUG ) std::cout << " -- looping photons : filling ids" << std::endl;
        //hist1d[107]->Fill( (*phoIsPixelSeed)[it] );
        //hist1d[108]->Fill( (*phoHadOverEM)[it] );
        hist1d[109]->Fill( (*selPhoR9)[it] );
        hist1d[128]->Fill( (*selPhoR9)[it] );//
        //hist1d[110]->Fill( (*phoMipTotEnergy)[it] );
        //hist1d[111]->Fill( (*phoEcalRHSumEtConeDR04)[it] );
        //hist1d[112]->Fill( (*phoHcalTwrSumEtConeDR04)[it] );

		//find intial cut values 
		//----------------------------------------------------------------
		if( DEBUG ) std::cout << " -- pho rechits" << std::endl;
		auto nrh = (*selPhoNrh)[it];
        auto phoClstrR9 = (*selPhoClstrRn)[it]; //clstrR9( (*phoRhIds)[it] );

        //auto isEB = (*phoIsEB)[it];
		auto isTightEB = std::abs((*selPhoEta)[it]) < 1.45;

		auto isTight = phoClass == 3;
        auto isLoose = phoClass > 1;
        auto isFake = phoClass == 1;

        auto goodRhCnt = nrh > 14;
        auto minPhoPt = (*selPhoPt)[it] > 30;

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

        auto usePho = isLoose && isTightEB;
        //auto usePho = isTight && isEB && isTightEB && isCmb && isSigPho;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isSigPho;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && not isSigPho;


		//if( DEBUG ) std::cout << " -- photons # " << it << " isSUSY " << isSUSY << " isOOT " << isOOT << " isCmb " << isCmb << std::endl;
        if( DEBUG ) std::cout << " -- sel photons # " << it << std::endl;
		if( usePho && goodRhCnt && minPhoPt ) { //----------------------------------- 
	
	        hist1d[132]->Fill(1);
	        //if( isSUSY ) hist1d[132]->Fill(2);
	        //if( isOOT ) hist1d[132]->Fill(3);
	        //if( isCmb ) hist1d[132]->Fill(4);
	
			if( DEBUG ) std::cout << " -- pho time" << std::endl;
            hist1d[104]->Fill( (*selPhoPt)[it] );
            hist1d[105]->Fill( (*selPhoEnergy)[it] );
			hist2d[205]->Fill( (*selPhoEnergy)[it], (*selPhoPt)[it] );

	        hist1d[119]->Fill( (*selPhoTime)[it] );//c mean
	        //hist1d[120]->Fill( (*phoSeedTOFTime)[it] );//lead time 
	        //hist1d[121]->Fill( (*selPhoTime)[it] - (*phoSeedTOFTime)[it] );//diff
			//for( auto rhid : (*phoRhIds)[it] ){ hist1d[115]->Fill((*rhEnergy)[getRhIdx(rhid)]);}

            //auto isPho = (phoGenIdx >= 0)?(((*genPdgId)[phoGenIdx] == 22)?1:0):-1;
            //hist2d[204]->Fill(isPho,isSUSY);

			if( DEBUG ) std::cout << " Not Finding Eigans for Photon W/ " << nrh << " rechits." << std::endl;

        }//<<>>if( usePho ) { //----------------------------------- 

    }//<<>>for( int it = 0; it < nPhotons; it++ )

	//-------- electrons --------------------------------------

/*
    for( int it = 0; it < nElectrons; it++ ){

        hist1d[350]->Fill( (*eleCMeanTime)[it] );//c mean
        hist1d[351]->Fill( (*eleSeedTOFTime)[it] );//lead time 
        hist1d[352]->Fill( (*eleCMeanTime)[it] - (*eleSeedTOFTime)[it] );//diff

    }//<<>>for( int it = 0; it < nElectrons; it++ )
*/

	//--------- jets --------------------------
    if( DEBUG ) std::cout << "Finding Jets with " << nSelJets << " selected. "<< std::endl;

	//// ---  base info  ----

	hist1d[9]->Fill(nSelJets);

    for( int it = 0; it < nSelJets; it++ ){

        //const auto jetepafrac   = (*jetPHEF)[it] + (*jetELEF)[it];;
        //const auto jetepe       = (*jetPHE)[it] + (*jetELE)[it];
        //const auto jeteme       = (*jetCEMF)[it] + (*jetNEMF)[it];
        //const auto jetemfrac    = jeteme/(*jetE)[it];
        //const auto jetepfrac    = jetepe/(*jetE)[it];
		
        //hist2d[61]->Fill(jetepafrac,jetepfrac);
        //hist2d[62]->Fill(jetepfrac,jetemfrac);

        fillTH1((*selJetPt)[it],hist1d[0]);//hist1d[12]->Fill(jet.pt());
        fillTH1((*selJetPhi)[it],hist1d[1]);//hist1d[13]->Fill(jet.phi());
        fillTH1((*selJetEta)[it],hist1d[2]);//hist1d[14]->Fill(jet.eta());

	}//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding genjet" << std::endl;
	//// --- genjet info -----------------------------------
/*
    for( int it = 0; it < nJets; it++ ){

        hist2d[161]->Fill((*jetGenTime)[it],(*jetGenTOF)[it]);
        hist1d[50]->Fill((*jetGenTime)[it]);
    	hist1d[51]->Fill((*jetGenTOF)[it]);

    }//<<>>for( int it = 0; it < nJets; it++ )
*/

}//<<>>void HistMaker::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

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

    //thresDivTH2D( hist2d[206], hist2d[207], 0 );
    //thresDivTH2D( hist2d[208], hist2d[209], 0 );

    //profileTH2D( hist2d[215], hist1d[116], hist1d[118] );


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


	//------- jets 0 - 99
    hist1d[0] = new TH1D("jetPt", "jetPt", 500, 0, 5000);
    hist1d[1] = new TH1D("jetPhi", "jetPhi", 700, -3.5, 3.5);
    hist1d[2] = new TH1D("jetEta", "jetEta", 700, -3.5, 3.5);
    hist1d[9] = new TH1D("nJet", "nJets", 21, -0.5, 20.5);

	//----- photons 100 - 249

    //hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);
    //hist1d[101] = new TH1D("genPhoOriginCode", "genPhoOriginCode;Code",1000,0,1000);

    //hist1d[102] = new TH1D("phoClstrR9", "phoClstrR9", 100, 0, 1);

    //hist1d[103] = new TH1D("phoNumRecHits", "phoNumRecHits", 80, 0, 80);
    hist1d[104] = new TH1D("phoPtPreCut", "phoPt PreCut", 1000, 0, 1000);
    hist1d[105] = new TH1D("phoEnergyPreCut", "phoEnergy PreCut", 2000, 0, 2000);

    //hist1d[106] = new TH1D("pho3DRotAngle", "pho3DRotAngle", 140, -0.5, 6.5);

    //hist1d[107] = new TH1D("phoIsPixelSeed", "phoIsPixelSeed", 3, 0, 2);
    //hist1d[108] = new TH1D("phoHadOverEM", "phoHadOverEM", 250, 0, 10);
    hist1d[109] = new TH1D("phoR9", "phoR9", 100, 0, 1);
    //hist1d[110] = new TH1D("phoMipTotEnergy", "phoMipTotEnergy", 250, 0, 250 );
    //hist1d[111] = new TH1D("phoEcalRHSumEtConeDR04", "phoEcalRHSumEtConeDR04", 750, 0, 750);
    //hist1d[112] = new TH1D("phoHcalTwrSumEtConeDR04", "phoHcalTwrSumEtConeDR04", 750, 0, 750);

    //hist1d[113] = new TH1D("phoEginValueSph", "phoEginValueSph", 150, 0.4, 1.1);
    //hist1d[114] = new TH1D("phoEginSlopeAP", "phoEginSlopeSph 1Way", sldiv, slmin, slmax);
    //hist1d[115] = new TH1D("clRhEnergy", "clRhEnergy", 3000, 0, 1000);

    //hist1d[116] = new TH1D("pho_etprofile", "Photon rEta Time Profile Sph;MajAxis [cm]", clsphdiv, -1*clsphtrn, clsphtrn);

    //hist1d[117] = new TH1D("phoEginSlopeSphTA", "phoEginSlopeSph TA", sldiv, slmin, slmax);
	//hist1d[118] = new TH1D("pho_proFitEtavChi", "Profile Fit Eta v Chi2Prob Sph", cwdiv, -1*cwtrn, cwtrn );

    hist1d[119] = new TH1D("phoClTimePreCut", "phoTime PreCut", jtdiv, -1*jtran, jtran);
    //hist1d[120] = new TH1D("phoSeedRhTime", "phoLeadRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[121] = new TH1D("phoSeedTimeDiff", "phoLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[122] = new TH1D("phoClass", "phoClass", 4, 0, 4);

    //hist1d[123] = new TH1D("phoClRhTime","phoClRhTime",1000, -50, 50);
    //hist1d[124] = new TH1D("phoClSMaj","phoClSMaj",50,0,5);
    //hist1d[125] = new TH1D("phoClSMin","phoClSMin",50,0,0.5);
    //hist1d[126] = new TH1D("phoSigIEtaIEta","phoSigIEtaIEta",300,0,0.03);

    //hist1d[127] = new TH1D("phoRotAngle","phoRotAngle",140,-0.5,6.5);

    hist1d[128] = new TH1D("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
    //hist1d[129] = new TH1D("phoGeoEginValueSph", "phoGeoEginValueSph", 150, 0.4, 1.1);
    //hist1d[130] = new TH1D("phoGenDr", "phoGenDr", 800, 0, 2.0 );
    //hist1d[131] = new TH1D("phoGenPdgId", "phoGenPdgId", 50, 0, 50 );
    hist1d[132] = new TH1D("phoType", "phoType(1-pho 2-susy 3-oot 4-cmb)", 6,0,5);

    //hist1d[133] = new TH1D("phoGeoValue", "phoGeoValue", 150, 0.4, 1.1);
    //hist1d[134] = new TH1D("phoGeoSlope", "phoGeoSlope", sldiv, slmin, slmax);
    //hist1d[135] = new TH1D("phoGeoAngle","phoGeoAngle",70,0,7);

    //hist1d[136] = new TH1D("pho2DEiganV0","pho2DEiganV0",200,-1,1);
    //hist1d[137] = new TH1D("pho2DEiganV1","pho2DEiganV1",200,-1,1);
    //hist1d[138] = new TH1D("phoG2DEiganV0","phoG2DEiganV0",200,-1,1);
    //hist1d[139] = new TH1D("phoG2DEiganV1","phoG2DEiganV1",200,-1,1);

    //hist1d[140] = new TH1D("phoMinDr","phoMinDr",500,0,5);

    //hist1d[141] = new TH1D("pho3DEiganV0","pho3DEiganV0",200,-1,1);
    //hist1d[142] = new TH1D("pho3DEiganV1","pho3DEiganV1",200,-1,1);
    //hist1d[143] = new TH1D("pho3DEiganV2","pho3DEiganV2",200,-1,1);
    //hist1d[144] = new TH1D("pho3d4Value","pho3d4Value",100,0,1);
    //hist1d[145] = new TH1D("pho3dTValue","pho3dTValue",400,0,0.4);

    //hist1d[146] = new TH1D("phoPtPostCut", "phoPt PostCut", 1000, 0, 1000);
    //hist1d[147] = new TH1D("phoEnergyPostCut", "phoEnergy PostCut", 2000, 0, 2000);
    //hist1d[148] = new TH1D("phoClTimePostCut", "phoTime PostCut", jtdiv, -1*jtran, jtran);

    //hist1d[149] = new TH1D("pho3dEV1","pho3dEV1",20,0,20);
    //hist1d[150] = new TH1D("pho3dEV2","pho3dEV2",20,0,20);
    //hist1d[151] = new TH1D("pho3dEV3","pho3dEV3",20,0,20);

    //hist1d[152] = new TH1D("phoSMaj", "phoSMaj", 300,0,3);
    //hist1d[153] = new TH1D("phoSMin", "phoSMin", 60,0,.6);
    //hist1d[154] = new TH1D("phoSAlp", "phoSAlp", 350,-1.75,1.75);
    //hist1d[155] = new TH1D("phoCovEtaEta", "phoCovEtaEtaSqrt", 100,0,0.001);
    //hist1d[156] = new TH1D("phoCovEtaPhi", "phoCovEtaPhiSqrt", 800,-0.0004,0.0004);
    //hist1d[157] = new TH1D("phoCovPhiPhi", "phoCovPhiPhiSqrt", 100,0,0.001);

	//------  genparticles 200 - 249

    //hist1d[200] = new TH1D("genPdgId", "genPdgId;PdgId",220,0,220);
    //hist1d[201] = new TH1D("nGenLLPPho", "nGenLLPPho", 5, 0, 5 );

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
    //hist1d[300] = new TH1D("rhRadius", "rhRadius", 300, 128, 131);

    //hist1d[301] = new TH1D("ebRhTime", "ebRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    //hist1d[302] = new TH1D("ebRhEnergy", "ebRhEnergy", 1000, 0, 1000);
    //hist1d[303] = new TH1D("eeRhTime", "eeRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    //hist1d[304] = new TH1D("eeRhEnergy", "eeRhEnergy", 1000, 0, 1000);
    //hist1d[305] = new TH1D("ebRhkOOT", "ebRhkOOT", 3, 0, 1);
    //hist1d[306] = new TH1D("eeRhkOOT", "eeRhkOOT", 3, 0, 1);

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

	//---jet id stuff 50 - 99 ---------------------------------------------------

*/

	//--- Photons 200 - 349 -------------------------------------------

    //hist2d[204] = new TH2D("phoID_SUSY", "phoIDvSUSY (-1 not matched);isPho;isSUSY", 3, -1, 2, 2, 0, 2);
    hist2d[205] = new TH2D("phoPreE_PrePt", "phoPreEvPrePt;Energy [GeV];Pt [GeV]", 2000, 0, 2000, 1000, 0, 1000);

    //hist2d[206] = new TH2D("pho_tmap_rot", "Photon Time Map Rotated;MajAxis [cm];MinAxis [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    //hist2d[207] = new TH2D("pho_occmap_rot", "Photon Occ Map Rotated;MajAxis [cm];MinAxis [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    //hist2d[208] = new TH2D("pho_tmap", "Photon Time Map;Z [cm];C [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    //hist2d[209] = new TH2D("pho_occmap", "Photon Occ Map;Z [cm];C [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);

    //hist2d[210] = new TH2D("phoSphEgnxy", "phoSphEgn z,c;z;c", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[211] = new TH2D("pho3DEgnxy", "pho3DEgn x,y;x;y", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[212] = new TH2D("pho3DEgnxz", "pho3DEgn x,z;x;z", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[213] = new TH2D("pho3DEgnyz", "pho3DEgn y,z;y;z", 220, -1.1, 1.1, 220, -1.1, 1.1);

    //hist2d[214] = new TH2D("pho_ptwtmap", "Photon Phi(x) Time(y) WtMap Sph;MinAxis [cm];Time [ns]", clsphdiv, -1*clsphtrn, clsphtrn, cwdiv, -1*cwtrn, cwtrn );
    //hist2d[215] = new TH2D("pho_etwtmap", "Photon Eta(x) Time(y) WtMap Sph;MajAxis [cm];Time [ns]", clsphdiv, -1*clsphtrn, clsphtrn, cwdiv, -1*cwtrn, cwtrn );

    //hist2d[216] = new TH2D("pho2dSlope_3dTValue","2D Slope v 3dTValue;Slope;3dTValue",200,-2,2,400,0,0.4);
    //hist2d[217] = new TH2D("phoPostE_PostPt", "phoPostEvPostPt;Energy [GeV];Pt [GeV]", 2000, 0, 2000, 1000, 0, 1000);

    //hist2d[218] = new TH2D("phoNRH_ClR9","Photon nClRecHits v clR9;nRecHits;ClusterR9",80,0,80,100,0,1);
    //hist2d[222] = new TH2D("phoNRH_phoClSMaj","phoNRH_phoClSMaj;#rh;SMaj",80,0,80,200,0,2);
    //hist2d[223] = new TH2D("phoNRH_phoClSMin","phoNRH_phoClSMin;#rh;SMin",80,0,80,80,0,0.8);
    //hist2d[224] = new TH2D("phoNRH_phoClSigIEtaIEta","phoNRH_phoClSigIEtaIEta;#rh;sieie",80,0,80,250,0,0.025);
    //hist2d[219] = new TH2D("phoNRH_phoE","Photon nClRecHits v phoE;nRecHits;Energy",80,0,80,2000,0,2000);
    //hist2d[220] = new TH2D("phoNRH_phoClTime","Photon nClRecHits v ClTime;nRecHits;ClTime [ns]",80,0,80,jtdiv, -1*jtran, jtran);
    //hist2d[229] = new TH2D("phoNRH_geoValue","Photon nClRecHits v geoValue;nRecHits;geoValue",80,0,80,50,0.5,1);
    //hist2d[228] = new TH2D("phoNRH_phoSMaj","Photon nClRecHits v phoSMaj;nRecHits;phoSMaj",80,0,80,250,0,2.5);
    //hist2d[230] = new TH2D("phoNRH_phoSMin","Photon nClRecHits v phoSMin;nRecHits;phoSMin",80,0,80,750,0,0.75);
    //hist2d[231] = new TH2D("phoNRH_phoSAlp","Photon nClRecHits v phoSAlp;nRecHits;phoSAlp",80,0,80,350,-1.75,1.75);
    //hist2d[232] = new TH2D("phoNRH_phoCovEtaEta","Photon nClRecHits v phoCovEtaEta;nRecHits;phoCovEtaEta",80,0,80,100,0,0.001);
    //hist2d[233] = new TH2D("phoNRH_phoCovEtaPhi","Photon nClRecHits v phoCovEtaPhi;nRecHits;phoCovEtaPhi",80,0,80,100,-0.0005,0.0005);
    //hist2d[234] = new TH2D("phoNRH_phoCovPhiPhi","Photon nClRecHits v phoCovPhiPhi;nRecHits;phoCovPhiPhi",80,0,80,100,0,0.001);

    //hist2d[247] = new TH2D("phoClR9_geoValue","Photon ClR9 v geoValue;ClR9;geoValue",100,0,1,50,0.5,1);
    //hist2d[252] = new TH2D("phoClSMaj_geoValue","phoClSMaj v geoValue;ClSMaj;geoValue",200,0,2,50,0.5,1);
    //hist2d[255] = new TH2D("phoClSMin_geoValue","phoClSMin v geoValue;ClSMin;geoValue",80,0,0.8,50,0.5,1);
    //hist2d[258] = new TH2D("phoClSigIEtaIEta_geoValue","phoClSigIEtaIEta v geoValue;ClSigIEtaIEta;geoValue",250,0,0.025,50,0.5,1);
    //hist2d[235] = new TH2D("phoE_geoValue","Photon phoE v geoValue;Energy;geoValue",2000,0,2000,50,0.5,1);
    //hist2d[236] = new TH2D("phoClTime_geoValue","Photon ClTime v geoValue;ClTime [ns];geoValue",jtdiv,-1*jtran,jtran,50,0.5,1);
    //hist2d[237] = new TH2D("phoSMaj_geoValue","Photon SMaj v geoValue;phoSMaj;geoValue",250,0,2.5,50,0.5,1);
    //hist2d[238] = new TH2D("phoSMin_geoValue","Photon SMin v geoValue;phoSMin;geoValue",750,0,0.75,50,0.5,1);
    //hist2d[239] = new TH2D("phoSAlp_geoValue","Photon SAlp v geoValue;phoSAlp;geoValue",350,-1.75,1.75,50,0.5,1);
    //hist2d[240] = new TH2D("phoCovEtaEta_geoValue","Photon CovEtaEta v geoValue;phoCovEtaEta;geoValue",100,0,0.001,50,0.5,1);
    //hist2d[241] = new TH2D("phoCovEtaPhi_geoValue","Photon CovEtaPhi v geoValue;phoCovEtaPhi;geoValue",100,-0.0005,0.0005,50,0.5,1);
    //hist2d[242] = new TH2D("phoCovPhiPhi_geoValue","Photon CovPhiPhi v geoValue;phoCovPhiPhi;geoValue",100,0,0.001,50,0.5,1);

    //hist2d[225] = new TH2D("phoClR9_phoClSMaj","phoClR9_phoClSMaj;clr9;SMaj",100,0,1,200,0,2);
    //hist2d[226] = new TH2D("phoClR9_phoClSMin","phoClR9_phoClSMin;clr9;SMin",100,0,1,80,0,0.8);
    //hist2d[227] = new TH2D("phoClR9_phoClSigIEtaIEta","phoClR9_phoClSigIEtaIEta;clr9;sieie",100,0,1,250,0,0.025);
    //hist2d[245] = new TH2D("phoClR9_phoE","Photon phoClR9 v phoE;clr9;Energy",100,0,1,2000,0,2000);
    //hist2d[246] = new TH2D("phoClR9_phoClTime","Photon phoClR9 v ClTime;clr9;ClTime [ns]",100,0,1,jtdiv, -1*jtran, jtran);
    //hist2d[248] = new TH2D("phoClR9_phoSMaj","Photon phoClR9 v phoSMaj;clr9;phoSMaj",100,0,1,250,0,2.5);
    //hist2d[249] = new TH2D("phoClR9_phoSMin","Photon phoClR9 v phoSMin;clr9;phoSMin",100,0,1,750,0,0.75);
    //hist2d[250] = new TH2D("phoClR9_phoSAlp","Photon phoClR9 v phoSAlp;clr9;phoSAlp",100,0,1,350,-1.75,1.75);
    //hist2d[251] = new TH2D("phoClR9_phoCovEtaEta","Photon phoClR9 v phoCovEtaEta;clr9;phoCovEtaEta",100,0,1,100,0,0.001);
    //hist2d[252] = new TH2D("phoClR9_phoCovEtaPhi","Photon phoClR9 v phoCovEtaPhi;clr9;phoCovEtaPhi",100,0,1,100,-0.0005,0.0005);
    //hist2d[253] = new TH2D("phoClR9_phoCovPhiPhi","Photon phoClR9 v phoCovPhiPhi;clr9;phoCovPhiPhi",100,0,1,100,0,0.001);

    //hist2d[254] = new TH2D("phoClSMaj_phoClSMin","phoClSMaj v phoClSMin;SMaj;SMin",200,0,2,80,0,0.8);
    //hist2d[255] = new TH2D("phoClSMaj_phoClSigIEtaIEta","phoClSMaj v phoClSigIEtaIEta;SMaj;sieie",200,0,2,250,0,0.025);
    //hist2d[256] = new TH2D("phoClSMaj_phoE","Photon phoClSMaj v phoE;SMaj;Energy",200,0,2,2000,0,2000);
    //hist2d[257] = new TH2D("phoClSMaj_phoClTime","Photon phoClSMaj v ClTime;SMaj;ClTime [ns]",200,0,2,jtdiv, -1*jtran, jtran);
    //hist2d[259] = new TH2D("phoClSMaj_phoSMaj","Photon phoClSMaj v phoSMaj;SMaj;phoSMaj",200,0,2,250,0,2.5);
    //hist2d[260] = new TH2D("phoClSMaj_phoSMin","Photon phoClSMaj v phoSMin;SMaj;phoSMin",200,0,2,750,0,0.75);
    //hist2d[261] = new TH2D("phoClSMaj_phoSAlp","Photon phoClSMaj v phoSAlp;SMaj;phoSAlp",200,0,2,350,-1.75,1.75);
    //hist2d[262] = new TH2D("phoClSMaj_phoCovEtaEta","Photon phoClSMaj v phoCovEtaEta;SMaj;phoCovEtaEta",200,0,2,100,0,0.001);
    //hist2d[263] = new TH2D("phoClSMaj_phoCovEtaPhi","Photon phoClSMaj v phoCovEtaPhi;SMaj;phoCovEtaPhi",200,0,2,100,-0.0005,0.0005);
    //hist2d[264] = new TH2D("phoClSMaj_phoCovPhiPhi","Photon phoClSMaj v phoCovPhiPhi;SMaj;phoCovPhiPhi",200,0,2,100,0,0.001);

    //hist2d[265] = new TH2D("phoClSMin_phoClSigIEtaIEta","phoClSMin v phoClSigIEtaIEta;SMin;sieie",80,0,0.8,250,0,0.025);
    //hist2d[266] = new TH2D("phoClSMin_phoE","Photon phoClSMin v phoE;SMin;Energy",80,0,0.8,2000,0,2000);
    //hist2d[267] = new TH2D("phoClSMin_phoClTime","Photon phoClSMin v ClTime;SMin;ClTime [ns]",80,0,0.8,jtdiv, -1*jtran, jtran);
    //hist2d[269] = new TH2D("phoClSMin_phoSMaj","Photon phoClSMin v phoSMaj;SMin;phoSMaj",80,0,0.8,250,0,2.5);
    //hist2d[270] = new TH2D("phoClSMin_phoSMin","Photon phoClSMin v phoSMin;SMin;phoSMin",80,0,0.8,750,0,0.75);
    //hist2d[271] = new TH2D("phoClSMin_phoSAlp","Photon phoClSMin v phoSAlp;SMin;phoSAlp",80,0,0.8,350,-1.75,1.75);
    //hist2d[272] = new TH2D("phoClSMin_phoCovEtaEta","Photon phoClSMin v phoCovEtaEta;SMin;phoCovEtaEta",80,0,0.8,100,0,0.001);
    //hist2d[273] = new TH2D("phoClSMin_phoCovEtaPhi","Photon phoClSMin v phoCovEtaPhi;SMin;phoCovEtaPhi",80,0,0.8,100,-0.0005,0.0005);
    //hist2d[274] = new TH2D("phoClSMin_phoCovPhiPhi","Photon phoClSMin v phoCovPhiPhi;SMin;phoCovPhiPhi",80,0,0.8,100,0,0.001);

    //hist2d[275] = new TH2D("phoClSigIEtaIEta_phoE","Photon phoClSigIEtaIEta v phoE;sieie;Energy",250,0,0.025,2000,0,2000);
    //hist2d[276] = new TH2D("phoClSigIEtaIEta_phoClTime","Photon phoClSigIEtaIEta v ClTime;sieie;ClTime [ns]",250,0,0.025,jtdiv, -1*jtran, jtran);
    //hist2d[278] = new TH2D("phoClSigIEtaIEta_phoSMaj","Photon phoClSigIEtaIEta v phoSMaj;sieie;phoSMaj",250,0,0.025,250,0,2.5);
    //hist2d[279] = new TH2D("phoClSigIEtaIEta_phoSMin","Photon phoClSigIEtaIEta v phoSMin;sieie;phoSMin",250,0,0.025,750,0,0.75);
    //hist2d[280] = new TH2D("phoClSigIEtaIEta_phoSAlp","Photon phoClSigIEtaIEta v phoSAlp;sieie;phoSAlp",250,0,0.025,350,-1.75,1.75);
    //hist2d[281] = new TH2D("phoClSigIEtaIEta_phoCovEtaEta","Photon phoClSigIEtaIEta v phoCovEtaEta;sieie;phoCovEtaEta",250,0,0.025,100,0,0.001);
    //hist2d[282] = new TH2D("phoClSigIEtaIEta_phoCovEtaPhi","Photon phoClSigIEtaIEta v phoCovEtaPhi;sieie;phoCovEtaPhi",250,0,0.025,100,-0.0005,0.0005);
    //hist2d[283] = new TH2D("phoClSigIEtaIEta_phoCovPhiPhi","Photon phoClSigIEtaIEta v phoCovPhiPhi;sieie;phoCovPhiPhi",250,0,0.025,100,0,0.001);

    //hist2d[284] = new TH2D("phoE_phoClTime","Photon phoE v ClTime;Energy;ClTime [ns]",2000,0,2000,jtdiv, -1*jtran, jtran);
    //hist2d[286] = new TH2D("phoE_phoSMaj","Photon phoE v phoSMaj;Energy;phoSMaj",2000,0,2000,250,0,2.5);
    //hist2d[287] = new TH2D("phoE_phoSMin","Photon phoE v phoSMin;Energy;phoSMin",2000,0,2000,750,0,0.75);
    //hist2d[288] = new TH2D("phoE_phoSAlp","Photon phoE v phoSAlp;Energy;phoSAlp",2000,0,2000,350,-1.75,1.75);
    //hist2d[289] = new TH2D("phoE_phoCovEtaEta","Photon phoE v phoCovEtaEta;Energy;phoCovEtaEta",2000,0,2000,100,0,0.001);
    //hist2d[290] = new TH2D("phoE_phoCovEtaPhi","Photon phoE v phoCovEtaPhi;Energy;phoCovEtaPhi",2000,0,2000,100,-0.0005,0.0005);
    //hist2d[291] = new TH2D("phoE_phoCovPhiPhi","Photon phoE v phoCovPhiPhi;Energy;phoCovPhiPhi",2000,0,2000,100,0,0.001);

    //hist2d[293] = new TH2D("phoClTime_phoSMaj","Photon phoClTime v phoSMaj;ClTime [ns];phoSMaj",jtdiv, -1*jtran, jtran,250,0,2.5);
    //hist2d[294] = new TH2D("phoClTime_phoSMin","Photon phoClTime v phoSMin;ClTime [ns];phoSMin",jtdiv, -1*jtran, jtran,750,0,0.75);
    //hist2d[295] = new TH2D("phoClTime_phoSAlp","Photon phoClTime v phoSAlp;ClTime [ns];phoSAlp",jtdiv, -1*jtran, jtran,350,-1.75,1.75);
    //hist2d[296] = new TH2D("phoClTime_phoCovEtaEta","Photon phoClTime v phoCovEtaEta;ClTime [ns];phoCovEtaEta",jtdiv, -1*jtran, jtran,100,0,0.001);
    //hist2d[297] = new TH2D("phoClTime_phoCovEtaPhi","Photon phoClTime v phoCovEtaPhi;ClTime [ns];phoCovEtaPhi",jtdiv, -1*jtran, jtran,100,-0.0005,0.0005);
    //hist2d[298] = new TH2D("phoClTime_phoCovPhiPhi","Photon phoClTime v phoCovPhiPhi;ClTime [ns];phoCovPhiPhi",jtdiv, -1*jtran, jtran,100,0,0.001);

/* avalible numbers -------
	237
	247
	258
	268
	277
	285
	292
------------------------*/

	//60 - 63

	//--- rechit collections 350 - 399 -------------------------------------------------
    //hist2d[350] = new TH2D("ebRhTime_Energy", "ebRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );
    //hist2d[351] = new TH2D("eeRhTime_Energy", "eeRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );

    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------
	nMaps = 0;
    std::string label("baseHists");
	if( doEBEEmaps ){ for(int it=0; it<nEBEEMaps; it++){
		fMap[it] = false;
		std::string label(";iEta;iPhi");
        std::string stt1("ebeeMapPhoCluster_"+std::to_string(it));
        ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
        std::string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
        ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
		std::string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
        ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
	}}//<<>>for(int it=0; it<nEBEEMaps; it++)

}//<<>>void HistMaker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "skim_v6_files/";

                auto infilename = "KUCMS_GMSB_Skim_List.txt";

                auto outfilename = "KUCMS_GMSB_Skim_BaseHists.root";

                HistMaker base;
                base.histMaker( listdir, infilename, outfilename );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

