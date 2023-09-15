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

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle, int cut ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");

    cutselection = cut;

	std::cout << "Producing Histograms for : " << outfilename << std::endl;
    std::ifstream infile(listdir+infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    auto fConfigTree = new TChain(configtreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
        if( str[0] == '#' ) continue;
        auto tfilename = eosdir + indir + str;
        std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
        fConfigTree->Add(tfilename.c_str());
		cnt++;
    }//<<>>while (std::getline(infile,str))

    std::cout << "Setting up For Main Loop." << std::endl;

	Init(fInTree);
	initHists(htitle);

    std::cout << "Filling Config Map." << std::endl;

    UInt_t          nEvents;
    UInt_t          nSelectedEvents;
    string          *sKey;
    Float_t         sCrossSection;
    Float_t         sGMSBGravMass;
    Float_t         sGMSBChi1Mass;
    Float_t         sMCWgt;
    Int_t           sMCType;

    TBranch        *b_nEvents;   //!
    TBranch        *b_nSelectedEvents;   //!
    TBranch        *b_sKey;   //!
    TBranch        *b_sCrossSection;   //!
    TBranch        *b_sGMSBGravMass;   //!
    TBranch        *b_sGMSBChi1Mass;   //!
    TBranch        *b_sMCWgt;   //!
    TBranch        *b_sMCType;   //!

    sKey = 0;

    fConfigTree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
    fConfigTree->SetBranchAddress("nSelectedEvents", &nSelectedEvents, &b_nSelectedEvents);
    fConfigTree->SetBranchAddress("sKey", &sKey, &b_sKey);
    fConfigTree->SetBranchAddress("sCrossSection", &sCrossSection, &b_sCrossSection);
    fConfigTree->SetBranchAddress("sGMSBGravMass", &sGMSBGravMass, &b_sGMSBGravMass);
    fConfigTree->SetBranchAddress("sGMSBChi1Mass", &sGMSBChi1Mass, &b_sGMSBChi1Mass);
    fConfigTree->SetBranchAddress("sMCWgt", &sMCWgt, &b_sMCWgt);
    fConfigTree->SetBranchAddress("sMCType", &sMCType, &b_sMCType);

    auto nConfigEntries = fConfigTree->GetEntries();
    std::cout << "Proccessing " << nConfigEntries << " config entries." << std::endl;
    for (Long64_t centry = 0; centry < nConfigEntries; centry++){
		
		auto entry = fConfigTree->LoadTree(centry);

        if(debug) std::cout << " - Getting Branches. " << std::endl;
    	b_nEvents->GetEntry(entry);   //!
    	b_nSelectedEvents->GetEntry(entry);   //!
    	b_sKey->GetEntry(entry);   //!  
    	b_sCrossSection->GetEntry(entry);   //!
    	b_sGMSBGravMass->GetEntry(entry);   //!
    	b_sGMSBChi1Mass->GetEntry(entry);   //!
    	b_sMCWgt->GetEntry(entry);   //!
    	b_sMCType->GetEntry(entry);   //!
        std::string configKey(*sKey);

		if( not configInfo.count(configKey) ){
        	if(debug) std::cout << " - Filling configValues. " << std::endl;
			std::map< std::string, float > configValues;
        	configValues["nEvents"] = nEvents;
        	configValues["nSelectedEvents"] = nSelectedEvents;
        	configValues["sCrossSection"] = sCrossSection;
        	configValues["sGMSBGravMass"] = sGMSBGravMass;
        	configValues["sGMSBChi1Mass"] = sGMSBChi1Mass;
        	configValues["sMCWgt"] = sMCWgt;
        	configValues["sMCType"] = sMCType;
        	if(debug) std::cout << " - Filling configInfo. " << std::endl;
        	configInfo[configKey] = configValues;
		} else {
			auto & configValues = configInfo[configKey];
			configValues["nEvents"] += nEvents;
			configValues["nSelectedEvents"] += nSelectedEvents;
		}//<<>>if( not configInfo.count(configKey) )

	}//<<>>for (Long64_t centry = 0; centry < nConfigEntries; centry++)

    for( auto item : configInfo ){ 
		std::cout << item.first << " ( ";  
		for( auto line : item.second ){ 
			std::cout << line.first <<  " " << line.second << " ";
		}//<<>>for( auto line : item.second )
		std::cout << ")" << std::endl;
	}//<<>>for( auto item : configInfo )

    std::cout << "<<<<<<<< Processing Event Loop <<<<<<<<<<<<<< " << std::endl;
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

	bool isgnsusy = ( cutselection == 1 ) ? true : false;
    bool iswzsusy = ( cutselection == 2 ) ? true : false;
    bool isOthersusy = ( cutselection == 3 ) ? true : false;
    bool isOthermatch = ( cutselection == 4 ) ? true : false;
    bool isXfsrsusy = ( cutselection == 5 ) ? true : false;
    bool isfsrsusy = ( cutselection == 6 ) ? true : false;

    bool isSig = ( cutselection == 7 ) ? true : false;
    bool isXfsr = ( cutselection == 8 ) ? true : false;
    bool isOther = ( cutselection == 9 ) ? true : false;
    bool isUnMatched = ( cutselection == 10 ) ? true : false;

    auto dskey  = *DataSetKey;
    //Float_t         selCMet;
    //Float_t         selCMetPx;
    //Float_t         selCMetPy;
    auto metPt = std::sqrt(rad2(selCMetPx,selCMetPy));
    if( metPt < 150 ) return;
    hist1d[400]->Fill( metPt );

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

    float phoeff = nPhotons ? static_cast<float>(nSelPhotons)/static_cast<float>(nPhotons) : 1.1;
    if( nPhotons > 0 ) hist1d[100]->Fill(phoeff); //("phoEff", "phoEff", 1000, 0, 1000);
    hist1d[118]->Fill(nSelPhotons);
    bool usepho = true;
    int nSusPhoGN1(0), nSusPhoGN1d(0), nSusPhoGN2(0), nSusPhoGN2d(0), nSusPhoWN(0), nSusPhoZC(0), nPhoNoGen(0);
    if( DEBUG ) std::cout << " - Looping over " << nSelPhotons << " photons" << std::endl;
    for( int it = 0; it < nSelPhotons; it++ ){

		// determine pog id class
		// -----------------------------------------------------

    //UInt_t          leadSelPho;
    //UInt_t          nPhotons;
    //UInt_t          nSelPhotons;
    //std::vector<float>   *selPhoClstrRn;
    //std::vector<float>   *selPhoEnergy;
    //std::vector<float>   *selPhoEta;
    //std::vector<float>   *selPhoGeoEgnVal;
    //std::vector<float>   *selPhoGeoSMaj;
    //std::vector<float>   *selPhoGeoSMin;
    //std::vector<unsigned int> *selPhoNrh;
    //std::vector<float>   *selPhoPhi;
    //std::vector<float>   *selPhoPt;
    //std::vector<int>     *selPhoQuality;
    //std::vector<float>   *selPhoR9;
    //std::vector<float>   *selPhoSMaj;
    //std::vector<float>   *selPhoSMin;
    //std::vector<float>   *selPhoSieie;
    //std::vector<float>   *selPhoTime;
    //UInt_t          subLeadSelPho;
//   vector<float>   *selPhoGenDp;
//   vector<float>   *selPhoGenDr;

        if( DEBUG ) std::cout << " -- looping photons : Start getting dr dp info" << std::endl;
		auto phoClass = (*selPhoQuality)[it];
        //if( phoClass < 1 ) continue;
        auto phoSusId = (*selPhoSusyId)[it];
        auto phoGenMatDr = (*selPhoGenDr)[it];
        auto phoGenMatDp = (*selPhoGenDp)[it];
        auto phoTime = (*selPhoTime)[it];

/*
        bool isgoodmatch = phoGenMatDr < 0.12 && phoGenMatDp < 1.0 && phoGenMatDr >= 0 && phoGenMatDp >= 0;
        bool isnogen = phoSusId < 20;
        if( isOthermatch && ( ( not isnogen ) && isgoodmatch )){ usepho = false; continue; } // no match
		if( ( not isOthermatch ) && ( isnogen || ( not isgoodmatch ) )){ usepho = false; continue; } // good match
*/

        bool isNino = phoSusId == 22;
        bool isXinoWZ = phoSusId == 23 || phoSusId == 24 || phoSusId == 33 || phoSusId == 34;
        bool isXfsr = phoSusId == 36 || phoSusId == 32 || phoSusId == 35;
        bool isPrompt = phoSusId == 29;
        bool isSfsr = phoSusId == 37 || phoSusId == 20 || phoSusId == 21 || phoSusId == 30 || phoSusId == 31;

		bool isSigType = isNino;
		bool isXfsrType = isXinoWZ || isXfsr;
		bool isOtherType = isPrompt || isSfsr; // || phoSusId < 1;

        if( isgnsusy && not isNino ){ usepho = false; continue; }
        if( iswzsusy && not isXinoWZ ){ usepho = false; continue; }
        if( isXfsrsusy && not isXfsr ){ usepho = false; continue; }
        if( isSfsrsusy && not isSfsr ){ usepho = false; continue; }
        if( isOthersusy && not isPrompt ){ usepho = false; continue; }

        if( isSig && not isSigType ){ usepho = false; continue; }
        if( isXfsr && not isXfsrType ){ usepho = false; continue; }
        if( isOther && not isOtherType ){ usepho = false; continue; }
        if( isUnMatched && phoSusId > 0 ){ usepho = false; continue; }


        hist2d[200]->Fill( phoGenMatDr, phoTime );
        hist2d[201]->Fill( phoGenMatDp, phoTime );
        hist2d[202]->Fill( phoGenMatDr, phoGenMatDp );
        hist1d[136]->Fill( phoGenMatDp );
        hist1d[137]->Fill( phoGenMatDr );

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;        

        if( isNino ) nSusPhoGN1++;
        if( isXfsr ) nSusPhoGN1d++;
        if( isPrompt ) nSusPhoGN2++;
        if( isSfsr ) nSusPhoGN2d++;        
        if( isXinoWZ ) nSusPhoWN++;
        if( false ) nSusPhoZC++;
        //if( isnogen ) nPhoNoGen++;
        hist1d[135]->Fill( phoSusId );

        hist2d[205]->Fill( (*selPhoEnergy)[it], (*selPhoPt)[it] );

        if( DEBUG ) std::cout << " -- looping photons : filling ids 0" << std::endl;

        hist1d[101]->Fill((*selPhoClstrRn)[it]); //("phoClstrRn", "phoClstrRn", 1000, 0, 1000);
        hist1d[102]->Fill((*selPhoEta)[it]); //("phoEta", "phoEta", 700, -3.5, 3.5);
        hist1d[103]->Fill((*selPhoPhi)[it]); //("phoPhi", "phoPhi", 700, -3.5, 3.5);
        hist1d[104]->Fill((*selPhoPt)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[105]->Fill((*selPhoEnergy)[it]); //("phoEnergy", "phoEnergy", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling ids 1" << std::endl;
        //hist1d[107]->Fill((*selPhoGeoEgnVal)[it]); //("phoGeoEgnVal", "phoGeoEgnVal", 100, 0, 1);
        //hist1d[108]->Fill((*selPhoGeoSMaj)[it]); //("phoGeoSMaj", "phoGeoSMaj", 100, 0, 1);
        hist1d[109]->Fill((*selPhoR9)[it]); //("phoR9", "phoR9", 100, 0, 1);
        hist1d[110]->Fill((*selPhoNrh)[it]); //("phoNrh", "phoNrh", 100, 0, 100);
        hist1d[111]->Fill((*selPhoSMaj)[it]); //("phoSMaj", "phoSMaj", 100, 0, 1);
        hist1d[112]->Fill((*selPhoSMin)[it]); //("phoSMin", "phoSMin", 100, 0, 1);
        hist1d[113]->Fill((*selPhoSieie)[it]); //("phoSieie", "phoSieie", 100, 0, 1);
        hist1d[114]->Fill((*selPhoTime)[it]); //("phoTime", "phoTime", 500, -25, 25);
        hist1d[115]->Fill((*selPhoQuality)[it]); //("phoQuality", "phoQuality", 4, 0, 4);
        hist1d[116]->Fill((*selPhoR9)[it]); //("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
        if( DEBUG ) std::cout << " -- looping photons : filling ids 2" << std::endl;
        //hist1d[117]->Fill((*selPhoGeoSMin)[it]);
        hist1d[125]->Fill((*selPhoHcalTowerSumEtBcConeDR04)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[126]->Fill((*selPhoTrkSumPtHollowConeDR03)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[127]->Fill((*selPhoTrkSumPtHollowConeDR04)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[128]->Fill((*selPhoTrkSumPtSolidConeDR04)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        //hist1d[129]->Fill((*selPhoPixelSeed)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling var hist 5" << std::endl;
        hist1d[130]->Fill((*selPhoEcalRHSumEtConeDR04)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[131]->Fill((*selPhoHadTowOverEM)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling var hist 6" << std::endl;
        hist1d[132]->Fill((*selPhoSieip)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        //if( DEBUG ) std::cout << " -- looping photons : selPhoSipip = " << (*selPhoSipip)[it] << std::endl;
        //hist1d[134]->Fill((*selPhoSipip)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling var hist done" << std::endl;

        if( DEBUG ) std::cout << " -- looping photons : gen counts" << std::endl; 
    	hist1d[119]->Fill(nSusPhoGN1);
    	hist1d[120]->Fill(nSusPhoGN1d);
    	hist1d[121]->Fill(nSusPhoGN2);
    	hist1d[122]->Fill(nSusPhoGN2d);
    	hist1d[123]->Fill(nSusPhoWN);
    	hist1d[124]->Fill(nSusPhoZC);
    	hist1d[125]->Fill(nPhoNoGen);

/*
		//find intial cut values 
		//----------------------------------------------------------------
		if( DEBUG ) std::cout << " -- pho rechits" << std::endl;
		auto nrh = (*selPhoNrh)[it];
        auto phoClstrR9 = (*selPhoClstrRn)[it]; //clstrR9( (*phoRhIds)[it] );

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


        if( DEBUG ) std::cout << " -- sel photons # " << it << std::endl;
		if( usePho && goodRhCnt && minPhoPt ) { //----------------------------------- 
	
	        //hist1d[132]->Fill(1);
	
			if( DEBUG ) std::cout << " -- pho time" << std::endl;


        }//<<>>if( usePho ) { //----------------------------------- 
*/

    }//<<>>for( int it = 0; it < nPhotons; it++ )
 
    if( DEBUG ) std::cout << " -- looping photons : filling lead phio" << std::endl;
    if( leadSelPho >= 0 && usepho ){

        float phoeffl = nPhotons ? static_cast<float>(nSelPhotons)/static_cast<float>(nPhotons) : 0;
        hist1d[150]->Fill(phoeffl); //("phoEff", "phoEff", 1000, 0, 1000);
        hist1d[151]->Fill((*selPhoClstrRn)[leadSelPho]); //("phoClstrRn", "phoClstrRn", 1000, 0, 1000);
        hist1d[152]->Fill((*selPhoEta)[leadSelPho]); //("phoEta", "phoEta", 700, -3.5, 3.5);
        hist1d[153]->Fill((*selPhoPhi)[leadSelPho]); //("phoPhi", "phoPhi", 700, -3.5, 3.5);
        hist1d[154]->Fill((*selPhoPt)[leadSelPho]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[155]->Fill((*selPhoEnergy)[leadSelPho]); //("phoEnergy", "phoEnergy", 500, 0, 5000);
        //hist1d[157]->Fill((*selPhoGeoEgnVal)[leadSelPho]); //("phoGeoEgnVal", "phoGeoEgnVal", 100, 0, 1);
        //hist1d[158]->Fill((*selPhoGeoSMaj)[leadSelPho]); //("phoGeoSMaj", "phoGeoSMaj", 100, 0, 1);
        hist1d[159]->Fill((*selPhoR9)[leadSelPho]); //("phoR9", "phoR9", 100, 0, 1);
        hist1d[160]->Fill((*selPhoNrh)[leadSelPho]); //("phoNrh", "phoNrh", 100, 0, 100);
        hist1d[161]->Fill((*selPhoSMaj)[leadSelPho]); //("phoSMaj", "phoSMaj", 100, 0, 1);
        hist1d[162]->Fill((*selPhoSMin)[leadSelPho]); //("phoSMin", "phoSMin", 100, 0, 1);
        hist1d[163]->Fill((*selPhoSieie)[leadSelPho]); //("phoSieie", "phoSieie", 100, 0, 1);
        hist1d[164]->Fill((*selPhoTime)[leadSelPho]); //("phoTime", "phoTime", 500, -25, 25);
        hist1d[165]->Fill((*selPhoQuality)[leadSelPho]); //("phoQuality", "phoQuality", 4, 0, 4);
        hist1d[166]->Fill((*selPhoR9)[leadSelPho]); //("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
        //hist1d[167]->Fill((*selPhoGeoSMin)[leadSelPho]);

	}//<<>>if( leadSelPho >= 0 )

    if( DEBUG ) std::cout << " -- looping photons : filling sub lead phio" << std::endl;
    if( subLeadSelPho >= 0 && usepho ){

        float phoeffsl = nPhotons ? static_cast<float>(nSelPhotons)/static_cast<float>(nPhotons) : 0;
        if( nPhotons > 0 ) hist1d[200]->Fill(phoeffsl); //("phoEff", "phoEff", 1000, 0, 1000);
        hist1d[201]->Fill((*selPhoClstrRn)[subLeadSelPho]); //("phoClstrRn", "phoClstrRn", 1000, 0, 1000);
        hist1d[202]->Fill((*selPhoEta)[subLeadSelPho]); //("phoEta", "phoEta", 700, -3.5, 3.5);
        hist1d[203]->Fill((*selPhoPhi)[subLeadSelPho]); //("phoPhi", "phoPhi", 700, -3.5, 3.5);
        hist1d[204]->Fill((*selPhoPt)[subLeadSelPho]); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[205]->Fill((*selPhoEnergy)[subLeadSelPho]); //("phoEnergy", "phoEnergy", 500, 0, 5000);
        //hist1d[207]->Fill((*selPhoGeoEgnVal)[subLeadSelPho]); //("phoGeoEgnVal", "phoGeoEgnVal", 100, 0, 1);
        //hist1d[208]->Fill((*selPhoGeoSMaj)[subLeadSelPho]); //("phoGeoSMaj", "phoGeoSMaj", 100, 0, 1);
        hist1d[209]->Fill((*selPhoR9)[subLeadSelPho]); //("phoR9", "phoR9", 100, 0, 1);
        hist1d[210]->Fill((*selPhoNrh)[subLeadSelPho]); //("phoNrh", "phoNrh", 100, 0, 100);
        hist1d[211]->Fill((*selPhoSMaj)[subLeadSelPho]); //("phoSMaj", "phoSMaj", 100, 0, 1);
        hist1d[212]->Fill((*selPhoSMin)[subLeadSelPho]); //("phoSMin", "phoSMin", 100, 0, 1);
        hist1d[213]->Fill((*selPhoSieie)[subLeadSelPho]); //("phoSieie", "phoSieie", 100, 0, 1);
        hist1d[214]->Fill((*selPhoTime)[subLeadSelPho]); //("phoTime", "phoTime", 500, -25, 25);
        hist1d[215]->Fill((*selPhoQuality)[subLeadSelPho]); //("phoQuality", "phoQuality", 4, 0, 4);
        hist1d[216]->Fill((*selPhoR9)[subLeadSelPho]); //("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
        //hist1d[217]->Fill((*selPhoGeoSMin)[subLeadSelPho]);

    }//<<>>if( subLeadSelPho >= 0 )

	//-------- electrons --------------------------------------

	//--------- jets --------------------------
    if( DEBUG ) std::cout << "Finding Jets with " << nSelJets << " selected. "<< std::endl;

    //UInt_t          nJets;
    //UInt_t          nSelJets;
    //std::vector<float>   *selJetEnergy;
    //std::vector<float>   *selJetEta;
    //std::vector<float>   *selJetMass;
    //std::vector<float>   *selJetPhi;
    //std::vector<float>   *selJetPt;
    //std::vector<int>     *selJetQuality;
    //std::vector<float>   *selJetTime;


	//// ---  base info  ----
    hist1d[9]->Fill(nSelJets); //("nJet", mkht(ht,"nJets"), 100, 0, 100);
    
    int nSusJetGino(0), nSusJetGinod(0), nSusJetSqrk(0), nSusJetSqrkd(0), nSusJetSusy(0);
    for( int it = 0; it < nSelJets; it++ ){

		if( (*selJetSusyId)[it] == 21 ) nSusJetGino++;
        if( (*selJetSusyId)[it] == 31 ) nSusJetGinod++;
        if( (*selJetSusyId)[it] == 20 ) nSusJetSqrk++;
        if( (*selJetSusyId)[it] == 30 ) nSusJetSqrkd++;
        if( (*selJetSusyId)[it] == 21 || (*selJetSusyId)[it] == 20 ) nSusJetSusy++;

/*
        bool isgoodmatch = jetGenMatDr < 0.3 && jetGenMatDp < 1.0 && jetGenMatDr >= 0 && jetGenMatDp >= 0;
        bool isnogen = jetSusId < 19;
        if( isOthermatch && ( ( not isnogen ) && isgoodmatch )){ continue; } // no match
        if( ( not isOthermatch ) && ( isnogen || ( not isgoodmatch ) )){ continue; } // good match
*/

        auto jetSusId = (*selJetSusyId)[it];
        auto jetGenMatDr = (*selJetLlpDr)[it];
        auto jetGenMatDp = (*selJetLlpDp)[it];

        bool issqg = jetSusId == 21 || jetSusId == 20 || jetSusId == 31 || jetSusId == 30;
        bool iswzino = jetSusId == 23 || jetSusId == 24 || jetSusId == 33 || jetSusId == 34;
        bool isxrad = jetSusId == 32 || jetSusId == 35 || jetSusId == 22 || jetSusId == 25;
        bool isproton = jetSusId == 29;
        bool issusfsr = jetSusId == 36 || jetSusId == 37;

        bool isSigt = issqg;
        bool isXfsrt = iswzino || issusfsr || isxrad;
        bool isOthert = isproton; // || jetSusId < 1;

        if( isgnsusy && not issqg ){ continue; }
        if( iswzsusy && not iswzino ){ continue; }
        if( isXfsrsusy && not isxrad ){ continue; }
        if( isfsrsusy && not issusfsr ){ continue; }
        if( isOthersusy && not isproton ){ continue; }

        if( isSig && not isSigt ){ continue; }
        if( isXfsr && not isXfsrt ){ continue; }
        if( isOther && not isOthert ){ continue; }
        if( isUnMatched && jetSusId > 0 ){ continue; }

        hist1d[0]->Fill((*selJetPt)[it]); //("jetPt", mkht(ht,"jetPt"), 500, 0, 5000);
        hist1d[1]->Fill((*selJetPhi)[it]); //("jetPhi", mkht(ht,"jetPhi"), 700, -3.5, 3.5);
        hist1d[2]->Fill((*selJetEta)[it]); //("jetEta", mkht(ht,"jetEta"), 700, -3.5, 3.5);
        hist1d[3]->Fill((*selJetEnergy)[it]); //("jetEnergy", mkht(ht,"jetEnergy"), 500, 0, 5000);
        hist1d[4]->Fill((*selJetMass)[it]); //("jetMass", mkht(ht,"jetMass"), 500, 0, 5000);
        hist1d[5]->Fill((*selJetQuality)[it]); //("jetQuality", mkht(ht,"jetQuality"), 5, 0, 5);
        hist1d[6]->Fill((*selJetTime)[it]); //("jetTime", mkht(ht,"jetTime"), 500, -25, 25);

        hist1d[30]->Fill((*selGenJetDpt)[it]); //("genJetDpt", mkht(ht,"genJetDpt").c_str(), 500, 0, 5);
        hist1d[31]->Fill((*selGenJetEnergy)[it]); //("genJetEnergy", mkht(ht,"genJetEnergy").c_str(), 5000, 0, 5000);
        hist1d[32]->Fill((*selGenJetImpAng)[it]); //("genJetImpAng", mkht(ht,"genJetImpAng").c_str(), 200, -10, 10);
        hist1d[33]->Fill((*selGenJetLlpTime)[it]); //("genJetLlpTime", mkht(ht,"genJetLlpTime").c_str(), 100, -5, 5);
        hist1d[34]->Fill((*selGenJetPt)[it]); //("genJetPt", mkht(ht,"genJetPt").c_str(), 3000, 0, 3000);
        hist1d[35]->Fill((*selGenJetTime)[it]); //("genJetTime", mkht(ht,"genJetTime").c_str(), 250, 0, 25);
        hist1d[36]->Fill((*selGenJetTof)[it]); //("genJetTof", mkht(ht,"genJetTof").c_str(), 400, 0, 400);
        hist1d[37]->Fill((*selGenJetdr)[it]); //("genJetdr", mkht(ht,"genJetdr").c_str(), 600, -1, 5);
        hist1d[38]->Fill((*selGenJeteta)[it]); //("genJeteta", mkht(ht,"genJeteta").c_str(), 160, -10, 6);

        hist1d[39]->Fill((*selJetSusyId)[it]); //("jetSusyId", mkht(ht,"jetSusyId").c_str(), 45, -5, 40);
        hist1d[40]->Fill((*selJetLlpDp)[it]); //("jetLlpDp", mkht(ht,"jetLlpDp").c_str(), 500, 0, 5);
        hist1d[41]->Fill((*selJetLlpDr)[it]); //("jetLlpDr", mkht(ht,"jetLlpDr").c_str(), 500, 0, 5);
        hist2d[203]->Fill((*selJetLlpDr)[it],(*selJetLlpDp)[it]);

        hist1d[42]->Fill((*selJetArea)[it]); //("jetArea", mkht(ht,"jetArea").c_str(), 100, 0, 1);
        hist1d[43]->Fill((*selJetChEmEF)[it]); //("jetChEmEF", mkht(ht,"jetChEmEF").c_str(), 100, 0, 1);
        hist1d[44]->Fill((*selJetChHM)[it]); //("jetChHM", mkht(ht,"jetChHM").c_str(), 100, 0, 100);
        hist1d[45]->Fill((*selJetMuEF)[it]); //("jetMuEF", mkht(ht,"jetMuEF").c_str(), 100, 0, 1);
        hist1d[46]->Fill((*selJetNeEmEF)[it]); //("jetNeEmEF", mkht(ht,"jetNeEmEF").c_str(), 100, 0, 1);
        hist1d[47]->Fill((*selJetNeHEF)[it]); //("jetNeHEF", mkht(ht,"jetNeHEF").c_str(), 100, 0, 1);
        hist1d[48]->Fill((*selJetNeHM)[it]); //("jetNeHM", mkht(ht,"jetNeHM").c_str(), 75, 0, 75);
        hist1d[49]->Fill((*selJetchHEF)[it]); //("jetchHEF", mkht(ht,"jetchHEF").c_str(), 100, 0, 1);


	}//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Filling jet stats" << std::endl;

    float jeteff = nSelJets ? static_cast<float>(nSelJets)/static_cast<float>(nJets) : 1.1;
    hist1d[7]->Fill(jeteff);
    hist1d[10]->Fill(nSusJetGino);
    hist1d[11]->Fill(nSusJetGinod);
    hist1d[12]->Fill(nSusJetSqrk);
    hist1d[13]->Fill(nSusJetSqrkd);
    hist1d[14]->Fill(nSusJetSusy);

    if( DEBUG ) std::cout << "Finding genjet" << std::endl;
	//// --- genjet info -----------------------------------

    if( DEBUG ) std::cout << "Finding rjr" << std::endl;

    //std::vector<float>   *SCosA;
    //std::vector<float>   *SMass;
    //std::vector<float>   *X1aCosA;
    //std::vector<float>   *X1aMass;
    //std::vector<float>   *X1bCosA;
    //std::vector<float>   *X1bMass;
    //std::vector<float>   *X2aCosA;
    //std::vector<float>   *X2aMass;
    //std::vector<float>   *X2bCosA;
    //std::vector<float>   *X2bMass;

    if( SCosA->size() > 0 ){

        hist1d[250]->Fill((*SCosA)[0]); //("SCosA", mkht(ht,"SCosA"), 70, -3.5, 3.5);
        hist1d[251]->Fill((*SMass)[0]); //("SMass", mkht(ht,"SMass"), 500, 0, 5000);
        hist1d[252]->Fill((*X1aCosA)[0]); //("X1aCosA", mkht(ht,"X1aCosA"), 70, -3.5, 3.5);
        hist1d[253]->Fill((*X1aMass)[0]); //("X1aMass", mkht(ht,"X1aMass"), 500, 0, 5000);
        hist1d[254]->Fill((*X1bCosA)[0]); //("X1bCosA", mkht(ht,"X1bCosA"), 70, -3.5, 3.5);
        hist1d[255]->Fill((*X1bMass)[0]); //("X1bMass", mkht(ht,"X1bMass"), 500, 0, 5000);
        hist1d[256]->Fill((*X2aCosA)[0]); //("X2aCosA", mkht(ht,"X2aCosA"), 70, -3.5, 3.5);
        hist1d[257]->Fill((*X2aMass)[0]); //("X2aMass", mkht(ht,"X2aMass"), 500, 0, 5000);
        hist1d[258]->Fill((*X2bCosA)[0]); //("X2bCosA", mkht(ht,"X2bCosA"), 70, -3.5, 3.5);
        hist1d[259]->Fill((*X2bMass)[0]); //("X2bMass", mkht(ht,"X2bMass"), 500, 0, 5000);

    }//<<>>if( SCosA->size() > 0 )

}//<<>>void HistMaker::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------


	//------- jets 0 - 99
    //UInt_t          nJets;
    //UInt_t          nSelJets;
    //std::vector<float>   *selJetEnergy;
    //std::vector<float>   *selJetEta;
    //std::vector<float>   *selJetMass;
    //std::vector<float>   *selJetPhi;
    //std::vector<float>   *selJetPt;
    //std::vector<int>     *selJetQuality;
    //std::vector<float>   *selJetTime;

    std::cout << " title test : " << mkht(ht,"jetPt") << std::endl;

    hist1d[0] = new TH1D("jetPt", mkht(ht,"jetPt").c_str(), 500, 0, 5000);
    hist1d[1] = new TH1D("jetPhi", mkht(ht,"jetPhi").c_str(), 700, -3.5, 3.5);
    hist1d[2] = new TH1D("jetEta", mkht(ht,"jetEta").c_str(), 700, -3.5, 3.5);
    hist1d[3] = new TH1D("jetEnergy", mkht(ht,"jetEnergy").c_str(), 75, 0, 750);
    hist1d[4] = new TH1D("jetMass", mkht(ht,"jetMass").c_str(), 500, 0, 5000);
    hist1d[5] = new TH1D("jetQuality", mkht(ht,"jetQuality").c_str(), 5, 0, 5);
    hist1d[6] = new TH1D("jetTime", mkht(ht,"jetTime").c_str(), 500, -25, 25);
    hist1d[7] = new TH1D("jetEff", mkht(ht,"jetEff").c_str(), 11, 0, 1.1);
    hist1d[9] = new TH1D("nJet", mkht(ht,"nJets").c_str(), 25, 0, 25);

    hist1d[10] = new TH1D("nGinoJets", mkht(ht,"nGinoJets").c_str(), 10, 0, 10);
    hist1d[11] = new TH1D("nGinodJets", mkht(ht,"nGinodJets").c_str(), 10, 0, 10);
    hist1d[12] = new TH1D("nSqrkJets", mkht(ht,"nSqrkJets").c_str(), 10, 0, 10);
    hist1d[13] = new TH1D("nSqrkdJets", mkht(ht,"nSqrkdJets").c_str(), 10, 0, 10);
    hist1d[14] = new TH1D("nSusyJets", mkht(ht,"nSusyJets").c_str(), 10, 0, 10);

    hist1d[30] = new TH1D("genJetRe", mkht(ht,"genJetRe").c_str(), 500, 0, 5);
    hist1d[31] = new TH1D("genJetEnergy", mkht(ht,"genJetEnergy").c_str(), 5000, 0, 5000);
    hist1d[32] = new TH1D("genJetImpAng", mkht(ht,"genJetImpAng").c_str(), 200, -10, 10);
    hist1d[33] = new TH1D("genJetLlpTime", mkht(ht,"genJetLlpTime").c_str(), 100, -5, 5);
    hist1d[34] = new TH1D("genJetPt", mkht(ht,"genJetPt").c_str(), 3000, 0, 3000);
    hist1d[35] = new TH1D("genJetTime", mkht(ht,"genJetTime").c_str(), 250, 0, 25);
    hist1d[36] = new TH1D("genJetTof", mkht(ht,"genJetTof").c_str(), 400, 0, 400);
    hist1d[37] = new TH1D("genJetDr", mkht(ht,"genJetDr").c_str(), 600, -1, 5);
    hist1d[38] = new TH1D("genJeteta", mkht(ht,"genJeteta").c_str(), 160, -10, 6);

    hist1d[39] = new TH1D("jetSusyId", mkht(ht,"jetSusyId").c_str(), 45, -5, 40);
    hist1d[40] = new TH1D("jetLlpRe", mkht(ht,"jetLlpRe").c_str(), 500, 0, 5);
    hist1d[41] = new TH1D("jetLlpDr", mkht(ht,"jetLlpDr").c_str(), 500, 0, 5);

    hist1d[42] = new TH1D("jetArea", mkht(ht,"jetArea").c_str(), 100, 0, 1);
    hist1d[43] = new TH1D("jetChEmEF", mkht(ht,"jetChEmEF").c_str(), 100, 0, 1);
    hist1d[44] = new TH1D("jetChHM", mkht(ht,"jetChHM").c_str(), 100, 0, 100);
    hist1d[45] = new TH1D("jetMuEF", mkht(ht,"jetMuEF").c_str(), 100, 0, 1);
    hist1d[46] = new TH1D("jetNeEmEF", mkht(ht,"jetNeEmEF").c_str(), 100, 0, 1);
    hist1d[47] = new TH1D("jetNeHEF", mkht(ht,"jetNeHEF").c_str(), 100, 0, 1);
    hist1d[48] = new TH1D("jetNeHM", mkht(ht,"jetNeHM").c_str(), 75, 0, 75);
    hist1d[49] = new TH1D("jetChHEF", mkht(ht,"jetChHEF").c_str(), 100, 0, 1);


	//----- photons 100 - 249
    //UInt_t          leadSelPho;
    //UInt_t          nPhotons;
    //UInt_t          nSelPhotons;
    //std::vector<float>   *selPhoClstrRn;
    //std::vector<float>   *selPhoEnergy;
    //std::vector<float>   *selPhoEta;
    //std::vector<float>   *selPhoGeoEgnVal;
    //std::vector<float>   *selPhoGeoSMaj;
    //std::vector<float>   *selPhoGeoSMin;
    //std::vector<unsigned int> *selPhoNrh;
    //std::vector<float>   *selPhoPhi;
    //std::vector<float>   *selPhoPt;
    //std::vector<int>     *selPhoQuality;
    //std::vector<float>   *selPhoR9;
    //std::vector<float>   *selPhoSMaj;
    //std::vector<float>   *selPhoSMin;
    //std::vector<float>   *selPhoSieie;
    //std::vector<float>   *selPhoTime;
    //UInt_t          subLeadSelPho;

    hist1d[100] = new TH1D("phoEff", mkht(ht,"phoEff").c_str(), 11, 0, 1.1);
    hist1d[101] = new TH1D("phoClstrRn", mkht(ht,"phoClstrRn").c_str(), 100, 0, 1);
    hist1d[102] = new TH1D("phoEta", mkht(ht,"phoEta").c_str(), 700, -3.5, 3.5);
    hist1d[103] = new TH1D("phoPhi", mkht(ht,"phoPhi").c_str(), 700, -3.5, 3.5);
    hist1d[104] = new TH1D("phoPt", mkht(ht,"phoPt").c_str(), 750, 0, 750);
    hist1d[105] = new TH1D("phoEnergy", mkht(ht,"phoEnergy").c_str(), 500, 0, 5000);
    //hist1d[107] = new TH1D("phoGeoEgnVal", mkht(ht,"phoGeoEgnVal").c_str(), 100, 0, 1);
    //hist1d[108] = new TH1D("phoGeoSMaj", mkht(ht,"phoGeoSMaj").c_str(), 100, 0, 1);
    hist1d[109] = new TH1D("phoR9", mkht(ht,"phoR9").c_str(), 100, 0, 1);
    hist1d[110] = new TH1D("phoNrh", mkht(ht,"phoNrh").c_str(), 100, 0, 100); 
    hist1d[111] = new TH1D("phoSMaj", mkht(ht,"phoSMaj").c_str(), 20, 0, 2); 
    hist1d[112] = new TH1D("phoSMin", mkht(ht,"phoSMin").c_str(), 100, 0, 1); 
    hist1d[113] = new TH1D("phoSieie", mkht(ht,"phoSieie").c_str(), 100, 0, 0.02); 
    hist1d[114] = new TH1D("phoTime", mkht(ht,"phoTime").c_str(), 500, -25, 25);
    hist1d[115] = new TH1D("phoQuality", mkht(ht,"phoQuality").c_str(), 4, 0, 4);
    hist1d[116] = new TH1D("phoR9zoom", mkht(ht,"phoR9zoom").c_str(), 200, 0.8, 1);
    //hist1d[117] = new TH1D("phoGeoSMin", mkht(ht,"phoGeoSMin").c_str(), 100, 0, 1);
    hist1d[118] = new TH1D("nPhotons", mkht(ht,"nPhotons").c_str(), 20, 0, 20);

    hist1d[119] = new TH1D("nGN1Photons", mkht(ht,"nGN1Photons").c_str(), 20, 0, 20);
    hist1d[120] = new TH1D("nGN1dPhotons", mkht(ht,"nGN1dPhotons").c_str(), 20, 0, 20);
    hist1d[121] = new TH1D("nGN2Photons", mkht(ht,"nGN2Photons").c_str(), 20, 0, 20);
    hist1d[122] = new TH1D("nGN2dPhotons", mkht(ht,"nGN2dPhotons").c_str(), 20, 0, 20);
    hist1d[123] = new TH1D("nWNPhotons", mkht(ht,"nWNPhotons").c_str(), 20, 0, 20);
    hist1d[124] = new TH1D("nZCPhotons", mkht(ht,"nZCPhotons").c_str(), 20, 0, 20);

    hist1d[125] = new TH1D("selPhoHcalTowerSumEtBcConeDR04", mkht(ht,"selPhoHcalTowerSumEtBcConeDR04").c_str(), 200, 0, 20);
    hist1d[126] = new TH1D("selPhoTrkSumPtHollowConeDR03", mkht(ht,"selPhoTrkSumPtHollowConeDR03").c_str(), 200, 0, 20);
    hist1d[127] = new TH1D("selPhoTrkSumPtHollowConeDR04", mkht(ht,"selPhoTrkSumPtHollowConeDR04").c_str(), 200, 0, 20);
    hist1d[128] = new TH1D("selPhoTrkSumPtSolidConeDR04", mkht(ht,"selPhoTrkSumPtSolidConeDR04").c_str(), 200, 0, 20);
    hist1d[129] = new TH1D("selPhoPixelSeed", mkht(ht,"selPhoPixelSeed").c_str(), 2, 0, 1);
    hist1d[130] = new TH1D("selPhoEcalRHSumEtConeDR04", mkht(ht,"selPhoEcalRHSumEtConeDR04").c_str(), 200, 0, 20);
    hist1d[131] = new TH1D("selPhoHadTowOverEM", mkht(ht,"selPhoHadTowOverEM").c_str(), 100, 0, 0.1);
    hist1d[132] = new TH1D("selPhoSieip", mkht(ht,"selPhoSieip").c_str(), 100, 0, 0.01);
    //hist1d[134] = new TH1D("selPhoSipip", mkht(ht,"selPhoSipip").c_str(), 100, 0, 0.01);

    hist1d[134] = new TH1D("nNoGenMatch", mkht(ht,"nNoGenMatch").c_str(), 20, 0, 20);
    hist1d[135] = new TH1D("phoGenType", mkht(ht,"phoGenType").c_str(), 45, -5, 40);
    hist1d[136] = new TH1D("phoGenRe", mkht(ht,"phoGenRe").c_str(), 500, 0, 5);
    hist1d[137] = new TH1D("phoGenDr", mkht(ht,"phoGenDr").c_str(), 50, 0, 0.5);

    hist1d[150] = new TH1D("phoLdEff", mkht(ht,"phoLdEff").c_str(), 10, 0, 1);
    hist1d[151] = new TH1D("phoLdClstrRn", mkht(ht,"phoLdClstrRn").c_str(), 100, 0, 1);
    hist1d[152] = new TH1D("phoLdEta", mkht(ht,"phoLdEta").c_str(), 700, -3.5, 3.5);
    hist1d[153] = new TH1D("phoLdPhi", mkht(ht,"phoLdPhi").c_str(), 700, -3.5, 3.5);
    hist1d[154] = new TH1D("phoLdPt", mkht(ht,"phoLdPt").c_str(), 200, 0, 2000);
    hist1d[155] = new TH1D("phoLdEnergy", mkht(ht,"phoLdEnergy").c_str(), 500, 0, 5000);
    //hist1d[157] = new TH1D("phoLdGeoEgnVal", mkht(ht,"phoLdGeoEgnVal").c_str(), 100, 0, 1);
    //hist1d[158] = new TH1D("phoLdGeoSMaj", mkht(ht,"phoLdGeoSMaj").c_str(), 100, 0, 1);
    hist1d[159] = new TH1D("phoLdR9", mkht(ht,"phoLdR9").c_str(), 100, 0, 1);
    hist1d[160] = new TH1D("phoLdNrh", mkht(ht,"phoLdNrh").c_str(), 100, 0, 100);
    hist1d[161] = new TH1D("phoLdSMaj", mkht(ht,"phoLdSMaj").c_str(), 20, 0, 2);
    hist1d[162] = new TH1D("phoLdSMin", mkht(ht,"phoLdSMin").c_str(), 100, 0, 1);
    hist1d[163] = new TH1D("phoLdSieie", mkht(ht,"phoLdSieie").c_str(), 100, 0, 0.01);
    hist1d[164] = new TH1D("phoLdTime", mkht(ht,"phoLdTime").c_str(), 500, -25, 25);
    hist1d[165] = new TH1D("phoLdQuality", mkht(ht,"phoLdQuality").c_str(), 4, 0, 4);
    hist1d[166] = new TH1D("phoLdR9_zoom", mkht(ht,"phoLdR9_zoom").c_str(), 200, 0.8, 1);
    //hist1d[167] = new TH1D("phoLdGeoSMin", mkht(ht,"phoLdGeoSMin").c_str(), 100, 0, 1);

    hist1d[200] = new TH1D("phoSLdEff", mkht(ht,"phoSLdEff").c_str(), 10, 0, 1);
    hist1d[201] = new TH1D("phoSLdClstrRn", mkht(ht,"phoSLdClstrRn").c_str(), 100, 0, 1);
    hist1d[202] = new TH1D("phoSLdEta", mkht(ht,"phoSLdEta").c_str(), 700, -3.5, 3.5);
    hist1d[203] = new TH1D("phoSLdPhi", mkht(ht,"phoSLdPhi").c_str(), 700, -3.5, 3.5);
    hist1d[204] = new TH1D("phoSLdPt", mkht(ht,"phoSLdPt").c_str(), 200, 0, 2000);
    hist1d[205] = new TH1D("phoSLdEnergy", mkht(ht,"phoSLdEnergy").c_str(), 500, 0, 5000);
    //hist1d[207] = new TH1D("phoSLdGeoEgnVal", mkht(ht,"phoSLdGeoEgnVal").c_str(), 100, 0, 1);
    //hist1d[208] = new TH1D("phoSLdGeoSMaj", mkht(ht,"phoSLdGeoSMaj").c_str(), 100, 0, 1);
    hist1d[209] = new TH1D("phoSLdR9", mkht(ht,"phoSLdR9").c_str(), 100, 0, 1);
    hist1d[210] = new TH1D("phoSLdNrh", mkht(ht,"phoSLdNrh").c_str(), 100, 0, 100);
    hist1d[211] = new TH1D("phoSLdSMaj", mkht(ht,"phoSLdSMaj").c_str(), 20, 0, 2);
    hist1d[212] = new TH1D("phoSLdSMin", mkht(ht,"phoSLdSMin").c_str(), 100, 0, 1);
    hist1d[213] = new TH1D("phoSLdSieie", mkht(ht,"phoSLdSieie").c_str(), 100, 0, 0.01);
    hist1d[214] = new TH1D("phoSLdTime", mkht(ht,"phoSLdTime").c_str(), 500, -25, 25);
    hist1d[215] = new TH1D("phoSLdQuality", mkht(ht,"phoSLdQuality").c_str(), 4, 0, 4);
    hist1d[216] = new TH1D("phoSLdR9_zoom", mkht(ht,"phoSLdR9_zoom").c_str(), 200, 0.8, 1);
    //hist1d[217] = new TH1D("phoSLdGeoSMin", mkht(ht,"phoSLdGeoSMin").c_str(), 100, 0, 1);

    //hist1d[132] = new TH1D("phoType", mkht(ht,"phoType(1-pho 2-susy 3-oot 4-cmb)", 6,0,5);

	//------  genparticles 200 - 249

	//------  rjr 250 - 299
    //std::vector<float>   *SCosA;
    //std::vector<float>   *SMass;
    //std::vector<float>   *X1aCosA;
    //std::vector<float>   *X1aMass;
    //std::vector<float>   *X1bCosA;
    //std::vector<float>   *X1bMass;
    //std::vector<float>   *X2aCosA;
    //std::vector<float>   *X2aMass;
    //std::vector<float>   *X2bCosA;
    //std::vector<float>   *X2bMass;

    hist1d[250] = new TH1D("SCosA", mkht(ht,"SCosA").c_str(), 70, -3.5, 3.5);
    hist1d[251] = new TH1D("SMass", mkht(ht,"SMass").c_str(), 500, 0, 5000);
    hist1d[252] = new TH1D("X1aCosA", mkht(ht,"X1aCosA").c_str(), 70, -3.5, 3.5);
    hist1d[253] = new TH1D("X1aMass", mkht(ht,"X1aMass").c_str(), 500, 0, 5000);
    hist1d[254] = new TH1D("X1bCosA", mkht(ht,"X1bCosA").c_str(), 70, -3.5, 3.5);   
    hist1d[255] = new TH1D("X1bMass", mkht(ht,"X1bMass").c_str(), 500, 0, 5000);
    hist1d[256] = new TH1D("X2aCosA", mkht(ht,"X2aCosA").c_str(), 70, -3.5, 3.5);
    hist1d[257] = new TH1D("X2aMass", mkht(ht,"X2aMass").c_str(), 500, 0, 5000);
    hist1d[258] = new TH1D("X2bCosA", mkht(ht,"X2bCosA").c_str(), 70, -3.5, 3.5);
    hist1d[259] = new TH1D("X2bMass", mkht(ht,"X2bMass").c_str(), 500, 0, 5000);

    //------ ecal rechits 300 - 349

    //------  electrons 350 - 400

	//-------- event vars 400 - 450
    hist1d[400] = new TH1D("evtMetPt", mkht(ht,"evtMetPt").c_str(), 500, 0, 5000 );
	
    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	//------- jets ( time ) 0-49 ------------------------------

	//---jet id stuff 50 - 99 ---------------------------------------------------

	//--- Photons 200 - 349 -------------------------------------------
	
    hist2d[200] = new TH2D("phoGenDrVTime", mkht(ht,"pho genDr v phoTime; dR; phoTime [ns]").c_str(), 30, 0, 0.3, 350, -10, 25);
    hist2d[201] = new TH2D("phoGenDeVTime", mkht(ht,"pho genDe v phoTime; dE; phoTime [ns]").c_str(), 600,0, 3, 350, -10, 25);
    hist2d[202] = new TH2D("phoGenDrVphoGenDe", mkht(ht,"pho genDr v genDe; dR; dE").c_str(), 30, 0, 0.3, 600, 0, 3);
    // ....
    hist2d[203] = new TH2D("jetGenDrVjetGenDe", mkht(ht,"jet genDr v genDe; dR; dE").c_str(), 30, 0, 0.3, 600, 0, 3);

    // ....
    hist2d[205] = new TH2D("phoEvPt", mkht(ht,"pho E v Pt;Energy [GeV];Pt [GeV]").c_str(), 2000, 0, 2000, 1000, 0, 1000);



	//60 - 63

	//--- rechit collections 350 - 399 -------------------------------------------------

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
                const std::string listdir = "skims_files/";

                auto infilename = "KUCMS_GMSB_Skim_List.txt";

				int cut = 9;

                //auto outfilename = "KUCMS_GMSB_L100t400_Met150_jet1m_v6_Skim_BaseHists.root"; //0   
                //auto outfilename = "KUCMS_GMSB_L100t400_Met150_jet1m_gn_v6_Skim_BaseHists.root"; //1
                //auto outfilename = "KUCMS_GMSB_L100t400_Met150_jet1m_wz_v6_Skim_BaseHists.root"; //2
                //auto outfilename = "KUCMS_GMSB_L100t400_Met150_jet1m_prt_v6_Skim_BaseHists.root"; //3
                //auto outfilename = "KUCMS_GMSB_L100t400_Met150_jet1m_not_v6_Skim_BaseHists.root"; //4 
                //auto outfilename = "KUCMS_GMSB_L100t400_Met150_jet1m_xsr_v6_Skim_BaseHists.root"; //5
                auto outfilename = "KUCMS_GMSB_L100t400_Met150_Not2_v9_Skim_BaseHists.root"; //6

                //auto htitle = "GMSB_met150_1t4_jet1m ";
                //auto htitle = "GMSB_met150_1t4_pho1m_ng";
                //auto htitle = "GMSB_met150_1t4_jet1m_sqj";     
                //auto htitle = "GMSB_met150_1t4_jet1m_wz ";    
                //auto htitle = "GMSB_met150_1t4_jet1m_prt ";   
                //auto htitle = "GMSB_met150_1t4_jet1m_not ";  
                //auto htitle = "GMSB_met150_1t4_jet1m_xsr "; 
				auto htitle = "GMSB_met150_1t4_v9_Not2_";

                HistMaker base;
                base.histMaker( listdir, infilename, outfilename, htitle, cut );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

