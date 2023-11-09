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

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle, int cut, float va, float vb, float vc ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");

    cutselection = cut;
	cutva = va;
    cutvb = vb;
    cutvc = vc;
    preCutNPhotons = 0;
    preCut30NPhotons = 0;
    preCut100NPhotons = 0;
    postCutNPhotons = 0;
    postCut30NPhotons = 0;
    postCut100NPhotons = 0;

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
    TBranch        *b_sumEvtGenWgt;

    sKey = 0;

    fConfigTree->SetBranchAddress("nEvents", &nEvents, &b_nEvents);
    fConfigTree->SetBranchAddress("nSelectedEvents", &nSelectedEvents, &b_nSelectedEvents);
    fConfigTree->SetBranchAddress("sKey", &sKey, &b_sKey);
    fConfigTree->SetBranchAddress("sCrossSection", &sCrossSection, &b_sCrossSection);
    fConfigTree->SetBranchAddress("sGMSBGravMass", &sGMSBGravMass, &b_sGMSBGravMass);
    fConfigTree->SetBranchAddress("sGMSBChi1Mass", &sGMSBChi1Mass, &b_sGMSBChi1Mass);
    fConfigTree->SetBranchAddress("sMCWgt", &sMCWgt, &b_sMCWgt);
    fConfigTree->SetBranchAddress("sMCType", &sMCType, &b_sMCType);
    fConfigTree->SetBranchAddress("sumEvtGenWgt", &sumEvtGenWgt, &b_sumEvtGenWgt);

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
	int loopCounter(100000);
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
    bool isNotSig = ( cutselection == 11 ) ? true : false;

    //auto dskey  = *DataSetKey
    //float evtgwt = evtGenWgt;
    float evtgwt = 1;
    float scale = 10;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    float segwt = (configInfo[configKey])["nEvents"];
    auto fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
    //std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
    //auto fillwt = (configInfo[configKey])["sCrossSection"] * evtwt;
    //if( configKey != "GMSBL100" ) return;

    //Float_t         selCMet;
    //Float_t         selCMetPx;
    //Float_t         selCMetPy;
    auto metPt = std::sqrt(rad2(selCMetPx,selCMetPy));
    if( metPt < 150 ) return;
    hist1d[400]->Fill( metPt, fillwt );

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

	int totNSelPhotons = 0;
    int totNSel30Photons = 0;
    int totNSel100Photons = 0;
	int totNCutPhotons = 0;
    int totNCut30Photons = 0;
    int totNCut100Photons = 0;
    float phoeff = nPhotons ? static_cast<float>(nSelPhotons)/static_cast<float>(nPhotons) : 1.1;
    if( nPhotons > 0 ) hist1d[100]->Fill(phoeff,fillwt); //("phoEff", "phoEff", 1000, 0, 1000);
    hist1d[118]->Fill(nSelPhotons,fillwt);
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

        if( (*selPhoPt)[it] < 20 ) continue;
        if( std::abs((*selPhoEta)[it]) > 1.479 ) continue;

		bool jetphoton = false;
		for( int jit = 0; jit < nSelJets; jit++ ){

			float dpjeta = (*selJetEta)[jit] - (*selPhoPhi)[it];
			float dpjphi = dPhi( (*selJetPhi)[jit], (*selPhoPhi)[it] );	
			float dr = hypo( dpjeta, dpjphi );
			if( dr < 0.4 ) jetphoton = true;

		} // for( int jit = 0; jit < nSelJets; jit++ )
		if( jetphoton ) continue; //{ std::cout << "Jet-Photon Overlap!!!" << std::endl; continue; }

        if( DEBUG ) std::cout << " -- looping photons : Start getting dr dp info" << std::endl;
		auto phoClass = (*selPhoQuality)[it];
        //if( phoClass < 1 ) continue;
        auto phoSusId = (*selPhoSusyId)[it]; // selPhoSusyId->at(it);
        auto phoGenMatDr = (*selPhoGenDr)[it];
        auto phoGenMatDp = (*selPhoGenDp)[it];
        auto phoTime = (*selPhoTime)[it];
        auto phoOOT = (*selPhoOOT)[it];

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
        bool isNotSusyType = isPrompt || phoSusId < 9;

//        if( isgnsusy && not isNino ){ usepho = false; continue; }
//        if( iswzsusy && not isXinoWZ ){ usepho = false; continue; }
//        if( isXfsrsusy && not isXfsr ){ usepho = false; continue; }
//        if( isfsrsusy && not isSfsr ){ usepho = false; continue; }
//        if( isOthersusy && not isPrompt ){ usepho = false; continue; }

        if( isSig && not isSigType ){ usepho = false; continue; }
        if( isXfsr && not isXfsrType ){ usepho = false; continue; }
        if( isOther && not isOtherType ){ usepho = false; continue; }
        if( isUnMatched && phoSusId > 10 ){ usepho = false; continue; }
		if( isNotSig && not isNotSusyType ){ usepho = false; continue; }

        hist2d[200]->Fill( phoGenMatDr, phoTime );
        hist2d[201]->Fill( phoGenMatDp, phoTime );
        hist2d[202]->Fill( phoGenMatDr, phoGenMatDp );
        hist1d[136]->Fill( phoGenMatDp,fillwt );
        hist1d[137]->Fill( phoGenMatDr,fillwt );

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;        

        if( isNino ) nSusPhoGN1++;
        if( isXfsr ) nSusPhoGN1d++;
        if( isPrompt ) nSusPhoGN2++;
        if( isSfsr ) nSusPhoGN2d++;        
        if( isXinoWZ ) nSusPhoWN++;
        if( false ) nSusPhoZC++;
        //if( isnogen ) nPhoNoGen++;
        hist1d[135]->Fill( phoSusId,fillwt );

		totNSelPhotons++;
		if( (*selPhoPt)[it] > 30 ) totNSel30Photons++;
        if( (*selPhoPt)[it] > 100 ) totNSel100Photons++;

        bool hcalsum = true;
        //bool hcalsum = (*selPhoHcalTowerSumEtBcConeDR04)[it] < 2.45;
        bool tsptscdr4 = (*selPhoTrkSumPtSolidConeDR04)[it] < cutvb; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
		bool ecalrhsum = (*selPhoEcalRHSumEtConeDR04)[it] < cutvc;
        bool htoem = (*selPhoHadTowOverEM)[it] < cutva;
        bool isoskip = not ( htoem && tsptscdr4 && ecalrhsum && hcalsum );
        if( isoskip ) continue;
        totNCutPhotons++;

        hist2d[205]->Fill( (*selPhoEnergy)[it], (*selPhoPt)[it] );

        auto pgpt = (*selPhoGenPt)[it];

        if( DEBUG ) std::cout << " -- looping photons : filling ids 0" << std::endl;

        hist1d[101]->Fill((*selPhoClstrRn)[it],fillwt); //("phoClstrRn", "phoClstrRn", 1000, 0, 1000);
        hist1d[102]->Fill((*selPhoEta)[it],fillwt); //("phoEta", "phoEta", 700, -3.5, 3.5);
        hist1d[103]->Fill((*selPhoPhi)[it],fillwt); //("phoPhi", "phoPhi", 700, -3.5, 3.5);
        hist1d[104]->Fill((*selPhoPt)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist1d[105]->Fill((*selPhoEnergy)[it],fillwt); //("phoEnergy", "phoEnergy", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling ids 0a" << std::endl;
        hist1d[106]->Fill((*selPhoGenPt)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[204]->Fill((*selPhoGenPt)[it],(*selPhoPt)[it]);
        if( DEBUG ) std::cout << " -- looping photons : filling ids 1" << std::endl;
        //hist1d[108]->Fill((*selPhoPhoIsoDr)[it],fillwt); //("phoGeoSMaj", "phoGeoSMaj", 100, 0, 1);
        //hist2d[212]->Fill((*selPhoPhoIsoDr)[it],(*selPhoPt)[it]);
        hist1d[109]->Fill((*selPhoR9)[it],fillwt); //("phoR9", "phoR9", 100, 0, 1);
        hist1d[110]->Fill((*selPhoNrh)[it],fillwt); //("phoNrh", "phoNrh", 100, 0, 100);
        hist1d[111]->Fill((*selPhoSMaj)[it],fillwt); //("phoSMaj", "phoSMaj", 100, 0, 1);
        hist1d[112]->Fill((*selPhoSMin)[it],fillwt); //("phoSMin", "phoSMin", 100, 0, 1);
        hist1d[113]->Fill((*selPhoSieie)[it],fillwt); //("phoSieie", "phoSieie", 100, 0, 1);
        hist1d[114]->Fill((*selPhoTime)[it],fillwt); //("phoTime", "phoTime", 500, -25, 25);
        hist1d[115]->Fill((*selPhoQuality)[it],fillwt); //("phoQuality", "phoQuality", 4, 0, 4);
        hist1d[116]->Fill((*selPhoR9)[it],fillwt); //("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
        if( DEBUG ) std::cout << " -- looping photons : filling ids 2" << std::endl;
        //hist1d[117]->Fill((*selPhoGeoSMin)[it]);
        hist1d[125]->Fill((*selPhoHcalTowerSumEtBcConeDR04)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[206]->Fill((*selPhoHcalTowerSumEtBcConeDR04)[it],(*selPhoPt)[it]);
        hist1d[126]->Fill((*selPhoTrkSumPtHollowConeDR03)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[207]->Fill((*selPhoTrkSumPtHollowConeDR03)[it],(*selPhoPt)[it]);
        hist1d[127]->Fill((*selPhoTrkSumPtHollowConeDR04)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[208]->Fill((*selPhoTrkSumPtHollowConeDR04)[it],(*selPhoPt)[it]);
        hist1d[128]->Fill((*selPhoTrkSumPtSolidConeDR04)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[209]->Fill((*selPhoTrkSumPtSolidConeDR04)[it],(*selPhoPt)[it]);
        //hist1d[129]->Fill((*selPhoPixelSeed)[it]); //("phoPt", "phoPt", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling var hist 5" << std::endl;
        hist1d[130]->Fill((*selPhoEcalRHSumEtConeDR04)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[210]->Fill((*selPhoEcalRHSumEtConeDR04)[it],(*selPhoPt)[it]);
        hist1d[131]->Fill((*selPhoHadTowOverEM)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        hist2d[211]->Fill((*selPhoHadTowOverEM)[it],(*selPhoPt)[it]);
        if( DEBUG ) std::cout << " -- looping photons : filling var hist 6" << std::endl;
        hist1d[132]->Fill((*selPhoSieip)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        //if( DEBUG ) std::cout << " -- looping photons : selPhoSipip = " << (*selPhoSipip)[it] << std::endl;
        //hist1d[134]->Fill((*selPhoSipip)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
        if( DEBUG ) std::cout << " -- looping photons : filling var hist done" << std::endl;

        //if( DEBUG ) std::cout << " -- looping photons : filling lead phio" << std::endl;
        if( (*selPhoPt)[it] >= 30  ){

			totNCut30Photons++;
            float phoeffl = nPhotons ? static_cast<float>(nSelPhotons)/static_cast<float>(nPhotons) : 0;
            hist1d[150]->Fill(phoeffl); //("phoEff", "phoEff", 1000, 0, 1000);
            hist1d[151]->Fill((*selPhoClstrRn)[it],fillwt); //("phoClstrRn", "phoClstrRn", 1000, 0, 1000);
            hist1d[152]->Fill((*selPhoEta)[it],fillwt); //("phoEta", "phoEta", 700, -3.5, 3.5);
            hist1d[153]->Fill((*selPhoPhi)[it],fillwt); //("phoPhi", "phoPhi", 700, -3.5, 3.5);
            hist1d[154]->Fill((*selPhoPt)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
            hist1d[155]->Fill((*selPhoEnergy)[it],fillwt); //("phoEnergy", "phoEnergy", 500, 0, 5000);
            hist1d[157]->Fill((*selPhoHcalTowerSumEtBcConeDR04)[it],fillwt); //("phoGeoEgnVal", "phoGeoEgnVal", 100, 0, 1);
            hist1d[158]->Fill((*selPhoTrkSumPtSolidConeDR04)[it],fillwt); //("phoGeoSMaj", "phoGeoSMaj", 100, 0, 1);
            hist1d[159]->Fill((*selPhoR9)[it],fillwt); //("phoR9", "phoR9", 100, 0, 1);
            hist1d[160]->Fill((*selPhoNrh)[it],fillwt); //("phoNrh", "phoNrh", 100, 0, 100);
            hist1d[161]->Fill((*selPhoSMaj)[it],fillwt); //("phoSMaj", "phoSMaj", 100, 0, 1);
            hist1d[162]->Fill((*selPhoSMin)[it],fillwt); //("phoSMin", "phoSMin", 100, 0, 1);
            hist1d[163]->Fill((*selPhoSieie)[it],fillwt); //("phoSieie", "phoSieie", 100, 0, 1);
            hist1d[164]->Fill((*selPhoTime)[it],fillwt); //("phoTime", "phoTime", 500, -25, 25);
            hist1d[165]->Fill((*selPhoQuality)[it],fillwt); //("phoQuality", "phoQuality", 4, 0, 4);
            hist1d[166]->Fill((*selPhoR9)[it],fillwt); //("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
            hist1d[167]->Fill((*selPhoEcalRHSumEtConeDR04)[it],fillwt);
            hist1d[168]->Fill((*selPhoHadTowOverEM)[it],fillwt);

		}//<<>>if( (*selPhoPt)[it] >= 30 )

        //if( DEBUG ) std::cout << " -- looping photons : filling sub lead phio" << std::endl;
        if( (*selPhoPt)[it] >= 100 ){

			totNCut100Photons++;
            float phoeffsl = nPhotons ? static_cast<float>(nSelPhotons)/static_cast<float>(nPhotons) : 0;
            if( nPhotons > 0 ) hist1d[200]->Fill(phoeffsl); //("phoEff", "phoEff", 1000, 0, 1000);
            hist1d[201]->Fill((*selPhoClstrRn)[it],fillwt); //("phoClstrRn", "phoClstrRn", 1000, 0, 1000);
            hist1d[202]->Fill((*selPhoEta)[it],fillwt); //("phoEta", "phoEta", 700, -3.5, 3.5);
            hist1d[203]->Fill((*selPhoPhi)[it],fillwt); //("phoPhi", "phoPhi", 700, -3.5, 3.5);
            hist1d[204]->Fill((*selPhoPt)[it],fillwt); //("phoPt", "phoPt", 500, 0, 5000);
            hist1d[205]->Fill((*selPhoEnergy)[it],fillwt); //("phoEnergy", "phoEnergy", 500, 0, 5000);
            //hist1d[206]->Fill((*selPhoGeoEgnVal)[it],fillwt); //("phoGeoEgnVal", "phoGeoEgnVal", 100, 0, 1);
            hist1d[207]->Fill((*selPhoHcalTowerSumEtBcConeDR04)[it],fillwt); //("phoGeoEgnVal", "phoGeoEgnVal", 100, 0, 1);
            hist1d[208]->Fill((*selPhoTrkSumPtSolidConeDR04)[it],fillwt); //("phoGeoSMaj", "phoGeoSMaj", 100, 0, 1);
            hist1d[209]->Fill((*selPhoR9)[it],fillwt); //("phoR9", "phoR9", 100, 0, 1);
            hist1d[210]->Fill((*selPhoNrh)[it],fillwt); //("phoNrh", "phoNrh", 100, 0, 100);
            hist1d[211]->Fill((*selPhoSMaj)[it],fillwt); //("phoSMaj", "phoSMaj", 100, 0, 1);
            hist1d[212]->Fill((*selPhoSMin)[it],fillwt); //("phoSMin", "phoSMin", 100, 0, 1);
            hist1d[213]->Fill((*selPhoSieie)[it],fillwt); //("phoSieie", "phoSieie", 100, 0, 1);
            hist1d[214]->Fill((*selPhoTime)[it],fillwt); //("phoTime", "phoTime", 500, -25, 25);
            hist1d[215]->Fill((*selPhoQuality)[it],fillwt); //("phoQuality", "phoQuality", 4, 0, 4);
            hist1d[216]->Fill((*selPhoR9)[it],fillwt); //("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
            hist1d[217]->Fill((*selPhoEcalRHSumEtConeDR04)[it],fillwt);
            hist1d[218]->Fill((*selPhoHadTowOverEM)[it],fillwt);

        }//<<>>if( (*selPhoPt)[it] >= 100 )

    }//<<>>for( int it = 0; it < nPhotons; it++ )

    preCutNPhotons += totNSelPhotons;
    preCut30NPhotons += totNSel30Photons;
    preCut100NPhotons += totNSel100Photons;
    postCutNPhotons += totNCutPhotons;
    postCut30NPhotons += totNCut30Photons;
    postCut100NPhotons += totNCut100Photons;

	//-------- electrons --------------------------------------

/*
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


//        bool isgoodmatch = jetGenMatDr < 0.3 && jetGenMatDp < 1.0 && jetGenMatDr >= 0 && jetGenMatDp >= 0;
//        bool isnogen = jetSusId < 19;
//        if( isOthermatch && ( ( not isnogen ) && isgoodmatch )){ continue; } // no match
//        if( ( not isOthermatch ) && ( isnogen || ( not isgoodmatch ) )){ continue; } // good match


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

        hist1d[0]->Fill((*selJetPt)[it],fillwt); //("jetPt", mkht(ht,"jetPt"), 500, 0, 5000);
        hist1d[1]->Fill((*selJetPhi)[it],fillwt); //("jetPhi", mkht(ht,"jetPhi"), 700, -3.5, 3.5);
        hist1d[2]->Fill((*selJetEta)[it],fillwt); //("jetEta", mkht(ht,"jetEta"), 700, -3.5, 3.5);
        hist1d[3]->Fill((*selJetEnergy)[it],fillwt); //("jetEnergy", mkht(ht,"jetEnergy"), 500, 0, 5000);
        hist1d[4]->Fill((*selJetMass)[it],fillwt); //("jetMass", mkht(ht,"jetMass"), 500, 0, 5000);
        hist1d[5]->Fill((*selJetQuality)[it],fillwt); //("jetQuality", mkht(ht,"jetQuality"), 5, 0, 5);
        hist1d[6]->Fill((*selJetTime)[it],fillwt); //("jetTime", mkht(ht,"jetTime"), 500, -25, 25);

        hist1d[30]->Fill((*selGenJetDpt)[it],fillwt); //("genJetDpt", mkht(ht,"genJetDpt").c_str(), 500, 0, 5);
        hist1d[31]->Fill((*selGenJetEnergy)[it],fillwt); //("genJetEnergy", mkht(ht,"genJetEnergy").c_str(), 5000, 0, 5000);
        hist1d[32]->Fill((*selGenJetImpAng)[it],fillwt); //("genJetImpAng", mkht(ht,"genJetImpAng").c_str(), 200, -10, 10);
        hist1d[33]->Fill((*selGenJetLlpTime)[it],fillwt); //("genJetLlpTime", mkht(ht,"genJetLlpTime").c_str(), 100, -5, 5);
        hist1d[34]->Fill((*selGenJetPt)[it],fillwt); //("genJetPt", mkht(ht,"genJetPt").c_str(), 3000, 0, 3000);
        hist1d[35]->Fill((*selGenJetTime)[it],fillwt); //("genJetTime", mkht(ht,"genJetTime").c_str(), 250, 0, 25);
        hist1d[36]->Fill((*selGenJetTof)[it],fillwt); //("genJetTof", mkht(ht,"genJetTof").c_str(), 400, 0, 400);
        hist1d[37]->Fill((*selGenJetdr)[it],fillwt); //("genJetdr", mkht(ht,"genJetdr").c_str(), 600, -1, 5);
        hist1d[38]->Fill((*selGenJeteta)[it],fillwt); //("genJeteta", mkht(ht,"genJeteta").c_str(), 160, -10, 6);

        hist1d[39]->Fill((*selJetSusyId)[it],fillwt); //("jetSusyId", mkht(ht,"jetSusyId").c_str(), 45, -5, 40);
        hist1d[40]->Fill((*selJetLlpDp)[it],fillwt); //("jetLlpDp", mkht(ht,"jetLlpDp").c_str(), 500, 0, 5);
        hist1d[41]->Fill((*selJetLlpDr)[it],fillwt); //("jetLlpDr", mkht(ht,"jetLlpDr").c_str(), 500, 0, 5);
        hist2d[203]->Fill((*selJetLlpDr)[it],(*selJetLlpDp)[it]);

        hist1d[42]->Fill((*selJetArea)[it],fillwt); //("jetArea", mkht(ht,"jetArea").c_str(), 100, 0, 1);
        hist1d[43]->Fill((*selJetChEmEF)[it],fillwt); //("jetChEmEF", mkht(ht,"jetChEmEF").c_str(), 100, 0, 1);
        hist1d[44]->Fill((*selJetChHM)[it],fillwt); //("jetChHM", mkht(ht,"jetChHM").c_str(), 100, 0, 100);
        hist1d[45]->Fill((*selJetMuEF)[it],fillwt); //("jetMuEF", mkht(ht,"jetMuEF").c_str(), 100, 0, 1);
        hist1d[46]->Fill((*selJetNeEmEF)[it],fillwt); //("jetNeEmEF", mkht(ht,"jetNeEmEF").c_str(), 100, 0, 1);
        hist1d[47]->Fill((*selJetNeHEF)[it],fillwt); //("jetNeHEF", mkht(ht,"jetNeHEF").c_str(), 100, 0, 1);
        hist1d[48]->Fill((*selJetNeHM)[it],fillwt); //("jetNeHM", mkht(ht,"jetNeHM").c_str(), 75, 0, 75);
        hist1d[49]->Fill((*selJetchHEF)[it],fillwt); //("jetchHEF", mkht(ht,"jetchHEF").c_str(), 100, 0, 1);


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
*/

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

	float eff =  static_cast<float>(postCutNPhotons)/static_cast<float>(preCutNPhotons); 
    float eff30 =  static_cast<float>(postCut30NPhotons)/static_cast<float>(preCut30NPhotons);
    float eff100 =  static_cast<float>(postCut100NPhotons)/static_cast<float>(preCut100NPhotons);
    std::cout << " Pho Eff : " << eff << " = " << postCutNPhotons << " / " << preCutNPhotons << std::endl;
    std::cout << " Pho Eff 30 : " << eff30 << " = " << postCut30NPhotons << " / " << preCut30NPhotons << std::endl;
    std::cout << " Pho Eff 100 : " << eff100 << " = " << postCut100NPhotons << " / " << preCut100NPhotons << std::endl;
    hist1d[197]->Fill(1,eff);
    hist1d[197]->Fill(2,eff30);
    hist1d[197]->Fill(3,eff100);
    hist1d[198]->Fill(1,preCutNPhotons);
    hist1d[199]->Fill(1,postCutNPhotons);
    hist1d[199]->Fill(2,postCut30NPhotons);
    hist1d[199]->Fill(3,postCut100NPhotons);

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

/*
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
*/

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
    hist1d[106] = new TH1D("phoGenPt", mkht(ht,"phoGenPt").c_str(), 750, 0, 750);
    hist1d[108] = new TH1D("phoPhoIsoDr", mkht(ht,"phoPhoIsoDr").c_str(), 100, 0, 10);
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

    //hist1d[119] = new TH1D("nGN1Photons", mkht(ht,"nGN1Photons").c_str(), 20, 0, 20);
    //hist1d[120] = new TH1D("nGN1dPhotons", mkht(ht,"nGN1dPhotons").c_str(), 20, 0, 20);
    //hist1d[121] = new TH1D("nGN2Photons", mkht(ht,"nGN2Photons").c_str(), 20, 0, 20);
    //hist1d[122] = new TH1D("nGN2dPhotons", mkht(ht,"nGN2dPhotons").c_str(), 20, 0, 20);
    //hist1d[123] = new TH1D("nWNPhotons", mkht(ht,"nWNPhotons").c_str(), 20, 0, 20);
    //hist1d[124] = new TH1D("nZCPhotons", mkht(ht,"nZCPhotons").c_str(), 20, 0, 20);

    hist1d[125] = new TH1D("selPhoHcalTowerSumEtBcConeDR04", mkht(ht,"selPhoHcalTowerSumEtBcConeDR04").c_str(), 200, 0, 20);
    hist1d[126] = new TH1D("selPhoTrkSumPtHollowConeDR03", mkht(ht,"selPhoTrkSumPtHollowConeDR03").c_str(), 200, 0, 20);
    hist1d[127] = new TH1D("selPhoTrkSumPtHollowConeDR04", mkht(ht,"selPhoTrkSumPtHollowConeDR04").c_str(), 200, 0, 20);
    hist1d[128] = new TH1D("selPhoTrkSumPtSolidConeDR04", mkht(ht,"selPhoTrkSumPtSolidConeDR04").c_str(), 200, 0, 20);
    hist1d[129] = new TH1D("selPhoPixelSeed", mkht(ht,"selPhoPixelSeed").c_str(), 2, 0, 1);
    hist1d[130] = new TH1D("selPhoEcalRHSumEtConeDR04", mkht(ht,"selPhoEcalRHSumEtConeDR04").c_str(), 200, 0, 20);
    hist1d[131] = new TH1D("selPhoHadTowOverEM", mkht(ht,"selPhoHadTowOverEM").c_str(), 200, 0, 0.2);
    hist1d[132] = new TH1D("selPhoSieip", mkht(ht,"selPhoSieip").c_str(), 100, 0, 0.01);
    //hist1d[134] = new TH1D("selPhoSipip", mkht(ht,"selPhoSipip").c_str(), 100, 0, 0.01);

    hist1d[134] = new TH1D("nNoGenMatch", mkht(ht,"nNoGenMatch").c_str(), 20, 0, 20);
    hist1d[135] = new TH1D("phoGenType", mkht(ht,"phoGenType").c_str(), 45, -5, 40);
    hist1d[136] = new TH1D("phoGenRe", mkht(ht,"phoGenRe").c_str(), 500, 0, 5);
    hist1d[137] = new TH1D("phoGenDr", mkht(ht,"phoGenDr").c_str(), 50, 0, 0.5);

    hist1d[150] = new TH1D("pho30ptEff", mkht(ht,"pho30ptEff").c_str(), 10, 0, 1);
    hist1d[151] = new TH1D("pho30ptClstrRn", mkht(ht,"pho30ptClstrRn").c_str(), 100, 0, 1);
    hist1d[152] = new TH1D("pho30ptEta", mkht(ht,"pho30ptEta").c_str(), 700, -3.5, 3.5);
    hist1d[153] = new TH1D("pho30ptPhi", mkht(ht,"pho30ptPhi").c_str(), 700, -3.5, 3.5);
    hist1d[154] = new TH1D("pho30ptPt", mkht(ht,"pho30ptPt").c_str(), 200, 0, 2000);
    hist1d[155] = new TH1D("pho30ptEnergy", mkht(ht,"pho30ptEnergy").c_str(), 500, 0, 5000);
    hist1d[157] = new TH1D("pho30ptHcalTowerSumEtBcConeDR04", mkht(ht,"pho30ptHcalTowerSumEtBcConeDR04").c_str(), 200, 0, 20);
    hist1d[158] = new TH1D("pho30ptTrkSumPtSolidConeDR04", mkht(ht,"pho30ptTrkSumPtSolidConeDR04").c_str(), 200, 0, 20);
    hist1d[159] = new TH1D("pho30ptR9", mkht(ht,"pho30ptR9").c_str(), 100, 0, 1);
    hist1d[160] = new TH1D("pho30ptNrh", mkht(ht,"pho30ptNrh").c_str(), 100, 0, 100);
    hist1d[161] = new TH1D("pho30ptSMaj", mkht(ht,"pho30ptSMaj").c_str(), 20, 0, 2);
    hist1d[162] = new TH1D("pho30ptSMin", mkht(ht,"pho30ptSMin").c_str(), 100, 0, 1);
    hist1d[163] = new TH1D("pho30ptSieie", mkht(ht,"pho30ptSieie").c_str(), 100, 0, 0.01);
    hist1d[164] = new TH1D("pho30ptTime", mkht(ht,"pho30ptTime").c_str(), 500, -25, 25);
    hist1d[165] = new TH1D("pho30ptQuality", mkht(ht,"pho30ptQuality").c_str(), 4, 0, 4);
    hist1d[166] = new TH1D("pho30ptR9_zoom", mkht(ht,"pho30ptR9_zoom").c_str(), 200, 0.8, 1);
    hist1d[167] = new TH1D("pho30ptEcalRHSumEtConeDR04", mkht(ht,"pho30ptEcalRHSumEtConeDR04").c_str(), 200, 0, 20);
    hist1d[168] = new TH1D("pho30ptHadTowOverEM", mkht(ht,"pho30ptHadTowOverEM").c_str(), 200, 0, 0.2);

    hist1d[200] = new TH1D("pho100ptEff", mkht(ht,"pho100ptEff").c_str(), 10, 0, 1);
    hist1d[201] = new TH1D("pho100ptClstrRn", mkht(ht,"pho100ptClstrRn").c_str(), 100, 0, 1);
    hist1d[202] = new TH1D("pho100ptEta", mkht(ht,"pho100ptEta").c_str(), 700, -3.5, 3.5);
    hist1d[203] = new TH1D("pho100ptPhi", mkht(ht,"pho100ptPhi").c_str(), 700, -3.5, 3.5);
    hist1d[204] = new TH1D("pho100ptPt", mkht(ht,"pho100ptPt").c_str(), 200, 0, 2000);
    hist1d[205] = new TH1D("pho100ptEnergy", mkht(ht,"pho100ptEnergy").c_str(), 500, 0, 5000);
    hist1d[207] = new TH1D("pho100ptHTSumEtBcConeDR04", mkht(ht,"pho100ptHcalTowerSumEtBcConeDR04").c_str(), 200, 0, 20);
    hist1d[208] = new TH1D("pho100ptTrkSumPtSolidConeDR04", mkht(ht,"pho100ptTrkSumPtSolidConeDR04").c_str(), 200, 0, 20);
    hist1d[209] = new TH1D("pho100ptR9", mkht(ht,"pho100ptR9").c_str(), 100, 0, 1);
    hist1d[210] = new TH1D("pho100ptNrh", mkht(ht,"pho100ptNrh").c_str(), 100, 0, 100);
    hist1d[211] = new TH1D("pho100ptSMaj", mkht(ht,"pho100ptSMaj").c_str(), 20, 0, 2);
    hist1d[212] = new TH1D("pho100ptSMin", mkht(ht,"pho100ptSMin").c_str(), 100, 0, 1);
    hist1d[213] = new TH1D("pho100ptSieie", mkht(ht,"pho100ptSieie").c_str(), 100, 0, 0.01);
    hist1d[214] = new TH1D("pho100ptTime", mkht(ht,"pho100ptTime").c_str(), 500, -25, 25);
    hist1d[215] = new TH1D("pho100ptQuality", mkht(ht,"pho100ptQuality").c_str(), 4, 0, 4);
    hist1d[216] = new TH1D("pho100ptR9_zoom", mkht(ht,"pho100ptR9_zoom").c_str(), 200, 0.8, 1);
    hist1d[217] = new TH1D("pho100ptEcalRHSumEtConeDR04", mkht(ht,"pho100ptEcalRHSumEtConeDR04").c_str(), 200, 0, 20);
    hist1d[218] = new TH1D("pho100ptHadTowOverEM", mkht(ht,"pho100ptHadTowOverEM").c_str(), 200, 0, 0.2);


	float cstart(0), cend(20), cdiv(20);
    hist1d[197] = new TH1D("phoCEffPhoton", mkht(ht,"phoCEffPhoton").c_str(), cdiv, cstart, cend);
    hist1d[198] = new TH1D("phoNPreCutPhoton", mkht(ht,"phoNPreCutPhoton").c_str(), cdiv, cstart, cend);
    hist1d[199] = new TH1D("phoNPostCutPhoton", mkht(ht,"phoNPostCutPhoton").c_str(), cdiv, cstart, cend);
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
    hist2d[204] = new TH2D("phoGenPtVPt", mkht(ht,"pho genPt v Pt; genPt [GeV]; Pt [GeV]").c_str(), 1000, 0, 1000, 1000, 0, 1000);
    // ....
    hist2d[205] = new TH2D("phoEvPt", mkht(ht,"pho E v Pt;Energy [GeV];Pt [GeV]").c_str(), 2000, 0, 2000, 1000, 0, 1000);

    hist2d[206] = new TH2D("phoEvHtsebcdr04", mkht(ht,"pho Htsebcdr04 v Pt;HcalTowerSumEtBcConeDR04;Pt [GeV]").c_str(), 200, 0, 20, 1000, 0, 1000);
    hist2d[207] = new TH2D("phoEvTspthcdr03", mkht(ht,"pho Tspthcdr03 v Pt;TrkSumPtHollowConeDR03;Pt [GeV]").c_str(), 200, 0, 20, 1000, 0, 1000);
    hist2d[208] = new TH2D("phoEvTspthcdr04", mkht(ht,"pho Tspthcdr04 v Pt;TrkSumPtHollowConeDR04;Pt [GeV]").c_str(), 200, 0, 20, 1000, 0, 1000);
    hist2d[209] = new TH2D("phoEvTsptscdr04", mkht(ht,"pho Tsptscdr04 v Pt;TrkSumPtSolidConeDR04;Pt [GeV]").c_str(), 200, 0, 20, 1000, 0, 1000);
    hist2d[210] = new TH2D("phoEvErhsecdr04", mkht(ht,"pho Erhsecdr04 v Pt;EcalRHSumEtConeDR04;Pt [GeV]").c_str(), 200, 0, 20, 1000, 0, 1000);
    hist2d[211] = new TH2D("phoEvHtoem", mkht(ht,"pho Htoem v Pt;HadTowOverEM;Pt [GeV]").c_str(), 100, 0, 0.1, 1000, 0, 1000);
    hist2d[212] = new TH2D("phoPhoIsovPt", mkht(ht,"pho Iso Dr v Pt;pho Iso Dr; Pt [GeV}").c_str(), 100, 0, 10, 1000, 0, 1000);

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

/*
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

	iso1 = selPhoTrkSumPtSolidConeDR04 a = 10, b = 12.5, c = 15

*/

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "skims_files/";

                //auto infilenameGMSB = "KUCMS_GMSB_Skim_List.txt";
                auto infilenameGMSB = "KUCMS_GMSB_Skim_List2.txt";
                auto infilenameGJets = "KUCMS_GJets_Skim_List.txt";
                auto infilenameJetHT = "KUCMS_JetHT_Skim_List.txt";

                int cuts = 7;
				int cutb = 11;
                float value0 = 10000.0;
                // -- ecalRHSumEtConeDR04 (1)
                float vcP = 3.6;// 0.001*(*Photon_pt)[it] + 3.5(L)2.0(T)
				std::list<float> cutsvc = { 8.0, 9.0, 10.0 };
                //float value1 = 7.0;
                //float value2 = 6.0;
                //float value3 = 5.0;
				// -- selPhoTrkSumPtSolidConeDR04 (2) 
                float vbP = 4.8;// 0.001*(*Photon_pt)[it] + 3.5(L)2.0(T)
				std::list<float> cutsvb = { 5.0, 6.0, 7.0 };
                //float value1 = 13.0;
                //float value2 = 15.0;
                //float value3 = 17.0;
				// -- selPhoHadTowOverEM (1)
                float vaP = 0.05;// 0.001*(*Photon_pt)[it] + 3.5(L)2.0(T)
                std::list<float> cutsva = { 0.01, 0.02, 0.03 };
                //float value1 = 0.03;
                //float value2 = 0.02;
                //float value3 = 0.01;

                //auto outfilenamePs = "KUCMS_GMSB_L100_Met150_Signal_v15_iso1_Skim_BaseHists.root"; //7
                //auto outfilename1s = "KUCMS_GMSB_L100_Met150_Signal_v15_iso2_Skim_BaseHists.root"; //7
                //auto outfilename2s = "KUCMS_GMSB_L100_Met150_Signal_v15_iso3_Skim_BaseHists.root"; //7
                //auto outfilename3s = "KUCMS_GMSB_L100_Met150_Signal_v15_iso4_Skim_BaseHists.root"; //7

                //auto outfilenamePb = "KUCMS_GMSB_L100_Met150_notSig_v15_iso1_Skim_BaseHists.root"; //11
                //auto outfilename1b = "KUCMS_GMSB_L100_Met150_notSig_v15_iso2_Skim_BaseHists.root"; //11
                //auto outfilename2b = "KUCMS_GMSB_L100_Met150_notSig_v15_iso3_Skim_BaseHists.root"; //11
                //auto outfilename3b = "KUCMS_GMSB_L100_Met150_notSig_v15_iso4_Skim_BaseHists.root"; //11
                //auto outfilenamePd = "KUCMS_JetHt_18D_Met150_notSig_v15_iso1_Skim_BaseHists.root"; //11
                //auto outfilename1d = "KUCMS_JetHt_18D_Met150_notSig_v15_iso2_Skim_BaseHists.root"; //11
                //auto outfilename2d = "KUCMS_JetHt_18D_Met150_notSig_v15_iso3_Skim_BaseHists.root"; //11
                //auto outfilename3d = "KUCMS_JetHt_18D_Met150_notSig_v15_iso4_Skim_BaseHists.root"; //11


                //auto outfilename = "KUCMS_GMSB_L100_Met150_XFSR_v13_Skim_BaseHists.root"; //8
                //auto outfilename = "KUCMS_GMSB_L100_Met150_Other_v15_Skim_BaseHists.root"; //9
                //auto outfilename = "KUCMS_GMSB_L100_Met150_UnMatched_v13_Skim_BaseHists.root"; //10
                //auto outfilename0b = "KUCMS_GMSB_L100_Met150_notSig_v15_Skim_BaseHists.root"; //11

                //auto outfilename = "KUCMS_GJets_HT40tInf_Met150_Other_v10a_Skim_BaseHists.root"; //9
                //auto outfilename = "KUCMS_GJets_HT40tInf_Met150_UnMatched_v10a_Skim_BaseHists.root"; //10

                //auto outfilename = "KUCMS_JetHT_Met150_Other_v15_Skim_BaseHists.root"; //9

				//auto htitle = "GMSB_met150_1t4_v9_Not2_";

                //auto htitle0s = "GMSB_met150_L100_v15_iso0_Signal_";
                //auto htitlePs = "GMSB_met150_L100_v15_iso1_Signal_";
                //auto htitle1s = "GMSB_met150_L100_v15_iso2_Signal_";
                //auto htitle2s = "GMSB_met150_L100_v15_iso3_Signal_";
                //auto htitle3s = "GMSB_met150_L100_v15_iso4_Signal_";

                //auto htitle0b = "GMSB_met150_L100_v15_iso0_notSig_";
                //auto htitlePb = "GMSB_met150_L100_v15_iso1_notSig_";
                //auto htitle1b = "GMSB_met150_L100_v15_iso2_notSig_";
                //auto htitle2b = "GMSB_met150_L100_v15_iso3_notSig_";
                //auto htitle3b = "GMSB_met150_L100_v15_iso4_notSig_";

                //auto htitle0d = "JetHT_met150_18D_v15_iso0_notSig_";
                //auto htitlePd = "JetHT_met150_18D_v15_iso1_notSig_";
                //auto htitle1d = "JetHT_met150_18D_v15_iso2_notSig_";
                //auto htitle2d = "JetHT_met150_18D_v15_iso3_notSig_";
                //auto htitle3d = "JetHT_met150_18D_v15_iso4_notSig_";

                //auto htitle = "GJets_met150_40tInf_Other_v10_";
                //auto htitle = "GJets_met150_40tInf_UnMatched_v10_";

                //auto htitle = "JetHT_met150_L100_v15_Other_";

                std::string outfilename0s = "KUCMS_GMSB_L100_Met150_Signal_v15_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename0b = "KUCMS_GMSB_L100_Met150_notSig_v15_"; //iso0_Skim_BaseHists.root"; //11
                std::string outfilename0d = "KUCMS_JetHt_18D_Met150_notSig_v15_"; //iso0_Skim_BaseHists.root"; //11                
                std::string ofnending = "Skim_BaseHists.root";

                std::string htitle0s = "GMSB_met150_L100_v15_Signal_";
                std::string htitle0b = "GMSB_met150_L100_v15_notSig_";
                std::string htitle0d = "JetHT_met150_18D_v15_notSig_";

                HistMaker base;

				//std::string outfilenames = outfilename0s + "iso_005_360_480_245_" + ofnending;
                //std::string outfilenameb = outfilename0b + "iso_005_360_480_245_" + ofnending;
                //std::string outfilenamed = outfilename0d + "iso_005_360_480_245_" + ofnending; 				

                //std::string htitles =  htitle0s + "iso_005_360_480_245_";
                //std::string htitleb =  htitle0b + "iso_005_360_480_245_";
                //std::string htitled =  htitle0d + "iso_005_360_480_245_";

                //base.histMaker( listdir, infilenameGMSB, outfilenames, htitles, cuts, vaP, vbP, vcP );
                //base.histMaker( listdir, infilenameGMSB, outfilenameb, htitleb, cutb, vaP, vbP, vcP );
                //base.histMaker( listdir, infilenameJetHT, outfilenamed, htitled, cutb, vaP, vbP, vcP );

				for( auto va : cutsva ){
					for( auto vb : cutsvb ){
						for( auto vc : cutsvc ){

							std::string astr = std::to_string( int(va * 100) );
                            std::string bstr = std::to_string( int(vb * 10) );
                            std::string cstr = std::to_string( int(vc * 10) );
							//std::cout << " - " << astr << " " << bstr << " " << cstr << std::endl;
							std::string isostr = "iso_" + astr + "_" + bstr + "_" + cstr + "_";
							//std::cout << " - " << isostr << std::endl;
							std::string htitles =  htitle0s + isostr;
							std::string htitleb =  htitle0b + isostr;
							std::string htitled =  htitle0d + isostr;
							//std::cout << " - " << htitles << std::endl;
                			std::string outfilenames = outfilename0s + isostr + ofnending;
            			    std::string outfilenameb = outfilename0b + isostr + ofnending;
			                std::string outfilenamed = outfilename0d + isostr + ofnending;
                            //std::cout << " - " << outfilenames << std::endl;
                			base.histMaker( listdir, infilenameGMSB, outfilenames, htitles, cuts, va, vb, vc );
                			base.histMaker( listdir, infilenameGMSB, outfilenameb, htitleb, cutb, va, vb, vc );
                			base.histMaker( listdir, infilenameJetHT, outfilenamed, htitled, cutb, va, vb, vc );

						}//for( vc : cutsvc )
					}//for( vb : cutsvb )
				}//for( va : cutsva )
			

                //base.histMaker( listdir, infilenameGMSB, outfilenamePs, htitlePs, cuts, valueP );
                //base.histMaker( listdir, infilenameGMSB, outfilename1s, htitle1s, cuts, value1 );
                //base.histMaker( listdir, infilenameGMSB, outfilename2s, htitle2s, cuts, value2 );
                //base.histMaker( listdir, infilenameGMSB, outfilename3s, htitle3s, cuts, value3 );

                //base.histMaker( listdir, infilenameGMSB, outfilenamePb, htitlePb, cutb, valueP );
                //base.histMaker( listdir, infilenameGMSB, outfilename1b, htitle1b, cutb, value1 );
                //base.histMaker( listdir, infilenameGMSB, outfilename2b, htitle2b, cutb, value2 );
                //base.histMaker( listdir, infilenameGMSB, outfilename3b, htitle3b, cutb, value3 );

                //base.histMaker( listdir, infilenameJetHT, outfilenamePd, htitlePd, cutb, valueP );
                //base.histMaker( listdir, infilenameJetHT, outfilename1d, htitle1d, cutb, value1 );
                //base.histMaker( listdir, infilenameJetHT, outfilename2d, htitle2d, cutb, value2 );
                //base.histMaker( listdir, infilenameJetHT, outfilename3d, htitle3d, cutb, value3 );

    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

