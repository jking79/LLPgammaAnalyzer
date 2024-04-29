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

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle, int cut, float va, float vb, float vc, float vd ){

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
    cutvd = vd;
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
        if( str == " " ) continue;
		if( str == "" ) continue;
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

    //auto dskey  = *DataSetKey
    //float evtgwt = evtGenWgt;
    float evtgwt = 1;
    float scale = 10;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    float segwt = (configInfo[configKey])["nEvents"];
    auto fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
	if( cutselection == 1 ) fillwt = 1;
    //std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
    //auto fillwt = (configInfo[configKey])["sCrossSection"] * evtwt;
    //if( configKey != "GMSBL"+lambda+"" ) return;

/*
    //Float_t         selCMet;
    //Float_t         selCMetPx;
    //Float_t         selCMetPy;
   Float_t         selCMet;
   Float_t         selCMetPx;
   Float_t         selCMetPy;
   Float_t         selMet;
   Float_t         selMetPx;
   Float_t         selMetPy;
*/

	for( int i = 0; i < rjrNRjrPhotons->size(); i++ ){

	bool genSigTypeCut( genSigType == cutva );
	bool nRjrPhoCut( rjrNRjrPhotons->at(0) == cutvb );

	if( cutva == 0 ) genSigTypeCut = true;
	if( cutvb == 0 ) nRjrPhoCut = true;
	if( not genSigTypeCut || not nRjrPhoCut ) continue;

    auto metCPt = std::sqrt(rad2(selCMetPx,selCMetPy));
    auto metPt = std::sqrt(rad2(selMetPx,selMetPy));

    hist1d[400]->Fill( selCMet, fillwt );
    hist1d[401]->Fill( metCPt, fillwt );
    hist1d[402]->Fill( rjrMET->at(0), fillwt );
    hist1d[403]->Fill( metPt, fillwt );

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

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
/*

   vector<int>     *genCharge;
   vector<float>   *genMass;
   vector<float>   *genPartEnergy;
   vector<float>   *genPartEta;
   vector<unsigned int> *genPartPdgId;
   vector<float>   *genPartPhi;
   vector<float>   *genPartPt;
   vector<int>     *genPartSusId;
   vector<float>   *genPx;
   vector<float>   *genPy;
   vector<float>   *genPz;
   vector<float>   *genSGMass;
   vector<float>   *genSQMass;
   Int_t           genSigType;
   Bool_t          genStatus;
   vector<float>   *genVx;
   vector<float>   *genVy;
   vector<float>   *genVz;

*/

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


    if( DEBUG ) std::cout << "Finding genjet" << std::endl;
	//// --- genjet info -----------------------------------
*/

    if( DEBUG ) std::cout << "Finding rjr" << std::endl;

/*

   vector<float>   *rjrDphiMETV;
   vector<float>   *rjrEVa;
   vector<float>   *rjrEVb;
   vector<float>   *rjrMET;
   vector<float>   *rjrMV;
   vector<float>   *rjrMVa;
   vector<float>   *rjrMVb;
   Int_t           rjrNJetsJa;
   Int_t           rjrNJetsJb;
   vector<int>     *rjrNRjrJets;
   vector<int>     *rjrNRjrPhotons;
   vector<float>   *rjrPV;
   vector<float>   *rjrPVa;
   vector<float>   *rjrPVb;
   vector<float>   *rjrPVlab;
   vector<float>   *rjrSCosA;
   vector<float>   *rjrSMass;
   vector<float>   *rjrSdphiBV;
   vector<float>   *rjrSdphiDA;
   vector<float>   *rjrX1aCosA;
   vector<float>   *rjrX1aMass;
   vector<float>   *rjrX1bCosA;
   vector<float>   *rjrX1bMass;
   vector<float>   *rjrX2Ave;
   vector<float>   *rjrX2Diff;
   vector<float>   *rjrX2GMean;
   vector<float>   *rjrX2NAve;
   vector<float>   *rjrX2NGMean;
   vector<float>   *rjrX2NVSum;
   vector<float>   *rjrX2VSum;
   vector<float>   *rjrX2aCosA;
   vector<float>   *rjrX2aMass;
   vector<float>   *rjrX2aPtS;
   vector<float>   *rjrX2bCosA;
   vector<float>   *rjrX2bMass;
   vector<float>   *rjrX2bPtS;

//rjrDphiMETV
////rjrMET
////rjrNJetsJa
////rjrNJetsJb
////rjrNRjrPhotons
////rjrSdphiBV
////rjrSdphiDA
////rjrX2Ave
////rjrX2Diff
////rjrX2GMean
////rjrX2NAve
////rjrX2NGMean
////rjrX2NVSum
////rjrX2VSum
////rjrX2aMass
////rjrX2bMass
////rjrX2aPtS
////rjrX2bPtS
////
////genSigType
////

*/

	hist1d[250]->Fill((*rjrSCosA)[0], fillwt ); //("SCosA", mkht(ht,"SCosA"), 70, -3.5, 3.5);
    hist1d[251]->Fill((*rjrSMass)[0], fillwt ); //("SMass", mkht(ht,"SMass"), 500, 0, 5000);
    hist1d[252]->Fill((*rjrX2aMass)[0], fillwt ); //("X2aMass", mkht(ht,"X2aMass"), 500, 0, 5000);
    hist1d[253]->Fill((*rjrX2bMass)[0], fillwt ); //("X2bMass", mkht(ht,"X2bMass"), 500, 0, 5000);

    hist1d[254]->Fill((*rjrDphiMETV)[0], fillwt ); 
    hist1d[255]->Fill(rjrNJetsJa, fillwt ); 
    hist1d[256]->Fill(rjrNJetsJb, fillwt ); 
    hist1d[257]->Fill((*rjrSdphiBV)[0], fillwt ); 
    hist1d[258]->Fill((*rjrSdphiDA)[0], fillwt ); 
    hist1d[259]->Fill((*rjrX2aPtS)[0], fillwt ); 
    hist1d[260]->Fill((*rjrX2bPtS)[0], fillwt ); 

    hist1d[261]->Fill((*rjrX2Ave)[0], fillwt ); 
    hist1d[262]->Fill((*rjrX2Diff)[0], fillwt ); 
    hist1d[263]->Fill((*rjrX2GMean)[0], fillwt ); 
    hist1d[264]->Fill((*rjrX2NAve)[0], fillwt ); 
    hist1d[265]->Fill((*rjrX2NGMean)[0], fillwt ); 
    hist1d[266]->Fill((*rjrX2NVSum)[0], fillwt ); 
    hist1d[267]->Fill((*rjrX2VSum)[0], fillwt ); 

	}//for( int i = 0; i < nRjr; i++ )

}//<<>>void HistMaker::eventLoop(` Long64_t entry )

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

	//------  genparticles 200 - 249

	//------  rjr 250 - 299
//rjrDphiMETV
//rjrMET
//rjrNJetsJa
//rjrNJetsJb
//rjrNRjrPhotons
//rjrSdphiBV
//rjrSdphiDA

//rjrX2Ave
//rjrX2Diff
//rjrX2GMean
//rjrX2NAve
//rjrX2NGMean
//rjrX2NVSum
//rjrX2VSum
//
//rjrX2aMass
//rjrX2bMass
//rjrX2aPtS
//rjrX2bPtS
//
//genSigType
//


    hist1d[250] = new TH1D("SCosA", mkht(ht,"SCosA").c_str(), 70, -3.5, 3.5);
    hist1d[251] = new TH1D("SMass", mkht(ht,"SMass").c_str(), 100, 0, 10000);
    hist1d[252] = new TH1D("X2aMass", mkht(ht,"X2aMass").c_str(), 90, 0, 3600);
    hist1d[253] = new TH1D("X2bMass", mkht(ht,"X2bMass").c_str(), 90, 0, 3600);

    hist1d[254] = new TH1D("rjrDphiMETV", mkht(ht,"rjrDphiMETV").c_str(), 70, -3.5, 3.5);
    hist1d[255] = new TH1D("rjrNJetsJa", mkht(ht,"rjrNJetsJa").c_str(), 20, 0, 20);
    hist1d[256] = new TH1D("rjrNJetsJb", mkht(ht,"rjrNJetsJb").c_str(), 20, 0, 20);
    hist1d[257] = new TH1D("rjrSdphiBV", mkht(ht,"rjrSdphiBV").c_str(), 35, 0, 3.5);
    hist1d[258] = new TH1D("rjrSdphiDA", mkht(ht,"rjrSdphiDA").c_str(), 35, 0, 3.5);
    hist1d[259] = new TH1D("rjrX2aPtS", mkht(ht,"rjrX2aPtS").c_str(), 100, 0, 6000);
    hist1d[260] = new TH1D("rjrX2bPtS", mkht(ht,"rjrX2bPtS").c_str(), 100, 0, 6000);

    hist1d[261] = new TH1D("rjrX2Ave", mkht(ht,"rjrX2Ave").c_str(), 90, 0, 3600);
    hist1d[262] = new TH1D("rjrX2Diff", mkht(ht,"rjrX2Diff").c_str(), 50, 0, 1);
    hist1d[263] = new TH1D("rjrX2GMean", mkht(ht,"rjrX2GMean").c_str(), 90, 0, 3600);
    hist1d[264] = new TH1D("rjrX2NAve", mkht(ht,"rjrX2NAve").c_str(), 50, 0, 1);
    hist1d[265] = new TH1D("rjrX2NGMean", mkht(ht,"rjrX2NGMean").c_str(), 100, 0, 1);
    hist1d[266] = new TH1D("rjrX2NVSum", mkht(ht,"rjrX2NVSum").c_str(), 50, 0, 1);
    hist1d[267] = new TH1D("rjrX2VSum", mkht(ht,"rjrX2VSum").c_str(), 90, 0, 3600);

    //------ ecal rechits 300 - 349

    //------  electrons 350 - 400

	//-------- event vars 400 - 450
    hist1d[400] = new TH1D("evtMetE", mkht(ht,"evtMetE").c_str(), 150, 0, 1500 );
    hist1d[401] = new TH1D("evtMetPt", mkht(ht,"evtMetPt ").c_str(), 150, 0, 1500 );	
    hist1d[402] = new TH1D("evtRjrMetPt", mkht(ht,"evtRjrMetPt ").c_str(), 150, 0, 1500 );
    hist1d[403] = new TH1D("evtRawMetPt", mkht(ht,"evtRawMetPt ").c_str(), 150, 0, 1500 );

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	//------- jets ( time ) 0-49 ------------------------------

	//---jet id stuff 50 - 99 ---------------------------------------------------

	//--- Photons 200 - 349 -------------------------------------------
	

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
                //const std::string listdir = "skims_files/";
                const std::string listdir = "rjr_study_files/";

                auto infilename = "KUCMS_RJR_Skim_List.txt";

                //std::list<float> cuts = { 10000000 }; 

			 	std::string lambda = "300";

                std::string outfilename = "KUCMS_GJets_v19_"; //iso0_Skim_BaseHists.root"; //7
                //std::string outfilename = "KUCMS_QCD_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename0 = "KUCMS_GMSB_L"+lambda+"_T22_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename1 = "KUCMS_GMSB_L"+lambda+"_T11_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename2 = "KUCMS_GMSB_L"+lambda+"_T44_v19_"; //iso0_Skim_BaseHists.root"; //7
                std::string outfilename3 = "KUCMS_GMSB_L"+lambda+"_T33_v19_"; //iso0_Skim_BaseHists.root"; //7

                std::string ofnending = "RjrSkim_multiHists.root";

                std::string htitle = "GJets_v19_";
                //std::string htitle = "QCD_v19_";
                std::string htitle0 = "GMSB_L"+lambda+"_v19_T22_";
                std::string htitle1 = "GMSB_L"+lambda+"_v19_T11_";
                std::string htitle2 = "GMSB_L"+lambda+"_v19_T33_";
                std::string htitle3 = "GMSB_L"+lambda+"_v19_T44_";

                HistMaker base;

				float nphos = 2;
                std::string isoline = "genSigPerfect_nSigPho2_";
                //std::string isoline = "BG_nSigPho1_";			
	
                std::string outfilenames = outfilename + ofnending;
				std::string outfilenames0 = outfilename0 + isoline + ofnending;
                std::string outfilenames1 = outfilename1 + isoline + ofnending;
                std::string outfilenames2 = outfilename2 + isoline + ofnending;
                std::string outfilenames3 = outfilename3 + isoline + ofnending;

                std::string htitles =  htitle + isoline;
                std::string htitle0s =  htitle0 + isoline;
                std::string htitle1s =  htitle1 + isoline;
                std::string htitle2s =  htitle2 + isoline;
                std::string htitle3s =  htitle3 + isoline;

                base.histMaker( listdir, infilename, outfilenames0, htitle0s, 1, 22, nphos, 0, 0 );
                base.histMaker( listdir, infilename, outfilenames1, htitle1s, 1, 11, nphos, 0, 0 );
                base.histMaker( listdir, infilename, outfilenames2, htitle2s, 1, 33, nphos, 0, 0 );
                base.histMaker( listdir, infilename, outfilenames3, htitle3s, 1, 44, nphos, 0, 0 );

				//base.histMaker( listdir, infilename, outfilenames, htitles, 1, 0, 0, 0, 0 );

				//base.histMaker( listdir, infilename, outfilenames, htitles, fillwt?, gensigtype?, #phos?, n/a, n/a );

    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

