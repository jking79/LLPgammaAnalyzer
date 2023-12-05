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

void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");

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
    float scale = 1;
    std::string configKey(*DataSetKey);
    float xsec = (configInfo[configKey])["sCrossSection"];
    float segwt = (configInfo[configKey])["nEvents"];
    auto fillwt = scale * ( xsec * 1000 ) * ( evtgwt / segwt );
    //std::cout << " evtfillwt :( " << xsec << " * 1000 ) * ( " << evtgwt << " / " << segwt << " ) = " << fillwt << std::endl;
    //auto fillwt = (configInfo[configKey])["sCrossSection"] * evtwt;
	// JetHT2018D
	// GMSBL100 - GMSBL400
    //if( configKey != "GMSBL100" ) return;

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

	int totNSelPhotons = 0;
    int totNSel30Photons = 0;
    int totNSel100Photons = 0;
	int totNCutPhotons = 0;
    int totNCut30Photons = 0;
    int totNCut100Photons = 0;
	std::map<float,int> phoOrderIndx;
	//vector<int> phoOrderIndx;
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

        bool isNino = phoSusId == 22;
        bool isXinoWZ = phoSusId == 23 || phoSusId == 24 || phoSusId == 33 || phoSusId == 34;
        bool isXfsr = phoSusId == 36 || phoSusId == 32 || phoSusId == 35;
        bool isPrompt = phoSusId == 29;
        bool isSfsr = phoSusId == 37 || phoSusId == 20 || phoSusId == 21 || phoSusId == 30 || phoSusId == 31;

		bool isSigType = isNino;
		bool isXfsrType = isXinoWZ || isXfsr;
		bool isOtherType = isPrompt || isSfsr; // || phoSusId < 1;
        bool isNotSusyType = isPrompt || phoSusId < 9;

		//if( isSigType ) continue;
		//if( isNotSusyType ) continue;

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;        

		totNSelPhotons++;
		if( (*selPhoPt)[it] > 30 ) totNSel30Photons++;
        if( (*selPhoPt)[it] > 100 ) totNSel100Photons++;

        bool hcalsum = true;
        bool tsptscdr4 = (*selPhoTrkSumPtSolidConeDR04)[it] < 6.0; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
		bool ecalrhsum = (*selPhoEcalRHSumEtConeDR04)[it] < 10.0;
        bool htoem = (*selPhoHadTowOverEM)[it] < 0.02;
        bool isoskip = not ( htoem && tsptscdr4 && ecalrhsum && hcalsum );
        if( isoskip ) continue;
        totNCutPhotons++;
		if( phoOrderIndx.count((*selPhoPt)[it]) > 0 ) std::cout << "Duplicate Pho Pt ---------------------- " << std::endl;
        else phoOrderIndx[(*selPhoPt)[it]] = it;


    }//<<>>for( int it = 0; it < nPhotons; it++ )
	bool pho1sig(false), pho2sig(false), pho12sig(false), pho3sig(false);
	int phocount = 0;
	if( phoOrderIndx.size() > 1 ){
		//std::cout << " Pho Pts : ";
		for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ ){ 
			auto index = phoptit->second;
			bool isGMSB = (*selPhoSusyId)[index] == 22;
			//std::cout << phoptit->first << " " << (*selPhoPt)[index] << " " << isGMSB << " : ";
			phocount++;
			if( phocount == 1 && isGMSB ){ pho1sig = true; pho12sig = true; }
			if( phocount == 2 && isGMSB ){ pho2sig = true; } 
			if( phocount == 2 && not isGMSB ){ pho12sig = pho12sig && false; }
			if( phocount > 2 && isGMSB ){ pho3sig = true; }
		}
		//std::cout << std::endl;
	}
	if( not( pho1sig || pho2sig || pho12sig || pho3sig ) ) hist1d[0]->Fill(0);
	if( pho12sig )  hist1d[0]->Fill(1);
    if( pho1sig )  hist1d[0]->Fill(2);
    if( pho2sig )  hist1d[0]->Fill(3);
    if( pho3sig )  hist1d[0]->Fill(4);

    preCutNPhotons += totNSelPhotons;
    preCut30NPhotons += totNSel30Photons;
    preCut100NPhotons += totNSel100Photons;
    postCutNPhotons += totNCutPhotons;
    postCut30NPhotons += totNCut30Photons;
    postCut100NPhotons += totNCut100Photons;

	//-------- electrons --------------------------------------


    //-------- jets  --------------------------------------
    //UInt_t          nJets;
    //UInt_t          nSelJets;
    //std::vector<float>   *selJetEnergy;
    //std::vector<float>   *selJetEta;
    //std::vector<float>   *selJetMass;
    //std::vector<float>   *selJetPhi;
    //std::vector<float>   *selJetPt;
    //std::vector<int>     *selJetQuality;
    //std::vector<float>   *selJetTime;


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

    hist1d[0] = new TH1D("phoSigMatchEff", "phoSigMatchEEff", 5, 0, 5); 

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

    //------  electrons 350 - 400

	//-------- event vars 400 - 450
	
    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------


    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

	//------------------------------------------------------------------------------------

}//<<>>void HistMaker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "skims_files/";
				auto infilename = "KUCMS_Master_Skim_List.txt";

                auto outfilename = "KUCMS_This_Skim_BaseHists.root"; //9

				auto htitle = "KUCMS_This_BaseHists_";

                HistMaker base;
                base.histMaker( listdir, infilename, outfilename, htitle );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

