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
    float evtgwt = evtGenWgt;
    //float evtgwt = 1;
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


    if( DEBUG ) std::cout << "Finding gen signal photons" << std::endl;
	///////////////////////    finding gen signal photons    ////////////////////////

	std::vector<uInt> genSigPhoAccIdx;
    std::vector<uInt> genSigPhoIdx;
	uInt nGenPart = genPartPdgId->size();
	for( uInt it = 0; it < nGenPart; it++ ){

		auto partType = (*genPartPdgId)[it];
		if( partType == 22 ){

            if( std::abs((*genPartEta)[it]) > 1.442 ) continue;
            if( (*genPartPt)[it] < 30 ) continue;

			auto phoSusType = (*genPartSusId)[it];
            if( phoSusType == 22 ){ genSigPhoIdx.push_back(it); }

            //if( std::abs((*genPartEta)[it]) > 1.442 ) continue;
            //if( (*genPartPt)[it] < 30 ) continue;

        	bool jetphoton = false;
        	for( int jit = 0; jit < nSelJets; jit++ ){

				bool minPt = (*selJetPt)[jit] < 75.0;
				bool minQual = (*selJetQuality)[jit]  < 2;
				if( minPt || minQual ) continue;

            	float dpjeta = (*selJetEta)[jit] - (*genPartEta)[it];
            	float dpjphi = dPhi( (*selJetPhi)[jit], (*genPartPhi)[it] );
            	float dr = hypo( dpjeta, dpjphi );
            	if( dr < 0.4 ) jetphoton = true;

        	} // for( int jit = 0; jit < nSelJets; jit++ )
        	if( jetphoton ) continue; //{ std::cout << "Jet-Photon Overlap!!!" << std::endl; continue; }

			if( phoSusType == 22 ){ genSigPhoAccIdx.push_back(it); }

		}//<<>>if( partType == 22 )

	}//<<>>for( int it = 0; it < nGenPart; it++ )


/*
	//////////////////////////////  gen matching exclusive   ///////////////////////////////////////////////////

	std::vector<int> genPhoPartIndx;
    std::vector<float> genPhoMDR;
    std::vector<int> genPhoRecoPhoIndx;
    //std::vector<int> prevPhoIndx;
	float drthres = 0.4;
	for( int idx = 0; idx < nGenPart; idx++ ){ 
	 	if( (*genPartPdgId)[idx] == 22 ){ 
			genPhoPartIndx.push_back(idx); 
			genPhoMDR.push_back(drthres); 
			genPhoRecoPhoIndx.push_back(-9); 
			//prevPhoIndx.push_back(-9);
		}//<<>>if( (*genPartPdgId)[idx] == 22 )
	}//<<>>for( int idx = 0; idx < nGenPart; idx++ )
    std::vector<int> recoPhoPartIndx;
	for( int it = 0; it < nSelPhotons; it++ ){ recoPhoPartIndx.push_back(-1); }
	bool done = false;
	int phoadvance = 0;
	while( not done ){
		
		for( int it = phoadvance; it < nSelPhotons; it++ ){
	    	for( int pgidx = 0; pgidx < genPhoRecoPhoIndx.size(); pgidx++ ){
				int partindx = genPhoPartIndx[pgidx];	
	            float dpgeta = (*genPartEta)[partindx] - (*selPhoEta)[it];
	            float dpgphi = dPhi( (*genPartPhi)[partindx], (*selPhoPhi)[it] );
	            float dr = hypo( dpgeta, dpgphi );
	            if( dr < genPhoMDR[pgidx] ){ genPhoMDR[pgidx] = dr; genPhoRecoPhoIndx[pgidx] = it; } 
			}//<<>>for( int pgidx = 0; pgidx < genSigPhoIdx.size(); pgidx++ )
		}//<<>>for( int it = 0; it < nSelPhotons; it++ )
		bool unassigned = false;
	    for( int pgidx = 0; pgidx < genPhoMDR.size(); pgidx++ ){
			for( int pgidx2 = pgidx+1; pgidx2 < genPhoMDR.size(); pgidx2++ ){
				if( ( genPhoRecoPhoIndx[pgidx] != -9 ) && ( genPhoRecoPhoIndx[pgidx] == genPhoRecoPhoIndx[pgidx2] ) ){
					unassigned = true;
					if( genPhoMDR[pgidx2] >= genPhoMDR[pgidx] ){ genPhoRecoPhoIndx[pgidx2] = -9; genPhoMDR[pgidx2] = drthres; }
					else { genPhoRecoPhoIndx[pgidx] = -9; genPhoMDR[pgidx] = drthres; }
				}//<<<>>if( genPhoIndx[pgidx] == genPhoIndx[pgidx2] )
			}//<<>>for( int pgidx2 = pgidx; pgidx2 < genGenPhoMDR.size(); pgidx2++ )
		}//<<>>for( int pgidx = 0; pgidx < genGenPhoMDR.size(); pgidx++ )
		//for( int iter = 0; iter < genPhoRecoPhoIndx.size(); iter++ ){
		//	if( genPhoRecoPhoIndx[iter] != prevPhoIndx[iter] ){ prevPhoIndx[iter] = genPhoRecoPhoIndx[iter]; }
		//}//<<>>for( int iter = 0; iter < genPhoRecoPhoIndx.size(); iter++ )
		if( unassigned ) phoadvance++;
		else done = true;

	}//<<>>while()
    for( int iter = 0; iter < genPhoRecoPhoIndx.size(); iter++ ){
        if( genPhoRecoPhoIndx[iter] != -9 ) recoPhoPartIndx[genPhoRecoPhoIndx[iter]] = genPhoPartIndx[iter];
    }//<<>>for( int iter = 0; iter < genPhoRecoPhoIndx.size(); iter++ )

*/

    if( DEBUG ) std::cout << "Finding reeco signal photons" << std::endl;
	////////////////////// find reco signal photons  /////////////////////////////////////////////////////////////////////////////////////////

	std::vector<int> recoPhoOrderIndx;
    std::vector<int> recoPhoAccOrderIndx;
    if( DEBUG ) std::cout << " - Looping over " << nSelPhotons << " photons" << std::endl;
    for( int it = 0; it < nSelPhotons; it++ ){

		if( std::abs((*selPhoEta)[it]) > 1.442 ) continue;
        if( (*selPhoPt)[it] < 20 ) continue;

        recoPhoOrderIndx.push_back(it);

		bool jetphoton = false;
		for( int jit = 0; jit < nSelJets; jit++ ){

            bool minPt = (*selJetPt)[jit] < 75.0;
            bool minQual = (*selJetQuality)[jit]  < 2;
            if( minPt || minQual ) continue;

			float dpjeta = (*selJetEta)[jit] - (*selPhoEta)[it];
			float dpjphi = dPhi( (*selJetPhi)[jit], (*selPhoPhi)[it] );	
			float dr = hypo( dpjeta, dpjphi );
			if( dr < 0.4 ) jetphoton = true;

		} // for( int jit = 0; jit < nSelJets; jit++ )
		if( jetphoton ) continue; //{ std::cout << "Jet-Photon Overlap!!!" << std::endl; continue; }

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;        
        bool hcalsum = true;
        bool tsptscdr4 = (*selPhoTrkSumPtHollowConeDR04)[it] < 6.0; //(*selPhoTrkSumPtHollowConeDR04)[it] < cutvalue;
		bool ecalrhsum = (*selPhoEcalRHSumEtConeDR04)[it] < 10.0;
        bool htoem = (*selPhoHadTowOverEM)[it] < 0.02;
        bool isoskip = not ( htoem && tsptscdr4 && ecalrhsum && hcalsum );
        if( isoskip ) continue;

        if( DEBUG ) std::cout << " -- looping photons : Filling recoPhoOrderIndx " << it << std::endl;

        recoPhoAccOrderIndx.push_back(it);

    }//<<>>for( int it = 0; it < nPhotons; it++ )

	/////////////////// find gen-reco signal photon match any ////////////////////////////////////////////////////////////////////////////////

    if( DEBUG ) std::cout << "Finding Susy-Gen Matched for Reco: " << recoPhoOrderIndx.size() << " Gen: " << genSigPhoIdx.size() << std::endl;

    std::vector<int> genPhoSigMatch;
	for( auto genIndx : genSigPhoIdx ){
		int genSusMatch = -1;
		if( DEBUG ) std::cout << " - for genIndx: " << genIndx << std::endl;
		for( auto phoptit : recoPhoOrderIndx ){
			if( DEBUG ) std::cout << " -- looping photons : getting gen index for " << phoptit << std::endl;
			auto partIndx = (*selPhoGenIdx)[phoptit];
            if( DEBUG ) std::cout << " -- looping photons : hass partIndx: " << partIndx << std::endl;
            //auto susyid =  ( partIndx >= 0 ) ? (*genPartSusId)[partIndx] : -1;
            //bool isGMSB = ( susyid == 22 );
            //if( not isGMSB ) continue;
			if( partIndx == genIndx ){ genSusMatch = genIndx; break; }
		}//<<>>for( auto phoptit = recoPhoOrderIndx.crbegin(); phoptit != recoPhoOrderIndx.crend(); phoptit++ )
		if( genSusMatch != -1 ){ genPhoSigMatch.push_back(genSusMatch); }
	}//<<>>for( auto genIndx : genSigPhoIdx )

    if( DEBUG ) std::cout << "Finding reco-gen signal photon match" << std::endl;
    /////////////////// find reco-gen signal photon match ////////////////////////////////////////////////////////////////////////////////

	int phocount = 0;
	std::vector<int> sigPhoIndx;
    std::vector<int> sigPhoOrder;
    std::vector<int> sigPhoGenMatch;
    std::vector<int> sigPhoGenOrder;
	if( recoPhoOrderIndx.size() > 0 ){
		//std::cout << " Pho Pts : ";
		for( auto phoptit : recoPhoOrderIndx ){ 
            phocount++;
            auto partIndx = (*selPhoGenIdx)[phoptit];
			auto susyid =  ( partIndx >= 0 ) ? (*genPartSusId)[partIndx] : -1;
            bool isGMSB = ( susyid == 22 );
            if( DEBUG ) std::cout << " --- looping phoOrder : getGenindx :  "  << partIndx << std::endl;
			//bool isGMSB = (*selPhoSusyId)[index] == 22;
            if( DEBUG ) std::cout << " --- looping phoOrder : find is susy :  "  << (*genPartSusId)[partIndx] << std::endl;
			//std::cout << phoptit->first << " " << (*selPhoPt)[index] << " " << isGMSB << " : ";
			if( not isGMSB ) continue;

			int genSusMatch = -1;
			int genSusOrder = 0;
			int genSusCnt = 0;

			for( auto genIndx : genSigPhoIdx ){ genSusCnt++; if( partIndx == genIndx ){ genSusMatch = genIndx; genSusOrder = phocount; break; }}
			sigPhoIndx.push_back(partIndx);
			sigPhoOrder.push_back(phocount);
			sigPhoGenMatch.push_back(genSusMatch);
            sigPhoGenOrder.push_back(genSusOrder);
		}//<<>>for( auto phoptit = recoPhoOrderIndx.crbegin(); phoptit != recoPhoOrderIndx.crend(); phoptit++ )
		//std::cout << std::endl;
	}//<<>>if( recoPhoOrderIndx.size() > 0 )

    if( DEBUG ) std::cout << "Finding reco-gen signal photon match acc" << std::endl;
    phocount = 0;
    std::vector<int> sigPhoAccIndx;
    std::vector<int> sigPhoAccOrder;
    std::vector<int> sigPhoGenAccMatch;
    std::vector<int> sigPhoGenAccOrder;
    if( recoPhoAccOrderIndx.size() > 0 ){
        for( auto phoptit : recoPhoAccOrderIndx ){
            phocount++;
            auto partIndx = (*selPhoGenIdx)[phoptit];
            auto susyid =  ( partIndx >= 0 ) ? (*genPartSusId)[partIndx] : -1;
            bool isGMSB = ( susyid == 22 );
          
            if( DEBUG ) std::cout << " --- looping phoOrder acc : getGenindx :  "  << partIndx << std::endl;  
            //bool isGMSB = (*selPhoSusyId)[index] == 22;
            if( DEBUG ) std::cout << " --- looping phoOrder acc : find is susy :  "  << (*genPartSusId)[partIndx] << std::endl;
            //std::cout << phoptit->first << " " << (*selPhoPt)[index] << " " << isGMSB << " : ";
            if( isGMSB ){
            
            	int genSusMatch = -1;
            	int genSusOrder = 0;
            	int genSusCnt = 0; 
            	for( auto genIndx : genSigPhoAccIdx ){ genSusCnt++; if( partIndx == genIndx ){genSusMatch = genIndx; genSusOrder = phocount; break;}}
            	sigPhoAccIndx.push_back(partIndx);
            	sigPhoAccOrder.push_back(phocount);
            	sigPhoGenAccMatch.push_back(genSusMatch);
            	sigPhoGenAccOrder.push_back(genSusOrder);

					hist1d[20]->Fill((*selPhoPt)[phoptit]);
                    hist1d[21]->Fill((*selPhoHadTowOverEM)[phoptit]);
                    hist1d[22]->Fill((*selPhoTrkSumPtHollowConeDR04)[phoptit]);
                    hist1d[23]->Fill((*selPhoEcalRHSumEtConeDR04)[phoptit]);
			
			} else { //<<>>if( isGMSB )
				
				//if( phocount < 3 && genSigPhoAccIdx.size() > 0 ){
                if( phocount < 3 ){

					hist1d[25]->Fill((*selPhoPt)[phoptit]);
                    hist1d[26]->Fill((*selPhoHadTowOverEM)[phoptit]);
                    hist1d[27]->Fill((*selPhoTrkSumPtHollowConeDR04)[phoptit]);
                    hist1d[28]->Fill((*selPhoEcalRHSumEtConeDR04)[phoptit]);

				}//if( phocount < 3 )
			}//<<>> else { //<<>>if( isGMSB )        

        }//<<>>for( auto phoptit = recoPhoAccOrderIndx.crbegin(); phoptit != recoPhoAccOrderIndx.crend(); phoptit++ )
    }//<<>>if( recoPhoAccOrderIndx.size() > 0 )

    if( DEBUG ) std::cout << "Filling Hists Acc" << std::endl;
    /////////////////// fill hists Acc ////////////////////////////////////////////////////////////////////////////////

    if( DEBUG ) std::cout << " -- looping photons : Filling gen cnt Acc hists "  << std::endl;
    auto nGenSigPhoAcc = genSigPhoAccIdx.size();
    if( nGenSigPhoAcc == 0 ){ hist1d[1]->Fill(0); hist1d[1]->Fill(1);}
    //if( nGenSigPho == 1 ){ hist1d[1]->Fill(2); hist1d[1]->Fill(3); hist1d[1]->Fill(4); hist1d[1]->Fill(5); hist1d[1]->Fill(6);}
    if( nGenSigPhoAcc == 2 ){ hist1d[1]->Fill(2); hist1d[1]->Fill(5); hist1d[1]->Fill(8); }//<<>>if( nGenSigPho == 2 )
    if( nGenSigPhoAcc == 1 ){ hist1d[1]->Fill(9); hist1d[1]->Fill(11); hist1d[1]->Fill(12);}

    auto nRecoAccPhos = recoPhoAccOrderIndx.size();
    if( DEBUG ) std::cout << " -- looping photons : Filling reco cnt Acc hists "  << std::endl;
	// sigPhoGenAccOrder[x] == 0 -> unmatched   sigPhoGenAccOrder[x] > 0 -> matched
    int nSigPhoAcc = sigPhoAccIndx.size();
	if( nSigPhoAcc == 0 && nGenSigPhoAcc == 0 ){
        hist1d[0]->Fill(0); hist1d[1]->Fill(13);
        if( nRecoAccPhos > 0 ) hist1d[0]->Fill(13);
    }//<<>>if( nSigPho == 0 && nGenSigPho == 0 )
    if( nSigPhoAcc > 0 && nGenSigPhoAcc == 0 ) hist1d[0]->Fill(1);
	if( nSigPhoAcc == 2 && nGenSigPhoAcc == 2 ){
		hist1d[0]->Fill(2); hist1d[1]->Fill(3); hist1d[1]->Fill(4);
		if( sigPhoGenAccMatch[0] > -1 && sigPhoGenAccMatch[1] > -1 ){
			if( sigPhoGenAccOrder[0] == 1 && sigPhoGenAccOrder[1] == 2 ) hist1d[0]->Fill(3);
            if( sigPhoGenAccOrder[0] == 1 && sigPhoGenAccOrder[1] > 2 ) hist1d[0]->Fill(4);
		}//<<>>if( sigPhoGenAccMatch[0] > -1 && sigPhoGenAccMatch[1] > -1 )
	}//<<>>if( nSigPhoAcc == 2 && nGenSigPhoAcc == 2 )
    if( nSigPhoAcc == 1 && nGenSigPhoAcc == 2 ){
        hist1d[0]->Fill(5); hist1d[1]->Fill(6); hist1d[1]->Fill(7); hist1d[1]->Fill(14);
        if( nRecoAccPhos > 1 ) hist1d[0]->Fill(14);
        if( sigPhoGenAccMatch[0] > -1 ){
            if( sigPhoGenAccOrder[0] == 1 ) hist1d[0]->Fill(6);
            if( sigPhoGenAccOrder[0] == 2 ) hist1d[0]->Fill(7);
    	}//<<>>if( sigPhoGenAccMatch[0] > -1 )
    }//<<>>if( nSigPhoAcc == 1 && nGenSigPhoAcc == 2 )
    if( nSigPhoAcc == 0 && nGenSigPhoAcc == 2 ){
        hist1d[0]->Fill(8); hist1d[1]->Fill(15); hist1d[1]->Fill(16);
        if( nRecoAccPhos > 0 ) hist1d[0]->Fill(15);
        if( nRecoAccPhos > 1 ) hist1d[0]->Fill(16);
    }//<<>>if( nSigPhoAcc == 0 && nGenSigPhoAcc == 2 )
    if( nSigPhoAcc == 1 && nGenSigPhoAcc == 1 ){
        hist1d[0]->Fill(9); hist1d[1]->Fill(10); hist1d[1]->Fill(17);
        if( nRecoAccPhos > 1 ) hist1d[0]->Fill(17);
        if( sigPhoGenAccMatch[0] > -1 ){
            if( sigPhoGenAccOrder[0] == 1 ) hist1d[0]->Fill(10);
        }//<<>>if( sigPhoGenAccMatch[0] > -1 )
    }//<<>>if( nSigPhoAcc == 1 && nGenSigPhoAcc == 2 )
    if( nSigPhoAcc == 0 && nGenSigPhoAcc == 1 ){
        hist1d[0]->Fill(11); hist1d[1]->Fill(18);
        if( nRecoAccPhos > 0 ) hist1d[0]->Fill(18);
    }//<<>>if( nSigPho == 0 && nGenSigPho == 1 )
    if( nSigPhoAcc > 1 && nGenSigPhoAcc == 1 ) hist1d[0]->Fill(12);


    /////////////////// fill hists no Acc ////////////////////////////////////////////////////////////////////////////////
    if( DEBUG ) std::cout << " -- looping photons : Filling gen cnt hists "  << std::endl;

    auto nGenSigPho = genSigPhoIdx.size();
    if( nGenSigPho == 0 ){ hist1d[7]->Fill(0); hist1d[7]->Fill(1);}
    //if( nGenSigPho == 1 ){ hist1d[7]->Fill(2); hist1d[7]->Fill(3); hist1d[7]->Fill(4); hist1d[7]->Fill(5); hist1d[7]->Fill(6);}
    if( nGenSigPho == 2 ){ hist1d[7]->Fill(2); hist1d[7]->Fill(5); hist1d[7]->Fill(8); }//<<>>if( nGenSigPho == 2 )
    if( nGenSigPho == 1 ){ hist1d[7]->Fill(9); hist1d[7]->Fill(11); hist1d[7]->Fill(12);}

	auto nRecoPhos = recoPhoOrderIndx.size();
    if( DEBUG ) std::cout << " -- looping photons : Filling reco cnt hists "  << std::endl;
    // sigPhoGenOrder[x] == 0 -> unmatched   sigPhoGenOrder[x] > 0 -> matched
    int nSigPho = sigPhoIndx.size();
    if( nSigPho == 0 && nGenSigPho == 0 ){ 
		hist1d[6]->Fill(0); hist1d[7]->Fill(13); 
		if( nRecoPhos > 0 ) hist1d[6]->Fill(13);
	}//<<>>if( nSigPho == 0 && nGenSigPho == 0 )
    if( nSigPho > 0 && nGenSigPho == 0 ) hist1d[6]->Fill(1);
    if( nSigPho == 2 && nGenSigPho == 2 ){
        hist1d[6]->Fill(2); hist1d[7]->Fill(3); hist1d[7]->Fill(4);
        if( sigPhoGenMatch[0] > -1 && sigPhoGenMatch[1] > -1 ){
            if( sigPhoGenOrder[0] == 1 && sigPhoGenOrder[1] == 2 ) hist1d[6]->Fill(3);
            if( sigPhoGenOrder[0] == 1 && sigPhoGenOrder[1] > 2 ) hist1d[6]->Fill(4);
        }//<<>>if( sigPhoGenMatch[0] > -1 && sigPhoGenMatch[1] > -1 )
    }//<<>>if( nSigPho == 2 && nGenSigPho == 2 )
    if( nSigPho == 1 && nGenSigPho == 2 ){
        hist1d[6]->Fill(5); hist1d[7]->Fill(6); hist1d[7]->Fill(7); hist1d[7]->Fill(14);
		if( nRecoPhos > 1 ) hist1d[6]->Fill(14);
        if( sigPhoGenMatch[0] > -1 ){
            if( sigPhoGenOrder[0] == 1 ) hist1d[6]->Fill(6);
            if( sigPhoGenOrder[0] == 2 ) hist1d[6]->Fill(7);
        }//<<>>if( sigPhoGenMatch[0] > -1 )
    }//<<>>if( nSigPho == 1 && nGenSigPho == 2 )
    if( nSigPho == 0 && nGenSigPho == 2 ){ 
		hist1d[6]->Fill(8); hist1d[7]->Fill(15); hist1d[7]->Fill(16); 
		if( nRecoPhos > 0 ) hist1d[6]->Fill(15);
		if( nRecoPhos > 1 ) hist1d[6]->Fill(16);
	}//<<>>if( nSigPho == 0 && nGenSigPho == 2 )
    if( nSigPho == 1 && nGenSigPho == 1 ){
        hist1d[6]->Fill(9); hist1d[7]->Fill(10); hist1d[7]->Fill(17);
		if( nRecoPhos > 1 ) hist1d[6]->Fill(17);
        if( sigPhoGenMatch[0] > -1 ){
            if( sigPhoGenOrder[0] == 1 ) hist1d[6]->Fill(10);
        }//<<>>if( sigPhoGenMatch[0] > -1 )
    }//<<>>if( nSigPho == 1 && nGenSigPho == 2 )
    if( nSigPho == 0 && nGenSigPho == 1 ){ 
		hist1d[6]->Fill(11); hist1d[7]->Fill(18); 
		if( nRecoPhos > 0 ) hist1d[6]->Fill(18);
	}//<<>>if( nSigPho == 0 && nGenSigPho == 1 )
    if( nSigPho > 1 && nGenSigPho == 1 ) hist1d[6]->Fill(12);

    /////////////////// fill hists other ////////////////////////////////////////////////////////////////////////////////

	if( DEBUG ) std::cout << " -- looping photons : Filling reco cnt other hists "  << std::endl;

	double nGenSigMatch = genPhoSigMatch.size();

	//hist1d[3]->Fill(0, argggghhhhh!!!!  no worky
    hist1d[3]->Fill(2, nGenSigPhoAcc);		hist1d[4]->Fill(2, nGenSigPho); //	hist1d[3]->Fill(1, nGenSigPhoAcc);

    hist1d[3]->Fill(3, nSigPho);			hist1d[4]->Fill(3, nGenSigPho); //		hist1d[3]->Fill(2, nGenSigPho);
    hist1d[3]->Fill(4, nSigPhoAcc);			hist1d[4]->Fill(4, nSigPho); //	hist1d[3]->Fill(3, nSigPhoAcc);
    hist1d[3]->Fill(5, nSigPhoAcc);			hist1d[4]->Fill(5, nGenSigPhoAcc); //    hist1d[3]->Fill(4, nSigPho);

    hist1d[3]->Fill(1, nGenSigMatch);       hist1d[4]->Fill(1, nGenSigPho);//   hist1d[3]->Fill(0, nGenPhoSigMatch);

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

	hist1d[2]->Divide(hist1d[0],hist1d[1]);
    hist1d[5]->Divide(hist1d[3],hist1d[4]);
    hist1d[8]->Divide(hist1d[6],hist1d[7]);

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------
	int bins = 20;
    hist1d[0] = new TH1D("phoSigMatchCntAcc", "phoSigMatchCntAcc", bins, 0, bins); 
    hist1d[1] = new TH1D("phoGenMatchCntAcc", "phoGenMatchCntAcc", bins, 0, bins);
    hist1d[2] = new TH1D("phoSigGenMatchEffAcc", "phoSigGenMatchEffAcc", bins, 0, bins);
    hist1d[3] = new TH1D("phoSigMatchCntOther", "phoSigMatchCntOther", bins, 0, bins);
    hist1d[4] = new TH1D("phoGenMatchCntOther", "phoGenMatchCntOther", bins, 0, bins);
    hist1d[5] = new TH1D("phoSigGenMatchEffOther", "phoSigGenMatchEffOther", bins, 0, bins);
    hist1d[6] = new TH1D("phoSigMatchCnt", "phoSigMatchCnt", bins, 0, bins);
    hist1d[7] = new TH1D("phoGenMatchCnt", "phoGenMatchCnt", bins, 0, bins);
    hist1d[8] = new TH1D("phoSigGenMatchEff", "phoSigGenMatchEff", bins, 0, bins);

	hist1d[20] = new TH1D("sigPhoPt", "sigPhoPt", 250, 0, 250);
    hist1d[21] = new TH1D("sigPhoHoe", "sigPhoHoe", 250, 0, 0.25);
    hist1d[22] = new TH1D("sigPhoTrk", "sigPhoTrk", 200, 0, 20);
    hist1d[23] = new TH1D("sigPhoEcal", "sigPhoEcal", 200, 0, 20);

    hist1d[25] = new TH1D("fakePhoPt", "fakePhoPt", 250, 0, 250);
    hist1d[26] = new TH1D("fakePhoHoe", "fakePhoHoe", 250, 0, 0.25);
    hist1d[27] = new TH1D("fakePhoTrk", "fakePhoTrk", 200, 0, 20);
    hist1d[28] = new TH1D("fakePhoEcal", "fakePhoEcal", 200, 0, 20);


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

                auto outfilename = "KUCMS_phoRecoGenSigEff_Hists_t3.root"; //9

				auto htitle = "KUCMS_phoRecoGenSigEff_Hists_";

                HistMaker base;
                base.histMaker( listdir, infilename, outfilename, htitle );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

