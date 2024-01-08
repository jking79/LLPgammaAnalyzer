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
// --   profile function ----------------
// ------------------------------------------------------------------------------------------------------------------------------------------

void profileThisTH2D(TH2D* nhist, TH1D* prof, TH1D* fithist, float range ){

    std::cout << "Profile " << " hist: " << nhist->GetName() << std::endl;

    const auto nXBins = nhist->GetNbinsX();
    TH1D* prohists[nXBins];
    for (auto ibinX = 1; ibinX <= nXBins; ibinX++){

		std::string temp_title = "profile_bin" + std::to_string( ibinX );
        prohists[ibinX-1] = (TH1D*)nhist->ProjectionY(temp_title.c_str(),ibinX,ibinX);
        auto mean = prohists[ibinX-1]->GetMean();
        auto stdv = prohists[ibinX-1]->GetStdDev();
        auto norm = prohists[ibinX-1]->GetBinContent(prohists[ibinX-1]->GetMaximumBin());
        auto high = mean + range*stdv;
        auto low = mean - range*stdv;
        if( abs(stdv) > 0.0 && abs(norm) > 1 ){
            auto tmp_form = new TFormula("tmp_formula","[0]*exp(-0.5*((x-[1])/[2])**2)");
            auto tmp_fit  = new TF1("tmp_fit",tmp_form->GetName(),low,high);
            tmp_fit->SetParameter(0,norm); //tmp_fit->SetParLimits(0,norm/2,norm*2);
            tmp_fit->SetParameter(1,mean); //tmp_fit->SetParLimits(1,-2,2);
            tmp_fit->SetParameter(2,stdv); //tmp_fit->SetParLimits(2,0,10);
            prohists[ibinX-1]->Fit(tmp_fit->GetName(),"RBQ0");
            auto fmean = tmp_fit->GetParameter(1);
            auto error = tmp_fit->GetParError(1);
            auto fNdf = tmp_fit->GetNDF();
            auto fProb = tmp_fit->GetProb();
            auto fChi2 = tmp_fit->GetChisquare();
            // set new contents
            //if( fNdf > 0 && fProb > 0.001 && error < fmean ){
            if( fNdf > 0 && error < stdv ){
                auto fChi2Ndf = fChi2/fNdf;
                fithist->SetBinContent( ibinX, fChi2Ndf );
                fithist->SetBinError( ibinX, 0 );
                prof->SetBinContent( ibinX, fmean );
                prof->SetBinError( ibinX, error );
            }//<<>>if( fmean < 1 && error < 0.1 )
            //delete tmp_form;
            delete tmp_fit;
        }//<<>>if( stdv > 0.01 )

    }//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)
	//for( int it = 0; it < nXBins; it++ ){ if(prohists[it]) prohists[it]->Write(); }

}//<<>>void profileTH2D(TH2D* hist, TH1D* prof)

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
	bool isGMSB = (configKey.find("GMSB") != std::string::npos) ? 1 : 0; 

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------

	int totNSelPhotons = 0;
    int totNSel30Photons = 0;
    int totNSel100Photons = 0;
	int totNCutPhotons = 0;
    int totNCut30Photons = 0;
    int totNCut100Photons = 0;
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

	///////////  pho selection ////////////////////////////////////////////////////////////////////
		if( DEBUG ) std::cout << " - Starting Pho Loop" << std::endl;

        //if( (*selPhoPt)[it] < 30 ) continue;
        if( std::abs((*selPhoEta)[it]) > 1.442 ) continue;

		bool jetphoton = false;
		for( int jit = 0; jit < nSelJets; jit++ ){

            bool minPt = (*selJetPt)[jit] < 75.0;
            bool minQual = (*selJetQuality)[jit]  < 2;
            if( minPt || minQual ) continue;

			float dpjeta = (*selJetEta)[jit] - (*selPhoPhi)[it];
			float dpjphi = dPhi( (*selJetPhi)[jit], (*selPhoPhi)[it] );	
			float dr = hypo( dpjeta, dpjphi );
			if( dr < 0.4 ) jetphoton = true;

		} // for( int jit = 0; jit < nSelJets; jit++ )
		if( jetphoton ) continue; //{ std::cout << "Jet-Photon Overlap!!!" << std::endl; continue; }

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;
        bool hcalsum = true;
        bool tsptscdr4 = (*selPhoTrkSumPtHollowConeDR04)[it] < 6.0; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
        bool ecalrhsum = (*selPhoEcalRHSumEtConeDR04)[it] < 10.0;
        bool htoem = (*selPhoHadTowOverEM)[it] < 0.02;
        bool isoskip = not ( htoem && tsptscdr4 && ecalrhsum && hcalsum );
        if( isoskip ) continue;

	///////////  pho selection ////////////////////////////////////////////////////////////////////


        if( DEBUG ) std::cout << " -- looping photons : Start getting dr dp info" << std::endl;
		auto phoClass = (*selPhoQuality)[it];
        //if( phoClass < 1 ) continue;
        auto phoSusId = (*selPhoSusyId)[it]; // selPhoSusyId->at(it);
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

		//if( not isSigType ) continue; //  select only SigType
		//if( not isNotSusyType ) continue; // select only NotSigType

	////////////   pho type selection /////////////////////////////////////////////////////////////

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;        

		totNSelPhotons++;

		if( DEBUG ) std::cout << " -- Filling gen/reco info 2ds" << std::endl;

        auto rvx = (*selPhoSCx)[it];
        auto rvy = (*selPhoSCy)[it];
        auto rvz = (*selPhoSCz)[it];
        auto rpt = (*selPhoPt)[it];
        auto re = (*selPhoEnergy)[it];
        auto rtime = (*selPhoTime)[it];
        auto peta = (*selPhoEta)[it];


        float cms000 = hypo(rvx,rvy,rvz);
        float calcor = cms000/SOL;
		float pvtof = hypo(rvx-PVx,rvy-PVy,rvz-PVz)/SOL;
		rtime = rtime + calcor - pvtof;

		if( rtime < 0.04 ) continue;

        float rt2 = sq2(rtime);
        float s1 = 0.99736-0.0035*rtime-0.00319*rt2;
        float s2 = 1.205-0.079*rtime+0.00093*rt2;//<100, <3 // 
		float s3 = std::exp( -1*sq2(rtime-0.2037)/144.942338);
        //float core = ( rtime < -0.549 ) ? 1 : ( rtime < 3.37 ) ? s1 : s2;
        float core = ( rtime < -0.549 ) ? 1 : s1;
		if( re < 325 && rtime > 7.0 ) core = s2;
        if( re < 200 && rtime > 6.0 ) core = s2;
        if( re < 150 && rtime > 5.0 ) core = s2;
        if( re < 100 && rtime > 4.0 ) core = s2;
		if( rtime > 16.0 ) core = s3;
        ////float core = ( rtime < 0.204 ) ? 1 : std::exp(-1*sq2(rtime-0.204)/144.942338);
        float ce = re/core;
        auto cpt = rpt/core;

		if( cpt < 30 ) continue;

		if( (*selPhoGenIdx)[it] > -1 ){

			if( DEBUG ) std::cout << " -- Filling gen/reco info 2ds in loop" << std::endl;
			auto gidx = (*selPhoGenIdx)[it];
			if( DEBUG ) std::cout << " --- gidx: " << gidx << std::endl;
			auto geta = (*genPartEta)[gidx];
			auto gphi = (*genPartPhi)[gidx];
			auto ge = (*genPartEnergy)[gidx];
            auto gpt = (*genPartPt)[gidx];
			auto gvx = (*genVx)[gidx];
            auto gvy = (*genVy)[gidx];
            auto gvz = (*genVz)[gidx];

			//if( ge < 20 ) continue;
	

            if( DEBUG ) std::cout << " -- doing reco calcs" << std::endl;
			float ceta = std::asinh((rvz-gvz)/hypo(rvx-gvx,rvy-gvy));
			float cphi = std::atan2(rvy-gvy,rvx-gvx);
			auto dr = dR1(geta,gphi,ceta,cphi);
            auto rege = re/ge;
            auto cege = ce/ge;

            float tofPVtoRH = hypo(rvx-PVx,rvy-PVy,rvz-PVz);
            float tmeasured = rtime*SOL+tofPVtoRH;

            float m1 = std::sqrt(sq2(tofPVtoRH)+8*sq2(tmeasured));
            float m2 = (tofPVtoRH+m1)*(tofPVtoRH/(4*sq2(tmeasured)));
            float m3 = tofPVtoRH/tmeasured;
			float m2_phys = ( m2 > 1 ) ? 1 : m2;
            float m3_phys = ( m3 > 1 ) ? 1 : m3;
            float MBetaEqual = 2*ce*std::sqrt(1-sq2(m2_phys));
            //float MBetaEqual = 2*ce*m2_phys;
            float MBetaPrompt = 2*ce*std::sqrt(1-sq2(m3_phys));
            //float MBetaPrompt = 2*ce*m3_phys;
			float bgbeta = gpt/ge;

			hist2d[66]->Fill(ce,ge,fillwt);

			if( not isGMSB ){

                hist1d[63]->Fill(MBetaEqual,fillwt);
                hist1d[64]->Fill(MBetaPrompt,fillwt);
                hist1d[67]->Fill(ce,fillwt);
                hist1d[68]->Fill(cpt,fillwt);

                hist2d[62]->Fill(MBetaEqual,rtime);
                hist2d[63]->Fill(MBetaPrompt,rtime);
                hist2d[65]->Fill(MBetaEqual,MBetaPrompt);
				hist2d[68]->Fill(ce,ge);

				hist2d[12]->Fill(rege,rtime);

                hist2d[70]->Fill(ce,rtime);
                hist2d[72]->Fill(cpt,rtime);

                hist1d[74]->Fill(bgbeta,fillwt);
                hist1d[75]->Fill(m2,fillwt);
                hist1d[76]->Fill(m3,fillwt);
                hist1d[77]->Fill(tmeasured,fillwt);
                hist1d[78]->Fill(m1,fillwt);
                hist1d[80]->Fill(rtime,fillwt);


			}//<<>>if( (*selPhoSusyId)[it] != 22 )


			if( DEBUG ) std::cout << " -- Filling reco 2d hists" << std::endl;

			if( isGMSB ){

				hist2d[0]->Fill(dr,rtime);
            	hist2d[2]->Fill(rege,rtime);
            	hist2d[13]->Fill(cege,rtime);
            	hist2d[4]->Fill(peta,rtime);
				//float space = ce/25.0;
				//for( int it = 0; it < 40; it++ ){ if( space < it && space > (it-1) ) hist2d[100+it]->Fill(rege,rtime);}
			}//<<>>if( isGMSB )

			if( isSigType && isGMSB ){

                //std::cout << " -- Filling gen 2d hists susID: " << (*selPhoSusyId)[it] << std::endl;
				auto nvx = (*selPhoGenSigMomVx)[it];
                auto nvy = (*selPhoGenSigMomVy)[it];
                auto nvz = (*selPhoGenSigMomVz)[it];
            	auto npx = (*selPhoGenSigMomPx)[it];
            	auto npy = (*selPhoGenSigMomPy)[it];
            	auto npz = (*selPhoGenSigMomPz)[it];
				auto ne = (*selPhoGenSigMomEnergy)[it];
                //std::cout << " -- Filling gen 2d hists nvx: " << nvx << " nvy: " << nvy << " nvz: " << nvz << std::endl;
				float np = hypo( npx, npy, npz );
				//float gp = (*selPhoGenSigMomPt)[it];
				float beta = np/ne;
                MBetaEqual = 2*ce*std::sqrt(1-sq2(beta));
                MBetaPrompt = 2*ce*std::sqrt(1-sq2(beta));
				//std::cout << " -- Filling gen 2d hists smpt: " << (*selPhoGenSigMomPt)[it] << " sme: " << (*selPhoGenSigMomEnergy)[it] << std::endl;
				float pl1 = hypo( gvx-nvx, gvy-nvy, gvz-nvz ); 
				float pl2 = hypo( rvx-gvx, rvy-gvy, rvz-gvz );
				//float tofPVtoRH = hypo( rvx-nvx, rvy-nvy, rvz-nvz );
				//std::cout << " -- Filling gen 2d hists gp: " << gp << " ge: " << ge << " SOL: " << SOL << " beta: " << beta << std::endl;
				float t1 = pl1/(beta*SOL);
				float t2 = pl2/SOL;
				float gtime = t1+t2-calcor;
                //float tmeasured = (t1+t2)*SOL;
				float rlalb = (pl1-pl2)/(pl1+pl2);

				//std::cout << " -- Filling gen 2d hists gtime: " << gtime << std::endl;
				if( DEBUG ) std::cout << " -- Filling gen 2d hists" << std::endl;

                hist1d[60]->Fill(rlalb,fillwt);
                hist1d[61]->Fill(MBetaEqual,fillwt);
                hist1d[62]->Fill(MBetaPrompt,fillwt);
                hist1d[65]->Fill(ce,fillwt);
                hist1d[66]->Fill(rpt,fillwt);

                hist1d[69]->Fill(beta,fillwt);
                hist1d[70]->Fill(m2,fillwt);
                hist1d[71]->Fill(m3,fillwt);
                hist1d[72]->Fill(tmeasured,fillwt);
                hist1d[73]->Fill(m1,fillwt);				
                hist1d[79]->Fill(rtime,fillwt);

				hist2d[60]->Fill(MBetaEqual,rtime);
				hist2d[61]->Fill(MBetaPrompt,rtime);

				hist2d[64]->Fill(MBetaEqual,MBetaPrompt);
				hist2d[67]->Fill(ce,ge);

				hist1d[50]->Fill(rtime-gtime,fillwt);
				hist1d[0]->Fill(gtime,fillwt);
           	 	hist2d[1]->Fill(rtime,gtime,fillwt);
            	hist2d[6]->Fill(dr,gtime);
            	hist2d[3]->Fill(rege,gtime);
            	hist2d[5]->Fill(peta,gtime);
            	hist2d[7]->Fill(dr,rtime);
                hist2d[8]->Fill(dr,rege);
            	hist2d[9]->Fill(rege,rtime);
            	hist2d[10]->Fill(peta,rtime);
                hist2d[11]->Fill(rtime,rege);

                hist2d[69]->Fill(ce,rtime);
                hist2d[71]->Fill(cpt,rtime);
                hist2d[73]->Fill(ce,MBetaEqual);
                hist2d[74]->Fill(cpt,MBetaEqual);

			}//<<>>if( (*selPhoSusyId)[it] != 22 ) 

		}//<<>>if( (*selPhoGenIdx) > -1 )

		//if( isGMSB ) continue;
        if( not isGMSB ) continue;

		if( DEBUG ) std::cout << " -- Filling pt > 20 sv hists" << std::endl;		

		hist1d[1]->Fill((*selPhoSigmaIEtaIEta)[it],fillwt);
        hist1d[2]->Fill((*selPhoClstrRn)[it],fillwt);
        hist1d[3]->Fill((*selPhoCovEtaEta)[it],fillwt);
        hist1d[4]->Fill((*selPhoCovEtaPhi)[it],fillwt);
        hist1d[5]->Fill((*selPhoCovPhiPhi)[it],fillwt);
        hist1d[6]->Fill((*selPhoNrh)[it],fillwt);
        hist1d[7]->Fill((*selPhoEtaWidth)[it],fillwt);
        hist1d[8]->Fill((*selPhoPhiWidth)[it],fillwt);
        hist1d[9]->Fill((*selPhoR9)[it],fillwt);
        hist1d[10]->Fill((*selPhoS4)[it],fillwt);
        hist1d[11]->Fill((*selPhoSAlp)[it],fillwt);
        hist1d[12]->Fill((*selPhoSMaj)[it],fillwt);
        hist1d[13]->Fill((*selPhoSMin)[it],fillwt);
        hist1d[14]->Fill((*selPhoSieie)[it],fillwt);
        hist1d[15]->Fill((*selPhoSieip)[it],fillwt);
        hist1d[16]->Fill((*selPhoSipip)[it],fillwt);

		hist2d[50]->Fill((*selPhoCovEtaEta)[it],sq2((*selPhoSieie)[it]));
        hist2d[51]->Fill((*selPhoSAlp)[it],(*selPhoSieie)[it]);

		if( (*selPhoPt)[it] < 30 ) continue;
		totNSel30Photons++;

        if( DEBUG ) std::cout << " -- Filling pt > 30 sv hists" << std::endl;

        hist1d[17]->Fill((*selPhoSigmaIEtaIEta)[it],fillwt);
        hist1d[18]->Fill((*selPhoClstrRn)[it],fillwt);
        hist1d[19]->Fill((*selPhoCovEtaEta)[it],fillwt);
        hist1d[20]->Fill((*selPhoCovEtaPhi)[it],fillwt);
        hist1d[21]->Fill((*selPhoCovPhiPhi)[it],fillwt);
        hist1d[22]->Fill((*selPhoNrh)[it],fillwt);
        hist1d[23]->Fill((*selPhoEtaWidth)[it],fillwt);
        hist1d[24]->Fill((*selPhoPhiWidth)[it],fillwt);
        hist1d[25]->Fill((*selPhoR9)[it],fillwt);
        hist1d[26]->Fill((*selPhoS4)[it],fillwt);
        hist1d[27]->Fill((*selPhoSAlp)[it],fillwt);
        hist1d[28]->Fill((*selPhoSMaj)[it],fillwt);
        hist1d[29]->Fill((*selPhoSMin)[it],fillwt);
        hist1d[30]->Fill((*selPhoSieie)[it],fillwt);
        hist1d[31]->Fill((*selPhoSieip)[it],fillwt);
        hist1d[32]->Fill((*selPhoSipip)[it],fillwt);

        if( (*selPhoPt)[it] < 100 ) continue;
		totNSel100Photons++;

        if( DEBUG ) std::cout << " -- Filling pt > 100 sv hists" << std::endl;

        hist1d[33]->Fill((*selPhoSigmaIEtaIEta)[it],fillwt);
        hist1d[34]->Fill((*selPhoClstrRn)[it],fillwt);
        hist1d[35]->Fill((*selPhoCovEtaEta)[it],fillwt);
        hist1d[36]->Fill((*selPhoCovEtaPhi)[it],fillwt);
        hist1d[37]->Fill((*selPhoCovPhiPhi)[it],fillwt);
        hist1d[38]->Fill((*selPhoNrh)[it],fillwt);
        hist1d[39]->Fill((*selPhoEtaWidth)[it],fillwt);
        hist1d[40]->Fill((*selPhoPhiWidth)[it],fillwt);
        hist1d[41]->Fill((*selPhoR9)[it],fillwt);
        hist1d[42]->Fill((*selPhoS4)[it],fillwt);
        hist1d[43]->Fill((*selPhoSAlp)[it],fillwt);
        hist1d[44]->Fill((*selPhoSMaj)[it],fillwt);
        hist1d[45]->Fill((*selPhoSMin)[it],fillwt);
        hist1d[46]->Fill((*selPhoSieie)[it],fillwt);
        hist1d[47]->Fill((*selPhoSieip)[it],fillwt);
        hist1d[48]->Fill((*selPhoSipip)[it],fillwt);


    }//<<>>for( int it = 0; it < nPhotons; it++ )


    preCutNPhotons += totNSelPhotons;
    preCut30NPhotons += totNSel30Photons;
    preCut100NPhotons += totNSel100Photons;
    postCutNPhotons += totNCutPhotons;
    postCut30NPhotons += totNCut30Photons;
    postCut100NPhotons += totNCut100Photons;

}//<<>>void HistMaker::eventLoop( Long64_t entry )

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

void HistMaker::endJobs(){

	profileThisTH2D( hist2d[11], hist1d[51], hist1d[52], 1.0 );
	//normTH1D(hist1d[52]);
	hist1d[400] = (TH1D*)hist1d[63]->Clone("ratio_mllp_b/s");
	hist1d[400]->Divide(hist1d[61]);
    hist1d[401] = (TH1D*)hist1d[67]->Clone("ratio_ce_b/s");
    hist1d[401]->Divide(hist1d[65]);
    hist1d[402] = (TH1D*)hist1d[68]->Clone("ratio_cpt_b/s");
    hist1d[402]->Divide(hist1d[66]);

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht ){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

    hist1d[0] = new TH1D("phoSigGenTime", "phoSigGenTime", 100, -5, 20); 
	//hist1d[0] = new TH1D("", "", , 0, );
    hist1d[1] = new TH1D("selPhoSigmaIEtaIEta_pt20","selPhoSigmaIEtaIEta_pt20",50, 0, 0.025 );
    hist1d[2] = new TH1D("selPhoClstrRn_pt20","selPhoClstrRn_pt20",100, 0, 1 );
    hist1d[3] = new TH1D("selPhoCovEtaEta_pt20","selPhoCovEtaEta_pt20",100, 0, 0.001 );
    hist1d[4] = new TH1D("selPhoCovEtaPhi_pt20","selPhoCovEtaPhi_pt20",200, -0.001, 0.001 );
    hist1d[5] = new TH1D("selPhoCovPhiPhi_pt20","selPhoCovPhiPhi_pt20",100, 0, 0.001 );
    hist1d[6] = new TH1D("selPhoNrh_pt20","selPhoNrh_pt20",150, 0, 150 );
    hist1d[7] = new TH1D("selPhoEtaWidth_pt20","selPhoEtaWidth_pt20",50, 0, 0.05 );
    hist1d[8] = new TH1D("selPhoPhiWidth_pt20","selPhoPhiWidth_pt20",250, 0, 0.25);
    hist1d[9] = new TH1D("selPhoR9_pt20","selPhoR9_pt20",100, 0, 1 );
    hist1d[10] = new TH1D("selPhoS4_pt20","selPhoS4_pt20",100, 0, 1 );
    hist1d[11] = new TH1D("selPhoSAlp_pt20","selPhoSAlp_pt20",60, -1.5, 1.5 );
    hist1d[12] = new TH1D("selPhoSMaj_pt20","selPhoSMaj_pt20",100, 0, 5 );
    hist1d[13] = new TH1D("selPhoSMin_pt20","selPhoSMin_pt20",50, 0, 2.5 );
    hist1d[14] = new TH1D("selPhoSieie_pt20","selPhoSieie_pt20",50, 0, 0.025 );
    hist1d[15] = new TH1D("selPhoSieip_pt20","selPhoSieip_pt20",20, -0.0005, 0.0005 );
    hist1d[16] = new TH1D("selPhoSipip_pt20","selPhoSipip_pt20",50, 0, 0.025 );

    hist1d[17] = new TH1D("selPhoSigmaIEtaIEta_pt30","selPhoSigmaIEtaIEta_pt30",50, 0, 0.025 );
    hist1d[18] = new TH1D("selPhoClstrRn_pt30","selPhoClstrRn_pt30",100, 0, 1 );
    hist1d[19] = new TH1D("selPhoCovEtaEta_pt30","selPhoCovEtaEta_pt30",100, 0, 0.001 );
    hist1d[20] = new TH1D("selPhoCovEtaPhi_pt30","selPhoCovEtaPhi_pt30",200, -0.001, 0.001 );
    hist1d[21] = new TH1D("selPhoCovPhiPhi_pt30","selPhoCovPhiPhi_pt30",100, 0, 0.001 );
    hist1d[22] = new TH1D("selPhoNrh_pt30","selPhoNrh_pt30",150, 0, 150 );
    hist1d[23] = new TH1D("selPhoEtaWidth_pt30","selPhoEtaWidth_pt30",50, 0, 0.05 );
    hist1d[24] = new TH1D("selPhoPhiWidth_pt30","selPhoPhiWidth_pt30",250, 0, 0.25);
    hist1d[25] = new TH1D("selPhoR9_pt30","selPhoR9_pt30",100, 0, 1 );
    hist1d[26] = new TH1D("selPhoS4_pt30","selPhoS4_pt30",100, 0, 1 );
    hist1d[27] = new TH1D("selPhoSAlp_pt30","selPhoSAlp_pt30",60, -1.5, 1.5 );
    hist1d[28] = new TH1D("selPhoSMaj_pt30","selPhoSMaj_pt30",100, 0, 5 );
    hist1d[29] = new TH1D("selPhoSMin_pt30","selPhoSMin_pt30",50, 0, 2.5 );
    hist1d[30] = new TH1D("selPhoSieie_pt30","selPhoSieie_pt30",50, 0, 0.025 );
    hist1d[31] = new TH1D("selPhoSieip_pt30","selPhoSieip_pt30",20, -0.0005, 0.0005 );
    hist1d[32] = new TH1D("selPhoSipip_pt30","selPhoSipip_pt30",50, 0, 0.025 );

    hist1d[33] = new TH1D("selPhoSigmaIEtaIEta_pt100","selPhoSigmaIEtaIEta_pt100",50, 0, 0.025 );
    hist1d[34] = new TH1D("selPhoClstrRn_pt100","selPhoClstrRn_pt100",100, 0, 1 );
    hist1d[35] = new TH1D("selPhoCovEtaEta_pt100","selPhoCovEtaEta_pt100",100, 0, 0.001 );
    hist1d[36] = new TH1D("selPhoCovEtaPhi_pt100","selPhoCovEtaPhi_pt100",200, -0.001, 0.001 );
    hist1d[37] = new TH1D("selPhoCovPhiPhi_pt100","selPhoCovPhiPhi_pt100",100, 0, 0.001 );
    hist1d[38] = new TH1D("selPhoNrh_pt100","selPhoNrh_pt100",150, 0, 150 );
    hist1d[39] = new TH1D("selPhoEtaWidth_pt100","selPhoEtaWidth_pt100",50, 0, 0.05 );
    hist1d[40] = new TH1D("selPhoPhiWidth_pt100","selPhoPhiWidth_pt100",250, 0, 0.25);
    hist1d[41] = new TH1D("selPhoR9_pt100","selPhoR9_pt100",100, 0, 1 );
    hist1d[42] = new TH1D("selPhoS4_pt100","selPhoS4_pt100",100, 0, 1 );
    hist1d[43] = new TH1D("selPhoSAlp_pt100","selPhoSAlp_pt100",60, -1.5, 1.5 );
    hist1d[44] = new TH1D("selPhoSMaj_pt100","selPhoSMaj_pt100",100, 0, 5 );
    hist1d[45] = new TH1D("selPhoSMin_pt100","selPhoSMin_pt100",50, 0, 2.5 );
    hist1d[46] = new TH1D("selPhoSieie_pt100","selPhoSieie_pt100",50, 0, 0.025 );
    hist1d[47] = new TH1D("selPhoSieip_pt100","selPhoSieip_pt100",20, -0.0005, 0.0005 );
    hist1d[48] = new TH1D("selPhoSipip_pt100","selPhoSipip_pt100",50, 0, 0.025 );

	hist1d[50] = new TH1D("rtgtdif","Diff( RecoTime - GenTime )",500,-5,5 );
    hist1d[51] = new TH1D("profile_srtvrege","profile_srtvrege;SRecoTime [ns]",100,-5,20);
    hist1d[52] = new TH1D("fit_srtvrege","fit_srtvrege;SRecoTime [ns]",100,-5,20);


	hist1d[60] = new TH1D("rlalb","(la-lb)/(la+lb);(la-lb)/(la+lb)",200,-1,1);
    hist1d[61] = new TH1D("mllpe","M_{llp}^{la=lb}; M [GeV]",260,-100,2500);
    hist1d[62] = new TH1D("mllp_pm","M_{llp}^{promt}; M [GeV]",260,-100,2500);
    hist1d[63] = new TH1D("mllpeb","M_{llp}^{la=lb} notSig;M [GeV]",260,-100,2500);
    hist1d[64] = new TH1D("mllp_pmb","M_{llp}^{promt} notSig; M [GeV]",260,-100,2500);
    hist1d[65] = new TH1D("phoe","Cor E_{pho};E [GeV]",50,0,1000);
    hist1d[66] = new TH1D("phopt","Cor Pt_{pho};Pt [GeV]",50,0,1000);
    hist1d[67] = new TH1D("phoeb","Cor E_{pho} notSig;E [GeV]",50,0,1000);
    hist1d[68] = new TH1D("phoptb","Cor Pt_{pho} notSig;Pt [GeV]",50,0,1000);

    hist1d[69] = new TH1D("beta","Beta;Beta",120,-1,5);
    hist1d[70] = new TH1D("e1beta","Beta(m2) La=Lb;Beta",120,-1,5);
    hist1d[71] = new TH1D("e2beta","Beta(m3) Prompt;Beta",120,-1,5);
    hist1d[72] = new TH1D("tmeas","Tmeas;Tmeas [cm]",110,-100,1000);
    hist1d[73] = new TH1D("m1","M1;M1",260,-100,2500);

    hist1d[74] = new TH1D("betabg","Beta BG;Beta",120,-1,5);
    hist1d[75] = new TH1D("e1betabg","Beta(m2) La=Lb BG;Beta",120,-1,5);
    hist1d[76] = new TH1D("e2betabg","Beta(m3) Prompt BG;Beta",120,-1,5);
    hist1d[77] = new TH1D("tmeasbg","Tmeas BG;Tmeas [cm]",110,-100,1000);
    hist1d[78] = new TH1D("m1bg","M1 BG;M1 [GeV]",260,-100,2500);

    hist1d[79] = new TH1D("rtime", "rtime", 100, -5, 20);
    hist1d[80] = new TH1D("rtimeb", "rtime BG", 100, -5, 20);

	//hisst1d 400 + for special use ( see end jobs )
	
    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	hist2d[0] = new TH2D("drvrt","dR vs RecoTime;dR;RecoTime [ns]",100,0,0.5,100,-5,20);
    hist2d[1] = new TH2D("srtvgt","SRecoTime vs GenTime;SRecoTime [ns];GenTime [ns]",100,-5,20,100,-5,20);
    hist2d[2] = new TH2D("regevrt","Ratio(RecoE/GenE) vs RecoTime; reco E / gen E; RecoTime [ns]",200,0.2,2.2,100,-5,20);
    hist2d[3] = new TH2D("regevgt","Ratio(RecoE/GenE) vs GenTime; reco E / gen E; GenTime [ns]",200,0.2,2.2,100,-5,20);
    hist2d[4] = new TH2D("evrt","Pho Eta vs RecoTime;eta;RecoTime [ns]",300,-1.5,1.5,100,-5,20);
    hist2d[5] = new TH2D("evgt","Pho Eta vs GenTime;eta;GenTime [ns]",300,-1.5,1.5,100,-5,20);
    hist2d[6] = new TH2D("drvgt","dR vs GenTime;dR;GenTime [ns]",100,0,0.5,100,-5,20);
    hist2d[7] = new TH2D("drvsrt","dR vs SReco Time;dR;SRecoTime [ns]",100,0,0.5,100,-5,20);
    hist2d[8] = new TH2D("drvrege","dR vs Ratio(RecoE/GenE);dR;reco E / gen E",100,0,0.5,200,0.2,2.2);
    hist2d[9] = new TH2D("sregevrt","Ratio(SRecoE/GenE) vs SRecoTime; sreco E / gen E; SRecoTime [ns]",200,0.2,2.2,100,-5,20);
    hist2d[10] = new TH2D("evsrt","Pho Eta vs SReco Time;eta;SRecoTime [ns]",300,-1.5,1.5,100,-5,20);
    hist2d[11] = new TH2D("srtvrege","SRecoTime v Ratio(SRecoE/GenE);SRecoTime [ns];Ratio(SRecoE/GenE)",100,-5,20,200,0.2,2.2);
    hist2d[12] = new TH2D("sregevrtbg","Ratio(SRecoE/GenE) vs SRecoTime BG; sreco E / gen E; SRecoTime [ns]",200,0.2,2.2,100,-5,20);
    hist2d[13] = new TH2D("cegevrt","Ratio(CRecoE/GenE) vs RecoTime; cor-reco E / gen E; RecoTime [ns]",200,0.2,2.2,100,-5,20);

    hist2d[50] = new TH2D("covievsie_pt20","CovEtaEtaVsSieie_pt20;CovEtaEta;Sieie^{2}", 100, 0, 0.001, 100, 0, 0.001 );
    hist2d[51] = new TH2D("salpvsie_pt20","SAlpVsSieie_pt20;SAlp;Sieie", 100, -2.5, 2.5, 50, 0, 0.025 );

    hist2d[60] = new TH2D("mllpevrt","M_{llp}^{la=lb} v dt_{mean-nom};M [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[61] = new TH2D("mllppvrt","M_{llp}^{promt} v dt_{mean-nom};M [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[62] = new TH2D("mllpevrtb","M_{llp}^{la=lb} v dt_{mean-nom} notSig;M [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[63] = new TH2D("mllppvrtb","M_{llp}^{promt} v dt_{mean-nom} notSig;M [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);

    hist2d[64] = new TH2D("mllpevmllpp","M_{llp}^{la=lb} v M_{llp}^{promt};M_{llp}^{la=lb} [GeV];M_{llp}^{promt} [GeV]",200,0,500,200,0,500);
    hist2d[65] = new TH2D("mllpeebvmllppb","M_{llp}^{la=lb} v M_{llp}^{promt} notSig;M_{llp}^{la=lb} [GeV];M_{llp}^{promt} [GeV]",200,0,500,200,0,500);
    hist2d[66] = new TH2D("cevge","Corr Energy Vs Gen Energy; E_{cor} [GeV];E_{gen} [GeV]",500,0,1000,500,0,1000);
    hist2d[67] = new TH2D("cevges","Corr Energy Vs Gen Energy Sig; E_{cor} [GeV];E_{gen} [GeV]",500,0,1000,500,0,1000);
    hist2d[68] = new TH2D("cevgeb","Corr Energy Vs Gen Energy notSig; E_{cor} [GeV];E_{gen} [GeV]",500,0,1000,500,0,1000);

    hist2d[69] = new TH2D("crevrt","Corr Reco E v dt_{mean-nom};E [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[70] = new TH2D("crevrtb","Corr Reco E v dt_{mean-nom} notSig;E [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[71] = new TH2D("crptvrt","Corr Reco Pt v dt_{mean-nom};Pt [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[72] = new TH2D("crptvrtb","Corr Reco Pt v dt_{mean-nom} notSig;Pt [GeV];dt_{mean-nom} [ns]",100,0,1000,210,-1,20);
    hist2d[73] = new TH2D("crevmllpe","Corr Reco E v M_{llp}^{la=lb};E [GeV];M_{llp}^{la=lb} [GeV]",100,0,1000,100,0,1000);
    hist2d[74] = new TH2D("crptvmllpe","Corr Reco Pt M_{llp}^{la=lb};Pt [GeV];M_{llp}^{la=lb} [GeV]",100,0,1000,100,0,1000);

    //for( int it = 0; it < 40; it++ ){
    //    std::string temp_title = "rregvre" + std::to_string( it );
    //    std::string temp_name = "Ratio(RecoE/GenE) V RTime by RecoE"+ std::to_string( it*25 )+" ;Ratio(SRecoE/GenE);RTime [ns]";
    //    hist2d[100+it] = new TH2D(temp_title.c_str(),temp_name.c_str(),100,0,1,250,-5,20);
    //}//<<>>for( int it = 100; it < 10; it++ )

    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

	//------------------------------------------------------------------------------------

}//<<>>void HistMaker::initHists()

//void HistMaker::histMaker( std::string indir, std::string infilelist, std::string outfilename )

int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                //const std::string listdir = "skims_files/";
                const std::string listdir = "root://cmseos.fnal.gov//store/user/jaking/llpSkims/";
				auto infilename = "KUCMS_Master_Skim_List.txt";

                //auto outfilename = "KUCMS_ShapeVar_Sig_L400_SkimHists_v4.root"; //9
                //auto outfilename = "KUCMS_ShapeVar_QCD_SkimHists_v3.root"; //9
                //auto outfilename = "KUCMS_ShapeVar_GJets_SkimHists_v4.root"; //9
                //auto outfilename = "KUCMS_ShapeVar_JetHt_Bkg_SkimHists.root"; //9
				auto outfilename = "KUCMS_ShapeVar_Sig_L400_GJets_SkimHists_v5.root";

				auto htitle = "KUCMS_ShapeVar_Hists_";

                HistMaker base;
                base.histMaker( listdir, infilename, outfilename, htitle );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )


