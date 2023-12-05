//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////


#include "skimHistMaker.hh"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TLatex.h"

//#define DEBUG true
#define DEBUG false
#define doEBEEmaps false

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// HistMaker class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

void HistMaker::histMaker( std::string indir, std::vector<std::string> infilelists, std::string outfilename, std::string htitle ){

    //bool debug = true;
    bool debug = false;

    std::cout << "Setting up For Main Loop." << std::endl;
    Nsample = infilelists.size();   
	initHists(htitle,Nsample);

    const std::string disphotreename("kuSkimTree");
    const std::string configtreename("kuSkimConfigTree");
    const std::string eosdir("");
    const std::string listdir("");

	std::cout << "Producing Histograms for : " << outfilename << std::endl;
	int curSample = 0;
	for( auto infilelist : infilelists ){

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

    	Init(fInTree);

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
			eventLoop(entry,curSample);
	    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop
		curSample++;
        delete fInTree; // = new TChain(disphotreename.c_str());
        delete fConfigTree; // = new TChain(configtreename.c_str());

	}//for( infilelist : infilelists )
   
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
void HistMaker::eventLoop( Long64_t entry, int chist ){

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
        if( std::abs((*selPhoEta)[it]) > 1.442 ) continue;

		bool jetphoton = false;
		for( int jit = 0; jit < nSelJets; jit++ ){

			float dpjeta = (*selJetEta)[jit] - (*selPhoEta)[it];
			float dpjphi = dPhi( (*selJetPhi)[jit], (*selPhoPhi)[it] );	
			float dr = hypo( dpjeta, dpjphi );
			if( dr < 0.4 ) jetphoton = true;

		} // for( int jit = 0; jit < nSelJets; jit++ )
		//if( jetphoton ) continue; //{ std::cout << "Jet-Photon Overlap!!!" << std::endl; continue; }

        if( DEBUG ) std::cout << " -- looping photons : Start getting dr dp info" << std::endl;
		auto phoClass = (*selPhoQuality)[it];
        //if( phoClass < 1 ) continue;
        auto phoSusId = (*selPhoSusyId)[it]; // selPhoSusyId->at(it);
        auto phoGenMatDr = (*selPhoGenDr)[it];
        auto phoGenMatDp = (*selPhoGenDp)[it];
        auto phoTime = (*selPhoTime)[it];
        auto phoOOT = (*selPhoOOT)[it];

		//if( phoSusId != 22 ) continue;

        if( DEBUG ) std::cout << " -- looping photons : getting susy ids for " << selPhoSusyId->size() << std::endl;        

        bool hcalsum = true;
        bool tsptscdr4 = (*selPhoTrkSumPtSolidConeDR04)[it] < 6.0; //(*selPhoTrkSumPtSolidConeDR04)[it] < cutvalue;
		bool ecalrhsum = (*selPhoEcalRHSumEtConeDR04)[it] < 10.0;
        bool htoem = (*selPhoHadTowOverEM)[it] < 0.02;
        bool isoskip = not ( htoem && tsptscdr4 && ecalrhsum && hcalsum );
        if( isoskip ) continue;

		hist1d[chist]->Fill((*selPhoClstrRn)[it]);


    }//<<>>for( int it = 0; it < nPhotons; it++ )

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

	std::string g_Xname("ClstrRn");
    std::string g_PlotTitle("Signal Photon ClstrRn");
    std::vector<std::string> g_Samples;
	g_Samples.push_back("GMSB");
    //g_Samples.push_back("JetHt");
    g_Samples.push_back("GJets");

	double max = -1.;
	int imax = -1;
	for(int i = 0; i < Nsample; i++){
	  	//std::cout << g_Samples[i]->GetTitle().c_str() << " " << hist1d[i]->Integral()*137 << " total" << endl;
	  	hist1d[i]->Scale(1./hist1d[i]->Integral());
	  	//hist1d[i]->Scale(1./hist[i]->GetMaximum());
	  	if(hist1d[i]->GetMaximum() > max){
	    	max = hist1d[i]->GetMaximum();
	    	imax = i;
	  	}//if(hist1d[i]->GetMaximum() > max)
	}//for(int i = 0; i < Nsample; i++)

	float width = hist1d[0]->GetBinWidth(1);
	char *yaxis = new char[100];
	sprintf(yaxis,"Events / %f", width);

	//auto gStyle = new TStyle("Default","Default Style");
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(11111111);
	TCanvas* can = (TCanvas*) new TCanvas("can","can",600.,500);

	can->SetLeftMargin(0.13);
	can->SetRightMargin(0.04);
	can->SetBottomMargin(0.16);
	can->SetTopMargin(0.085);
	can->SetGridx();
	can->SetGridy();
	can->Draw();
	can->cd();
	// for(int i = 0; i < hist1d[imax]->GetNbinsX(); i++){
	//   char* sbin = new char[20];
	//   sprintf(sbin,"#geq %d", i);
	//   hist1d[imax]->GetXaxis()->SetBinLabel(i+1,sbin);
	// }

	hist1d[imax]->Draw("hist");
	hist1d[imax]->GetXaxis()->CenterTitle();
	hist1d[imax]->GetXaxis()->SetTitleFont(42);
	hist1d[imax]->GetXaxis()->SetTitleSize(0.06);
	hist1d[imax]->GetXaxis()->SetTitleOffset(1.06);
	hist1d[imax]->GetXaxis()->SetLabelFont(42);
	hist1d[imax]->GetXaxis()->SetLabelSize(0.05);
	hist1d[imax]->GetXaxis()->SetTitle(g_Xname.c_str());
	hist1d[imax]->GetYaxis()->CenterTitle();
	hist1d[imax]->GetYaxis()->SetTitleFont(42);
	hist1d[imax]->GetYaxis()->SetTitleSize(0.06);
	hist1d[imax]->GetYaxis()->SetTitleOffset(1.1);
	hist1d[imax]->GetYaxis()->SetLabelFont(42);
	hist1d[imax]->GetYaxis()->SetLabelSize(0.05);
	hist1d[imax]->GetYaxis()->SetTitle("a. u.");
	hist1d[imax]->GetYaxis()->SetRangeUser(0., hist1d[imax]->GetMaximum()*1.1);
	//hist1d[imax]->GetYaxis()->SetTitle(yaxis);
	//hist1d[imax]->GetYaxis()->SetTitle("N_{evt} / fb^{-1}");
	int Ntype[3];

	int mycolor[8];
	mycolor[0] = kBlue+2;
	mycolor[1] = kGreen+3;
	mycolor[2] = kRed+1;
	mycolor[3] = kYellow+2;
	mycolor[4] = kMagenta+1;
	mycolor[5] = kMagenta+2;
	mycolor[6] = kCyan+2;
	mycolor[7] = kCyan+3;

	Ntype[0] = 0;
	for(int i = Nsample-1; i >= 0; i--){
		hist1d[i]->SetLineColor(mycolor[i]);
		hist1d[i]->SetLineWidth(3);
		if(i >= 4 && false){
	    	if(i%2 == 0) hist1d[i]->SetLineStyle(7);
	    	if(i%2 == 1) hist1d[i]->SetLineStyle(9);
	    	hist1d[i]->SetLineWidth(4);
	  	}// if(i >= 4 && false)
	  	hist1d[i]->SetMarkerColor(mycolor[i]);
	  	hist1d[i]->SetMarkerSize(0);
	  	hist1d[i]->SetFillColor(kWhite);
	  	Ntype[0]++;
	  	hist1d[i]->Draw("hist SAME");
	}//for(int i = Nsample-1; i >= 0; i--)

	TLegend* leg = new TLegend(0.688,0.22,0.93,0.42);
	leg->SetTextFont(42);
	leg->SetTextSize(0.04);
	leg->SetFillColor(kWhite);
	leg->SetLineColor(kWhite);
	leg->SetShadowColor(kWhite);
	for(int i = 0; i < Nsample; i++){ leg->AddEntry(hist1d[i],g_Samples[i].c_str()); }
	leg->SetLineColor(kWhite);
	leg->SetFillColor(kWhite);
	leg->SetShadowColor(kWhite);
	leg->Draw("SAME");

	TLatex l;
	l.SetTextFont(132);
	l.SetNDC();
	l.SetTextSize(0.045);
	l.SetTextFont(42);
	// l.DrawLatex(0.17,0.855,g_PlotTitle.c_str());
	l.DrawLatex(0.6,0.943,g_PlotTitle.c_str());
	l.SetTextSize(0.045);
	l.SetTextFont(42);
	l.DrawLatex(0.135,0.943,"#bf{CMS} Simulation Preliminary");
	l.SetTextSize(0.04);
	l.DrawLatex(0.69,0.85,"#tilde{#chi}_{2}^{0} #tilde{#chi}_{1}^{#pm} #rightarrow Z #tilde{#chi}_{1}^{0} W #tilde{#chi}_{1}^{0}");
	l.DrawLatex(0.69,0.85,"#tilde{#chi}_{1}^{#pm} #tilde{#chi}_{1}^{#pm} #rightarrow W* #tilde{#chi}_{1}^{0} W* #tilde{#chi}_{1}^{0}");
	l.DrawLatex(0.69,0.85,"#tilde{t} #tilde{t} #rightarrow b W* #tilde{#chi}_{1}^{0} b W* #tilde{#chi}_{1}^{0}");
	l.DrawLatex(0.69,0.85,"#tilde{t} #tilde{t} #rightarrow b #tilde{#chi}^{#pm}_{1}(f #bar{f} #tilde{#chi}^{0}_{1}) b #tilde{#chi}^{#pm}_{1}(f #bar{f} #tilde{#chi}^{0}_{1})");
	l.DrawLatex(0.69,0.85,"#tilde{t} #tilde{t} #rightarrow b f #bar{f} #tilde{#chi}^{0}_{1} b f #bar{f} #tilde{#chi}^{0}_{1}");
	l.DrawLatex(0.69,0.85,"#tilde{l} #tilde{l} #rightarrow l #tilde{#chi}^{0}_{1} l #tilde{#chi}^{0}_{1}");

	// l.SetTextSize(0.045);
	// l.SetTextFont(132);
	// string bla = "#scale[0.6]{#int} #it{L dt} = "+to_string(int(g_lumi))+" fb^{-1},  #Delta_{N#scale[0.8]{bkg}} = ";
	// bla += to_string(int(g_deltaNbkg))+" %";
	// l.DrawLatex(0.61,0.943,bla.c_str());

	//std::cout << "Press 'q' to continue..." << std::endl;
	//while( std::cin.get() != 'q' ){};

}//<<>>void HistMaker::endJobs()

void HistMaker::initHists( std::string ht, int nHists){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------
	
	int nbins = 100;
	double str = 0;
	double end = 1;
	for( int it = 0; it < nHists; it++ ){

		std::string sname("crStyleSkimHist_");
		sname += std::to_string(it);
    	hist1d[it] = new TH1D(sname.c_str(),sname.c_str(), nbins, str, end ); 

	}//<<>>for( int it = 0; it < nHists; it++ )

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

//int main ( int argc, char *argv[] ){
void crStyleSkimHistMaker(){ 

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "skims_files/";
				std::vector<std::string> infileLists;
				infileLists.push_back("KUCMS_GMSB_Skim_List2.txt");
                //infileLists.push_back("KUCMS_JetHT_Skim_List.txt");
                infileLists.push_back("KUCMS_GJets_Skim_List.txt");

                auto outfilename = "KUCMS_crStyle_SkimHists.root"; //9

				auto htitle = "KUCMS_crStyle_SkimHists_";

                HistMaker base;
                base.histMaker( listdir, infileLists, outfilename, htitle );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )

