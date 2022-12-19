#include "Skimmer.hh"
#include "TROOT.h"
#include <string> 
#include <iostream>
#include <sstream> 
#include <TRandom.h>

enum ECAL {EB, EM, EP, NONE};
//std::string ecal_config_path("/home/t3-ku/jaking/ecaltiming/CMSSW_10_2_5/src/Timing/TimingAnalyzer/macros/ecal_config/");
std::string ecal_config_path("/uscms/home/jaking/nobackup/ecaltiming/CMSSW_11_3_0_pre6/src/Timing/TimingAnalyzer/macros/ecal_config/");

struct DetIDStruct
{ 
  DetIDStruct() {}
  DetIDStruct(const Int_t inI1, const Int_t inI2, const Int_t inTT, const Int_t & inEcal) : i1(inI1), i2(inI2), TT(inTT), ecal(inEcal)  {}

  Int_t i1; // EB: iphi, EE: ix
  Int_t i2; // EB: ieta, EE: iy
  Int_t TT; // trigger tower
  Int_t ecal; // EB, EM, EP
};

void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )
{
    const std::string detIDConfigEB(ecal_config_path+"fullinfo_detids_EB.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    TString pos;

    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM)
    {
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << EB << std::endl;
        DetIDMap[cmsswId] = {iphi,ieta,TT25,EB};
        auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )
{
    const std::string detIDConfigEE(ecal_config_path+"fullinfo_detids_EE.txt");
    std::ifstream infile( detIDConfigEE, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    TString EE;

    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant)
    {
        ECAL ec = EM;
        if( side > 0 ) ec = EP;
        //std::cout << "DetID Input Line: " << cmsswId << " " << ix << " "  << iy << " " << ec << std::endl; 
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec};
        auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

void NormX(TH2F *& hist){ 	//, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY){
 
	std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;
    
    for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++){

    	//const auto binwidthX = hist->GetXaxis()->GetBinWidth(ibinX);
        const auto norm = hist->Integral(ibinX,ibinX,1,hist->GetYaxis()->GetNbins());
		if( norm == 0.0 ) continue;
        for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++){
        	//const auto binwidthY = hist->GetYaxis()->GetBinWidth(ibinY);
               
            //get multiplier/divisor
            //auto multiplier = 1.f;      
            //if (varBinsX) multiplier *= binwidthX;
            //if (varBinsY) multiplier *= binwidthY;
               
            // get content/error
            auto content = hist->GetBinContent(ibinX,ibinY);
            auto error   = hist->GetBinError  (ibinX,ibinY);
               
            // scale it
            //if (isUp){ 
            //   content *= multiplier;
            //   error   *= multiplier;
            //} else {
            content /= norm;
            error   /= norm;
            //}
               
            // set new contents
            hist->SetBinContent(ibinX,ibinY,content);
            hist->SetBinError  (ibinX,ibinY,error);
        }
	}
}

void setBins(std::string & str, std::vector<Double_t> & bins, Bool_t & var_bins){

  	if(str.find("CONSTANT") != std::string::npos){

    	var_bins = false;
    	str = Common::RemoveDelim(str,"CONSTANT");
    	Int_t nbins = 0; Double_t low = 0.f, high = 0.f;
    	std::stringstream ss(str);
    	ss >> nbins >> low >> high;
    	Double_t bin_width = (high-low)/nbins;
      	//std::cout << "Setting Const bins : ";
    	for (Int_t ibin = 0; ibin <= nbins; ibin++){
    		//std::cout << low+ibin*bin_width << " ";
    		bins.push_back(low+ibin*bin_width);
    	}//<<>.for (Int_t ibin = 0; ibin <= nbins; ibin++)
    	//std::cout << std::endl;

	} else if(str.find("VARIABLE") != std::string::npos) {

    	var_bins = true;
    	str = Common::RemoveDelim(str,"VARIABLE");
    	Float_t bin_edge;
    	std::stringstream ss(str);
    	//std::cout << "Setting Var bins : ";
    	while(ss >> bin_edge){
    		//std::cout << bin_edge << " ";
    		bins.push_back(bin_edge);
    	}//<<>>while (ss >> bin_edge)
    	//std::cout << std::endl;

  	} else {

    	std::cerr << "Aye... bins are either VARIABLE or CONSTANT! Exiting..." << std::endl;
    	exit(1);

	}//<<>>if      (str.find("CONSTANT") != std::string::npos)

}//<<>>void setBins(std::string & str, std::vector<Double_t> & bins, Bool_t & var_bins)

void scaleHist(TH2F *& hist, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY){

	std::cout << "Scaling " << (isUp?"up":"down") << " hist: " << hist->GetName() << std::endl;

    for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++){

    	const auto binwidthX = hist->GetXaxis()->GetBinWidth(ibinX);
        for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++){

        	const auto binwidthY = hist->GetYaxis()->GetBinWidth(ibinY);
            // get multiplier/divisor
            auto multiplier = 1.f;
            if( varBinsX ) multiplier *= binwidthX;
            if( varBinsY ) multiplier *= binwidthY;
			auto scale = ( not isUp ) ? multiplier : 1/multiplier; 
            // get content/error
            auto content = hist->GetBinContent(ibinX,ibinY)*scale;
            auto error = hist->GetBinError(ibinX,ibinY)*scale;
			// set new contents
            hist->SetBinContent(ibinX,ibinY,content);
            hist->SetBinError  (ibinX,ibinY,error);

        }//<<>>for (auto ibinY = 1; ibinY <= hist->GetYaxis()->GetNbins(); ibinY++)
	}//<<>>for (auto ibinX = 1; ibinX <= hist->GetXaxis()->GetNbins(); ibinX++)

}//<<>>void scaleHist(TH2D *& hist, const Bool_t isUp, const Bool_t varBinsX, const Bool_t varBinsY)

//void plot2dResolution( string indir, string infilelistname, string outfilename, string tvarname, string calimapname, string isd_type, int brun, int erun, int leta, int heta  ){
void plot2dResolution( string indir, string infilelistname, string outfilename, string tvarname, string calimapname, string isd_type ){

	bool debug = false;
	TRandom* getRandom = new TRandom();
	getRandom->SetSeed(0);
	unsigned int nRand(100);

    std::cout << "opening output file" << std::endl;
    string histoutfilename(outfilename+".root");
    TFile* fOutFile = new TFile( histoutfilename.c_str(), "RECREATE" );
    std::cout << "fOutFile : " << fOutFile << std::endl;

    double phoseedtimeCaliIc_0(0.0);
    double phoseedtimeCaliIcErr_0(0.0);
    double phoseedtimeCaliIc_1(0.0);
    double phoseedtimeCaliIcErr_1(0.0);

    // Declaration of leaf types
    UInt_t          run;
    //UInt_t          lumi;
    //vector<unsigned int> *rhCaliID;
    //vector<float>   *rhCaliRtTime;
    //vector<float>   *rhCaliCCTime;
    vector<unsigned int> *resResVecRhID;
    vector<float>   *resAmp;
    vector<float>   *resE;
    vector<float>   *resRtTime;
    vector<float>   *resCCtime;
    vector<float>   *resTOF;

    // List of branches
    TBranch        *b_run;   //!
    //TBranch        *b_lumi;   //!
    //TBranch        *b_rhCaliID;   //!
    ///TBranch        *b_rhCaliRtTime;   //!
    //TBranch        *b_rhCaliCCTime;   //!
    TBranch        *b_resResVecRhID;   //!
    TBranch        *b_resAmp;   //!
    TBranch        *b_resE;   //!
    TBranch        *b_resRtTime;   //!
    TBranch        *b_resCCtime;   //!
    TBranch        *b_resTOF;   //!


    // >> calcs  <<
    std::cout << "Setting up 2D plot" << std::endl;

	string locsname("SRO_");
    string locdname("DRO_");
    string globname("ZEE_");

    string histname("Data_Hist");
    string histnameg("Data_GHist");
    string histnames("Data_SHist");
    string histname0("Data_OthHist");
    string histname1("td_Hist");
    string histname2("tdc_Hist");
    string histname3("teffa_Hist");
    string histname4("tdfc_Hist");
    string histname5("td_eta_Hist");
    string histname6("td_phi_Hist");
    string fTitle("#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)");
    string fTitle0("#Delta(Photon Seed Time) [ns] vs. A_{0}/#sigma_{n} (EBEB)");
    string fTitle1("Photon Seed Time [ns]");
    string fTitle2("Photon Seed Time Calibrated [ns]");
    string fTitle3("Amplitude 0 vs Amplitude 1");
    string fTitle4("Photon Seed Time Filtered Calibrated [ns]");
    string fTitle5("#Delta(Photon Seed Time) [ns] vs iEta");
    string fTitle6("#Delta(Photon Seed Time) [ns] vs iPhi");
    string fXTitle("A_{eff}/#sigma_{n} (EBEB)");
    std::vector<Double_t> fXBins;
    Bool_t fXVarBins = false;
    //string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 750 950 1275 1700 2250");
    string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250");
    int nMyBins = 13;
    //std::vector<TString> fXLabels;
    string fYTitle("#Delta(Photon Seed Time) [ns] (EBEB)");
    std::vector<Double_t> fYBins;
    Bool_t fYVarBins = false;
    //string ybinstr("CONSTANT 1920 -3 3");
    string ybinstr("CONSTANT 600 -3 3");
    //std::vector<TString> fYLabels;
    string fZTitle("");
    //string yetabinstr("CONSTANT 800 -4 4");
    //string xetabinstr("CONSTANT 181 -90.5 90.5");
    //string xphibinstr("CONSTANT 362 -0.5 361.5");
    std::vector<Double_t> fXetaBins;
    std::vector<Double_t> fYetaBins;


    setBins(xbinstr,fXBins,fXVarBins);
    setBins(ybinstr,fYBins,fYVarBins);
    //Common::SetupBins(xetabinstr,fXetaBins,fXVarBins);
    //Common::SetupBins(yetabinstr,fYetaBins,fYVarBins);

    const auto xbins = &fXBins[0];
    const auto ybins = &fYBins[0];

    int tdiv = 600;
    float tstart = -3.0;
    float tend = 3.0;

    std::string tehistname = "timeErrHist";
	auto timeErrBinHist = new TH1F(tehistname.c_str(), tehistname.c_str(),1000,0,1);
    TH1F * caliErrBinHist[nMyBins+1];
    std::string eahistname = "effAHist";
    caliErrBinHist[0] = new TH1F(eahistname.c_str(),eahistname.c_str(),225,0,2250);
    caliErrBinHist[0]->Sumw2();
    for( auto ibin = 1; ibin < nMyBins+1; ibin++ ){
		auto binhistname = "bin" + to_string(ibin) + "ErrHist";
		caliErrBinHist[ibin] = new TH1F(binhistname.c_str(),binhistname.c_str(),500,0,0.5);
        caliErrBinHist[ibin]->Sumw2();
    }

	// local same
    auto theHistLS = new TH2F((locsname+histname)).c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto theGHistLS = new TH2F((locsname+histnameg).c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto theSHistLS = new TH2F((locsname+histnames).c_str(),fTitle.c_str(), 2250, 0, 2250, tdiv, tstart, tend);
    auto theOthHistLS = new TH2F((locsname+histname0).c_str(),fTitle0.c_str(), 2250, 0, 2250, tdiv, tstart, tend);
    auto thetdHistLS = new TH1F((locsname+histname1).c_str(),fTitle1.c_str(),tdiv,tstart,tend);
    auto thetdcHistLS = new TH1F((locsname+histname2).c_str(),fTitle2.c_str(),tdiv,tstart,tend);

    //auto thetdfcHist = new TH1F(histname4.c_str(),fTitle4.c_str(),tdiv,tstart,tend);
    auto theEffaHistLS = new TH2F((locsname+histname3).c_str(),fTitle3.c_str(),2250,0,2250,2250,0,2250);

    auto theetaHistLS = new TH2F((locsname+histname5).c_str(),fTitle5.c_str(),182,-90.5,90.5,tdiv,tstart,tend);
    auto thephiHistLS = new TH2F((locsname+histname6).c_str(),fTitle6.c_str(),363,-1.5,361.5,tdiv,tstart,tend);

    theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
    theGHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    theGHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    theGHistLS->GetZaxis()->SetTitle(fZTitle.c_str());

	// local diffrnent
    auto theHistLD = new TH2F((locdname+histname).c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto theGHistLD = new TH2F((locdname+histnameg).c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto theSHistLD = new TH2F((locdname+histnames).c_str(),fTitle.c_str(), 2250, 0, 2250, tdiv, tstart, tend);
    auto theOthHistLD = new TH2F((locdname+histname0).c_str(),fTitle0.c_str(), 2250, 0, 2250, tdiv, tstart, tend);
    auto thetdHistLD = new TH1F((locdname+histname1).c_str(),fTitle1.c_str(),tdiv,tstart,tend);
    auto thetdcHistLD = new TH1F((locdname+histname2).c_str(),fTitle2.c_str(),tdiv,tstart,tend);
    //auto thetdfcHist = new TH1F(histname4.c_str(),fTitle4.c_str(),tdiv,tstart,tend);
    auto theEffaHistLD = new TH2F((locdname+histname3).c_str(),fTitle3.c_str(),2250,0,2250,2250,0,2250);

    auto theetaHistLD = new TH2F((locdname+histname5).c_str(),fTitle5.c_str(),182,-90.5,90.5,tdiv,tstart,tend);
    auto thephiHistLD = new TH2F((locdname+histname6).c_str(),fTitle6.c_str(),363,-1.5,361.5,tdiv,tstart,tend);

    theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());
    theGHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
    theGHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
    theGHistLD->GetZaxis()->SetTitle(fZTitle.c_str());

	// global
    auto theHistGB = new TH2F((globname+histname).c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto theGHistGB = new TH2F((globname+histnameg).c_str(),fTitle.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    auto theSHistGB = new TH2F((globname+histnames).c_str(),fTitle.c_str(), 2250, 0, 2250, tdiv, tstart, tend);
    auto theOthHistGB = new TH2F((globname+histname0).c_str(),fTitle0.c_str(), 2250, 0, 2250, tdiv, tstart, tend);
    auto thetdHistGB = new TH1F((globname+histname1).c_str(),fTitle1.c_str(),tdiv,tstart,tend);
    auto thetdcHistGB = new TH1F((globname+histname2).c_str(),fTitle2.c_str(),tdiv,tstart,tend);

    //auto thetdfcHist = new TH1F(histname4.c_str(),fTitle4.c_str(),tdiv,tstart,tend);
    auto theEffaHistGB = new TH2F((globname+histname3).c_str(),fTitle3.c_str(),2250,0,2250,2250,0,2250);

    auto theetaHistGB = new TH2F((globname+histname5).c_str(),fTitle5.c_str(),182,-90.5,90.5,tdiv,tstart,tend);
    auto thephiHistGB = new TH2F((globname+histname6).c_str(),fTitle6.c_str(),363,-1.5,361.5,tdiv,tstart,tend);

    theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
    theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
    theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());
    theGHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
    theGHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
    theGHistGB->GetZaxis()->SetTitle(fZTitle.c_str());

    std::cout << "Setting up DetIDs." << std::endl;
    std::map<UInt_t,DetIDStruct> DetIDMap;
    SetupDetIDsEB( DetIDMap );
    SetupDetIDsEE( DetIDMap );

    std::cout << "open input files list : " << infilelistname << std::endl;
    string disphotreename("tree/llpgtree");

    std::ifstream infilelist(infilelistname);
    std::string infiles;
    while (std::getline(infilelist,infiles)){
        std::stringstream ss(infiles);
        std::string infilename;
        //std::string califilename;
        ss >> infilename; // >> califilename;
        std::cout << "open input file : " << infilename << std::endl;
        //std::cout << "open input cali : " << califilename << std::endl;

        auto fInFile = TFile::Open((indir+infilename).c_str(), "update");
        fInFile->cd();
        auto fInTree = (TTree*)fInFile->Get(disphotreename.c_str());
		//auto calidir = "/home/t3-ku/jaking/ecaltiming/skimmed_trees/local_chain/";
		//auto calidir = "cali_root_files/";
        //auto fCaliFile = TFile::Open((calidir+califilename).c_str(), "update");
        //std::cout << "fInFile : " << fInFile  << " fInTree : " << fInTree << " fCaliFile : " << fCaliFile << std::endl;

        std::cout << "set branches to get from fInFile : fInTree" << std::endl;

		run = 0;
		resResVecRhID = 0;
		resAmp = 0;
		resE = 0;
		resRtTime = 0;
		resCCtime = 0;
		resTOF = 0;

        fInTree->SetBranchAddress( "run", &run, &b_run );   //!
        fInTree->SetBranchAddress( "resResVecRhID", &resResVecRhID, &b_resResVecRhID );   //!
        fInTree->SetBranchAddress( "resAmp", &resAmp, &b_resAmp);   //!
        fInTree->SetBranchAddress( "resE", &resE, &b_resE);   //!
        fInTree->SetBranchAddress( "resRtTime", &resRtTime, &b_resRtTime);   //!
        fInTree->SetBranchAddress( "resCCtime" , &resCCtime, &b_resCCtime);   //!
        fInTree->SetBranchAddress( "resTOF" , &resTOF, &b_resTOF);   //!

//         get maps from fCaliFile

		vector<TH2F*> calimaps;

/*
        std::cout << "get maps from fCaliFile" << std::endl;
        fCaliFile->cd();
        //string cmbs("AveXtalRtOOTStcPhoIcRecTime");
        //string cmbs("AveXtalRtStcRecTimeE5");
        string cmbs(calimapname);
        //string itcnt("_i95");
        string itcnt("");
        string ebmapstring(cmbs+"EBMap"+itcnt);
        string eberrmapstring(cmbs+"ErrEBMap"+itcnt);
        string epmapstring(cmbs+"EPMap"+itcnt);
        string emmapstring(cmbs+"EMMap"+itcnt);
        auto ebmapic = (TH2F*)fCaliFile->Get(ebmapstring.c_str());
        auto eberrmapic = (TH2F*)fCaliFile->Get(eberrmapstring.c_str());
        auto epmapic = (TH2F*)fCaliFile->Get(epmapstring.c_str());
        auto emmapic = (TH2F*)fCaliFile->Get(emmapstring.c_str());
        std::cout << " Ic hists for " << cmbs << std::endl;
        std::cout << "  " << ebmapstring << " : " << ebmapic << std::endl;
        std::cout << "  " << eberrmapstring << " : " << eberrmapic << std::endl;
        std::cout << "  " << epmapstring << " : " << epmapic << std::endl;
        std::cout << "  " << emmapstring << " : " << emmapic << std::endl;
*/

        std::cout << "Getting calibration values and plotting" << std::endl;

        fInFile->cd();
		//const int bump = 1;
		unsigned int lastRun(0);
        const auto nEntries = fInTree->GetEntries();
        for (auto entry = 0U; entry < nEntries; entry++){
        	// if( entry%int(nEntries*0.1) == 0 ) std::cout << "Proccessed " << entry/nEntries << "\% of " << nEntries << " entries." << std::endl;
	        if( entry%(100000) == 0 ) std::cout << "Mf2d Proccessed " << entry << " of " << nEntries << " entries." << std::endl;        
			//if( entry%(1000*bump) == 0 ) std::cout << "Mf2d Proccessed " << entry << " of " << nEntries << " entries." << std::endl;
            //if( entry%bump != 0 ) continue;

			if(debug) std::cout << " - Start loop " << std::endl;
	        fInFile->cd();

            b_run->GetEntry(entry);   //!
            b_resResVecRhID->GetEntry(entry);   //!
            b_resAmp->GetEntry(entry);   //!
            b_resE->GetEntry(entry);   //!
            b_resRtTime->GetEntry(entry);   //!
            b_resCCtime->GetEntry(entry);   //!
            b_resTOF->GetEntry(entry);   //!

			if( run != lastRun ){
				lastRun = run;
				calimaps = getCaliMaps(run,calimapname);  ///  <<<<<  create this  !!!!!!!!
			}//<<>>if( run != lastRun )					

			if(debug) std::cout << " - Finshed Get Entry " << std::endl;

			auto idinfoL0 = DetIDMap[resResVecRhID[0]];
            auto idinfoL1 = DetIDMap[resResVecRhID[1]];
            auto idinfoG0 = DetIDMap[resResVecRhID[2]];
            auto idinfoG1 = DetIDMap[resResVecRhID[3]];

            auto i1L0 = idinfoL0.i1;
            auto i1L1 = idinfoL1.i1;
            auto i2L0 = idinfoL0.i2;     
            auto i2L1 = idinfoL1.i2;

            auto i1G0 = idinfoG0.i1;
            auto i1G1 = idinfoG1.i1;
            auto i2G0 = idinfoG0.i2;
            auto i2G1 = idinfoG1.i2;

			auto L0EB = idinfoL0.ecal;
            auto L1EB = idinfoL1.ecal;
            auto G0EB = idinfoG0.ecal;
            auto G1EB = idinfoG1.ecal;

	        int bin_offset = 86;
	        int adjust = 0.0;

	        if( calimapname != "none" ){ //and califilename != "none" ){

            	if ( L0EB == ECAL::EB ){
            		lSeedTimeIC0 = calimaps[0]->GetBinContent( i2L0 + bin_offset, i1L0 ) - adjust;
					lSeedTimeICE0 = calimaps[1]->GetBinContent( i2L0 + bin_offset, i1L0 );
            	}else if ( L0EB == ECAL::EP ){
            	    lSeedTimeIC0 = calimaps[2]->GetBinContent( i2L0, i1L0 ) - adjust;
            	}else if (L0EB == ECAL::EM ){
            	    lSeedTimeIC0 = calimaps[3]->GetBinContent( i2L0, i1L0 ) - adjust;
            	}//<<>>if ( L0EB == ECAL::EB )
            
        		if ( L1EB == ECAL::EB ){
            		lSeedTimeIC1 = calimaps[0]->GetBinContent( i2L1 + bin_offset, i1L1 ) - adjust;
                	lSeedTimeICE1 = calimaps[1]->GetBinContent( i2L1 + bin_offset, i1L1 );
            	}else if ( L1EB == ECAL::EP ){
                	lSeedTimeIC1 = calimaps[2]->GetBinContent( i2L1, i1L1 ) - adjust;
            	}else if (L1EB == ECAL::EM ){
                	lSeedTimeIC1 = calimaps[3]->GetBinContent( i2L1, i1L1 ) - adjust;
				}//<<>>if ( L1EB == ECAL::EB )

                if ( G0EB == ECAL::EB ){
                    gSeedTimeIC0 = calimaps[0]->GetBinContent( i2G0 + bin_offset, i1G0 ) - adjust;
                    gSeedTimeICE0 = calimaps[1]->GetBinContent( i2G0 + bin_offset, i1G0 );
                }else if ( G0EB == ECAL::EP ){
                    gSeedTimeIC0 = calimaps[2]->GetBinContent( i2G0, i1G0 ) - adjust;
                }else if (G0EB == ECAL::EM ){
                    gSeedTimeIC0 = calimaps[3]->GetBinContent( i2G0, i1G0 ) - adjust;
                }//<<>>if ( G0EB == ECAL::EB )

                if ( G1EB == ECAL::EB ){
                    gSeedTimeIC1 = calimaps[0]->GetBinContent( i2G1 + bin_offset, i1G1 ) - adjust;
                    gSeedTimeICE1 = calimaps[1]->GetBinContent( i2G1 + bin_offset, i1G1 );
                }else if ( G1EB == ECAL::EP ){
                    gSeedTimeIC1 = calimaps[2]->GetBinContent( i2G1, i1G1 ) - adjust;
                }else if (G1EB == ECAL::EM ){
                    gSeedTimeIC1 = calimaps[3]->GetBinContent( i2G1, i1G1 ) - adjust;
                }//<<>>if ( G1EB == ECAL::EB )

	        }//<<>>if( calimapname != "none" )


//-------------------set for local, repo calcs for global --------------------------------

            if(debug) std::cout << " - Calc 2D Hist" << std::endl;

			double dTOF = resTOF[0]-resTOF[1]; //phoseedTOF_0-phoseedTOF_1;
            double yfill = (resRtTime[0]-lSeedTimeIC0)-(resRtTime[1]-lSeedTimeIC1)+resTOF[0]-resTOF[1];
			if(debug) std::cout << "delta seedTOF: " << to_string(dTOF) 
								<< " 0: " << to_string(phoseedTOF_0) << " 1: "  << to_string(phoseedTOF_1) 
								<< std::endl;
            double effa0 = resAmp[0]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
            double effa1 = resAmp[1]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
            double xfill = (effa0*effa1)/sqrt(pow(effa0,2)+pow(effa1,2));
            double timeErr0 = 0; //phoseedtimeErr_0/25.0;
            double timeErr1 = 0; //phoseedtimeErr_1/25.0;
			double dtserr0 = timeErr0*timeErr0+lSeedTimeICE0*lSeedTimeICE0;
            double dtserr1 = timeErr1*timeErr1+lSeedTimeICE1*lSeedTimeICE1;
			double dterr = sqrt(dtserr0+dtserr1);

		    auto e_cut = (resE[0]>=10)&&(resE[0]<=120)&&(resE[1]>=10)&&(resE[1]<=120);
	        auto eta_cut = (L0EB == ECAL::EB)&&(L1EB == ECAL::EB);

			auto isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, fasle = different
	        auto event_good = e_cut && eta_cut;

//-----------------   set up fills for local same diffrent and global -----------------------

			if(debug) std::cout << " - Fill 2D Hist" << std::endl;
	        if( event_good ){
				theHist->Fill(xfill,yfill);
				theSHist->Fill(xfill,yfill);
				theOthHist->Fill(effa0,yfill);
				for( auto i = 0; i < nRand; i++ ){ 
					//virtual Double_t	Gaus(Double_t mean = 0, Double_t sigma = 1)
					auto gfill = getRandom->Gaus(yfill,dterr);
                  	theGHist->Fill(xfill,gfill);
				}//<<>>for( int i = 0; i < 100; i++ )	
				if(debug) std::cout << " - Fill effA dist Hist" << std::endl;
				//string xbinstr("VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250");
				caliErrBinHist[0]->Fill(xfill);
				timeErrBinHist->Fill(timeErr0); timeErrBinHist->Fill(timeErr1);
               	if(debug) std::cout << " - Fill bin Hists" << std::endl;
				if( xfill < 75 ) caliErrBinHist[1]->Fill(dterr);
				else if( xfill < 100 ) caliErrBinHist[2]->Fill(dterr);
               	else if( xfill < 125 ) caliErrBinHist[3]->Fill(dterr);
               	else if( xfill < 150 ) caliErrBinHist[4]->Fill(dterr);
               	else if( xfill < 175 ) caliErrBinHist[5]->Fill(dterr);
               	else if( xfill < 225 ) caliErrBinHist[6]->Fill(dterr);
               	else if( xfill < 275 ) caliErrBinHist[7]->Fill(dterr);
               	else if( xfill < 325 ) caliErrBinHist[8]->Fill(dterr);
               	else if( xfill < 375 ) caliErrBinHist[9]->Fill(dterr);
               	else if( xfill < 475 ) caliErrBinHist[10]->Fill(dterr);
               	else if( xfill < 600 ) caliErrBinHist[11]->Fill(dterr);
               	else if( xfill < 950 ) caliErrBinHist[12]->Fill(dterr);
               	else if( xfill < 2250 ) caliErrBinHist[13]->Fill(dterr);
				else std::cout << "Over 2250" << std::endl;
            }
            if(debug) std::cout << " - Fill 1D Hist" << std::endl;
            if( event_good ){
				theEffaHist->Fill(effa0,effa1); 
				thetdHist->Fill(phoseedtime_0); thetdHist->Fill(phoseedtime_1);
                thetdcHist->Fill(phoseedtime_0-lSeedTimeIC0); thetdcHist->Fill(phoseedtime_1-lSeedTimeIC1);
				theetaHist->Fill(i20,yfill); theetaHist->Fill(i21,yfill);
                thephiHist->Fill(i10,yfill); thephiHist->Fill(i11,yfill);
			}
			if(debug) std::cout << " - Fill hists done" << std::endl;
        } // end of for loop
	 	delete fInFile;

        //delete fCaliFile;  <<<<<<<<<<<<<<   delete califiles ????????????????

    } // end of while loop

// --------------  process new histos ----------------------------------------------

	if(debug) std::cout << " - Scale 2D Hist" << std::endl;
    scaleHist(theHist,false,fXVarBins,fYVarBins);
    scaleHist(theGHist,false,fXVarBins,fYVarBins);
    NormX(theetaHist);
    NormX(thephiHist);

    fOutFile->cd();
    theHist->Write();
    theGHist->Write();
    theSHist->Write();
    theOthHist->Write();
    theEffaHist->Write();
    thetdHist->Write();
    thetdcHist->Write();
    theetaHist->Write();
    thephiHist->Write();
    timeErrBinHist->Write();
    for( auto ibin = 0; ibin < nMyBins; ibin++ ){ caliErrBinHist[ibin]->Write(); }

    delete theHist;
    delete theGHist;
    delete theSHist;
    delete theOthHist;
    delete theEffaHist;
    delete thetdHist;
    delete thetdcHist;
    delete theetaHist;
    delete thephiHist;
    delete timeErrBinHist;
    for( auto ibin = 0; ibin < nMyBins; ibin++ ){ delete caliErrBinHist[ibin]; }
    delete fOutFile;
	delete getRandom;

    std::cout << "Thats all Folks!" << std::endl;
}


int main ( int argc, char *argv[] ){

        //if( argc != 7 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {
        	auto indir = argv[1];
            auto infilelistname = argv[2];
            auto outfilename = argv[3];
            auto tvarname = argv[4];
			auto calimapname = argv[5];
            auto isd_type = argv[6];
            //auto brun = std::stoi(argv[7]);
            //auto erun = std::stoi(argv[8]);
            //auto leta = std::stoi(argv[9]);
            //auto heta = std::stoi(argv[10]);
      		plot2dResolution( indir, infilelistname, outfilename, tvarname, calimapname, isd_type ); //, brun, erun, leta, heta );
        //}
        return 1;
}
