// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TText.h"
#include "TChain.h"

// STL includes
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <sys/stat.h>

using namespace std;

typedef unsigned int uInt;

template <typename T> std::string to_string(T value)
{  
   //create an output string stream
   std::ostringstream os ;
   
   //throw the value into the string stream
   os << value ;
   
   //convert the string stream into a string and return
   return os.str() ;

//you can also do this
//std::string output;
//os >> output;  //throw whats in the string stream into the string
}

const auto splitString( std::string str, const char* separator ) {

    std::vector < std::string > strings;
    int startIndex(0), endIndex(0);
    for ( int i = 0; i <= str.size(); i++ ){
        if ( str[i] == *separator || i == str.size() ){

            endIndex = i;
            std::string temp;
            temp.append(str, startIndex, endIndex - startIndex);
            strings.push_back(temp);
            startIndex = endIndex + 1;

        }//<<>>if (str[i] == separator || i == str.size())
    }//<<>>for (int i = 0; i <= str.size(); i++)
    return strings;

}//<<>>const auto  splitString(string str, char separator)

enum ECAL {EB, EM, EP, NONE};
//std::string ecal_config_path("/home/t3-ku/jaking/ecaltiming/CMSSW_10_2_5/src/Timing/TimingAnalyzer/macros/ecal_config/");
//std::string ecal_config_path("/uscms/home/jaking/nobackup/ecaltiming/CMSSW_11_3_0_pre6/src/Timing/TimingAnalyzer/macros/ecal_config/");
std::string ecal_config_path("ecal_config/");

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

const auto sq2(const float x){return x*x;}
const auto sq2(const double x){return x*x;}
std::string addstr( std::string current, std::string input ){ return (current+input); }

void fillRatioHist(TH1D* numi, TH1D* denom, TH1D* result ){

    const auto nbins = numi->GetNbinsX();
    for (auto ibin = 0; ibin <= nbins; ibin++){
        auto nc = numi->GetBinContent(ibin);
        auto ncer = numi->GetBinError(ibin);
        auto dc = denom->GetBinContent(ibin);
        auto dcer = denom->GetBinError(ibin);
        auto ratio(0.0);
        auto rerr(0.0);
        if( dc > 20 ){
            ratio = nc/dc;
			rerr = std::sqrt( sq2(ncer/dc)-(sq2(ratio)/dc) );
			//rerr = std::sqrt((sq2(ncer/dc)+sq2((nc/sq2(dc))*dcer))/dc);
        }//<<>>if( dc > 0 )
        result->SetBinContent(ibin,ratio);
        result->SetBinError(ibin,rerr);
    }//<<>>for (auto ibinX = 1; ibinX <= nXbins; ibinX++)

}//<<>>fillRatioHist(TH1F* numi, TH1F* denom, TH1F* result )

void wc_ku_InterCali_aveRecHit_mini( string indir, string infilelistname, string outfilename, float minE ){

    const int  nAlgos = 1; // Mini, MfootCCStc
    //const int  nPhotons = 4;
    const double offset = 0.0;
    const int bin_offset = 86;
	const float minRhEnergy = minE;
    //const float minRhEnergy = 10.0;

	//const bool debug = true;
    const bool debug = false;
	const bool useEnergy = true;

    const string treename("tree/llpgtree");

//    auto fInFile = TFile::Open(infilename.c_str(), "update");
//    auto fInTree = (TTree*)fInFile->Get(disphotreename.c_str());

         std::cout << "Opening Outfile : " << outfilename << std::endl;
         TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
         fOutFile->cd();

         TH2F * IcMapEB[nAlgos];
         TH2F * IcMapErrEB[nAlgos];
         TH2F * IcMapOccEB[nAlgos];
         TH1F * IcDistEB[nAlgos];
         TH1F * IcDistErrEB[nAlgos];
         TH2F * IcMapEP[nAlgos];
         TH2F * IcMapEM[nAlgos];

         TH1F * IcDistMeanEBEta[nAlgos];
         TH1F * IcDistMeanErrEBEta[nAlgos];
         TH1F * IcDistMeanEBPhi[nAlgos];
         TH1F * IcDistMeanErrEBPhi[nAlgos];
         TH1F * IcMapETEta[nAlgos];
         TH1F * IcMapETPhi[nAlgos];
         TH1F * IcMapETErrEta[nAlgos];
         TH1F * IcMapETErrPhi[nAlgos];

         //string algostring[2] = { "Mini", "MfootCCStc" };
         //string algostring[3] = { "Mini", "RtStc", "MfootCCStc" };
         string algostring[2] = { "Ratio", "Kucc" };
         //string algostring[6] = { "RtStc", "RtOOTStc", "WtStc", "WtOOTStc", "MfootStc", "MfootCCStc" };
         for( auto i = 0; i < nAlgos; i++){
             string hnameEB( "AveXtal"+algostring[i]+"RecTimeEBMap");
             string htitleEB( "AveXtal "+algostring[i]+" RecTimeEBMap EB ");
             IcMapEB[i] = new TH2F(hnameEB.c_str(),htitleEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
             IcMapEB[i]->Sumw2();
             string hnameEP( "AveXtal"+algostring[i]+"RecTimeEPMap");
             string htitleEP( "AveXtal "+algostring[i]+" RecTimeEPMap EP ");
             IcMapEP[i] = new TH2F(hnameEP.c_str(),htitleEP.c_str(),100,0.5,100.5,100,0.5,100.5);
             IcMapEP[i]->Sumw2();
             string hnameEM( "AveXtal"+algostring[i]+"RecTimeEMMap");
             string htitleEM( "AveXtal "+algostring[i]+" RecTimeEBMap EM ");
             IcMapEM[i] = new TH2F(hnameEM.c_str(),htitleEM.c_str(),100,0.5,100.5,100,0.5,100.5);
             IcMapEM[i]->Sumw2();
             string hnameOccEB( "AveXtal"+algostring[i]+"OccEBMap");
             string htitleOccEB( "AveXtal "+algostring[i]+" OccEBMap EB ");
             IcMapOccEB[i] = new TH2F(hnameOccEB.c_str(),htitleOccEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
             IcMapOccEB[i]->Sumw2();
             string hnameErrEB( "AveXtal"+algostring[i]+"RecTimeErrEBMap");
             string htitleErrEB( "AveXtal "+algostring[i]+" RecTimeErrEBMap EB ");
             IcMapErrEB[i] = new TH2F(hnameErrEB.c_str(),htitleErrEB.c_str(),171,-85.5,85.5,360,0.5,360.5);
             IcMapErrEB[i]->Sumw2();
             string hnameDistEB( "AveXtal"+algostring[i]+"EBdist");
             string htitleDistEB( "AveXtal "+algostring[i]+" EBMeanDist EB ");
             IcDistEB[i] = new TH1F(hnameDistEB.c_str(),htitleDistEB.c_str(),240,-0.6,0.6);
             IcDistEB[i]->Sumw2();
             string hnameErrDistEB( "AveXtal"+algostring[i]+"EBErrdist");
             string htitleErrDistEB( "AveXtal "+algostring[i]+" EBMeanErrDist EB ");
             IcDistErrEB[i] = new TH1F(hnameErrDistEB.c_str(),htitleErrDistEB.c_str(),200,0,0.2);
             IcDistErrEB[i]->Sumw2();
		}

         	 IcDistMeanEBEta[0] = new TH1F("AveEtaRecTimeEta","AveEtaRecTimeEta;iEta;MeanTime [ns]",171,-85.5,85.5);
             IcDistMeanErrEBEta[0] = new TH1F("AvePhiRecTimeEtaDist","AvePhiRecTimeEtaDist;MeanTime [ns]",120,-0.3,0.3);
             IcDistMeanEBEta[0]->Sumw2();
             IcDistMeanErrEBEta[0]->Sumw2();

         	 IcDistMeanEBPhi[0] = new TH1F("AveEtaRecTimePhi","AveEtaRecTimePhi;iPhi;MeanTime [ns]",360,0.5,360.5);
             IcDistMeanErrEBPhi[0] = new TH1F("AveEtaRecTimePhiDist","AveEtaRecTimePhiDist;MeanTime [ns]",120,-0.3,0.3);
             IcDistMeanEBPhi[0]->Sumw2();
             IcDistMeanErrEBPhi[0]->Sumw2();

         	 IcMapETEta[0] = new TH1F("AveEtaRecMMTimeEta","AveEtaRecMeanMeanTimeEta;iEta;Mean MeanTime [ns]",171,-85.5,85.5);
             IcMapETPhi[0] = new TH1F("AvePhiRecMMTimePhi","AvePhiRecMeanMeanTimePhi;iPhi;Mean MeanTime [ns]",360,0.5,360.5);
             IcMapETPhi[0]->Sumw2();
             IcMapETEta[0]->Sumw2();

             IcMapETErrEta[0] = new TH1F("AveEtaRecMMTimeEtaDist","AveEtaRecMMTimeEtaDist;Mean MeanTime [ns]",120,-0.3,0.3);
             IcMapETErrPhi[0] = new TH1F("AvePhiRecMMTimePhiDist","AvePhiRecMMTimePhiDist;Mean MeanTime [ns]",120,-0.3,0.3);
             IcMapETErrEta[0]->Sumw2();
             IcMapETErrPhi[0]->Sumw2();

		 //auto IcDistCompEB = new TH2F("ratiovcccomp","Ratio v CC Comp; Rt [ns]; CC [ns]",240,-3,3,240,-3,3);
		 //IcDistCompEB->Sumw2();

/*
		std::vector<std::vector<int>> crystals{{22,130},{-19,207},{59,115},{-68,64}};
        std::vector<std::string> fHTitle{"22x130","-19x207","59x115","-68x64"};
		auto nCrys = fHTitle.size();
		TH1D *hist1d136[nCrys], *hist1d140[nCrys], *hist1d137[nCrys], *hist1d138[nCrys], *hist1d139[nCrys], *hist1d143[nCrys], *hist1d144[nCrys];;

		for( int i = 0; i < nCrys; i++ ){
    	//hist1d136[i] = new TH1D(addstr("rhcalCCTimeby100k_",fHTitle[i]).c_str(),addstr(fHTitle[i]," cry CC Time by 100k").c_str(),2000,0,2000);
        //hist1d140[i] = new TH1D(addstr("rhcalRtTimeby100k_",fHTitle[i]).c_str(),addstr(fHTitle[i]," cry Rt Time by 100k").c_str(),2000,0,2000);
    	//hist1d137[i] = new TH1D(addstr("rhcalCntby100k_",fHTitle[i]).c_str(),addstr(fHTitle[i]," cry Cnt by 100k").c_str(),2000,0,2000);
    	hist1d138[i] = new TH1D(addstr("rhcalRtTimeDist_",fHTitle[i]).c_str(),addstr(fHTitle[i]," cry rhcalRtTimeDist").c_str(),500,-25,25);
    	//hist1d139[i] = new TH1D(addstr("rhcalCCTimeDist_",fHTitle[i]).c_str(),addstr(fHTitle[i],"cry rhcalCCTimeDist").c_str(),500,-25,25);
    	//hist1d143[i] = new TH1D(addstr("rhcalCCAveby100k_",fHTitle[i]).c_str(),addstr(fHTitle[i]," cry CC Ave by 100k").c_str(),2000,0,2000);
        hist1d144[i] = new TH1D(addstr("rhcalRtAveby100k_",fHTitle[i]).c_str(),addstr(fHTitle[i]," cry Rt Ave by 100k").c_str(),2000,0,2000);
		//hist1d136[i]->Sumw2();
        //hist1d140[i]->Sumw2();
        //hist1d137[i]->Sumw2();
       // hist1d138[i]->Sumw2();
        //hist1d139[i]->Sumw2();
        //hist1d143[ic]->Sumw2();
        //hist1d144[i]->Sumw2();
		}//<<>>for( int i = 0; i < nCrys; i++ )
*/

         std::cout << "Setting up DetIDs." << std::endl;
         std::map<UInt_t,DetIDStruct> DetIDMap;
         SetupDetIDsEB( DetIDMap );
         SetupDetIDsEE( DetIDMap );

   	 	 std::map<UInt_t,Float_t> sumXtalMiniRecTime;
    	 std::map<int,Float_t> sumXtalEtaRecTime;
         std::map<int,Float_t> sumXtalPhiRecTime;

         std::map<UInt_t,Float_t> sumXtal2MiniRecTime;
         std::map<int,Float_t> sumXtal2EtaRecTime;
         std::map<int,Float_t> sumXtal2PhiRecTime;

    	 std::map<UInt_t,UInt_t> numXtalMiniRecTime;
    	 std::map<int,UInt_t> numXtalEtaRecTime;
         std::map<int,UInt_t> numXtalPhiRecTime;

         std::map<int,Float_t> sumXtalEtaRecMMTime;
         std::map<int,Float_t> sumXtalPhiRecMMTime;
         std::map<int,Float_t> sumXtal2EtaRecMMTime;
         std::map<int,Float_t> sumXtal2PhiRecMMTime;
         std::map<int,UInt_t> numXtalEtaRecMMTime;
         std::map<int,UInt_t> numXtalPhiRecMMTime;


		 //std::map<UInt_t,Float_t> xtalMeanRtRecTime;
         //std::map<int,Float_t> xtalMeanEtaRecTime;
         //std::map<int,Float_t> xtalMeanPhiRecTime;

         // Declaration of leaf types
         uInt          	   run;
		 uInt       	   event;
         //UInt_t          lumi;
         vector<uInt> 	 *rhCaliID;
         vector<float>   *rhCaliRtTime;
         //vector<float>   *rhCaliCCTime;
         //vector<uInt>    *resResVecRhID;
         //vector<float>   *resAmp;
         vector<float>   *rhEnergy;
         vector<float>   *pvTOF;
         vector<float>   *cms0TOF;
         //vector<float>   *resRtTime;
         //vector<float>   *resCCtime;
         //vector<float>   *resTOF;

         vector<float>   *htsecdr4;
         vector<float>   *tspscdr4;
         vector<float>   *erhsecdr4;
         vector<float>   *htoem;

		 vector<int> *scIndex;
		 vector<vector<unsigned int>> *scRhIds;

//Photon_scIndex vector<int>
//SuperCluster_rhIds vector<vector<unsigned int>>
		

		 TBranch        *b_scIndex;
         TBranch        *b_scRhIds;

         // List of branches
         TBranch        *b_run;   //!
         TBranch        *b_event;   //!
         //TBranch        *b_lumi;   //!
         TBranch        *b_rhCaliID;   //!
         TBranch        *b_rhCaliRtTime;   //!
         //TBranch        *b_rhCaliCCTime;   //!
         //TBranch        *b_resResVecRhID;   //!
         //TBranch        *b_resAmp;   //!
         TBranch        *b_rhEnergy;   //!
         TBranch        *b_pvTOF;
         TBranch        *b_cms0TOF;
         //TBranch        *b_resRtTime;   //!
         //TBranch        *b_resCCtime;   //!
         //TBranch        *b_resTOF;   //!

         TBranch        *b_htsecdr4;
         TBranch        *b_tspscdr4;
         TBranch        *b_erhsecdr4;
         TBranch        *b_htoem;

    	std::ifstream infilelist(infilelistname);
    	std::string infilestr;
    	while (std::getline(infilelist,infilestr)){

			//std::cout << "open input string 1: " << infilestr << " ; " << infilestr.size() << std::endl;
			if( infilestr[0] == '#' ) continue;
            if( infilestr.size() == 0 ) continue;
            //std::cout << "open input string 2: " << infilestr << std::endl;
			auto vinstrs = splitString( infilestr, " " );
        	std::string infilename = vinstrs[0];
        	std::string srunstr = vinstrs[1];
        	std::string erunstr = vinstrs[2];
        	std::cout << "open input file : " << infilename << std::endl;
        	std::cout << "For Run " << srunstr << " to Run " << erunstr << std::endl;
        	auto srun = std::stoi(srunstr);
        	auto erun = std::stoi(erunstr);



    	 std::ifstream infile(infilename);
    	 std::string instr;
         auto fInTree = new TChain(treename.c_str());
         std::cout << "Adding files to TChain" << std::endl;
		 //const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
		 const std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");	
         while (std::getline(infile,instr)){
			auto tfilename = eosdir + indir + instr;
         	//auto tfilename = indir + "/" + str;
         	//std::cout << "--  adding file: " << tfilename << std::endl;
         	std::cout << ".";
         	fInTree->Add(tfilename.c_str());
			//break;
         }//<<>>while (std::getline(infile,str))
		std::cout <<std::endl;

         run = 0;
         event = 0;
         rhCaliID = 0;
         rhCaliRtTime = 0;
         //rhCaliCCTime = 0;
   		 rhEnergy = 0;
		 pvTOF = 0;
         cms0TOF = 0;

		 scIndex = 0;
		 scRhIds = 0;
		 htsecdr4 = 0;
		 tspscdr4 = 0;
		 erhsecdr4 = 0;
		 htoem = 0;

         fInTree->SetBranchAddress("Evt_run", &run, &b_run);
         fInTree->SetBranchAddress("Evt_event", &event, &b_event);
         fInTree->SetBranchAddress("ECALRecHit_ID", &rhCaliID, &b_rhCaliID);
         fInTree->SetBranchAddress("ECALRecHit_time", &rhCaliRtTime, &b_rhCaliRtTime);
         fInTree->SetBranchAddress("ECALRecHit_pvTOF", &pvTOF, &b_pvTOF);
         fInTree->SetBranchAddress("ECALRecHit_0TOF", &cms0TOF, &b_cms0TOF);
         //fInTree->SetBranchAddress("rhCaliCCTime", &rhCaliCCTime, &b_rhCaliCCTime);
         //if( useEnergy ) 
		 fInTree->SetBranchAddress("ECALRecHit_energy", &rhEnergy, &b_rhEnergy);
    
         fInTree->SetBranchAddress("Photon_scIndex", &scIndex, &b_scIndex);
         fInTree->SetBranchAddress("SuperCluster_rhIds", &scRhIds, &b_scRhIds);
         fInTree->SetBranchAddress("Photon_hcalTowerSumEtConeDR04", &htsecdr4, &b_htsecdr4);
         fInTree->SetBranchAddress("Photon_trkSumPtSolidConeDR04", &tspscdr4, &b_tspscdr4);
         fInTree->SetBranchAddress("Photon_ecalRHSumEtConeDR04", &erhsecdr4, &b_erhsecdr4);
         fInTree->SetBranchAddress("Photon_hadTowOverEM", &htoem, &b_htoem);

         // >> calcs  <<
     
         std::cout << "Starting entry loops "<< std::endl;
         auto nEntries = fInTree->GetEntries();
		 if( debug ) nEntries = 100;
		 //nEntries = 1000000;
         for (Long64_t centry = 0; centry < nEntries; centry++){
			  
     	     if( centry%100000 == 0 or centry == 0) std::cout << "Proccessed " << centry << " of " << nEntries 
																<< " (" << static_cast<float>((10000*centry)/nEntries)/(100) << "%)" << std::endl;
             //if( centry > 100 ) break;
     	     //if( centry%100 != 0 ) continue;   //{  std::cout << "Continuing on : " << centry << endl; continue;	}

			 auto entry = fInTree->LoadTree(centry);

             if( debug ) std::cout << " - GetEntry : run/event "  << std::endl;
             b_run->GetEntry(entry);   //!
             b_event->GetEntry(entry);   //!
             if( debug ) std::cout << " - GetEntry : id/time "  << std::endl;
             //b_lumi->GetEntry(entry);   //!
             b_rhCaliID->GetEntry(entry);   //!
             b_rhCaliRtTime->GetEntry(entry);   //!
			 b_pvTOF->GetEntry(entry);   //!
			 b_cms0TOF->GetEntry(entry);   //!
             //b_rhCaliCCTime->GetEntry(entry);   //!
             //b_resResVecRhID->GetEntry(entry);   //!
             //b_resAmp->GetEntry(entry);   //!
             //if( useEnergy ) 
             if( debug ) std::cout << " - GetEntry : Energy "  << std::endl;
			 if( useEnergy ) b_rhEnergy->GetEntry(entry);   //!
             //b_resRtTime->GetEntry(entry);   //!
             //b_resCCtime->GetEntry(entry);   //!
             //b_resTOF->GetEntry(entry);   //!
			 b_scIndex->GetEntry(entry);
			 b_scRhIds->GetEntry(entry);

             b_htsecdr4->GetEntry(entry);
             b_tspscdr4->GetEntry(entry);
             b_erhsecdr4->GetEntry(entry);
             b_htoem->GetEntry(entry);

			 //if( run < srun || run > erun ) continue;

         	 //vector<int> scIndex;
         	 //vector<vector<unsigned int>> scRhIds;

             const auto nRecHits1 = rhCaliID->size(); //(cluster[ipho0])->size();
			 const auto nPhos = scIndex->size();
			 for( int pindx = 0; pindx < nPhos; pindx++ ){

				bool passTrkSum = (*tspscdr4)[pindx] < 6.0; 
				bool passsEcalRhSum = (*erhsecdr4)[pindx] < 10.0;
				bool passHOE = (*htoem)[pindx] < 0.02;
				bool passHcalSum = true;
				bool failPhoIso = not ( passHOE && passsEcalRhSum && passTrkSum && passHcalSum );			
				if( failPhoIso ) continue;

				 auto rhlist = (*scRhIds)[pindx];
				 const int nScRh = rhlist.size();
				 for( int rhindx = 0; rhindx < nScRh; rhindx++){
					 uInt scrhid = rhlist[rhindx]; 			
		
		             if( debug ) std::cout << "Looping over first recHits"  << std::endl;
		             for ( auto rh_i = 0U; rh_i < nRecHits1; rh_i++ ){
		
						  auto id_i = (*rhCaliID)[rh_i];
						  if( scrhid == id_i ){
		
						  	auto rhe = useEnergy ? (*rhEnergy)[rh_i] : 999;	
						  	if( rhe < minRhEnergy ) continue;
						  	const auto & fill_idinfo = DetIDMap[id_i];
		                  	auto iEta = fill_idinfo.i2;
		                  	auto iPhi = fill_idinfo.i1;
		                  	auto Mini_t_i = (*rhCaliRtTime)[rh_i] + (*cms0TOF)[rh_i] - (*pvTOF)[rh_i];
		                  	//auto CCStc_t_i = (*rhCaliCCTime)[rh_i];
		     	     	  	if( debug ) std::cout << "Getting maps " << std::endl;
		                  	sumXtalMiniRecTime[id_i] += Mini_t_i; 
						  	numXtalMiniRecTime[id_i] += 1;
		                  	sumXtal2MiniRecTime[id_i] += Mini_t_i*Mini_t_i;
		
						  	if( fill_idinfo.ecal == ECAL::EB && Mini_t_i != 0.0 ){
		                  		sumXtalEtaRecTime[iEta] += Mini_t_i;
		                  		numXtalEtaRecTime[iEta] += 1;
		                  		sumXtal2EtaRecTime[iEta] += Mini_t_i*Mini_t_i;
		
		                  		sumXtalPhiRecTime[iPhi] += Mini_t_i;
		                  		numXtalPhiRecTime[iPhi] += 1;
		                  		sumXtal2PhiRecTime[iPhi] += Mini_t_i*Mini_t_i;
		                  		//sumXtalCCStcRecTime[id_i] += CCStc_t_i; 
						  		//numXtalCCStcRecTime[id_i] += 1;
		                  		//sumXtal2CCStcRecTime[id_i] += CCStc_t_i*CCStc_t_i;
						  	}//<<>>if( fill_idinfo.ecal == ECAL::EB )
						  break;
						  }//if( scrhid == id_i )
		
		             }//<<>>for (auto i = 0U; i < nRecHits1; i++) // end loop over rechits
	
	             }//for( int rhindx = 0; rhindex < nScRh; rhidx++){
             }//for( int pindx = 0; pindx < nPhos; pindx++ ){

             if( debug ) std::cout << "RecHits Loop done "<< std::endl;
         }  //  end entry loop

	delete fInTree;

    } // while (std::getline(infilelist,infiles))

    //double  norm[nAlgos] = { normRtStc, normRtOOTStc };
    std::map<UInt_t,Float_t> *  icmaps[nAlgos] = {&sumXtalMiniRecTime}; //, &sumXtalCCStcRecTime };
    std::map<UInt_t,UInt_t> *  nicmaps[nAlgos] = {&numXtalMiniRecTime}; //, &numXtalCCStcRecTime };
    std::map<UInt_t,Float_t> *  ic2maps[nAlgos] = {&sumXtal2MiniRecTime}; //, &sumXtal2CCStcRecTime };
	//std::map<UInt_t,Float_t> *  meanMaps[nAlgos] = {&xtalMeanRtRecTime}; //, &xtalMeanCCRecTime };
    for( auto ai = 0; ai < nAlgos; ai++ ){
         for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it){
            const auto & fill_idinfo = DetIDMap[it->first];
            const auto & map_time = (((*icmaps[ai])[it->first])/((*nicmaps[ai])[it->first])) + offset; // - (drift/(icmaps[ai]->size()))) + offset;
            const auto & map_occ = (*nicmaps[ai])[it->first];
            const auto & map_err = sqrt((((*ic2maps[ai])[it->first])/map_occ - map_time*map_time)/map_occ);
			if( debug ) std::cout << "Fill hist for Algo " << ai << " at " << fill_idinfo.i2; 
			if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << " for iter " << std::endl;
            if( fill_idinfo.ecal == ECAL::EB ){
		   		if( debug ) std::cout << "Fill EB hist for Algo " << ai << " at " << fill_idinfo.i2 << " "; 
                if( debug ) std::cout << fill_idinfo.i1 << " with " << map_time << std::endl;
            	(IcMapEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
                (IcMapOccEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_occ );
                (IcMapErrEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_err );
                (IcDistEB[ai])->Fill(map_time);
                (IcDistErrEB[ai])->Fill(map_err);
				//(*meanMaps[ai])[it->first] = map_time;
				if( map_time != 0.0 ){
					sumXtalEtaRecMMTime[fill_idinfo.i2] += map_time;
					sumXtalPhiRecMMTime[fill_idinfo.i1] += map_time;
					sumXtal2EtaRecMMTime[fill_idinfo.i2] += map_time*map_time;
					sumXtal2PhiRecMMTime[fill_idinfo.i1] += map_time*map_time;
					numXtalEtaRecMMTime[fill_idinfo.i2] += 1;
					numXtalPhiRecMMTime[fill_idinfo.i1] += 1;
				}//<<>>if( map_time != 0.0 )
            } else if( fill_idinfo.ecal == ECAL::EP ){
                if( debug ) std::cout << "Fill EP hist for Algo " << ai << " at " << fill_idinfo.i2; 
                if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                (IcMapEP[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            } else if( fill_idinfo.ecal == ECAL::EM ){
                if( debug ) std::cout << "Fill EM hist for Algo " << ai << " at " << fill_idinfo.i2;
                if( debug ) std::cout << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                (IcMapEM[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            }//<<>>if( fill_idinfo.ecal == ECAL::EB )
         }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)
    }//<<>>for( auto ai = 0; ai < nAlgos; ai++ )

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecTime.begin(); it!=sumXtalEtaRecTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecTime[iEta]/numXtalEtaRecTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcDistMeanEBEta[0]->SetBinContent(iEta+86,map_time);
           IcDistMeanEBEta[0]->SetBinError(iEta+86,map_err);
           IcDistMeanErrEBEta[0]->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalPhiRecTime.begin(); it!=sumXtalPhiRecTime.end(); ++it){
       const auto iPhi = it->first;
       const auto & map_time = sumXtalPhiRecTime[iPhi]/numXtalPhiRecTime[iPhi]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalPhiRecTime[iPhi];
       const auto & map_err = sqrt((sumXtal2PhiRecTime[iPhi]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcDistMeanEBPhi[0]->SetBinContent(iPhi,map_time);
           IcDistMeanEBPhi[0]->SetBinError(iPhi,map_err);
           IcDistMeanErrEBPhi[0]->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalEtaRecMMTime.begin(); it!=sumXtalEtaRecMMTime.end(); ++it){
       const auto iEta = it->first;
       const auto & map_time = sumXtalEtaRecMMTime[iEta]/numXtalEtaRecMMTime[iEta]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalEtaRecMMTime[iEta];
       const auto & map_err = sqrt((sumXtal2EtaRecMMTime[iEta]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           //IcMapETEta[0]->Fill( iEta, map_time );
           IcMapETErrEta[0]->Fill( map_time );
		   IcMapETEta[0]->SetBinContent(iEta+86,map_time);
		   IcMapETEta[0]->SetBinError(iEta+86,map_err);
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

    for( std::map<int,Float_t>::iterator it=sumXtalPhiRecMMTime.begin(); it!=sumXtalPhiRecMMTime.end(); ++it){
       const auto iPhi = it->first;
       const auto & map_time = sumXtalPhiRecMMTime[iPhi]/numXtalPhiRecMMTime[iPhi]; // - (drift/(icmaps[ai]->size()))) + offset;
       const auto & map_occ = numXtalPhiRecMMTime[iPhi];
       const auto & map_err = sqrt((sumXtal2PhiRecMMTime[iPhi]/map_occ - map_time*map_time)/map_occ);
       //if( fill_idinfo.ecal == ECAL::EB ){
           IcMapETPhi[0]->SetBinContent(iPhi,map_time);
           IcMapETPhi[0]->SetBinError(iPhi,map_err);
           IcMapETErrPhi[0]->Fill( map_time );
       //}//<<>>if( fill_idinfo.ecal == ECAL::EB )
    }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)

/*
	for( std::map<UInt_t,Float_t>::iterator it0=(*meanMaps[0]).begin(); it0!=(*meanMaps[0]).end(); ++it0){
		bool found( false );
		for( std::map<UInt_t,Float_t>::iterator it1=(*meanMaps[1]).begin(); it1!=(*meanMaps[1]).end(); ++it1){
			if( it0->first == it1->first ){ IcDistCompEB->Fill(it0->second,it1->second); found = true; continue; }
		}//<<>>for( std::map<UInt_t,Float_t>::iterator it1=(*IcMapEM[1]).begin(); it1!=(*meanMaps[1]).end(); ++it1)
		if( not found ) IcDistCompEB->Fill(it0->second,-3);
	}//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*IcMapEM[0]).begin(); it!=(*meanMaps[0]).end(); ++it)
*/

    //for( auto ai = 0; ai < nAlgos; ai++ ){ (*icmaps[ai]).clear(); (*nicmaps[ai]).clear(); (*meanMaps[ai]).clear(); }

    fOutFile->cd();

/*
    std::cout << "Write AveXtal Rechit Time Maps" << std::endl;

	for( int i = 0; i < nCrys; i++ ){
	hist1d136[i]->Write(); 
    hist1d140[i]->Write();
	hist1d137[i]->Write();
	hist1d138[i]->Write();
	hist1d139[i]->Write();
	fillRatioHist(hist1d136[i],hist1d137[i],hist1d143[i]);
	hist1d143[i]->Write();
    fillRatioHist(hist1d140[i],hist1d137[i],hist1d144[i]);
    hist1d144[i]->Write();
	}//<<>>for( int i = 0; i < nCrys; i++ )
*/

	IcMapETErrPhi[0]->Write();
	IcMapETPhi[0]->Write();
	IcMapETErrEta[0]->Write();
	IcMapETEta[0]->Write();
	IcDistMeanErrEBPhi[0]->Write();
	IcDistMeanEBPhi[0]->Write();
	IcDistMeanErrEBEta[0]->Write();
	IcDistMeanEBEta[0]->Write();


    for( auto i = 0; i < nAlgos; i++){

         IcMapEB[i]->Write();
         IcMapErrEB[i]->Write();
         IcMapOccEB[i]->Write();
         IcDistEB[i]->Write();
         IcDistErrEB[i]->Write();
         IcMapEP[i]->Write();
         IcMapEM[i]->Write();

         delete IcMapEB[i];
         delete IcMapErrEB[i];
         delete IcMapOccEB[i];
         delete IcDistEB[i];
         delete IcDistErrEB[i];
         delete IcMapEP[i];
         delete IcMapEM[i];

    }//<<>>for( auto i = 0; i < nAlgos; i++)

/*
	for( int i = 0; i < nCrys; i++ ){
    delete hist1d136[i];
    delete hist1d140[i];
    delete hist1d137[i];
    delete hist1d138[i];
    delete hist1d139[i];
    delete hist1d143[i];
    delete hist1d144[i];
	}//<<>>for( int i = 0; i < nCrys; i++ ){
*/
	//IcDistCompEB->Write();
	//delete IcDistCompEB;

    //delete fInTree;
    delete fOutFile;
	std::cout << "Thats All Folkss !!!!!" << std::endl;

}

int main ( int argc, char *argv[] ){

        //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {

				float minE = 5.0;

                //auto indir = "KUCMSNtuple/kucmsntuple_QCD_AOD_v14B/";
                auto indir = "KUCMSNtuple/kucmsntuple_GJETS_AOD_v14B/";
                //auto indir = "KUCMSNtuple/kucmsntuple_GMSB_AOD_v14/";
                //auto indir = "KUCMSNtuple/kucmsntuple_JetHT_Met150_AOD_v14/";

                //auto infilelist = "cali_list_files/KUCMS_QCD_v14_califilelist.txt";
                auto infilelist = "cali_list_files/KUCMS_GJets_v14_califilelist.txt";
                //auto infilelist = "cali_list_files/KUCMS_GMSB_v14_califilelist.txt";
				//auto infilelist = "cali_list_files/KUCMS_JetHT18D_v14_califilelist.txt";

                //auto outfilename = "KUCMS_QCD_v14B3_phorhE5_Cali.root";
                auto outfilename = "KUCMS_GJets_v14B23_phorhE5_Cali.root";
                //auto outfilename = "KUCMS_GMSB_ct200_L350_v14_rhE0_Cali.root";
                //auto outfilename = "KUCMS_JetHT18D_v14_rhE2_Cali.root";

                wc_ku_InterCali_aveRecHit_mini( indir, infilelist, outfilename, minE );
        //}
        return 1;
}

