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


void wc_ku_InterCali_aveRecHit_mini( string indir, string infilelistname, string outfilename ){

    const int  nAlgos = 2; // Mini, MfootCCStc
    //const int  nPhotons = 4;
    const double offset = 0.0;
    const int bin_offset = 86;
    const double ri_ecut = 5.0;

	//const bool debug = true;
    const bool debug = false;

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
             string htitleDistEB( "AveXtal "+algostring[i]+" EBDist EB ");
             IcDistEB[i] = new TH1F(hnameDistEB.c_str(),htitleDistEB.c_str(),320,-4,4);
             IcDistEB[i]->Sumw2();
             string hnameErrDistEB( "AveXtal"+algostring[i]+"EBErrdist");
             string htitleErrDistEB( "AveXtal "+algostring[i]+" EBErrDist EB ");
             IcDistErrEB[i] = new TH1F(hnameErrDistEB.c_str(),htitleErrDistEB.c_str(),200,0,0.2);
             IcDistErrEB[i]->Sumw2();
         }

		 auto IcDistCompEB = new TH2F("ratiovcccomp","Ratio v CC Comp; Rt [ns]; CC [ns]",240,-3,3,240,-3,3);
		 IcDistCompEB->Sumw2();

         std::cout << "Setting up DetIDs." << std::endl;
         std::map<UInt_t,DetIDStruct> DetIDMap;
         SetupDetIDsEB( DetIDMap );
         SetupDetIDsEE( DetIDMap );

    	 std::map<UInt_t,Float_t> sumXtalMiniRecTime;
    	 std::map<UInt_t,Float_t> sumXtalCCStcRecTime;

         std::map<UInt_t,Float_t> sumXtal2MiniRecTime;
         std::map<UInt_t,Float_t> sumXtal2CCStcRecTime;

    	 std::map<UInt_t,UInt_t> numXtalMiniRecTime;
    	 std::map<UInt_t,UInt_t> numXtalCCStcRecTime;

		 std::map<UInt_t,Float_t> xtalMeanRtRecTime;
         std::map<UInt_t,Float_t> xtalMeanCCRecTime;

         // Declaration of leaf types
         UInt_t          run;
         //UInt_t          lumi;
         vector<uInt> 	 *rhCaliID;
         vector<float>   *rhCaliRtTime;
         vector<float>   *rhCaliCCTime;
         //vector<uInt>    *resResVecRhID;
         //vector<float>   *resAmp;
         //vector<float>   *resE;
         //vector<float>   *resRtTime;
         //vector<float>   *resCCtime;
         //vector<float>   *resTOF;

         // List of branches
         TBranch        *b_run;   //!
         //TBranch        *b_lumi;   //!
         TBranch        *b_rhCaliID;   //!
         TBranch        *b_rhCaliRtTime;   //!
         TBranch        *b_rhCaliCCTime;   //!
         //TBranch        *b_resResVecRhID;   //!
         //TBranch        *b_resAmp;   //!
         //TBranch        *b_resE;   //!
         //TBranch        *b_resRtTime;   //!
         //TBranch        *b_resCCtime;   //!
         //TBranch        *b_resTOF;   //!


    	std::ifstream infilelist(infilelistname);
    	std::string infilestr;
    	while (std::getline(infilelist,infilestr)){

        	std::stringstream ss(infilestr);
        	std::string infilename;
        	std::string srunstr;
        	std::string erunstr;
        	ss >> infilename >> srunstr >> erunstr;
        	std::cout << "open input file : " << infilename << std::endl;
        	std::cout << "For Run " << srunstr << " to Run " << erunstr << std::endl;
        	auto srun = std::stoi(srunstr);
        	auto erun = std::stoi(erunstr);



    	 std::ifstream infile(infilename);
    	 std::string instr;
         auto fInTree = new TChain(treename.c_str());
         std::cout << "Adding files to TChain." << std::endl;
         while (std::getline(infile,instr)){
			const std::string eosdir("root://cmseos.fnal.gov//store/user/");
			auto tfilename = eosdir + indir + instr;
         	//auto tfilename = indir + "/" + str;
         	std::cout << "--  adding file: " << tfilename << std::endl;
         	fInTree->Add(tfilename.c_str());
         }//<<>>while (std::getline(infile,str))


         run = 0;
         rhCaliID = 0;
         rhCaliRtTime = 0;
         rhCaliCCTime = 0;

         fInTree->SetBranchAddress("run", &run, &b_run);
         fInTree->SetBranchAddress("rhCaliID", &rhCaliID, &b_rhCaliID);
         fInTree->SetBranchAddress("rhCaliRtTime", &rhCaliRtTime, &b_rhCaliRtTime);
         fInTree->SetBranchAddress("rhCaliCCTime", &rhCaliCCTime, &b_rhCaliCCTime);
     
         // >> calcs  <<
     
         std::cout << "Starting entry loops "<< std::endl;
         auto nEntries = fInTree->GetEntries();
		 if( debug ) nEntries = 100;
         for (Long64_t centry = 0; centry < nEntries; centry++){
			  
     	     if( centry%100000 == 0 or centry == 0) std::cout << "Proccessed " << centry << " of " << nEntries 
																<< " (" << static_cast<float>((10000*centry)/nEntries)/(100) << "%)" << std::endl;
             //if( centry > 100 ) break;
     	     //if( centry%100 != 0 ) continue;   //{  std::cout << "Continuing on : " << centry << endl; continue;	}

			 auto entry = fInTree->LoadTree(centry);

             b_run->GetEntry(entry);   //!
             //b_lumi->GetEntry(entry);   //!
             b_rhCaliID->GetEntry(entry);   //!
             b_rhCaliRtTime->GetEntry(entry);   //!
             b_rhCaliCCTime->GetEntry(entry);   //!
             //b_resResVecRhID->GetEntry(entry);   //!
             //b_resAmp->GetEntry(entry);   //!
             //b_resE->GetEntry(entry);   //!
             //b_resRtTime->GetEntry(entry);   //!
             //b_resCCtime->GetEntry(entry);   //!
             //b_resTOF->GetEntry(entry);   //!

			 if( run < srun || run > erun ) continue;

             const auto nRecHits1 = (*rhCaliID).size(); //(cluster[ipho0])->size();
             if( debug ) std::cout << "Looping over first recHits"  << std::endl;
             for ( auto rh_i = 0U; rh_i < nRecHits1; rh_i++ ){
                  auto id_i = (*rhCaliID)[rh_i];
                  auto Mini_t_i = (*rhCaliRtTime)[rh_i];
                  auto CCStc_t_i = (*rhCaliCCTime)[rh_i];
     	     	  if( debug ) std::cout << "Getting maps " << std::endl;
                  sumXtalMiniRecTime[id_i] += Mini_t_i; 
				  numXtalMiniRecTime[id_i] += 1;
                  sumXtal2MiniRecTime[id_i] += Mini_t_i*Mini_t_i;
                  sumXtalCCStcRecTime[id_i] += CCStc_t_i; 
				  numXtalCCStcRecTime[id_i] += 1;
                  sumXtal2CCStcRecTime[id_i] += CCStc_t_i*CCStc_t_i;

             }//<<>>for (auto i = 0U; i < nRecHits1; i++) // end loop over rechits
             if( debug ) std::cout << "RecHits Loop done "<< std::endl;
         }  //  end entry loop
	delete fInTree;

    } // while (std::getline(infilelist,infiles))

    //double  norm[nAlgos] = { normRtStc, normRtOOTStc };
    std::map<UInt_t,Float_t> *  icmaps[nAlgos] = {&sumXtalMiniRecTime, &sumXtalCCStcRecTime };
    std::map<UInt_t,UInt_t> *  nicmaps[nAlgos] = {&numXtalMiniRecTime, &numXtalCCStcRecTime };
    std::map<UInt_t,Float_t> *  ic2maps[nAlgos] = {&sumXtal2MiniRecTime, &sumXtal2CCStcRecTime };
	std::map<UInt_t,Float_t> *  meanMaps[nAlgos] = {&xtalMeanRtRecTime, &xtalMeanCCRecTime };
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
				(*meanMaps[ai])[it->first] = map_time;
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

	for( std::map<UInt_t,Float_t>::iterator it0=(*meanMaps[0]).begin(); it0!=(*meanMaps[0]).end(); ++it0){
		bool found( false );
		for( std::map<UInt_t,Float_t>::iterator it1=(*meanMaps[1]).begin(); it1!=(*meanMaps[1]).end(); ++it1){
			if( it0->first == it1->first ){ IcDistCompEB->Fill(it0->second,it1->second); found = true; continue; }
		}//<<>>for( std::map<UInt_t,Float_t>::iterator it1=(*IcMapEM[1]).begin(); it1!=(*meanMaps[1]).end(); ++it1)
		if( not found ) IcDistCompEB->Fill(it0->second,-3);
	}//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*IcMapEM[0]).begin(); it!=(*meanMaps[0]).end(); ++it)

    for( auto ai = 0; ai < nAlgos; ai++ ){ (*icmaps[ai]).clear(); (*nicmaps[ai]).clear(); (*meanMaps[ai]).clear(); }

    fOutFile->cd();

    std::cout << "Write AveXtal Rechit Time Maps" << std::endl;
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

	IcDistCompEB->Write();
	delete IcDistCompEB;

    //delete fInTree;
    delete fOutFile;

}

int main ( int argc, char *argv[] ){

        //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {

                //auto indir = "lpcsusylep/jaking/ecalTiming/tt_KUCCRes_126_Test/EGamma/"; //argv[1];
                auto indir = "jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/";
                //auto indir = "jaking/ecalTiming/gammares_tt_kucc_126_v4_flipped/EGamma/";
                //auto indir = "jaking/ecalTiming/gammares_tt_kucc_126_v3/EGamma/";

                ///auto infilelist = "egamma_miniaod_run2018A_316241-316245.txt"; //argv[2];
                ///auto infilelist = "egamma_miniaod_run2022A_352400-358400.txt";
                //auto infilelist = "egamma_run352400-run358400_califilelist.txt";
                //auto infilelist = "egamma_run22C_partial_126_gammares_v2a_califilelist.txt";
                //auto infilelist = "egamma_run18A_316000-316499_126_gammares_v2a_califilelist.txt";
                //auto infilelist = "egamma_run22A_352400-358400_126_gammares_v2a_califilelist.txt";

                //auto infilelist = "egamma_run3_prompt_352319_356513_126_gammares_v2b_califilelist.txt";//1
                //auto infilelist = "egamma_run3_prompt_356514_357289_126_gammares_v2b_califilelist.txt";//2
                //auto infilelist = "egamma_run3_prompt_357290_358883_126_gammares_v2b_califilelist.txt";//3
                //auto infilelist = "egamma_run3_prompt_358884_359420_126_gammares_v2b_califilelist.txt";//4
                //auto infilelist = "egamma_run3_prompt_359421_360089_126_gammares_v2b_califilelist.txt";//5
                //auto infilelist = "egamma_run3_prompt_360090_360981_126_gammares_v2b_califilelist.txt";//6
                //auto infilelist = "egamma_run3_prompt_361417_362522_126_gammares_v2b_califilelist.txt";//8
                //auto infilelist = "egamma_run3_prompt_362523_362760_126_gammares_v2b_califilelist.txt";//9
				auto infilelist = "egamma_run3_prompt_2022C_126_gammares_v2b_califilelist.txt";

                //auto infilelist = "egamma_run3_prompt_357487_357733_126_gammares_v4Flip_califilelist.txt"; //D1
                //auto infilelist = "egamma_run3_iov5_359421_360089_126_gammares_v4Flip_califilelist.txt"; //5
                //auto infilelist = "egamma_run3_prompt_362523_362760_126_gammares_v4Flip_califilelist.txt"; //9
                //auto infilelist = "egamma_run3_iov2_356514_357289_126_gammares_v4Flip_califilelist.txt";//2
                //auto infilelist = "egamma_run3_prompt_357290_358883_126_gammares_v4Flip_califilelist.txt";//3
                //auto infilelist = "egamma_run3_iov4_358884_359420_126_gammares_v4Flip_califilelist.txt";//4
                //auto infilelist = "egamma_run3_22E_359022_360331_126_gammares_v4Flip_califilelist.txt";
                //auto infilelist = "egamma_run3_iov6_360090_360981_126_gammares_v4Flip_califilelist.txt";
                //auto infilelist = "egamma_run3_iov7_360982_361416_126_gammares_v4Flip_califilelist.txt";
                //auto infilelist = "egamma_run3_iov8_361417_362522_126_gammares_v4Flip_califilelist.txt";

                //auto infilelist = "egamma_run3_prompt_357734_358219_126_gammares_v3_califilelist.txt";

                //auto infilelist = "egamma_run3_iov2_356514_357289_126_gammares_v3_califilelist.txt"; 
                //auto infilelist = "egamma_run3_iov3_357290_358883_126_gammares_v3_califilelist.txt";
                //auto infilelist = "egamma_run3_iov4_358884_359420_126_gammares_v3_califilelist.txt";
                //auto infilelist = "egamma_run3_iov5_359421_360089_126_gammares_v3_califilelist.txt";
                //auto infilelist = "egamma_run3_iov6_360090_360981_126_gammares_v3_califilelist.txt"; 
                //auto infilelist = "egamma_run3_iov7_360982_361416_126_gammares_v3_califilelist.txt";
                //auto infilelist = "egamma_run3_iov8_361417_362522_126_gammares_v3_califilelist.txt";
                //auto infilelist = "egamma_run3_iov9_362523_362760_126_gammares_v3_califilelist.txt";

                //auto outfilename = "tt_KUCCRes_126_run2022A_352400-358400_Cali.root"; //argv[3];
                //auto outfilename = "tt_KUCCRes_126_v2a_run2022C_partial_Cali.root"; //argv[3];
                //auto outfilename = "st_RatioRes_126_v2a_run2018A_316000-316499_Cali.root"; //argv[3];
 
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022IOV1_352319_356513_Cali.root"; //argv[3];
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022IOV2_356514_357289_Cali.root"; //argv[3];
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022IOV3_357290_358883_Cali.root"; //argv[3]; 
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022_358884_359420_Cali.root";
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022IOV5_359421_360089_Cali.root"; 
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022_360090_360981_Cali.root";
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022_361417_362522_Cali.root";
                //auto outfilename = "tt_KUCCRes_126_v2b_run3_2022IOV9_362523_362760_Cali.root";
                auto outfilename = "tt_KUCCRes_126_v2b_run3_2022C_355794_357486_Cali.root";

                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022E_359022_360331_Cali.root"; //E
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022_357487_357733_Cali.root";//D1 
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022_357734_358219_Cali.root";

                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022IOV2_356514_357289_Cali2.root"; //2
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022_357290_358883_Cali.root"; //3
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022IOV4_358884_359420_Cali.root"; //4
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022IOV5_359421_360089_Cali.root"; //5
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022IOV6_360090_360981_Cali.root"; //6
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022IOV7_360982_361416_Cali.root"; //7
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022IOV8_361417_362522_Cali.root"; //8
                //auto outfilename = "tt_KUCCRes_126_v4Flip_run3_2022_362523_362760_Cali.root"; //9

                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV2_356514_357289_Cali.root"; //2
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV3_357290_358883_Cali.root"; //3
        		//auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV4_358884_359420_Cali.root"; //4
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV5_359421_360089_Cali.root"; //5
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV6_360090_360981_Cali.root"; //6
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV7_360982_361416_Cali.root"; //7
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV8_361417_362522_Cali.root"; //8
                //auto outfilename = "tt_KUCCRes_126_v3_run3_2022IOV9_362523_362760_Cali.root"; //9

                wc_ku_InterCali_aveRecHit_mini( indir, infilelist, outfilename );
        //}
        return 1;
}

