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
             IcDistEB[i] = new TH1F(hnameDistEB.c_str(),htitleDistEB.c_str(),800,-4,4);
             IcDistEB[i]->Sumw2();
             string hnameErrDistEB( "AveXtal"+algostring[i]+"EBErrdist");
             string htitleErrDistEB( "AveXtal "+algostring[i]+" EBErrDist EB ");
             IcDistErrEB[i] = new TH1F(hnameErrDistEB.c_str(),htitleErrDistEB.c_str(),200,0,0.2);
             IcDistErrEB[i]->Sumw2();
         }

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
    for( auto ai = 0; ai < nAlgos; ai++ ){
         for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it){
            const auto & fill_idinfo = DetIDMap[it->first];
            const auto & map_time = (((*icmaps[ai])[it->first])/((*nicmaps[ai])[it->first])) + offset; // - (drift/(icmaps[ai]->size()))) + offset;
            const auto & map_occ = (*nicmaps[ai])[it->first];
            const auto & map_err = sqrt((((*ic2maps[ai])[it->first])/map_occ - map_time*map_time)/map_occ);
		    if( debug ) std::cout << "Fill hist for Algo " << ai << " at " << fill_idinfo.i2 << " " << fill_idinfo.i1 << " with " << map_time << " for iter " << std::endl;
            if( fill_idinfo.ecal == ECAL::EB ){
		   		if( debug ) std::cout << "Fill EB hist for Algo " << ai << " at " << fill_idinfo.i2 << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
            	(IcMapEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
                (IcMapOccEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_occ );
                (IcMapErrEB[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_err );
                (IcDistEB[ai])->Fill(map_time);
                (IcDistErrEB[ai])->Fill(map_err);
            } else if( fill_idinfo.ecal == ECAL::EP ){
                if( debug ) std::cout << "Fill EP hist for Algo " << ai << " at " << fill_idinfo.i2 << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                (IcMapEP[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            } else if( fill_idinfo.ecal == ECAL::EM ){
                if( debug ) std::cout << "Fill EM hist for Algo " << ai << " at " << fill_idinfo.i2 << " " << fill_idinfo.i1 << " with " << map_time << std::endl;
                (IcMapEM[ai])->Fill( fill_idinfo.i2, fill_idinfo.i1, map_time );
            }//<<>>if( fill_idinfo.ecal == ECAL::EB )
         }//<<>>for( std::map<UInt_t,Float_t>::iterator it=(*icmaps[ai]).begin(); it!=(*icmaps[ai]).end(); ++it)
    }//<<>>for( auto ai = 0; ai < nAlgos; ai++ )

    for( auto ai = 0; ai < nAlgos; ai++ ){ (*icmaps[ai]).clear(); (*nicmaps[ai]).clear();}

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

    }

    //delete fInTree;
    delete fOutFile;

}

int main ( int argc, char *argv[] ){

        //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {
                //auto indir = "lpcsusylep/jaking/ecalTiming/tt_KUCCRes_126_Test/EGamma/"; //argv[1];
                auto indir = "jaking/ecalTiming/gammares_tt_kucc_126_v2b/EGamma/";
                ///auto infilelist = "egamma_miniaod_run2018A_316241-316245.txt"; //argv[2];
                ///auto infilelist = "egamma_miniaod_run2022A_352400-358400.txt";
                //auto infilelist = "egamma_run352400-run358400_califilelist.txt";
                //auto infilelist = "egamma_run22C_partial_126_gammares_v2a_califilelist.txt";
                //auto infilelist = "egamma_run18A_316000-316499_126_gammares_v2a_califilelist.txt";
                auto infilelist = "egamma_run22A_352400-358400_126_gammares_v2a_califilelist.txt";
                //auto outfilename = "tt_KUCCRes_126_run2022A_352400-358400_Cali.root"; //argv[3];
                //auto outfilename = "tt_KUCCRes_126_v2a_run2022C_partial_Cali.root"; //argv[3];
                //auto outfilename = "st_RatioRes_126_v2a_run2018A_316000-316499_Cali.root"; //argv[3];
                auto outfilename = "tt_KUCCRes_126_v2a_run2022A_352400-358400_Cali.root"; //argv[3];

                wc_ku_InterCali_aveRecHit_mini( indir, infilelist, outfilename );
        //}
        return 1;
}

