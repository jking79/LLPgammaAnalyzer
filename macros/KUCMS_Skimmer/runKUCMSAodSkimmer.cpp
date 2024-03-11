//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

//#include "KUCMSAodSkimmer_cc.hh"
//#include "KUCMSAodSkimmer_rh_cc.hh"
#include "KUCMSAodSkimmer_v16_rh_cc.hh"
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "ntuple_master_lists/";
                //const string KUCMSAodSkimmer::eosdir = "root://cmseos.fnal.gov//store/user/jaking/";
                //const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/";
				const std::string eosdir = "";				

                //const std::string infilename = "KUCMS_Ntuple_Master_BG_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_GMSB_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_JetHT_Files_List.txt";
                ////const std::string infilename = "KUCMS_Ntuple_Master_GMSB_Test_Files_List.txt";
                const std::string infilename = "KUCMS_Ntuple_Test_Files_List.txt";

				//const std::string outfilename = "_Ntuple_v14_LLPgama_Skim_v15.root";
				//const std::string outfilename = "_LLPgama_Skim_v15b_rhe1k.root";
                const std::string outfilename = "_LLPgama_Skim_v16_rhe0p2.root";
				//const std::string outfilename = "_LLPgama_Skim_v15_genSigPerfect.root";

				//bool hasGenInfo = true;
                bool hasGenInfo = false;
				//bool genSigPerfect = true;
                bool genSigPerfect = false;
                KUCMSAodSkimmer llpgana;
                llpgana.kucmsAodSkimmer( listdir, eosdir, infilename, outfilename, hasGenInfo, genSigPerfect );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )
