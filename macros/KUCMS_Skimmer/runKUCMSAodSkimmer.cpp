//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSAodSkimmer_cc.hh"

// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "ntuple_master_lists/";
                //const string KUCMSAodSkimmer::eosdir = "root://cmseos.fnal.gov//store/user/jaking/";
                const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/";

                auto infilename = "KUCMS_Ntuple_Master_BG_Files_List.txt";
                //auto infilename = "KUCMS_Ntuple_Master_GMSB_Files_List.txt";
                //auto infilename = "KUCMS_Ntuple_Master_JetHT_Files_List.txt";

				auto outfilename = "_Ntuple_v12_LLPgama_Skim_v13.root";

				bool hasGenInfo = true;
                //bool hasGenInfo = false;
                KUCMSAodSkimmer llpgana;
                llpgana.kucmsAodSkimmer( listdir, eosdir, infilename, outfilename, hasGenInfo );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )
