//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSAodSkimmer_code.hh"

// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
				auto indir = "KUCMSNtuple/GMSB_AOD_v1/";
                //auto indir = "KUCMSNtuple/GJETS_AOD_v1/";
                //auto indir = "KUCMSNtuple/WJETS_AOD_v1/";

				auto infilename = "test_GMSB_Files_List.txt";
                //auto infilename = "GMSB_Files_List.txt";
                //auto infilename = "GJets_Files_List.txt";
                //auto infilename = "WJets_Files_List.txt";

				auto outfilename = "_AODSIM_Ntuple_v1_Skim_v4Test.root";

                KUCMSAodSkimmer llpgana;
                llpgana.kucmsAodSkimmer( indir, infilename, outfilename );
    //}
    return 1;

}//<<>>int main ( int argc, char *argv[] )
