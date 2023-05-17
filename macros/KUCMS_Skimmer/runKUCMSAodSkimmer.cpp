//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include"KUCMSAodSkimmer_code.hh"

// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                //auto indir = "LLPGamma/llpga_GMSB_AOD_v60/"; //argv[1];
                auto indir = "LLPGamma/llpga_GJets_AOD_v60/";

                //auto infilename = "llpgana_mc_AODSIM_GMSB_AOD_v60_Full.txt"; //argv[2];
                //auto infilename = "llpgana_mc_AODSIM_GMSB_AOD_v60_L100.txt";
                auto infilename = "llpgana_mc_AODSIM_GJets_AOD_v60_Full.txt";

                //auto outfilename = "lpgana_mc_test12_tight_max3d4v_sig.root";


                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_GV85Pt30logE_10th_Tight_tEB_nRH15_isSig_v3_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_GV85Pt30logE_10th_Tight_tEB_nRH15_isNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV59_mhV30_GV85Pt30logE_20th_Tight_tEB_nRH15_isCmb_v3_Hists.root";
                //v3 has TAngle 33->27

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_10th_loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_SMM_10th_loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_GEO_10th_loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_SMMGEO_10th_loose_tEB_nRH15_isSig_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_50th_loose_tEB_nRH15_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_SMM_50th_loose_tEB_nRH15_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_GEO_50th_loose_tEB_nRH15_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_SMMGEO_20th_loose_tEB_nRH15_Hists.root";

				auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_skim_v1.root";  // photon only diphton selection
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_L100_AOD_aV60_skim_v1.root"; 

                int pct = 10;
                KUCMSAodSkimmer llpgana;
                llpgana.kucmsAodSkimmer( indir, infilename, outfilename, pct );
    //}
    return 1;

}//<<>>int main ( int argc, char *argv[] )
