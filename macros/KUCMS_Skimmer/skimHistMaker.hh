//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"

#include "kuSkimTree_v3.h"

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36

//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSSkimmer class -----------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class HistMaker : public kuSkimTree {

	public:

	HistMaker(){};
	//~HistMaker();

	void histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle, int cut );	
	void initHists( std::string htitle );
	//void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
 	void endJobs();	

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];
    int cutselection;

    std::map< std::string, std::map< std::string, float > > configInfo;

	int nMaps;
	bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
    void makeEBEEMaps( int phoit );
    void makeEBEEMaps( std::vector<unsigned int> rhcol );

};

//--------------------------------------------------------------------------------------------------------------------------------
