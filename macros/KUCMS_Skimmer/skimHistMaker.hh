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

#include "kuSkimTree_Mod.h"

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

    void histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle );
	void histMaker( std::string indir, std::string infilelist, std::string outfilename, std::string htitle,int cut,float va,float vb,float vc,float vd );
	void histMaker( std::string indir, std::vector<std::string> infilelists, std::string outfilename, std::string htitle );
    void histMaker( std::string indir, std::string infilelist, std::string outfilename, std::vector<std::vector<std::string>> deflist, std::vector<float> params );

	void initHists( std::string htitle );
    void initHists( std::string htitle, int nhists );
	//void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
	void eventLoop( Long64_t entry, int chist );
 	void endJobs();	

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];
    TH1D *sigHist[n1dHists];
    TH1D *bkgHist[n1dHists];
    TH1D *dataHist[n1dHists];

    int cutselection;
	int preCutNPhotons, preCut30NPhotons, preCut100NPhotons; 
    int postCutNPhotons, postCut30NPhotons, postCut100NPhotons;
    float cutva, cutvb, cutvc, cutvd;
    float sumEvtGenWgt;

    std::map< std::string, std::map< std::string, float > > configInfo;

	std::vector<std::string> bkglist, siglist, datalist, bkgleg, sigleg, dataleg, title, varsel;
	float lumi, maxy, miny, maxr;

	int nMaps, Nsample;
	bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
    void makeEBEEMaps( int phoit );
    void makeEBEEMaps( std::vector<unsigned int> rhcol );

};

//--------------------------------------------------------------------------------------------------------------------------------
