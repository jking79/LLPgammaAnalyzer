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

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36

//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSSkimmer class -----------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class KUCMSSkimmer {

    public:

    void kucmsSkimmer( std::string indir, std::string infilelist, std::string outfilename, int pct );
    void initHists();
    void getBranches( Long64_t entry );
    bool eventLoop( Long64_t entry );
	void startJobs();
    void endJobs( TTree* fOutTree );
    void setOutputBranches( TTree* fOutTree );

	void processEvntVars();
	void processRechits();
	void processGenParticles();
	void processCalojets();
	void processPhotons();
	void processElectrons();
	void processMuons();
    void processJets();
	void processMet();

	bool eventSelection();

    std::map<UInt_t,DetIDStruct> DetIDMap;

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

    int nMaps;
    bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
    void makeEBEEMaps( int phoit );
    void makeEBEEMaps( std::vector<unsigned int> rhcol );

};
