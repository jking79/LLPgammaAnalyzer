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

#include "llpgana_hist_rebase_v18.hh"

//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSAodSkimmer class ---------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

class KUCMSAodSkimmer : public root_rebase {

    public:

	// constructor varibles
    std::string disphotreename;
    std::string eosdir;
    std::string listdir;

    KUCMSAodSkimmer();
	//~KUCMSAodSkimmer();

	// tchian processing functions
    void kucmsAodSkimmer( std::string indir, std::string infilelist, std::string outfilename, int pct );
    void initHists();
    void getBranches( Long64_t entry );
    bool eventLoop( Long64_t entry );
	void startJobs();
    void endJobs( TTree* fOutTree );
    void setOutputBranches( TTree* fOutTree );

	// object processing & selection
	void processEvntVars();
	void processRechits();
	void processGenParticles();
	void processCalojets();
	void processPhotons();
	void processElectrons();
	void processMuons();
    void processJets();
	void processMet();

	// event processing and selection
	bool eventSelection();

	// histogram and recdhit map varibles & functions
    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

    int nMaps;
    bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
    void makeEBEEMaps( int phoit );
    void makeEBEEMaps( std::vector<unsigned int> rhcol );

    // aod skimmer helper functions & varibles

    std::map<UInt_t,DetIDStruct> DetIDMap;

    int getRhIdx( uInt rhDetID );
    uInt getLeadRhID( vector<uInt> recHitIds );
    float clstrR9( vector<uInt> recHitIds );
    vector<float> getLeadTofRhTime( vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
    vector<float> getRhGrpEigen_sph( vector<float> times, vector<uInt> rechitids );

    // Output Branch variables

    uInt RunNumber;
    uInt nEvents, nSelectedEvents;
	float Met;
	uInt nSelPhotons; 
	vector<uInt> selPhoID;
	vector<float> selPhoTime, selPhoGEgnVal, selPhoEta, selPhoPt, selPhoSMaj, selPhoSMin;
	uint nSelJets;
	vector<uInt> selJetID;
	float JetHt;
	vector<float> selJetPt, selJetEta, selJetTime;


};
