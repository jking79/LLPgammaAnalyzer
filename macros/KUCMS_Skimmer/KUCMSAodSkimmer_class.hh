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
#include "RestFrames/RestFrames.hh"

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36

#include "kucmsntuple_root_rebase_v1.hh"

//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSAodSkimmer class ---------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

using namespace RestFrames;

class KUCMSAodSkimmer : public root_rebase {

    public:

	// constructor varibles
    std::string disphotreename;
    std::string eosdir;
    std::string listdir;

    KUCMSAodSkimmer();
	~KUCMSAodSkimmer();

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
	void processRJR();

	int getPhoQuality( int it );
    int getJetQuality( int it );	

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
    uInt getLeadRhID( std::vector<uInt> recHitIds );
    float clstrR9( std::vector<uInt> recHitIds );
    std::vector<float> getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
    std::vector<float> getRhGrpEigen_sph( std::vector<float> times, std::vector<uInt> rechitids );

  	// RestFrames frames and friends

  	LabRecoFrame*     LAB;
  	DecayRecoFrame*   S;
  	DecayRecoFrame*   X2a;
  	DecayRecoFrame*   X2b;
  	VisibleRecoFrame*   Ja;
  	VisibleRecoFrame*   Jb;
  	InvisibleRecoFrame* X1a;
  	InvisibleRecoFrame* X1b;

  	InvisibleGroup*       INV;
  	SetMassInvJigsaw*     InvM;
  	SetRapidityInvJigsaw* InvEta;
  	MinMassesSqInvJigsaw* InvSplit;

  	CombinatoricGroup*   COMB_J;
  	MinMassesSqCombJigsaw* CombSplit_J;

    // Output Branch variables

    uInt RunNumber;
    uInt nEvents, nSelectedEvents;
	float Met;
	uInt nSelPhotons, leadSelPho, subLeadSelPho; 
	std::vector<uInt> selPhoQuality;
	std::vector<float> selPhoTime, selPhoGEgnVal, selPhoEta, selPhoPhi, selPhoPt, selPhoSMaj, selPhoSMin, selPhoClstrRn;
	uint nSelJets;
	std::vector<uInt> selJetQuality;
	float JetHt;
	std::vector<float> selJetMass, selJetPt, selJetEnergy, selJetEta, selJetPhi, selJetTime;

};
