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

#include "KUCMSItemManager.hh"
#include "KUCMSBranchManager2.hh"

#include "RestFrames/RestFrames.hh"
#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"

#include "kucmsntuple_root_rebase_v2.hh"


#ifndef KUCMSAodSkimmer_header
#define KUCMSAodSkimmer_header
//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSAodSkimmer class ---------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

using namespace RestFrames;

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36

class KUCMSAodSkimmer : public root_base {

    public:

    KUCMSAodSkimmer();
	~KUCMSAodSkimmer();

	// tchian processing functions
    void kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename );
    void initHists();
    bool eventLoop( Long64_t entry );
	void startJobs();
    void endJobs();
	void fillConfigTree( TTree* fOutTree );
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
	std::vector<float> getRhGrpTimes( std::vector<uInt> rechitids );
    std::vector<float> getRhGrpEnergies( std::vector<uInt> rechitids );
	std::vector<float> getRhGrpEigenFromAngles( std::vector<uInt> rechitids );
    std::vector<float> getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
    std::vector<float> getRhGrpEigen_sph( std::vector<float> times, std::vector<uInt> rechitids );

  	// RestFrames frames and friends

  	LabRecoFrame* LAB;
  	DecayRecoFrame* S;
	DecayRecoFrame*	X2a; 
	DecayRecoFrame* X2b;
  	VisibleRecoFrame* Ja; 
    VisibleRecoFrame* Jb;
  	InvisibleRecoFrame* X1a; 
    InvisibleRecoFrame* X1b;

  	InvisibleGroup* INV;
  	SetMassInvJigsaw* InvM;
  	SetRapidityInvJigsaw* InvEta;
  	MinMassesSqInvJigsaw* InvSplit;

  	CombinatoricGroup* COMB_J;
  	MinMassesSqCombJigsaw* CombSplit_J;

    // config vars
    std::string dataSetKey;
    float xsctn;
    float lam;
	float ctau;
	float mcwgt;
	int mctype;

	// event varibles

    ItemManager<std::vector<float>> geVects;	
    ItemManager<uInt> geCnts;
	ItemManager<float> geVars;	
	uInt nEvents, nSelectedEvents;

    // Output Branch variables

	KUCMSBranchManager selEvtVars;
	KUCMSBranchManager selMet;
	KUCMSBranchManager selPhotons;
	KUCMSBranchManager selJets;
	KUCMSBranchManager selRjrVars;

};

#endif
