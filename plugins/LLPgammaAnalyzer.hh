// -*- C++ -*-
//
// Package:    LLPgamma/LLPgammaAnalyzer
// Class:      LLPgammaAnalyzer
//
/**\class LLPgammaAnalyzer LLPgammaAnalyzer.cc LLPgamma/LLPgammaAnalyzer/plugins/LLPgammaAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// system include files
#include <memory>

// basic C++ types
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include <utility>
#include <map>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <tuple>
#include <random>
#include <sys/stat.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// HLT + Trigger info
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

// Gen Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// DataFormats
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// DetIds 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

// Ecal RecHits
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"

// Supercluster info
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// EGamma Tools
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

// Geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// ECAL Record info (Laser Constants)
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

// ECAL Record info (Intercalibration Constants)
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

// ECAL Record info (ADCToGeV)
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"

// ECAL Record info (ADCToGeV)
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"

// ECAL Record info (Pedestals)
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

// JECS
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JERs
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// TOOLS
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TFormula.h"
#include "TF1.h"
#include "TTree.h"
#include "Math/PositionVector3D.h"

using namespace std;
using namespace edm;

//
// In class declaration :
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//


using reco::TrackCollection;

//
//  constants, enums and typedefs
//

typedef edm::View<reco::Candidate> CandidateView;
typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> recHitCol;
typedef vector<EcalRecHit> rhGroup;
typedef vector<reco::SuperCluster> scGroup;
typedef vector<reco::CaloCluster> bcGroup;
typedef unsigned int uInt;

#define nEBEEMaps 32
#define nHists 64
//const float sol = 29.9792458; // speed of light in cm/ns
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.1415926535 // pie ...   

enum class ECAL {EB, EM, EP, EE, NONE};
#define ecal_config_path "/home/t3-ku/jaking/ecaltiming/CMSSW_10_2_5/src/Timing/TimingAnalyzer/macros/ecal_config/"
struct DetIDStruct{
   DetIDStruct() {}
   DetIDStruct(const int inI1, const int inI2, const int inTT, const ECAL & inEcal) : i1(inI1), i2(inI2), TT(inTT), ecal(inEcal)  {}
   int i1; // EB: iphi, EE: ix
   int i2; // EB: ieta, EE: iy
   int TT; // trigger tower
   ECAL ecal; // EB, EM, EP
};//>>>>struct DetIDStruct Def
typedef std::map<UInt_t,DetIDStruct> detIdMap;

//
//  Class Declaration
//

class LLPgammaAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {


   public:
      explicit LLPgammaAnalyzer(const edm::ParameterSet&);
      ~LLPgammaAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		detIdMap SetupDetIDs();
      int getPFJetID(const pat::Jet & jet);
		rhGroup getRHGroup( float eta, float phi, float drmin, float minenr );
      rhGroup getRHGroup( const scGroup superClusterGroup, float minenr );
      rhGroup getRHGroup( const reco::CaloCluster basicCluster, float minenr );
      rhGroup getRHGroup( uInt detid );
		rhGroup getRHGroup();
		EcalRecHit getLeadRh( rhGroup recHits );
      vector<float> getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
		vector<float> getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
		vector<float> getTimeDistStats( vector<float> input );
		vector<float> getTimeDistStats( vector<float> times, vector<float> wts );
      vector<float> getTimeDistStats( vector<float> input, rhGroup rechits );
		float getdt( float t1, float t2 );
		void mrgRhGrp( rhGroup & x, rhGroup & y);
		bool reduceRhGrps( vector<rhGroup> & x);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      TTree *outTree;

      TH1D *jetTimeHist, *jetRHTimeHist;
		TH1D *hist1d[nHists];
      TH2D *hist2d[nHists];

		TH2D *ebeeMapSc[nEBEEMaps];
      TH2D *ebeeMapBc[nEBEEMaps];
      TH2D *ebeeMapDr[nEBEEMaps];

      TH2D *ebeeMapT[nEBEEMaps];
      TH2D *ebeeMapE[nEBEEMaps];

		uInt nGoodJetEvents;
		detIdMap DetIDMap;

      // Event
      unsigned long int event; // technically unsigned long long in Event.h...
      uInt run, lumi; 

      // Tracks
      const edm::InputTag tracksTag;
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::Handle<std::vector<reco::Track> > tracks_;

      // PF Candidates
      const edm::InputTag pfcandTag;
      edm::EDGetTokenT<edm::View<reco::Candidate>> pfcand_token_;
      edm::Handle<edm::View<reco::Candidate>> pfcands_;

      // vertices
      const edm::InputTag verticesTag;
      edm::EDGetTokenT<std::vector<reco::Vertex> > verticesToken_;
      edm::Handle<std::vector<reco::Vertex> > vertices_;
      int nVtx;
      float vtxX, vtxY, vtxZ;

      // mets
      const edm::InputTag metsTag;
      edm::EDGetTokenT<std::vector<pat::MET> > metsToken_;
      edm::Handle<std::vector<pat::MET> > mets_;

	   // supercluster
	   const edm::InputTag superClusterCollectionTag;
      edm::EDGetTokenT<reco::SuperClusterCollection> scToken_;
      edm::Handle<reco::SuperClusterCollection> superCluster_;

      const edm::InputTag ootSuperClusterCollectionTag;
      edm::EDGetTokenT<reco::SuperClusterCollection> ootScToken_;
		edm::Handle<reco::SuperClusterCollection> ootSuperCluster_;

		// calocluster
      const edm::InputTag caloClusterTag;
      edm::EDGetTokenT<vector<reco::CaloCluster>> ccToken_;
      edm::Handle<vector<reco::CaloCluster>> caloCluster_;

      // jets
      const edm::InputTag jetsTag;
      edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken_;
      edm::Handle<std::vector<pat::Jet> > jets_;

      uInt nJets;
      std::vector<float> jetE, jetPt, jetPhi, jetEta; 
      std::vector<float> jetMuTime, jetTimeError, jetTimeRMS, jetMedTime, jetCMuTime, jetCMedTime;
      std::vector<float> jetSCMuTime, jetSCMedTime, jetCSCMuTime, jetCSCMedTime, jetCBCMuTime, jetCBCMedTime;
      std::vector<int>   jetID, njetKids, jetKidOfJet, njetSubs, njetRecHits, jetRecHitOfJet;
      std::vector<int>   jetKidPdgID, jetKidCharge, jetKid3Charge, jetPHM, jetELM;
      std::vector<uInt> jetRecHitId;
		std::vector<bool> jetKidLLP;
      std::vector<double> jetKidMass, jetKidVx, jetKidVy, jetKidVz;
      std::vector<float> jetKidE, jetKidPt, jetKidPhi, jetKidEta, jetKidTime, jetKidMedTime;
      std::vector<float> jetNHF, jetNEMF, jetCHF, jetCEMF, jetMUF, jetNHM, jetCHM, jetC, jetPHE, jetPHEF;
      std::vector<float> jetELE, jetELEF, jetMUE, jetCharge;

      // electrons
      const edm::InputTag electronsTag;
      edm::EDGetTokenT<std::vector<pat::Electron> > electronsToken_;
      edm::Handle<std::vector<pat::Electron> > electrons_;
      std::vector<pat::Electron> electrons;
      
      // muons
      const edm::InputTag muonsTag;
      edm::EDGetTokenT<std::vector<pat::Muon> > muonsToken_;
      edm::Handle<std::vector<pat::Muon> > muons_;
      std::vector<pat::Muon> muons;
      
      // RecHits
      const edm::InputTag recHitsEBTag;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEBToken_;
      edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEB_;
      const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > * recHitsEB;

      const edm::InputTag recHitsEETag;
      edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEEToken_;
      edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEE_;
      const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > * recHitsEE;

      int nRecHits;
      std::vector<float> rhX, rhY, rhZ, rhE, rhtime, rhtimeErr, rhTOF;
      std::vector<uInt> rhID;
      std::vector<bool> rhisOOT, rhisGS6, rhisGS1;
      std::vector<float> rhadcToGeV;
      std::vector<float> rhped12, rhped6, rhped1;
      std::vector<float> rhpedrms12, rhpedrms6, rhpedrms1;

      // gedPhotons
      const edm::InputTag gedPhotonsTag;
      edm::EDGetTokenT<std::vector<pat::Photon> > gedPhotonsToken_;
      edm::Handle<std::vector<pat::Photon> > gedPhotons_;

      // ootPhotons
      const edm::InputTag ootPhotonsTag;
      edm::EDGetTokenT<std::vector<pat::Photon> > ootPhotonsToken_;
      edm::Handle<std::vector<pat::Photon> > ootPhotons_;

     
      // geometry CaloSubdetectorGeometry
      edm::ESHandle<CaloGeometry> caloGeo_;
      const CaloSubdetectorGeometry * barrelGeometry;
      const CaloSubdetectorGeometry * endcapGeometry;  

};

//
// Helper functions ( single line function defs, mostly )
//

// sort functions

#define CAuto const auto
#define CFlt  const float
#define CVFlt const vector<float>

CAuto sortByPt = [](CAuto & obj1, CAuto & obj2) {return obj1.pt() > obj2.pt();};

// math functions
CAuto sq2			(CFlt x){return x*x;}
CAuto rad2  		(CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
CAuto hypo  		(CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
CAuto phi   		(CFlt x, CFlt y){return std::atan2(y,x);}
CAuto theta 		(CFlt r, CFlt z){return std::atan2(r,z);}
CAuto eta   		(CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
CAuto effMean   	(CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}

// stats functions
CAuto mean			(CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f)/x.size();}
CAuto mean  		(CVFlt x, CFlt w){return std::accumulate(x.begin(),x.end(),0.0f)/w;}
CAuto mean			(CVFlt x, CVFlt wv){float sum(0.0), wt(0.0); int it(0); for( auto ix : x ){ sum+=ix*wv[it]; wt+=wv[it]; it++; } return sum/wt;}
CAuto stdev			(CVFlt x, float m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return std::sqrt(sum/x.size());}	
CAuto stdev			(CVFlt x, float m, CVFlt wv, CFlt w)
							{float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return std::sqrt(sum/(((it-1)/it)*w));}
CAuto rms			(CVFlt x){float sum(0.0); for(auto ix : x ){ sum += sq2(ix); } return std::sqrt(sum/x.size());}
CAuto sum         (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f);}

// rh group functions
CAuto getRawID		(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.rawId();}
CAuto rhMatch		(const EcalRecHit rhx, const EcalRecHit rhy){ return getRawID(rhx) == getRawID(rhy);}
CAuto dupRhFnd		(const rhGroup x, const rhGroup y){for(CAuto rhx : x ){for(CAuto rhy : y ){if(rhMatch(rhx,rhy)){ return true;}}} return false;}
CAuto isRhGrpEx 	(const rhGroup x){int s=x.size();for( int i=0;i<s;i++){for( int j=i+1;j<s;j++){if(rhMatch(x[i],x[j])) return false;}} return true;}
CAuto getRhGrpEnr	(const rhGroup x){float e(0.0);for( auto ix : x ){e+=ix.energy();} return e;}
CAuto getDupCnt	(const vector<rhGroup> x){int c=0; int s=x.size();for( int a=0;a<s;a++){for( int b=a+1;b<s;b++){if(dupRhFnd(x[a],x[b]))c++;}} return c;}

//
// static data member definitions
//


//-------------------------------------------------------------------------------------------------------------------
