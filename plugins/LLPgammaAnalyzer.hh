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
#include "DataFormats/Common/interface/ValueMap.h"
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

// ROOT
#include "TH1F.h"
#include "TTree.h"
#include "Math/PositionVector3D.h"

using namespace std;
using namespace edm;

//
// class declaration
//
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class LLPgammaAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit LLPgammaAnalyzer(const edm::ParameterSet&);
      ~LLPgammaAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      int GetPFJetID(const pat::Jet & jet);

      // ----------member data ---------------------------
     
      TTree * outTree;
 
      // Event
      unsigned long int event; // technically unsigned long long in Event.h...
      unsigned int run, lumi; 

      // Tracks
      const edm::InputTag tracksTag;
      edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      edm::Handle<std::vector<reco::Track> > tracks_;

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

      // jets
      const edm::InputTag jetsTag;
      edm::EDGetTokenT<std::vector<pat::Jet> > jetsToken_;
      edm::Handle<std::vector<pat::Jet> > jets_;
      std::vector<pat::Jet> jets;
      int nJets;
      std::vector<float> jetE, jetpt, jetphi, jeteta;
      std::vector<int>   jetID;
      std::vector<float> jetNHF, jetNEMF, jetCHF, jetCEMF, jetMUF, jetNHM, jetCHM;

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
      std::vector<unsigned int> rhID;
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

};

//
// constants, enums and typedefs
//

      const auto sortByPt = [](const auto & obj1, const auto & obj2) {return obj1.pt() > obj2.pt();};

//
// static data member definitions
//


//-------------------------------------------------------------------------------------------------------------------
