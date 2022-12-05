// -*- C++ -*-
//
// Package:    GammaResTool
// Class:      GammaResTool
//
/**\class GammaResTool GammaResTool.cc LLPgammaAnalyzer/plugins/GammaResTool.cc

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
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TGraph.h"
#include "TMathBase.h"

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
#define nHists 256
//const float sol = 29.9792458; // speed of light in cm/ns
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.1415926535 // pie ...  

enum class ECAL {EB, EM, EP, EE, NONE};
#define ecal_config_path "/uscms/home/jaking/nobackup/llpa/CMSSW_10_6_20/src/LLPGamma/GammaResTool/macros/ecal_config/"
//#define ecal_config_path "/LLPGamma/GammaResTool/macros/ecal_config/"

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

class GammaResTool : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

	public:

    	explicit GammaResTool(const edm::ParameterSet&);
      	~GammaResTool();

      	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		detIdMap SetupDetIDs();
		float getPhotonSeedTime( pat::Photon );
		int getRhIdx( uInt rhDetID );

	private:

    	virtual void beginJob() override;
      	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      	virtual void endJob() override;

      	// ----------member data ---------------------------

		// oputput tree
      	TTree *outTree;

		// histograms
		TH1D *hist1d[nHists];
      	TH2D *hist2d[nHists];

		// Flags
		const bool hasGenInfo;

      	// Event
      	unsigned long int event; // technically unsigned long long in Event.h...
      	uInt run, lumi; 

      	// Tracks
      	const edm::InputTag tracksTag;
      	edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      	edm::Handle<std::vector<reco::Track>> tracks_;

        // vertices
        const edm::InputTag verticesTag;
        edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
        edm::Handle<std::vector<reco::Vertex>> vertices_;
        int nVtx;
        float vtxX, vtxY, vtxZ;

      	// electrons
      	const edm::InputTag electronsTag;
      	edm::EDGetTokenT<std::vector<pat::Electron> > electronsToken_;
      	edm::Handle<std::vector<pat::Electron> > electrons_;
      	std::vector<pat::Electron> electrons;
      
      	// RecHits
      	const edm::InputTag recHitsEBTag;
      	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEBToken_;
      	edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEB_;
      	const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > * recHitsEB;

      	const edm::InputTag recHitsEETag;
      	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEEToken_;
      	edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEE_;
      	const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > * recHitsEE;

        detIdMap DetIDMap;
      	int nRecHits;
      	std::vector<float> rhX, rhY, rhZ, rhE, rhtime, rhtimeErr, rhTOF;
      	std::vector<uInt>  rhID;
        std::vector<int>   rhTT;
      	std::vector<bool>  rhisOOT, rhisGS6, rhisGS1;
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

		// skimmer output ::
		// Local : 0,1 ; Z->ee : 2,3
        std::vector<int> 	rhResVecI1;
        std::vector<int> 	rhResVecI2;
        std::vector<int> 	rhResVecEcal;
        std::vector<int> 	rhResVecTT;
        std::vector<float> 	rhResVecE;
        std::vector<float> 	rhResVecadcToGeV;
        std::vector<float> 	rhResVecpedrms12;
        std::vector<float> 	rhResVecTOF;
        std::vector<float> 	rhResVectime;


};//<<>>class GammaResTool : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

//
// Helper functions ( single line function defs, mostly )
//

// #include "TMath.h"
#include <cmath>

void fillTH1( float val, TH1F* hist ){

   	auto nBins = hist->GetNbinsX();
   	auto low = hist->GetBinCenter(1);
   	auto high = hist->GetBinCenter(nBins);
   	if( val < low ) hist->Fill( low );
   	else if ( val > high ) hist->Fill( high );
   	else hist->Fill( val );

}//<<>>void fillTH1( float val, TH1F* hist )

void fillTH1( float val, TH1D* hist ){

   	auto nBins = hist->GetNbinsX();
   	auto low = hist->GetBinCenter(1);
   	auto high = hist->GetBinCenter(nBins);
   	if( val < low ) hist->Fill( low );
   	else if ( val > high ) hist->Fill( high );
   	else hist->Fill( val );

}//<<>>void fillTH1( float val, TH1D* hist )

void normTH2D(TH2D* hist){

    std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;

	const auto nXbins = hist->GetNbinsX();
    const auto nYbins = hist->GetNbinsY();

    for (auto ibinX = 1; ibinX <= nXbins; ibinX++){

        const auto norm = hist->Integral(ibinX,ibinX,1,nYbins);
        if( norm == 0.0 ) continue;
        for (auto ibinY = 1; ibinY <= nYbins; ibinY++){

            // get content/error
            auto content = hist->GetBinContent(ibinX,ibinY);
            auto error   = hist->GetBinError  (ibinX,ibinY);
            // set new contents
            content /= norm;
            error /= norm;
            hist->SetBinContent(ibinX,ibinY,content);
            hist->SetBinError  (ibinX,ibinY,error);

        }//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
	}//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

}//<<>>void NormTH2D(TH2D* hist){

void normTH1D(TH1D* hist){

    std::cout << "Normilizing " << " hist: " << hist->GetName() << std::endl;

    const auto nBins = hist->GetNbinsX();
    const auto norm = hist->Integral();
    for (auto ibinX = 1; ibinX <= nBins; ibinX++){

        if( norm == 0.0 ) continue;
		// get content/error
		auto content = hist->GetBinContent(ibinX);
		auto error   = hist->GetBinError(ibinX);
		// set new contents
        content /= norm;
		error /= norm;
		hist->SetBinContent(ibinX,content);
		hist->SetBinError  (ibinX,error);

	}//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)

}//<<>>void NormTH1D(TH1D* hist)

void profileTH2D(TH2D* nhist, TH1D* prof, TH1D* fithist){

    std::cout << "Profile " << " hist: " << nhist->GetName() << std::endl;

    const auto nXBins = nhist->GetNbinsX();
    //const auto nYBins = nhist->GetNbinsY();
    for (auto ibinX = 1; ibinX <= nXBins; ibinX++){

		auto phist = (TH1F*)nhist->ProjectionY("temp",ibinX,ibinX);
		//double error;
        //auto content = hist->IntegralAndError(ibinX,ibinX,1,nYBins,error);

        auto mean = phist->GetMean();
        //auto mean = 0.f;
        auto stdv = phist->GetStdDev();
		//auto error = phist->GetMeanError();
		auto norm = phist->GetBinContent(phist->GetMaximumBin());
		auto high = mean + 0.2*stdv;
        //auto high = 0.2;
        auto low = mean - 0.2*stdv;
        //auto low = -0.2;
		//std::cout << " - Profile: m " << mean << " s " << stdv << " h " << high << " l " << low << " n " << norm << std::endl;
		if( abs(stdv) > 0.01 && abs(norm) > 1 ){
			auto tmp_form = new TFormula("tmp_formula","[0]*exp(-0.5*((x-[1])/[2])**2)");
			auto tmp_fit  = new TF1("tmp_fit",tmp_form->GetName(),low,high);
            //auto tmp_fit  = new TF1("tmp_fit","crystalball",low,high);
			//auto tmp_fit = new TF1("crystalball", twosided_crystalball_function, low, high, 7);
            //auto tmp_fit  = new TF1("tmp_fit","gaus",high,low);
    		tmp_fit->SetParameter(0,norm); //tmp_fit->SetParLimits(0,norm/2,norm*2);
    		tmp_fit->SetParameter(1,mean); //tmp_fit->SetParLimits(1,-2,2);
    		tmp_fit->SetParameter(2,stdv); //tmp_fit->SetParLimits(2,0,10);
            //tmp_fit->SetParameter(3,0.25);
            //tmp_fit->SetParameter(4,0.25);
            //tmp_fit->SetParameter(5,1);
            //tmp_fit->SetParameter(6,1);
			phist->Fit(tmp_fit->GetName(),"RBQ0");
            //phist->Fit(tmp_fit->GetName(),"R");
			//auto fnorm = tmp_fit->GetParameter(0);
			auto fmean = tmp_fit->GetParameter(1);
            //auto fstd = tmp_fit->GetParameter(2);
			auto error = tmp_fit->GetParError(1);
            //auto fChi2 = tmp_fit->GetChisquare();
            auto fNdf = tmp_fit->GetNDF();
            auto fProb = tmp_fit->GetProb();
			//auto error = fstd/std::sqrt(fnorm);
			//std::cout << " - Profile: fm " << fmean << " fChi2 " << fProb  << " e " << error << " fNdf " << fNdf << std::endl;
	
        	// set new contents
        	if( fNdf > 0 && fProb > 0.05 && error < 1.0 ){
				//auto fChi2Ndf = fChi2/fNdf;
				fithist->SetBinContent( ibinX, fProb );
                fithist->SetBinError( ibinX, 0 );
				prof->SetBinContent( ibinX, fmean );
        		prof->SetBinError( ibinX, error );
			}//<<>>if( fmean < 1 && error < 0.1 )

			//delete tmp_form;
        	delete tmp_fit;
		}//<<>>if( stdv > 0.01 )

    }//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)

}//<<>>void profileTH2D(TH2D* hist, TH1D* prof)

void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

    std::cout << "Threshold Division - " << " hist: " << numi->GetName() << std::endl;

    const auto nXbins = numi->GetNbinsX();
    const auto nYbins = numi->GetNbinsY();

    for (auto ibinX = 1; ibinX <= nXbins; ibinX++){
        for (auto ibinY = 1; ibinY <= nYbins; ibinY++){
 
            // get content/error
            auto ncontent = numi->GetBinContent(ibinX,ibinY);
            auto nerror   = numi->GetBinError  (ibinX,ibinY);
            auto dcontent = denom->GetBinContent(ibinX,ibinY);
            auto derror   = denom->GetBinError  (ibinX,ibinY);
            // set new contents
            auto content(0.0);
            auto error(0.0);
			if( dcontent > thres ){ content = ncontent/dcontent; error = nerror/derror; } 
            numi->SetBinContent(ibinX,ibinY,content);
            numi->SetBinError  (ibinX,ibinY,error);
 
        }//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
    }//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){
 
}//<<>>void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

const float getATan2 ( const float x, const float y){

    if( x == 0 && y == 0) return 6.39;
    else return std::atan2(y,x);

}//<<>> const float getAngle (CFlt x, CFlt y) with atan2

//  --------------------------------------------------------------------------
//  ---------   move all chase functions into class   !!!!!!!!!!!!!!!!!!!!!!!!
//  --------------------------------------------------------------------------

// sort functions

#define CAuto const auto
#define CFlt  const float
#define CDbl  const double
#define CVFlt const vector<float>
#define CVDbl const vector<double>

CAuto sortByPt = [](CAuto & obj1, CAuto & obj2) {return obj1.pt() > obj2.pt();};

// math functions
CAuto sq2			(CFlt x){return x*x;}
CAuto sq2           (CDbl x){return x*x;}
CAuto rad2  		(CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
CAuto hypo  		(CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
CAuto phi   		(CFlt x, CFlt y){return std::atan2(y,x);}
CAuto theta 		(CFlt r, CFlt z){return std::atan2(r,z);}
CAuto eta   		(CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
CAuto effMean   	(CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}
CAuto dPhi			(CFlt x, CFlt y){auto dp(x-y); if( dp > 180 ){dp-=360.0;} else if( dp < -180 ){ dp+=360.0;} return dp;}
CAuto vfsum         (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f);}
CAuto max           (CVFlt x){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

// stats functions
CAuto mean			(CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f)/x.size();}
CAuto mean  		(CVFlt x, CFlt w){return std::accumulate(x.begin(),x.end(),0.0f)/w;}
CAuto mean			(CVFlt x, CVFlt wv){float sum(0.0), wt(0.0); int it(0); for( auto ix : x ){ sum+=ix*wv[it]; wt+=wv[it]; it++; } return sum/wt;}
CAuto wnum			(CFlt it, CFlt w){return (((it-1)*w)/it);}
CAuto stdev			(CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return std::sqrt(sum/(x.size()-1));}	
CAuto var           (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return sum/(x.size()-1);}

CAuto stdev			(CVFlt x, CFlt m, CVFlt wv, CFlt w){
						float sum(0.0); int it(0); 
						for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } 
						return std::sqrt(sum/wnum(it,w));
					}//CAuto stdev

CAuto var           (CVFlt x, CFlt m, CVFlt wv, CFlt w){float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,w);}
CAuto var           (CVFlt x, CFlt m, CVFlt wv){float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,vfsum(wv));}

CAuto cvar          (CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w){
						float sum(0.0); int it(0); 
						for(auto ix : x ){ sum += wv[it]*(ix-mx)*(y[it]-my); it++; } 
						return sum/wnum(it,w);
					}//<<>> CAuto cvar(CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w)

CAuto cvar          (CVFlt x, CFlt mx, CVFlt y, CFlt my){
						float sum(0.0); int it(0); 
						for( auto ix : x ){ sum += (ix-mx)*(y[it]-my); it++; } 
						return sum/(x.size()-1);
					}//CAuto cvar

CAuto rms			(CVFlt x){float sum(0.0); for(auto ix : x ){ sum += sq2(ix); } return std::sqrt(sum/x.size());}
CAuto chisq			(CVFlt x, CFlt m, CFlt v){ float chi(0); for(auto ix : x ){ chi += sq2(ix-m)/v; } return chi; }

CAuto meanPhi		(CVFlt x){
						float sum(0.0);
                        auto maxphi = max(x); 
						for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360); else sum+=ix; } 
						auto rslt = sum/x.size(); 
						if( rslt > 360 ) rslt-=360; 
						return rslt;
					}//<<>> CAuto meanPhi(CVFlt x)

CAuto meanPhi       (CVFlt x, CVFlt wv){
						float wt(0.0), sum(0.0); int it(0); 
						auto maxphi = max(x); 
						for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360)*wv[it]; else sum+=ix*wv[it]; wt+=wv[it]; it++; }
                        auto rslt = sum/wt; 
						if( rslt > 360 ) rslt-=360; 
						return rslt;
					 }//<<>> CAuto meanPhi(CVFlt x, CVFlt wv)

CAuto wsin2         (CVFlt x, CVFlt wv){
						double sum(0.0), wt(0.0); int it(0); 
						for(auto ix : x ){ 
							sum += wv[it]*sq2(sin(ix)); 
							wt += wv[it];
							//std::cout << " ---- wsin2 : " << it << " x: " << ix << " sin^2: " << sq2(sin(ix)) << " sum: " << sum << " wt: " << wt << std::endl;
							it++;
						}//for(auto ix : x )
						return sum/wt;
					}//CAuto wsin2(CVFlt x, CVFlt wv)

CAuto wcos2         (CVFlt x, CVFlt wv){
						double sum(0.0), wt(0.0); int it(0); 
						for(auto ix : x ){ sum += wv[it]*sq2(cos(ix)); wt += wv[it]; it++;} 
						return sum/wt;
					}//CAuto wcos2

CAuto wsincos       (CVFlt x, CVFlt wv){
						double sum(0.0), wt(0.0); int it(0); 
						for(auto ix : x ){ sum += wv[it]*sin(ix)*cos(ix); wt += wv[it]; it++;} 
						return sum/wt;
					}//CAuto wsincos

//
// static data member definitions
//


//-------------------------------------------------------------------------------------------------------------------
