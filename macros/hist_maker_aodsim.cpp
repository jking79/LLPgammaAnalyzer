// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TText.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFormula.h"
#include "Math/PositionVector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TGraph.h"
#include "TMathBase.h"

// STL includes
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <utility>
#include <algorithm>
#include <sys/stat.h>

//#include "llpgana_hist_base_v2.hh" 
//#include "llpgana_ntuple_base.hh"
//#include "llpgana_hist_base_v10.hh"
//#include "llpgana_hist_rebase_v11.hh"
#include "llpgana_hist_rebase_v17.hh"

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.14159265358979323846 // pie ... 

#define CFlt  const float
#define CDbl  const double
#define CVFlt const vector<float>
#define CVDbl const vector<double>

#define DEBUG false
//#define DEBUG true

typedef unsigned int uInt;

enum ECAL {EB, EM, EP, NONE};

struct DetIDStruct { 
  	DetIDStruct() {}
  	DetIDStruct(const Int_t ni1, const Int_t ni2, const Int_t nTT, const Int_t & necal) : i1(ni1), i2(ni2), TT(nTT), ecal(necal){}
  	Int_t i1; // EB: iphi, EE: ix
  	Int_t i2; // EB: ieta, EE: iy
  	Int_t TT; // trigger tower
  	Int_t ecal; // EB, EM, EP
};//<<>>struct DetIDStruct

//------------------------------------------------------------------------------------------------------------------------------------
//  ------------  helper functions --------------
//------------------------------------------------------------------------------------------------------------------------------------

void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap ){   

    const std::string detIDConfigEB("../test/fullinfo_detids_EB.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);
    
    UInt_t cmsswId, dbID; 
    Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    TString pos;
    
    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM){   
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << 0 << std::endl;
        DetIDMap[cmsswId] = {iphi,ieta,TT25,0};
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }//while (infile >>

}//<<>>void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )

void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap ){

    const std::string detIDConfigEE("../test/fullinfo_detids_EE.txt");
    std::ifstream infile( detIDConfigEE, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    TString EE;

    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant){
        int ec = 1;
        if( side > 0 ) ec = 2;
        //std::cout << "DetID Input Line: " << cmsswId << " " << ix << " "  << iy << " " << ec << std::endl; 
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec};
        //auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )

//------------------------------------------------------------

const float getAngle ( const float x, const float y){

    if( x == 0 && y == 0) return 6.39;
    auto a = std::atan2(y,x);
    if( a < 0 ) a = 2*PI+a;
    return a;

}//<<>> const float getAngle (CFlt x, CFlt y) with atan2

float getdt( float t1, float t2 ){

    auto dt = t1 - t2;
    if( dt == 0.0 ) dt = -5.5;// to be commented out in final version
    if( t1 == 0.0 || t2 == 0.0 ) dt = -5.0;// to be copmmented out in final version
    if( t1 < -28.9 || t2 < -28.9 ) dt = -3.75;

    return dt;
}//<<>>getdt( float t1, float t2 )

// math functions -------------------------

const auto sq2      (CFlt x){return x*x;}
const auto sq2      (CDbl x){return x*x;}
const auto rad2     (CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
const auto hypo     (CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
const auto phi      (CFlt x, CFlt y){return std::atan2(y,x);}
const auto theta    (CFlt r, CFlt z){return std::atan2(r,z);}
const auto eta      (CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
const auto effMean  (CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}
const auto dltIPhi  (CFlt x, CFlt y){auto dp(x-y); if( dp > 180 ){dp-=360.0;} else if( dp < -180 ){ dp+=360.0;} return dp;}
const auto dltPhi   (CFlt x, CFlt y){auto dp(x-y);if(dp>PI) dp-=2*PI; else if(dp<=-PI) dp+=2*PI; return dp;}
const auto max      (CVFlt x){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

const auto deltaR2	(CDbl e0, CDbl e1, CDbl p0, CDbl p1 ){ auto dp(p1-p0); if(dp>PI) dp-=2*PI; else if(dp<=-PI) dp+=2*PI; return sq2(dp)+sq2(e1-e0);}
const auto deltaR	(CDbl e0, CDbl e1, CDbl p0, CDbl p1 ){ return std::sqrt(deltaR2(e0,e1,p0,p1));} 

// Histogram functions ------------------------------

/*
void fillOUHist1F( float val, float low, float high, float div, TH1F * hist ){

    auto step = ((high-low)/div)/2;
    if( val < low ) hist->Fill( low+step );
    else if ( val > high ) hist->Fill( high-step );
    else hist->Fill( val );

}//<<>>void fillOUHist1F( float val, float low, float high, TH1F & hist )
*/

void fillTH1( float val, TH1F *& hist ){

   auto nBins = hist->GetNbinsX();
   auto low = hist->GetBinCenter(1);
   auto high = hist->GetBinCenter(nBins);
   if( val < low ) hist->Fill( low );
   else if ( val > high ) hist->Fill( high );
   else hist->Fill( val );

}//<<>>void fillTH1( float val, TH1F *& hist )

void fillTH1( float val, TH1D *& hist ){

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
            if( dcontent > 0 ){
            	auto content(0.0);
            	auto error(0.0);
            	if( dcontent > thres ){ content = 100+(ncontent/dcontent); error = nerror/derror; }
            	numi->SetBinContent(ibinX,ibinY,content);
            	numi->SetBinError  (ibinX,ibinY,error);
			}//<<>>if( dcontent > 0 ){
        }//<<>>for (auto ibinY = 1; ibinY <= nXbins; ibinY++){
    }//<<>>for (auto ibinX = 1; ibinX <+ nYbins; ibinX++){

}//<<>>void thresDivTH2D(TH2D* numi, TH2D* denom, float thres){

// stats functions --------------------------------------------------------

const auto accum	(CVFlt x){float sum(0.0); for( auto ix : x ){ sum += ix; } return sum; } 
const auto accuminv (CVFlt x){float sum(0.0); for( auto ix : x ){ sum += 1/ix; } return sum; }
const auto mean     (CVFlt x){return accum(x)/x.size();}
const auto mean     (CVFlt x, CFlt w){return accum(x)/w;}
const auto mean     (CVFlt x, CVFlt wv){float sum(0.0), wt(0.0); int it(0); for( auto ix : x ){ sum+=ix*wv[it]; wt+=wv[it]; it++; } return sum/wt;}
const auto wnum     (CFlt it, CFlt w){return (((it-1)*w)/it);}
const auto stdev    (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return std::sqrt(sum/(x.size()-1));}
const auto var      (CVFlt x, CFlt m){float sum(0.0); for( auto ix : x ){ sum += sq2(ix-m); } return sum/(x.size()-1);}

const auto stdev    (CVFlt x, CFlt m, CVFlt wv, CFlt w){
                        float sum(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; }
                        return std::sqrt(sum/wnum(it,w));
                    }//const auto stdev

const auto var      (CVFlt x, CFlt m, CVFlt wv, CFlt w){float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,w);}
const auto var      (CVFlt x, CFlt m, CVFlt wv){float sum(0.0); int it(0); for(auto ix : x ){ sum += wv[it]*sq2(ix-m); it++; } return sum/wnum(it,accum(wv));}

const auto cvar     (CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w){
                        float sum(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*(ix-mx)*(y[it]-my); it++; }
                        return sum/wnum(it,w);
                    }//<<>> const auto cvar(CVFlt x, CFlt mx, CVFlt y, CFlt my, CVFlt wv, CFlt w)

const auto cvar     (CVFlt x, CFlt mx, CVFlt y, CFlt my){
                        float sum(0.0); int it(0);
                        for( auto ix : x ){ sum += (ix-mx)*(y[it]-my); it++; }
                        return sum/(x.size()-1);
                    }//const auto cvar

const auto rms      (CVFlt x){float sum(0.0); for(auto ix : x ){ sum += sq2(ix); } return std::sqrt(sum/x.size());}
const auto chisq    (CVFlt x, CFlt m, CFlt v){ float chi(0); for(auto ix : x ){ chi += sq2(ix-m)/v; } return chi; }

const auto chisqv   (CVFlt x, CFlt m, CVFlt wv, CFlt mv ){  //  !!!!! assumes wv's are inverted + squared && mv is squared
						float sum(0.0), var(0); int it(0); 
						for(auto ix : x ){ sum += sq2(ix-m)/std::sqrt(1/wv[it]+mv); it++; } 
						return sum;
					}//<<>>const auto chisqv 

const auto meanIPhi  (CVFlt x){
                        float sum(0.0);
                        auto maxphi = max(x);
                        for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360); else sum+=ix; }
                        auto rslt = sum/x.size();
                        if( rslt > 360 ) rslt-=360;
                        return rslt;
                    }//<<>> const auto meanIPhi(CVFlt x)

const auto meanIPhi  (CVFlt x, CVFlt wv){
                        float wt(0.0), sum(0.0); int it(0);
                        auto maxphi = max(x);
                        for(auto ix : x ){ if( (maxphi-ix) > 180 ) sum+=(ix+360)*wv[it]; else sum+=ix*wv[it]; wt+=wv[it]; it++; }
                        auto rslt = sum/wt;
                        if( rslt > 360 ) rslt-=360;
                        return rslt;
                     }//<<>> const auto meanIPhi(CVFlt x, CVFlt wv)

const auto wsin2    (CVFlt x, CVFlt wv){
                        double sum(0.0), wt(0.0); int it(0);
                        for(auto ix : x ){
                            sum += wv[it]*sq2(sin(ix));
                            wt += wv[it];
                    //std::cout << " ---- wsin2 : " << it << " x: " << ix << " sin^2: " << sq2(sin(ix)) << " sum: " << sum << " wt: " << wt << std::endl;
                            it++;
                        }//for(auto ix : x )
                        return sum/wt;
                    }//const auto wsin2(CVFlt x, CVFlt wv)


const auto wcos2    (CVFlt x, CVFlt wv){
                        double sum(0.0), wt(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*sq2(cos(ix)); wt += wv[it]; it++;}
                        return sum/wt;
                    }//const auto wcos2

const auto wsincos  (CVFlt x, CVFlt wv){
                        double sum(0.0), wt(0.0); int it(0);
                        for(auto ix : x ){ sum += wv[it]*sin(ix)*cos(ix); wt += wv[it]; it++;}
                        return sum/wt;
                    }//const auto wsincos

vector<float> getRhGrpEigen( vector<float> xs, vector<float> wts ){
//spherical

    vector<float> results;

    auto ts2 = wsin2( xs, wts );
    auto tc2 = wcos2( xs, wts );
    auto tsc = wsincos( xs, wts );
    double array[] = { ts2, tsc, tsc, tc2 };
    TMatrixDSym mat(2,array);
    TMatrixDSymEigen eigen(mat);
    const TVectorD& eginVal = eigen.GetEigenValues();
    const TMatrixD& eginVec = eigen.GetEigenVectors();
    TVectorD evaules(eginVal);

    int index(1);
    if( (eginVal(0) >= eginVal(1)) ){
        index = 0;
        //if( (eginVal(0) == eginVal(1)) ) std::cout << " -- rhGrp is Spherical" << std::endl;
    }//<<>>if( (eginVal(0) >= eginVal(1)) )
    //std::cout << "Largest Eigin Value Found" << std::endl;

    auto ex = eginVec(index,0);
    auto ey = eginVec(index,1);
    auto ev = eginVal(index);
    //std::cout << "Eigin Angles Found" << std::endl;

    results.push_back(ex);
    results.push_back(ey);
    results.push_back(ev);

    return results;
}//<<>>vector<float> getRhGrpEigen2D( vector<float> xs, vector<float> wts )

//--------------------------------------------------------------------------------------------------------------------------------------
// Make hist class -----------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

class makehists : llpgana_hist_rebase {

	public:

	void llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename, int pct );	
	void initHists();
	void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
 	void endJobs();	

	int getRhIdx( uInt rhDetID );
	uInt getLeadRhID( vector<uInt> recHitIds );
	float clstrR9( vector<uInt> recHitIds );
	vector<float> getLeadTofRhTime( vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
	vector<float> getRhGrpEigen_sph( vector<float> times, vector<uInt> rechitids, TH2D* hist2d73, TH2D* hist2d74, TH2D* hist2d75,
                                    	TH2D* hist2d76, TH2D* hist2d86, TH2D* hist2d87 );

    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

	int nMaps;
	bool fMap[nEBEEMaps];
    std::map<UInt_t,DetIDStruct> DetIDMap;
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];

};

//------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------
//// Make hist class ----------------------------------------------------------------------------------------------------------------
////-----------------------------------------------------------------------------------------------------------------------------

//makehists::makehists(){}

//makehists::~makehists(){}


int makehists::getRhIdx( uInt rhDetID ){

	//b_rhID->GetEntry(entry);
	for( int idx = 0; idx < rhID->size(); idx++ ){ if( rhDetID == (*rhID)[idx] ) return idx; }
	//std::cout << " -- !! no rhDetID to (*rhID)[idx] match !! ---------------------- " << std::endl;
	return -1;

}//<<>>int makehists::getRhIdx( int rhDetID )

uInt makehists::getLeadRhID( vector<uInt> recHitIds ){

    uInt result;
    float enr(0.0);
	//b_rhEnergy->GetEntry(entry);
    for( auto id : recHitIds ){
        auto rhenr = (*rhEnergy)[getRhIdx(id)];
        if( rhenr > enr ){ enr = rhenr; result = id; }
    }//<<>>for (const auto recHit : recHits )

    return result;

}//>>>>EcalRecHit LLPgammaAnalyzer_AOD::getLeadRh( rhGroup recHits

float makehists::clstrR9( vector<uInt> recHitIds ){

	auto leadRhID = getLeadRhID( recHitIds );
	auto leadRhEn = (*rhEnergy)[getRhIdx(leadRhID)];
	float sumRhEn(0);
	for ( auto id : recHitIds ){ sumRhEn +=  (*rhEnergy)[getRhIdx(id)]; }
	return sumRhEn > 0 ? leadRhEn/sumRhEn  : 1.2;

}//<<>>float makehists::clstrR9( vector<uInt> recHitIds )

vector<float> makehists::getLeadTofRhTime( vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ ){

	//b_rhPosX->GetEntry(entry);
	//b_rhPosY->GetEntry(entry);
	//b_rhPosZ->GetEntry(entry);
	//b_rhTime->GetEntry(entry);

    vector<float> result;
    if( recHitIds.size() < 1 ){ result.push_back(-99); return result; }
    auto lrhid = getLeadRhID(recHitIds);
	auto lrhidx = getRhIdx(lrhid);
	auto X = (*rhPosX)[lrhidx];
	auto Y = (*rhPosY)[lrhidx];
	auto Z = (*rhPosZ)[lrhidx];
    const auto d_rh = hypo( X, Y, Z);
    const auto d_pv = hypo( X-vtxX, Y-vtxY, Z-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    for( int idx = 0; idx < rhTime->size(); idx++ ){result.push_back((*rhTime)[idx]-tof);}
    return result;

}//>>>>vector<float>  LLPgammaAnalyzer_AOD::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

vector<float> makehists::getRhGrpEigen_sph( vector<float> times, vector<uInt> rechitids, TH2D* hist2d73, TH2D* hist2d74, TH2D* hist2d75, 
												TH2D* hist2d76, TH2D* hist2d86, TH2D* hist2d87 ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);
    //float twoPI(2*PI);

    vector<float> egwts;
    vector<float> etas;
    vector<float> phis;
    vector<float> xs;
    vector<float> ys;
    vector<float> zs;
    vector<float> ebtimes;
    vector<float> angles;
    vector<float> zcAngles;
    //vector<float> energies;
    vector<float> reswts;
	vector<float> tvar;
    //float sumenergy = 0;
    vector<float> emepReturn(8,-9);

	// --------- prepare inputs for eigan calcs --------------------------

    auto nRecHits = rechitids.size();
	if( nRecHits < 5 ) return emepReturn;
    for( uInt it(0); it < nRecHits; it++ ){

        const auto rhIDX = getRhIdx(rechitids[it]);
		auto idinfo = DetIDMap[rechitids[it]];
		//std::cout << "In getRhGrpEigen_sph w/ idx : " << rhIDX << std::endl;
		if( rhIDX == -1 ){ return emepReturn; std::cout << " -- Bad idx !!!!! -- In getRhGrpEigen_sph ---- " << std::endl; }
        if( idinfo.ecal == ECAL::EB ){
            const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
            etas.push_back(rhEtaPos);
            const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
            phis.push_back(rhPhiPos);
            const auto rhXPos = (*rhPosX)[rhIDX];
            xs.push_back(rhXPos);
            const auto rhYPos = (*rhPosY)[rhIDX];
            ys.push_back(rhYPos);
            const auto rhZPos = (*rhPosZ)[rhIDX];
            zs.push_back(rhZPos);
            ebtimes.push_back(times[it]);
            auto rhenergy = (*rhEnergy)[rhIDX];
            auto resolution = sq2(N/rhenergy)+2*C*C;
            reswts.push_back(1/resolution);
			tvar.push_back(resolution);
            //std::cout << "In getRhGrpEigen_sph w/ rheta " << rhEtaPos << " : rhphi " << rhPhiPos << " : rht " << times[it] << std::endl;      
        } else { return emepReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
        }//<<>>if( idinfo.ecal == ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    auto meta = mean(etas,reswts);
    auto mphi = meanIPhi(phis,reswts);
    auto mtime = mean(ebtimes,reswts);
    float totRes = accum(reswts);
	auto mtvar = 1/totRes;

	auto mx = mean(xs,reswts);
    auto my = mean(ys,reswts);
    auto mz = mean(zs,reswts);
	auto mr = hypo(mx,my);
	auto ma = std::atan2(my,mx);

    vector<float> letas;
    vector<float> lphis;
    vector<float> lzs;
    vector<float> lcs;
    vector<float> lts;

    //std::cout << "In getRhGrpEigen_sph w/ meta " << meta << " : mphi " << mphi << " : mt " << mtime << std::endl;

    for( uInt it(0); it < ebtimes.size(); it++ ){

        float leta = etas[it]-meta;
		letas.push_back(leta);
        float lphi = dltIPhi(phis[it],mphi);
        lphis.push_back(lphi);
        if( leta == 0 && lphi == 0 ) continue;
        float angle = getAngle( leta, lphi );
        angles.push_back(angle);

		float lz = zs[it]-mz;
        lzs.push_back(lz);
		float lc = mr*(std::atan2(ys[it],xs[it])-ma);
        lcs.push_back(lc);
		//float zcAngle = getAngle( lz, lc );
        float zcAngle = std::atan2(lc,lz);
		zcAngles.push_back(zcAngle);

        float ltim = ebtimes[it]-mtime;
        lts.push_back(ltim);
        egwts.push_back(ltim*ltim*reswts[it]);

    }//<<>>for( uInt it(0); it < ebtimes.size(); it++ )

	// --------------  get eigan values and vectors ---------------------------------------

    auto epEigens =  getRhGrpEigen( angles, egwts );//0 x, 1 y, 2 values
    auto eigens =  getRhGrpEigen( zcAngles, egwts );//0 x, 1 y, 2 values
    auto geoeigens =  getRhGrpEigen( zcAngles, reswts );//0 x, 1 y, 2 values

	// --------------  get eigan vector angles ------------------------------------------ 

    auto rotangle = getAngle(eigens[0], eigens[1]);
    //float eignsin = std::sin(rotangle);
    float eignsin = eigens[1];
    //float eigncos = std::cos(rotangle);
    float eigncos = eigens[0];

    auto orgrota = rotangle;
    auto orgsin = eignsin;
    auto orgcos = eigncos;

    // ----------------------
    // aligning slope along eigan vector to point same direction
    // ----------------------
    float ltsum(0.0);
    int slopeCorr(1);
	float egflip(1);
    for( uInt it(0); it < egwts.size(); it++ ){

        //float leta = etas[it] - meta;
        //float lphi = dltIPhi(phis[it],mphi);
        //float lz = zs[it]-mz;
        //float lc = mr*(std::atan2(ys[it],xs[it])-ma);
        //float ltime = ebtimes[it]-mtime;
        auto epxcor = orgcos*(letas[it]) - orgsin*(lphis[it]);
        auto xcor = orgcos*(lzs[it]) - orgsin*(lcs[it]);
        //if( xcor > 0 ) ltsum += ltime*reswts[it];
        ltsum += lts[it]*(reswts[it])/xcor;
        //ltsum += lts[it]/xcor;

    }//<<>>for( uInt it(0); it < egwts.size(); it++ )

    //if( ltsum/totRes < 0 ){
    if( ltsum < 0 ){

		egflip = -1;
        rotangle = getAngle(-1*eigens[0], -1*eigens[1]);
        //eignsin = std::sin(rotangle);
        eignsin = -1*eigens[1];
        //eigncos = std::cos(rotangle);
        eigncos = -1*eigens[0];
        //slopeCorr = -1;

    }//if( ebp && ebn && ltsum < 0 )


    // -----------------------------------------
    // finding nemo ( slope )
    // -----------------------------------------

    auto nWts = reswts.size();

    vector<float> xs1;
    vector<float> xs2;
    vector<float> xsgeo;
    vector<float> slvars1;
    vector<float> slvars2;
    vector<float> slvarsgeo;
    float xsum1(0.0);
    float xsum2(0.0);
    float xsumgeo(0.0);
    //float adj(2.2);
	float adj(1.0);

    auto dsxcor = eigncos*(2.2) - eignsin*(2.2);
    auto dxcor = orgcos*(2.2) - orgsin*(2.2);
    auto dgxcor = (geoeigens[0])*(2.2) - (geoeigens[1])*(2.2);
    auto sxcorvar = sq2(dsxcor)/12;
    auto xcorvar = sq2(dxcor)/12;
    auto gxcorvar = sq2(dgxcor)/12;

	// for pairs method
    vector<float> plzs;
    vector<float> plcs;
    vector<float> plws;
    vector<float> plgzs;


    for( uInt it(0); it < nWts; it++ ){

        auto epsxcor = eigncos*(letas[it]*adj) - eignsin*(lphis[it]*adj);
        auto epsycor = eignsin*(letas[it]*adj) + eigncos*(lphis[it]*adj);
        auto epxcor = (orgcos*(letas[it]*adj) - orgsin*(lphis[it]*adj));
        auto epycor = orgsin*(letas[it]*adj) + orgcos*(lphis[it]*adj);

        auto sxcor = eigncos*(lzs[it]*adj) - eignsin*(lcs[it]*adj);
        auto sycor = eignsin*(lzs[it]*adj) + eigncos*(lcs[it]*adj);
        auto xcor = (orgcos*(lzs[it]*adj) - orgsin*(lcs[it]*adj));
        auto ycor = orgsin*(lzs[it]*adj) + orgcos*(lcs[it]*adj);
        auto gxcor = ((geoeigens[0])*(lzs[it]*adj) - (geoeigens[1])*(lcs[it]*adj));
        auto gycor = (geoeigens[1])*(lzs[it]*adj) + (geoeigens[0])*(lcs[it]*adj);

		plzs.push_back(xcor);
		plcs.push_back(ycor);
        plgzs.push_back(gxcor);
		//plws.push_back(
		//plts.push_back(ebtimes[it1]);

        if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << letas[it] << " : lphi " << lphis[it]
                                << " : xcor " << xcor << " : ycor " << ycor << " : dt " << egwts[it] << std::endl;
        if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << letas[it] << " : lphi " << lphis[it]
                                << " : sxcor " << sxcor << " : sycor " << sycor << " : dt " << egwts[it] << std::endl;

        auto fill = lts[it]*reswts[it];
        hist2d73->Fill(sxcor,sycor,fill);
        hist2d74->Fill(sxcor,sycor,reswts[it]);
        hist2d75->Fill(lzs[it],lcs[it],fill);
        hist2d76->Fill(lzs[it],lcs[it],reswts[it]);
        hist2d86->Fill(sycor,lts[it],reswts[it]);
        hist2d87->Fill(sxcor,lts[it],reswts[it]);

        auto sl1 = (lts[it])/(sxcor);//*slopeCorr;
        auto sl2 = (lts[it])/(xcor);//*slopeCorr;
        auto slg = (lts[it])/(gxcor);//*slopeCorr;
        xs1.push_back(sl1);
        xs2.push_back(sl2);
        xsgeo.push_back(slg);

		slvars1.push_back(1/((tvar[it]+mtvar+sq2(sl1)*sxcorvar*(1.0+(1.0/nWts)))/sq2(sxcor)));
        slvars2.push_back(1/((tvar[it]+mtvar+sq2(sl2)*xcorvar*(1.0+(1.0/nWts)))/sq2(xcor)));
        slvarsgeo.push_back(1/((tvar[it]+mtvar+sq2(slg)*gxcorvar*(1.0+(1.0/nWts)))/sq2(gxcor)));
		xsum1 += sl1*slvars1[it];
        xsum2 += sl2*slvars2[it];
        xsumgeo += slg*slvarsgeo[it];

    }//<<>>for( uInt it(0); it < wts.size(); it++ )

//===================================================================

	vector<float> pzs1;
    vector<float> plrs;
    float plsum(0.0);
    vector<float> pgzs1;
    vector<float> plgrs;
    float plgsum(0.0);

    //=============================================================
    // pairs method set up
    //============================================================

    for( uInt it1(0); it1 < plzs.size(); it1++ ){
        for( uInt it2(it1); it2 < plzs.size(); it2++ ){

            auto plz = plzs[it1]-plzs[it2];
			if( plz < 2.2 ) continue;
        	//auto plc = plcs[it1]-plcs[it2]; 
            auto pltime = ebtimes[it1]-ebtimes[it2];
 
			auto psl = pltime/plz;
			pzs1.push_back(psl);
			
			auto plr1 = 1/((tvar[it1]+mtvar+sq2(psl)*sxcorvar*(1.0+(1.0/nWts)))/sq2(plzs[it1]));	
            auto plr2 = 1/((tvar[it2]+mtvar+sq2(psl)*sxcorvar*(1.0+(1.0/nWts)))/sq2(plzs[it2]));
			auto gaplr = std::sqrt(plr1*plr2);

			plrs.push_back( gaplr );
			plsum += psl*gaplr;

           //egwts.push_back(ltim*ltim*reswts[it]);i

        }//<<>>for( uInt it2(it1); it2 < ebtimes.size(); it2++ )
    }//<<>>for( uInt it1(0); it1 < ebtimes.size(); it1++ )

    //============================================================
    // piars geo method set up
    //============================================================

    for( uInt it1(0); it1 < plgzs.size(); it1++ ){
        for( uInt it2(it1); it2 < plgzs.size(); it2++ ){

            auto plz = plgzs[it1]-plgzs[it2];
            if( plz < 2.2 ) continue;
            //auto plc = plcs[it1]-plcs[it2]; 
            auto pltime = ebtimes[it1]-ebtimes[it2];

            auto psl = pltime/plz;
            pgzs1.push_back(psl);

            auto plr1 = 1/((tvar[it1]+mtvar+sq2(psl)*sxcorvar*(1.0+(1.0/nWts)))/sq2(plgzs[it1]));
            auto plr2 = 1/((tvar[it2]+mtvar+sq2(psl)*sxcorvar*(1.0+(1.0/nWts)))/sq2(plgzs[it2]));
            auto gaplr = std::sqrt(plr1*plr2);

            plgrs.push_back( gaplr );
            plgsum += psl*gaplr;

           //egwts.push_back(ltim*ltim*reswts[it]);i

        }//<<>>for( uInt it2(it1); it2 < ebtimes.size(); it2++ )
    }//<<>>for( uInt it1(0); it1 < ebtimes.size(); it1++ )

//--------------------------------------------------------------------

	// -----------   compute final outputs --------------------------------

    auto ebside = (mz > 0)?1:-1;
    auto taflip = ((eigens[0] > 0)?1:-1)*ebside;

    auto nXSum = xs1.size();
	auto totSloRes1 = accum(slvars1);
    //auto slope1 = xsum1/totRes;
    auto slope1 = xsum1/totSloRes1;
    auto slope1err = std::sqrt(1/totSloRes1);
    //auto varsl = var(xs1,slope1,reswts,totRes);
    auto varsl = var(xs1,slope1,slvars1,totSloRes1);
    auto chi1 = chisq(xs1,slope1,varsl);
	auto slchi2v1 = chisqv(xs1,slope1,slvars1,varsl);//?????????????????
    //chisqr probs
    auto chi2pf1 = 1 - TMath::Prob(chi1, nWts);

    //auto slope2 = xsum2/totRes;
    auto totSloRes2 = accum(slvars2);
    auto slope2 = taflip*xsum2/totSloRes2;
    auto slope2err = std::sqrt(1/totSloRes2);
    //auto vars2 = var(xs2,slope2,reswts,totRes);
    auto vars2 = var(xs2,slope2,slvars2,totSloRes2);
    auto chi2 = chisq(xs2,slope2,vars2);
    auto slchi2v2 = chisqv(xs2,slope2,slvars2,vars2);//?????????????????
	//chisqr probs
    auto chi2pf2 = 1 - TMath::Prob(chi2, nWts);

	//geo slope
    auto totSloResGeo = accum(plgrs);
    auto slopeg = plgsum/totSloResGeo;
    auto slopegerr = std::sqrt(1/totSloResGeo);

	//pairs slope
    auto totSloResPrs = accum(plrs);
    auto slopeprs = plsum/totSloResPrs;
    auto slopeprserr = std::sqrt(1/totSloResPrs);

	// eigens 0 = vector x, 1 = vector y, 2 = vector mag
    eigens.push_back(slope1);//3  aligned slope
    eigens.push_back(chi2pf1);//4 aligned slope chi sqr prob
    eigens.push_back(slope2);//5 unaligned slope ( uslope )
    eigens.push_back(chi2pf2);//6 uslope chi sqr prob
    eigens.push_back(rotangle);//7 aligned rotation angle
    eigens.push_back(nXSum);//8 # of entries ( rechits )
    eigens.push_back(orgrota);//9 unaligned rotation angle
    eigens.push_back(std::sqrt(varsl));//10 stdev aligned slope
    eigens.push_back(slope1err);//11 err aligned slope
    eigens.push_back(std::sqrt(vars2));//12 stdev unaligned slope
    eigens.push_back(slope2err);//13 errr unaligned slope
    eigens.push_back(slchi2v1);//14 chisqr like gof aligned slope
    eigens.push_back(slchi2v2);//15 chisqr like gof unaligned slope
    eigens.push_back(geoeigens[0]);//16 geoeigan x vec
    eigens.push_back(geoeigens[1]);//17 geoeigan y vec
    eigens.push_back(geoeigens[2]);//18 geoeigan mag vec
    eigens.push_back(ebside);//19 EB side
    eigens.push_back(taflip);//20 towards(+)/away(-) slope sign flip
    eigens.push_back(egflip);//21 eigan vector slope aligment flip
    eigens.push_back(slopeg);//22 geoeigan slope
    eigens.push_back(slopegerr);//23 geoeigan slope error
    eigens.push_back(slopeprs);//24 pairs slope method
    eigens.push_back(slopeprserr);//25 pairs slope method error
	//std::cout << "Slope egin : " << slope1 << " " << chi2pf1 << " " << rotangle << " " << std::sqrt(varsl) << " " << slope1err << std::endl;

    return eigens;
}//>>>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen_sph( vector<float> times, rhGroup rechits ){

//-----------------------------------------------------------------------------------------------------------------------------
//----------------------------------------  Make Hists function call  ---------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------------
void makehists::llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename, int pct ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("tree/llpgtree");
    //const std::string eosdir("root://cmseos.fnal.gov//store/user/jaking/");
    const std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");

    std::ifstream infile(infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
		std::cout << "--  for Fine #" << cnt << " moduls " << cnt%pct << " ";
		if( cnt%pct == 0 ){ 
        auto tfilename = eosdir + indir + str;
        std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
		}//<<>>if( cnt%4 == 0 ){
		else std::cout << " do not add file" << std::endl;
		cnt++;
    }//<<>>while (std::getline(infile,str))

	Init(fInTree);
	initHists();

	SetupDetIDsEB(DetIDMap);
	SetupDetIDsEE(DetIDMap);

    std::cout << "Setting up For Main Loop." << std::endl;
	int loopCounter(10000);
    auto nEntries = fInTree->GetEntries();
	//nEntries = 100;
	//loopCounter = 10;
    if(debug){ nEntries = 4; loopCounter = 1; }
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%loopCounter == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
        if(debug) std::cout << "*****************************************************************************" << std::endl;
        auto entry = fInTree->LoadTree(centry);
		if(debug) std::cout << " - getBranches " << std::endl;
		getBranches(entry);
		if(debug) std::cout << " - eventLoop " << std::endl;
		eventLoop(entry);
    }//<<>>for (Long64_t centry = 0; centry < nEntries; centry++)  end entry loop

    TFile* fOutFile = new TFile( outfilename.c_str(), "RECREATE" );
    fOutFile->cd();

    std::cout << "<<<<<<<< Write Output Maps and Hists <<<<<<<<<<<<<< " << std::endl;

	endJobs();
	for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]){ hist1d[it]->Write(); delete hist1d[it]; } }
    for( int it = 0; it < n2dHists; it++ ){ if(hist2d[it]){ hist2d[it]->Write(); delete hist2d[it]; } }
    for( int it = 0; it < n3dHists; it++ ){ if(hist3d[it]){ hist3d[it]->Write(); delete hist3d[it]; } }

	for( int it = 0; it < nEBEEMaps; it++ ){ 
		ebeeMapP[it]->Write(); delete ebeeMapP[it]; 								 
		ebeeMapT[it]->Write(); delete ebeeMapT[it]; 
		ebeeMapR[it]->Write(); delete ebeeMapR[it];
	}//<<>>for( int it = 0; it < nEBEEMaps; it++ )

    fOutFile->Close();
    std::cout << "llpgana_hist_maker : Thats all Folks!!" << std::endl;
}//<<>>void llpgana_hist_maker
//------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------

// ----------------------------------------------- event loop -------------------------------------------------------------------------
void makehists::eventLoop( Long64_t entry ){

    float   jetPTmin        = 100.0;// for energy/time comp
    int     jetIDmin        = 2; //3;
    float   jetETAmax       = 1.5; //1.5;
    int     minRHcnt        = 15; //32;
    float   minRHenr        = 2.0;
    float   bcMinEnergy     = 0.667;
    int     bcMinRHGrpSize  = 3;
    float   minEmf          = 0.0;//0.2

    //if( DEBUG ) std::cout << "Finding Events" << std::endl;
	//------------ event varibles ------------------

    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

    if( DEBUG ) std::cout << " -- Looping over " << nRecHits << " rechits" << std::endl;

    for( int it = 0; it < nRecHits; it++ ){

		auto id = (*rhID)[it];
		auto idinfo = DetIDMap[id];
		if( idinfo.ecal == ECAL::EB ){
			auto radius = hypo( (*rhPosX)[it], (*rhPosY)[it] );
			fillTH1( radius, hist1d[300]);
			fillTH1( (*rhTime)[it], hist1d[301]); 
			fillTH1( (*rhEnergy)[it], hist1d[302]); 
			hist2d[350]->Fill( (*rhTime)[it], (*rhEnergy)[it]);
			fillTH1( (*rhisOOT)[it], hist1d[305]);
		} else {
			fillTH1( (*rhTime)[it], hist1d[303]); 
			fillTH1( (*rhEnergy)[it], hist1d[304]); 
			hist2d[351]->Fill( (*rhTime)[it], (*rhEnergy)[it]);
			fillTH1( (*rhisOOT)[it], hist1d[306]);
		}//<<>>if( (*rhSubdet)[it] == 0 )

	}//<<>>for( int it = 0; it < nRecHits; it++ )

    if( true ) { // genparts lock
    if( DEBUG ) std::cout << "Finding genParticles" << std::endl;
    //------------  genparts ------------------------

	int nLlpPho(0);
    for( int it = 0; it < nGenParts; it++ ){

		if( (*genPdgId)[it] == 22 && (*genLLP)[it] > 0 ) nLlpPho++;
        auto genID = (*genPdgId)[it];
        if( abs((*genPdgId)[it]) > 2000000 ) genID = abs((*genPdgId)[it]) - 2000000 + 200;
        else if( abs((*genPdgId)[it]) > 1000000 ) genID = abs((*genPdgId)[it]) - 1000000 + 100;
        hist1d[200]->Fill( genID );	

    }//<<>>for( int it = 0; it < nGenParts; it++ )
	hist1d[201]->Fill( nLlpPho );

    }//<<>>if( false ) { // genparts lock


	if( false ) { // calojet lock
/*    if( DEBUG ) std::cout << "Finding calojets" << std::endl;

	//------------  calojets ------------------------

	for( int it = 0; it < nCaloJets; it++ ){

    	hist1d[250]->Fill( (*cljCMeanTime)[it] );//c mean
        hist1d[253]->Fill( (*cljCDrMeanTime)[it] );//c dr mean
    	hist1d[251]->Fill( (*cljSeedTOFTime)[it] );//lead time 
    	hist1d[252]->Fill( (*cljCMeanTime)[it] - (*cljSeedTOFTime)[it] );//diff

	}//<<>>for( int it = 0; it < nCaloJets; it++ )
*/	}//<<>>if( false ) { // calojet lock


    if( DEBUG ) std::cout << "Finding photons" << std::endl;
    //----------------- photons ------------------
	
    if( DEBUG ) std::cout << " - Looping over for " << nPhotons << " photons" << std::endl;
    for( int it = 0; it < nPhotons; it++ ){

		auto isCmb = not (*phoExcluded)[it];
        //auto isCmb = phoExcluded[it];
        auto phoGenIdx = (*genPhoIdx)[it];
        auto isLLP = (phoGenIdx >= 0)?((*genLLP)[phoGenIdx] > 0):false;
		auto isOOT = (*phoIsOotPho)[it];

		if( DEBUG ) std::cout << " -- pho id" << std::endl;
		auto rhIso = (*phoEcalRHSumEtConeDR04)[it] < ( 0.006*(*phoPt)[it] + 4.2 );
		auto hcalTowIso = (*phoHcalTowerSumEtBcConeDR04)[it] < ( 0.0025*(*phoPt)[it] + 2.2 );
        if( DEBUG ) std::cout << " -- pho id 1" << std::endl;
		auto hcTrkIsoL = (*phoTrkSumPtSolidConeDR04)[it] < ( 0.001*(*phoPt)[it] + 3.5 ); //hallow cone track iso
        auto hcTrkIsoT = (*phoTrkSumPtSolidConeDR04)[it] < ( 0.001*(*phoPt)[it] + 2 ); //hallow cone track iso
        if( DEBUG ) std::cout << " -- pho id 2" << std::endl;
        auto hadOverE = (*phohadTowOverEM)[it] < 0.05;
        auto sigmaIeieEE = (*phoSigmaIEtaIEta)[it] < 0.03; // tight only EE
        auto sigmaIeieEB = (*phoSigmaIEtaIEta)[it] < 0.013; // tight only EB

        if( DEBUG ) std::cout << " -- pho id set cuts" << std::endl;
		auto baseCut = rhIso && hcalTowIso && hadOverE;
		auto looseCut = baseCut && hcTrkIsoL;
		auto tightCut = baseCut && hcTrkIsoT;
		auto tightEB = tightCut && sigmaIeieEB;
		auto tightEE = tightCut && sigmaIeieEE;

		auto phoClass = tightCut?3:looseCut?2:1;
		hist1d[122]->Fill( phoClass );//122
		//hist2d[256]->Fill(looseCut,tightCut);

		if( DEBUG ) std::cout << " -- pho rechits" << std::endl;
		auto nrh = ((*phoRhIds)[it]).size();
        auto phoClstrR9 = clstrR9( (*phoRhIds)[it] );

        auto isEB = (*phoIsEB)[it];
		auto isTightEB = std::abs((*phoEta)[it]) < 1.45;

		auto isTight = phoClass == 3;
        auto isLoose = phoClass > 1;
        auto isFake = phoClass == 1;

		auto isClR9r46 = phoClstrR9 > 0.4 && phoClstrR9 < 0.6;
        auto isClR9r68 = phoClstrR9 > 0.6 && phoClstrR9 < 0.8;
        auto isClR9r26 = phoClstrR9 > 0.2 && phoClstrR9 < 0.6;

        auto goodRhCnt = nrh > 14;

        //auto usePho = true;
		//auto usePho = not isLLP && isFake;
        //auto usePho = not isLLP && isLoose;
        //auto usePho = not isLLP && isLoose && isEB && isTightEB;
        //auto usePho = isLLP && isFake;
        //auto usePho = isLLP && isLoose && isEB && isTightEB;
        //auto usePho = not isLLP && isLoose && phoExcluded[it];
        //auto usePho = isLoose && isEB && isTightEB;
        
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isClR9r26;
        //auto usePho = isLoose && isEB && isTightEB && isClR9r68;

        //auto usePho = isLoose && isEB && isTightEB && isCmb;
        //auto usePho = isLoose && isEB && isTightEB && not isCmb;
		//auto usePho = isLoose && isEB && isTightEB && isCmb && isLLP;
		auto usePho = isLoose && isEB && isTightEB && isCmb && isOOT;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isOOT && isLLP;
		//auto usePho = isLoose && isEB && isTightEB && not isCmb;


		if( DEBUG ) std::cout << " -- photons # " << it << " isLLP " << isLLP << " isOOT " << isOOT << " isCmb " << isCmb << std::endl;
		if( usePho && goodRhCnt ) { //----------------------------------- 

        hist1d[132]->Fill(1);
        if( isLLP ) hist1d[132]->Fill(2);
        if( isOOT ) hist1d[132]->Fill(3);
        if( isCmb ) hist1d[132]->Fill(4);

		//hist3d[0]->Fill( ((*phoRhIds)[it]).size(), phoClstrR9, phoClass ); 		
        //if( phoClass == 1 ){ hist2d[216]->Fill( nrh, phoClstrR9 ); hist2d[224]->Fill( nrh, (*phoR9)[it]);}// tightCut?3:looseCut?2:1 );
        //if( phoClass > 1 ) hist2d[217]->Fill( nrh, phoClstrR9 ); // tightCut?3:looseCut?2:1 );
        //if( phoClass > 2 ){ }//tightCut?3:looseCut?2:1
        hist2d[218]->Fill( nrh, phoClstrR9 ); 
		hist2d[225]->Fill( nrh, (*phoR9)[it]);


		//if( phoClass < 3 ) continue;

		hist1d[102]->Fill( phoClstrR9 );
		hist2d[223]->Fill( (*phoR9)[it], phoClstrR9 );

		if( DEBUG ) std::cout << " -- pho time" << std::endl;
        hist1d[119]->Fill( (*phoCMeanTime)[it] );//c mean
        hist1d[120]->Fill( (*phoSeedTOFTime)[it] );//lead time 
        hist1d[121]->Fill( (*phoCMeanTime)[it] - (*phoSeedTOFTime)[it] );//diff

		// count only rhs above a certian energy threshold?
		auto goodEiganRhCnt = true;
		//auto goodEiganRhCnt = nrh > 14;
		//auto goodEiganRhCnt = nrh && nrh < 50;

/*
		if( DEBUG ) std::cout << " - looping photons : filling Eigans" << std::endl;
        if( (*phoSc3dEx)[it] != -9 ){
            auto epanlge = getAngle( (*phoSc3dEx)[it], (*phoSc3dEy)[it] );
            hist1d[104]->Fill(epanlge);//etaphi angle
            auto ephypo3D = hypo( (*phoSc3dEx)[it], (*phoSc3dEy)[it] );
            auto etanlge = getAngle( ephypo3D, (*phoSc3dEz)[it] );
            hist1d[105]->Fill(etanlge);//etatim angle
            hist2d[211]->Fill( (*phoSc3dEx)[it], (*phoSc3dEy)[it] );
            hist2d[212]->Fill( (*phoSc3dEx)[it], (*phoSc3dEz)[it] );
            hist1d[115]->Fill( (*phoSc3dEv)[it] );
        }//<<>>if( phoSCEigen3D[0] != -999 )
        //if( (*phoSc2dEx)[it] != -9 ){
        //    auto sphanlge = getAngle( (*phoSc2dEx)[it], (*phoSc2dEy)[it] );
        //    hist1d[106]->Fill( sphanlge );//eliptical angle
        //    hist2d[210]->Fill( (*phoSc2dEx)[it], (*phoSc2dEy)[it] );
        //    hist1d[113]->Fill( (*phoSc2dEv)[it] );
        //    //hist2d[262]->Fill( sphanlge, (*phoSc2dEv)[it] );
        //}//<<>>if( phoSCEigen2D[0] != -999 )i
*/

		if( DEBUG ) std::cout << " - looping photons : calc 2D Eigans" << std::endl;

        hist1d[103]->Fill(nrh);

		auto isPho = (phoGenIdx >= 0)?(((*genPdgId)[phoGenIdx] == 22)?1:0):-1; 
        hist2d[204]->Fill(isPho,isLLP);

        if( DEBUG ) std::cout << " - looping photons : making ebee maps" << std::endl;
		auto doEPMap = nMaps < nEBEEMaps && isEB;
		if( doEPMap ){
			if( DEBUG ) std::cout << " -- looping over ebee maps:" << std::endl;
			bool fill(false);
			auto rhcol = (*phoRhIds)[it];
			if( nMaps < 6 ){ if( nrh < 25 && not isOOT ) fill = true; }
            else if( nMaps < 12 ){ if( nrh < 25 && isOOT ) fill = true; }
			else if( nMaps < 18 ){ if( nrh >= 25 && nrh < 50 && not isOOT ) fill = true; }
            else if( nMaps < 24 ){ if( nrh >= 25 && nrh < 50 && isOOT ) fill = true; }
            else if( nMaps < 30 ){ if( nrh >= 50 && not isOOT ) fill = true; }
			else { if( nrh >= 50 && isOOT ) fill = true; }
            if( DEBUG ) std::cout << " -- start fill of ebee maps:" << std::endl;
			if(fill){
				if( DEBUG ) std::cout << " -- Filling ebeeMapT : " << nMaps << " with nRH : " << nRecHits << std::endl;
                vector<float> times;
                for( int idx = 0; idx < nrh; idx++ ) times.push_back((*rhTime)[getRhIdx(rhcol[idx])]);
                auto mtime = mean(times);
            	for( int idx = 0; idx < nrh; idx++ ){
        			const auto rhIDX = getRhIdx(rhcol[idx]);
					auto idinfo = DetIDMap[rhcol[idx]];
            		const auto rhEtaPos = idinfo.i2;//recHitPos.eta();
            		const auto rhPhiPos = idinfo.i1;//recHitPos.phi();
					auto res = std::sqrt(sq2(3.64/(*rhEnergy)[rhIDX])+0.18);
					ebeeMapP[nMaps]->Fill( rhEtaPos, rhPhiPos,1);
                	ebeeMapP[nMaps]->Fill( nMaps, 1,1);
                    ebeeMapT[nMaps]->Fill( rhEtaPos, rhPhiPos,100+(*rhTime)[rhIDX] - mtime);
                    ebeeMapR[nMaps]->Fill( rhEtaPos, rhPhiPos,res);
            	}//<<>>for( idx = 0; idx < nRecHits; idx++ )
                nMaps++;
			}//<<>>if(fill)	
		}//<<>>if( nMaps < 36 )
        if( DEBUG ) std::cout << " - looping photons : Finished making ebee maps" << std::endl;

		if( isEB && goodEiganRhCnt ){
			if( DEBUG ) std::cout << " Finding Eigans for Photon W/ " << nrh << " rechits." << std::endl;
			auto tofTimes = getLeadTofRhTime( (*phoRhIds)[it], vtxX, vtxY, vtxZ );
			auto phoEigens2D = getRhGrpEigen_sph( tofTimes, (*phoRhIds)[it], hist2d[206], hist2d[207], hist2d[208], hist2d[209], hist2d[214],  hist2d[215] );

            auto eig1way = phoEigens2D[21];
            auto eigtaf = phoEigens2D[20];
			auto rotAngle = getAngle( eig1way*phoEigens2D[0], eig1way*phoEigens2D[1] );

            auto notAnglePeak1 = rotAngle > 0.2;
            auto notAnglePeak2 = rotAngle < 1.4 || rotAngle > 1.75;
            auto notAnglePeak3 = rotAngle < 2.9 || rotAngle > 3.35;
            auto notAnglePeak4 = rotAngle < 4.4 || rotAngle > 4.85;
            auto notAnglePeak5 = rotAngle < 6.1;
			auto notAnglePeak = notAnglePeak1 && notAnglePeak2 && notAnglePeak3 && notAnglePeak4 && notAnglePeak5;

			auto goodEigan = phoEigens2D[24] != -9;
            //auto useEigan = goodEigan;
			auto useEigan = goodEigan && notAnglePeak;

	        if( useEigan ){

                if( DEBUG ) std::cout << " -- start eigan hists fill p0:" << std::endl;
                auto sphanlge = getAngle( phoEigens2D[0], phoEigens2D[1] );
                hist1d[127]->Fill( sphanlge );//eliptical angle  Eigan Rotation Angle
                hist1d[106]->Fill( rotAngle );//eliptical angle  Eta Phi Angle 2D : algined
            	hist2d[210]->Fill( eig1way*phoEigens2D[0], eig1way*phoEigens2D[1] );
            	hist1d[113]->Fill( phoEigens2D[2] ); 
                if( DEBUG ) std::cout << " -- start eigan hists fill p0a:" << std::endl;
                //hist1d[129]->Fill( phoEigens2D[18] );
				hist1d[114]->Fill(eig1way*phoEigens2D[24]); // aligned slope
                hist1d[117]->Fill(eigtaf*phoEigens2D[24]); // TA slope
				hist2d[200]->Fill((*phoEta)[it],eig1way*phoEigens2D[24]); 
                hist2d[202]->Fill((*phoEta)[it],eigtaf*phoEigens2D[24]);
				hist2d[201]->Fill((*phoEta)[it],phoEigens2D[2]);
                if( DEBUG ) std::cout << " -- start eigan hists fill p0b:" << std::endl;
				hist1d[123]->Fill(phoEigens2D[10]);
				hist1d[124]->Fill(phoEigens2D[25]);
				hist1d[125]->Fill(phoEigens2D[12]);
				hist1d[126]->Fill(phoEigens2D[25]);
                if( DEBUG ) std::cout << " -- start eigan hists fill p0c:" << std::endl;
				hist2d[219]->Fill( phoEigens2D[10], phoEigens2D[25] );
                hist2d[220]->Fill( eigtaf*phoEigens2D[24], phoEigens2D[10] );
                hist2d[221]->Fill( eigtaf*phoEigens2D[24], phoEigens2D[25] );
                if( DEBUG ) std::cout << " -- start eigan hists fill p0d:" << std::endl;
                hist2d[222]->Fill( phoEigens2D[25], phoClstrR9 );
				//hist2d[226]->Fill( (*phoEta)[it], rotAngle ); 
                //hist2d[227]->Fill( (*phoPhi)[it], rotAngle );
                //hist2d[228]->Fill( rotAngle, phoEigens2D[2] );
                //hist2d[229]->Fill( phoEigens2D[9], eigtaf*phoEigens2D[24] );
                if( DEBUG ) std::cout << " -- eigan hists fill p1:" << std::endl;

                hist2d[232]->Fill( phoClstrR9, rotAngle );
                hist2d[233]->Fill( phoEigens2D[2], eigtaf*phoEigens2D[24] );
                hist2d[235]->Fill( rotAngle, phoEigens2D[2] );
                hist2d[234]->Fill( rotAngle, eigtaf*phoEigens2D[24] );
                hist2d[230]->Fill( phoClstrR9, phoEigens2D[2] );
                hist2d[231]->Fill( phoClstrR9, eigtaf*phoEigens2D[24] );

                hist1d[133]->Fill( phoEigens2D[18] );
                hist1d[134]->Fill( phoEigens2D[22] );
				auto geoRotAngle = getAngle( eig1way*phoEigens2D[16], eig1way*phoEigens2D[17] );
                hist1d[135]->Fill( geoRotAngle );
                if( DEBUG ) std::cout << " -- eigan hists fill p2:" << std::endl;

                hist2d[236]->Fill( phoEigens2D[18], rotAngle );
                hist2d[237]->Fill( eigtaf*phoEigens2D[22], rotAngle );
                hist2d[238]->Fill( eigtaf*phoEigens2D[22], phoEigens2D[2] );
                hist2d[239]->Fill( geoRotAngle, rotAngle );
                hist2d[240]->Fill( eigtaf*phoEigens2D[24], eigtaf*phoEigens2D[22] );
                hist2d[241]->Fill( geoRotAngle, eigtaf*phoEigens2D[24] );
                hist2d[242]->Fill( geoRotAngle, phoEigens2D[2] );
                hist2d[243]->Fill( phoEigens2D[18], phoEigens2D[2] );
                hist2d[244]->Fill( phoEigens2D[18], eigtaf*phoEigens2D[24] );

                if( DEBUG ) std::cout << " -- eigan hists fill p3:" << std::endl;
                hist2d[245]->Fill( phoClstrR9, geoRotAngle );
                hist2d[246]->Fill( phoEigens2D[18], eigtaf*phoEigens2D[22] );
                hist2d[249]->Fill( geoRotAngle, phoEigens2D[18] );
                hist2d[250]->Fill( geoRotAngle, eigtaf*phoEigens2D[22] );
                hist2d[247]->Fill( phoClstrR9, phoEigens2D[18] );
                hist2d[248]->Fill( phoClstrR9, eigtaf*phoEigens2D[22] );

			}//<<>>if( phoEigens2D[3] != -9 )
		}//<<>>if( phoIsEB && tightEB )
	
		if( DEBUG ) std::cout << " -- looping photons : filling ids" << std::endl;
    	hist1d[107]->Fill( (*phoIsPixelSeed)[it] );
        hist1d[108]->Fill( (*phoHadOverEM)[it] );
        hist1d[109]->Fill( (*phoR9)[it] );
		hist1d[128]->Fill( (*phoR9)[it] );
        //hist1d[110]->Fill( (*phoMipTotEnergy)[it] );
        hist1d[111]->Fill( (*phoEcalRHSumEtConeDR04)[it] );
        hist1d[112]->Fill( (*phoHcalTwrSumEtConeDR04)[it] );

        if( DEBUG ) std::cout << " -- looping photons : filling gen" << std::endl;
		
		if(phoGenIdx >= 0){

        	hist1d[100]->Fill( (*genPt)[phoGenIdx] );
			fillTH1( (*genPhoDr)[it], hist1d[130]); //->Fill( (*genPhoDr)[it] );
        	fillTH1( abs((*genPdgId)[phoGenIdx]), hist1d[131]);
			auto genParentExoID = 98;
			if( abs((*genLLP)[phoGenIdx]) > 2000000 ) genParentExoID = abs((*genLLP)[phoGenIdx]) - 2000000 + 200;
			else if( abs((*genLLP)[phoGenIdx]) > 1000000 ) genParentExoID = abs((*genLLP)[phoGenIdx]) - 1000000 + 100;
			if( abs((*genPdgId)[phoGenIdx]) == 0 ) genParentExoID = 95;
        	hist1d[101]->Fill( genParentExoID );

		}//<<>>if(phoGenIdx >= 0)

        }//<<>>if( usePho ) { //----------------------------------- 

    }//<<>>for( int it = 0; it < nPhotons; it++ )

	if( false ) { // electrons lock
/*  if( DEBUG ) std::cout << "Finding electrons" << std::endl;
	//-------- electrons --------------------------------------

	
    for( int it = 0; it < nElectrons; it++ ){

        hist1d[350]->Fill( (*eleCMeanTime)[it] );//c mean
        hist1d[351]->Fill( (*eleSeedTOFTime)[it] );//lead time 
        hist1d[352]->Fill( (*eleCMeanTime)[it] - (*eleSeedTOFTime)[it] );//diff

    }//<<>>for( int it = 0; it < nElectrons; it++ )
*/	}//<<>>if( false ) { // electrons lock


	//if( true ) { // JETS lock
	if( false ) { // JETS lock
/*    if( DEBUG ) std::cout << "Finding pfjets" << std::endl;
	//--------- jets --------------------------


	//// ---  base info  ----

	hist1d[9]->Fill(nJets);

    for( int it = 0; it < nJets; it++ ){

        const auto jetepafrac   = (*jetPHEF)[it] + (*jetELEF)[it];;
        const auto jetepe       = (*jetPHE)[it] + (*jetELE)[it];
        const auto jeteme       = (*jetCEMF)[it] + (*jetNEMF)[it];
        const auto jetemfrac    = jeteme/(*jetE)[it];
        const auto jetepfrac    = jetepe/(*jetE)[it];
		
        hist2d[61]->Fill(jetepafrac,jetepfrac);
        hist2d[62]->Fill(jetepfrac,jetemfrac);

        fillTH1((*jetPt)[it],hist1d[0]);//hist1d[12]->Fill(jet.pt());
        fillTH1((*jetPhi)[it],hist1d[1]);//hist1d[13]->Fill(jet.phi());
        fillTH1((*jetEta)[it],hist1d[2]);//hist1d[14]->Fill(jet.eta());

	}//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding genjet" << std::endl;
	//// --- genjet info -----------------------------------

    for( int it = 0; it < nJets; it++ ){

        hist2d[161]->Fill((*jetGenTime)[it],(*jetGenTOF)[it]);
        hist1d[50]->Fill((*jetGenTime)[it]);
    	hist1d[51]->Fill((*jetGenTOF)[it]);

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet dr times" << std::endl;
	//// ---  jet time dr method --------------------------- (*)[it]

	for( int it = 0; it < nJets; it++ ){

        fillTH1((*jetDRMuTime)[it],hist1d[22]);//hist1d[29]->Fill(jmutime);
        fillTH1((*jetDrRhCnt)[it],hist1d[19]);//hist1d[1]->Fill(rhCount);
        fillTH1((*jetDRTimeError)[it],hist1d[20]);//hist1d[2]->Fill(jterr);
        fillTH1((*jetDRTimeRMS)[it],hist1d[21]);//hist1d[3]->Fill(jtrms);
        fillTH1((*jetDRMedTime)[it],hist1d[23]);//hist1d[4]->Fill(jmedtime);
        fillTH1((*jetCDRMuTime)[it],hist1d[24]);//hist1d[6]->Fill(jcmutime);
        fillTH1((*jetCDRMedTime)[it],hist1d[25]);//hist1d[7]->Fill(jcmedtime);

		if( DEBUG ) std::cout << " - Finding jet dr times 0" << std::endl;
        hist2d[1]->Fill((*jetDRMuTime)[it],(*jetPt)[it]);
        hist2d[2]->Fill((*jetDRMuTime)[it],(*jetID)[it]);
        hist2d[3]->Fill((*jetDRMuTime)[it],(*jetNHF)[it]);//jetNHF
        hist2d[4]->Fill((*jetDRMuTime)[it],(*jetCHF)[it]);//jetCHF
        hist2d[5]->Fill((*jetDRMuTime)[it],(*jetNEMF)[it]);//jetNEMF
        hist2d[6]->Fill((*jetDRMuTime)[it],(*jetCEMF)[it]);//jetCEMF
        hist2d[7]->Fill((*jetDRMuTime)[it],(*jetMUF)[it]);//jetMUF
        hist2d[8]->Fill((*jetDRMuTime)[it],(*jetNHM)[it]);//jetNHM
        hist2d[9]->Fill((*jetDRMuTime)[it],(*jetCHM)[it]);//jetCHM
    
        if( DEBUG ) std::cout << " - Finding jet dr times 1" << std::endl;
        hist2d[10]->Fill((*jetDRMuTime)[it],(*jetDRMedTime)[it]);
        hist2d[24]->Fill((*jetDRMuTime)[it],(*jetDrRhCnt)[it]);
        hist2d[25]->Fill((*jetDRMedTime)[it],(*jetDrRhCnt)[it]);
        hist2d[11]->Fill((*jetDRMuTime)[it],(*jetDRTimeRMS)[it]);
        hist2d[12]->Fill((*jetDRMuTime)[it],(*jetDRTimeError)[it]);

        if( DEBUG ) std::cout << " - Finding jet dr times 2" << std::endl;
        //hist2d[32]->Fill((*jetDRMuTime)[it],(*jetDrLeadEta)[it]);
        if( DEBUG ) std::cout << " - Finding jet dr times 2a" << std::endl;
        //hist2d[33]->Fill((*jetDRMuTime)[it],(*jetDrLeadPhi)[it]);
        if( DEBUG ) std::cout << " - Finding jet dr times 2b" << std::endl;
        //hist2d[34]->Fill((*jetDRMuTime)[it],(*jetDrLeadEnr)[it]);
        if( DEBUG ) std::cout << " - Finding jet dr times 2c" << std::endl;
        //hist2d[35]->Fill((*jetDRMedTime)[it],(*jetDrLeadEta)[it]);
        //hist2d[36]->Fill((*jetDRMedTime)[it],(*jetDrLeadPhi)[it]);
        //hist2d[37]->Fill((*jetDRMedTime)[it],(*jetDrLeadEnr)[it]);
    
        if( DEBUG ) std::cout << " - Finding jet dr times 3" << std::endl;
        hist2d[13]->Fill((*jetDRMedTime)[it],(*jetPt)[it]);
        hist2d[14]->Fill((*jetDRMedTime)[it],(*jetID)[it]);
        hist2d[15]->Fill((*jetDRMedTime)[it],(*jetNHF)[it]);
        hist2d[16]->Fill((*jetDRMedTime)[it],(*jetCHF)[it]);
        hist2d[17]->Fill((*jetDRMedTime)[it],(*jetNEMF)[it]);
        hist2d[18]->Fill((*jetDRMedTime)[it],(*jetCEMF)[it]);
        hist2d[19]->Fill((*jetDRMedTime)[it],(*jetMUF)[it]);
        hist2d[20]->Fill((*jetDRMedTime)[it],(*jetNHM)[it]);
        hist2d[21]->Fill((*jetDRMedTime)[it],(*jetCHM)[it]);
		
    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet sc times" << std::endl;
	//// --- jet time SC method --------------------------

    if( DEBUG ) std::cout << " - Finding jet sc times " << std::endl;
    for( int it = 0; it < nJets; it++ ){
*/

/*
        if( (*jetSc3dEy)[it] != -9 ){
                auto epanlge = getAngle( (*jetSc3dEx)[it], (*jetSc3dEy)[it] );
                hist1d[58]->Fill(epanlge);//etaphi angle
                auto ephypo3D = hypo( (*jetSc3dEx)[it], (*jetSc3dEy)[it] );
                auto etanlge = getAngle( ephypo3D, (*jetSc3dEz)[it] );
                fillTH1(etanlge,hist1d[59]);//hist1d[64]->Fill(etanlge);//etatim angle
                hist2d[173]->Fill( (*jetSc3dEx)[it], (*jetSc3dEy)[it] );
                hist2d[174]->Fill( (*jetSc3dEx)[it], (*jetSc3dEz)[it] );
                hist1d[57]->Fill((*jetSc3dEv)[it]);
                if( (*jetSc3dEchisp)[it] > 0.95 && (*jetSc3dEv)[it] < 0.9 && (*jetSc3dEv)[it] > 0.7 ){
                    hist2d[102]->Fill( (*jetImpactAngle)[it], (*jetSc3dEslope)[it] );
                    hist2d[126]->Fill( (*jetEta)[it], (*jetSc3dEslope)[it] );
                    hist2d[127]->Fill( (*jetPhi)[it], (*jetSc3dEslope)[it] );
                }//<<>>if( (*jetSc3dEchisp)[it] > 0.95 ):
            }//<<>>if( (*jetSc3dEx)[it] != -999 )
            if( (*jetSc2dEx)[it] != -9 ){
                auto sphanlge = getAngle( (*jetSc2dEx)[it], (*jetSc2dEy)[it] );
                hist1d[60]->Fill(sphanlge);//eliptical angle
                hist2d[128]->Fill( (*jetSc2dEx)[it], (*jetSc2dEy)[it] );
                hist1d[61]->Fill((*jetSc2dEv)[it]);
                //hist2d[262]->Fill( sphanlge, (*jetSc2dEv)[it] );
                if( (*jetSc2dEchisp)[it] > 0.95 && (*jetSc2dEv)[it] < 0.9 && (*jetSc2dEv)[it] > 0.7 ){
                    if( (*jetEta)[it] > 1.0 ) hist1d[55]->Fill( (*jetSc2dEslope)[it] );
                    if( (*jetEta)[it] < 0.5 ) hist1d[56]->Fill( (*jetSc2dEslope)[it] );
                    hist2d[100]->Fill( (*jetEta)[it], (*jetSc2dEslope)[it] );
                    hist2d[121]->Fill( (*jetPhi)[it], (*jetSc2dEslope)[it] );
                    hist2d[101]->Fill( (*jetImpactAngle)[it], (*jetSc2dEslope)[it] );
                }//<<>>if( (*jetSc2dEchisp)[it] < 0.1 )
                //hist2d[263]->Fill( (*jetSc2dEslope)[it], (*jetSc2dEchisp)[it]);
        }//<<>>if( (*jetSc2dEx)[it] != -999 )
*/

/*
		if( DEBUG ) std::cout << " -- Filling jet sc times " << std::endl;

        hist2d[39]->Fill( (*jetDrRhCnt)[it], (*jetScRhCnt)[it] );
        hist1d[27]->Fill((*jetSCMedTime)[it]);//median
        hist1d[26]->Fill((*jetSCMuTime)[it]);//mean
        //hist1d[46]->Fill(jetSCTimeStats[4]);//rms
        //hist1d[50]->Fill(jetSCTimeStats[5]);//skew
        hist1d[28]->Fill((*jetCSCMuTime)[it]);//c mean
        hist1d[29]->Fill((*jetCSCMedTime)[it]);//c med   

        hist1d[42]->Fill((*jetCSCMuTime)[it]-(*jetGenTime)[it]);

        hist2d[53]->Fill( (*jetScEMF)[it], (*jetEMEnrFrac)[it] );
        hist2d[56]->Fill( (*sJetScRhEnergy)[it], (*jetEMEnergy)[it] );
        hist2d[59]->Fill( (*jetDrEMF)[it], (*jetScEMF)[it] );

        //hist2d[63]->Fill( (*jetEta)[it], jetSCTimeStats[7] );
        //hist2d[64]->Fill( (*jetEta)[it]etaMmt, jetSCTimeStats[7] );
        //hist2d[65]->Fill( (*jetPhi)[it]phiMnt, jetSCTimeStats[7] );
        //hist2d[66]->Fill( (*jetEta)[it]phiMnt, jetSCTimeStats[7] );
        //hist2d[67]->Fill( jetMaxD, jetSCTimeStats[7] );
        //hist2d[68]->Fill( jetConPtDis, jetSCTimeStats[7] );
        //hist2d[204]->Fill( jetConEtaPhiSprd, jetSCTimeStats[7] );
        //hist2d[205]->Fill( jetArea, jetSCTimeStats[7] );
        //hist2d[255]->Fill( jetNCarry, jetSCTimeStats[7] );
        //hist2d[256]->Fill( jetNConst, jetSCTimeStats[7] );

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet bc times" << std::endl;
	//// --- jet time BC method ---------------------------

    for( int it = 0; it < nJets; it++ ){

        hist2d[38]->Fill((*jetDrRhCnt)[it],(*jetBcRhCnt)[it]);
        hist1d[41]->Fill((*jetBcGrpCnt)[it]);
        hist1d[30]->Fill((*jetCBCMedTime)[it]);//c med
        hist1d[31]->Fill((*jetCBCMuTime)[it]);//c mu
        hist2d[54]->Fill( (*jetBcEMFr)[it], (*jetEMEnrFrac)[it] );
        hist2d[57]->Fill( (*jetBcSumRHEnr)[it], (*jetEMEnergy)[it] );
        hist2d[58]->Fill( (*jetBcSumRHEnr)[it], (*sJetScRhEnergy)[it] );
        hist2d[60]->Fill( (*jetScEMF)[it], (*jetBcEMFr)[it] );

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet genjet comps" << std::endl;
    // ****************************  photon/electron to kid pfcand -> SC matcher ********************************

    for( int it = 0; it < nJets; it++ ){

        auto jetERatio = (*jetE)[it]/(*jetGenEnergy)[it];
        auto difSCTime = std::abs((*jetGenTime)[it]-(*jetCSCMuTime)[it]);
        auto difDrTime = std::abs((*jetGenTime)[it]-(*jetCDRMuTime)[it]);
        //if( (*jetCSCMuTime)[it] < -27.9 ) (*jetCSCMuTime)[it] = -25.0;
        //if( (*jetCDRMuTime)[it] < -27.9 ) (*jetCDRMuTime)[it] = -25.0;
        if( (*jetCSCMuTime)[it] < -25.0 || (*jetGenTime)[it] < -25.0 ) difSCTime = 100.0;
        if( (*jetCDRMuTime)[it] < -25.0 || (*jetGenTime)[it] < -25.0 ) difDrTime = 100.0;
        auto hasGoodGenTime = (*jetGenTime)[it] > -25.0;
        auto etaCut = std::abs((*jetGenEta)[it]) < 1.5 && std::abs((*jetEta)[it]) < 1.5;
        auto genEnergyCut = (*jetGenEnergy)[it] > 0.0;
        auto genVarCut = (*jetGenLLPPurity)[it] > 0.0 && (*jetGenTimeVar)[it] < 100;
        //auto hasGoodGenSCMatch = difSCTime < 0.8;
        auto hasGoodGenSCMatch = difSCTime < 20.0;
        auto genNoCut = true;// no cut
        //auto genCutTime = 9.0 - 4.0 * jetERatio;
        //auto genPhSpaceCut = (*jetGenTime)[it] > genCutTime;
        //auto genCutDr = (*jetGenDrMatch)[it] > 0.04;
        //auto genCutSCdiff = difSCTime > 4;
        //auto genDrPhSpaceCut = genCutSCdiff && genCutDr;
        auto hasGoodSCTime = (*jetCSCMuTime)[it] > -14.0;

        if( genNoCut ){

            hist1d[52]->Fill((*jetGenTimeLLP)[it]);
            hist1d[53]->Fill((*jetGenLLPPurity)[it]);
            hist2d[166]->Fill((*jetGenLLPPurity)[it],(*jetGenTimeVar)[it]);
            hist1d[54]->Fill((*jetGenNKids)[it]);
            hist2d[167]->Fill((*jetGenTimeVar)[it],(*jetGenNKids)[it]);
            hist2d[168]->Fill((*jetGenLLPPurity)[it],(*jetGenNKids)[it]);
            hist1d[47]->Fill((*jetGenTime)[it]);
            hist1d[43]->Fill((*jetGenImpactAngle)[it]);
            hist1d[44]->Fill((*jetGenDrMatch)[it]);
            hist1d[48]->Fill((*jetGenTimeVar)[it]);
            hist1d[49]->Fill((*jetGenNextBX)[it]);
            hist1d[45]->Fill( difSCTime );
            hist1d[46]->Fill( difDrTime );
            hist2d[153]->Fill( jetERatio, (*jetGenTime)[it] );
            hist2d[169]->Fill( difSCTime, (*jetGenDrMatch)[it] );
            hist2d[170]->Fill( (*jetGenTime)[it], (*jetEMEnrFrac)[it] );
            hist2d[171]->Fill( (*jetEMEnrFrac)[it], (*jetGenDrMatch)[it] );
            hist2d[172]->Fill( jetERatio, (*jetGenDrMatch)[it] );

        //}//<<>>if( genSpaceCut )

        //if( hasGoodGenSCMatch && etaCut && genEnergyCut && genVarCut && hasGoodGenTime && hasGoodSCTime ){

            hist2d[150]->Fill( (*jetE)[it], (*jetGenEnergy)[it] );
            hist2d[151]->Fill( jetERatio, (*jetGenTime)[it] );
            hist2d[162]->Fill( jetERatio, (*jetCSCMuTime)[it] );
            hist2d[163]->Fill( jetERatio, (*jetCDRMuTime)[it] );
            hist2d[152]->Fill( (*jetEMEnrFrac)[it], (*jetCSCMuTime)[it] );
            hist2d[155]->Fill( jetERatio, difSCTime );
            hist2d[156]->Fill( (*jetCDRMuTime)[it], (*jetCSCMuTime)[it] );
            hist2d[154]->Fill( (*jetGenTime)[it], (*jetCDRMuTime)[it] );
            hist2d[157]->Fill( (*jetGenTime)[it], (*jetCSCMuTime)[it] );
            hist2d[158]->Fill( (*jetGenTime)[it], (*jetGenEnergy)[it] );
            hist2d[159]->Fill( (*jetGenDrMatch)[it], difSCTime );
            hist2d[160]->Fill( (*jetGenDrMatch)[it], difDrTime );
            hist2d[164]->Fill( (*jetGenTimeVar)[it], difSCTime );
            hist2d[165]->Fill( (*jetGenLLPPurity)[it], difSCTime );

        }//<<>>if( (*jetCSCMuTime)[it] > -28.0 && (*jetGenTime)[it] > -28.0 )

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jets info" << std::endl;
	//// --- jet time resolution -------------------------

    hist1d[5]->Fill(jetHt);
    hist1d[10]->Fill(nGoodDrJets);
    hist1d[11]->Fill(nGoodScJets);
    hist1d[12]->Fill(nGoodBcJets);
    hist1d[13]->Fill(nUnJets);
    if( nUnJets != 0 ) hist1d[14]->Fill(float(nJets)/nUnJets);
    if( nGoodScJets != 0 ) hist1d[18]->Fill(float(nGoodBcJets)/nGoodScJets);
    if( nJets != 0 ){
        hist1d[15]->Fill(float(nGoodDrJets)/nJets);
        hist1d[16]->Fill(float(nGoodScJets)/nJets);
        hist1d[17]->Fill(float(nGoodBcJets)/nJets);
    }//<<>>if( nUnJets != 0 )

    if( DEBUG ) std::cout << "Finding jet dt pairs" << std::endl;
    //-----------------------------------------------------------------------------------------------------
    // ***************************** d jetTime for back-to-back high pt jets  *****************************
    auto dijetIdCut = 1;
    auto dijetPtMin = 200.0;
    auto difPtLmt = 0.8;
    auto htPctLmt = 0.8;
    auto dPhiLmt = 2.8;


    for( int q = 0; q < nJets; q++ ) hist1d[40]->Fill( getdt((*jetCSCMuTime)[q],(*jetCBCMuTime)[q]) );

    if( DEBUG ) std::cout << "Finding jet dt pairs" << std::endl;
    for ( int q = 0; q < nJets; q++ ){
        for ( int p = q+1; p < nJets; p++ ){

            //if( DEBUG ) std::cout << " - filter jet pairs" << std::endl;
            if( (*jetPt)[q] < dijetPtMin ) continue;
            auto diffPt = (*jetPt)[p]/(*jetPt)[q];
            hist1d[6]->Fill(diffPt);
            if( diffPt < difPtLmt ) continue;
            auto htPct= ((*jetPt)[q]+(*jetPt)[p])/jetHt;
            hist1d[7]->Fill(htPct);
            if( htPct < htPctLmt ) continue;
            //auto dPhi = reco::deltaPhi((*jetPhi)[q],(*jetPhi)[p]);
            auto dPhi = std::abs(dltPhi((*jetPhi)[q],(*jetPhi)[p])); 
            hist1d[8]->Fill(dPhi);
            //if( dPhi > dPhiLmt ) continue;

            if( DEBUG ) std::cout << " - get jet pair dt" << std::endl;
            auto dTmu = getdt( (*jetDRMuTime)[q], (*jetDRMuTime)[p] );
            auto dTmed = getdt( (*jetDRMedTime)[q], (*jetDRMedTime)[p] );
            auto dTcmu = getdt( (*jetCDRMuTime)[q], (*jetCDRMuTime)[p] );
            auto dTcmed = getdt( (*jetCDRMedTime)[q], (*jetCDRMedTime)[p] );
            if( DEBUG ) std::cout << "dT dR      : " << dTmu <<  " " << (*jetDRMuTime)[q] << " " << (*jetDRMuTime)[p] << std::endl;
            auto dTmusc = getdt( (*jetCSCMuTime)[q], (*jetCSCMuTime)[p] );
            auto dTmedsc = getdt( (*jetSCMedTime)[q], (*jetSCMedTime)[p] );
            auto dTcmusc = getdt( (*jetCSCMuTime)[q], (*jetCSCMuTime)[p] );
            auto dTcmedsc = getdt( (*jetCSCMedTime)[q], (*jetCSCMedTime)[p] );
            if( DEBUG ) std::cout << "dT SC      : " << dTmusc <<  " " << (*jetCSCMuTime)[q] << " " << (*jetCSCMuTime)[p] << std::endl;
            auto dTcmubc = getdt( (*jetCBCMuTime)[q], (*jetCBCMuTime)[p] );
            if( DEBUG ) std::cout << "dT cMu BC  : " << dTcmubc <<  " " << (*jetCBCMuTime)[q] << " " << (*jetCBCMuTime)[p] << std::endl;
            auto dTcmedbc = getdt( (*jetCBCMedTime)[q], (*jetCBCMedTime)[p] );
            if( DEBUG ) std::cout << "dT cMed BC : " << dTcmedbc <<  " " << (*jetCBCMedTime)[q] << " " << (*jetCBCMedTime)[p] << std::endl;
            if( DEBUG ) std::cout << " - fill hists" << std::endl;

            auto dtThrs = -2.5;//removes default dt values from getdt
            if( dTmu > dtThrs ) hist1d[3]->Fill(dTmu);
            if( dTmed > dtThrs ) hist1d[4]->Fill(dTmed);
            if( dTcmu > dtThrs ) hist1d[34]->Fill(dTcmu);
            if( dTcmed > dtThrs ) hist1d[35]->Fill(dTcmed);

            if( dTmedsc > dtThrs ) hist1d[38]->Fill(dTmedsc);
            if( dTmusc > dtThrs ) hist1d[39]->Fill(dTmusc);
            if( dTcmusc > dtThrs ) hist1d[36]->Fill(dTcmusc);
            if( dTcmedsc > dtThrs ) hist1d[37]->Fill(dTcmedsc);

            if( dTcmubc > dtThrs ) hist1d[32]->Fill(dTcmubc);
            if( dTcmedbc > dtThrs ) hist1d[33]->Fill(dTcmedbc);

            hist2d[22]->Fill(dTmu,nJets);
            hist2d[26]->Fill(dTmu,diffPt);
            hist2d[27]->Fill(dTmu,htPct);
            hist2d[28]->Fill(dTmu,dPhi);
            hist2d[23]->Fill(dTmed,nJets);
            hist2d[29]->Fill(dTmu,diffPt);
            hist2d[30]->Fill(dTmu,htPct);
            hist2d[31]->Fill(dTmu,dPhi);

            if( DEBUG ) std::cout << " - fill dt vs eff e hists" << std::endl;
            auto effje = effMean((*jetPHE)[p],(*jetPHE)[q]);
            hist2d[30]->Fill(dTmusc,effje);
            hist2d[31]->Fill(dTmu,effje);

        }//<<>>for ( int p = q+1; p < nJets; p++ )
    }//<<>>for ( int q = 0; q < nJets; q++ )

    if( DEBUG ) std::cout << "Finding calojet dt pairs" << std::endl;
    for ( int q = 0; q < nCaloJets; q++ ){
        for ( int p = q+1; p < nCaloJets; p++ ){

            //if( DEBUG ) std::cout << " - filter calojet pairs" << std::endl;
            if( (*cljPt)[q] < dijetPtMin ) continue;
            auto diffPt = (*cljPt)[p]/(*cljPt)[q];
            hist1d[254]->Fill(diffPt);
            if( diffPt < difPtLmt ) continue;
            //auto htPct= ((*cljPt)[q]+(*cljPt)[p])/jetHt;
            //hist1d[7]->Fill(htPct);
            //if( htPct < htPctLmt ) continue;
            //auto dPhi = reco::deltaPhi((*jetPhi)[q],(*jetPhi)[p]);
            auto dPhi = dltPhi((*cljPhi)[q],(*cljPhi)[p]);
            hist1d[255]->Fill(dPhi);
            if( dPhi < dPhiLmt ) continue;

            if( DEBUG ) std::cout << " - get calojet pair dt" << std::endl;
            auto dTmu = getdt( (*cljCMeanTime)[q], (*cljCMeanTime)[p] );

            auto dtThrs = -2.5;//removes default dt values from getdt
            if( dTmu > dtThrs ) hist1d[256]->Fill(dTmu);

        }//<<>>for ( int p = q+1; p < nJets; p++ )
    }//<<>>for ( int q = 0; q < nJets; q++ )
*/	}//<<>>if( false ) { // JETS lock


}//<<>>void makehists::eventLoop( Long64_t entry )

void makehists::getBranches( Long64_t entry ){

	//fChain->GetEntry(entry);

   //b_run->GetEntry(entry);   //!
   //b_lumi->GetEntry(entry);   //!
   //b_event->GetEntry(entry);   //!

   b_nVtx->GetEntry(entry);   //!
   b_vtxX->GetEntry(entry);   //!
   b_vtxY->GetEntry(entry);   //!
   b_vtxZ->GetEntry(entry);   //!

   //b_metSumEt->GetEntry(entry);   //!
   //b_metPt->GetEntry(entry);   //!
   //b_metPx->GetEntry(entry);   //!
   //b_metPy->GetEntry(entry);   //!
   //b_metPhi->GetEntry(entry);   //!
   //b_metCSumEt->GetEntry(entry);   //!
   //b_metCPx->GetEntry(entry);   //!
   //b_metCPy->GetEntry(entry);   //!

   //b_jetHt->GetEntry(entry);   //!
   //b_nJets->GetEntry(entry);   //!
   //b_nGoodDrJets->GetEntry(entry);   //!
   //b_nGoodScJets->GetEntry(entry);   //!
   //b_nGoodBcJets->GetEntry(entry);   //!
   //b_nUnJets->GetEntry(entry);   //!
   //b_jetE->GetEntry(entry);   //!

   b_nPhotons->GetEntry(entry);   //!
   b_phoIsOotPho->GetEntry(entry);   //!
   b_phoExcluded->GetEntry(entry);   //!
   b_phoSeedTOFTime->GetEntry(entry);   //!
   b_phoCMeanTime->GetEntry(entry);   //!
   b_phoSc2dEv->GetEntry(entry);   //!
   b_phoPt->GetEntry(entry);   //!
   b_phoEnergy->GetEntry(entry);   //!
   b_phoPhi->GetEntry(entry);   //!
   b_phoEta->GetEntry(entry);   //!
   b_phoPx->GetEntry(entry);   //!
   b_phoPy->GetEntry(entry);   //!
   b_phoPz->GetEntry(entry);   //!
   b_phoRhIds->GetEntry(entry);   //!
   b_phoIsPFPhoton->GetEntry(entry);   //!
   b_phoIsStdPhoton->GetEntry(entry);   //!
   b_phoHasConTracks->GetEntry(entry);   //!
   b_phoIsPixelSeed->GetEntry(entry);   //!
   b_phoIsPhoton->GetEntry(entry);   //!
   b_phoIsEB->GetEntry(entry);   //!
   b_phoIsEE->GetEntry(entry);   //!
   b_phoHadOverEM->GetEntry(entry);   //!
   b_phoHadD1OverEM->GetEntry(entry);   //!
   b_phoHadD2OverEM->GetEntry(entry);   //!
   b_phoHadOverEMVaid->GetEntry(entry);   //!
   b_phohadTowOverEM->GetEntry(entry);   //!
   b_phohadTowD10OverEM->GetEntry(entry);   //!
   b_phohadTowD20OverEM->GetEntry(entry);   //!
   b_phohadTowOverEMValid->GetEntry(entry);   //!
   b_phoE1x5->GetEntry(entry);   //!
   b_phoE2x5->GetEntry(entry);   //!
   b_phoE3x3->GetEntry(entry);   //!
   b_phoE5x5->GetEntry(entry);   //!
   b_phoMaxEnergyXtal->GetEntry(entry);   //!
   b_phoSigmaEtaEta->GetEntry(entry);   //!
   b_phoSigmaIEtaIEta->GetEntry(entry);   //!
   b_phoR1x5->GetEntry(entry);   //!
   b_phoR2x5->GetEntry(entry);   //!
   b_phoR9->GetEntry(entry);   //!
   b_phoFull5x5_e1x5->GetEntry(entry);   //!
   b_phoFull5x5_e2x5->GetEntry(entry);   //!
   b_phoFull5x5_e3x3->GetEntry(entry);   //!
   b_phoFull5x5_e5x5->GetEntry(entry);   //!
   b_phoFull5x5_maxEnergyXtal->GetEntry(entry);   //!
   b_phoFull5x5_sigmaEtaEta->GetEntry(entry);   //!
   b_phoFull5x5_sigmaIEtaIEta->GetEntry(entry);   //!
   b_phoFull5x5_r9->GetEntry(entry);   //!
   b_phoEcalRHSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalTwrSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalDepth1TowerSumEtConeDR04->GetEntry(entry);   //!
   b_phoCalDepth2TowerSumEtConeDR04->GetEntry(entry);   //!
   b_phoHcalTowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoHcalDepth1TowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoHcalDepth2TowerSumEtBcConeDR04->GetEntry(entry);   //!
   b_phoTrkSumPtSolidConeDR04->GetEntry(entry);   //!
   b_phoTrkSumPtHollowConeDR04->GetEntry(entry);   //!
   b_phoNTrkSolidConeDR04->GetEntry(entry);   //!
   b_phoNTrkHollowConeDR04->GetEntry(entry);   //!
   b_genPhoIdx->GetEntry(entry);   //!
   b_genPhoDr->GetEntry(entry);   //!

   b_nGenParts->GetEntry(entry);   //!
   b_genPt->GetEntry(entry);   //!
   //b_genEnergy->GetEntry(entry);   //!
   //b_genPhi->GetEntry(entry);   //!
   //b_genEta->GetEntry(entry);   //!
   //b_genPx->GetEntry(entry);   //!
   //b_genPy->GetEntry(entry);   //!
   //b_genPz->GetEntry(entry);   //!
   b_genPdgId->GetEntry(entry);   //!
   b_genLLP->GetEntry(entry);   //!

   b_nRecHits->GetEntry(entry);   //!
   b_rhPosX->GetEntry(entry);   //!
   b_rhPosY->GetEntry(entry);   //!
   b_rhPosZ->GetEntry(entry);   //!
   b_rhPosEta->GetEntry(entry);   //!
   b_rhPosPhi->GetEntry(entry);   //!
   b_rhEnergy->GetEntry(entry);   //!
   b_rhTime->GetEntry(entry);   //!
   b_rhTOF->GetEntry(entry);   //!
   b_rhID->GetEntry(entry);   //!
   b_rhisOOT->GetEntry(entry);   //!

}//<<>>void makehists::getBranches( Long64_t entry )

void makehists::endJobs(){

    normTH1D(hist1d[3]);
    normTH1D(hist1d[4]);
    normTH1D(hist1d[34]);
    normTH1D(hist1d[35]);
    normTH1D(hist1d[256]);
    normTH1D(hist1d[38]);
    normTH1D(hist1d[39]);
    normTH1D(hist1d[36]);
    normTH1D(hist1d[37]);
    normTH1D(hist1d[32]);
    normTH1D(hist1d[33]);

    thresDivTH2D( hist2d[206], hist2d[207], 0 );
    thresDivTH2D( hist2d[208], hist2d[209], 0 );

    profileTH2D( hist2d[215], hist1d[116], hist1d[118] );

    normTH2D(hist2d[100]);
    normTH2D(hist2d[121]);
    normTH2D(hist2d[126]);
    normTH2D(hist2d[127]);
    normTH2D(hist2d[202]);
    normTH2D(hist2d[201]);
    normTH2D(hist2d[200]);
    //normTH2D(hist2d[226]);
    //normTH2D(hist2d[227]);

}//<<>>void makehists::endJobs()

void makehists::initHists(){

	for( int it = 0; it < n1dHists; it++ ){ hist1d[it] = NULL; }
    for( int it = 0; it < n2dHists; it++ ){ hist2d[it] = NULL; }
    for( int it = 0; it < n3dHists; it++ ){ hist3d[it] = NULL; }
    //for( int it = 0; it < n1dHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

	// jet time
    //int jtdiv(400);
    //float jtran(8);
    int jtdiv(1000);
    float jtran(25);
    int jdtdiv(160);
    float jdtran(4);
    int rhcnt(160);

	// jet id stuff
    //auto stddiv = 120;
    //auto stdtran = 3;

	// eigan cluster 2d maps
    auto cldiv = 1200;
    auto cltrn = 30;

	// eigan cluster 3d maps
    auto cl3ddiv = 160;
    auto cl3dtrn = 4;
    auto cl3ddiv1 = 160;
    auto cl3dtrn1 = 4;

	// photon time phi - eta map
    auto clsphdiv = 80;
    auto clsphtrn = 4;

	// time profile maps
    auto cwdiv = 1000;
    auto cwtrn = 5;

	// time cl slope plots
    auto slmax = 2;
    auto slmin = -2;
    auto sldiv = 400;

    auto sl3max = 2;
    auto sl3min = -2;
    auto sl3div = 400;

    //auto chimax = 1.01;
    //auto chimin = 0.91;
    //auto chidiv = 2880;


	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

	//------- jets 0 - 99

    hist1d[0] = new TH1D("jetPt", "jetPt", 500, 0, 5000);
    hist1d[1] = new TH1D("jetPhi", "jetPhi", 700, -3.5, 3.5);
    hist1d[2] = new TH1D("jetEta", "jetEta", 700, -3.5, 3.5);
    hist1d[3] = new TH1D("jetdtmu", "jetdtmu", jdtdiv, -1*jdtran, jdtran);
    hist1d[4] = new TH1D("jetdtmed", "jetdtmed", jdtdiv, -1*jdtran, jdtran);
    hist1d[5] = new TH1D("jetHt", "jetHt", 1000, 0, 5000);
    hist1d[6] = new TH1D("diffPt", "diffPt", 1000, 0, 10);
    hist1d[7] = new TH1D("htPct", "htPct", 100, 0, 1);
    hist1d[8] = new TH1D("dPhi", "dPhi", 70, -3.5, 3.5);

    hist1d[9] = new TH1D("nJet", "nJets", 21, -0.5, 20.5);
    hist1d[10] = new TH1D("nGoodDrJets", "nGoodDrJets", 21, -0.5, 20.5);
    hist1d[11] = new TH1D("nGoodScJets", "nGoodScJets", 21, -0.5, 20.5);
    hist1d[12] = new TH1D("nGoodBcJets", "nGoodBcJets", 21, -0.5, 20.5);
    hist1d[13] = new TH1D("nUnJets", "nUnJets", 101, -0.5, 100.5);
    hist1d[14] = new TH1D("pJets", "pJets", 110, 0, 1.1);
    hist1d[15] = new TH1D("pGoodDrJets", "pGoodDrJets", 110, 0, 1.1);
    hist1d[16] = new TH1D("pGoodScJets", "pGoodScJets", 110, 0, 1.1);
    hist1d[17] = new TH1D("pGoodBcJets", "pGoodBcJets", 110, 0, 1.1);
    hist1d[18] = new TH1D("pGoodBcToScJets", "pGoodBcToScJets", 110, 0, 1.1);

    hist1d[19] = new TH1D("jetDRRHMulti", "jetDRRHMulti", rhcnt, 0, rhcnt);
    hist1d[20] = new TH1D("jetDRTimeError", "jetDRTimeError", 300, 0, 3);
    hist1d[21] = new TH1D("jetDRTimeRMS", "jetDRTimeRMS", 200, 0, 20);
    hist1d[22] = new TH1D("jetDrMuTime", "jetDrMuTime", jtdiv, -1*jtran, jtran);
    hist1d[23] = new TH1D("jetDRMedTime", "jetMedDRTime", jtdiv, -1*jtran, jtran);
    hist1d[24] = new TH1D("jetCDRMuTime", "jetCDRMuTime", jtdiv, -1*jtran, jtran);
    hist1d[25] = new TH1D("jetCDRMedTime", "jetCDRMedTime", jtdiv, -1*jtran, jtran);

    hist1d[26] = new TH1D("jetSCmuTime", "jetSCmuTime", jtdiv, -1*jtran, jtran);
    hist1d[27] = new TH1D("jetSCmedTime", "jetSCmedTime", jtdiv, -1*jtran, jtran);
    hist1d[28] = new TH1D("jetCSCMuTime", "jetCSCMuTime", jtdiv, -1*jtran, jtran);
    hist1d[29] = new TH1D("jetCSCMedTime", "jetCSCMedTime", jtdiv, -1*jtran, jtran);

    hist1d[30] = new TH1D("jetCBCMedTime", "jetCBCMedTime", jtdiv, -1*jtran, jtran);
    hist1d[31] = new TH1D("jetCBCMuTime", "jetCBCMuTime", jtdiv, -1*jtran, jtran);


    hist1d[32] = new TH1D("jetcmudtbc", "jetcmudtbc", jdtdiv, -1*jdtran, jdtran);
    hist1d[33] = new TH1D("jetcmeddtbc", "jetcmeddtbc", jdtdiv, -1*jdtran, jdtran);
    hist1d[34] = new TH1D("jetcmudtdr", "jetcmudtdr", jdtdiv, -1*jdtran, jdtran);
    hist1d[35] = new TH1D("jetcmeddtdr", "jetcmeddtdr", jdtdiv, -1*jdtran, jdtran);
    hist1d[36] = new TH1D("jetcmudtsc", "jetcmudtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[37] = new TH1D("jetcmeddtsc", "jetcmeddtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[38] = new TH1D("jetmeddtsc", "jetmeddtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[39] = new TH1D("jetmudtsc", "jetmudtsc", jdtdiv, -1*jdtran, jdtran);

    hist1d[40] = new TH1D("scbcdt", "scbcdt", jdtdiv, -1*jdtran, jdtran);
    hist1d[41] = new TH1D("nBCinJet", "nBCinJet", 11, -0.5, 10.5);

    hist1d[42] = new TH1D("jetmudtgen", "jetmudtgen", jdtdiv, -1*jdtran, jdtran);
    hist1d[43] = new TH1D("genJetImpactAngle", "genJetImpactAngle", 6600, -0.2, 6.4);

    hist1d[44] = new TH1D("genJetDrMatchJet", "genJetDrMatchJet", 100, 0, 0.1);
    hist1d[45] = new TH1D("genJetSCTimeDiff", "genJetSCTimeDiff", 300, 0, 30);
    hist1d[46] = new TH1D("genJetDrTimeDiff", "genJetSCTimeDiff", 300, 0, 30);

    hist1d[47] = new TH1D("jetGenTimeCut", "jetGenTimeCut", jtdiv, -1*jtran, jtran);
    hist1d[48] = new TH1D("jetGenTimeVar", "jetGenTimeVar", 408, -2, 100);
    hist1d[49] = new TH1D("jetGenTimeNextBX", "jetGenTimeNextBX", 3, -1, 2);
    hist1d[50] = new TH1D("jetGenTime", "jetGenTime", jtdiv, -1*jtran, jtran);
    hist1d[51] = new TH1D("jetGenTOF", "jetGenTOF", 300, 0, 30);
    hist1d[52] = new TH1D("jetGenTimeIsLLP", "jetGenTimeIsLLP", 3, -1, 2);
    hist1d[53] = new TH1D("jetGenTimeLLPPurity", "jetGenTimeLLPPurity", 100, 0, 1);
    hist1d[54] = new TH1D("jetGenNKids", "jetGenNKids", 100, 0, 100);

    hist1d[55] = new TH1D("jetClEtaTimeSlopeRangeA", "Jet Cluster Eta Time Slope Eta > 1.0", 480, -240, 240);
    hist1d[56] = new TH1D("jerClEtaTimeSlopeRangeB", "Jet Cluster Eta Time Slope Eta < 0.5", 480, -240, 240);

    hist1d[57] = new TH1D("jetEginValue3D", "jetEginValue3D", 110, 0, 1.1);
    hist1d[58] = new TH1D("jetEtaPhiAngle3D", "jetEtaPhiAngle3D", 660, -0.2, 6.4);
    hist1d[59] = new TH1D("jetEtaTimAngle3D", "jetEtaTimAngle3D", 660, -0.2, 6.4);
    hist1d[60] = new TH1D("jetEtaPhiAngle2D", "jetEtaTimAngle2D", 660, -0.2, 6.4);
    hist1d[61] = new TH1D("jetEginValueSph", "jetEginValueSph", 150, 0.4, 1.1);


	//----- photons 100 - 249

    hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);
    hist1d[101] = new TH1D("genPhoParentPdgId", "genPhoParent:PdgId;PdgId",130,90,220);

    hist1d[102] = new TH1D("phoClstrR9", "phoClstrR9", 100, 0, 1);

    hist1d[103] = new TH1D("phoNumRecHits", "phoNumRecHits", 100, 0, 100);

    //hist1d[104] = new TH1D("phoEtaPhiAngle3D", "phoEtaPhiAngle3D", 660, -0.2, 6.4);
    //hist1d[105] = new TH1D("phoEtaTimAngle3D", "phoEtaTimAngle3D", 660, -0.2, 6.4);
    hist1d[106] = new TH1D("phoRotAngle", "phoRotAngle", 660, -0.2, 6.4);

    hist1d[107] = new TH1D("phoIsPixelSeed", "phoIsPixelSeed", 3, 0, 2);
    hist1d[108] = new TH1D("phoHadOverEM", "phoHadOverEM", 250, 0, 10);
    hist1d[109] = new TH1D("phoR9", "phoR9", 100, 0, 1);
    hist1d[110] = new TH1D("phoMipTotEnergy", "phoMipTotEnergy", 250, 0, 250 );
    hist1d[111] = new TH1D("phoEcalRHSumEtConeDR04", "phoEcalRHSumEtConeDR04", 750, 0, 750);
    hist1d[112] = new TH1D("phoHcalTwrSumEtConeDR04", "phoHcalTwrSumEtConeDR04", 750, 0, 750);

    hist1d[113] = new TH1D("phoEginValueSph", "phoEginValueSph", 150, 0.4, 1.1);
    hist1d[114] = new TH1D("phoEginSlopeSph", "phoEginSlopeSph", 20000, -100, 100);
    //hist1d[115] = new TH1D("phoEginValue3D", "phoEginValue3D", 110, 0, 1.1);

    hist1d[116] = new TH1D("pho_etprofile", "Photon rEta Time Profile Sph", cwdiv, -1*cwtrn, cwtrn);

    hist1d[117] = new TH1D("phoEginUnCorSlopeSph", "phoEginUnCorSlopeSph", 20000, -100, 100);
    hist1d[118] = new TH1D("pho_proFitEtavChi", "Profile Fit Eta v Chi2Prob Sph", cwdiv, -1*cwtrn, cwtrn );

    hist1d[119] = new TH1D("phoClTime", "phoTime", jtdiv, -1*jtran, jtran);
    hist1d[120] = new TH1D("phoSeedRhTime", "phoLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[121] = new TH1D("phoSeedTimeDiff", "phoLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[122] = new TH1D("phoClass", "phoClass", 4, 0, 4);

    hist1d[123] = new TH1D("phoSlopeVar1","phoSlopeVar1",200,0,20);
    hist1d[124] = new TH1D("phoSlopeErr1","phoSlopeErr1",100,0,10);
    hist1d[125] = new TH1D("phoSlopeVar2","phoSlopeVar2",200,0,20);
    hist1d[126] = new TH1D("phoSlopeErr2","phoSlopeErr2",100,0,10);

    hist1d[127] = new TH1D("phoEiganRotAngle","phoEiganRotAngle",70,0,7);

    hist1d[128] = new TH1D("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
    //hist1d[129] = new TH1D("phoGeoEginValueSph", "phoGeoEginValueSph", 150, 0.4, 1.1);
    hist1d[130] = new TH1D("phoGenDr", "phoGenDr", 150, 0, 0.15 );
    hist1d[131] = new TH1D("phoGenPdgId", "phoGenPdgId", 50, 0, 50 );
    hist1d[132] = new TH1D("phoType", "phoType( 1 pho; 2 llp; 3 oot; 4 cmb )", 6,0,5);

    hist1d[133] = new TH1D("phoGeoEginValueSph", "phoGeoEginValueSph", 150, 0.4, 1.1);
    hist1d[134] = new TH1D("phoGeoEginSlopeSph", "phoGeoEginSlopeSph", 20000, -100, 100);
    hist1d[135] = new TH1D("phoGeoEiganRotAngle","phoGeoEiganRotAngle",70,0,7);

	//------  genparticles 200 - 249

    hist1d[200] = new TH1D("genPdgId", "genPdgId;PdgId",220,0,220);
    hist1d[201] = new TH1D("nGenLLPPho", "nGenLLPPho", 20, 0, 20 );

	//------  cluster jet 250 - 299

    hist1d[250] = new TH1D("cljTime", "cljTime", jtdiv, -1*jtran, jtran);
    hist1d[251] = new TH1D("cljLeadRhTime", "cljLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[252] = new TH1D("cljClstLeadTimeDiff", "cljClstLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[253] = new TH1D("cljDrTime", "cljDrTime", jtdiv, -1*jtran, jtran);
    hist1d[254] = new TH1D("cljDiffPt", "cljDiffPt", 1000, 0, 10);
    hist1d[255] = new TH1D("cljdPhi", "cljdPhi", 70, -3.5, 3.5);
    hist1d[256] = new TH1D("cljdtmu", "cljdtmu", jdtdiv, -1*jdtran, jdtran);

    //------ ecal rechits 300 - 349
    hist1d[300] = new TH1D("rhRadius", "rhRadius", 300, 128, 131);

    hist1d[301] = new TH1D("ebRhTime", "ebRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[302] = new TH1D("ebRhEnergy", "ebRhEnergy", 1000, 0, 1000);
    hist1d[303] = new TH1D("eeRhTime", "eeRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[304] = new TH1D("eeRhEnergy", "eeRhEnergy", 1000, 0, 1000);
    hist1d[305] = new TH1D("ebRhkOOT", "ebRhkOOT", 3, 0, 1);
    hist1d[306] = new TH1D("eeRhkOOT", "eeRhkOOT", 3, 0, 1);

    //------  electrons 350 - 400

    hist1d[350] = new TH1D("eleClTime", "eleClTime", jtdiv, -1*jtran, jtran);
    hist1d[351] = new TH1D("eleSeedRhTime", "eleLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[352] = new TH1D("eleClSeedTimeDiff", "eleClLeadTimeDiff", jtdiv, -1*jtran, jtran);

	//-------- event vars 400 - 450
	
	hist1d[300] = new TH1D("nphoexclusions", " Number Pho (0) Excl (1) OOT (2) Excl (3) ", 4, 0, 4);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	//------- jets ( time ) 0-49 ------------------------------

	//hist2d[0]
    hist2d[1] = new TH2D("jetDrMuTime_pt", "jetDrMuTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[2] = new TH2D("jetDrMuTime_id", "jetDrMuTime_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
    hist2d[3] = new TH2D("jetDrMuTime_nhf", "jetDrMuTime_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[4] = new TH2D("jetDrMuTime_chf", "jetDrMuTime_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[5] = new TH2D("jetDrMuTime_nemf", "jetDrMuTime_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[6] = new TH2D("jetDrMuTime_cemf", "jetDrMuTime_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[7] = new TH2D("jetDrMuTime_muf", "jetDrMuTime_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[8] = new TH2D("jetDrMuTime_nhm", "jetDrMuTime_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
    hist2d[9] = new TH2D("jetDrMuTime_chm", "jetDrMuTime_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

    hist2d[10] = new TH2D("jetDrMuTime_jetDrMedTime", "jetDrMuTime_jetDrMedTime", jtdiv, -1*jtran, jtran, 200, -10, 10);
    hist2d[11] = new TH2D("jetDrMuTime_rms", "jetDrMuTime_rms", jtdiv, -1*jtran, jtran, 200, 0, 20);
    hist2d[12] = new TH2D("jetDrMuTime_err", "jetDrMuTime_err", jtdiv, -1*jtran, jtran, 300, 0, 3);

    hist2d[13] = new TH2D("jetDrMedTime_pt", "jetDrMedTime_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[14] = new TH2D("jetDrMedTime_id", "jetDrMedTime_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
    hist2d[15] = new TH2D("jetDrMedTime_nhf", "jetDrMedTime_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[16] = new TH2D("jetDrMedTime_chf", "jetDrMedTime_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[17] = new TH2D("jetDrMedTime_nemf", "jetDrMedTime_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[18] = new TH2D("jetDrMedTime_cemf", "jetDrMedTime_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[19] = new TH2D("jetDrMedTime_muf", "jetDrMedTime_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[20] = new TH2D("jetDrMedTime_nhm", "jetDrMedTime_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
    hist2d[21] = new TH2D("jetDrMedTime_chm", "jetDrMedTime_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

    hist2d[22] = new TH2D("jetDrDtMu_nJets", "jetDrDtMu_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);
    hist2d[23] = new TH2D("jetDrDtMed_nJets", "jetDrDtMed_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);

    hist2d[24] = new TH2D("jetDrMu_nrh", "jetDrMu_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);
    hist2d[25] = new TH2D("jetDrMedTime_nrh", "jetDrMedTime_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);

    hist2d[26] = new TH2D("jetDrDtMu_diffPt", "jetDrDtMu_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[27] = new TH2D("jetDrDtMu_htPct", "jetDrDtMu_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[28] = new TH2D("jetDrDtMu_dPhi", "jetDrDtMu_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);
    hist2d[29] = new TH2D("jetDrDtMed_diffPt", "jetDrDtMed_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[30] = new TH2D("jetDrDtMed_htPct", "jetDrDtMed_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[31] = new TH2D("jetDrDtMed_dPhi", "jetDrDtMed_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);

    //hist2d[32] = new TH2D("jetDrMuTime_sceta", "jetDrMuTime_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[33] = new TH2D("jetDrMuTime_scphi","jetDrMuTime_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[34] = new TH2D("jetDrMuTime_scenr", "jetDrMuTime_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    //hist2d[35] = new TH2D("jetDrMedTime_sceta", "jetDrMedTime_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[36] = new TH2D("jetDrMedTime_scphi","jetDrMedTime_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    //hist2d[37] = new TH2D("jetDrMedTime_scenr", "jetDrMedTime_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    hist2d[38] = new TH2D("nRhDrJet_nRhBcJet", "nRhDrJet_nRhBcJet", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);
    hist2d[39] = new TH2D("nRhDrJet_nRhScJet", "nRhDrJet_nRhScJet", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);

    hist2d[30] = new TH2D("jetScDtMu_effje", "jetScDtMu_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);
    hist2d[31] = new TH2D("jetDrDtMu_effje", "jetDrDtMu_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);


	//---jet id stuff 50 - 99 ---------------------------------------------------
	//hist2d[50]
	//hist2d[51]
	//hist2d[52]
    hist2d[53] = new TH2D("scemf_emf", "scemf_emf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[54] = new TH2D("bcemf_emf", "bcemf_emf", 110, 0, 1.1, 110, 0, 1.1);
    //hist2d[55] = new TH2D("dreme_eme", "dreme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[56] = new TH2D("sceme_eme", "sceme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[57] = new TH2D("bceme_eme", "bceme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[58] = new TH2D("sce_bce", "sce_bce", 500, 0, 500, 500, 0, 500);
    hist2d[59] = new TH2D("dremf_scemf", "dremf_scemf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[60] = new TH2D("scemf_bcemf", "scemf_bcemf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[61] = new TH2D("epaf_epf", "epaf_epf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[62] = new TH2D("epf_emf", "epf_emf", 110, 0, 1.1, 110, 0, 1.1);

	//--- jet cluster slope 100 - 149 --------------------------------------------------
    hist2d[100] = new TH2D("jetEta_Slope", "Jet Eta v Slope;Eta;Slope ps/cm", 120, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[101] = new TH2D("jetImpAngle_Slope", "Jet ImpactAngle v Slope;ImpactAngle;Slope ps/cm", 120, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[102] = new TH2D("jetImpAngle_Slope3D", "Jet ImpactAngle v Slope 3D;ImpactAngle;Slope", 150, 0, 1.5, sl3div, sl3min, sl3max);

    //hist2d[119] = new TH2D("jetSlope_RotSlope", "Jet Slope v rotated Slope;Slope;rotated Slope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    //hist2d[120] = new TH2D("jetSlope_DifSlope", "Jet Slope v dif w/ rotSlope;Slope;difSlope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    hist2d[121] = new TH2D("jetPhi_Slope", "Jet Phi v Slope;Phi;Slope ps/cm", 70, -3.5, 3.5, sldiv, slmin, slmax);

    hist2d[126] = new TH2D("jetEta_Slope3D", "Jet Eta v Slope 3D;Eta;Slope ps/cm", 60, -1.5, 1.5, sl3div, sl3min, sl3max);
    hist2d[127] = new TH2D("jetPhi_Slope3D", "Jet Phi v Slope 3D;Phi;Slope ps/cm", 140, -3.5, 3.5, sl3div, sl3min, sl3max);

	hist2d[128] = new TH2D("jetSphEgnxy", "jetSphEgn x,y;x;y", 200, -5, 5, 200, -5, 5);

	// jet gen stuff 150 - 199 --------------------------------------------------------------
    hist2d[150] = new TH2D("jetE_GenE", "Jet Energy v GenEnergy;JetEnergy;GenEnergy", 100, 0, 1000, 100, 0, 1000 );
    hist2d[151] = new TH2D("jetEGenERatio_GenTime", "Jet E/GenE v GenTime;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[152] = new TH2D("jetEMFrac_SCTime", "Jet EMFrac v SCTime;EMFrac;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[153] = new TH2D("jetGenRatio_GenTimePre", "Jet E/GenE v GenTime Pre;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[154] = new TH2D("jetGenTime_DrJetTime", "GenTime v DrJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[155] = new TH2D("jetEGenERatio_SCTimeDiff", "Jet E/GenE v JetSC GenJet TimeDif;E/GenE;SCTimeDif", 80, 0, 2, 300, 0, 30.0 );
    hist2d[156] = new TH2D("jetSCTime_DrTime", "JetSCTime v JetDrTime;JetSCTime;JetDrTime", 280, -15, 25, 280, -15, 25 );
    hist2d[157] = new TH2D("jetGenTime_SCJetTime", "GenTime v SCJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[158] = new TH2D("jetGenTime_GenEnergy", "GenTime v GenEnergy;GenTime;GenEnergy", 280, -15, 25, 100, 0, 1000 );
    hist2d[159] = new TH2D("jetGenjetDr_SCTimeDiff", "Jet GenJet Dr v SCJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[160] = new TH2D("jetGenjetDr_DRTimeDiff", "Jet GenJet Dr v DRJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[161] = new TH2D("jetGenTime_TOFcorr", "GenTime v TOFcorr;GenTime;TOFcorr", 250, 0, 25, 250, 0, 25 );
    hist2d[162] = new TH2D("jetEGenERatio_SCTime", "Jet E/GenE v SCTime;E/GenE;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[163] = new TH2D("jetEGenERatio_DRTime", "Jet E/GenE v DRTime;E/GenE;DRTime", 80, 0, 2, 40, -15, 25 );
    hist2d[164] = new TH2D("jetGenVar_SCTimeDiff", "Jet GenJet Var v SCJet GenJet TimeDif;Var;TimeDif", 270, -2, 25, 200, 0, 20.0 );
    hist2d[165] = new TH2D("jetGenPurity_SCTimeDiff", "Jet GenJet Purity v SCJet GenJet TimeDif;Purity;TimeDif", 100, 0, 1, 200, 0, 20.0 );
    hist2d[166] = new TH2D("jetGenPurity_GenJet_ar", "Jet GenJet Purity v GenJet Var;Purity;Var", 100, 0, 1, 270, -2, 25 );
    hist2d[167] = new TH2D("jetGenVar_GenJetNKids", "Jet GenJet Var v GenJet nKids;Var;nKids", 270, -2, 25, 100, 0, 100 );
    hist2d[168] = new TH2D("jetGenPurity_GenJetNKids", "Jet GenJet Purity v GenJet nKids;Purity;nKids", 100, 0, 1, 100, 0, 100 );
    hist2d[169] = new TH2D("jetGenSCTimeDiff_DrMatchJet", "genJet SCTimeDiff v DrMatchJet;SCTimeDiff;DrMatchJet", 300, 0, 30, 320, 0, 3.2 );
    hist2d[170] = new TH2D("jetGenTime_JetEMFrac", "GenTime v JetEMFrac;GenTime;JetEMFrac", 280, -15, 25, 150, 0, 1.5 );

    hist2d[171] = new TH2D("jetEmFrac_GenjetDr", "Jet emFrac v GenJet Dr;emFrac;jetGenDr", 40, 0, 1, 200, 0, 0.5 );
    hist2d[172] = new TH2D("jetEGenERatio_GenjetDr", "Jet E/GenE v GenJet Dr;E/GenE;jetGenDr", 80, 0, 2, 200, 0, 0.5 );

    hist2d[173] = new TH2D("jet3DEgnxy", "jet3DEgn x,y;x;y", 200, -5, 5, 200, -5, 5);
    hist2d[174] = new TH2D("jet3DEgnxz", "jet3DEgn x,z;x;z", 200, -5, 5, 200, -5, 5);

	//--- Photons 200 - 349 -------------------------------------------
    hist2d[200] = new TH2D( "phoEta_Slope", "Photon Eta Vs Slope;Eta;Slope [ps/cm]", 60, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[201] = new TH2D( "phoEta_EV2D", "Photon Eta Vs EiganValue2D;Eta;EiganValue",60, -1.5, 1.5, 150, 0.4, 1.1);
    hist2d[202] = new TH2D( "phoEta_UnCorSlope", "Photon Eta Vs Slope;Eta;UnCorSlope [ps/cm]", 60, -1.5, 1.5, sldiv, slmin, slmax);
    //hist2d[203]

    hist2d[204] = new TH2D("phoID_LLP", "phoIDvLLP;isPho;isLLP", 3, -1, 2, 2, 0, 2);
    hist2d[205] = new TH2D("phoIDL_T", "phoIDLvT;isLoose;isTight", 2, 0, 2, 2, 0, 2);

    hist2d[206] = new TH2D("pho_tmap_rot", "Photon Time Map Rotated;rEta;rPhi", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[207] = new TH2D("pho_occmap_rot", "Photon Occ Map Rotated;rEta;rPhi", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[208] = new TH2D("pho_tmap", "Photon Time Map;iEta;iPhi", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[209] = new TH2D("pho_occmap", "Photon Occ Map;iEta;iPhi", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);

    hist2d[210] = new TH2D("phoSphEgnxy", "phoSphEgn x,y;x;y", 200, -5, 5, 200, -5, 5);
    //hist2d[211] = new TH2D("pho3DEgnxy", "pho3DEgn x,y;x;y", 200, -5, 5, 200, -5, 5);
    //hist2d[212] = new TH2D("pho3DEgnxz", "pho3DEgn x,z;x;z", 200, -5, 5, 200, -5, 5);
    //hist2d[213] = new TH2D("phoNRH_ClR9vID","Photon nClRecHits v clR9 v phoId;nRecHits;ClusterR9;PhotonID(fake1,loose2,tight3)",100,0,100,100,0,100);

    hist2d[214] = new TH2D("pho_ptwtmap", "Photon Phi(x) Time(y) Wt Map Sph;rPhi;Time", cwdiv, -1*cwtrn, cwtrn, clsphdiv, -1*clsphtrn, clsphtrn);
    hist2d[215] = new TH2D("pho_etwtmap", "Photon Eta(x) Time(y) Wt Map Sph;rEta;Time", cwdiv, -1*cwtrn, cwtrn, clsphdiv, -1*clsphtrn, clsphtrn);

    //hist2d[216] = new TH2D("phoNRH_ClR9_fake","Photon Fake nClRecHits v clR9;nRecHits;ClusterR9",200,0,200,100,0,1);
    //hist2d[217] = new TH2D("phoNRH_ClR9_loose","Photon Loose nClRecHits v clR9;nRecHits;ClusterR9",200,0,200,100,0,1);
    hist2d[218] = new TH2D("phoNRH_ClR9","Photon nClRecHits v clR9;nRecHits;ClusterR9",200,0,200,100,0,1);

    hist2d[219] = new TH2D("phoSlopeVar1_Err1","phoSlope Var1 vs Err1;Varience;Error",250,0,25,500,0,5);
    hist2d[220] = new TH2D("phoSlope_Var1","phoSlope Slope vs Var1;Slope;Varience",sldiv, slmin, slmax,250,0,25);
    hist2d[221] = new TH2D("phoSlope_Err1","phoSlope Slope vs Err1;Slope;Error",sldiv, slmin, slmax,500,0,5);

    hist2d[222] = new TH2D("phoSlErr_ClR9","Photon Slope Err v clR9;error;ClusterR9",500,0,5,100,0,1);
    hist2d[223] = new TH2D("phoR9_ClR9","phoR9_ClR9;R9;clR9",100,0,1,100,0,1);

    //hist2d[224] = new TH2D("phoNRH_R9_fake","Photon Fake nClRecHits v R9;nRecHits;R9",200,0,200,100,0,1);
    hist2d[225] = new TH2D("phoNRH_R9","Photon nClRecHits v R9;nRecHits;R9",200,0,200,100,0,1);

    //hist2d[226] = new TH2D("phoEignAngleUncr_eta","Photon EignAngleNoCor v Eta;Eta;Angle",60,-1.5,1.5,70,0,7);
    //hist2d[227] = new TH2D("phoEignAngleUncr_phi","Photon EignAngleNoCor v Phi;Phi;Angle",140,-3.5,3.5,70,0,7);
    //hist2d[228] = new TH2D("phoEignAngleUncr_value","Photon EignAngleNoCor v eiganValue;Angle;Value",70,0,7,50,0.5,1);
    //hist2d[229] = new TH2D("phoEignAngleUncr_slope","Photon EignAngleNoCor v eiganSlope;Angle;Slope",70,0,7,sldiv, slmin, slmax);

	// ---   pho v/s/a/c comps  ----

    //hist2d[] = new TH2D("phoGeoEignAngleUncr_eta","Photon GeoEignAngleNoCor v Eta;Eta;Angle",60,-1.5,1.5,70,0,7);
    //hist2d[] = new TH2D("phoGeoEignAngleUncr_phi","Photon GeoEignAngleNoCor v Phi;Phi;Angle",140,-3.5,3.5,70,0,7);
    //hist2d[] = new TH2D("phoGeoEignAngleUncr_value","Photon GeoEignAngleNoCor v eiganValue;Angle;Value",70,0,7,50,0.5,1);
    //hist2d[] = new TH2D("phoGeoEignAngleUncr_slope","Photon GeoEignAngleNoCor v eiganSlope;Angle;Slope",70,0,7,sldiv, slmin, slmax);

	// eigan plots value/slope/angle/clr9

    hist2d[235] = new TH2D("phoRotAngle_value","Photon rotAngle v Value;rotAngle;Value",70,0,7,50,0.5,1);
    hist2d[234] = new TH2D("phoRotAngle_rotslope","Photon rotAngle v Slope;rotAngle;Slope",70,0,7,sldiv, slmin, slmax);
    hist2d[230] = new TH2D("phoClR9_value","Photon ClR9 v Value;ClR9;Value",100,0,1,50,0.5,1);
    hist2d[231] = new TH2D("phoClR9_slope","Photon ClR9 v Slope;ClR9;Slope",100,0,1,sldiv, slmin, slmax);
    hist2d[232] = new TH2D("phoClR9_rotAngle","Photon ClR9 v rotAngle;ClR9;rotAngle",100,0,1,70,0,7);
    hist2d[233] = new TH2D("phoValue_slope","Photon Value v Slope;Value;Slope",50,0.5,1,sldiv, slmin, slmax);

    hist2d[236] = new TH2D("phoGeoValue_angle","Photon geoValue v angle;geoValue;angle",50,0.5,1,70,0,7);
    hist2d[237] = new TH2D("phoGeoSlope_angle","Photon geoSlope v angle;geoSlope;angle",sldiv, slmin, slmax,70,0,7);
    hist2d[238] = new TH2D("phoGeoSlope_value","Photon geoSlope v value;geoSlope;value",sldiv, slmin, slmax,50,0.5,1);

    hist2d[239] = new TH2D("phoGeoRotAngle_Angle","Photon geoRotAngle v Angle;geoRotAngle;rotAngle",70,0,7,70,0,7);
    hist2d[242] = new TH2D("phoGeoRotAngle_Value","Photon geoRotAngle v Value;geoRotAngle;Value",70,0,7,50,0.5,1);
    hist2d[241] = new TH2D("phoGeoRotAngle_slope","Photon geoRotAngle v Slope;geoRotAngle;slope",70,0,7,sldiv, slmin, slmax);
    hist2d[243] = new TH2D("phoGeoValue_value","Photon geoValue v value;geoValue;value",50,0.5,1,50,0.5,1);
    hist2d[244] = new TH2D("phoGeoValue_slope","Photon geoValue v slope;geoValue;slope",50,0.5,1,sldiv, slmin, slmax);
    hist2d[240] = new TH2D("phoGeoSlope_slope","Photon geoSlope v slope;geoSlope;slope",sldiv, slmin, slmax,sldiv, slmin, slmax);

    hist2d[245] = new TH2D("phoGeoRotAngle_geoValue","Photon geoRotAngle v geoValue;geoRotAngle;geoValue",70,0,7,50,0.5,1);
    hist2d[246] = new TH2D("phoGeoRotAngle_geoRotslope","Photon geoRotAngle v geoSlope;geoRotAngle;geoSlope",70,0,7,sldiv, slmin, slmax);
    hist2d[247] = new TH2D("phoClR9_geoValue","Photon ClR9 v geoValue;ClR9;geoValue",100,0,1,50,0.5,1);
    hist2d[248] = new TH2D("phoClR9_geoSlope","Photon ClR9 v geoSlope;ClR9;geoSlope",100,0,1,sldiv, slmin, slmax);
    hist2d[249] = new TH2D("phoClR9_georotAngle","Photon ClR9 v geoRotAngle;ClR9;geoRotAngle",100,0,1,70,0,7);
    hist2d[250] = new TH2D("phoGeoValue_geoSlope","Photon geoValue v geoSlope;geoValue;geoSlope",50,0.5,1,sldiv, slmin, slmax);


	//--- rechit collections 350 - 399 -------------------------------------------------
    hist2d[350] = new TH2D("ebRhTime_Energy", "ebRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );
    hist2d[351] = new TH2D("eeRhTime_Energy", "eeRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );

    //------------------------------------------------------------------------------------------
    //------ 3D Hists --------------------------------------------------------------------------

    //hist3d[0] = new TH3D("phoNRH_ClR9_phoID","Photon nClRecHits v clR9 v phoId;nRecHits;ClusterR9;PhotonID(fake1,loose2,tight3)",200,0,200,100,0,1,8,0,4);

	//------------------------------------------------------------------------------------
    // Cluster maps -----------------------------------------------------------------------
	nMaps = 0;
	for(int it=0; it<nEBEEMaps; it++){
		fMap[it] = false;
		string label(";iEta;iPhi");
        string stt1("ebeeMapPhoCluster_"+std::to_string(it));
        ebeeMapP[it] = new TH2D( stt1.c_str(), (stt1+label).c_str(), 361, -90, 90, 721, 0, 360);
        string stt2("ebeeMapPhoClusterTime_"+std::to_string(it));
        ebeeMapT[it] = new TH2D( stt2.c_str(), (stt2+label).c_str(), 361, -90, 90, 721, 0, 360);
		string stt3("ebeeMapPhoClusterRes_"+std::to_string(it));
        ebeeMapR[it] = new TH2D( stt3.c_str(), (stt3+label).c_str(), 361, -90, 90, 721, 0, 360);
	}//<<>>for(int it=0; it<nEBEEMaps; it++)

}//<<>>void makehists::initHists()

// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                auto indir = "LLPGamma/llpga_GMSB_AOD_v57/"; //argv[1];
                //auto indir = "LLPGamma/llpga_GJets_AOD_v57/";

                auto infilename = "llpgana_mc_AODSIM_GMSB_AOD_v57_Full.txt"; //argv[2];
                //auto infilename = "llpgana_mc_AODSIM_GJets_AOD_v57_Full.txt";

				auto outfilename = "lpgana_mc_test2.root";
                ///auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v53_V12_3th_notLLP_Tight_Excluded_Hists.root"; //argv[3];
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v53_V22_10th_Loose_tEB_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v53_V22_10th_Loose_tEB_noPeak_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v53_V22_10th_Loose_tEB_isLLP_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v53_V22_3th_Loose_tEB_notLLP_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v53_V22_10th_Loose_tEB_ClR9r68_Hists.root";
            	///auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v53_V12_3th_notLLP_Tight_Excluded_Hists.root"; //argv[3];
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v55_V23_5th_Loose_tEB_notCmb_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v53_V22_3th_Jets2_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v57_V24_20th_Loose_tEB_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v57_V24_20th_Loose_tEB_noPeak_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v57_V24_20th_Loose_tEB_isOOT_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v57_V24_20th_Loose_tEB_ClR9r26_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v57_V24_50th_Loose_tEB_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v57_V24_50th_Loose_tEB_ClR9r26_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v57_V24_50th_Loose_tEB_noPeak_Hists.root";

				int pct = 50;
				makehists base;
                base.llpgana_hist_maker( indir, infilename, outfilename, pct );
    //}
    return 1;

}//<<>>int main ( int argc, char *argv[] )



