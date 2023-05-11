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
#include "llpgana_hist_rebase_v18.hh"

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
const auto dltAngle (CFlt x, CFlt y){auto dp(x-y);if(dp>=2*PI) dp-=2*PI; else if(dp<=-2*PI) dp+=2*PI; return dp;}
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
		auto width = 2.0;
        auto high = mean + width*stdv;
        auto low = mean - width*stdv;
        //std::cout << " - Profile: m " << mean << " s " << stdv << " h " << high << " l " << low << " n " << norm << std::endl;
        if( std::abs(stdv) < 5.0 && std::abs(norm) > 0 ){
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
            if( error < 2.0 && std::abs(fmean) < 10.0 ){
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
                        for(auto ix : x ){ if( (maxphi-ix) > 179 ) sum+=(ix+360); else sum+=ix; }
                        auto rslt = sum/x.size();
                        if( rslt > 359 ) rslt-=360;
                        return rslt;
                    }//<<>> const auto meanIPhi(CVFlt x)

const auto meanIPhi  (CVFlt x, CVFlt wv){
                        float wt(0.0), sum(0.0); int it(0);
                        auto maxphi = max(x);
                        for(auto ix : x ){ if( (maxphi-ix) > 179 ) sum+=(ix+360)*wv[it]; else sum+=ix*wv[it]; wt+=wv[it]; it++; }
                        auto rslt = sum/wt;
                        if( rslt > 359 ) rslt-=360;
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

vector<float> getRhGrpEigen( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts ){
// ieipt

    vector<float> results;

    auto mean_x = mean( xs, wts );
    auto mean_y = mean( ys, wts );
    auto mean_z = mean( zs, wts );
    auto swts = accum( wts );
    auto var_x = var( xs, mean_x, wts, swts );
    auto var_y = var( ys, mean_y, wts, swts );
    auto var_z = var( zs, mean_z, wts, swts );
    auto var_xy = cvar( xs, mean_x, ys, mean_y, wts, swts );
    auto var_xz = cvar( xs, mean_x, zs, mean_z, wts, swts );
    auto var_yz = cvar( ys, mean_y, zs, mean_z, wts, swts );

    //TMatrixDSym mat(3,3);
    double array[] = { var_x, var_xy, var_xz, var_xy, var_y, var_yz, var_xz, var_yz, var_z };
    TMatrixDSym mat(3,array);
    //std::cout << "Input matrix created" << std::endl;
    //mat.Print();

    TMatrixDSymEigen eigen(mat);
    const TVectorD& eginVal = eigen.GetEigenValues();
    const TMatrixD& eginVec = eigen.GetEigenVectors();
    TVectorD evaules(eginVal);

    int index;
    int index2;
    int index3;
    if( (eginVal[0] >= eginVal[1]) && (eginVal[0] >= eginVal[2]) ){
        index = 0;
        if( eginVal[1] >= eginVal[2] ){ index2 = 1; index3 = 2; }
		else { index2 = 2; index3 = 1; }
        //if( (eginVal[0] == eginVal[1]) && (eginVal[0] == eginVal[2]) ); // std::cout << " -- rhGrp is Spherical" << std::endl;
    } else if( eginVal[1] >= eginVal[2] ) {
        index = 1;
        if( eginVal[0] >= eginVal[2] ){ index2 = 0; index3 = 2; }
        else { index2 = 2; index3 = 0; }
        //if( eginVal[1] == eginVal[2] ) ; //std::cout << " -- rhGrp is a Flatend Sphere" << std::endl;
    } else { 
		index = 2;
		if( eginVal[0] >= eginVal[1] ){ index2 = 0; index3 = 1; }
		else { index2 = 1; index3 = 0; }
	} //<<>>if( (eginVal[0] >= eginVal[1]) && (eginVal[0] >= eginVal[2]) )
    //std::cout << "Largest Eigin Value Found" << std::endl; 

    auto ex  = eginVec(index,0);
    auto ey  = eginVec(index,1);
    auto ez  = eginVec(index,2);
    auto ev  = eginVal[index]/(eginVal[index]+eginVal[index2]); // z vs c+z oval
    auto ev2 = eginVal[index2]/(eginVal[index2]+hypo(eginVal[index], eginVal[index3])); // c vs c+zt oval
	auto ev3 = eginVal[index3]/(eginVal[index3]+hypo(eginVal[index], eginVal[index2])); // t vs t+zc oval
    auto ev4 = eginVal[index2]/(eginVal[index3]+eginVal[index2]); // c vs c+t
    auto ev5 = eginVal[index]/(eginVal[index]+hypo(eginVal[index2], eginVal[index3])); // z vs z+ct oval
    //std::cout << "Eigin Angles Found" << std::endl;
    
    results.push_back(ex);//0
    results.push_back(ey);//1
    results.push_back(ez);//2

    results.push_back(ev);//3
    results.push_back(ev2);//4
    results.push_back(ev3);//5
    results.push_back(ev4);//6
    results.push_back(ev5);//7

    results.push_back(eginVal[index]);//8
    results.push_back(eginVal[index2]);//9
    results.push_back(eginVal[index3]);//10
    
    return results;
}//<<>>vector<float> LLPgammaAnalyzer_AOD::getRhGrpEigen3D( vector<float> xs, vector<float> ys, vector<float> zs, vector<float> wts )

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
	void makeEBEEMaps( int phoit );
	void makeEBEEMaps( vector<unsigned int> rhcol );
	vector<float> getLeadTofRhTime( vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
	vector<float> getRhGrpEigen_sph( vector<float> times, vector<uInt> rechitids );

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

void makehists::makeEBEEMaps( vector<unsigned int> rhcol ){

    if( DEBUG ) std::cout << " -- looping over ebee maps:" << std::endl;
	if( nMaps >= nEBEEMaps ) return;
    bool fill(false);
	auto nrh = rhcol.size();
    if( nMaps < 6 ){ if( nrh < 25 ) fill = true; }
    else if( nMaps < 12 ){ if( nrh < 35 ) fill = true; }
    else if( nMaps < 18 ){ if( nrh < 45 ) fill = true; }
    else if( nMaps < 24 ){ if( nrh < 55 ) fill = true; }
    else if( nMaps < 30 ){ if( nrh < 65 ) fill = true; }
    else { if( nrh >= 65 ) fill = true; }
    if( DEBUG ) std::cout << " -- start fill of ebee maps:" << std::endl;
    if(fill){
        if( DEBUG ) std::cout << " -- Filling ebeeMapT : " << nMaps << " with nRH : " << nRecHits << std::endl;
        vector<float> times;
        for( int idx = 0; idx < nrh; idx++ ) times.push_back((*rhTime)[getRhIdx(rhcol[idx])]);
		// sum energy
        auto mtime = mean(times);
        for( int idx = 0; idx < nrh; idx++ ){
            const auto rhIDX = getRhIdx(rhcol[idx]);
            auto idinfo = DetIDMap[rhcol[idx]];
            const auto rhEtaPos = idinfo.i2;//recHitPos.eta();
            const auto rhPhiPos = idinfo.i1;//recHitPos.phi();
            //auto res = 1/(sq2(3.64/(*rhEnergy)[rhIDX])+0.18);
			//auto logwt = std::max(0.0, 4.2 + log(rhenergy/sumRhEn));
            ebeeMapP[nMaps]->Fill( rhEtaPos, rhPhiPos,1 );
            ebeeMapP[nMaps]->Fill( nMaps, 1,1 );
            ebeeMapT[nMaps]->Fill( rhEtaPos, rhPhiPos,100+(*rhTime)[rhIDX] - mtime );
            ebeeMapR[nMaps]->Fill( rhEtaPos, rhPhiPos );
        }//<<>>for( idx = 0; idx < nRecHits; idx++ )
        nMaps++;
    }//<<>>if(fill) 

}//<<>>void makehists::makeEBEEMaps( vector<unsigned int> )


void makehists::makeEBEEMaps( int phoit ){

    if( DEBUG ) std::cout << " - looping photons : making ebee maps" << std::endl;
	auto isEB = (*phoIsEB)[phoit];
	auto rhcol = (*phoRhIds)[phoit];
    if( isEB ) makeEBEEMaps(rhcol);
    if( DEBUG ) std::cout << " - looping photons : Finished making ebee maps" << std::endl;
	return;

}//int makehists::makeEBEEMaps( int phoit )

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

vector<float> makehists::getRhGrpEigen_sph( vector<float> times, vector<uInt> rechitids ){

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2

    float N(3.64);
    float C(0.3000);

    vector<float> eg2wts;
    vector<float> rhetas, rhphis;;
    vector<float> xs, ys, zs;
    vector<float> rhxs, rhys, rhzs;
    vector<float> rhtimes;
    vector<float> grhtimes;
    vector<float> angles;
    vector<float> zcAngles;
    vector<float> invtresvec, rhinvtresvec;
	vector<float> rhlogwtvec, rhtresvec;
	vector<float> logwtvec, tresvec;
    vector<float> emptyReturn(8,-9);

	// --------- prepare inputs for eigan calcs --------------------------
	//std::cout << " getRhGrpEigen_sph 1, ";

    auto nRecHits = rechitids.size();
	if( nRecHits < 16 ) return emptyReturn;
    float sumRhEn(0);
    for ( auto id : rechitids ){ sumRhEn +=  (*rhEnergy)[getRhIdx(id)]; }
	if( sumRhEn <= 0 ) return emptyReturn;
    for( uInt it(0); it < nRecHits; it++ ){

        const auto rhIDX = getRhIdx(rechitids[it]);
		auto idinfo = DetIDMap[rechitids[it]];
        auto isEB = idinfo.ecal == ECAL::EB;
		if( isEB ) hist1d[123]->Fill(times[it]);
		//std::cout << "In getRhGrpEigen_sph w/ idx : " << rhIDX << std::endl;
		if( rhIDX == -1 ){ return emptyReturn; std::cout << " -- Bad idx !!!!! -- In getRhGrpEigen_sph ---- " << std::endl; }
        if( isEB ){
            const auto rhEtaPos = idinfo.i2;//recHitPos.ieta();
            rhetas.push_back((rhEtaPos>0)?rhEtaPos+84.5:rhEtaPos+85.5);
            const auto rhPhiPos = idinfo.i1;//recHitPos.iphi();
            rhphis.push_back(rhPhiPos-0.5);
            const auto rhXPos = (*rhPosX)[rhIDX];
            rhxs.push_back(rhXPos);
            const auto rhYPos = (*rhPosY)[rhIDX];
            rhys.push_back(rhYPos);
            const auto rhZPos = (*rhPosZ)[rhIDX];
            rhzs.push_back(rhZPos);
            rhtimes.push_back(times[it]);
            auto rhenergy = (*rhEnergy)[rhIDX];
            auto resolution = sq2(N/rhenergy)+2*C*C;
            auto logwt = std::max(0.0, 4.2 + log(rhenergy/sumRhEn));// cut at rh energy < 1.5% of cluster
			rhlogwtvec.push_back(logwt);
			rhinvtresvec.push_back(1/resolution);
			rhtresvec.push_back(resolution);
            //std::cout << "In getRhGrpEigen_sph w/ rheta " << rhEtaPos << " : rhphi " << rhPhiPos << " : rht " << times[it] << std::endl;      
        } else { 
			return emptyReturn; //std::cout << "In getRhGrpEigen_sph : NOT EB !!!!!!" << std::endl; }
        }//<<>> else - if( idinfo.ecal == ECAL::EB )
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )
	if( rhtimes.size() < 9 ) return emptyReturn;

    bool uselog(true);

    auto rhewgt = (uselog)?rhlogwtvec:rhinvtresvec;
    auto tmtime = mean(rhtimes,rhewgt);
	for( uInt it(0); it < rhtimes.size(); it++ ){ 
		auto rht = rhtimes[it];
		auto goodDifTime = (rhtimes[it]-tmtime) < 10.0;
		auto isInTime = ( rhtimes[it] > -50.0 ) && ( rhtimes[it] < 50.0 );
		if( isInTime && goodDifTime ){ 
			grhtimes.push_back(rhtimes[it]);
			xs.push_back(rhxs[it]);
            ys.push_back(rhys[it]);
            zs.push_back(rhzs[it]);		
			logwtvec.push_back(rhlogwtvec[it]);
			tresvec.push_back(rhtresvec[it]);
			invtresvec.push_back(rhinvtresvec[it]);
		}//<<>>if( isInTime && goodDifTime )
	}//<<>>for( auto rht : rhtimes )
	if( grhtimes.size() < 9 ) return emptyReturn;

    //std::cout << "2, ";
    vector<float> letas;
    vector<float> lphis;
    auto meta = mean(rhetas,rhewgt);
    auto mphi = meanIPhi(rhphis,rhewgt);
    //auto meta = mean(rhetas);
    //auto mphi = meanIPhi(rhphis);
    for( uInt it(0); it < rhetas.size(); it++ ){
        float leta = rhetas[it]-meta;
        letas.push_back(leta);
        float lphi = dltIPhi(rhphis[it],mphi);
        lphis.push_back(lphi);
        float angle = getAngle( leta, lphi );
        angles.push_back(angle);
    }//<<>>for( uInt it(0); it < etas.size(); it++ )

    auto ewgt = (uselog)?logwtvec:invtresvec;
    auto mtime = mean(grhtimes,ewgt);
	auto mx = mean(xs,ewgt);
    auto my = mean(ys,ewgt);
    auto mz = mean(zs,ewgt);
	auto mr = hypo(mx,my);
	auto ma = std::atan2(my,mx);

    //std::cout << "In getRhGrpEigen_sph w/ meta " << meta << " : mphi " << mphi << " : mt " << mtime << std::endl;
    vector<float> lzs;
    vector<float> lcs;
    vector<float> lts;
    vector<float> ltds;
    vector<float> nolts;
    vector<float> invltres;
	float minDr(100.0);
	float vslz(0.0);
	float vslc(0.0);
    float vs3lz(0.0);
    float vs3lc(0.0);
    float vs3lt(0.0);
    for( uInt it(0); it < grhtimes.size(); it++ ){

        float ltim = grhtimes[it]-mtime;
        lts.push_back(ltim);
		auto ltd = ltim*SOL;
        ltds.push_back(ltd);
		invltres.push_back(ewgt[it]);
		auto lt2reswt = ltim*ltim*ewgt[it];
        eg2wts.push_back(lt2reswt);
        nolts.push_back(0.0);

		float lz = zs[it]-mz;
        lzs.push_back(lz);
		float lc = mr*(std::atan2(ys[it],xs[it])-ma);
        lcs.push_back(lc);
		//float zcAngle = getAngle( lz, lc );
        float zcAngle = std::atan2(lc,lz);
		zcAngles.push_back(zcAngle);

		auto dr = hypo(lz,lc);
		auto drt = hypo(ltd,dr);
		if( dr < minDr ) minDr = dr;
		// do dr calc in ieta,iphi cross check

		auto sqrtreswt = (uselog)?ewgt[it]:std::sqrt(ewgt[it]);
		auto ltdsqrtreswt = std::abs(ltd)*sqrtreswt;
        vs3lz += lz*ltdsqrtreswt/drt;
        vs3lc += lc*ltdsqrtreswt/drt;
        vs3lt += ltd*sqrtreswt/drt;
        vslz += lz*ltdsqrtreswt/dr;
        vslc += lc*ltdsqrtreswt/dr;
		//if( (std::abs(lz) < minDl) && (std::abs(lc) < minDl) ) ewgt[it] = 0;

    }//<<>>for( uInt it(0); it < grhtimes.size(); it++ )

	// --------------  get eigan values and vectors ---------------------------------------
    //std::cout << "3, ";

    auto eigens =  getRhGrpEigen( zcAngles, eg2wts );//0 x, 1 y, 2 values
	auto d2dot = eigens[0]*vslz + eigens[1]*vslc;
	if( d2dot < 0 ){ eigens[0] *= -1; eigens[1] *= -1; }
    auto eigens3 =  getRhGrpEigen( lzs, lcs, ltds, invltres );//0 x, 1 y, 2 values
    auto d3dot = eigens3[0]*vs3lz + eigens3[1]*vs3lc + eigens3[2]*vs3lt;
    if( d3dot < 0 ){ eigens3[0] *= -1; eigens3[1] *= -1; eigens3[2] *= -1; }
    auto geoeigens =  getRhGrpEigen( zcAngles, invltres );//0 x, 1 y, 2 values
    auto geoddot = geoeigens[0]*vslz + geoeigens[1]*vslc;
    if( geoddot < 0 ){ geoeigens[0] *= -1; geoeigens[1] *= -1; }
    ////auto geoeigens3 =  getRhGrpEigen( lzs, lcs, nolts, invltres );//0 x, 1 y, 2 values

	// --------------  get eigan vector angles ------------------------------------------ 
    //std::cout << "4, ";

    float rotangle = getAngle(eigens[0], eigens[1]);
    float e2sin = std::sin(rotangle); //eigens[1];
    float e2cos = std::cos(rotangle); //eigens[0];
    //float rot3angle = getAngle(eigens3[0], eigens3[1]);
    float rot3angle = getAngle(eigens3[0], eigens3[1]);
    float e3sin = std::sin(rot3angle);
    float e3cos = std::cos(rot3angle);
    float rotgangle = getAngle(geoeigens[0], geoeigens[1]);
    float egsin = std::sin(rotgangle);
    float egcos = std::cos(rotgangle);

    // -----------------------------------------
    // finding nemo ( slope )
    // -----------------------------------------
    //std::cout << "6, ";

    auto nWts = invltres.size();

    vector<float> xs1;
    vector<float> xs3;
    vector<float> xsgeo;
    vector<float> slvars1;
    vector<float> slvars3;
    vector<float> slvarsgeo;
    float xsum1(0.0);
    float xsum3(0.0);
    float xsumgeo(0.0);

    auto dsxcor = e2cos*(2.2) - e2sin*(2.2);
    auto d3xcor = e3cos*(2.2) - e3sin*(2.2);
    auto dgxcor = (geoeigens[0])*(2.2) - (geoeigens[1])*(2.2);
    auto xscorvar = sq2(dsxcor)/12;
    auto x3corvar = sq2(d3xcor)/12;
    auto xgcorvar = sq2(dgxcor)/12;

	// for pairs method
    vector<float> plzs;
    vector<float> plcs;
    //vector<float> plws;
    vector<float> pl3zs;
    vector<float> pl3cs;
    vector<float> plgzs;
    vector<float> plgcs;

    //std::cout << "7, ";
    for( uInt it(0); it < lzs.size(); it++ ){

        auto xscor = e2cos*lzs[it] - e2sin*lcs[it];
        auto yscor = e2sin*lzs[it] + e2cos*lcs[it];
        auto x3cor = e3cos*lzs[it] - e3sin*lcs[it];
        auto y3cor = e3sin*lzs[it] + e3cos*lcs[it];
        auto xgcor = egcos*lzs[it] - egsin*lcs[it];
        auto ygcor = egsin*lzs[it] + egcos*lcs[it];

		// for pairs method && histogram maps
		plzs.push_back(xscor);
		plcs.push_back(yscor);
        pl3zs.push_back(x3cor);
        pl3cs.push_back(y3cor);
        plgzs.push_back(xgcor);
        plgcs.push_back(ygcor);

        //if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << letas[it] << " : lphi " << lphis[it]
        //                        << " : xsor " << xscor << " : ycor " << yscor << " : dt " << eg2wts[it] << std::endl;
        //if( false ) std::cout << "In getRhGrpEigen_sph w/2 leta " << letas[it] << " : lphi " << lphis[it]
        //                        << " : sxcor " << x3cor << " : sycor " << y3cor << " : dt " << eg2wts[it] << std::endl;

		// calc slope info
        auto sl1 = (lts[it])/(xscor);//*slopeCorr;
        auto sl3 = (lts[it])/(x3cor);//*slopeCorr;
        auto slg = (lts[it])/(xgcor);//*slopeCorr;
        xs1.push_back(sl1);
        xs3.push_back(sl3);
        xsgeo.push_back(slg);
		//slvars1.push_back(1/((notresvec[it]+tottresvec+sq2(sl1)*xscorvar*(1.0+(1.0/nWts)))/sq2(xscor)));
		slvars1.push_back(invltres[it]);
        //slvars3.push_back(1/((notresvec[it]+tottresvec+sq2(sl3)*x3corvar*(1.0+(1.0/nWts)))/sq2(x3cor)));
		slvars3.push_back(invltres[it]);
        //slvarsgeo.push_back(1/((notresvec[it]+tottresvec+sq2(slg)*xgcorvar*(1.0+(1.0/nWts)))/sq2(xgcor)));
		slvarsgeo.push_back(invltres[it]);
		xsum1 += sl1*slvars1[it];
        xsum3 += sl3*slvars3[it];
        xsumgeo += slg*slvarsgeo[it];

    }//<<>>for( uInt it(0); it < wts.size(); it++ )

//===================================================================
    //std::cout << "8, ";

	vector<float> pszs1;
    vector<float> plsrs;
    float plssum(0.0);
    vector<float> p3zs1;
    vector<float> pl3rs;
    float pl3sum(0.0);
    vector<float> pgzs1;
    vector<float> plgrs;
    float plgsum(0.0);

    //=============================================================
    // pairs method set up
    //============================================================

    for( uInt it1(0); it1 < plzs.size(); it1++ ){
        for( uInt it2(it1); it2 < plzs.size(); it2++ ){

			auto minDz = 3.3;
            auto plz = plzs[it1]-plzs[it2];
			if( std::abs(plz) < minDz ) plz = 999;
            auto pl3z = pl3zs[it1]-pl3zs[it2];
            if( std::abs(pl3z) < minDz ) pl3z = 999;
            auto plgz = plgzs[it1]-plgzs[it2];
            if( std::abs(plgz) < minDz ) plgz = 999;
        	//auto plc = plcs[it1]-plcs[it2]; 
            auto pltime = grhtimes[it1]-grhtimes[it2];

			//auto gaplr = std::sqrt(invltres[it1]*invltres[it2]); 
			auto gaplr = std::sqrt(std::sqrt(invltres[it1])*std::sqrt(invltres[it2]));

			if( plz != 999 ){
            	auto psl = pltime/plz;
            	pszs1.push_back(psl);
				//auto gaplr = (slvars1[it1] + slvars1[it2]);
				plsrs.push_back( gaplr );
				plssum += psl*gaplr;
			}//<<>>if( plz != 999 )

            if( pl3z != 999 ){
            	auto p3sl = pltime/pl3z;
            	p3zs1.push_back(p3sl);
            	//auto gaplr = (slvars3[it1] + slvars3[it2]);
            	pl3rs.push_back( gaplr );
            	pl3sum += p3sl*gaplr;
            }//<<>>if( plz != 999 )

            if( plgz != 999 ){
            	auto pgsl = pltime/plgz;
            	pgzs1.push_back(pgsl);
            	//auto gaplr = (slvarsgeo[it1] + slvarsgeo[it2]);
            	plgrs.push_back( gaplr );
            	plgsum += pgsl*gaplr;
            }//<<>>if( plz != 999 )

        }//<<>>for( uInt it2(it1); it2 < grhtimes.size(); it2++ )
    }//<<>>for( uInt it1(0); it1 < grhtimes.size(); it1++ )

//--------------------------------------------------------------------
    //std::cout << "9, ";


    //find eigan vecgtor aligment
    //float eigensOld2d0(eigens[0]), eigensOld2d1(eigens[1]);
    //if( xsum1 < 0 ){ eigensOld2d0 *= -1; eigensOld2d1 *= -1; xsum1 = std::abs(xsum1);}
    //if( plssum < 0 ){ eigens[0] *= -1; eigens[1] *= -1; plssum = std::abs(plssum);}
    //if( pl3sum < 0 ){ eigens3[0] *= -1; eigens3[1] *= -1; eigens3[2] *= -1; pl3sum = std::abs(pl3sum);}
    //if( plgsum < 0 ){ geoeigens[0] *= -1; geoeigens[1] *= -1; plgsum = std::abs(plgsum);}


	// -----------   compute final outputs --------------------------------

	auto phiCorrFactor = 0.8;
	auto sxx = var( letas, 0., rhlogwtvec );
	auto syy = var( lphis, 0., rhlogwtvec, accum(rhlogwtvec)/phiCorrFactor );
	auto sxy = cvar( letas, 0., lphis, 0., rhlogwtvec, accum(rhlogwtvec)/std::sqrt(phiCorrFactor) );
	auto smaj = (sxx + syy + std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
    auto smin = (sxx + syy - std::sqrt(sq2(sxx - syy) + 4.*sq2(sxy)))/2.;
	auto sang = std::atan((sxx-syy+std::sqrt(sq2(syy-sxx)+4.*sq2(sxy)))/(2.*sxy));

    auto ebside = ( mz > 0 ) ? 1 : -1;
    auto taflip = ( ((std::abs(mz) < std::abs(vtxZ)) && (mz*vtxZ > 0) ) ? -1 : 1 )*ebside;

	//2d egians slope
    auto nXSum = xs1.size();
	auto totSloRes1 = accum(slvars1);
    auto slope1 = xsum1/totSloRes1;
    auto slope1err = std::sqrt(1/totSloRes1);
    auto varsl = var(xs1,slope1,slvars1,totSloRes1);
    auto chi1 = chisq(xs1,slope1,varsl);
	auto slchi2v1 = chisqv(xs1,slope1,slvars1,varsl);//?????????????????
    auto chi2pf1 = 1 - TMath::Prob(chi1, nWts);

	// 3d eigan slope
    auto totSloRes3 = accum(pl3rs);
    auto slope3 = pl3sum/totSloRes3;
    auto slope3err = std::sqrt(1/totSloRes3);

	//geo eigan slope
    auto totSloResGeo = accum(plgrs);
    auto slopeg = plgsum/totSloResGeo;
    auto slopegerr = std::sqrt(1/totSloResGeo);

	//pairs slope
    auto totSloResPrs = accum(plsrs);
    auto slopeprs = plssum/totSloResPrs;
    auto slopeprserr = std::sqrt(1/totSloResPrs);

    // eigan3d angle slope
    auto angle3d = getAngle(eigens3[1],eigens3[2]);
    auto hypo3d02 = hypo(eigens3[0],eigens3[2]);
    auto hypo3d01 = hypo(eigens3[0],eigens3[1]);
    auto hypo3d12 = hypo(eigens3[1],eigens3[2]);
    auto slope3d = std::atan(eigens3[2]/hypo3d01);
    //auto slope3d = 100*(eigens3[0]/hypo3d12)/SOL;

	// Fill Histograms
	for( uInt it(0); it < lzs.size(); it++ ){

        //plzs.push_back(xscor);
        //plcs.push_back(yscor);
        //pl3zs.push_back(x3cor);
        //pl3cs.push_back(y3cor);
        //plgzs.push_back(xgcor);
        //plgcs.push_back(ygcor);
        if( invltres[it] == 0 ) continue;
		auto sqrtwt = std::sqrt(invltres[it]);
		auto wichsum = slopeprs;
        //auto wichsum = pl3sum;
        //auto wichsum = xsum1;
        //auto wichsum = 1;
        auto xcor = ( wichsum < 0 ) ? -1*plzs[it] : plzs[it];
        auto ycor = ( wichsum < 0 ) ? -1*plcs[it] : plcs[it];
        auto fill = lts[it]*sqrtwt;
		hist2d[206]->Fill(xcor,ycor,fill);
		hist2d[207]->Fill(xcor,ycor,sqrtwt);
		hist2d[208]->Fill(lzs[it],lcs[it],fill);
		hist2d[209]->Fill(lzs[it],lcs[it],sqrtwt);
		hist2d[214]->Fill(ycor,lts[it],sqrtwt);
		hist2d[215]->Fill(xcor,lts[it],sqrtwt);

	}//<<>>for( uInt it(0); it < lzs.size(); it++ )

	// Fill results vector
    // std::cout << "10, ";
	// eigens 0 = vector x, 1 = vector y, 2 = vector mag
    eigens.push_back(slope1);//3  aligned slope
    eigens.push_back(chi2pf1);//4 aligned slope chi sqr prob
    eigens.push_back(slope3);//5 3d pairs slope
    eigens.push_back(0);//6 
    eigens.push_back(rotangle);//7 aligned rotation angle
    eigens.push_back(nXSum);//8 # of entries ( rechits )
    eigens.push_back(rot3angle);//9 3d rotation angle
    eigens.push_back(std::sqrt(varsl));//10 stdev aligned slope
    eigens.push_back(slope1err);//11 err aligned slope
    eigens.push_back(0);//12 
    eigens.push_back(slope3err);//13 errr unaligned slope
    eigens.push_back(slchi2v1);//14 chisqr like gof aligned slope
    eigens.push_back(minDr);//15 
    eigens.push_back(geoeigens[0]);//16 geoeigan x vec
    eigens.push_back(geoeigens[1]);//17 geoeigan y vec
    eigens.push_back(geoeigens[2]);//18 geoeigan mag vec
    eigens.push_back(ebside);//19 EB side
    eigens.push_back(taflip);//20 towards(+)/away(-) slope sign flip
    eigens.push_back(0);//21
    eigens.push_back(slopeg);//22 geoeigan slope
    eigens.push_back(slopegerr);//23 geoeigan slope error
    eigens.push_back(slopeprs);//24 pairs slope method
    eigens.push_back(slopeprserr);//25 pairs slope method error
    eigens.push_back(0);//26
    eigens.push_back(slope3d);//27 slope from 3d major eiganvector in "time" deminsion
    eigens.push_back(eigens3[3]);//28 "scaled magnitude" of 3d eiganvalue3
    eigens.push_back(angle3d);//29 rotation angle for eigan 3d
    eigens.push_back(eigens3[0]);//30 ev x	"z"
    eigens.push_back(eigens3[1]);//31 ev y	"c"
    eigens.push_back(eigens3[2]);//32 ev z  "time"
    eigens.push_back(eigens3[6]);//33 3dEV 6 = c vs c+t , ? 3d4Value
    eigens.push_back(eigens3[5]);//34 3d Time EV 5 = t vs t+zc oval ? 3dTvalue
	//std::cout << "Slope egin : " << slope1 << " " << chi2pf1 << " " << rotangle << " " << std::sqrt(varsl) << " " << slope1err << std::endl;
    eigens.push_back(smaj);//35
    eigens.push_back(smin);//36
    eigens.push_back(sang);//37
    eigens.push_back(eigens3[8]);//38
    eigens.push_back(eigens3[9]);//39
    eigens.push_back(eigens3[10]);//40

    //std::cout << " Done" << std::endl;;
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
    const std::string listdir("llpgana_list_files/");

	std::cout << "Producing Histograms for : " << outfilename << std::endl;
    std::ifstream infile(listdir+infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::cout << " - With : " << infilelist << " >> " << fInTree << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
		//std::cout << "--  for Fine #" << cnt << " moduls " << cnt%pct << " ";
		if( cnt%pct == 0 ){ 
        	auto tfilename = eosdir + indir + str;
        	std::cout << "--  adding file: " << tfilename << std::endl;
        	fInTree->Add(tfilename.c_str());
		}//<<>>if( cnt%4 == 0 ){
		//else std::cout << " do not add file" << std::endl;
		cnt++;
    }//<<>>while (std::getline(infile,str))

	Init(fInTree);
	initHists();

	SetupDetIDsEB(DetIDMap);
	SetupDetIDsEE(DetIDMap);

    std::cout << "Setting up For Main Loop." << std::endl;
	int loopCounter(10000);
    auto nEntries = fInTree->GetEntries();
	//nEntries = 10000;
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

	nMaps = 0;
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

		if( (*genPdgId)[it] == 22 && (*genLLP)[it] == 1 ) nLlpPho++;
        auto genID = (*genPdgId)[it];
        //if( abs((*genPdgId)[it]) > 2000000 ) genID = abs((*genPdgId)[it]) - 2000000 + 200;
        //else if( abs((*genPdgId)[it]) > 1000000 ) genID = abs((*genPdgId)[it]) - 1000000 + 100;
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


		// detrimine photon classification ( sig/susy/ect... )
		//---------------------------------------------------
		auto isCmb = not (*phoExcluded)[it];
        //auto isCmb = phoExcluded[it];
        auto phoGenIdx = (*genPhoIdx)[it];
        auto isSUSY = (phoGenIdx >= 0)?((*genLLP)[phoGenIdx] < 700 ):false;
		auto isOOT = (*phoIsOotPho)[it];
		auto isSigPho = (phoGenIdx >= 0)?((*genLLP)[phoGenIdx] == 1):false;

        if( DEBUG ) std::cout << " -- looping photons : filling gen" << std::endl;
    
        if(phoGenIdx >= 0){
    
            hist1d[100]->Fill( (*genPt)[phoGenIdx] );
            fillTH1( (*genPhoDr)[it], hist1d[130]); //->Fill( (*genPhoDr)[it] );
            fillTH1( abs((*genPdgId)[phoGenIdx]), hist1d[131]);
            //auto genParentExoID = 98;
            //if( abs((*genLLP)[phoGenIdx]) > 2000000 ) genParentExoID = abs((*genLLP)[phoGenIdx]) - 2000000 + 200;
            //else if( abs((*genLLP)[phoGenIdx]) > 1000000 ) genParentExoID = abs((*genLLP)[phoGenIdx]) - 1000000 + 100;
            //if( abs((*genPdgId)[phoGenIdx]) == 0 ) genParentExoID = 95;
            hist1d[101]->Fill( (*genLLP)[phoGenIdx] );
    
        }//<<>>if(phoGenIdx >= 0)


		// determine pog id class
		// -----------------------------------------------------
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

        if( DEBUG ) std::cout << " -- looping photons : filling ids" << std::endl;
        hist1d[107]->Fill( (*phoIsPixelSeed)[it] );
        hist1d[108]->Fill( (*phoHadOverEM)[it] );
        hist1d[109]->Fill( (*phoR9)[it] );
        hist1d[128]->Fill( (*phoR9)[it] );//
        //hist1d[110]->Fill( (*phoMipTotEnergy)[it] );
        hist1d[111]->Fill( (*phoEcalRHSumEtConeDR04)[it] );
        hist1d[112]->Fill( (*phoHcalTwrSumEtConeDR04)[it] );

		//find intial cut values 
		//----------------------------------------------------------------
		if( DEBUG ) std::cout << " -- pho rechits" << std::endl;
		auto nrh = ((*phoRhIds)[it]).size();
        auto phoClstrR9 = clstrR9( (*phoRhIds)[it] );

        auto isEB = (*phoIsEB)[it];
		auto isTightEB = std::abs((*phoEta)[it]) < 1.45;

		auto isTight = phoClass == 3;
        auto isLoose = phoClass > 1;
        auto isFake = phoClass == 1;

        auto goodRhCnt = nrh > 14;
        auto minPhoPt = (*phoPt)[it] > 30;

        //auto usePho = true;
		//auto usePho = not isSUSY && isFake;
        //auto usePho = not isSUSY && isLoose;
        //auto usePho = not isSUSY && isLoose && isEB && isTightEB;
        //auto usePho = isSUSY && isFake;
        //auto usePho = isSUSY && isLoose && isEB && isTightEB;
        //auto usePho = not isSUSY && isLoose && phoExcluded[it];
        //auto usePho = isLoose && isEB && isTightEB;
        
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isClR9r26;
        //auto usePho = isLoose && isEB && isTightEB && isClR9r68;

        auto usePho = isLoose && isEB && isTightEB && isCmb;
        //auto usePho = isTight && isEB && isTightEB && isCmb && isSigPho;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && isSigPho;
        //auto usePho = isLoose && isEB && isTightEB && isCmb && not isSigPho;


		if( DEBUG ) std::cout << " -- photons # " << it << " isSUSY " << isSUSY << " isOOT " << isOOT << " isCmb " << isCmb << std::endl;
		if( usePho && goodRhCnt && minPhoPt ) { //----------------------------------- 
	
	        hist1d[132]->Fill(1);
	        if( isSUSY ) hist1d[132]->Fill(2);
	        if( isOOT ) hist1d[132]->Fill(3);
	        if( isCmb ) hist1d[132]->Fill(4);
	
			if( DEBUG ) std::cout << " -- pho time" << std::endl;
            hist1d[104]->Fill( (*phoPt)[it] );
            hist1d[105]->Fill( (*phoEnergy)[it] );
			hist2d[205]->Fill( (*phoEnergy)[it], (*phoPt)[it] );

	        hist1d[119]->Fill( (*phoCMeanTime)[it] );//c mean
	        hist1d[120]->Fill( (*phoSeedTOFTime)[it] );//lead time 
	        hist1d[121]->Fill( (*phoCMeanTime)[it] - (*phoSeedTOFTime)[it] );//diff
			for( auto rhid : (*phoRhIds)[it] ){ hist1d[115]->Fill((*rhEnergy)[getRhIdx(rhid)]);}

            auto isPho = (phoGenIdx >= 0)?(((*genPdgId)[phoGenIdx] == 22)?1:0):-1;
            hist2d[204]->Fill(isPho,isSUSY);

			if( DEBUG ) std::cout << " Finding Eigans for Photon W/ " << nrh << " rechits." << std::endl;
			auto tofTimes = getLeadTofRhTime( (*phoRhIds)[it], vtxX, vtxY, vtxZ );
			auto phoEigens2D = getRhGrpEigen_sph( tofTimes,(*phoRhIds)[it] );
			auto isEiganGood = phoEigens2D[0] != -9;	
			//if( true ){
			if( isEiganGood ){

	            //auto eigtaf = phoEigens2D[20];
	            //auto eigtaf = 1;
				//auto rotAngle = getAngle( phoEigens2D[0], phoEigens2D[1] );
				//auto rot3DAnlge = phoEigens2D[9];// : (((phoEigens2D[9] + PI) < 2*PI) ? phoEigens2D[9] + PI : phoEigens2D[9] - PI);
				//auto geoRotAngle = getAngle( phoEigens2D[16], phoEigens2D[17] );	
/*
				//float dAngleMin = 0.1;
	            float dAngleMin = 0.2;
	            //float dAngleMin = 0.3;
				auto dGeoAngle1 = dltAngle(geoRotAngle,PI);
	            auto dRotAngle1 = dltAngle(rotAngle,PI);
				auto dRAngle1 = hypo(dGeoAngle1,dRotAngle1);
	            auto notAnglePeak1 = dRAngle1 > dAngleMin;
				auto angle2 = PI/2;
	            auto dGeoAngle2 = dltAngle(geoRotAngle,angle2);
	            auto dRotAngle2 = dltAngle(rotAngle,angle2);
	            auto dRAngle2 = hypo(dGeoAngle2,dRotAngle2);
	            auto notAnglePeak2 = dRAngle2 > dAngleMin;
	            auto angle3 = 3*PI/2;
	            auto dGeoAngle3 = dltAngle(geoRotAngle,angle3);
	            auto dRotAngle3 = dltAngle(rotAngle,angle3);
	            auto dRAngle3 = hypo(dGeoAngle3,dRotAngle3);
	            auto notAnglePeak3 = dRAngle3 > dAngleMin;
	            auto angle4 = 2*PI;
	            auto dGeoAngle4 = dltAngle(geoRotAngle,angle4);
	            auto dRotAngle4 = dltAngle(rotAngle,angle4);
	            auto dRAngle4 = hypo(dGeoAngle4,dRotAngle4);
	            auto notAnglePeak4 = dRAngle4 > dAngleMin;
	            auto angle5 = 0;
	            auto dGeoAngle5 = dltAngle(geoRotAngle,angle5);
	            auto dRotAngle5 = dltAngle(rotAngle,angle5);
	            auto dRAngle5 = hypo(dGeoAngle5,dRotAngle5);
	            auto notAnglePeak5 = dRAngle5 > dAngleMin;
	            //auto notAnglePeak = true;
	            //auto notAnglePeak = notAnglePeak1;
				auto notAnglePeak = notAnglePeak1 && notAnglePeak2 && notAnglePeak3 && notAnglePeak4 && notAnglePeak5;
		
				//if( not notAnglePeak ) std::cout << "dAngle cuts 1: " << dRAngle3 << " 2: " << angle3 << " 3: " << dGeoAngle3 << " 4: " << dRotAngle3 << std::endl;
*/
		
				//auto evalue2d = phoEigens2D[2];
	        	auto evaluegeo = phoEigens2D[18];
	            //auto geoslope = phoEigens2D[22];
				//auto slope = phoEigens2D[27];// 24 -> pairs 2d slope, 3 -> old slope, 5 -> pairs 3d slope
	            //auto slope = phoEigens2D[24]; //default;
	            //auto slope = phoEigens2D[5];

				auto smmcut = (*phoSMaj)[it] < 1.3 && (*phoSMin)[it] < 0.4;
	
	        	auto isClR9r46 = phoClstrR9 > 0.4 && phoClstrR9 < 0.6;
	        	auto isClR9r68 = phoClstrR9 > 0.6 && phoClstrR9 < 0.8;
	        	auto isClR9r26 = phoClstrR9 > 0.2 && phoClstrR9 < 0.6;
	        	auto isClR9lt0p7 = phoClstrR9 < 0.7;
	
				auto maxPhoGeoEV = evaluegeo < 0.85;
				auto minRHCnt = nrh > 25;
				//auto maxEvRh = evalue2d < ( 0.005*nrh + 0.8 );
				//auto min3dTval = phoEigens2D[34] > 0.02;
                //auto max3d4val = phoEigens2D[33] < 0.8;	
				//auto minSlope = std::abs(slope) > 0.05;

				//auto passCuts = maxPhoGeoEV;
				//auto passCuts = maxPhoGeoEV && notAnglePeak;
	            //auto passCuts = maxPhoGeoEV && minRHCnt;
				//auto passCuts = maxPhoGeoEV && maxEvRh;
	            //auto passCuts = maxPhoGeoEV && notAnglePeak && minRHCnt;
	            //auto passCuts = min3dTval;
	            //auto passCuts = max3d4val;
	            //auto passCuts = minSlope;

				auto passCuts = true;
				//auto passCuts = maxPhoGeoEV;
                //auto passCuts = smmcut;
                //auto passCuts = smmcut && maxPhoGeoEV;	

				if( passCuts ){

                    makeEBEEMaps( it );
	
	            	hist1d[102]->Fill( phoClstrR9 );//<<   
					hist1d[124]->Fill( phoEigens2D[35] );//smaj //<<
					hist1d[125]->Fill( phoEigens2D[36] );//smin //<<
					hist1d[126]->Fill( (*phoSigmaIEtaIEta)[it] );
	            	hist1d[146]->Fill( (*phoPt)[it] );
	            	hist1d[147]->Fill( (*phoEnergy)[it] );
					hist2d[217]->Fill( (*phoEnergy)[it], (*phoPt)[it] );
	            	hist1d[148]->Fill( (*phoCMeanTime)[it] );//c mean
                    hist1d[103]->Fill( nrh );//<<
                    hist1d[133]->Fill( evaluegeo );//<<
                    hist1d[152]->Fill( (*phoSMaj)[it] );
                    hist1d[153]->Fill( (*phoSMin)[it] );
                    hist1d[154]->Fill( (*phoSAlp)[it] );
                    hist1d[155]->Fill( (*phoCovEtaEta)[it] );
                    hist1d[156]->Fill( (*phoCovEtaPhi)[it] );
                    hist1d[157]->Fill( (*phoCovPhiPhi)[it] );

	            	hist2d[218]->Fill( nrh, phoClstrR9 );
                    hist2d[222]->Fill( nrh, phoEigens2D[35] );
                    hist2d[223]->Fill( nrh, phoEigens2D[36] );
                    hist2d[224]->Fill( nrh, (*phoSigmaIEtaIEta)[it] );
	                hist2d[219]->Fill( nrh, (*phoEnergy)[it] );
	                hist2d[220]->Fill( nrh, (*phoCMeanTime)[it] );
                    hist2d[229]->Fill( nrh, evaluegeo );
                    hist2d[228]->Fill( nrh, (*phoSMaj)[it] );	
                    hist2d[230]->Fill( nrh, (*phoSMin)[it] );
                    hist2d[231]->Fill( nrh, (*phoSAlp)[it] );
                    hist2d[232]->Fill( nrh, (*phoCovEtaEta)[it] );
                    hist2d[233]->Fill( nrh, (*phoCovEtaPhi)[it] );
                    hist2d[234]->Fill( nrh, (*phoCovPhiPhi)[it] );

                    hist2d[247]->Fill( phoClstrR9, evaluegeo );
                    hist2d[252]->Fill( phoEigens2D[35], evaluegeo );
                    hist2d[255]->Fill( phoEigens2D[36], evaluegeo );
                    hist2d[258]->Fill( (*phoSigmaIEtaIEta)[it], evaluegeo );
                    hist2d[235]->Fill( (*phoEnergy)[it], evaluegeo );
                    hist2d[236]->Fill( (*phoCMeanTime)[it], evaluegeo );
                    hist2d[237]->Fill( (*phoSMaj)[it], evaluegeo ); 
                    hist2d[238]->Fill( (*phoSMin)[it], evaluegeo );
                    hist2d[239]->Fill( (*phoSAlp)[it], evaluegeo );
                    hist2d[240]->Fill( (*phoCovEtaEta)[it], evaluegeo );
                    hist2d[241]->Fill( (*phoCovEtaPhi)[it], evaluegeo );
                    hist2d[242]->Fill( (*phoCovPhiPhi)[it], evaluegeo );

                    hist2d[225]->Fill( phoClstrR9, phoEigens2D[35] );
                    hist2d[226]->Fill( phoClstrR9, phoEigens2D[36] );
                    hist2d[227]->Fill( phoClstrR9, (*phoSigmaIEtaIEta)[it] );
                    hist2d[245]->Fill( phoClstrR9, (*phoEnergy)[it] );
                    hist2d[246]->Fill( phoClstrR9, (*phoCMeanTime)[it] );
                    hist2d[248]->Fill( phoClstrR9, (*phoSMaj)[it] ); 
                    hist2d[249]->Fill( phoClstrR9, (*phoSMin)[it] );
                    hist2d[250]->Fill( phoClstrR9, (*phoSAlp)[it] );
                    hist2d[251]->Fill( phoClstrR9, (*phoCovEtaEta)[it] );
                    hist2d[252]->Fill( phoClstrR9, (*phoCovEtaPhi)[it] );
                    hist2d[253]->Fill( phoClstrR9, (*phoCovPhiPhi)[it] );

                    hist2d[254]->Fill( phoEigens2D[35], phoEigens2D[36] );
                    hist2d[255]->Fill( phoEigens2D[35], (*phoSigmaIEtaIEta)[it] );
                    hist2d[256]->Fill( phoEigens2D[35], (*phoEnergy)[it] );
                    hist2d[257]->Fill( phoEigens2D[35], (*phoCMeanTime)[it] );
                    hist2d[259]->Fill( phoEigens2D[35], (*phoSMaj)[it] );
                    hist2d[260]->Fill( phoEigens2D[35], (*phoSMin)[it] );
                    hist2d[261]->Fill( phoEigens2D[35], (*phoSAlp)[it] );
                    hist2d[262]->Fill( phoEigens2D[35], (*phoCovEtaEta)[it] );
                    hist2d[263]->Fill( phoEigens2D[35], (*phoCovEtaPhi)[it] );
                    hist2d[264]->Fill( phoEigens2D[35], (*phoCovPhiPhi)[it] );

                    hist2d[265]->Fill( phoEigens2D[36], (*phoSigmaIEtaIEta)[it] );
                    hist2d[266]->Fill( phoEigens2D[36], (*phoEnergy)[it] );
                    hist2d[267]->Fill( phoEigens2D[36], (*phoCMeanTime)[it] );
                    hist2d[269]->Fill( phoEigens2D[36], (*phoSMaj)[it] );
                    hist2d[270]->Fill( phoEigens2D[36], (*phoSMin)[it] );
                    hist2d[271]->Fill( phoEigens2D[36], (*phoSAlp)[it] );
                    hist2d[272]->Fill( phoEigens2D[36], (*phoCovEtaEta)[it] );
                    hist2d[273]->Fill( phoEigens2D[36], (*phoCovEtaPhi)[it] );
                    hist2d[274]->Fill( phoEigens2D[36], (*phoCovPhiPhi)[it] );

                    hist2d[275]->Fill( (*phoSigmaIEtaIEta)[it], (*phoEnergy)[it] );
                    hist2d[276]->Fill( (*phoSigmaIEtaIEta)[it], (*phoCMeanTime)[it] );
                    hist2d[278]->Fill( (*phoSigmaIEtaIEta)[it], (*phoSMaj)[it] );
                    hist2d[279]->Fill( (*phoSigmaIEtaIEta)[it], (*phoSMin)[it] );
                    hist2d[280]->Fill( (*phoSigmaIEtaIEta)[it], (*phoSAlp)[it] );
                    hist2d[281]->Fill( (*phoSigmaIEtaIEta)[it], (*phoCovEtaEta)[it] );
                    hist2d[282]->Fill( (*phoSigmaIEtaIEta)[it], (*phoCovEtaPhi)[it] );
                    hist2d[283]->Fill( (*phoSigmaIEtaIEta)[it], (*phoCovPhiPhi)[it] );

                    hist2d[284]->Fill( (*phoEnergy)[it], (*phoCMeanTime)[it] );
                    hist2d[286]->Fill( (*phoEnergy)[it], (*phoSMaj)[it] );
                    hist2d[287]->Fill( (*phoEnergy)[it], (*phoSMin)[it] );
                    hist2d[288]->Fill( (*phoEnergy)[it], (*phoSAlp)[it] );
                    hist2d[289]->Fill( (*phoEnergy)[it], (*phoCovEtaEta)[it] );
                    hist2d[290]->Fill( (*phoEnergy)[it], (*phoCovEtaPhi)[it] );
                    hist2d[291]->Fill( (*phoEnergy)[it], (*phoCovPhiPhi)[it] );

                    hist2d[293]->Fill( (*phoCMeanTime)[it], (*phoSMaj)[it] );
                    hist2d[294]->Fill( (*phoCMeanTime)[it], (*phoSMin)[it] );
                    hist2d[295]->Fill( (*phoCMeanTime)[it], (*phoSAlp)[it] );
                    hist2d[296]->Fill( (*phoCMeanTime)[it], (*phoCovEtaEta)[it] );
                    hist2d[297]->Fill( (*phoCMeanTime)[it], (*phoCovEtaPhi)[it] );
                    hist2d[298]->Fill( (*phoCMeanTime)[it], (*phoCovPhiPhi)[it] );
	
					//hist1d[136]->Fill( phoEigens2D[0] );
		            //hist1d[137]->Fill( phoEigens2D[1] );
		            //hist1d[138]->Fill( phoEigens2D[16] );
		            //hist1d[139]->Fill( phoEigens2D[17] );
		            //hist1d[140]->Fill( phoEigens2D[15] );
		
		            //hist1d[141]->Fill( phoEigens2D[30] );//ev1 : x  "z"
		            //hist1d[142]->Fill( phoEigens2D[31] );//ev1 : y  "c"
		            //hist1d[143]->Fill( phoEigens2D[32] );//ev1 : z  "time"
		            //hist1d[144]->Fill( phoEigens2D[27] );//c vs c+t , ? 3d4Value//33
		            //hist1d[145]->Fill( phoEigens2D[34] );//t vs t+zc oval ? 3dTvalue
                    //hist2d[216]->Fill( eigtaf*slope, phoEigens2D[34] );


                    //hist1d[149]->Fill( phoEigens2D[38] );//3dev1
                    //hist1d[150]->Fill( phoEigens2D[39] );//3dev2
                    //hist1d[151]->Fill( phoEigens2D[40] );//3dev3

					////hist2d[222]->Fill( phoClstrR9, phoEigens2D[15] );
		
		            //if( DEBUG ) std::cout << " -- start eigan hists fill p0:" << std::endl;
		            //hist1d[127]->Fill( rotAngle );//eliptical angle  Eigan Rotation Angle
		            //hist1d[106]->Fill( rot3DAnlge );//eliptical angle  Eta Phi Angle 2D : algined
		            //hist2d[203]->Fill( rotAngle, rot3DAnlge );
		            //hist2d[210]->Fill( phoEigens2D[0], phoEigens2D[1] );
		            //hist2d[211]->Fill( phoEigens2D[30], phoEigens2D[31] );
		            //hist2d[212]->Fill( phoEigens2D[30], phoEigens2D[32] );
		            //hist2d[213]->Fill( phoEigens2D[31], phoEigens2D[32] );
		            //hist1d[113]->Fill( evalue2d ); 
		            //if( DEBUG ) std::cout << " -- start eigan hists fill p0a:" << std::endl;
					//hist1d[114]->Fill(slope); // rot way slope aligned
		            //hist1d[117]->Fill(eigtaf*slope); // TA slope
					//hist2d[200]->Fill((*phoEta)[it],eigtaf*slope); 
		            //hist2d[202]->Fill((*phoEta)[it],phoEigens2D[27]);//3d vector slope
					//hist2d[201]->Fill((*phoEta)[it],evalue2d);
		            //if( DEBUG ) std::cout << " -- start eigan hists fill p0b:" << std::endl;
					//hist1d[]->Fill(phoEigens2D[10]);
					//hist1d[]->Fill(phoEigens2D[25]);
					//hist1d[125]->Fill(phoEigens2D[12]);
					//hist1d[126]->Fill(phoEigens2D[25]);
		            //if( DEBUG ) std::cout << " -- start eigan hists fill p0c:" << std::endl;
					//hist2d[]->Fill( phoEigens2D[10], phoEigens2D[25] );
		            //hist2d[]->Fill( eigtaf*phoEigens2D[24], phoEigens2D[10] );
		            //if( DEBUG ) std::cout << " -- start eigan hists fill p0d:" << std::endl;
		            ////hist2d[222]->Fill( phoEigens2D[25], phoClstrR9 );
					//hist2d[]->Fill( (*phoEta)[it], rotAngle ); 
		            //hist2d[]->Fill( (*phoPhi)[it], rotAngle );
		            //hist2d[]->Fill( rotAngle, phoEigens2D[2] );
		            //hist2d[]->Fill( phoEigens2D[9], eigtaf*slope );
		            //if( DEBUG ) std::cout << " -- eigan hists fill p1:" << std::endl;
		
		            //hist2d[232]->Fill( phoClstrR9, rotAngle );
		            //hist2d[233]->Fill( evalue2d, eigtaf*slope );
		            //hist2d[235]->Fill( rotAngle, evalue2d );
		            //hist2d[234]->Fill( rotAngle, eigtaf*slope );
		            //hist2d[230]->Fill( phoClstrR9, evalue2d );
		            //hist2d[231]->Fill( phoClstrR9, eigtaf*slope );
		
		            //hist1d[134]->Fill( geoslope );
		            //hist1d[135]->Fill( geoRotAngle );
		            //if( DEBUG ) std::cout << " -- eigan hists fill p2:" << std::endl;
		
		            //hist2d[236]->Fill( evaluegeo, rotAngle );
		            //hist2d[237]->Fill( eigtaf*geoslope, rotAngle );
		            //hist2d[238]->Fill( eigtaf*geoslope, evalue2d );
		            //hist2d[239]->Fill( geoRotAngle, rotAngle );
		            //hist2d[240]->Fill( eigtaf*geoslope, eigtaf*slope );
		            //hist2d[241]->Fill( geoRotAngle, eigtaf*slope );
		            //hist2d[242]->Fill( geoRotAngle, evalue2d );
		            //hist2d[243]->Fill( evaluegeo, evalue2d );
		            //hist2d[244]->Fill( evaluegeo, eigtaf*slope );
		
		            //if( DEBUG ) std::cout << " -- eigan hists fill p3:" << std::endl;
		            //hist2d[245]->Fill( geoRotAngle, evaluegeo );
		            //hist2d[246]->Fill( geoRotAngle, eigtaf*geoslope );
		            //hist2d[248]->Fill( phoClstrR9, eigtaf*geoslope );
		            //hist2d[249]->Fill( phoClstrR9, geoRotAngle );
		            //hist2d[250]->Fill( evaluegeo, eigtaf*geoslope );
		
					//hist2d[224]->Fill( nrh, rotAngle );
					//hist2d[226]->Fill( nrh, evalue2d );
		            //hist2d[227]->Fill( nrh, eigtaf*slope );
		            //hist2d[228]->Fill( nrh, eigtaf*geoslope );
		            //hist2d[221]->Fill( nrh, geoRotAngle );

					//hist2d[251]->Fill( phoEigens2D[35], evalue2d );
                    //hist2d[253]->Fill( phoEigens2D[35], phoEigens2D[27] );//33->27
                    //hist2d[254]->Fill( phoEigens2D[36], evalue2d );
                    //hist2d[256]->Fill( phoEigens2D[36], phoEigens2D[27] );//33->27
                    //hist2d[257]->Fill( (*phoSigmaIEtaIEta)[it], evalue2d );
                    //hist2d[259]->Fill( (*phoSigmaIEtaIEta)[it], phoEigens2D[27] );//33->27

                    //hist2d[228]->Fill( rotAngle, phoEigens2D[35] );
                    //hist2d[]->Fill( rotAngle, phoEigens2D[36] );
                    //hist2d[260]->Fill( rotAngle, (*phoSigmaIEtaIEta)[it] );

                    //hist2d[261]->Fill( geoRotAngle, phoEigens2D[35] );
                    //hist2d[262]->Fill( geoRotAngle, phoEigens2D[36] );
                    //hist2d[263]->Fill( geoRotAngle, (*phoSigmaIEtaIEta)[it] );


				}//<<>>if( doCuts ) //--------------------------------------

			}//<<>if( isEiganGood ){

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

//------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------

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
   //b_phoSc2dEv->GetEntry(entry);   //!
   b_phoPt->GetEntry(entry);   //!
   b_phoEnergy->GetEntry(entry);   //!
   b_phoPhi->GetEntry(entry);   //!
   b_phoEta->GetEntry(entry);   //!
   b_phoPx->GetEntry(entry);   //!
   b_phoPy->GetEntry(entry);   //!
   b_phoPz->GetEntry(entry);   //!
   b_phoRhIds->GetEntry(entry);   //!
   //b_phoIsPFPhoton->GetEntry(entry);   //!
   //b_phoIsStdPhoton->GetEntry(entry);   //!
   b_phoHasConTracks->GetEntry(entry);   //!
   b_phoIsPixelSeed->GetEntry(entry);   //!
   //b_phoIsPhoton->GetEntry(entry);   //!
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

   b_phoSMaj->GetEntry(entry);   //!
   b_phoSMin->GetEntry(entry);   //!
   b_phoSAlp->GetEntry(entry);   //!
   b_phoCovEtaEta->GetEntry(entry);   //!
   b_phoCovEtaPhi->GetEntry(entry);   //!
   b_phoCovPhiPhi->GetEntry(entry);   //!

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

/* 1D jets 0-99

    normTH1D(hist1d[3]);
    normTH1D(hist1d[4]);
    normTH1D(hist1d[34]);
    normTH1D(hist1d[35]);
    normTH1D(hist1d[38]);
    normTH1D(hist1d[39]);
    normTH1D(hist1d[36]);
    normTH1D(hist1d[37]);
    normTH1D(hist1d[32]);
    normTH1D(hist1d[33]);
*/

/* 1D cljets 250 - 299

    normTH1D(hist1d[256]);
*/

    thresDivTH2D( hist2d[206], hist2d[207], 0 );
    thresDivTH2D( hist2d[208], hist2d[209], 0 );

    profileTH2D( hist2d[215], hist1d[116], hist1d[118] );


/* jets 0 - 199

    normTH2D(hist2d[100]);
    normTH2D(hist2d[121]);
    normTH2D(hist2d[126]);
    normTH2D(hist2d[127]);
*/

    //1normTH2D(hist2d[202]);
    //normTH2D(hist2d[201]);
    //normTH2D(hist2d[200]);
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
    int jtdiv(300);
    float jtran(15);
    int jdtdiv(160);
    float jdtran(4);
    int rhcnt(160);

	// jet id stuff
    //auto stddiv = 120;
    //auto stdtran = 3;

	// eigan cluster 2d maps
    //auto cldiv = 1200;
    //auto cltrn = 30;
    auto cldiv = 36;
    auto cltrn = 19.8;

	// eigan cluster 3d maps
    //auto cl3ddiv = 160;
    //auto cl3dtrn = 4;
    //auto cl3ddiv1 = 160;
    //auto cl3dtrn1 = 4;

	// photon time phi - eta map
    auto clsphdiv = 72/2;
    auto clsphtrn = 19.8/2;

	// time profile maps
    auto cwdiv = 100;
    auto cwtrn = 5/2;

	// time cl slope plots
    auto slmax = 2.0;
    auto slmin = -2.0;
    auto sldiv = 200;

    auto sl3max = 2.0;
    auto sl3min = -2.0;
    auto sl3div = 200;

    //auto chimax = 1.01;
    //auto chimin = 0.91;
    //auto chidiv = 2880;


	//------------------------------------------------------------------------------------------
    //------ 1D Hists --------------------------------------------------------------------------

/*
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
*/

	//----- photons 100 - 249

    hist1d[100] = new TH1D("genPhoPt", "genPhoPt;Pt [GeV]",500,0,1000);
    hist1d[101] = new TH1D("genPhoOriginCode", "genPhoOriginCode;Code",1000,0,1000);

    hist1d[102] = new TH1D("phoClstrR9", "phoClstrR9", 100, 0, 1);

    hist1d[103] = new TH1D("phoNumRecHits", "phoNumRecHits", 80, 0, 80);
    hist1d[104] = new TH1D("phoPtPreCut", "phoPt PreCut", 1000, 0, 1000);
    hist1d[105] = new TH1D("phoEnergyPreCut", "phoEnergy PreCut", 2000, 0, 2000);

    //hist1d[106] = new TH1D("pho3DRotAngle", "pho3DRotAngle", 140, -0.5, 6.5);

    hist1d[107] = new TH1D("phoIsPixelSeed", "phoIsPixelSeed", 3, 0, 2);
    hist1d[108] = new TH1D("phoHadOverEM", "phoHadOverEM", 250, 0, 10);
    hist1d[109] = new TH1D("phoR9", "phoR9", 100, 0, 1);
    //hist1d[110] = new TH1D("phoMipTotEnergy", "phoMipTotEnergy", 250, 0, 250 );
    hist1d[111] = new TH1D("phoEcalRHSumEtConeDR04", "phoEcalRHSumEtConeDR04", 750, 0, 750);
    hist1d[112] = new TH1D("phoHcalTwrSumEtConeDR04", "phoHcalTwrSumEtConeDR04", 750, 0, 750);

    //hist1d[113] = new TH1D("phoEginValueSph", "phoEginValueSph", 150, 0.4, 1.1);
    //hist1d[114] = new TH1D("phoEginSlopeAP", "phoEginSlopeSph 1Way", sldiv, slmin, slmax);
    hist1d[115] = new TH1D("clRhEnergy", "clRhEnergy", 3000, 0, 1000);

    hist1d[116] = new TH1D("pho_etprofile", "Photon rEta Time Profile Sph;MajAxis [cm]", clsphdiv, -1*clsphtrn, clsphtrn);

    //hist1d[117] = new TH1D("phoEginSlopeSphTA", "phoEginSlopeSph TA", sldiv, slmin, slmax);
	hist1d[118] = new TH1D("pho_proFitEtavChi", "Profile Fit Eta v Chi2Prob Sph", cwdiv, -1*cwtrn, cwtrn );

    hist1d[119] = new TH1D("phoClTimePreCut", "phoTime PreCut", jtdiv, -1*jtran, jtran);
    hist1d[120] = new TH1D("phoSeedRhTime", "phoLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[121] = new TH1D("phoSeedTimeDiff", "phoLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[122] = new TH1D("phoClass", "phoClass", 4, 0, 4);

    hist1d[123] = new TH1D("phoClRhTime","phoClRhTime",1000, -50, 50);
    hist1d[124] = new TH1D("phoClSMaj","phoClSMaj",50,0,5);
    hist1d[125] = new TH1D("phoClSMin","phoClSMin",50,0,0.5);
    hist1d[126] = new TH1D("phoSigIEtaIEta","phoSigIEtaIEta",300,0,0.03);

    //hist1d[127] = new TH1D("phoRotAngle","phoRotAngle",140,-0.5,6.5);

    hist1d[128] = new TH1D("phoR9_zoom", "phoR9_zoom", 200, 0.8, 1);
    //hist1d[129] = new TH1D("phoGeoEginValueSph", "phoGeoEginValueSph", 150, 0.4, 1.1);
    hist1d[130] = new TH1D("phoGenDr", "phoGenDr", 800, 0, 2.0 );
    hist1d[131] = new TH1D("phoGenPdgId", "phoGenPdgId", 50, 0, 50 );
    hist1d[132] = new TH1D("phoType", "phoType(1-pho 2-susy 3-oot 4-cmb)", 6,0,5);

    hist1d[133] = new TH1D("phoGeoValue", "phoGeoValue", 150, 0.4, 1.1);
    //hist1d[134] = new TH1D("phoGeoSlope", "phoGeoSlope", sldiv, slmin, slmax);
    //hist1d[135] = new TH1D("phoGeoAngle","phoGeoAngle",70,0,7);

    //hist1d[136] = new TH1D("pho2DEiganV0","pho2DEiganV0",200,-1,1);
    //hist1d[137] = new TH1D("pho2DEiganV1","pho2DEiganV1",200,-1,1);
    //hist1d[138] = new TH1D("phoG2DEiganV0","phoG2DEiganV0",200,-1,1);
    //hist1d[139] = new TH1D("phoG2DEiganV1","phoG2DEiganV1",200,-1,1);

    hist1d[140] = new TH1D("phoMinDr","phoMinDr",500,0,5);

    //hist1d[141] = new TH1D("pho3DEiganV0","pho3DEiganV0",200,-1,1);
    //hist1d[142] = new TH1D("pho3DEiganV1","pho3DEiganV1",200,-1,1);
    //hist1d[143] = new TH1D("pho3DEiganV2","pho3DEiganV2",200,-1,1);
    //hist1d[144] = new TH1D("pho3d4Value","pho3d4Value",100,0,1);
    //hist1d[145] = new TH1D("pho3dTValue","pho3dTValue",400,0,0.4);

    hist1d[146] = new TH1D("phoPtPostCut", "phoPt PostCut", 1000, 0, 1000);
    hist1d[147] = new TH1D("phoEnergyPostCut", "phoEnergy PostCut", 2000, 0, 2000);
    hist1d[148] = new TH1D("phoClTimePostCut", "phoTime PostCut", jtdiv, -1*jtran, jtran);

    //hist1d[149] = new TH1D("pho3dEV1","pho3dEV1",20,0,20);
    //hist1d[150] = new TH1D("pho3dEV2","pho3dEV2",20,0,20);
    //hist1d[151] = new TH1D("pho3dEV3","pho3dEV3",20,0,20);

    hist1d[152] = new TH1D("phoSMaj", "phoSMaj", 300,0,3);
    hist1d[153] = new TH1D("phoSMin", "phoSMin", 60,0,.6);
    hist1d[154] = new TH1D("phoSAlp", "phoSAlp", 350,-1.75,1.75);
    hist1d[155] = new TH1D("phoCovEtaEta", "phoCovEtaEtaSqrt", 100,0,0.001);
    hist1d[156] = new TH1D("phoCovEtaPhi", "phoCovEtaPhiSqrt", 800,-0.0004,0.0004);
    hist1d[157] = new TH1D("phoCovPhiPhi", "phoCovPhiPhiSqrt", 100,0,0.001);

	//------  genparticles 200 - 249

    hist1d[200] = new TH1D("genPdgId", "genPdgId;PdgId",220,0,220);
    hist1d[201] = new TH1D("nGenLLPPho", "nGenLLPPho", 5, 0, 5 );

	//------  cluster jet 250 - 299
/*
    hist1d[250] = new TH1D("cljTime", "cljTime", jtdiv, -1*jtran, jtran);
    hist1d[251] = new TH1D("cljLeadRhTime", "cljLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[252] = new TH1D("cljClstLeadTimeDiff", "cljClstLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[253] = new TH1D("cljDrTime", "cljDrTime", jtdiv, -1*jtran, jtran);
    hist1d[254] = new TH1D("cljDiffPt", "cljDiffPt", 1000, 0, 10);
    hist1d[255] = new TH1D("cljdPhi", "cljdPhi", 70, -3.5, 3.5);
    hist1d[256] = new TH1D("cljdtmu", "cljdtmu", jdtdiv, -1*jdtran, jdtran);
*/
    //------ ecal rechits 300 - 349
    hist1d[300] = new TH1D("rhRadius", "rhRadius", 300, 128, 131);

    hist1d[301] = new TH1D("ebRhTime", "ebRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[302] = new TH1D("ebRhEnergy", "ebRhEnergy", 1000, 0, 1000);
    hist1d[303] = new TH1D("eeRhTime", "eeRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[304] = new TH1D("eeRhEnergy", "eeRhEnergy", 1000, 0, 1000);
    hist1d[305] = new TH1D("ebRhkOOT", "ebRhkOOT", 3, 0, 1);
    hist1d[306] = new TH1D("eeRhkOOT", "eeRhkOOT", 3, 0, 1);

    //------  electrons 350 - 400

    //hist1d[350] = new TH1D("eleClTime", "eleClTime", jtdiv, -1*jtran, jtran);
    //hist1d[351] = new TH1D("eleSeedRhTime", "eleLeadRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[352] = new TH1D("eleClSeedTimeDiff", "eleClLeadTimeDiff", jtdiv, -1*jtran, jtran);

	//-------- event vars 400 - 450
	
	//hist1d[400] = new TH1D("nphoexclusions", " Number Pho (0) Excl (1) OOT (2) Excl (3) ", 4, 0, 4);

    for( int it = 0; it < n1dHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

	//------------------------------------------------------------------------------------------
    //------ 2D Hists --------------------------------------------------------------------------

	//------- jets ( time ) 0-49 ------------------------------
/*
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
*/

	//--- Photons 200 - 349 -------------------------------------------
    //hist2d[200] = new TH2D( "phoEta_TASpSlope", "Photon Eta Vs TA Slope;Eta;Slope [ns/cm]", 60, -1.5, 1.5, sldiv, slmin, slmax);
    //hist2d[201] = new TH2D( "phoEta_2dtValue", "Photon Eta Vs EiganValue2D;Eta;EiganValue",60, -1.5, 1.5, 25, 0.5, 1.0);
    //hist2d[202] = new TH2D( "phoEta_3DTAngle", "Photon Eta Vs 3D TAngle;Eta;TAngle", 60, -1.5, 1.5, 640, -3.2, 3.2);
	//hist2d[203] = new TH2D( "phoSAngle_pho3Dangle","phoRotAngle Vs pho3DAngle",140,-0.5,6.5,140,-0.5,6.5);

    hist2d[204] = new TH2D("phoID_SUSY", "phoIDvSUSY (-1 not matched);isPho;isSUSY", 3, -1, 2, 2, 0, 2);
    hist2d[205] = new TH2D("phoPreE_PrePt", "phoPreEvPrePt;Energy [GeV];Pt [GeV]", 2000, 0, 2000, 1000, 0, 1000);

    hist2d[206] = new TH2D("pho_tmap_rot", "Photon Time Map Rotated;MajAxis [cm];MinAxis [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[207] = new TH2D("pho_occmap_rot", "Photon Occ Map Rotated;MajAxis [cm];MinAxis [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[208] = new TH2D("pho_tmap", "Photon Time Map;Z [cm];C [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);
    hist2d[209] = new TH2D("pho_occmap", "Photon Occ Map;Z [cm];C [cm]", cldiv, -1*cltrn, cltrn, cldiv, -1*cltrn, cltrn);

    //hist2d[210] = new TH2D("phoSphEgnxy", "phoSphEgn z,c;z;c", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[211] = new TH2D("pho3DEgnxy", "pho3DEgn x,y;x;y", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[212] = new TH2D("pho3DEgnxz", "pho3DEgn x,z;x;z", 220, -1.1, 1.1, 220, -1.1, 1.1);
    //hist2d[213] = new TH2D("pho3DEgnyz", "pho3DEgn y,z;y;z", 220, -1.1, 1.1, 220, -1.1, 1.1);

    hist2d[214] = new TH2D("pho_ptwtmap", "Photon Phi(x) Time(y) WtMap Sph;MinAxis [cm];Time [ns]", clsphdiv, -1*clsphtrn, clsphtrn, cwdiv, -1*cwtrn, cwtrn );
    hist2d[215] = new TH2D("pho_etwtmap", "Photon Eta(x) Time(y) WtMap Sph;MajAxis [cm];Time [ns]", clsphdiv, -1*clsphtrn, clsphtrn, cwdiv, -1*cwtrn, cwtrn );

    //hist2d[216] = new TH2D("pho2dSlope_3dTValue","2D Slope v 3dTValue;Slope;3dTValue",200,-2,2,400,0,0.4);
    hist2d[217] = new TH2D("phoPostE_PostPt", "phoPostEvPostPt;Energy [GeV];Pt [GeV]", 2000, 0, 2000, 1000, 0, 1000);

    hist2d[218] = new TH2D("phoNRH_ClR9","Photon nClRecHits v clR9;nRecHits;ClusterR9",80,0,80,100,0,1);
    hist2d[222] = new TH2D("phoNRH_phoClSMaj","phoNRH_phoClSMaj;#rh;SMaj",80,0,80,200,0,2);
    hist2d[223] = new TH2D("phoNRH_phoClSMin","phoNRH_phoClSMin;#rh;SMin",80,0,80,80,0,0.8);
    hist2d[224] = new TH2D("phoNRH_phoClSigIEtaIEta","phoNRH_phoClSigIEtaIEta;#rh;sieie",80,0,80,250,0,0.025);
    hist2d[219] = new TH2D("phoNRH_phoE","Photon nClRecHits v phoE;nRecHits;Energy",80,0,80,2000,0,2000);
    hist2d[220] = new TH2D("phoNRH_phoClTime","Photon nClRecHits v ClTime;nRecHits;ClTime [ns]",80,0,80,jtdiv, -1*jtran, jtran);
    hist2d[229] = new TH2D("phoNRH_geoValue","Photon nClRecHits v geoValue;nRecHits;geoValue",80,0,80,50,0.5,1);
    hist2d[228] = new TH2D("phoNRH_phoSMaj","Photon nClRecHits v phoSMaj;nRecHits;phoSMaj",80,0,80,250,0,2.5);
    hist2d[230] = new TH2D("phoNRH_phoSMin","Photon nClRecHits v phoSMin;nRecHits;phoSMin",80,0,80,750,0,0.75);
    hist2d[231] = new TH2D("phoNRH_phoSAlp","Photon nClRecHits v phoSAlp;nRecHits;phoSAlp",80,0,80,350,-1.75,1.75);
    hist2d[232] = new TH2D("phoNRH_phoCovEtaEta","Photon nClRecHits v phoCovEtaEta;nRecHits;phoCovEtaEta",80,0,80,100,0,0.001);
    hist2d[233] = new TH2D("phoNRH_phoCovEtaPhi","Photon nClRecHits v phoCovEtaPhi;nRecHits;phoCovEtaPhi",80,0,80,100,-0.0005,0.0005);
    hist2d[234] = new TH2D("phoNRH_phoCovPhiPhi","Photon nClRecHits v phoCovPhiPhi;nRecHits;phoCovPhiPhi",80,0,80,100,0,0.001);

    hist2d[247] = new TH2D("phoClR9_geoValue","Photon ClR9 v geoValue;ClR9;geoValue",100,0,1,50,0.5,1);
    hist2d[252] = new TH2D("phoClSMaj_geoValue","phoClSMaj v geoValue;ClSMaj;geoValue",200,0,2,50,0.5,1);
    hist2d[255] = new TH2D("phoClSMin_geoValue","phoClSMin v geoValue;ClSMin;geoValue",80,0,0.8,50,0.5,1);
    hist2d[258] = new TH2D("phoClSigIEtaIEta_geoValue","phoClSigIEtaIEta v geoValue;ClSigIEtaIEta;geoValue",250,0,0.025,50,0.5,1);
    hist2d[235] = new TH2D("phoE_geoValue","Photon phoE v geoValue;Energy;geoValue",2000,0,2000,50,0.5,1);
    hist2d[236] = new TH2D("phoClTime_geoValue","Photon ClTime v geoValue;ClTime [ns];geoValue",jtdiv,-1*jtran,jtran,50,0.5,1);
    hist2d[237] = new TH2D("phoSMaj_geoValue","Photon SMaj v geoValue;phoSMaj;geoValue",250,0,2.5,50,0.5,1);
    hist2d[238] = new TH2D("phoSMin_geoValue","Photon SMin v geoValue;phoSMin;geoValue",750,0,0.75,50,0.5,1);
    hist2d[239] = new TH2D("phoSAlp_geoValue","Photon SAlp v geoValue;phoSAlp;geoValue",350,-1.75,1.75,50,0.5,1);
    hist2d[240] = new TH2D("phoCovEtaEta_geoValue","Photon CovEtaEta v geoValue;phoCovEtaEta;geoValue",100,0,0.001,50,0.5,1);
    hist2d[241] = new TH2D("phoCovEtaPhi_geoValue","Photon CovEtaPhi v geoValue;phoCovEtaPhi;geoValue",100,-0.0005,0.0005,50,0.5,1);
    hist2d[242] = new TH2D("phoCovPhiPhi_geoValue","Photon CovPhiPhi v geoValue;phoCovPhiPhi;geoValue",100,0,0.001,50,0.5,1);

    hist2d[225] = new TH2D("phoClR9_phoClSMaj","phoClR9_phoClSMaj;clr9;SMaj",100,0,1,200,0,2);
    hist2d[226] = new TH2D("phoClR9_phoClSMin","phoClR9_phoClSMin;clr9;SMin",100,0,1,80,0,0.8);
    hist2d[227] = new TH2D("phoClR9_phoClSigIEtaIEta","phoClR9_phoClSigIEtaIEta;clr9;sieie",100,0,1,250,0,0.025);
    hist2d[245] = new TH2D("phoClR9_phoE","Photon phoClR9 v phoE;clr9;Energy",100,0,1,2000,0,2000);
    hist2d[246] = new TH2D("phoClR9_phoClTime","Photon phoClR9 v ClTime;clr9;ClTime [ns]",100,0,1,jtdiv, -1*jtran, jtran);
    hist2d[248] = new TH2D("phoClR9_phoSMaj","Photon phoClR9 v phoSMaj;clr9;phoSMaj",100,0,1,250,0,2.5);
    hist2d[249] = new TH2D("phoClR9_phoSMin","Photon phoClR9 v phoSMin;clr9;phoSMin",100,0,1,750,0,0.75);
    hist2d[250] = new TH2D("phoClR9_phoSAlp","Photon phoClR9 v phoSAlp;clr9;phoSAlp",100,0,1,350,-1.75,1.75);
    hist2d[251] = new TH2D("phoClR9_phoCovEtaEta","Photon phoClR9 v phoCovEtaEta;clr9;phoCovEtaEta",100,0,1,100,0,0.001);
    hist2d[252] = new TH2D("phoClR9_phoCovEtaPhi","Photon phoClR9 v phoCovEtaPhi;clr9;phoCovEtaPhi",100,0,1,100,-0.0005,0.0005);
    hist2d[253] = new TH2D("phoClR9_phoCovPhiPhi","Photon phoClR9 v phoCovPhiPhi;clr9;phoCovPhiPhi",100,0,1,100,0,0.001);

    hist2d[254] = new TH2D("phoClSMaj_phoClSMin","phoClSMaj v phoClSMin;SMaj;SMin",200,0,2,80,0,0.8);
    hist2d[255] = new TH2D("phoClSMaj_phoClSigIEtaIEta","phoClSMaj v phoClSigIEtaIEta;SMaj;sieie",200,0,2,250,0,0.025);
    hist2d[256] = new TH2D("phoClSMaj_phoE","Photon phoClSMaj v phoE;SMaj;Energy",200,0,2,2000,0,2000);
    hist2d[257] = new TH2D("phoClSMaj_phoClTime","Photon phoClSMaj v ClTime;SMaj;ClTime [ns]",200,0,2,jtdiv, -1*jtran, jtran);
    hist2d[259] = new TH2D("phoClSMaj_phoSMaj","Photon phoClSMaj v phoSMaj;SMaj;phoSMaj",200,0,2,250,0,2.5);
    hist2d[260] = new TH2D("phoClSMaj_phoSMin","Photon phoClSMaj v phoSMin;SMaj;phoSMin",200,0,2,750,0,0.75);
    hist2d[261] = new TH2D("phoClSMaj_phoSAlp","Photon phoClSMaj v phoSAlp;SMaj;phoSAlp",200,0,2,350,-1.75,1.75);
    hist2d[262] = new TH2D("phoClSMaj_phoCovEtaEta","Photon phoClSMaj v phoCovEtaEta;SMaj;phoCovEtaEta",200,0,2,100,0,0.001);
    hist2d[263] = new TH2D("phoClSMaj_phoCovEtaPhi","Photon phoClSMaj v phoCovEtaPhi;SMaj;phoCovEtaPhi",200,0,2,100,-0.0005,0.0005);
    hist2d[264] = new TH2D("phoClSMaj_phoCovPhiPhi","Photon phoClSMaj v phoCovPhiPhi;SMaj;phoCovPhiPhi",200,0,2,100,0,0.001);

    hist2d[265] = new TH2D("phoClSMin_phoClSigIEtaIEta","phoClSMin v phoClSigIEtaIEta;SMin;sieie",80,0,0.8,250,0,0.025);
    hist2d[266] = new TH2D("phoClSMin_phoE","Photon phoClSMin v phoE;SMin;Energy",80,0,0.8,2000,0,2000);
    hist2d[267] = new TH2D("phoClSMin_phoClTime","Photon phoClSMin v ClTime;SMin;ClTime [ns]",80,0,0.8,jtdiv, -1*jtran, jtran);
    hist2d[269] = new TH2D("phoClSMin_phoSMaj","Photon phoClSMin v phoSMaj;SMin;phoSMaj",80,0,0.8,250,0,2.5);
    hist2d[270] = new TH2D("phoClSMin_phoSMin","Photon phoClSMin v phoSMin;SMin;phoSMin",80,0,0.8,750,0,0.75);
    hist2d[271] = new TH2D("phoClSMin_phoSAlp","Photon phoClSMin v phoSAlp;SMin;phoSAlp",80,0,0.8,350,-1.75,1.75);
    hist2d[272] = new TH2D("phoClSMin_phoCovEtaEta","Photon phoClSMin v phoCovEtaEta;SMin;phoCovEtaEta",80,0,0.8,100,0,0.001);
    hist2d[273] = new TH2D("phoClSMin_phoCovEtaPhi","Photon phoClSMin v phoCovEtaPhi;SMin;phoCovEtaPhi",80,0,0.8,100,-0.0005,0.0005);
    hist2d[274] = new TH2D("phoClSMin_phoCovPhiPhi","Photon phoClSMin v phoCovPhiPhi;SMin;phoCovPhiPhi",80,0,0.8,100,0,0.001);

    hist2d[275] = new TH2D("phoClSigIEtaIEta_phoE","Photon phoClSigIEtaIEta v phoE;sieie;Energy",250,0,0.025,2000,0,2000);
    hist2d[276] = new TH2D("phoClSigIEtaIEta_phoClTime","Photon phoClSigIEtaIEta v ClTime;sieie;ClTime [ns]",250,0,0.025,jtdiv, -1*jtran, jtran);
    hist2d[278] = new TH2D("phoClSigIEtaIEta_phoSMaj","Photon phoClSigIEtaIEta v phoSMaj;sieie;phoSMaj",250,0,0.025,250,0,2.5);
    hist2d[279] = new TH2D("phoClSigIEtaIEta_phoSMin","Photon phoClSigIEtaIEta v phoSMin;sieie;phoSMin",250,0,0.025,750,0,0.75);
    hist2d[280] = new TH2D("phoClSigIEtaIEta_phoSAlp","Photon phoClSigIEtaIEta v phoSAlp;sieie;phoSAlp",250,0,0.025,350,-1.75,1.75);
    hist2d[281] = new TH2D("phoClSigIEtaIEta_phoCovEtaEta","Photon phoClSigIEtaIEta v phoCovEtaEta;sieie;phoCovEtaEta",250,0,0.025,100,0,0.001);
    hist2d[282] = new TH2D("phoClSigIEtaIEta_phoCovEtaPhi","Photon phoClSigIEtaIEta v phoCovEtaPhi;sieie;phoCovEtaPhi",250,0,0.025,100,-0.0005,0.0005);
    hist2d[283] = new TH2D("phoClSigIEtaIEta_phoCovPhiPhi","Photon phoClSigIEtaIEta v phoCovPhiPhi;sieie;phoCovPhiPhi",250,0,0.025,100,0,0.001);

    hist2d[284] = new TH2D("phoE_phoClTime","Photon phoE v ClTime;Energy;ClTime [ns]",2000,0,2000,jtdiv, -1*jtran, jtran);
    hist2d[286] = new TH2D("phoE_phoSMaj","Photon phoE v phoSMaj;Energy;phoSMaj",2000,0,2000,250,0,2.5);
    hist2d[287] = new TH2D("phoE_phoSMin","Photon phoE v phoSMin;Energy;phoSMin",2000,0,2000,750,0,0.75);
    hist2d[288] = new TH2D("phoE_phoSAlp","Photon phoE v phoSAlp;Energy;phoSAlp",2000,0,2000,350,-1.75,1.75);
    hist2d[289] = new TH2D("phoE_phoCovEtaEta","Photon phoE v phoCovEtaEta;Energy;phoCovEtaEta",2000,0,2000,100,0,0.001);
    hist2d[290] = new TH2D("phoE_phoCovEtaPhi","Photon phoE v phoCovEtaPhi;Energy;phoCovEtaPhi",2000,0,2000,100,-0.0005,0.0005);
    hist2d[291] = new TH2D("phoE_phoCovPhiPhi","Photon phoE v phoCovPhiPhi;Energy;phoCovPhiPhi",2000,0,2000,100,0,0.001);

    hist2d[293] = new TH2D("phoClTime_phoSMaj","Photon phoClTime v phoSMaj;ClTime [ns];phoSMaj",jtdiv, -1*jtran, jtran,250,0,2.5);
    hist2d[294] = new TH2D("phoClTime_phoSMin","Photon phoClTime v phoSMin;ClTime [ns];phoSMin",jtdiv, -1*jtran, jtran,750,0,0.75);
    hist2d[295] = new TH2D("phoClTime_phoSAlp","Photon phoClTime v phoSAlp;ClTime [ns];phoSAlp",jtdiv, -1*jtran, jtran,350,-1.75,1.75);
    hist2d[296] = new TH2D("phoClTime_phoCovEtaEta","Photon phoClTime v phoCovEtaEta;ClTime [ns];phoCovEtaEta",jtdiv, -1*jtran, jtran,100,0,0.001);
    hist2d[297] = new TH2D("phoClTime_phoCovEtaPhi","Photon phoClTime v phoCovEtaPhi;ClTime [ns];phoCovEtaPhi",jtdiv, -1*jtran, jtran,100,-0.0005,0.0005);
    hist2d[298] = new TH2D("phoClTime_phoCovPhiPhi","Photon phoClTime v phoCovPhiPhi;ClTime [ns];phoCovPhiPhi",jtdiv, -1*jtran, jtran,100,0,0.001);

/* avalible numbers -------
	237
	247
	258
	268
	277
	285
	292
------------------------*/

    //hist2d[] = new TH2D("phoClR9_minDr","Photon clR9 Vs minDr;ClusterR9;minDr",100,0,1,500,0,5);
    //hist2d[] = new TH2D("phoR9_ClR9","phoR9_ClR9;R9;clR9",100,0,1,100,0,1);

    //hist2d[] = new TH2D("phoNRH_R9_fake","Photon Fake nClRecHits v R9;nRecHits;R9",200,0,200,100,0,1);
    //hist2d[] = new TH2D("phoNRH_angle","Photon nClRecHits v Angle;nRecHits;Angle",120,0,120,70,0,7);
    //hist2d[] = new TH2D("phoNRH_","Photon nClRecHits v R9;nRecHits;R9",120,0,120,100,0,1);
    //hist2d[] = new TH2D("phoNRH_value","Photon nClRecHits v 2dtValue;nRecHits;2dtValue",120,0,120,50,0.5,1);
    //hist2d[] = new TH2D("phoNRH_slope","Photon nClRecHits v Slope;nRecHits;Slope",120,0,120,sldiv, slmin, slmax);
    //hist2d[ = new TH2D("phoNRH_geoSlope","Photon nClRecHits v geoSlope;nRecHits;geoSlope",120,0,120,sldiv, slmin, slmax);
    //hist2d[] = new TH2D("phoNRH_geoAngle","Photon nClRecHits v geoAngle;nRecHits;geoAngle",120,0,120,70,0,7);

    //hist2d[228] = new TH2D("phoAngle_phoClSMaj","phoAngle_phoClSMaj;2dtAng;Smaj",70,0,7,500,0,5);
    //hist2d[] = new TH2D("phoAngle_phoClSMin","phoAngle_phoClSMin;2dtAng;SMin",70,0,7,200,0,2);
    //hist2d[260] = new TH2D("phoAngle_phoClSigIEtaIEta","phoAngle_phoClSigIEtaIEta;2dtAng;sieie",70,0,7,500,0,0.05);

    //hist2d[261] = new TH2D("phoGeoAngle_phoClSMaj","phoGeoAngle_phoClSMaj;geoAng;SMaj",70,0,7,500,0,5);
    //hist2d[262] = new TH2D("phoGeoAngle_phoClSMin","phoGeoAngle_phoClSMin;geoAng;SMin",70,0,7,200,0,2);
    //hist2d[263] = new TH2D("phoGeoAngle_phoClSigIEtaIEta","phoGeoAngle_phoClSigIEtaIEta;geoAng;sieie",70,0,7,500,0,0.05);

	// eigan plots value/slope/angle/clr9

    //hist2d[230] = new TH2D("phoClR9_value","Photon ClR9 v tValue;ClR9;tValue",100,0,1,50,0.5,1);
    //hist2d[231] = new TH2D("phoClR9_slope","Photon ClR9 v Slope;ClR9;Slope",100,0,1,sldiv, slmin, slmax);
    //hist2d[232] = new TH2D("phoClR9_Angle","Photon ClR9 v Angle;ClR9;Angle",100,0,1,70,0,7);
    //hist2d[233] = new TH2D("phoValue_slope","Photon Value v Slope;Value;Slope",50,0.5,1,sldiv, slmin, slmax);
    //hist2d[235] = new TH2D("phoAngle_value","Photon Angle v tValue;Angle;Value",70,0,7,50,0.5,1);
    //hist2d[234] = new TH2D("phoAngle_slope","Photon Angle v Slope;Angle;Slope",70,0,7,sldiv, slmin, slmax);

    //hist2d[236] = new TH2D("phoGeoValue_angle","Photon geoValue v angle;geoValue;angle",50,0.5,1,70,0,7);
    //hist2d[237] = new TH2D("phoGeoSlope_angle","Photon geoSlope v angle;geoSlope;angle",sldiv, slmin, slmax,70,0,7);
    //hist2d[238] = new TH2D("phoGeoSlope_value","Photon geoSlope v value;geoSlope;value",sldiv, slmin, slmax,50,0.5,1);

    //hist2d[239] = new TH2D("phoGeoAngle_Angle","Photon geoAngle v Angle;geoAngle;rotAngle",70,0,7,70,0,7);
    //hist2d[242] = new TH2D("phoGeoAngle_Value","Photon geoAngle v tValue;geoAngle;tValue",70,0,7,50,0.5,1);
    //hist2d[241] = new TH2D("phoGeoAngle_slope","Photon geoAngle v Slope;geoAngle;slope",70,0,7,sldiv, slmin, slmax);
    //hist2d[243] = new TH2D("phoGeoValue_value","Photon geoValue v tValue;geoValue;tValue",50,0.5,1,50,0.5,1);
    //hist2d[244] = new TH2D("phoGeoValue_slope","Photon geoValue v slope;geoValue;slope",50,0.5,1,sldiv, slmin, slmax);
    //hist2d[240] = new TH2D("phoGeoSlope_slope","Photon geoSlope v slope;geoSlope;slope",sldiv, slmin, slmax,sldiv, slmin, slmax);

    //hist2d[245] = new TH2D("phoGeoAngle_geoValue","Photon geoAngle v geoValue;geoAngle;geoValue",70,0,7,50,0.5,1);
    //hist2d[246] = new TH2D("phoGeoAngle_geoSlope","Photon geoAngle v geoSlope;geoAngle;geoSlope",70,0,7,sldiv, slmin, slmax);
    //hist2d[248] = new TH2D("phoClR9_geoSlope","Photon ClR9 v geoSlope;ClR9;geoSlope",100,0,1,sldiv, slmin, slmax);
    //hist2d[249] = new TH2D("phoClR9_geoAngle","Photon ClR9 v geoAngle;ClR9;geoAngle",100,0,1,70,0,7);
    //hist2d[250] = new TH2D("phoGeoValue_geoSlope","Photon geoValue v geoSlope;geoValue;geoSlope",50,0.5,1,sldiv, slmin, slmax);

    //hist2d[251] = new TH2D("phoClSMaj_2dtValue","phoClSMaj_2dtValue;ClSMaj;2dtValue",500,0,5,50,0.5,1);
    //hist2d[253] = new TH2D("phoClSMaj_3d4Value","phoClSMaj_3d4Value;ClSMaj;3d4Value",500,0,5,50,0.5,1);
    //hist2d[254] = new TH2D("phoClSMin_2dtValue","phoClSMin_2dtValue;ClSMin;2dtValue",200,0,2,50,0.5,1);
    //hist2d[256] = new TH2D("phoClSMin_3d4Value","phoClSMin_3d4Value;ClSMin;3d4Value",200,0,2,50,0.5,1);
    //hist2d[257] = new TH2D("phoClSigIEtaIEta_2dtValue","phoClSigIEtaIEta_2dtValue;ClSigIEtaIEta;2dtValue",500,0,0.05,50,0.5,1);
    //hist2d[259] = new TH2D("phoClSigIEtaIEta_3d4Value","phoClSigIEtaIEta_3d4Value;ClSigIEtaIEta;3d4Value",500,0,0.05,50,0.5,1);
	//60 - 63

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
                //auto indir = "LLPGamma/llpga_GMSB_AOD_v60/"; //argv[1];
                auto indir = "LLPGamma/llpga_GJets_AOD_v60/";

                //auto infilename = "llpgana_mc_AODSIM_GMSB_AOD_v60_Full.txt"; //argv[2];
                auto infilename = "llpgana_mc_AODSIM_GJets_AOD_v60_Full.txt";

				//auto outfilename = "lpgana_mc_test12_tight_max3d4v_sig.root";

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

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v58_V24_10th_Loose_tEB_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v58_V24_10th_Loose_tEB_noPeak_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v58_V24_10th_Loose_tEB_isOOT_Hists.root";				
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v58_V24_10th_Loose_tEB_ClR9r26_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v58_V24_20th_Loose_tEB_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v58_V24_20th_Loose_tEB_noPeak_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v58_V24_20th_Loose_tEB_ClR9r26_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v58_V24_5th_Loose_tEB_isOOT_Hists.root";

				//auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_v58_V24_3DEigan3_20th_Loose_tEB_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_v58_V24_3DEiganSlope_10th_Loose_tEB_Hists.root";
                
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_8th_Loose_tEB_CMB_isNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_8th_Loose_tEB_CMB_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_8th_Loose_tEB_CMB_isSusyNotSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_8th_Loose_tEB_CMB_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV58_mhV29_24th_Loose_tEB_CMB_Hists.root";


                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_Base_10th_Loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_Base_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_Base_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_Base_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30_10th_Loose_tEB_nRH15_isSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA101_10th_Loose_tEB_nRH15_isSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA101_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA101_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30DA101_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA102_10th_Loose_tEB_nRH15_isSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA102_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA102_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30DA102_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA301_10th_Loose_tEB_nRH15_isSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA301_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA301_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30DA301_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA302_10th_Loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA302_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA302_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30DA302_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA502_10th_Loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA502_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA502_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30DA502_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30nRH35_10th_Loose_tEB_nRH15_isSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30nRH35_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30nRH35_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30nRH35_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30mEVRH_10th_Loose_tEB_nRH15_isSig_Hists.root";
				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30mEVRH_10th_Loose_tEB_nRH15_isSusyNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30mEVRH_10th_Loose_tEB_nRH15_isNotSusy_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_Gjets_AOD_aV59_mhV29_GV85E30mEVRH_30th_Loose_tEB_nRH15_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV29_GV85E30DA502nRH35_10th_Loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV59_mhV29_GV85E30DA502nRH35_30th_Loose_tEB_nRH15_isSig_Hists.root";

				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_GV85Pt30logE_4th_Loose_tEB_nRH15_isSig_Hists.root";

				//auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_Pt30logE_10th_Tight_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_Pt30logE_10th_Tight_tEB_nRH15_isNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV59_mhV30_Pt30logE_20th_Tight_tEB_nRH15_isCmb_Hists.root";

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_GV85Pt30logE_10th_Tight_tEB_nRH15_isSig_v3_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV59_mhV30_GV85Pt30logE_10th_Tight_tEB_nRH15_isNotSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV59_mhV30_GV85Pt30logE_20th_Tight_tEB_nRH15_isCmb_v3_Hists.root";
				//v3 has TAngle 33->27

                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_10th_loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_SMM_10th_loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_GEO_10th_loose_tEB_nRH15_isSig_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GMSB_AOD_aV60_mhV31_Pt30logE_SMMGEO_10th_loose_tEB_nRH15_isSig_Hists.root";

                auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_50th_loose_tEB_nRH15_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_SMM_50th_loose_tEB_nRH15_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_GEO_50th_loose_tEB_nRH15_Hists.root";
                //auto outfilename = "llpgana_mc_AODSIM_GJets_AOD_aV60_mhV31_Pt30logE_SMMGEO_20th_loose_tEB_nRH15_Hists.root";

				int pct = 50;
				makehists base;
                base.llpgana_hist_maker( indir, infilename, outfilename, pct );
    //}
    return 1;

}//<<>>int main ( int argc, char *argv[] )

