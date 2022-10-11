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
#include "llpgana_ntuple_base.hh"

#define nHists 256
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.1415926535 // pie ... 

#define CFlt  const float
#define CDbl  const double
#define CVFlt const vector<float>
#define CVDbl const vector<double>

#define DEBUG false
//#define DEBUG true

enum ECAL {EB, EM, EP, NONE};

struct DetIDStruct { 
  	DetIDStruct() {}
  	DetIDStruct(const Int_t ni1, const Int_t ni2, const Int_t nTT, const Int_t & necal) : i1(ni1), i2(ni2), TT(nTT), ecal(necal){}
  	Int_t i1; // EB: iphi, EE: ix
  	Int_t i2; // EB: ieta, EE: iy
  	Int_t TT; // trigger tower
  	Int_t ecal; // EB, EM, EP
};

void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )
{   
    const std::string detIDConfigEB("ecal_config/fullinfo_detids_EB.txt");
    std::ifstream infile( detIDConfigEB, std::ios::in);
    
    UInt_t cmsswId, dbID; 
    Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    TString pos;
    
    while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM)
    {   
        //std::cout << "DetID Input Line: " << cmsswId << " " << iphi << " "  << ieta << " " << 0 << std::endl;
        DetIDMap[cmsswId] = {iphi,ieta,TT25,0};
        auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )
{
    const std::string detIDConfigEE("ecal_config/fullinfo_detids_EE.txt");
    std::ifstream infile( detIDConfigEE, std::ios::in);

    UInt_t cmsswId, dbID;
    Int_t hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    TString EE;

    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant)
    {
        int ec = 1;
        if( side > 0 ) ec = 2;
        //std::cout << "DetID Input Line: " << cmsswId << " " << ix << " "  << iy << " " << ec << std::endl; 
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec};
        auto idinfo = DetIDMap[cmsswId];
        //std::cout << "DetID set to : " << idinfo.i1 << " " << idinfo.i2 << " " << idinfo.ecal << std::endl;
    }
}

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

const auto sq2      (CFlt x){return x*x;}
const auto sq2      (CDbl x){return x*x;}
const auto rad2     (CFlt x, CFlt y, CFlt z = 0.f){return x*x+y*y+z*z;}
const auto hypo     (CFlt x, CFlt y, CFlt z = 0.f){return std::sqrt(rad2(x,y,z));}
const auto phi      (CFlt x, CFlt y){return std::atan2(y,x);}
const auto theta    (CFlt r, CFlt z){return std::atan2(r,z);}
const auto eta      (CFlt x, CFlt y, CFlt z){return -1.0f*std::log(std::tan(theta(hypo(x,y),z)/2.f));}
const auto effMean  (CFlt x, CFlt y){return (x*y)/sqrt(x*x+y*y);}
const auto dltPhi   (CFlt x, CFlt y){auto dp(x-y); if( dp > 180 ){dp-=360.0;} else if( dp < -180 ){ dp+=360.0;} return dp;}
//const auto vfsum    (CVFlt x){return std::accumulate(x.begin(),x.end(),0.0f);}
const auto max      (CVFlt x){float m(x[0]); for(auto ix : x ){ if( ix > m ) m = ix; } return m;}

void fillOUHist1F( float val, float low, float high, float div, TH1F * hist ){

    auto step = ((high-low)/div)/2;
    if( val < low ) hist->Fill( low+step );
    else if ( val > high ) hist->Fill( high-step );
    else hist->Fill( val );

}//<<>>void fillOUHist1F( float val, float low, float high, TH1F & hist )

void fillTH1F( float val, TH1F *& hist ){

   auto nBins = hist->GetNbinsX();
   auto low = hist->GetBinCenter(1);
   auto high = hist->GetBinCenter(nBins);
   if( val < low ) hist->Fill( low );
   else if ( val > high ) hist->Fill( high );
   else hist->Fill( val );

}//<<>>void fillTH1F( float val, TH1F *& hist )

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

class makehists : llpgana_ntuple_base {

	public:

	void llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename );	
	void initHists();
	void getBranches( Long64_t entry );
	void eventLoop( Long64_t entry );
 	void endJobs();	

    TH1D *hist1d[nHists];
    TH2D *hist2d[nHists];

};

//makehists::makehists(){}

//makehists::~makehists(){}

void makehists::llpgana_hist_maker( std::string indir, std::string infilelist, std::string outfilename ){

    //bool debug = true;
    bool debug = false;

    const std::string disphotreename("tree/llpgtree");

    std::ifstream infile(infilelist);
    auto fInTree = new TChain(disphotreename.c_str());
    std::cout << "Adding files to TChain." << std::endl;
    std::string str;
	int cnt = 1;
    while (std::getline(infile,str)){
		if( cnt%4 == 0 ){ 
        auto tfilename = "root://cmseos.fnal.gov//store/user/jaking/" + indir + str;
        std::cout << "--  adding file: " << tfilename << std::endl;
        fInTree->Add(tfilename.c_str());
		}//<<>>if( cnt%4 == 0 ){
		cnt++;
    }//<<>>while (std::getline(infile,str))

	Init(fInTree);
	initHists();

    std::cout << "Setting up For Main Loop." << std::endl;
    auto nEntries = fInTree->GetEntries();
    if(debug) nEntries = 1000;
    //nEntries = 60000;
    //nEntries = 1;
    std::cout << "Proccessing " << nEntries << " entries." << std::endl;
    for (Long64_t centry = 0; centry < nEntries; centry++){
        if( centry%10000 == 0 ) std::cout << "Proccessed " << centry << " of " << nEntries << " entries." << std::endl;
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
	for( int it = 0; it < nHists; it++ ){ if(hist1d[it]) hist1d[it]->Write(); }
    for( int it = 0; it < nHists; it++ ){ if(hist2d[it]) hist2d[it]->Write(); }
    for( int it = 0; it < nHists; it++ ){ if(hist1d[it]) delete hist1d[it]; if(hist2d[it]) delete hist2d[it]; }

    fOutFile->Close();
    std::cout << "llpgana_hist_maker : Thats all Folks!!" << std::endl;
}//<<>>void llpgana_hist_maker

void makehists::eventLoop( Long64_t entry ){

    float   jetPTmin        = 100.0;// for energy/time comp
    int     jetIDmin        = 2; //3;
    float   jetETAmax       = 1.5; //1.5;
    int    minRHcnt        = 5; //32;
    float   minRHenr        = 2.0;
    float   bcMinEnergy     = 0.667;
    int    bcMinRHGrpSize  = 3;
    float   minEmf          = 0.0;//0.2

    if( DEBUG ) std::cout << "Finding Events" << std::endl;
	//------------ event varibles ------------------

    b_run->GetEntry(entry);   //!
    b_lumi->GetEntry(entry);   //!
    b_event->GetEntry(entry);   //!
    b_nVtx->GetEntry(entry);   //!
    b_vtxX->GetEntry(entry);   //!
    b_vtxY->GetEntry(entry);   //!
    b_vtxZ->GetEntry(entry);   //!

    if( DEBUG ) std::cout << "Finding rechits" << std::endl;
	//------------ rechits -------------------------

    b_nRecHits->GetEntry(entry);   //!
    b_rhPosX->GetEntry(entry);   //!
    b_rhPosY->GetEntry(entry);   //!
    b_rhPosZ->GetEntry(entry);   //!
    b_rhPosEta->GetEntry(entry);   //!
    b_rhPosPhi->GetEntry(entry);   //!
    b_rhEnergy->GetEntry(entry);   //!
    b_rhTime->GetEntry(entry);   //!
    b_rhTimeErr->GetEntry(entry);   //!
    b_rhTOF->GetEntry(entry);   //!
    b_rhID->GetEntry(entry);   //!
    b_rhXtalI1->GetEntry(entry);   //!
    b_rhXtalI2->GetEntry(entry);   //!
    b_rhSubdet->GetEntry(entry);   //!
    b_rhisOOT->GetEntry(entry);   //!
    b_rhisGS6->GetEntry(entry);   //!
    b_rhisGS1->GetEntry(entry);   //!
    b_rhadcToGeV->GetEntry(entry);   //!
    b_rhped12->GetEntry(entry);   //!
    b_rhped6->GetEntry(entry);   //!
    b_rhped1->GetEntry(entry);   //!
    b_rhpedrms12->GetEntry(entry);   //!
    b_rhpedrms6->GetEntry(entry);   //!
    b_rhpedrms1->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " -- Finding rechits" << std::endl;

    for( int it = 0; it < nRecHits; it++ ){

		if( (*rhSubdet)[it] == 0 ){
		fillTH1( (*rhTime)[it], hist1d[131]); 
		fillTH1( (*rhEnergy)[it], hist1d[132]); 
		hist2d[122]->Fill( (*rhTime)[it], (*rhEnergy)[it]);
		fillTH1( (*rhisOOT)[it], hist1d[135]);
		} else {
		fillTH1( (*rhTime)[it], hist1d[133]); 
		fillTH1( (*rhEnergy)[it], hist1d[134]); 
		hist2d[123]->Fill( (*rhTime)[it], (*rhEnergy)[it]);
		fillTH1( (*rhisOOT)[it], hist1d[136]);
		}

	}//<<>>for( int it = 0; it < nRecHits; it++ )

    if( DEBUG ) std::cout << "Finding calojets" << std::endl;
	//------------  calojets ------------------------

    b_cljBcCnt->GetEntry(entry);   //!
    b_nCaloJets->GetEntry(entry);   //!
    b_cljSeedTOFTime->GetEntry(entry);   //!
    b_cljCMeanTime->GetEntry(entry);   //!
    b_cljBc3dEx->GetEntry(entry);   //!
    b_cljBc3dEy->GetEntry(entry);   //!
    b_cljBc3dEz->GetEntry(entry);   //!
    b_cljBc3dEv->GetEntry(entry);   //!
    b_cljBc3dEslope->GetEntry(entry);   //!
    b_cljBc3dEchisp->GetEntry(entry);   //!
    b_cljBc2dEx->GetEntry(entry);   //!
    b_cljBc2dEy->GetEntry(entry);   //!
    b_cljBc2dEv->GetEntry(entry);   //!
    b_cljBc2dEslope->GetEntry(entry);   //!
    b_cljBc2dEchisp->GetEntry(entry);   //!
    b_cljCDrMeanTime->GetEntry(entry);   //!
    b_cljPt->GetEntry(entry);   //!
    b_cljEnergy->GetEntry(entry);   //!
    b_cljPhi->GetEntry(entry);   //!
    b_cljEta->GetEntry(entry);   //!
    b_cljPx->GetEntry(entry);   //!
    b_cljPy->GetEntry(entry);   //!
    b_cljPz->GetEntry(entry);   //!

	for( int it = 0; it < nCaloJets; it++ ){

    	hist1d[146]->Fill( (*cljCMeanTime)[it] );//c mean
        hist1d[158]->Fill( (*cljCDrMeanTime)[it] );//c dr mean
    	hist1d[147]->Fill( (*cljSeedTOFTime)[it] );//lead time 
    	hist1d[148]->Fill( (*cljCMeanTime)[it] - (*cljSeedTOFTime)[it] );//diff

	}//<<>>for( int it = 0; it < nCaloJets; it++ )

    if( DEBUG ) std::cout << "Finding photons" << std::endl;
	//----------------- photons ------------------

    b_nPhotons->GetEntry(entry);   //!
    b_phoSeedTOFTime->GetEntry(entry);   //!
    b_phoCMeanTime->GetEntry(entry);   //!
    b_phoSc3dEx->GetEntry(entry);   //!
    b_phoSc3dEy->GetEntry(entry);   //!
    b_phoSc3dEz->GetEntry(entry);   //!
    b_phoSc3dEv->GetEntry(entry);   //!
    b_phoSc3dEslope->GetEntry(entry);   //!
    b_phoSc3dEchisp->GetEntry(entry);   //!
    b_phoSc2dEx->GetEntry(entry);   //!
    b_phoSc2dEy->GetEntry(entry);   //!
    b_phoSc2dEv->GetEntry(entry);   //!
    b_phoSc2dEslope->GetEntry(entry);   //!
    b_phoSc2dEchisp->GetEntry(entry);   //!
    b_phoPt->GetEntry(entry);   //!
    b_phoEnergy->GetEntry(entry);   //!
    b_phoPhi->GetEntry(entry);   //!
    b_phoEta->GetEntry(entry);   //!
    b_phoPx->GetEntry(entry);   //!
    b_phoPy->GetEntry(entry);   //!
    b_phoPz->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " - loopoing phtons" << std::endl;
    for( int it = 0; it < nPhotons; it++ ){

        hist1d[155]->Fill( (*phoCMeanTime)[it] );//c mean
        hist1d[156]->Fill( (*phoSeedTOFTime)[it] );//lead time 
        hist1d[157]->Fill( (*phoCMeanTime)[it] - (*phoSeedTOFTime)[it] );//diff

        if( (*phoSc3dEx)[it] != -9 ){
            auto epanlge = getAngle( (*phoSc3dEx)[it], (*phoSc3dEy)[it] );
            hist1d[63]->Fill(epanlge);//etaphi angle
            auto ephypo3D = hypo( (*phoSc3dEx)[it], (*phoSc3dEy)[it] );
            auto etanlge = getAngle( ephypo3D, (*phoSc3dEz)[it] );
            hist1d[64]->Fill(etanlge);//etatim angle
            hist2d[82]->Fill( (*phoSc3dEx)[it], (*phoSc3dEy)[it] );
            hist2d[83]->Fill( (*phoSc3dEx)[it], (*phoSc3dEz)[it] );
            hist1d[86]->Fill( (*phoSc3dEv)[it] );
        }//<<>>if( phoSCEigen3D[0] != -999 )
        if( (*phoSc2dEx)[it] != -9 ){
            auto sphanlge = getAngle( (*phoSc2dEx)[it], (*phoSc2dEy)[it] );
            hist1d[65]->Fill( sphanlge );//eliptical angle
            hist2d[81]->Fill( (*phoSc2dEx)[it], (*phoSc2dEy)[it] );
            hist1d[81]->Fill( (*phoSc2dEv)[it] );
            //hist2d[88]->Fill( sphanlge, (*phoSc2dEv)[it] );
        }//<<>>if( phoSCEigen2D[0] != -999 )i

    }//<<>>for( int it = 0; it < nPhotons; it++ )

    if( DEBUG ) std::cout << "Finding ootphotons" << std::endl;
	//-------------- OOTPhotons ----------------------------

    b_nOotPhotons->GetEntry(entry);   //!
    b_ootPhoSeedTOFTime->GetEntry(entry);   //!
    b_ootPhoCMeanTime->GetEntry(entry);   //!
    b_ootPhoSc3dEx->GetEntry(entry);   //!
    b_ootPhoSc3dEy->GetEntry(entry);   //!
    b_ootPhoSc3dEz->GetEntry(entry);   //!
    b_ootPhoSc3dEv->GetEntry(entry);   //!
    b_ootPhoSc3dEslope->GetEntry(entry);   //!
    b_ootPhoSc3dEchisp->GetEntry(entry);   //!
    b_ootPhoSc2dEx->GetEntry(entry);   //!
    b_ootPhoSc2dEy->GetEntry(entry);   //!
    b_ootPhoSc2dEv->GetEntry(entry);   //!
    b_ootPhoSc2dEslope->GetEntry(entry);   //!
    b_ootPhoSc2dEchisp->GetEntry(entry);   //!
    b_ootPhoPt->GetEntry(entry);   //!
    b_ootPhoEnergy->GetEntry(entry);   //!
    b_ootPhoPhi->GetEntry(entry);   //!
    b_ootPhoEta->GetEntry(entry);   //!
    b_ootPhoPx->GetEntry(entry);   //!
    b_ootPhoPy->GetEntry(entry);   //!
    b_ootPhoPz->GetEntry(entry);   //!

    for( int it = 0; it < nOotPhotons; it++ ){

        hist1d[149]->Fill( (*ootPhoCMeanTime)[it] );//c mean
        hist1d[150]->Fill( (*ootPhoSeedTOFTime)[it] );//lead time 
        hist1d[151]->Fill( (*ootPhoCMeanTime)[it] - (*ootPhoSeedTOFTime)[it] );//diff

	}//<<>>for( int it = 0; it < nOotPhotons; it++ )

    if( DEBUG ) std::cout << "Finding electrons" << std::endl;
	//-------- electrons --------------------------------------
	
    b_nElectrons->GetEntry(entry);   //!
    b_eleSeedTOFTime->GetEntry(entry);   //!
    b_eleCMeanTime->GetEntry(entry);   //!
    b_eleSc3dEx->GetEntry(entry);   //!
    b_eleSc3dEy->GetEntry(entry);   //!
    b_eleSc3dEz->GetEntry(entry);   //!
    b_eleSc3dEv->GetEntry(entry);   //!
    b_eleSc3dEslope->GetEntry(entry);   //!
    b_eleSc3dEchisp->GetEntry(entry);   //!
    b_eleSc2dEx->GetEntry(entry);   //!
    b_eleSc2dEy->GetEntry(entry);   //!
    b_eleSc2dEv->GetEntry(entry);   //!
    b_eleSc2dEslope->GetEntry(entry);   //!
    b_eleSc2dEchisp->GetEntry(entry);   //!
    b_elePt->GetEntry(entry);   //!
    b_eleEnergy->GetEntry(entry);   //!
    b_elePhi->GetEntry(entry);   //!
    b_eleEta->GetEntry(entry);   //!
    b_elePx->GetEntry(entry);   //!
    b_elePy->GetEntry(entry);   //!
    b_elePz->GetEntry(entry);   //!

    for( int it = 0; it < nElectrons; it++ ){

        hist1d[152]->Fill( (*eleCMeanTime)[it] );//c mean
        hist1d[153]->Fill( (*eleSeedTOFTime)[it] );//lead time 
        hist1d[154]->Fill( (*eleCMeanTime)[it] - (*eleSeedTOFTime)[it] );//diff

    }//<<>>for( int it = 0; it < nOotPhotons; it++ )

    if( DEBUG ) std::cout << "Finding pfjets" << std::endl;
	//--------- jets --------------------------
	
	//// ---  base info  ----

    b_jetHt->GetEntry(entry);   //!
    b_nJets->GetEntry(entry);   //!
    b_nGoodDrJets->GetEntry(entry);   //!
    b_nGoodScJets->GetEntry(entry);   //!
    b_nGoodBcJets->GetEntry(entry);   //!
    b_nUnJets->GetEntry(entry);   //!
    b_jetE->GetEntry(entry);   //!
    b_jetPt->GetEntry(entry);   //!
    b_jetEta->GetEntry(entry);   //!
    b_jetPhi->GetEntry(entry);   //!
    b_jetNHF->GetEntry(entry);   //!
    b_jetNEMF->GetEntry(entry);   //!
    b_jetCHF->GetEntry(entry);   //!
    b_jetCEMF->GetEntry(entry);   //!
    b_jetMUF->GetEntry(entry);   //!
    b_jetNHM->GetEntry(entry);   //!
    b_jetCHM->GetEntry(entry);   //!
    b_jetPHM->GetEntry(entry);   //!
    b_jetELM->GetEntry(entry);   //!
    b_jetPHE->GetEntry(entry);   //!
    b_jetPHEF->GetEntry(entry);   //!
    b_jetELE->GetEntry(entry);   //!
    b_jetELEF->GetEntry(entry);   //!
    b_jetMUE->GetEntry(entry);   //!
    b_jetID->GetEntry(entry);   //!

    b_jetSumEPFrac->GetEntry(entry);   //!
    b_jetEPEnergy->GetEntry(entry);   //!
    b_jetEMEnergy->GetEntry(entry);   //!
    b_jetEMEnrFrac->GetEntry(entry);   //!
    b_jetEPEnrFrac->GetEntry(entry);   //!

    for( int it = 0; it < nJets; it++ ){

        const auto jetepafrac   = (*jetPHEF)[it] + (*jetELEF)[it];;
        const auto jetepe       = (*jetPHE)[it] + (*jetELE)[it];
        const auto jeteme       = (*jetCEMF)[it] + (*jetNEMF)[it];
        const auto jetemfrac    = jeteme/(*jetE)[it];
        const auto jetepfrac    = jetepe/(*jetE)[it];
		
        hist2d[61]->Fill(jetepafrac,jetepfrac);
        hist2d[62]->Fill(jetepfrac,jetemfrac);

        fillTH1((*jetPt)[it],hist1d[12]);//hist1d[12]->Fill(jet.pt());
        fillTH1((*jetPhi)[it],hist1d[13]);//hist1d[13]->Fill(jet.phi());
        fillTH1((*jetEta)[it],hist1d[14]);//hist1d[14]->Fill(jet.eta());

	}//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding genjet" << std::endl;
	//// --- genjet info -----------------------------------

    b_jetGenImpactAngle->GetEntry(entry);   //!
    b_jetGenTime->GetEntry(entry);   //!
    b_jetGenPt->GetEntry(entry);   //!
    b_jetGenEta->GetEntry(entry);   //!
    b_jetGenEnergy->GetEntry(entry);   //!
    b_jetGenEMFrac->GetEntry(entry);   //!
    b_jetGenDrMatch->GetEntry(entry);   //!
    b_jetGenTimeVar->GetEntry(entry);   //!
    b_jetGenTimeLLP->GetEntry(entry);   //!
    b_jetGenLLPPurity->GetEntry(entry);   //!
    b_jetGenNextBX->GetEntry(entry);   //!
    b_jetGenNKids->GetEntry(entry);   //!
    b_jetGenTOF->GetEntry(entry);   //!

    for( int it = 0; it < nJets; it++ ){

        hist2d[109]->Fill((*jetGenTime)[it],(*jetGenTOF)[it]);
        hist1d[110]->Fill((*jetGenTime)[it]);
    	hist1d[111]->Fill((*jetGenTOF)[it]);

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet dr times" << std::endl;
	//// ---  jet time dr method --------------------------- (*)[it]

    b_jetDrLeadEta->GetEntry(entry);   //!
    b_jetDrLeadPhi->GetEntry(entry);   //!
    b_jetDrLeadEnr->GetEntry(entry);   //!
    b_sJetDrRHEnergy->GetEntry(entry);   //!
    b_jetDrEMF->GetEntry(entry);   //!
    b_jetDrRhCnt->GetEntry(entry);   //!
    b_jetDRMuTime->GetEntry(entry);   //!
    b_jetDRTimeError->GetEntry(entry);   //!
    b_jetDRTimeRMS->GetEntry(entry);   //!
    b_jetDRMedTime->GetEntry(entry);   //!
    b_jetCDRMuTime->GetEntry(entry);   //!
    b_jetCDRMedTime->GetEntry(entry);   //!

	for( int it = 0; it < nJets; it++ ){

        fillTH1((*jetDRMuTime)[it],hist1d[29]);//hist1d[29]->Fill(jmutime);
        fillTH1((*jetDrRhCnt)[it],hist1d[1]);//hist1d[1]->Fill(rhCount);
        fillTH1((*jetDRTimeError)[it],hist1d[2]);//hist1d[2]->Fill(jterr);
        fillTH1((*jetDRTimeRMS)[it],hist1d[3]);//hist1d[3]->Fill(jtrms);
        fillTH1((*jetDRMedTime)[it],hist1d[4]);//hist1d[4]->Fill(jmedtime);
        fillTH1((*jetCDRMuTime)[it],hist1d[6]);//hist1d[6]->Fill(jcmutime);
        fillTH1((*jetCDRMedTime)[it],hist1d[7]);//hist1d[7]->Fill(jcmedtime);

        hist2d[1]->Fill((*jetDRMuTime)[it],(*jetPt)[it]);
        hist2d[2]->Fill((*jetDRMuTime)[it],(*jetID)[it]);
        hist2d[3]->Fill((*jetDRMuTime)[it],(*jetNHF)[it]);//jetNHF
        hist2d[4]->Fill((*jetDRMuTime)[it],(*jetCHF)[it]);//jetCHF
        hist2d[5]->Fill((*jetDRMuTime)[it],(*jetNEMF)[it]);//jetNEMF
        hist2d[6]->Fill((*jetDRMuTime)[it],(*jetCEMF)[it]);//jetCEMF
        hist2d[7]->Fill((*jetDRMuTime)[it],(*jetMUF)[it]);//jetMUF
        hist2d[8]->Fill((*jetDRMuTime)[it],(*jetNHM)[it]);//jetNHM
        hist2d[9]->Fill((*jetDRMuTime)[it],(*jetCHM)[it]);//jetCHM
    
        hist2d[10]->Fill((*jetDRMuTime)[it],(*jetDRMedTime)[it]);
        hist2d[24]->Fill((*jetDRMuTime)[it],(*jetDrRhCnt)[it]);
        hist2d[25]->Fill((*jetDRMedTime)[it],(*jetDrRhCnt)[it]);
        hist2d[11]->Fill((*jetDRMuTime)[it],(*jetDRTimeRMS)[it]);
        hist2d[12]->Fill((*jetDRMuTime)[it],(*jetDRTimeError)[it]);

        hist2d[32]->Fill((*jetDRMuTime)[it],(*jetDrLeadEta)[it]);
        hist2d[33]->Fill((*jetDRMuTime)[it],(*jetDrLeadPhi)[it]);
        hist2d[34]->Fill((*jetDRMuTime)[it],(*jetDrLeadEnr)[it]);
        hist2d[35]->Fill((*jetDRMedTime)[it],(*jetDrLeadEta)[it]);
        hist2d[36]->Fill((*jetDRMedTime)[it],(*jetDrLeadPhi)[it]);
        hist2d[37]->Fill((*jetDRMedTime)[it],(*jetDrLeadEnr)[it]);
    

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

    b_nJetScMatch->GetEntry(entry);   //!
    b_jetScRhCnt->GetEntry(entry);   //!
    b_sJetScEnergy->GetEntry(entry);   //!
    b_sJetScPhEnergy->GetEntry(entry);   //!
    b_sJetScRhEnergy->GetEntry(entry);   //!
    b_jetScEMF->GetEntry(entry);   //!
    b_jetSCMuTime->GetEntry(entry);   //!
    b_jetSCMedTime->GetEntry(entry);   //!
    b_jetCSCMuTime->GetEntry(entry);   //!
    b_jetCSCMedTime->GetEntry(entry);   //!
	b_jetImpactAngle->GetEntry(entry);   //!

    b_jetSc3dEx->GetEntry(entry);   //!
    b_jetSc3dEy->GetEntry(entry);   //!
    b_jetSc3dEz->GetEntry(entry);   //!
    b_jetSc3dEv->GetEntry(entry);   //!
    b_jetSc3dEslope->GetEntry(entry);   //!
    b_jetSc3dEchisp->GetEntry(entry);   //!
    b_jetSc2dEx->GetEntry(entry);   //!
    b_jetSc2dEy->GetEntry(entry);   //!
    b_jetSc2dEv->GetEntry(entry);   //!
    b_jetSc2dEslope->GetEntry(entry);   //!
    b_jetSc2dEchisp->GetEntry(entry);   //!
    b_jetSc2dEslope2->GetEntry(entry);   //!
    b_jetSc2dEchisp2->GetEntry(entry);   //!
    b_jetSc2dErangle->GetEntry(entry);   //!
    b_jetSc2dEnxsum->GetEntry(entry);   //!

    if( DEBUG ) std::cout << " - Finding jet sc times " << std::endl;
    for( int it = 0; it < nJets; it++ ){

        if( (*jetSc3dEy)[it] != -9 ){
                auto epanlge = getAngle( (*jetSc3dEx)[it], (*jetSc3dEy)[it] );
                hist1d[63]->Fill(epanlge);//etaphi angle
                auto ephypo3D = hypo( (*jetSc3dEx)[it], (*jetSc3dEy)[it] );
                auto etanlge = getAngle( ephypo3D, (*jetSc3dEz)[it] );
                fillTH1(etanlge,hist1d[64]);//hist1d[64]->Fill(etanlge);//etatim angle
                hist2d[82]->Fill( (*jetSc3dEx)[it], (*jetSc3dEy)[it] );
                hist2d[83]->Fill( (*jetSc3dEx)[it], (*jetSc3dEz)[it] );
                hist1d[86]->Fill((*jetSc3dEv)[it]);
                if( (*jetSc3dEchisp)[it] > 0.95 && (*jetSc3dEv)[it] < 0.9 && (*jetSc3dEv)[it] > 0.7 ){
                    hist2d[96]->Fill( (*jetImpactAngle)[it], (*jetSc3dEslope)[it] );
                    hist2d[126]->Fill( (*jetEta)[it], (*jetSc3dEslope)[it] );
                    hist2d[127]->Fill( (*jetPhi)[it], (*jetSc3dEslope)[it] );
                }//<<>>if( (*jetSc3dEchisp)[it] > 0.95 ):
            }//<<>>if( (*jetSc3dEx)[it] != -999 )
            if( (*jetSc2dEx)[it] != -9 ){
                auto sphanlge = getAngle( (*jetSc2dEx)[it], (*jetSc2dEy)[it] );
                hist1d[65]->Fill(sphanlge);//eliptical angle
                hist2d[81]->Fill( (*jetSc2dEx)[it], (*jetSc2dEy)[it] );
                hist1d[81]->Fill((*jetSc2dEv)[it]);
                //hist2d[88]->Fill( sphanlge, (*jetSc2dEv)[it] );
                if( (*jetSc2dEchisp)[it] > 0.95 && (*jetSc2dEv)[it] < 0.9 && (*jetSc2dEv)[it] > 0.7 ){
                    if( (*jetEta)[it] > 1.0 ) hist1d[117]->Fill( (*jetSc2dEslope)[it] );
                    if( (*jetEta)[it] < 0.5 ) hist1d[118]->Fill( (*jetSc2dEslope)[it] );
                    hist2d[89]->Fill( (*jetEta)[it], (*jetSc2dEslope)[it] );
                    hist2d[121]->Fill( (*jetPhi)[it], (*jetSc2dEslope)[it] );
                    hist2d[90]->Fill( (*jetImpactAngle)[it], (*jetSc2dEslope)[it] );
                }//<<>>if( (*jetSc2dEchisp)[it] < 0.1 )
                //hist2d[91]->Fill( (*jetSc2dEslope)[it], (*jetSc2dEchisp)[it]);
        }//<<>>if( (*jetSc2dEx)[it] != -999 )

		if( DEBUG ) std::cout << " -- Filling jet sc times " << std::endl;

        hist2d[39]->Fill( (*jetDrRhCnt)[it], (*jetScRhCnt)[it] );
        hist1d[44]->Fill((*jetSCMedTime)[it]);//median
        hist1d[45]->Fill((*jetSCMuTime)[it]);//mean
        //hist1d[46]->Fill(jetSCTimeStats[4]);//rms
        //hist1d[50]->Fill(jetSCTimeStats[5]);//skew
        hist1d[8]->Fill((*jetCSCMuTime)[it]);//c mean
        hist1d[9]->Fill((*jetCSCMedTime)[it]);//c med   

        hist1d[89]->Fill((*jetCSCMuTime)[it]-(*jetGenTime)[it]);

        hist2d[53]->Fill( (*jetScEMF)[it], (*jetEMEnrFrac)[it] );
        hist2d[56]->Fill( (*sJetScRhEnergy)[it], (*jetEMEnergy)[it] );
        hist2d[59]->Fill( (*jetDrEMF)[it], (*jetScEMF)[it] );

        //hist2d[63]->Fill( (*jetEta)[it], jetSCTimeStats[7] );
        //hist2d[64]->Fill( (*jetEta)[it]etaMmt, jetSCTimeStats[7] );
        //hist2d[65]->Fill( (*jetPhi)[it]phiMnt, jetSCTimeStats[7] );
        //hist2d[66]->Fill( (*jetEta)[it]phiMnt, jetSCTimeStats[7] );
        //hist2d[67]->Fill( jetMaxD, jetSCTimeStats[7] );
        //hist2d[68]->Fill( jetConPtDis, jetSCTimeStats[7] );
        //hist2d[69]->Fill( jetConEtaPhiSprd, jetSCTimeStats[7] );
        //hist2d[70]->Fill( jetArea, jetSCTimeStats[7] );
        //hist2d[71]->Fill( jetNCarry, jetSCTimeStats[7] );
        //hist2d[72]->Fill( jetNConst, jetSCTimeStats[7] );

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet bc times" << std::endl;
	//// --- jet time BC method ---------------------------

    b_jetBcTimesCnt->GetEntry(entry);   //!
    b_jetBcSumRHEnr->GetEntry(entry);   //!
    b_jetBcEMFr->GetEntry(entry);   //!
    b_jetBcRhCnt->GetEntry(entry);   //!
    b_jetBcGrpCnt->GetEntry(entry);   //!
    b_jetCBCMuTime->GetEntry(entry);   //!
    b_jetCBCMedTime->GetEntry(entry);   //!

    for( int it = 0; it < nJets; it++ ){

        hist2d[51]->Fill((*jetDrRhCnt)[it],(*jetBcRhCnt)[it]);
        hist1d[55]->Fill((*jetBcGrpCnt)[it]);
        hist1d[10]->Fill((*jetCBCMedTime)[it]);//c med
        hist1d[19]->Fill((*jetCBCMuTime)[it]);//c mu
        hist2d[54]->Fill( (*jetBcEMFr)[it], (*jetEMEnrFrac)[it] );
        hist2d[57]->Fill( (*jetBcSumRHEnr)[it], (*jetEMEnergy)[it] );
        hist2d[58]->Fill( (*jetBcSumRHEnr)[it], (*sJetScRhEnergy)[it] );
        hist2d[60]->Fill( (*jetScEMF)[it], (*jetBcEMFr)[it] );

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jet genjet comps" << std::endl;
    //****************************  photon/electron to kid pfcand -> SC matcher ********************************

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

            hist1d[112]->Fill((*jetGenTimeLLP)[it]);
            hist1d[113]->Fill((*jetGenLLPPurity)[it]);
            hist2d[114]->Fill((*jetGenLLPPurity)[it],(*jetGenTimeVar)[it]);
            hist1d[116]->Fill((*jetGenNKids)[it]);
            hist2d[115]->Fill((*jetGenTimeVar)[it],(*jetGenNKids)[it]);
            hist2d[116]->Fill((*jetGenLLPPurity)[it],(*jetGenNKids)[it]);
            hist1d[93]->Fill((*jetGenTime)[it]);
            hist1d[90]->Fill((*jetGenImpactAngle)[it]);
            hist1d[105]->Fill((*jetGenDrMatch)[it]);
            hist1d[108]->Fill((*jetGenTimeVar)[it]);
            hist1d[109]->Fill((*jetGenNextBX)[it]);
            hist1d[106]->Fill( difSCTime );
            hist1d[107]->Fill( difDrTime );
            hist2d[101]->Fill( jetERatio, (*jetGenTime)[it] );
            hist2d[117]->Fill( difSCTime, (*jetGenDrMatch)[it] );
            hist2d[118]->Fill( (*jetGenTime)[it], (*jetEMEnrFrac)[it] );
            hist2d[124]->Fill( (*jetEMEnrFrac)[it], (*jetGenDrMatch)[it] );
            hist2d[125]->Fill( jetERatio, (*jetGenDrMatch)[it] );

        //}//<<>>if( genSpaceCut )

        //if( hasGoodGenSCMatch && etaCut && genEnergyCut && genVarCut && hasGoodGenTime && hasGoodSCTime ){

            hist2d[98]->Fill( (*jetE)[it], (*jetGenEnergy)[it] );
            hist2d[99]->Fill( jetERatio, (*jetGenTime)[it] );
            hist2d[110]->Fill( jetERatio, (*jetCSCMuTime)[it] );
            hist2d[111]->Fill( jetERatio, (*jetCDRMuTime)[it] );
            hist2d[100]->Fill( (*jetEMEnrFrac)[it], (*jetCSCMuTime)[it] );
            hist2d[103]->Fill( jetERatio, difSCTime );
            hist2d[104]->Fill( (*jetCDRMuTime)[it], (*jetCSCMuTime)[it] );
            hist2d[102]->Fill( (*jetGenTime)[it], (*jetCDRMuTime)[it] );
            hist2d[105]->Fill( (*jetGenTime)[it], (*jetCSCMuTime)[it] );
            hist2d[106]->Fill( (*jetGenTime)[it], (*jetGenEnergy)[it] );
            hist2d[107]->Fill( (*jetGenDrMatch)[it], difSCTime );
            hist2d[108]->Fill( (*jetGenDrMatch)[it], difDrTime );
            hist2d[112]->Fill( (*jetGenTimeVar)[it], difSCTime );
            hist2d[113]->Fill( (*jetGenLLPPurity)[it], difSCTime );

        }//<<>>if( (*jetCSCMuTime)[it] > -28.0 && (*jetGenTime)[it] > -28.0 )

    }//<<>>for( int it = 0; it < nJets; it++ )

    if( DEBUG ) std::cout << "Finding jets info" << std::endl;
	//// --- jet time resolution -------------------------

    hist1d[17]->Fill(jetHt);
    hist1d[30]->Fill(nGoodDrJets);
    hist1d[31]->Fill(nGoodScJets);
    hist1d[32]->Fill(nGoodBcJets);
    hist1d[33]->Fill(nUnJets);
    if( nUnJets != 0 ) hist1d[34]->Fill(float(nJets)/nUnJets);
    if( nGoodScJets != 0 ) hist1d[38]->Fill(float(nGoodBcJets)/nGoodScJets);
    if( nJets != 0 ){
        hist1d[35]->Fill(float(nGoodDrJets)/nJets);
        hist1d[36]->Fill(float(nGoodScJets)/nJets);
        hist1d[37]->Fill(float(nGoodBcJets)/nJets);
    }//<<>>if( nUnJets != 0 )

    if( DEBUG ) std::cout << "Finding jet dt pairs" << std::endl;
    //-----------------------------------------------------------------------------------------------------
    // ***************************** d jetTime for back-to-back high pt jets  *****************************
    auto dijetIdCut = 1;
    auto dijetPtMin = 200.0;
    auto difPtLmt = 0.8;
    auto htPctLmt = 0.8;
    auto dPhiLmt = 2.8;


    for( int q = 0; q < nJets; q++ ) hist1d[53]->Fill( getdt((*jetCSCMuTime)[q],(*jetCBCMuTime)[q]) );

    if( DEBUG ) std::cout << "Finding jet dt pairs" << std::endl;
    for ( int q = 0; q < nJets; q++ ){
        for ( int p = q+1; p < nJets; p++ ){

            if( DEBUG ) std::cout << " - filter jet pairs" << std::endl;
            if( (*jetPt)[q] < dijetPtMin ) continue;
            auto diffPt = (*jetPt)[p]/(*jetPt)[q];
            hist1d[24]->Fill(diffPt);
            if( diffPt < difPtLmt ) continue;
            auto htPct= ((*jetPt)[q]+(*jetPt)[p])/jetHt;
            hist1d[25]->Fill(htPct);
            if( htPct < htPctLmt ) continue;
            //auto dPhi = reco::deltaPhi((*jetPhi)[q],(*jetPhi)[p]);
            auto dPhi = dltPhi((*jetPhi)[q],(*jetPhi)[p]); 
            hist1d[26]->Fill(dPhi);
            if( dPhi < dPhiLmt ) continue;

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
            if( dTmu > dtThrs ) hist1d[15]->Fill(dTmu);
            if( dTmed > dtThrs ) hist1d[16]->Fill(dTmed);
            if( dTcmu > dtThrs ) hist1d[39]->Fill(dTcmu);
            if( dTcmed > dtThrs ) hist1d[40]->Fill(dTcmed);

            if( dTmedsc > dtThrs ) hist1d[48]->Fill(dTmedsc);
            if( dTmusc > dtThrs ) hist1d[49]->Fill(dTmusc);
            if( dTcmusc > dtThrs ) hist1d[41]->Fill(dTcmusc);
            if( dTcmedsc > dtThrs ) hist1d[42]->Fill(dTcmedsc);

            if( dTcmubc > dtThrs ) hist1d[27]->Fill(dTcmubc);
            if( dTcmedbc > dtThrs ) hist1d[28]->Fill(dTcmedbc);

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
            hist2d[43]->Fill(dTmusc,effje);
            hist2d[44]->Fill(dTmu,effje);

        }//<<>>for ( int p = q+1; p < nJets; p++ )
    }//<<>>for ( int q = 0; q < nJets; q++ )

    if( DEBUG ) std::cout << "Finding calojet dt pairs" << std::endl;
    for ( int q = 0; q < nCaloJets; q++ ){
        for ( int p = q+1; p < nCaloJets; p++ ){

            if( DEBUG ) std::cout << " - filter calojet pairs" << std::endl;
            if( (*cljPt)[q] < dijetPtMin ) continue;
            auto diffPt = (*cljPt)[p]/(*cljPt)[q];
            hist1d[159]->Fill(diffPt);
            if( diffPt < difPtLmt ) continue;
            //auto htPct= ((*cljPt)[q]+(*cljPt)[p])/jetHt;
            //hist1d[25]->Fill(htPct);
            //if( htPct < htPctLmt ) continue;
            //auto dPhi = reco::deltaPhi((*jetPhi)[q],(*jetPhi)[p]);
            auto dPhi = dltPhi((*cljPhi)[q],(*cljPhi)[p]);
            hist1d[160]->Fill(dPhi);
            if( dPhi < dPhiLmt ) continue;

            if( DEBUG ) std::cout << " - get calojet pair dt" << std::endl;
            auto dTmu = getdt( (*cljCMeanTime)[q], (*cljCMeanTime)[p] );

            auto dtThrs = -2.5;//removes default dt values from getdt
            if( dTmu > dtThrs ) hist1d[161]->Fill(dTmu);

        }//<<>>for ( int p = q+1; p < nJets; p++ )
    }//<<>>for ( int q = 0; q < nJets; q++ )

}

void makehists::getBranches( Long64_t entry ){

}

void makehists::endJobs(){

    normTH1D(hist1d[15]);
    normTH1D(hist1d[16]);
    normTH1D(hist1d[39]);
    normTH1D(hist1d[40]);
    normTH1D(hist1d[161]);

    normTH1D(hist1d[48]);
    normTH1D(hist1d[49]);
    normTH1D(hist1d[41]);
    normTH1D(hist1d[42]);

    normTH1D(hist1d[27]);
    normTH1D(hist1d[28]);

    //normTH1D(hist1d[59]);
    //normTH1D(hist1d[60]);

//    hist2d[73]->Divide(hist2d[74]);
//    hist2d[75]->Divide(hist2d[76]);
//    hist2d[77]->Divide(hist2d[78]);

}

void makehists::initHists(){

	for( int it = 0; it < nHists; it++ ){ hist1d[it] = NULL; hist2d[it] = NULL; }
    //for( int it = 0; it < nHists; it++ ){ std::cout << " hist set check: " << hist1d[it] << std::endl; }

    //int jtdiv(400);
    //float jtran(8);
    int jtdiv(625);
    float jtran(25);
    int jdtdiv(200);
    float jdtran(4);
    int jztdiv(100);
    float jztran(2);
    int rhcnt(80);

    auto stddiv = 120;
    auto stdtran = 3;

    auto cldiv = 80;
    auto cltrn = 20;

    auto cl3ddiv = 200;
    auto cl3dtrn = 4;
    auto cl3ddiv1 = 200;
    auto cl3dtrn1 = 4;

    auto clsphdiv = 400;
    auto clsphtrn = 4;
    auto cwdiv = 80;
    auto cwtrn = 40;
    auto slmax = 120;
    auto slmin = -120;
    //auto sldiv = 120;
    auto sl3max = 12;
    auto sl3min = -12;
    auto sl3div = 2400;
    auto sldiv = 320;
    auto chimax = 1.01;
    auto chimin = 0.91;
    auto chidiv = 400;



    //------ 1D Hists --------------------------------------------------------------------------

    //hist1d[0] = new TH1D("jetRHTime", "jetRHTime", 2000, -100, 100);
    hist1d[1] = new TH1D("jetDRRHMulti", "jetDRRHMulti", rhcnt, 0, rhcnt);
    hist1d[2] = new TH1D("jetDRTimeError", "jetDRTimeError", 300, 0, 3);
    hist1d[3] = new TH1D("jetDRTimeRMS", "jetDRTimeRMS", 200, 0, 20);
    hist1d[4] = new TH1D("jetDRMedTime", "jetMedDRTime", jtdiv, -1*jtran, jtran);
    //hist1d[5] = new TH1D("jetRawTime", "jetRawTime", jtdiv, -1*jtran, jtran);
    hist1d[6] = new TH1D("jetCDRMuTime", "jetCDRMuTime", jtdiv, -1*jtran, jtran);
    hist1d[7] = new TH1D("jetCDRMedTime", "jetCDRMedTime", jtdiv, -1*jtran, jtran);
    hist1d[8] = new TH1D("jetCSCMuTime", "jetCSCMuTime", jtdiv, -1*jtran, jtran);
    hist1d[9] = new TH1D("jetCSCMedTime", "jetCSCMedTime", jtdiv, -1*jtran, jtran);
    hist1d[10] = new TH1D("jetCBCMedTime", "jetCBCMedTime", jtdiv, -1*jtran, jtran);
    //hist1d[11]
    hist1d[12] = new TH1D("jetPt", "jetPt", 500, 0, 5000);
    hist1d[13] = new TH1D("jetPhi", "jetPhi", 700, -3.5, 3.5);
    hist1d[14] = new TH1D("jetEta", "jetEta", 700, -3.5, 3.5);
    hist1d[15] = new TH1D("jetdtmu", "jetdtmu", jdtdiv, -1*jdtran, jdtran);
    hist1d[16] = new TH1D("jetdtmed", "jetdtmed", jdtdiv, -1*jdtran, jdtran);
    hist1d[17] = new TH1D("jetHt", "jetHt", 1000, 0, 5000);
    hist1d[18] = new TH1D("nJet", "nJets", 21, -0.5, 20.5);
    hist1d[19] = new TH1D("jetCBCMuTime", "jetCBCMuTime", jtdiv, -1*jtran, jtran);
    //hist1d[20]
    //hist1d[21]
    //hist1d[22]
    //hist1d[23]
    hist1d[24] = new TH1D("diffPt", "diffPt", 1000, 0, 10);
    hist1d[25] = new TH1D("htPct", "htPct", 100, 0, 1);
    hist1d[26] = new TH1D("dPhi", "dPhi", 70, -3.5, 3.5);

    hist1d[27] = new TH1D("jetcmudtbc", "jetcmudtbc", jdtdiv, -1*jdtran, jdtran);
    hist1d[28] = new TH1D("jetcmeddtbc", "jetcmeddtbc", jdtdiv, -1*jdtran, jdtran);
    hist1d[29] = new TH1D("jetMuTime", "jetMuTime", jtdiv, -1*jtran, jtran);
    hist1d[30] = new TH1D("nGoodDrJets", "nGoodDrJets", 21, -0.5, 20.5);
    hist1d[31] = new TH1D("nGoodScJets", "nGoodScJets", 21, -0.5, 20.5);
    hist1d[32] = new TH1D("nGoodBcJets", "nGoodBcJets", 21, -0.5, 20.5);
    hist1d[33] = new TH1D("nUnJets", "nUnJets", 101, -0.5, 100.5);
    hist1d[34] = new TH1D("pJets", "pJets", 110, 0, 1.1);
    hist1d[35] = new TH1D("pGoodDrJets", "pGoodDrJets", 110, 0, 1.1);
    hist1d[36] = new TH1D("pGoodScJets", "pGoodScJets", 110, 0, 1.1);
    hist1d[37] = new TH1D("pGoodBcJets", "pGoodBcJets", 110, 0, 1.1);
    hist1d[38] = new TH1D("pGoodBcToScJets", "pGoodBcToScJets", 110, 0, 1.1);

    hist1d[39] = new TH1D("jetcmudt", "jetcmudt", jdtdiv, -1*jdtran, jdtran);
    hist1d[40] = new TH1D("jetcmeddt", "jetcmeddt", jdtdiv, -1*jdtran, jdtran);
    hist1d[41] = new TH1D("jetcmudtsc", "jetcmudtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[42] = new TH1D("jetcmeddtsc", "jetcmeddtsc", jdtdiv, -1*jdtran, jdtran);

    //hist1d[43] = new TH1D("nPhotonsPerJet","nPhotonsPerJet", 21, -0.5, 20.5);

    hist1d[44] = new TH1D("jetSCmedTime", "jetSCmedTime", jtdiv, -1*jtran, jtran);
    hist1d[45] = new TH1D("jetSCmuTime", "jetSCmuTime", jtdiv, -1*jtran, jtran);
    //hist1d[46] = new TH1D("jetSCTimeRms", "jetSCTimeRms", 60, 0, 6);
    //hist1d[47] = new TH1D("jetSCrawTime", "jetSCrawTime", jtdiv, -1*jtran, jtran);

    hist1d[48] = new TH1D("jetmeddtsc", "jetmeddtsc", jdtdiv, -1*jdtran, jdtran);
    hist1d[49] = new TH1D("jetmudtsc", "jetmudtsc", jdtdiv, -1*jdtran, jdtran);

    //hist1d[50] = new TH1D("jetSCTimeSkew", "jetSCTimeSkew", 80, -4.0, 4.0);
    //hist1d[51] = new TH1D("jetPhotons", "jetPhotons", 21, -0.5, 20.5);
    //hist1d[52] = new TH1D("jetElectrons", "jetElectrons", 21, -0.5, 20.5);

    hist1d[53] = new TH1D("scbcdt", "scbcdt", jdtdiv, -1*jdtran, jdtran);
    //hist1d[54] = new TH1D("bc1rhef", "bc1rhef", 110, 0, 1.1);
    hist1d[55] = new TH1D("nBCinJet", "nBCinJet", 11, -0.5, 10.5);
    //hist1d[56] = new TH1D("bcMrhef", "bcMrhef", 110, 0, 1.1);

    //hist1d[57] = new TH1D("jetCPhMuTime", "jetCPhMuTime", jtdiv, -1*jtran, jtran);
    //hist1d[58] = new TH1D("jetCPhMedTime", "jetCPhMedTime", jtdiv, -1*jtran, jtran);
    //hist1d[59] = new TH1D("jetmudtph", "jetmudtph", jdtdiv, -1*jdtran, jdtran);
    //hist1d[60] = new TH1D("jetmudtel", "jetmudtel", jdtdiv, -1*jdtran, jdtran);
    //hist1d[61] = new TH1D("jetCEleMuTime", "jetCEleMuTime", jtdiv, -1*jtran, jtran);
    //hist1d[62] = new TH1D("jetCEleMedTime", "jetCEleMedTime", jtdiv, -1*jtran, jtran);

    hist1d[63] = new TH1D("scEtaPhiAngle3D", "scEtaPhiAngle3D", 660, -0.2, 6.4);
    hist1d[64] = new TH1D("scEtaTimAngle3D", "scEtaTimAngle3D", 660, -0.2, 6.4);
    hist1d[65] = new TH1D("scEtaPhiAngle2D", "scEtaTimAngle2D", 660, -0.2, 6.4);

    //hist1d[66] = new TH1D("sciEta3D", "sciEta", 171, -85, 85);
    //hist1d[67] = new TH1D("sciPhi3D", "sciPhi", 361, 0, 360);
    //hist1d[68] = new TH1D("scTime3D", "scTime", 5000, -25, 25);
    //hist1d[69] = new TH1D("sciEta2D", "sciEtaDiff", 201, -100, 100 );
    //hist1d[70] = new TH1D("sciPhi2D", "sciPhiDiff", 201, -100, 100);
    //hist1d[71] = new TH1D("sciTim2D", "sciTimDiff", 400, -10, 10);
    //hist1d[72] = new TH1D("sciAngle2D", "sciAngleSph", 660, -0.2, 6.4);
    //hist1d[73] = new TH1D("scRotAngle", "scRotAngle", 660, -0.2, 6.4);
    //hist1d[74] = new TH1D("scAngleTest", "scAngleTest", 660, -0.2, 6.4);
    //hist1d[75] = new TH1D("sciEta3Diff", "sciEta3Diff", 400, -10, 10 );
    //hist1d[76] = new TH1D("sciPhi3Diff", "sciPhi3Diff", 400, -10, 10);
    //hist1d[77] = new TH1D("sciTim3Diff", "sciTim3Diff", 400, -10, 10);
    //hist1d[78] = new TH1D("scSphEgn0", "scSphEgn0", 400, -10, 10);
    //hist1d[79] = new TH1D("scSphEgn1", "scSphEgn1", 400, -10, 10);
    //hist1d[80] = new TH1D("scSinTest", "scSinTest", 200, -1, 1);
    hist1d[81] = new TH1D("eginValueSph", "eginValueSph", 150, 0.4, 1.1);
    //hist1d[82] = new TH1D("dIPhiTest", "dIPhiTest", 1440, -720, 720);
    //hist1d[83] = new TH1D("meanIPhiTest", "meanIPhiTest", 1440, -720, 720);
    //hist1d[84] = new TH1D("meanIEtaTest", "meanIEtaTest", 200, -100, 100);
    //hist1d[85] = new TH1D("meanTimeTest", "meanTimeTest", 2000, -25, 25);
    hist1d[86] = new TH1D("eginValue3D", "eginValue3D", 110, 0, 1.1);

    //hist1d[87] = new TH1D("cluster_etprofile", "Cluster Eta Time Profile Sph", cwdiv, -1*cwtrn, cwtrn);
    //hist1d[88] = new TH1D("cluster_et3Dprofl", "Cluster Eta Time Profile 3D", cl3ddiv1, -1*cl3dtrn1, cl3dtrn1);

    hist1d[89] = new TH1D("jetmudtgen", "jetmudtgen", jdtdiv, -1*jdtran, jdtran);
    hist1d[90] = new TH1D("genJetImpactAngle", "genJetImpactAngle", 660, -0.2, 6.4);

    //hist1d[91] = new TH1D("clEtaTimeSlope", "Cluster Eta Time Slope", 240, -120, 120);
    //hist1d[92] = new TH1D("clEtaTimeSlopeChi", "Cluster Eta Time Slope Chi2", 100, 0, 1);

    hist1d[93] = new TH1D("jetGenTimeCut", "jetGenTimeCut", jtdiv, -1*jtran, jtran);

    //hist1d[94] = new TH1D("clBinProfile_m1", "spCluster Bin -1 Profile", 200, -10, 10);
    //hist1d[95] = new TH1D("clBinProfile_p1", "spCluster Bin +1 Profile", 200, -10, 10);
    //hist1d[96] = new TH1D("clBinProfile_m2", "spCluster Bin -2 Profile", 200, -10, 10);
    //hist1d[97] = new TH1D("clBinProfile_p2", "spCluster Bin +2 Profile", 200, -10, 10);
    //hist1d[98] = new TH1D("clBinProfile_m3", "spCluster Bin -3 Profile", 200, -10, 10);
    //hist1d[99] = new TH1D("clBinProfile_p3", "spCluster Bin +3 Profile", 200, -10, 10);

    //hist1d[100] = new TH1D("clETSlopes", "Cluster ET Slopes", 240, -120, 120);

    //hist1d[101]
    //hist1d[102]

    //hist1d[103] = new TH1D("clEtaTimeSlopeInv", "Cluster Eta Time SlopeInv", 350, -0.1, 34.9);
    //hist1d[104] = new TH1D("clEtaTimeSlopeInv3D", "Cluster Eta Time SlopeInv 3D", 350, -0.1, 34.9);

    hist1d[105] = new TH1D("genJetDrMatchJet", "genJetDrMatchJet", 100, 0, 0.1);
    hist1d[106] = new TH1D("genJetSCTimeDiff", "genJetSCTimeDiff", 300, 0, 30);
    hist1d[107] = new TH1D("genJetDrTimeDiff", "genJetSCTimeDiff", 300, 0, 30);

    hist1d[108] = new TH1D("jetGenTimeVar", "jetGenTimeVar", 408, -2, 100);
    hist1d[109] = new TH1D("jetGenTimeNextBX", "jetGenTimeNextBX", 3, -1, 2);
    hist1d[110] = new TH1D("jetGenTime", "jetGenTime", jtdiv, -1*jtran, jtran);
    hist1d[111] = new TH1D("jetGenTOF", "jetGenTOF", 300, 0, 30);
    hist1d[112] = new TH1D("jetGenTimeIsLLP", "jetGenTimeIsLLP", 3, -1, 2);
    hist1d[113] = new TH1D("jetGenTimeLLPPurity", "jetGenTimeLLPPurity", 100, 0, 1);

    //hist1d[114] = new TH1D("clETSlope3D", "Cluster EtaTimeSlope 3D", 500, -250, 250);
    //hist1d[115] = new TH1D("clETSlopeChi3D", "Cluster EtaTimeSlope Chi2 3D", 100, 0, 1);

    hist1d[116] = new TH1D("jetGenNKids", "jetGenNKids", 100, 0, 100);
    hist1d[117] = new TH1D("clEtaTimeSlopeRangeA", "Cluster Eta Time Slope Eta > 1.0", 480, -240, 240);
    hist1d[118] = new TH1D("clEtaTimeSlopeRangeB", "Cluster Eta Time Slope Eta < 0.5", 480, -240, 240);

    //hist1d[119] = new TH1D("ootPhotonTime", "ootPhotonTime", jtdiv, -1*jtran, jtran);

    //hist1d[120] = new TH1D("jetCOOTPhMuTime", "jetCOOTPhMuTime", jtdiv, -1*jtran, jtran);
    //hist1d[121] = new TH1D("jetCOOTPhMedTime", "jetCOOTPhMedTime", jtdiv, -1*jtran, jtran);

    //hist1d[122] = new TH1D("jetPhClRhTime", "phClRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[123] = new TH1D("jetPhClRhPkOOT", "phClRhPkOOT", 120, -0.1, 1.1);
    //hist1d[124] = new TH1D("jetPhClRhPMatched", "phClRhPMatched", 120, -0.1, 1.1);
    //hist1d[125] = new TH1D("jetOOTPhClRhTime", "ootPhClRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[126] = new TH1D("jetOOTPhClRhPkOOT", "ootPhClRhPkOOT", 120, -0.1, 1.1);
    //hist1d[127] = new TH1D("jetOOTPhClRhPMatched", "ootPhClRhPMatched", 120, -0.1, 1.1);
    //hist1d[128] = new TH1D("jetEleClRhTime", "eleClRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[129] = new TH1D("jetEleClRhPkOOT", "eleClRhPkOOT", 120, -0.1, 1.1);
    //hist1d[130] = new TH1D("jetEleClRhPMatched", "eleClRhPMatched", 120, -0.1, 1.1);

    hist1d[131] = new TH1D("ebRhTime", "ebRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[132] = new TH1D("ebRhEnergy", "ebRhEnergy", 1000, 0, 1000);
    hist1d[133] = new TH1D("eeRhTime", "eeRhTime", jtdiv*2, -1*jtran*2, jtran*2);
    hist1d[134] = new TH1D("eeRhEnergy", "eeRhEnergy", 1000, 0, 1000);
    hist1d[135] = new TH1D("ebRhkOOT", "ebRhkOOT", 3, 0, 1);
    hist1d[136] = new TH1D("eeRhkOOT", "eeRhkOOT", 3, 0, 1);

    //hist1d[137] = new TH1D("phClRhTime", "phClRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[138] = new TH1D("phClRhPkOOT", "phClRhPkOOT", 120, -0.1, 1.1);
    //hist1d[139] = new TH1D("phClRhPMatched", "phClRhPMatched", 120, -0.1, 1.1);
    //hist1d[140] = new TH1D("ootPhClRhTime", "ootPhClRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[141] = new TH1D("ootPhClRhPkOOT", "ootPhClRhPkOOT", 120, -0.1, 1.1);
    //hist1d[142] = new TH1D("ootPhClRhPMatched", "ootPhClRhPMatched", 120, -0.1, 1.1);
    //hist1d[143] = new TH1D("eleClRhTime", "eleClRhTime", jtdiv, -1*jtran, jtran);
    //hist1d[144] = new TH1D("eleClRhPkOOT", "eleClRhPkOOT", 120, -0.1, 1.1);
    //hist1d[145] = new TH1D("eleClRhPMatched", "eleClRhPMatched", 120, -0.1, 1.1);

    hist1d[146] = new TH1D("cljTime", "cljTime", jtdiv, -1*jtran, jtran);
    hist1d[147] = new TH1D("cljLeadRhTime", "cljLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[148] = new TH1D("cljClstLeadTimeDiff", "cljClstLeadTimeDiff", jtdiv, -1*jtran, jtran);
    hist1d[149] = new TH1D("ootPhoTime", "ootPhoTime", jtdiv, -1*jtran, jtran);
	hist1d[150] = new TH1D("ootPhoSeedRhTime", "ootPhoLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[151] = new TH1D("ootPhoSeedTimeDiff", "ootPhoLeadTimeDiff", jtdiv, -1*jtran, jtran);
    hist1d[152] = new TH1D("eleClTime", "eleClTime", jtdiv, -1*jtran, jtran);
    hist1d[153] = new TH1D("eleSeedRhTime", "eleLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[154] = new TH1D("eleClSeedTimeDiff", "eleClLeadTimeDiff", jtdiv, -1*jtran, jtran);
    hist1d[155] = new TH1D("phoClTime", "phoTime", jtdiv, -1*jtran, jtran);
    hist1d[156] = new TH1D("phoSeedRhTime", "phoLeadRhTime", jtdiv, -1*jtran, jtran);
    hist1d[157] = new TH1D("phoSeedTimeDiff", "phoLeadTimeDiff", jtdiv, -1*jtran, jtran);

    hist1d[158] = new TH1D("cljDrTime", "cljDrTime", jtdiv, -1*jtran, jtran);
    hist1d[159] = new TH1D("cljDiffPt", "cljDiffPt", 1000, 0, 10);
    hist1d[160] = new TH1D("cljdPhi", "cljdPhi", 70, -3.5, 3.5);
    hist1d[161] = new TH1D("cljdtmu", "cljdtmu", jdtdiv, -1*jdtran, jdtran);

    for( int it = 0; it < nHists; it++ ){ if(hist1d[it]) hist1d[it]->Sumw2();}

    //------ 2D Hists --------------------------------------------------------------------------

    //hist2d[0] = new TH2D("jt", "jt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[1] = new TH2D("jt_pt", "jt_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[2] = new TH2D("jt_id", "jt_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
    hist2d[3] = new TH2D("jt_nhf", "jt_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[4] = new TH2D("jt_chf", "jt_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[5] = new TH2D("jt_nemf", "jt_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[6] = new TH2D("jt_cemf", "jt_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[7] = new TH2D("jt_muf", "jt_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[8] = new TH2D("jt_nhm", "jt_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
    hist2d[9] = new TH2D("jt_chm", "jt_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

    hist2d[10] = new TH2D("jt_medt", "jt_medt", jtdiv, -1*jtran, jtran, 200, -10, 10);
    hist2d[11] = new TH2D("jt_rms", "jt_rms", jtdiv, -1*jtran, jtran, 200, 0, 20);
    hist2d[12] = new TH2D("jt_err", "jt_err", jtdiv, -1*jtran, jtran, 300, 0, 3);

    hist2d[13] = new TH2D("medt_pt", "medt_pt", jtdiv, -1*jtran, jtran, 500, 0, 500);
    hist2d[14] = new TH2D("medt_id", "medt_id", jtdiv, -1*jtran, jtran, 5, 0, 5);
    hist2d[15] = new TH2D("medt_nhf", "medt_nhf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[16] = new TH2D("medt_chf", "medt_chf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[17] = new TH2D("medt_nemf", "medt_nemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[18] = new TH2D("medt_cemf", "medt_cemf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[19] = new TH2D("medt_muf", "medt_muf", jtdiv, -1*jtran, jtran, 100, 0, 1);
    hist2d[20] = new TH2D("medt_nhm", "medt_nhm", jtdiv, -1*jtran, jtran, 40, 0, 40);
    hist2d[21] = new TH2D("medt_chm", "medt_chm", jtdiv, -1*jtran, jtran, 40, 0, 40);

    hist2d[22] = new TH2D("jdtmu_nJets", "jdtmu_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);
    hist2d[23] = new TH2D("jdtmed_nJets", "jdtmed_nJets", jdtdiv, -1*jdtran, jdtran, 6, 2, 8);

    hist2d[24] = new TH2D("jt_nrh", "jt_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);
    hist2d[25] = new TH2D("medt_nrh", "medt_nrh", jtdiv, -1*jtran, jtran, 50, 0, 50);

    hist2d[26] = new TH2D("jdtmu_diffPt", "jdtmu_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[27] = new TH2D("jdtmu_htPct", "jdtmu_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[28] = new TH2D("jdtmu_dPhi", "jdtmu_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);
    hist2d[29] = new TH2D("jdtmed_diffPt", "jdtmed_diffPt", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[30] = new TH2D("jdtmed_htPct", "jdtmed_htPct", jdtdiv, -1*jdtran, jdtran, 200, 0.8, 1);
    hist2d[31] = new TH2D("jdtmed_dPhi", "jdtmed_dPhi", jdtdiv, -1*jdtran, jdtran, 400, 2.8, 3.2);

    hist2d[32] = new TH2D("jt_sceta", "jt_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    hist2d[33] = new TH2D("jt_scphi","jt_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    hist2d[34] = new TH2D("jt_scenr", "jt_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    hist2d[35] = new TH2D("medt_sceta", "medt_sceta", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    hist2d[36] = new TH2D("medt_scphi","medt_scphi", jtdiv, -1*jtran, jtran, 700, -3.5, 3.5);
    hist2d[37] = new TH2D("medt_scenr", "medt_scenr", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    hist2d[38] = new TH2D("rht_rhe", "rht_rhe", jtdiv, -1*jtran, jtran, 1000, 0, 1000);

    hist2d[39] = new TH2D("njrh_nscrh", "njrh_nscrh", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);
    //hist2d[40] = new TH2D("je_sce", "je_sce", 500, 0, 500, 500, 0, 500);
    //hist2d[41] = new TH2D("je_ege", "je_ege", 500, 0, 500, 500, 0, 500);
    //hist2d[42] = new TH2D("ege_sce", "ege_sce", 500, 0, 500, 500, 0, 500);
    hist2d[43] = new TH2D("scdt_effje", "scdt_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);
    hist2d[44] = new TH2D("jdt_effje", "jdt_effje", jdtdiv, -1*jdtran, jdtran, 250, 0, 500);

    //hist2d[45]
    //hist2d[46]
    //hist2d[47]
    //hist2d[48]
    //hist2d[49]
    //hist2d[50]

    hist2d[51] = new TH2D("njrh_nsbcrh", "njrh_nsbcrh", rhcnt, 0, rhcnt, rhcnt, 0, rhcnt);

    //hist2d[52] = new TH2D("dremf_emf", "dremf_emf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[53] = new TH2D("scemf_emf", "scemf_emf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[54] = new TH2D("bcemf_emf", "bcemf_emf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[55] = new TH2D("dreme_eme", "dreme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[56] = new TH2D("sceme_eme", "sceme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[57] = new TH2D("bceme_eme", "bceme_eme", 500, 0, 500, 500, 0, 500);
    hist2d[58] = new TH2D("sce_bce", "sce_bce", 500, 0, 500, 500, 0, 500);
    hist2d[59] = new TH2D("dremf_scemf", "dremf_scemf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[60] = new TH2D("scemf_bcemf", "scemf_bcemf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[61] = new TH2D("epaf_epf", "epaf_epf", 110, 0, 1.1, 110, 0, 1.1);
    hist2d[62] = new TH2D("epf_emf", "epf_emf", 110, 0, 1.1, 110, 0, 1.1);

    //hist2d[63] = new TH2D("jetEta_stdSCt", "jetEta_stdSCt", 200, -2.5, 2.5, stddiv, 0, stdtran);
    //hist2d[64] = new TH2D("jetEtaetaMmt_stdSCt", "jetEtaetaMmt_stdSCt", 800, 0, 0.02, stddiv, 0, stdtran);
    //hist2d[65] = new TH2D("jetPhiphiMnt_stdSCt", "jetPhiphiMnt_stdSCt", 800, 0, 0.02, stddiv, 0, stdtran);
    //hist2d[66] = new TH2D("jetEtaphiMnt_stdSCt", "jetEtaphiMnt_stdSCt", 1600, -0.02, 0.02, stddiv, 0, stdtran);
    //hist2d[67] = new TH2D("jetMaxD_stdSCt", "jetMaxD_stdSCt", 400, 0, 1, stddiv, 0, stdtran);
    //hist2d[68] = new TH2D("jetConPtDis_stdSCt", "jetConPtDis_stdSCt", 400, 0, 1, stddiv, 0, stdtran);
    //hist2d[69] = new TH2D("jetConEtaPhiSprd_stdSCt", "jetConEtaPhiSprd_stdSCt", 600, -0.005, 0.01, stddiv, 0, stdtran);
    //hist2d[70] = new TH2D("jetArea_stdSCt", "jetArea_stdSCt", 200, 0, 1, stddiv, 0, stdtran);
    //hist2d[71] = new TH2D("jetNCarry_stdSCt", "jetNCarry_stdSCt", 40, 0, 40, stddiv, 0, stdtran);
    //hist2d[72] = new TH2D("jetNConst_stdSCt", "jetNConst_stdSCt", 100, 0, 100, stddiv, 0, stdtran);

    //hist2d[73]
    //hist2d[74]
    //hist2d[75]
    //hist2d[76]
    //hist2d[77]
    //hist2d[78]
    //hist2d[79]
    //hist2d[80]

    hist2d[81] = new TH2D("scSphEgn01", "scSphEgn01", 200, -5, 5, 200, -5, 5);
    hist2d[82] = new TH2D("sc3DEgn01", "sc3DEgn01", 200, -5, 5, 200, -5, 5);
    hist2d[83] = new TH2D("sc3DEgn02", "sc3DEgn02", 200, -5, 5, 200, -5, 5);
    //hist2d[84] = new TH2D("cluster_3D_etwtmap", "Cluster Eta(x)Time(y) WtMap 3D;eta;time",cl3ddiv1,-1*cl3dtrn1,cl3dtrn1,cl3ddiv,-1*cl3dtrn,cl3dtrn);
    //hist2d[85] = new TH2D("cluster_3D_etwtprof","Clstr Eta(x)Time(y) WtMap 3DPrfle;eta;time",cl3ddiv1,-1*cl3dtrn1,cl3dtrn1,cl3ddiv,-1*cl3dtrn,cl3dtrn);

    //hist2d[86] = new TH2D("cluster_ptwtmap", "Cluster Phi(x)Time(y) Wt Map Sph;Phi;Time", cwdiv, -1*cwtrn, cwtrn, clsphdiv, -1*clsphtrn, clsphtrn);
    //hist2d[87] = new TH2D("cluster_etwtmap", "Cluster Eta(x)Time(y) Wt Map Sph;Eta;Time", cwdiv, -1*cwtrn, cwtrn, clsphdiv, -1*clsphtrn, clsphtrn);

    hist2d[89] = new TH2D("jetEtavSlope", "Jet Eta v Slope;Eta;Slope ps/cm", 60, -1.5, 1.5, sldiv, slmin, slmax);
    hist2d[90] = new TH2D("jetImpAnglevSlope", "Jet ImpactAngle v Slope;ImpactAngle;Slope ps/cm", 60, -1.5, 1.5, sldiv, slmin, slmax);
    //hist2d[91] = new TH2D("clEtaTimeSlvChi", "Cluster EtaTime Slope v Chi2;Slope;Chi2", sldiv, slmin, slmax, chidiv, chimin, chimax);
    //hist2d[92] = new TH2D("clEtaTimeSlvEVal", "Cluster EtaTime Slope v EigenValue;Slope;EigenValue", sldiv, slmin, slmax, 180, 0.55, 1.0);
    //hist2d[93] = new TH2D("clEtaTimeSlvRotAng", "Cluster EtaTime Slope v Rotation Angle;Slope;rotAngle", sldiv, slmin, slmax, 660, -0.2, 6.4);
    //hist2d[94] = new TH2D("clEtaTimeSlvNumClRHs", "Cluster EtaTime Slope v nClRecHits;Slope;nClRecHits", sldiv, slmin, slmax, 60, 0, 60);
    //hist2d[95] = new TH2D("clEtaTimeChi2vEVal", "Cluster EtaTime Chi2 v EigenValue;Chi2;EigenValue", chidiv, chimin, chimax, 60, 0.45, 1.05);
    hist2d[96] = new TH2D("jetImpAnglevSlope3D", "Jet ImpactAngle v Slope 3D;ImpactAngle;Slope", 150, 0, 1.5, sl3div, sl3min, sl3max);
    //hist2d[97] = new TH2D("clEtaTimeChi2vNumClRHs", "Cluster EtaTime Chi2 v nClRecHits;Chi2;nClRecHits", chidiv, chimin, chimax, 60, 0, 60);

    hist2d[98] = new TH2D("jetEvGenE", "Jet Energy v GenEnergy;JetEnergy;GenEnergy", 100, 0, 1000, 100, 0, 1000 );
    hist2d[99] = new TH2D("jetEGenERatiovGenTime", "Jet E/GenE v GenTime;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[100] = new TH2D("jetEMFracvSCTime", "Jet EMFrac v SCTime;EMFrac;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[101] = new TH2D("jetGenRatiovGenTimePre", "Jet E/GenE v GenTime Pre;E/GenE;GenTime", 80, 0, 2, 40, -15, 25 );
    hist2d[102] = new TH2D("jetGenTimevDrJetTime", "GenTime v DrJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[103] = new TH2D("jetEGenERatiovSCTimeDiff", "Jet E/GenE v JetSC GenJet TimeDif;E/GenE;SCTimeDif", 80, 0, 2, 300, 0, 30.0 );
    hist2d[104] = new TH2D("jetSCTimevDrTime", "JetSCTime v JetDrTime;JetSCTime;JetDrTime", 280, -15, 25, 280, -15, 25 );
    hist2d[105] = new TH2D("jetGenTimevSCJetTime", "GenTime v SCJetTime;GenTime;JetTime", 280, -15, 25, 280, -15, 25 );
    hist2d[106] = new TH2D("jetGenTimevGenEnergy", "GenTime v GenEnergy;GenTime;GenEnergy", 280, -15, 25, 100, 0, 1000 );
    hist2d[107] = new TH2D("jetGenjetDrvSCTimeDiff", "Jet GenJet Dr v SCJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[108] = new TH2D("jetGenjetDrvDRTimeDiff", "Jet GenJet Dr v DRJet GenJet TimeDif;jetGenDr;TimeDif", 200, 0, 0.5, 300, 0, 30.0 );
    hist2d[109] = new TH2D("jetGenTimevTOFcorr", "GenTime v TOFcorr;GenTime;TOFcorr", 250, 0, 25, 250, 0, 25 );
    hist2d[110] = new TH2D("jetEGenERatiovSCTime", "Jet E/GenE v SCTime;E/GenE;SCTime", 80, 0, 2, 40, -15, 25 );
    hist2d[111] = new TH2D("jetEGenERatiovDRTime", "Jet E/GenE v DRTime;E/GenE;DRTime", 80, 0, 2, 40, -15, 25 );
    hist2d[112] = new TH2D("jetGenVarvSCTimeDiff", "Jet GenJet Var v SCJet GenJet TimeDif;Var;TimeDif", 270, -2, 25, 200, 0, 20.0 );
    hist2d[113] = new TH2D("jetGenPurityvSCTimeDiff", "Jet GenJet Purity v SCJet GenJet TimeDif;Purity;TimeDif", 100, 0, 1, 200, 0, 20.0 );
    hist2d[114] = new TH2D("jetGenPurityvGenJetVar", "Jet GenJet Purity v GenJet Var;Purity;Var", 100, 0, 1, 270, -2, 25 );
    hist2d[115] = new TH2D("jetGenVarvGenJetNKids", "Jet GenJet Var v GenJet nKids;Var;nKids", 270, -2, 25, 100, 0, 100 );
    hist2d[116] = new TH2D("jetGenPurityvGenJetNKids", "Jet GenJet Purity v GenJet nKids;Purity;nKids", 100, 0, 1, 100, 0, 100 );
    hist2d[117] = new TH2D("genJetSCTimeDiffvDrMatchJet", "genJet SCTimeDiff v DrMatchJet;SCTimeDiff;DrMatchJet", 300, 0, 30, 320, 0, 3.2 );
    hist2d[118] = new TH2D("jetGenTimevJetEMFrac", "GenTime v JetEMFrac;GenTime;JetEMFrac", 280, -15, 25, 150, 0, 1.5 );

    //hist2d[119] = new TH2D("jetSlopevRotSlope", "Jet Slope v rotated Slope;Slope;rotated Slope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    //hist2d[120] = new TH2D("jetSlopevDifSlope", "Jet Slope v dif w/ rotSlope;Slope;difSlope", sldiv, slmin, slmax, sldiv, slmin, slmax);
    hist2d[121] = new TH2D("jetPhivSlope", "Jet Phi v Slope;Phi;Slope ps/cm", 140, -3.5, 3.5, sldiv, slmin, slmax);

    hist2d[122] = new TH2D("ebRhTimevEnergy", "ebRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );
    hist2d[123] = new TH2D("eeRhTimevEnergy", "eeRhTimevEnergy;Time [ns];Energy [GeV]", jtdiv*2, -1*jtran*2, jtran*2, 1000, 0, 1000 );

    hist2d[124] = new TH2D("jetEmFracvGenjetDr", "Jet emFrac v GenJet Dr;emFrac;jetGenDr", 40, 0, 1, 200, 0, 0.5 );
    hist2d[125] = new TH2D("jetEGenERatiovGenjetDr", "Jet E/GenE v GenJet Dr;E/GenE;jetGenDr", 80, 0, 2, 200, 0, 0.5 );

    hist2d[126] = new TH2D("jetEtavSlope3D", "Jet Eta v Slope 3D;Eta;Slope ps/cm", 60, -1.5, 1.5, sl3div, sl3min, sl3max);
    hist2d[127] = new TH2D("jetPhivSlope3D", "Jet Phi v Slope 3D;Phi;Slope ps/cm", 140, -3.5, 3.5, sl3div, sl3min, sl3max);

}//<<>>void makehists::initHists()


int main ( int argc, char *argv[] ){

    if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    else {
                auto indir = argv[1];
                auto infilename = argv[2];
                auto outfilename = argv[3];
				makehists base;				
                base.llpgana_hist_maker( indir, infilename, outfilename );
    }
    return 1;
}



