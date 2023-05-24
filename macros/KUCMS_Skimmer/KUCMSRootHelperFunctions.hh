// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

#include "KUCMSHelperFunctions.hh"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
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
#include "TF1.h"
#include "TFormula.h"
#include "Math/PositionVector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TGraph.h"
#include "TMathBase.h"

#ifndef KUCMSRootHelperFunctionsH
#define KUCMSRootHelperFunctionsH 

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
    for (auto ibinX = 1; ibinX <= nXBins; ibinX++){

        auto phist = (TH1F*)nhist->ProjectionY("temp",ibinX,ibinX);
        auto mean = phist->GetMean();
        auto stdv = phist->GetStdDev();
        auto norm = phist->GetBinContent(phist->GetMaximumBin());
        auto high = mean + 0.2*stdv;
        auto low = mean - 0.2*stdv;
        if( abs(stdv) > 0.01 && abs(norm) > 1 ){
            auto tmp_form = new TFormula("tmp_formula","[0]*exp(-0.5*((x-[1])/[2])**2)");
            auto tmp_fit  = new TF1("tmp_fit",tmp_form->GetName(),low,high);
            tmp_fit->SetParameter(0,norm); //tmp_fit->SetParLimits(0,norm/2,norm*2);
            tmp_fit->SetParameter(1,mean); //tmp_fit->SetParLimits(1,-2,2);
            tmp_fit->SetParameter(2,stdv); //tmp_fit->SetParLimits(2,0,10);
            phist->Fit(tmp_fit->GetName(),"RBQ0");
            auto fmean = tmp_fit->GetParameter(1);
            auto error = tmp_fit->GetParError(1);
            auto fNdf = tmp_fit->GetNDF();
            auto fProb = tmp_fit->GetProb();
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

std::vector<float> getDistStats( std::vector<float> values ){

    std::vector<float> results;
    int size = values.size();

    if( size == 0 ){
        std::cout << "value group is empty" << std::endl;
        results.push_back(-29.25);//0
        results.push_back(99.9);//1
        results.push_back(-29.25);//2
        results.push_back(99.9);//3
        results.push_back(99.9);//4
        results.push_back(9.9);//5
        results.push_back(-29.25);//6
        results.push_back(99.9);//7
        results.push_back(99.9);//8
        results.push_back(99.9);//9
        results.push_back(-29.25);//10
        results.push_back(99.9);//11
        return results;
    }//<<>>if( size == 0 )

    TH1D valueDist( "temp", "temp", 500, -25, 25 );
    for( const auto value : values ) valueDist.Fill(value);
    auto nEntries = valueDist.GetEntries();
    double medvalue(0.0);
    if( nEntries > 0 ){
        double quant(0.5);
        valueDist.ComputeIntegral();
        valueDist.GetQuantiles(1, &medvalue, &quant);
    } else medvalue = -9.9;


    auto error = valueDist.GetMeanError();
    results.push_back(valueDist.GetMean());//0
    results.push_back(error);//1
    results.push_back(medvalue);//2
    results.push_back(1.2533*error);//3 valid for noraml distibutions
    results.push_back(valueDist.GetRMS());//4
    results.push_back(valueDist.GetSkewness());//5

    auto rmu = mean(values);
    results.push_back(rmu);//6 mean
    auto rsd = stdev(values, rmu);
    results.push_back(rsd);//7 stdev
    results.push_back(rms(values));//8 rms
    auto err = rsd/std::sqrt(size);
    results.push_back(err);//9 error of mean

    std::sort(values.begin(),values.end());
    if( size%2 == 0 ) results.push_back((values[(size/2)-1] + values[(size/2)])/2);//10
    else results.push_back(values[size/2]);//10 median
    results.push_back(1.2533*err);//11 error of median

    return results;

}//>>>> std::vector<float> getDistStats( std::vector<float> values )

std::vector<float> getDistStats( std::vector<float> values, std::vector<float> wts ){

    std::vector<float> results;
    int size = values.size();
    int wtsize = wts.size();
    if( size != wtsize || size == 0 ){
        std::cout << "value & rechit groups not same or empty" << std::endl;
        results.push_back(-29.25);//0 mu
        results.push_back(99.9);//1
        results.push_back(-29.25);//2
        results.push_back(99.9);//3
        results.push_back(99.9);//4
        results.push_back(9.9);//5
        results.push_back(-29.25);//6
        results.push_back(99.9);//7
        results.push_back(99.9);//8
        results.push_back(99.9);//9
        results.push_back(-29.25);//10
        results.push_back(9.9);//11
        results.push_back(-10.0);//12
        return results;
    }//<<>>if( size != rechits.size() || size == 0 )

    TH1D valueDist( "temp", "temp", 500, -10, 10 );
    std::vector<float> wvalues;
    float wtot(0.0);
    for( int it(0); it < size; it++ ){
        valueDist.Fill(values[it], wts[it]);
        wvalues.push_back(values[it]*wts[it]);
        wtot += wts[it];
    }//<<>>for( uInt it(0); it < size; it++ )
    //auto nEntries = valueDist.GetEntries();
    double hmedvalue(9.9);
    //if( nEntries > 0 ){
    //  double quant(0.5);
    //  valueDist.ComputeIntegral();
    //  valueDist.GetQuantiles(1, &hmedvalue, &quant);
    //} else hmedvalue = -9.9; 

    auto herror = valueDist.GetMeanError();
    auto hmuvalue = valueDist.GetMean();
    auto hrms = valueDist.GetRMS();
    results.push_back(hmuvalue);//0
    results.push_back(herror);//1
    results.push_back(hmedvalue);//2
    results.push_back(1.2533*herror);//3 - valid for noraml distibutions
    results.push_back(hrms);//4
    results.push_back(valueDist.GetSkewness());//5

    auto rmu = mean(wvalues,wtot);
    results.push_back(rmu);//6 mean
    auto rsd = stdev(values, rmu, wts, wtot);
    results.push_back(rsd);//7 stdev
    results.push_back(rms(wvalues));//8 rms
    auto err = rsd/std::sqrt(size);
    results.push_back(err);//9 error of mean

    std::sort(values.begin(),values.end());
    if( size%2 == 0 ) results.push_back((values[(size/2)-1]+values[(size/2)])/2);//10
    else results.push_back(values[size/2]);//10 median
    results.push_back(1.2533*err);//11 error of median  
    results.push_back(wtot);//12 tot e of all rechits
    //results[2] = results[6] - results[10];

    return results;

}//>>>>std::vector<float> getDistStats( std::vector<float> values, std::vector<float> weights )

std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> wts ){
//spherical

    std::vector<float> results;

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

}//<<>>std::vector<float> getRhGrpEigen2D( std::vector<float> xs, std::vector<float> ys, std::vector<float> zs, std::vector<float> wts )

std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> ys, std::vector<float> wts ){

    std::vector<float> results;

    auto mean_x = mean( xs, wts );
    auto mean_y = mean( ys, wts );
    auto swts = accum( wts );
    auto var_x = var( xs, mean_x, wts, swts );
    auto var_y = var( ys, mean_y, wts, swts );
    auto var_xy = cvar( xs, mean_x, ys, mean_y, wts, swts );
    //std::cout << "Varencies and Means Calculated" << std::endl;

    //TMatrixDSym mat(3,3);
    double array[] = { var_x, var_xy, var_xy, var_y };
    TMatrixDSym mat(2,array);
    //std::cout << "Input matrix created" << std::endl;
    //mat.Print();

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

}//<<>>std::vector<float> = getRhGrpEigen2D( std::vector<float> xs, std::vector<float> ys, std::vector<float> zs, std::vector<float> wts )

std::vector<float> getRhGrpEigen( std::vector<float> xs, std::vector<float> ys, std::vector<float> zs, std::vector<float> wts ){
// ieipt

    std::vector<float> results;

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

    int index(2);
    int index2(2);
    if( (eginVal(0) >= eginVal(1)) && (eginVal(0) >= eginVal(2)) ){
        index = 0;
        if( eginVal(1) >= eginVal(2) ) index2 = 1;
        if( (eginVal(0) == eginVal(1)) && (eginVal(0) == eginVal(2)) ) std::cout << " -- rhGrp is Spherical" << std::endl;
    } else if( eginVal(1) >= eginVal(2) ) {
        index = 1;
        if( eginVal(0) >= eginVal(2) ) index2 = 0;
        if( eginVal(1) == eginVal(2) ) std::cout << " -- rhGrp is a Flatend Sphere" << std::endl;
    } else { if( eginVal(0) >= eginVal(1) ) index2 = 0; else index2 = 1; } //<<>>if( (eginVal(0) >= eginVal(1)) && (eginVal(0) >= eginVal(2)) )
    //std::cout << "Largest Eigin Value Found" << std::endl;

    auto index3 = 3 - index - index2;
    auto ev23hypo = hypo(eginVal[index2], eginVal[index3]);
    auto ev1Shypo = hypo(eginVal[index], ev23hypo);
    auto speginval = sq2(eginVal[index]/ev1Shypo);

    auto ex = eginVec(index,0);
    auto ey = eginVec(index,1);
    auto ez = eginVec(index,2);
    auto ev = speginval; //eginVal(index);
    //std::cout << "Eigin Angles Found" << std::endl;

    results.push_back(ex);//0
    results.push_back(ey);//1
    results.push_back(ez);//2
    results.push_back(ev);//3

    return results;

}//<<>>std::vector<float> getRhGrpEigen3D( std::vector<float> xs, std::vector<float> ys, std::vector<float> zs, std::vector<float> wts )

#endif
//----------------------------------------------------------------------------------------------------------------------



