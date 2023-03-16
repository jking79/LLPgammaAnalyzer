// ROOT inludes
#include "TStyle.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
//#include "TString.h"
#include "TColor.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
#include "TMath.h"
#include "TROOT.h"
#include "TPaveText.h"

// STL includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <cstdlib>
#include <utility>
#include <sys/stat.h>

// Common include(s)
//#include "Common.cpp"  --  still used for SetUpBins :(
//#include "CommonTimeFit.cpp"
#include "wc_ku_gausTimeFit.cpp"

// set bins helper function

std::string RemoveDelim(std::string tmp, const std::string & delim){return tmp.erase(tmp.find(delim),delim.length());}

void setBins(std::string & str, std::vector<Double_t> & bins, Bool_t & var_bins){

    if(str.find("CONSTANT") != std::string::npos){

        var_bins = false;
        str = RemoveDelim(str,"CONSTANT");
        Int_t nbins = 0; Double_t low = 0.f, high = 0.f;
        std::stringstream ss(str);
        ss >> nbins >> low >> high;
        Double_t bin_width = (high-low)/nbins;
        //std::cout << "Setting Const bins : ";
        for (Int_t ibin = 0; ibin <= nbins; ibin++){
            //std::cout << low+ibin*bin_width << " ";
            bins.push_back(low+ibin*bin_width);
        }//<<>.for (Int_t ibin = 0; ibin <= nbins; ibin++)
        //std::cout << std::endl;

    } else if(str.find("VARIABLE") != std::string::npos) {

        var_bins = true;
        str = RemoveDelim(str,"VARIABLE");
        Float_t bin_edge;
        std::stringstream ss(str);
        //std::cout << "Setting Var bins : ";
        while(ss >> bin_edge){
            //std::cout << bin_edge << " ";
            bins.push_back(bin_edge);
        }//<<>>while (ss >> bin_edge)
        //std::cout << std::endl;

    } else {

        std::cerr << "Aye... bins are either VARIABLE or CONSTANT! Exiting..." << std::endl;
        exit(1);

    }//<<>>if      (str.find("CONSTANT") != std::string::npos)

}//<<>>void setBins(std::string & str, std::vector<Double_t> & bins, Bool_t & var_bins)

// fit struct
struct FitStruct{
 
  FitStruct() {}
  FitStruct(const std::string & label, const std::string & inHistName) : label(label), inHistName(inHistName) {}
  
  const std::string label;
  const std::string inHistName;
  
  TH2D * Hist2D;
  std::map<Int_t,TimeFitStruct*> TimeFitStructMap;
  std::map<std::string,TH1F*> ResultsMap;
  
  TFormula * SigmaForm;
  TF1 * SigmaFit;

};//<<>>struct FitStruct

// sigma fit params
struct SigmaFitParams {
 
  SigmaFitParams() {}
  
  Float_t low;
  Float_t val;
  Float_t up;

};//<<>>struct SigmaFitParams

std::vector<Double_t> makeConstantBins( int nbins, float low, float high ){

	double bin_width = (high-low)/nbins;
	std::vector<Double_t> result;
	for (Int_t ibin = 0; ibin <= nbins; ibin++){ result.push_back(low+ibin*bin_width); }
	return result;

}//<<>>std::vector<Double_t> makeConstantBins

void SetupTimeFitType(const std::string & str, TimeFitType & type){

  if      (str.find("Gaus1")       != std::string::npos) type = TimeFitType::Gaus1;
  else if (str.find("Gaus1coretail")   != std::string::npos) type = TimeFitType::Gaus1coretail;
  else if (str.find("Gaus1core")   != std::string::npos) type = TimeFitType::Gaus1core;
  else if (str.find("Gaus2fm")     != std::string::npos) type = TimeFitType::Gaus2fm;
  else if (str.find("Gaus2fmcore") != std::string::npos) type = TimeFitType::Gaus2fmcore;
  else if (str.find("Gaus3fm")     != std::string::npos) type = TimeFitType::Gaus3fm;
  else if (str.find("Gaus3fmcore") != std::string::npos) type = TimeFitType::Gaus3fmcore;
  else { std::cerr << "Specified a non-supported fit type: " << str.c_str() << " ... Exiting..." << std::endl; exit(1);}

}//<<>>void SetupTimeFitType(const std::string & str, TimeFitType & type)


void DumpFitInfo(FitStruct & DataInfo, const std::string & outfiletext, std::vector<Double_t> fXBins ){

  std::cout << "Dumping fit info into text file..." << std::endl;
  // get histograms!
  const auto & data_mu_hist = DataInfo.ResultsMap["mu"];
  const auto & data_sigma_hist = DataInfo.ResultsMap["sigma"];

  // make dumpfile object 
  const std::string filename = outfiletext+"_FitInfo.log"; //Common::outFitText.Data()+"."+Common::outTextExt.Data();
  std::ofstream dumpfile(filename.c_str(),std::ios_base::out);
  auto & TimeFitStructMap = DataInfo.TimeFitStructMap;

  dumpfile << std::setw(5)  << "Bin |"
           << std::setw(18) << "  pT range          |"
           << std::setw(18) << "  Data mu           |"
           << std::setw(18) << "  Data sigma        |"
           << std::setw(18) << "  Data n            |"
           << std::setw(18) << "  Data n2           |"
           << std::setw(18) << "  Data sigma2       |"
           << std::setw(22) << "  Data chisqr/ndf   "
           << std::endl;

  std::string space = "";
  const auto nw = 158;// = 5+18+19+19+19+19+17+38+38;
  for (auto i = 0; i < nw; i++) space += "-";

  dumpfile << space.c_str() << std::endl;

  auto fNBinsX = fXBins.size() - 1;
  for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++){

    std::cout << "Dumping Info for " << ibinX << " of " << fNBinsX << " bins. " << std::endl;
    auto & TimeFit = TimeFitStructMap[ibinX];
    const auto pt_low = fXBins[ibinX-1];
    const auto pt_up  = fXBins[ibinX];

    const auto data_mu   = data_mu_hist->GetBinContent(ibinX);
    const auto data_mu_e = data_mu_hist->GetBinError  (ibinX);
    const auto data_sigma   = data_sigma_hist->GetBinContent(ibinX);
    const auto data_sigma_e = data_sigma_hist->GetBinError  (ibinX);

    auto param0 = 0.f;
    auto paramErr0 = 0.f;
    auto param1 = 0.f;
    auto paramErr1 = 0.f;
    auto param2 = 0.f;
    auto paramErr2 = 0.f;
    auto param3 = 0.f;
    auto paramErr3 = 0.f;
    auto param4 = 0.f;
    auto paramErr4 = 0.f;
    auto chisqr = 0.f;
    auto ndf = 0;

    if (not TimeFit->isEmpty()){
    std::cout << "Number of Parameters " << TimeFit->fit->GetNpar() << std::endl;
         if( (TimeFit->fit)->GetNpar() > 0 ) {
              param0 = TimeFit->fit->GetParameter(0);
              paramErr0 = TimeFit->fit->GetParError (0);
              param1 = TimeFit->fit->GetParameter(1);
              paramErr1 = TimeFit->fit->GetParError (1);
              param2 = TimeFit->fit->GetParameter(2);
              paramErr2 = TimeFit->fit->GetParError (2);
              chisqr = TimeFit->fit->GetChisquare();
              ndf = TimeFit->fit->GetNDF();
			}//<<>>if( (TimeFit->fit)->GetNpar() > 0 )
     		if( (TimeFit->fit)->GetNpar() > 3 ) {
              param3 = TimeFit->fit->GetParameter(3);
              paramErr3 = TimeFit->fit->GetParError (3);
              param4 = TimeFit->fit->GetParameter(4);
              paramErr4 = TimeFit->fit->GetParError (4);
			}//<<>>if( (TimeFit->fit)->GetNpar() > 3 )
    } else { std::cout << "Time Fit Empty " << std::endl; }//<<>>if (not TimeFit->isEmpty())
    dumpfile << std::setw(5)  << Form("%i |",ibinX)
        		 << std::setw(18) << Form(" %6.1f - %6.1f |",pt_low,pt_up)
             << std::setw(20) << Form(" %4.4f +/- %3.4f |",data_mu,data_mu_e)
             << std::setw(20) << Form(" %4.4f +/- %3.4f |",data_sigma,data_sigma_e)
             << std::setw(20) << Form(" %4.4f +/- %3.4f |",param0,paramErr0)
             << std::setw(20) << Form(" %4.4f +/- %3.4f |",param3,paramErr3)
             << std::setw(20) << Form(" %4.4f +/- %3.4f |",param4,paramErr4)
             << std::setw(22) << Form(" %6.4f/%i = %4.4f ",chisqr,ndf,chisqr/ndf)
             << std::endl;
    if (ibinX % 20 == 0) dumpfile << space.c_str() << std::endl;
  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

  // get fits!
  std::cout << "Do Sigma Fit Info Dump"<< std::endl;
  const auto & data_sigma_fit = DataInfo.SigmaFit;
  dumpfile << space.c_str() << std::endl << std::endl;
  dumpfile << "Sigma Fit Results" << std::endl;
  dumpfile << std::setw(12)  << "Par |" << std::setw(23) << "         Data         |" << std::endl;
  // loop over params and dump
  for (auto ipar = 0; ipar < data_sigma_fit->GetNpar(); ipar++){
    // get constants
    const auto data_sigma_fit_par   = data_sigma_fit->GetParameter(ipar);
    const auto data_sigma_fit_par_e = data_sigma_fit->GetParError (ipar);

    dumpfile << std::setw(12)  << Form("%s |",data_sigma_fit->GetParName(ipar)) 
				 << std::setw(18) << Form(" %8.4f +/- %7.4f |",data_sigma_fit_par,data_sigma_fit_par_e) << std::endl;

  }//<<>>for (auto ipar = 0; ipar < data_sigma_fit->GetNpar(); ipar++)
  const auto data_sigma_fit_chi   = data_sigma_fit->GetChisquare();
  const auto data_sigma_fit_ndf   = data_sigma_fit->GetNDF();
  dumpfile << std::setw(12)  << "ChiSqr/NDF |"  << std::setw(23)  
           << Form(" %8.4f/%i = %4.4f ", data_sigma_fit_chi, data_sigma_fit_ndf, data_sigma_fit_chi/data_sigma_fit_ndf ) << std::endl;

}//<<>>void DumpFitInfo(FitStruct & DataInfo, const std::string & outfiletext, std::vector<Double_t> fXBins )

//************** Primary function *************************************************************************************
void runTimeFitter(const std::string & infilename, const std::string & plotconfig, const std::string & miscconfig,
		   const std::string & timefitconfig, const std::string & era, const std::string & outfiletext, const std::string & inplotname, std::string xbinstr ){

//----header------------------------------------------

  // settings
  const std::string fInFileName(infilename);
  const std::string fPlotConfig(plotconfig);
  const std::string fMiscConfig(miscconfig);
  const std::string fTimeFitConfig(timefitconfig);
  const std::string fEra(era);
  const std::string fOutFileText(outfiletext);

  //bool doEstTOF = true;
  bool doEstTOF = false;
  float estTOF = 0.01;

  // style
  TStyle * fTDRStyle;
  std::string fTitle;
  std::string fXTitle;
  //std::vector<Double_t> fXBins;
  Int_t fNBinsX(0);
  //Bool_t fXVarBins;

  // var fit config
  TimeFitType fTimeFitType;
  Float_t fRangeLow;
  Float_t fRangeUp;
  std::string fTimeText;
  std::string f2DHistName;
  //std::string f2DGHistName;

  // sigma fit config
  std::string fSigmaVarText;
  std::string fSigmaVarUnit;
  SigmaFitParams fSigmaInitN;
  SigmaFitParams fSigmaInitC;

//--constructor----------------------------------------------------------------------
  std::cout << "Initializing TimeFitter..." << std::endl;
  // input
  TFile * fInFile;

  // output
  TFile * fOutFile;
  TPaveText * fConfigPave;

  fInFile = TFile::Open(infilename.c_str());
  fTDRStyle = new TStyle("TDRStyle","Style for P-TDR");
  std::string fOutFileName = fOutFileText + ".root";
  fOutFile = TFile::Open(fOutFileName.c_str(),"UPDATE");

  //std::cout << "Setting up Common..." << std::endl;
  //Common::SetupEras();
  //Common::SetupSamples();
  //Common::SetupGroups();
  //Common::SetupHistNames();

  fTitle = "#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)";
  fXTitle = "A_{eff}/#sigma_{n} (EBEB)";
  //std::string xbinstr  = "VARIABLE 0 75 100 125 150 175 225 275 325 375 475 600 950 2250";
  std::vector<Double_t> fXBins;
  Bool_t fXVarBins = false;//dummy not used
  setBins(xbinstr,fXBins,fXVarBins);
  //std::vector<Double_t> fXBins{75.0,100.0,125.0,150.0,175.0,225.0,275.0,325.0,375.0,475.0,600.0,950.0,2250.0};
  //std::vector<Double_t> fXBins{100.0,125.0,150.0,175.0,225.0,275.0,325.0,375.0,475.0,600.0,950.0,2250.0};
  //std::vector<Double_t> fXBins{75.0,150.0,175.0,225.0,275.0,325.0,375.0,475.0,600.0,950.0,2250.0,4500.0};
  //std::vector<Double_t> fXBins{75.0,150.0,175.0,225.0,275.0,325.0,375.0,475.0,600.0,950.0,2250.0,9000.0};
  //std::vector<Double_t> fXBins{0.0,75.0,125.0,175.0,225.0,275.0,325.0,375.0,425.0,500.0,700.0,2250.0,9000.0};
  fNBinsX = fXBins.size();
  std::string ftstr = "Gaus1core";
  SetupTimeFitType(ftstr,fTimeFitType);
  fRangeLow = 2.0f;
  fRangeUp = 2.0f;
  fTimeText = "t_{1}-t_{2}";
  fSigmaVarText = "A_{eff}/#sigma_{n}";
  fSigmaVarUnit = "";
  //TimeFitter::ReadInitParams("0 50 100",fSigmaInitN);
  fSigmaInitN.low = 0; fSigmaInitN.val = 50; fSigmaInitN.up = 100;
  //TimeFitter::ReadInitParams("0 .5 1",fSigmaInitC);
  fSigmaInitC.low = 0; fSigmaInitC.val = .5; fSigmaInitC.up = 1;
  f2DHistName=inplotname;
  //f2DGHistName="SRO_Data_GHist";

  // set fitter
  //TVirtualFitter::SetDefaultFitter("Minuit2");

//---------fitter.MakeTimeFits()------ run script substances class with <TimeFitter fitter;> then runs <fitter.MakeTimeFits();> -------------------
  std::cout << "Making time fits..." << std::endl;
  // Do data first
  FitStruct FitInfo("Data","Data"); //Common::HistNameMap["Data"].Data());
  //FitStruct GFitInfo("GData","GData");
  //MakeTimeFit(DataInfo>>FitInfo);
  //......................................
  const auto & label = FitInfo.label;
  //const auto & glabel = GFitInfo.label;
  std::cout << "Making time fits for: " << label << std::endl;

  // Get input hist
  //------------TimeFitter::GetInputHist(FitInfo);
  //const auto & label = FitInfo.label;
  std::cout << "Getting input hist: " << label << std::endl;
  // get input
  const auto & inHistName = FitInfo.inHistName;
  //const auto & inGHistName = GFitInfo.inHistName;
  auto & Hist2D = FitInfo.Hist2D;
  //auto & GHist2D = GFitInfo.Hist2D;
  // get the hist
  Hist2D = (TH2D*)fInFile->Get(f2DHistName.c_str());
  //GHist2D = (TH2D*)fInFile->Get(f2DGHistName.c_str());
  // save to output
  fOutFile->cd();
  //Hist2D->Write(Hist2D->GetName(),TObject::kWriteDelete);
  //GHist2D->Write(GHist2D->GetName(),TObject::kWriteDelete);

  // Init time fits
  //-------------TimeFitter::InitTimeFits(FitInfo);
  //const auto & label = FitInfo.label;
  std::cout << "Initializing TimeFitStructMap for: " << label << std::endl;
  // get inputs/outputs
  auto & TimeFitStructMap = FitInfo.TimeFitStructMap;
  //auto & GTimeFitStructMap = GFitInfo.TimeFitStructMap;
  // setup a time fit struct for each bin!
  for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++){ TimeFitStructMap[ibinX] = new TimeFitStruct(fTimeFitType,fRangeLow,fRangeUp);}
  //for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++){ GTimeFitStructMap[ibinX] = new TimeFitStruct(fTimeFitType,fRangeLow,fRangeUp);}

  // Project out 2D hists into map
  //-------------TimeFitter::Project2Dto1DHists(FitInfo);
  //const auto & label = FitInfo.label;
  std::cout << "Projecting to 1D from 2D plot: " << label << std::endl;
  // get inputs/outputs
  //const auto & Hist2D = FitInfo.Hist2D;
  //auto & TimeFitStructMap = FitInfo.TimeFitStructMap;
  const std::string histname1 = Hist2D->GetName();
  //const std::string ghistname1 = GHist2D->GetName();
  for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++){
    auto & hist = TimeFitStructMap[ibinX]->hist;
    //auto & ghist = GTimeFitStructMap[ibinX]->hist;
	std::string fullHistName = histname1 + "_ibin" + std::to_string(ibinX);
    //std::string fullGHistName = ghistname1 + "_ibin" + std::to_string(ibinX);
    hist = (TH1F*)Hist2D->ProjectionY(fullHistName.c_str(),ibinX,ibinX);
    //ghist = (TH1F*)GHist2D->ProjectionY(fullGHistName.c_str(),ibinX,ibinX);
	hist->SetTitle(fullHistName.c_str());
    //ghist->SetTitle(fullGHistName.c_str());
  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

  // Fit each 1D hist
  //--------------TimeFitter::Fit1DHists(FitInfo);
  //const auto & label = FitInfo.label;
  std::cout << "Fitting hists for: " << label << std::endl;
  // get inputs/outputs
  //auto & TimeFitStructMap = FitInfo.TimeFitStructMap;
  for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++){
    // get pair input
    auto & TimeFit = TimeFitStructMap[ibinX];
    //auto & GTimeFit = GTimeFitStructMap[ibinX];
    // get hist, and skip if no entries
    if (TimeFit->isEmpty()) continue;
    // Prep the fit
    TimeFit->PrepFit();
    //GTimeFit->PrepFit();
    // do the fit!
    TimeFit->DoFit();
    //GTimeFit->DoFit();
    // save output
    fOutFile->cd();
    TimeFit->hist->Write(TimeFit->hist->GetName(),TObject::kWriteDelete);
    //GTimeFit->hist->Write(GTimeFit->hist->GetName(),TObject::kWriteDelete);
    TimeFit->form->Write(TimeFit->form->GetName(),TObject::kWriteDelete);
    //GTimeFit->form->Write(GTimeFit->form->GetName(),TObject::kWriteDelete);
    TimeFit->fit->Write(TimeFit->fit->GetName(),TObject::kWriteDelete);
    //GTimeFit->fit->Write(GTimeFit->fit->GetName(),TObject::kWriteDelete);
  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

  // Extract mu and sigma into maps
  //---------------TimeFitter::ExtractFitResults(FitInfo);
  //const auto & label = FitInfo.label;
  std::cout << "Extracting results for: " << label << std::endl;
  // get inputs/outputs
  //auto & TimeFitStructMap = FitInfo.TimeFitStructMap;
  auto & ResultsMap = FitInfo.ResultsMap;
  //auto & GResultsMap = GFitInfo.ResultsMap;
  // setup hists
  const auto xbins = &fXBins[0];
  //TH1F * TimeFitter::SetupHist(const TString & ytitle, const TString & yextra, const TString & label) 
  //ResultsMap["chi2ndf"]  = TimeFitter::SetupHist("#chi^{2}/NDF","chi2ndf",label);
  auto chi2ndfName = inplotname+"_chi2ndf";
  //auto gchi2ndfName = glabel+"_chi2ndf";
  auto chi2ndfTitle = fTitle+" "+"#chi^{2}/NDF."+";"+fXTitle+";"+"#chi^{2}/NDF.";
  ResultsMap["chi2ndf"]  = new TH1F(chi2ndfName.c_str(),chi2ndfName.c_str(),fNBinsX,xbins);
  ResultsMap["chi2ndf"]->Sumw2();
  //GResultsMap["chi2ndf"]  = new TH1F(gchi2ndfName.c_str(),chi2ndfName.c_str(),fNBinsX,xbins);
  //GResultsMap["chi2ndf"]->Sumw2();
  //ResultsMap["chi2prob"] = TimeFitter::SetupHist("#chi^{2} Prob.","chi2prob",label);
  auto chi2probName = inplotname+"_chi2prob";
  //auto gchi2probName = glabel+"_chi2prob";
  auto chi2probTitle = fTitle+" "+"#chi^{2} Prob."+";"+fXTitle+";"+"#chi^{2} Prob.";
  ResultsMap["chi2prob"] = new TH1F(chi2probName.c_str(),chi2probTitle.c_str(),fNBinsX,xbins);
  ResultsMap["chi2prob"]->Sumw2();
  //GResultsMap["chi2prob"] = new TH1F(gchi2probName.c_str(),chi2probTitle.c_str(),fNBinsX,xbins);
  //GResultsMap["chi2prob"]->Sumw2();
  auto mutext = "#mu("+fTimeText+") [ns]";
  auto muName = inplotname+"_mu";
  //auto gmuName = glabel+"_mu";
  auto muTitle = fTitle+" "+mutext+";"+fXTitle+";"+mutext;
  //ResultsMap["mu"]       = TimeFitter::SetupHist(mutext,"mu",label);
  ResultsMap["mu"] = new TH1F(muName.c_str(),muTitle.c_str(),fNBinsX,xbins);
  ResultsMap["mu"]->Sumw2();
  //GResultsMap["mu"] = new TH1F(gmuName.c_str(),muTitle.c_str(),fNBinsX,xbins);
  //GResultsMap["mu"]->Sumw2();
  auto sigtext = "#sigma("+fTimeText+") [ns]";
  auto sigName = inplotname+"_sigma";
  //auto gsigName = glabel+"_sigma";
  auto sigTitle = fTitle+" "+sigtext+";"+fXTitle+";"+sigtext;
  //ResultsMap["sigma"]    = TimeFitter::SetupHist(sigtext,"sigma",label);
  ResultsMap["sigma"] = new TH1F(sigName.c_str(),sigTitle.c_str(),fNBinsX,xbins);
  ResultsMap["sigma"]->Sumw2();
  //GResultsMap["sigma"] = new TH1F(gsigName.c_str(),sigTitle.c_str(),fNBinsX,xbins);
  //GResultsMap["sigma"]->Sumw2();
  auto occtext = "#occ("+fTimeText+")";
  auto occName = inplotname+"_occ";
  //auto goccName = glabel+"_occ";
  auto occTitle = fTitle+" "+occtext+";"+fXTitle+";"+occtext;
  //ResultsMap["occ"]       = TimeFitter::SetupHist(occtext,"occ",label);
  ResultsMap["occ"] = new TH1F(occName.c_str(),occTitle.c_str(),fNBinsX,xbins);
  ResultsMap["occ"]->Sumw2();
  //GResultsMap["occ"] = new TH1F(goccName.c_str(),occTitle.c_str(),fNBinsX,xbins);
  //GResultsMap["occ"]->Sumw2();
  auto rmstext = "#rms("+fTimeText+")";
  auto rmsName = inplotname+"_rms";
  //auto grmsName = glabel+"_rms";
  auto rmsTitle = fTitle+" "+rmstext+";"+fXTitle+";"+rmstext;
  //ResultsMap["rms"]       = TimeFitter::SetupHist(occtext,"occ",label);
  ResultsMap["rms"] = new TH1F(rmsName.c_str(),rmsTitle.c_str(),fNBinsX,xbins);
  ResultsMap["rms"]->Sumw2();
  //GResultsMap["rms"] = new TH1F(grmsName.c_str(),rmsTitle.c_str(),fNBinsX,xbins);
  //GResultsMap["rms"]->Sumw2();

  // set bin content!
  for (auto ibinX = 1; ibinX < fNBinsX; ibinX++){ 
    // get time fit
    auto & TimeFit = TimeFitStructMap[ibinX];
    //auto & GTimeFit = GTimeFitStructMap[ibinX];
    // skip if fit not present
    if (TimeFit->isEmpty()) continue;
    // get results
    TimeFit->GetFitResult();
    //GTimeFit->GetFitResult();
    //const auto & gresult = GTimeFit->result;
    const auto & result = TimeFit->result;
	auto errHistName = "bin" + std::to_string(ibinX) + "ErrHist";
    auto errHist = (TH2F*)fInFile->Get(errHistName.c_str());
	auto merr = errHist->GetMean();
	auto nerr = errHist->Integral();
    auto sgerr = errHist->GetStdDev();
    auto err = sgerr/sqrt(nerr);
    //auto gerr = (0.5f)*(gresult.sigma - result.sigma);
    //auto gerr = (gresult.sigma - result.sigma);
	 //auto gserr = sqrt(gerr*gerr+result.esigma*result.esigma+gresult.esigma*gresult.esigma);
    //std::cout << "sigma: " << result.sigma ;
    //std::cout << "Sigma ( err hist: " << err << " sig: " << sgerr << " n: " << nerr << " m: " << merr << " )";
    std::cout << " Bin " << ibinX  <<  " : " << fXBins[ibinX-1] << "-" << fXBins[ibinX]; 
	std::cout << " : no guass: " << result.sigma-estTOF << " err: " << result.esigma << std::endl; 
	//std::cout << " :: w/ guass:" << gresult.sigma << " err:" << gresult.esigma << std::endl;
	//auto fserr = gresult.esigma * 8;
	//if( ibinX == 1 ) fserr = gresult.esigma * 12;
    //if( ibinX == 2 ) fserr *= 3;
    if( not doEstTOF ) estTOF = 0.0f; 
    //std::cout << " :: w/ guass:" << gresult.sigma-estTOF << " err:" << fserr << std::endl;
    //auto fserr = merr;
    // set bin content
    ResultsMap["chi2ndf"] ->SetBinContent(ibinX,result.chi2ndf);
    ResultsMap["chi2prob"]->SetBinContent(ibinX,result.chi2prob);
    //GResultsMap["chi2ndf"] ->SetBinContent(ibinX,gresult.chi2ndf);
    //GResultsMap["chi2prob"]->SetBinContent(ibinX,gresult.chi2prob);
    ResultsMap["mu"]      ->SetBinContent(ibinX,result.mu);
    ResultsMap["mu"]      ->SetBinError  (ibinX,result.emu);
    //GResultsMap["mu"]      ->SetBinContent(ibinX,gresult.mu);
    //GResultsMap["mu"]      ->SetBinError  (ibinX,gresult.emu);
    ResultsMap["sigma"]   ->SetBinContent(ibinX,result.sigma-estTOF);
    ResultsMap["sigma"]   ->SetBinError  (ibinX,result.esigma);
    //GResultsMap["sigma"]   ->SetBinContent(ibinX,gresult.sigma-estTOF);
    //GResultsMap["sigma"]   ->SetBinError  (ibinX,fserr);
    ResultsMap["occ"]   ->SetBinContent(ibinX,result.occ);
    ResultsMap["occ"]   ->SetBinError  (ibinX,sqrt(result.occ));
    //GResultsMap["occ"]   ->SetBinContent(ibinX,gresult.occ);
    //GResultsMap["occ"]   ->SetBinError  (ibinX,sqrt(gresult.occ));
    ResultsMap["rms"]   ->SetBinContent(ibinX,result.rms);
    ResultsMap["rms"]   ->SetBinError  (ibinX,err);
    //GResultsMap["rms"]   ->SetBinContent(ibinX,gresult.rms);
    //GResultsMap["rms"]   ->SetBinError  (ibinX,err);

  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)
  // save output
  fOutFile->cd();
  for (const auto & ResultsPair : ResultsMap) ResultsPair.second->Write(ResultsPair.second->GetName(),TObject::kWriteDelete);
  //for (const auto & ResultsPair : GResultsMap) ResultsPair.second->Write(ResultsPair.second->GetName(),TObject::kWriteDelete);

  //..........................................
  //MakeSigmaFit(DataInfo>>FitInfo);
  // Dump mu's and sigma's into text file
  //const auto & label = FitInfo.label;
  const auto & sflabel = FitInfo.label;
  //std::cout << "Making sigma fit for: " << sflabel << std::endl;

  // Prep sigma fit
  //----------------------TimeFitter::PrepSigmaFit(FitInfo);
  std::cout << "Prepping sigma fit for: " << sflabel << std::endl;
  // get input hist
  const auto & hist = FitInfo.ResultsMap["sigma"];
  // get range
  const auto x_low = hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetFirst());
  const auto x_up  = hist->GetXaxis()->GetBinUpEdge (hist->GetXaxis()->GetLast());
  // set names
  std::cout << "Prepping sigma form for: " << label << std::endl;
  const std::string histname = hist->GetName();
  const auto formname = histname+"_form";
  const auto fitname  = histname+"_fit";
  // get and set formula
  auto & form = FitInfo.SigmaForm;
  auto fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))";
  //auto fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1])+([2]*[2]/x))";
  form = new TFormula(formname.c_str(),fitformstr);
  // get and set fit
  std::cout << "Prepping sigma fit params for: " << label << std::endl;
  auto & fit = FitInfo.SigmaFit;
  fit = new TF1(fitname.c_str(),form->GetName(),x_low,x_up);
  // init params
  fit->SetParName(0,"N"); fit->SetParameter(0,fSigmaInitN.val); fit->SetParLimits(0,fSigmaInitN.low,fSigmaInitN.up);
  fit->SetParName(1,"C"); fit->SetParameter(1,fSigmaInitC.val); fit->SetParLimits(1,fSigmaInitC.low,fSigmaInitC.up);
  //fit->SetParName(2,"S"); fit->SetParameter(2,1.0); fit->SetParLimits(2,0.0,10.0);
  // set line color
  fit->SetLineColor(hist->GetLineColor());
  // save form
  std::cout << "Saving sigma form for: " << sflabel << std::endl;
  fOutFile->cd();
  form->Write(form->GetName(),TObject::kWriteDelete);

  // Fit sigma hist
  //---------------------TimeFitter::FitSigmaHist(FitInfo);
  //const auto & label = FitInfo.label;
  std::cout << "Fitting sigma hist for: " << label << std::endl;
  // get inputs/outputs
  //auto & hist = FitInfo.ResultsMap["sigma"];
  //auto & fit  = FitInfo.SigmaFit;
  // and fit it
  //auto c1 = new TCanvas( "c1", "canvas" , 200, 10, 700, 500 );
  //c1->SetBatch(true);
  //hist->Draw();
  //hist->Fit(fit->GetName(),"RBQ0");
  hist->Fit(fit->GetName(),"RBQM");
  //c1->Modified();
  //c1->Update();
  //c1->Draw();
  // save to output
  fOutFile->cd();
  fit->Write(fit->GetName(),TObject::kWriteDelete);
  //c1->Print( "test.png" );
  //.......................................
  std::string fOutFileLogText = fOutFileText + "_" + inplotname;
  DumpFitInfo(FitInfo, fOutFileLogText, fXBins);

  //.......................................
  //DeleteInfo(FitInfo);
  // Delete infos
  //const auto & label = FitInfo.label;
  std::cout << "Deleting info for: " << label << std::endl;
  delete FitInfo.SigmaFit;
  delete FitInfo.SigmaForm;
  //delete GFitInfo.SigmaFit;
  //delete GFitInfo.SigmaForm;
  //Common::DeleteMap(FitInfo.ResultsMap);
  for (auto & Pair : FitInfo.ResultsMap) delete Pair.second;
  //for (auto & Pair : GFitInfo.ResultsMap) delete Pair.second;
  FitInfo.ResultsMap.clear();
  //GFitInfo.ResultsMap.clear();
  // loop over good bins to delete things
  //auto & TimeFitStructMap = FitInfo.TimeFitStructMap;
  for (auto ibinX = 1; ibinX < fNBinsX; ibinX++){
    // get pair input
    auto & TimeFit = TimeFitStructMap[ibinX];
    //auto & GTimeFit = GTimeFitStructMap[ibinX];
    // delete internal members
    TimeFit->DeleteInternal();
    //GTimeFit->DeleteInternal();
    // finally, delete the object itself
    delete TimeFit;
    //delete GTimeFit;
  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)
  delete FitInfo.Hist2D;
  //delete GFitInfo.Hist2D;

//-----destructor--------------------------------------------------
  std::cout << "Tidying up in destructor..." << std::endl;
  //delete fConfigPave;
  delete fOutFile;
  delete fTDRStyle;
  delete fInFile;

}//<<>>void runTimeFitter(......

/*
//k********************************  the Main function *******************************************
int main ( int argc, char *argv[] ){

        //if( argc != 7 ) { std::cout << "Insufficent arguments." << std::endl; }
        //else {

               std::string plotconfig("");
               std::string miscconfig("");
               std::string timefitconfig("");
               std::string era("");

               std::string infilename = "gres_Run2022A_2dhists.root";
               std::string outfiletext = "test";
               std::string inplotname = "SRO_Data_Hist";

               runTimeFitter( infilename, plotconfig, miscconfig, timefitconfig, era, outfiletext, inplotname );

        //}
        return 1;

}//<<>>int main ( int argc, char *argv[] )
*/
