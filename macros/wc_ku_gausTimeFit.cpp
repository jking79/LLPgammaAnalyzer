// ROOT inludes
#include "TFormula.h"
#include "TF1.h"
#include "TH1F.h"

// STL includes
#include <iostream>
#include <cmath>
#include <string>

// Common include
//#include "Common.cpp"

// Time fit enum
enum TimeFitType {Gaus1, Gaus1core, Gaus1coretail, Gaus2fm, Gaus2fmcore, Gaus3fm, Gaus3fmcore};

// fit result struct
struct TimeFitResult {

//  TimeFitResult() {}

  Float_t chi2ndf;
  Float_t chi2prob;
  Float_t occ;
  Float_t rms;
  Float_t std;
  Float_t mu;
  Float_t emu;
  Float_t sigma;
  Float_t esigma;

};//<<>>struct TimeFitResult

class TimeFitStruct {

  public:

  TimeFitStruct() {}
  TimeFitStruct(const TimeFitType type, const Float_t rangeLow, const Float_t rangeUp) : type(type), rangeLow(rangeLow), rangeUp(rangeUp) {}

  // helper functions for making fits to variables
  Bool_t isEmpty() const {return (hist->GetEntries() == 0);}
  void PrepFit();
  void DoFit();
  void GetFitResult();
  void SetupTimeFitType(const std::string & str, TimeFitType & type);

  // cleanup
  void DeleteInternal();

  // internal data members
  TimeFitType type;
  Float_t rangeLow;
  Float_t rangeUp;
  TH1F * hist;
  Bool_t varBinsX;
  TFormula * form;
  TF1 * fit;
  TimeFitResult result;

};//<<>>struct TimeFitStruct

void TimeFitStruct::PrepFit() {

  constexpr float SqrtPI   = 1.77245385090551602729f;
  constexpr float Sqrt2PI   = 2.506628274631001f;
  std::string ftypes[] = {"Gaus1", "Gaus1core", "Gaus1coretail", "Gaus2fm", "Gaus2fmcore", "Gaus3fm", "Gaus3fmcore"};
  //std::cout << "PrepFit for " << ftypes[type] << std::endl;  

  // Word on fit notation
  // "GausN" == N Gaussians fit
  // "fm" == "fixed mean", i.e. for N Gaussian fit, all Gaussians share the same mu
  // "core" == mid point of range of fit is mean of the histogram, range is n times the std. dev of hist

  // set tmp init vals
  Float_t hsum  = hist->Integral(varBinsX?"width":"");
  Float_t mu    = hist->GetMean();
  Float_t sigma = hist->GetStdDev(); 
  result.rms = hist->GetRMS();  
  result.occ = hsum;
  result.std = sigma;
  auto norm = hsum/(sigma*Sqrt2PI);

  // range vars
  auto rangelow = 0.f;
  auto rangeup  = 0.f;

  // make tmp fit first if not gausNcore, set range
  if (type == TimeFitType::Gaus1 || type == TimeFitType::Gaus2fm || type == TimeFitType::Gaus3fm){

    // set range for tmp and main fit
    rangelow = rangeLow;
    rangeup  = rangeUp;

    auto tmp_form = new TFormula("tmp_formula","[0]*exp(-0.5*((x-[1])/[2])**2)");
    auto tmp_fit  = new TF1("tmp_fit",tmp_form->GetName(),rangelow,rangeup);

    tmp_fit->SetParameter(0,norm); //tmp_fit->SetParLimits(0,amp/2,norm*2);
    tmp_fit->SetParameter(1,mu); //tmp_fit->SetParLimits(1,-2,2);
    tmp_fit->SetParameter(2,sigma); tmp_fit->SetParLimits(2,0,10);

    // fit hist with tmp tf1
    hist->Fit(tmp_fit->GetName(),"RBQ0");

    norm  = tmp_fit->GetParameter(0); // constant
    mu    = tmp_fit->GetParameter(1); // mu
    sigma = tmp_fit->GetParameter(2); // sigma

    delete tmp_form;
    delete tmp_fit;

  } else { // "core" fits <<>>if (type == TimeFitType::Gaus1 || type == TimeFitType::Gaus2fm || type == TimeFitType::Gaus3fm)
    // set range for main fit
    rangelow = (mu-(rangeLow*sigma));
    rangeup  = (mu+(rangeUp *sigma));
  }//<<>>if (type == TimeFitType::Gaus1 || type == TimeFitType::Gaus2fm || type == TimeFitType::Gaus3fm) else
  
  // names for fits and formulas
  const std::string histname = hist->GetName();
  const auto formname = histname+"_formula";
  const auto fitname  = histname+"_fit";

  //--------------------------------------------------------------------
  if (type == TimeFitType::Gaus1 || type == TimeFitType::Gaus1core){

    auto rangelow = (mu-(rangeLow*sigma));
	 if( rangelow < -1.0 ) rangelow = -0.5;
    auto rangeup  = (mu+(rangeUp *sigma));
    if( rangeup > 1.0 ) rangeup = 0.5;
	 if( abs(mu) > 0.03 ) mu = -0.02;
    if( sigma > 1  ) sigma = 0.3;

    form = new TFormula(formname.c_str(),"[0]*exp(-0.5*((x-[1])/[2])**2)");
    fit  = new TF1(fitname.c_str(),form->GetName(),rangelow,rangeup);

    fit->SetParName(0,"N");      fit->SetParameter(0,norm); fit->SetParLimits(0,norm/2,norm*10);
    fit->SetParName(1,"#mu");    fit->SetParameter(1,mu); fit->SetParLimits(1,-0.03,0.03);
    fit->SetParName(2,"#sigma"); fit->SetParameter(2,sigma); fit->SetParLimits(2,0,1);

	 std::cout << "Prep " << ftypes[type] << " Fit: (" << rangelow << "t" << rangeup;   
	 std::cout << ") N: " << norm << "(" << norm/2 << "t" << norm*10 << ") Mu: " << mu << " (-0.03t0.03) s: " << sigma << " (0t1)" << std::endl; 

  } else if (type == TimeFitType::Gaus1coretail) { 

    form = new TFormula(formname.c_str(),"[0]*exp(-0.5*((x-[1])/[2])**2)");
    fit  = new TF1(fitname.c_str(),form->GetName(),rangelow,rangeup);

    mu = 0.0;
    fit->SetParName(0,"N");      fit->SetParameter(0,norm);
    fit->SetParName(1,"#mu");    fit->FixParameter(1,mu);
    fit->SetParName(2,"#sigma"); fit->SetParameter(2,sigma); fit->SetParLimits(2,0,10);

  } else if (type == TimeFitType::Gaus2fm || type == TimeFitType::Gaus2fmcore) { 

    form = new TFormula(formname.c_str(),"[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[1])/[4])**2)");
    fit  = new TF1(fitname.c_str(),form->GetName(),rangelow,rangeup);

    fit->SetParName(0,"N_{1}");      fit->SetParameter(0,norm);
    fit->SetParName(1,"#mu");        fit->SetParameter(1,mu);
    fit->SetParName(2,"#sigma_{1}"); fit->SetParameter(2,sigma);   fit->SetParLimits(2,0,10);
    fit->SetParName(3,"N_{2}");      fit->SetParameter(3,norm/10);
    fit->SetParName(4,"#sigma_{2}"); fit->SetParameter(4,sigma*4); fit->SetParLimits(4,0,10);

  } else if (type == TimeFitType::Gaus3fm || type == TimeFitType::Gaus3fmcore) {

    form = new TFormula(formname.c_str(),"[0]*exp(-0.5*((x-[1])/[2])**2)+[3]*exp(-0.5*((x-[1])/[4])**2)+[5]*exp(-0.5*((x-[1])/[6])**2)");
    fit  = new TF1(fitname.c_str(),form->GetName(),rangelow,rangeup);

    fit->SetParName(0,"N_{1}");      fit->SetParameter(0,norm*0.8);  fit->SetParLimits(0,norm*0.5,norm);
    fit->SetParName(1,"#mu");        fit->SetParameter(1,mu);
    fit->SetParName(2,"#sigma_{1}"); fit->SetParameter(2,sigma*0.7); fit->SetParLimits(2,sigma*0.5,sigma);
    fit->SetParName(3,"N_{2}");      fit->SetParameter(3,norm*0.3);  fit->SetParLimits(3,norm*0.1,norm*0.5);
    fit->SetParName(4,"#sigma_{2}"); fit->SetParameter(4,sigma*1.4); fit->SetParLimits(4,sigma,sigma*1.5);
    fit->SetParName(5,"N_{3}");      fit->SetParameter(5,norm*0.01); fit->SetParLimits(5,norm*0.005,norm*0.1);
    fit->SetParName(6,"#sigma_{3}"); fit->SetParameter(6,sigma*2.5); fit->SetParLimits(6,sigma*1.5,sigma*5.0);

  } else { std::cerr << "How did this happen?? Fit was not set for prepping fits! Exiting..." << std::endl; exit(1); }
  //--------------------------------------------------------------------------

}//<<>> void TimeFitStruct::PrepFit()

void TimeFitStruct::DoFit(){

  hist->Fit(fit->GetName(),"RBQM");

}//<<>>void TimeFitStruct::DoFit()

void TimeFitStruct::GetFitResult(){

  // get common result
  result.mu       = fit->GetParameter(1);
  result.emu      = fit->GetParError (1);
  result.chi2ndf  = fit->GetChisquare();
  result.chi2prob = fit->GetProb();

  if (type == TimeFitType::Gaus1 || type == TimeFitType::Gaus1core){

    result.sigma  = fit->GetParameter(2);
    result.esigma = fit->GetParError (2);

  } else if (type == TimeFitType::Gaus2fm || type == TimeFitType::Gaus2fmcore) {

    const Float_t const1 = fit->GetParameter(0); 
    const Float_t const2 = fit->GetParameter(3);
    const Float_t denom  = const1 + const2;

    result.sigma  = (const1*fit->GetParameter(2)+const2*fit->GetParameter(4))/denom;
    result.esigma = std::hypot(const1*fit->GetParError(2),const2*fit->GetParError(4))/denom;

  } else if (type == TimeFitType::Gaus3fm || type == TimeFitType::Gaus3fmcore) {

    const Double_t const1 = fit->GetParameter(0); 
    const Double_t const2 = fit->GetParameter(3);
    const Double_t const3 = fit->GetParameter(5);
    const Double_t denom  = const1 + const2 + const3;
    
    result.sigma  = (const1*fit->GetParameter(2) + const2*fit->GetParameter(4) + const3*fit->GetParameter(6))/denom;
    result.esigma = std::hypot(std::hypot(const1*fit->GetParError(2),const2*fit->GetParError(4)),const3*fit->GetParError(6))/denom; // need c++17...

  } else { std::cerr << "How did this happen?? Fit was not set for getting result! Exiting..." << std::endl; exit(1);}

}//<<>>void TimeFitStruct::GetFitResult()

void TimeFitStruct::DeleteInternal(){

  // skip dependent deletes if no entries
  if (!isEmpty()){ delete fit; delete form; }
  // always delete hist
  delete hist;

}//<<>>void TimeFitStruct::DeleteInternal()
