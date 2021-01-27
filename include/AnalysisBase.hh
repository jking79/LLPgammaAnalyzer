#ifndef AnalysisBase_h
#define AnalysisBase_h

#include <iostream>

#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1D.h>
#include <string>

#include "NeventTool.hh"
#include "XsecTool.hh"

using namespace std;

class ParticleList;

template <class Base>
class AnalysisBase : public Base {

public:
  AnalysisBase(TTree* tree = 0);
  virtual ~AnalysisBase();

  void AddLabels(const string& dataset, const string& filetag);
  void AddEventCountFile(const string& rootfile);
  void DoSMS(){ m_DoSMS = true; }

  string GetEntry(int entry);

  // analysis functions
  virtual TVector3 GetMET();
  virtual ParticleList GetJets();
  virtual ParticleList GetElectrons();
  virtual ParticleList GetMuons();

  virtual TVector3 GetGenMET();
  virtual ParticleList GetGenElectrons();
  virtual ParticleList GetGenMuons();
  virtual ParticleList GetGenNeutrinos();
  virtual ParticleList GetGenBosons();
  virtual ParticleList GetGenSparticles();
 
  double DeltaPhiMin(const vector<TLorentzVector>& JETs, const TVector3& MET, int N = -1);
  double DeltaPhiMin(const vector<pair<TLorentzVector, bool> >& JETs, const TVector3& MET, int N = -1);
  
  void MomTensorCalc(vector<TLorentzVector>& input, vector<double>& eigenvalues, double pow = 1., bool threeD = true); 

protected:
  bool m_DoSMS;
  
  virtual double GetEventWeight();
  virtual double GetXsec();


private:
  string m_DataSet;
  string m_FileTag;

  NeventTool m_NeventTool;
  XsecTool   m_XsecTool;

  int m_SampleIndex;
  virtual int GetSampleIndex();
  int m_Nsample;
  std::map<int,int>         m_HashToIndex;
  std::map<int,std::string> m_IndexToSample;
  std::map<int,double>      m_IndexToXsec;
  std::map<int,double>      m_IndexToNevent;
  std::map<int,double>      m_IndexToNweight;
  
};

#endif









