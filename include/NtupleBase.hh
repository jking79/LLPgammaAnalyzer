#ifndef NtupleBase_h
#define NtupleBase_h

#include "AnalysisBase.hh"

template <class Base>
class NtupleBase : public AnalysisBase<Base> {

public:
  NtupleBase(TTree* tree = 0);
  virtual ~NtupleBase();

  void WriteNtuple(const std::string& filename);

protected:
  std::vector<TTree*>     m_Trees;
  std::map<string,TTree*> m_Label2Tree;
 

private:
  virtual TTree* InitOutputTree(const std::string& sample) = 0;
  virtual void FillOutputTree(TTree* tree) = 0;
  
  

};

#endif
