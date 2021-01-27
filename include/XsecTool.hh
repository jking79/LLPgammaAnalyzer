#ifndef XsecTool_h
#define XsecTool_h

#include <iostream>
#include <string>
#include <map>
#include <vector>

class XsecTool {

public:
  XsecTool();
  virtual ~XsecTool();

  double GetXsec_BKG(const std::string& dataset) const;
  double GetXsec_SMS(const std::string& dataset, double MP) const;

private:
  static std::map<std::string,double> m_Label2Xsec_BKG;
  static std::map<std::string,double> InitMap_Xsec_BKG();

  static std::map<std::string,int> m_N_SMS;
  static std::map<std::string,int> InitMap_N_SMS();
  static std::map<std::string,std::vector<double> > m_Label2Mass_SMS;
  static std::map<std::string,std::vector<double> > InitMap_Mass_SMS();
  static std::map<std::string,std::vector<double> > m_Label2Xsec_SMS;
  static std::map<std::string,std::vector<double> > InitMap_Xsec_SMS();
  static std::map<std::string,std::vector<double> > m_Label2XsecUnc_SMS;
  static std::map<std::string,std::vector<double> > InitMap_XsecUnc_SMS();
  
};

#endif



