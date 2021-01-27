#ifndef SampleSet_HH
#define SampleSet_HH

#include <string>

using std::string;

vector<string> g_File;
vector<string> g_Tree;
vector<string> g_Title;
vector<bool> g_Bkg;
vector<int> g_Color;
vector<int> g_Hist;

class SampleSet {
public:
  SampleSet();
    
  virtual ~SampleSet();

  void SetBkg(bool is_bkg);
  bool GetBkg() const;
  
  void   AddFile(const string& filename);
  int    GetNFile() const;
  string GetFile(int n);
  
  void   SetTitle(const string& title);
  string GetTitle() const;

  void   SetTreeName(const string& treename);
  string GetTreeName() const;
  
  void SetColor(int icolor);
  int  GetColor() const;
  
  void SetSkip(int iskip);
  int  GetSkip() const;
  
  void   SetScale(double scale);
  double GetScale() const;
  
private:
  bool m_IsBkg;
  std::vector<string> m_FileNames;
  string m_TreeName;
  string m_Title;
  int m_Color;
  int m_Skip;
  double m_Scale;

};

#endif

inline SampleSet::SampleSet(){
  m_IsBkg = true;
  m_Title = "";
  m_TreeName = "KUAnalysis";
  m_Color = kBlue;
  m_Skip = 1;
  m_Scale = 1.;
}
    
inline SampleSet::~SampleSet() {}

inline void SampleSet::SetBkg(bool is_bkg){
  m_IsBkg = is_bkg;
}

inline bool SampleSet::GetBkg() const {
  return m_IsBkg;
}
  
inline void SampleSet::AddFile(const string& filename){
  m_FileNames.push_back(filename);
}

inline int SampleSet::GetNFile() const {
  return m_FileNames.size();
}

inline string SampleSet::GetFile(int n){
  int N = GetNFile();
  if(n < 0 || n >= N)
    return "NO FILE";
  return m_FileNames[n];
}
  
inline void SampleSet::SetTitle(const string& title){
  m_Title = title;
}

inline string SampleSet::GetTitle() const {
  return m_Title;
}

inline void SampleSet::SetTreeName(const string& treename){
  m_TreeName = treename;
}
inline string SampleSet::GetTreeName() const {
  return m_TreeName;
}
  
inline void SampleSet::SetColor(int icolor){
  m_Color = icolor;
}

inline int SampleSet::GetColor() const {
  return m_Color;
}
  
inline void SampleSet::SetSkip(int iskip){
  m_Skip = iskip;
}

inline int SampleSet::GetSkip() const {
  return m_Skip;
}
  
inline void SampleSet::SetScale(double scale){
  m_Scale = scale;
}

inline double SampleSet::GetScale() const {
  return m_Scale;
}
