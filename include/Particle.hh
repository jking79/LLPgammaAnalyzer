#ifndef Particle_HH
#define Particle_HH

#include "TLorentzVector.h"

class ParticleList;

/// Particle ID level
enum ParticleIDType { kNothing, kVeryLoose, kLoose, kMedium, kTight, kVeryTight };

class Particle : public TLorentzVector {
public:
  Particle();
    
  virtual ~Particle();

  int Charge() const;
  void SetCharge(int charge);

  int PDGID() const;
  void SetPDGID(int pdgid);

  int MomPDGID() const;
  void SetMomPDGID(int pdgid);

  ParticleIDType ParticleID() const;
  void SetParticleID(ParticleIDType id);

  double RelIso() const;
  double MiniIso() const;
  void SetRelIso(double iso);
  void SetMiniIso(double iso);

  double Btag() const;
  void SetBtag(double btag);
    
  operator ParticleList() const;
  ParticleList operator + (const Particle& part) const; 
  ParticleList operator + (const ParticleList& parts) const; 
    
private:
  int m_Charge;
  int m_PDGID;
  int m_MomPDGID;
  ParticleIDType m_ParticleID;

  double m_RelIso;
  double m_MiniIso;

  double m_Btag;

};

bool sortbypt(const Particle& p1, const Particle& p2);


#endif
