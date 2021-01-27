#ifndef ParticleList_HH
#define ParticleList_HH

#include <vector>
#include "TLorentzVector.h"
#include "Particle.hh"
  
class ParticleList : public std::vector<Particle> {
public:
  ParticleList();
  
  virtual ~ParticleList();
  
  operator TLorentzVector() const;

  int Charge() const;

  ParticleList PtEtaCut(double pt, double eta = -1) const;
  ParticleList ParticleIDCut(ParticleIDType id) const;
  ParticleList RemoveOverlap(const ParticleList& parts, double deltaR = 0.2) const;

  
  ParticleList& SortByPt();

  ParticleList operator + (const ParticleList& parts) const;
  void operator += (const ParticleList& parts);
};

  

#endif
