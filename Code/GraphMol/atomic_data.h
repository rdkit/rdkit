//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//

/*! \file atomic_data.h

  \brief No user-serviceable parts inside

  This stuff is used by the PeriodicTable interface

*/  
#ifndef __RD_ATOMIC_DATA_H
#define __RD_ATOMIC_DATA_H

#include <RDGeneral/types.h>

namespace RDKit {
  extern const std::string periodicTableAtomData;
  
  class atomicData {
  public :
    atomicData(const std::string &dataLine);
    ~atomicData() {};
    
    int AtomicNum() const { return anum;};
    
    int DefaultValence() const { return valence.front();};
    
    int NumValence() const { return static_cast<int>(valence.size());};
    
    const INT_VECT &ValenceList() const {
      return valence;
    };
    
    double Mass() const { return mass;};
    
    std::string Symbol() const { 
      return symb;
    }

    double Rcov() const { return rCov; }
    
    double Rb0() const {return rB0;}
    
    double Rvdw() const { return rVdw;}
    
    int NumOuterShellElec() const { return nVal;}
    
  private:
    int anum; //atomic number
    std::string symb; // atomic symbol
    double rCov, rB0, rVdw; //radii
    INT_VECT valence; //list of all valences, the first one is the default valence, -1 at the end signifies that any upper valence is tolerated
    double mass;  // atomic mass
    int nVal; // number of outer shell electrons
  };
  
};
#endif
