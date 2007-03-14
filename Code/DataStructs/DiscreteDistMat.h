//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//  @@ All Rights Reserved @@
//
#ifndef __RD_DISCRETEDISTMAT_H__
#define __RD_DISCRETEDISTMAT_H__
#include "DiscreteValueVect.h"

namespace RDKit{
  class DiscreteDistMat {
  public:
    DiscreteDistMat();
    ~DiscreteDistMat(){};
    unsigned int getDist(unsigned char v1, 
                         unsigned char v2, 
                         DiscreteValueVect::DiscreteValueType type);

  private:
    unsigned int d_oneBitTab[256*256];
    unsigned int d_twoBitTab[256*256];
    unsigned int d_fourBitTab[256*256];
  
  };
  extern DiscreteDistMat *getDiscreteDistMat();
}
#endif
