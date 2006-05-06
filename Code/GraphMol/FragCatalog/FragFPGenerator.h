//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_FRAG_FP_GENERATOR_H_
#define _RD_FRAG_FP_GENERATOR_H_

#include <vector>
#include <Catalogs/Catalog.h>
#include "FragCatalogEntry.h"
#include "FragCatParams.h"

class ExplicitBitVect;
namespace RDKit {
  class ROMol;
  typedef RDCatalog::HierarchCatalog<FragCatalogEntry, FragCatParams, int> FragCatalog;
  typedef std::vector< std::pair<int,int> > MatchVectType; 


  class FragFPGenerator {
  public:
    FragFPGenerator() {}

    ExplicitBitVect *getFPForMol(const ROMol &mol, const FragCatalog &fcat);

  private:
    void computeFP(const ROMol &mol, const FragCatalog &fcat,
		   const MatchVectType &aidToFid, ExplicitBitVect *fp);
  };
}

#endif
