//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef RD_FRAG_CAT_GENERATOR_H
#define RD_FRAG_CAT_GENERATOR_H

#include <Catalogs/Catalog.h>
#include "FragCatalogEntry.h"
#include "FragCatParams.h"
#include <GraphMol/Subgraphs/Subgraphs.h>

namespace RDKit {
  class ROMol;
  
  typedef RDCatalog::HierarchCatalog<FragCatalogEntry, FragCatParams, int> FragCatalog;


  class FragCatGenerator {
  public:
    
    FragCatGenerator() {}
   
    unsigned int addFragsFromMol(const ROMol &mol, FragCatalog *fcat);
 };
}

#endif    
