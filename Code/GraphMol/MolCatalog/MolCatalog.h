//
//  Copyright (C) 2006 Greg Landrum
//
//
#ifndef _RD_MOL_CATALOG_H_
#define _RD_MOL_CATALOG_H_

#include <Catalogs/Catalog.h>
#include <GraphMol/MolCatalog/MolCatalogEntry.h>
#include <GraphMol/MolCatalog/MolCatalogParams.h>

namespace RDKit {

  //! a hierarchical catalog for holding molecules
  typedef RDCatalog::HierarchCatalog<MolCatalogEntry, MolCatalogParams, int> MolCatalog;
}

#endif    
