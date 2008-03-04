//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SMARTSWRITE_H
#define _RD_SMARTSWRITE_H

#include <string>

namespace RDKit {
  class QueryAtom;
  class QueryBond;
  namespace SmartsWrite {
    //! returns the SMARTS for a QueryAtom
    std::string GetAtomSmarts(const QueryAtom *qatom);
    //! returns the SMARTS for a QueryBond
    std::string GetBondSmarts(const QueryBond *qbond);
  }

  class ROMol;
  //! returns the SMARTS for a molecule
  std::string MolToSmarts(ROMol &mol,bool doIsomericSmarts=false);
};

#endif
