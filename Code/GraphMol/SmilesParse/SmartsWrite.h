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
    std::string GetAtomSmarts(const QueryAtom *qatom);
    std::string GetBondSmarts(const QueryBond *qbond);
  }

  class ROMol;
  std::string MolToSmarts(ROMol &mol,bool doIsomericSmarts=false);
};

#endif
