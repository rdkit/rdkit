//
//  Copyright (C) 2004-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SMARTSWRITE_H
#define _RD_SMARTSWRITE_H

#include <string>

namespace RDKit {
  class ROMol;
  class QueryAtom;
  class QueryBond;

  std::string MolToSmarts(ROMol &mol,bool doIsomericSmarts=false);
  std::string GetAtomSmarts(const QueryAtom *qatom);
  std::string GetBondSmarts(const QueryBond *qbond);
};

#endif
