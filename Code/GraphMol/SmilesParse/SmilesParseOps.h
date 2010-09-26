//
//  Copyright (C) 2001-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_SMILESPARSEOPS_H
#define _RD_SMILESPARSEOPS_H
#include <GraphMol/Bond.h>

namespace RDKit{
  class RWMol;
  class Atom;
}
namespace SmilesParseOps {
  void ReportParseError(const char *message,bool throwIt=true);
  void CleanupAfterParseError(RDKit::RWMol *mol);
  void AddFragToMol(RDKit::RWMol *mol,RDKit::RWMol *frag,
		    RDKit::Bond::BondType bondOrder=RDKit::Bond::UNSPECIFIED,
		    RDKit::Bond::BondDir bondDir=RDKit::Bond::NONE,
		    bool closeRings=false,bool doingQuery=false);
  RDKit::Bond::BondType GetUnspecifiedBondType(const RDKit::RWMol *mol,
					       const RDKit::Atom *atom1,
					       const RDKit::Atom *atom2);
  void CloseMolRings(RDKit::RWMol *mol,bool toleratePartials);
  void AdjustAtomChiralityFlags(RDKit::RWMol *mol);
};

#endif
