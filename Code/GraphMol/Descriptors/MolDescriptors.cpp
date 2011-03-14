// $Id$
//
//  Copyright (C) 2005-2008 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include "MolDescriptors.h"

namespace RDKit{
  namespace Descriptors {

    static std::string _amwVersion="1.0.0";
    double calcAMW(const ROMol &mol,bool onlyHeavy){
      double res=0.0;
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt != mol.endAtoms();++atomIt){
	int atNum=(*atomIt)->getAtomicNum();
	if(atNum!=1 || !onlyHeavy){
	  res += (*atomIt)->getMass();
	}

	// add our implicit Hs if we need to:
	if(!onlyHeavy){
      const PeriodicTable *table=PeriodicTable::getTable();
	  res += (*atomIt)->getTotalNumHs()*table->getAtomicWeight(1);
	}
      }
      return res;
    }
  } // end of namespace Descriptors
} // end of namespace RDKit
