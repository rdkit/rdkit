// $Id$
//
//  Copyright (C) 2005-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#include <RDGeneral/Invariant.h>
#include <GraphMol/RDKitBase.h>
#include "MolDescriptors.h"

namespace RDKit{
  namespace Descriptors {

    static std::string _amwVersion="1.0.0";
    double CalcAMW(const ROMol &mol,bool onlyHeavy){
      double res=0.0;
      const PeriodicTable *table=PeriodicTable::getTable();
      for(ROMol::ConstAtomIterator atomIt=mol.beginAtoms();
	  atomIt != mol.endAtoms();++atomIt){
	int atNum=(*atomIt)->getAtomicNum();
	if(atNum!=1 || !onlyHeavy){
	  res += table->getAtomicWeight(atNum);
	}

	// add our implicit Hs if we need to:
	if(!onlyHeavy){
	  res += (*atomIt)->getImplicitValence()*table->getAtomicWeight(1);
	}
      }
      return res;
    }
  } // end of namespace Descriptors
} // end of namespace RDKit
