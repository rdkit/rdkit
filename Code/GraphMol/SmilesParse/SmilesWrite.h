//
//  Copyright (C) 2002-2006 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved  @@
//
#ifndef _RD_SMILESWRITE_H
#define _RD_SMILESWRITE_H

#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <string>
#include <vector>


namespace RDKit{
  namespace SmilesWrite {
    std::string GetAtomSmiles(const Atom *atom,bool doKekule=false);
    std::string GetBondSmiles(const Bond *bond,int atomToLeftIdx=-1,bool doKekule=false);
  } 
  
  std::string MolToSmiles(ROMol &mol,bool doIsomericSmiles=false,
			  bool doKekule=false,int rootedAtAtom=-1);
}
#endif
