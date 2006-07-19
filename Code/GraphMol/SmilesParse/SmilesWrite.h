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

namespace SmilesWrite {
  std::string GetAtomSmiles(const RDKit::Atom *atom,bool doKekule=false);
}

namespace RDKit{
  std::string MolToSmiles(ROMol &mol,bool doIsomericSmiles=false,
			  bool doKekule=false);
}
#endif
