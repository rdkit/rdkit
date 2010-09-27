//
//  Copyright (C) 2003-2006 Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#ifndef _RD_FRAG_CATALOG_UTILS_H_
#define _RD_FRAG_CATALOG_UTILS_H_

#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include "FragCatParams.h"
#include <iostream>

namespace RDKit {
  
  // get the functional groups from file or stream
  // each functional groups is read in as a molecule with queryatoms and 
  // querybonds
  MOL_SPTR_VECT readFuncGroups(std::string fileName);
  MOL_SPTR_VECT readFuncGroups(std::istream &inStream,int nToRead=-1);
  
  // REVIEW: should this return a vector of pairs or a map?
  // mark the functional groups of interest on the molecule
  // and return a vector os std::pair <aid, fid>
  //    aid - is the atom id in mol that connect to a functional (fid)
  //    fid - the functional groups in the list maintained in params
  // ARGUMENTS:
  //  mol - molecule of interest
  //  params - fragment catalog paramter object (contains a list of functional
  //             groups of interest
  //  fgBonds - container for bondIds in mol that are part of the functional groups 
  //            the connection bond is included. these need to be chopped from 
  //            the molecule later
  
  MatchVectType findFuncGroupsOnMol(const ROMol &mol, 
				 const FragCatParams *params,
				 INT_VECT &fgBonds);
 
  // This functions is called before either adding the fragments from a molecule
  // to a fragment catalog or generating the fincgerprint for this molecule
  // using a fragment catalog. These are the things this function does
  // - recognize the function groups (and their location) on the molecule
  // - chop these functional groups of the molecule to create a core molecule
  //   "coreMol"
  // - map the function group locations onto this "coreMol" (bacause the atom ids
  //   on coreMol are different from the original molecule
  // - return coreMol to the caller of this function and the enter the atom ids to func
  //   group ids mapping into aToFmap argument 
  ROMol *prepareMol(const ROMol &mol, const FragCatParams *fparams,
		    MatchVectType &aToFmap);

}

#endif
