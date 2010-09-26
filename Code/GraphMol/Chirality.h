//
//  Copyright (C) 2008 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Chirality.h

  \brief Not intended for client-code use.

*/
#ifndef _RD_CHIRALITY_20AUG2008_H_
#define _RD_CHIRALITY_20AUG2008_H_
#include <RDGeneral/types.h>

/// @cond
namespace RDKit{
  class ROMol;
  namespace Chirality {
    /*!
      \param mol the molecule to be altered
      \param ranks  used to return the set of ranks.  
                    Should be at least mol.getNumAtoms() long.
    
      <b>Notes:</b>
         - All atoms gain a property "_CIPRank" with their overall
           CIP ranking.
    
    */
    void assignAtomCIPRanks(const ROMol &mol,INT_VECT &ranks);
  }
}
/// @endcond
#endif
