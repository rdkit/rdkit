//
//  Copyright (C) 2023 Rocco Moretti and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
// This file implements the non-contiguous atom molecular structural similarity (NAMS)
// method of Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)
// for an RDKit context

#include <RDGeneral/export.h>
#ifndef __RD_NAMS_H__
#define __RD_NAMS_H__

namespace RDKit {
class ROMol;
namespace NAMS {

/*!
  The NAMS MolInfo object functions analogously to a "fingerprint" for a
  particular molecule, though it doesn't integrate with other fingerprint code,
  and is only used by the other NAMS functionalities
*/
class RDKIT_FINGERPRINTS_EXPORT NAMSMolInfo {

};

//! returns a NAMS MolInfo object for a molecule
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

  \param mol:    the molecule to be fingerprinted
   
  \return a pointer to the molinfo object. The client is
  responsible for calling delete on this.
*/
RDKIT_FINGERPRINTS_EXPORT NAMSMolInfo * getNAMSMolInfo(const ROMol &mol);

//! returns the NAMS similarity between two molecules (encoded in NAMS MolInfo objects)
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

  The similarity score should be invariant to the order of molecule listing.

  \param molinfo1:   one of the molinfos to use
  \param molinfo2:   one of the molinfos to use
   
  \return the similarity between the two molecules (on a scale between 0-1)
*/
RDKIT_FINGERPRINTS_EXPORT double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2);

} // namespace NAMS
} // namespace RDKit
#endif
