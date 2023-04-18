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

#include <vector>
#include <string>

namespace RDKit {
class ROMol;
namespace NAMS {

/*!
  Not exposed to Python, this is an internal implementation detail of the NAMSMolInfo type
*/
struct BondInfoType {
  BondInfoType() = default;
  BondInfoType(const ROMol &mol, unsigned int at1, unsigned int at2, bool do_isomerism = false);

  int ele1=0, ele2=0;
  unsigned int nring1=0, nring2=0;
  int chr1=0, chr2=0;
  bool ringbond=false, aromatic=false;
  double order=0;
  int dbstero12=0;

  bool operator<( const BondInfoType & r ) const;
};

/*!
  The NAMS MolInfo object functions analogously to a "fingerprint" for a
  particular molecule, though it doesn't integrate with other fingerprint code,
  and is only used by the other NAMS functionalities
*/
class RDKIT_FINGERPRINTS_EXPORT NAMSMolInfo {
public:
  NAMSMolInfo(const ROMol &mol);

  std::string smiles; // Needed?
  double molwt=99.9; // Needed?
  std::vector< BondInfoType > aba_types;
  // natoms by bond by aba_type index
  std::vector< std::vector< unsigned int > > mat_aba_types;
  // natoms by bond by level number
  std::vector< std::vector< int > > mat_levels;

  //! Primarily for debugging purposes.
  //! Dump data in NAMS database format (more or less)
  std::string dump(int cid=1) const; // For debugging.
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
