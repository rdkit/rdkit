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
// method of Teixeira & Falcao (http://dx.doi.org/10.1021/ci400324u)
// for an RDKit context

#include <RDGeneral/export.h>
#ifndef __RD_NAMS_H__
#define __RD_NAMS_H__

#include <vector>
#include <string>

namespace RDKit {
class ROMol;
namespace NAMS {

struct NAMSParameters {

  NAMSParameters() = default;
  ~NAMSParameters();

  // bool include_chirality = false; // Should we include the chirality in matching?
  float BS_ALPHA = 0.9f;
  float ANRINGS_FAC = 0.8f;	//  #number of rings an atom belongs to
  float ACHIR_FAC = 0.95f;	//  #chiral atom
  float DBSTEREO_FAC = 0.95f;	//  #double bond stereo
  float BRING_FAC = 0.9f;	//  #bond in ring
  float BAROM_FAC = 0.9f;	//  #bond aromaticity
  float BORDER_FAC = 0.9f;	//  #bond order
  float PEN = 1.0f;
  int ADM = 3;

  static constexpr int MAX_LEVELS = 128;

  float ELEMS_DISTS(int ele1, int ele2) const;

  // These are non-owning raw pointers (do not attempt to delete/free).
  const int* getBondLevelsMatrix() const;

  //! Get a default NAMSParameters object;
  static const NAMSParameters & getDefault();

private:
  void calcBondLevelsMatrix() const;

  mutable int * blev_mat = nullptr;
  mutable float blev_alpha = 0; // The BS_ALPHA alpha level the blev_mat is calculated for.
};

struct NAMSResult {
  double self_similarity1 = -1, self_similarity2 = -1;
  double similarity = -1;
  double jaccard = -1;
  std::vector< int > mapping1to2;
  std::vector< double > atom_scores;
};

/*!
  Not exposed to Python, this is an internal implementation detail of the NAMSMolInfo type
*/
struct BondInfoType {
  BondInfoType() = default;
  BondInfoType(const ROMol &mol, unsigned int at1, unsigned int at2, bool do_isomerism = false);

  int ele1=0, ele2=0;
  int nrings1=0, nrings2=0; // Needs to be int due to subtraction
  int chir1=0, chir2=0;
  bool inring=false, aromatic=false;
  double order=0;
  int dbcistrans=0;

  bool operator<( const BondInfoType & r ) const;
};

/*!
  The NAMS MolInfo object functions analogously to a "fingerprint" for a
  particular molecule, though it doesn't integrate with other fingerprint code,
  and is only used by the other NAMS functionalities
*/
class RDKIT_FINGERPRINTS_EXPORT NAMSMolInfo {
public:
  NAMSMolInfo(const ROMol &mol, bool include_chirality=false);
  NAMSMolInfo(const ROMol &mol, const NAMSParameters & parms);

  unsigned int natoms() const;
  unsigned int nbonds() const;
  unsigned int naba_types() const;

  std::vector< BondInfoType > aba_types;
  // natoms by bond by aba_type index
  std::vector< std::vector< unsigned int > > mat_aba_types;
  // natoms by bond by level number
  std::vector< std::vector< int > > mat_levels;

  double self_similarity = -1;
};

//! returns a NAMS MolInfo object for a molecule
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

  \param mol:    the molecule to be fingerprinted

  \return a pointer to the molinfo object. The client is
  responsible for calling delete on this.
*/
RDKIT_FINGERPRINTS_EXPORT NAMSMolInfo * getNAMSMolInfo(const ROMol &mol);

//! returns a NAMS MolInfo object for a molecule
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

  \param mol:    the molecule to be fingerprinted

  \return a pointer to the molinfo object. The client is
  responsible for calling delete on this.
*/
RDKIT_FINGERPRINTS_EXPORT NAMSMolInfo * getNAMSMolInfo(const ROMol &mol, const NAMSParameters & parms);

//! returns the NAMS similarity between two molecules (encoded in NAMS MolInfo objects)
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)
  Default parameters will be used

  The similarity score should be invariant to the order of molecule listing.

  \param molinfo1:   one of the molinfos to use
  \param molinfo2:   one of the molinfos to use

  \return the similarity between the two molecules (on a scale between 0-1)
*/
RDKIT_FINGERPRINTS_EXPORT double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2);

//! returns the NAMS similarity between two molecules (encoded in NAMS MolInfo objects)
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

  The similarity score should be invariant to the order of molecule listing.

  \param molinfo1:   one of the molinfos to use
  \param molinfo2:   one of the molinfos to use

  \return the similarity between the two molecules (on a scale between 0-1)
*/
RDKIT_FINGERPRINTS_EXPORT double getNAMSSimilarity(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2, const NAMSParameters & params);

//! returns the NAMS similarity result between two molecules (encoded in NAMS MolInfo objects)
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)
  Default parameters will be used

  The similarity score should be invariant to the order of molecule listing.

  \param molinfo1:   one of the molinfos to use
  \param molinfo2:   one of the molinfos to use
  \param params:     The NAMSParameters object to control how to run the calculation

  \return a pointer to the NAMSResult object which contains detailed information about the similarity calculation.
  The client is responsible for calling delete on this.
*/
RDKIT_FINGERPRINTS_EXPORT NAMSResult * getNAMSResult(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2);

//! returns the NAMS similarity result between two molecules (encoded in NAMS MolInfo objects)
/*!
  The NAMS algorithm is described by Teixeira & Falcao (https://pubs.acs.org/doi/abs/10.1021/ci400324u)

  The similarity score should be invariant to the order of molecule listing.

  \param molinfo1:   one of the molinfos to use
  \param molinfo2:   one of the molinfos to use
  \param params:     The NAMSParameters object to control how to run the calculation

  \return a pointer to the NAMSResult object which contains detailed information about the similarity calculation.
  The client is responsible for calling delete on this.
*/
RDKIT_FINGERPRINTS_EXPORT NAMSResult * getNAMSResult(const NAMSMolInfo & molinfo1, const NAMSMolInfo & molinfo2, const NAMSParameters & params);

} // namespace NAMS
} // namespace RDKit
#endif
