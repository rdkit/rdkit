/*==============================================*/
/* Copyright (C)  2019       NextMove Software  */
/* All rights reserved.                         */
/*                                              */
/* This file is part of molhash.                */
/*                                              */
/* The contents are covered by the terms of the */
/* BSD license, which is included in the file   */
/* license.txt.                                 */
/*==============================================*/
#include <RDGeneral/export.h>
#ifndef NMS_MOLHASH_H
#define NMS_MOLHASH_H

#include <string>
#include <vector>

namespace RDKit {
class RWMol;
namespace MolHash {
enum class RDKIT_NMMOLHASHLIB_EXPORT HashFunction {
  AnonymousGraph = 1,
  ElementGraph = 2,
  CanonicalSmiles = 3,
  MurckoScaffold = 4,
  ExtendedMurcko = 5,
  MolFormula = 6,
  AtomBondCounts = 7,
  DegreeVector = 8,
  Mesomer = 9,
  HetAtomTautomer = 10,
  HetAtomProtomer = 11,
  RedoxPair = 12,
  Regioisomer = 13,
  NetCharge = 14,
  SmallWorldIndexBR = 15,
  SmallWorldIndexBRL = 16,
  ArthorSubstructureOrder = 17
};

RDKIT_NMMOLHASHLIB_EXPORT std::string MolHash(RWMol *mol,
                                              enum HashFunction func);

enum class RDKIT_NMMOLHASHLIB_EXPORT StripType {
  AtomStereo = 1,
  BondStereo = 2,
  Isotope = 4,
  AtomMap = 8,
  Hydrogen = 16
};

RDKIT_NMMOLHASHLIB_EXPORT void Strip(RWMol *mol, unsigned int striptype);
RDKIT_NMMOLHASHLIB_EXPORT void SplitMolecule(RWMol *mol,
                                             std::vector<RWMol *> &molv);
}  // namespace MolHash
}  // namespace RDKit
#endif  // NMS_MOLHASH_H
