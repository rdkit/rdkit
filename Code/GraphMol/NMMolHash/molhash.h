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
struct RDKIT_NMMOLHASHLIB_EXPORT HashFunction {
  static const unsigned int AnonymousGraph = 1;
  static const unsigned int ElementGraph = 2;
  static const unsigned int CanonicalSmiles = 3;
  static const unsigned int MurckoScaffold = 4;
  static const unsigned int ExtendedMurcko = 5;
  static const unsigned int MolFormula = 6;
  static const unsigned int AtomBondCounts = 7;
  static const unsigned int DegreeVector = 8;
  static const unsigned int Mesomer = 9;
  static const unsigned int HetAtomTautomer = 10;
  static const unsigned int HetAtomProtomer = 11;
  static const unsigned int RedoxPair = 12;
  static const unsigned int Regioisomer = 13;
  static const unsigned int NetCharge = 14;
  static const unsigned int SmallWorldIndexBR = 15;
  static const unsigned int SmallWorldIndexBRL = 16;
  static const unsigned int ArthorSubstructureOrder = 17;
};

RDKIT_NMMOLHASHLIB_EXPORT std::string MolHash(RWMol *mol, unsigned int func);

struct RDKIT_NMMOLHASHLIB_EXPORT StripType {
  static const unsigned int AtomStereo = 1;
  static const unsigned int BondStereo = 2;
  static const unsigned int Isotope = 4;
  static const unsigned int AtomMap = 8;
  static const unsigned int Hydrogen = 16;
};

RDKIT_NMMOLHASHLIB_EXPORT void Strip(RWMol *mol, unsigned int striptype);
RDKIT_NMMOLHASHLIB_EXPORT void SplitMolecule(RWMol *mol,
                                             std::vector<RWMol *> &molv);
}  // namespace MolHash
}  // namespace RDKit
#endif  // NMS_MOLHASH_H
