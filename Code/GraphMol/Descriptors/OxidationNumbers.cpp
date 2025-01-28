//
//  Copyright (c) 2023, David Cosgrove, CozChemIx Limited
//  All rights reserved.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
//
// Calculate the oxidation numbers (states) of the atoms in a molecule.
// Based on the code at
// https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/src/linchemin/cheminfo/functions.py#L544
// and therefore also subject to the MIT licence as detailed at
// https://github.com/syngenta/linchemin/blob/f44fda38e856eaa876483c94284ee6788d2c27f4/LICENSE

#include <GraphMol/RDKitBase.h>
#include <GraphMol/Descriptors/OxidationNumbers.h>

namespace RDKit {
namespace Descriptors {

namespace {
int calcOxidationNumberByEN(const Atom *atom) {
  const static std::map<int, float> pauling_en_map = {
      {1, 2.2},   {3, 0.98},  {4, 1.57},  {5, 2.04},  {6, 2.55},  {7, 3.04},
      {8, 3.44},  {9, 3.98},  {11, 0.93}, {12, 1.31}, {13, 1.61}, {14, 1.9},
      {15, 2.19}, {16, 2.58}, {17, 3.16}, {19, 0.82}, {20, 1},    {21, 1.36},
      {22, 1.54}, {23, 1.63}, {24, 1.66}, {25, 1.55}, {26, 1.83}, {27, 1.88},
      {28, 1.91}, {29, 1.9},  {30, 1.65}, {31, 1.81}, {32, 2.01}, {33, 2.18},
      {34, 2.55}, {35, 2.96}, {36, 3.0},  {37, 0.82}, {38, 0.95}, {39, 1.22},
      {40, 1.33}, {41, 1.6},  {42, 2.16}, {43, 1.9},  {44, 2.2},  {45, 2.28},
      {46, 2.2},  {47, 1.93}, {48, 1.69}, {49, 1.78}, {50, 1.96}, {51, 2.05},
      {52, 2.1},  {53, 2.66}, {54, 2.6},  {55, 0.79}, {56, 0.89}, {57, 1.1},
      {58, 1.12}, {59, 1.13}, {60, 1.14}, {62, 1.17}, {64, 1.2},  {66, 1.22},
      {67, 1.23}, {68, 1.24}, {69, 1.25}, {71, 1.27}, {72, 1.3},  {73, 1.5},
      {74, 2.36}, {75, 1.9},  {76, 2.2},  {77, 2.2},  {78, 2.28}, {79, 2.54},
      {80, 2.0},  {81, 1.62}, {82, 2.33}, {83, 2.02}, {84, 2.0},  {85, 2.2},
      {88, 0.9},  {89, 1.1},  {90, 1.3},  {91, 1.5},  {92, 1.38}, {93, 1.36},
      {94, 1.28}, {95, 1.3},  {96, 1.3},  {97, 1.3},  {98, 1.3},  {99, 1.3},
      {100, 1.3}, {101, 1.3}, {102, 1.39}};
  PRECONDITION(atom != nullptr, "must have an atom");
  PRECONDITION(atom->hasOwningMol(), "atom must have owning molecule")
  auto get_en = [&](int atomicNum) -> float {
    auto res = pauling_en_map.find(atomicNum);
    return res == pauling_en_map.end() ? 0.0 : res->second;
  };
  auto sf = [](float enDiff) -> int {
    if (enDiff > 0.0) {
      return -1;
    } else if (enDiff < 0.0) {
      return 1;
    } else {
      return 0;
    }
  };

  int oxNum = 0;
  float parEN = get_en(atom->getAtomicNum());

  for (const auto &bond : atom->getOwningMol().atomBonds(atom)) {
    if (bond->getBondType() == Bond::DATIVE ||
        bond->getBondType() == Bond::DATIVEONE ||
        bond->getBondType() == Bond::DATIVEL ||
        bond->getBondType() == Bond::DATIVER ||
        bond->getBondType() == Bond::ZERO) {
      continue;
    }
    auto otherAtom = bond->getOtherAtom(atom);
    if (otherAtom->getAtomicNum() > 1) {
      float en_diff = parEN - get_en(otherAtom->getAtomicNum());
      double bondType = bond->getBondTypeAsDouble();
      // Make sure this is a kekulized mol i.e. no bond type of 1.5.  This
      // shouldn't happen if called from calcOxidationNumbers, but who knows
      // what might happen in future.
      if (bondType > 1.0 && bondType < 2.0) {
        throw ValueErrorException(
            "Molecule appears not to be Kekulized,"
            " oxidation number calculation fails.");
      }
      oxNum += bondType * sf(en_diff);
    }
  }
  oxNum += sf(parEN - get_en(1)) * atom->getTotalNumHs();
  oxNum += atom->getFormalCharge();
  return oxNum;
}
}  // namespace

void calcOxidationNumbers(const ROMol &mol) {
  RWMol molCp(mol);
  RDKit::MolOps::Kekulize(molCp);
  for (const auto &atom : mol.atoms()) {
    auto cpAtom = molCp.getAtomWithIdx(atom->getIdx());
    int oxNum = calcOxidationNumberByEN(cpAtom);
    atom->setProp<int>(common_properties::OxidationNumber, oxNum);
  }
}

}  // end of namespace Descriptors
}  // end of namespace RDKit
