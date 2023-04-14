//
//  Copyright (C) 2016-2022 Novartis Institutes for BioMedical Research and
//  other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "../../RDGeneral/types.h"
#include "../../Geometry/point.h"
#include "Utilites.h"
#include <algorithm>
#include <tuple>

namespace RDKit {
namespace StructureCheck {

void SetupNeighbourhood(const ROMol &mol,
                        std::vector<Neighbourhood> &neighbours) {
  neighbours.clear();
  neighbours.resize(mol.getNumAtoms());

  for (unsigned i = 0; i < mol.getNumBonds(); i++) {
    const Bond *bond = mol.getBondWithIdx(i);
    unsigned a1 = bond->getBeginAtomIdx();
    unsigned a2 = bond->getEndAtomIdx();

    neighbours[a1].Atoms.push_back(a2);
    neighbours[a1].Bonds.push_back(i);

    neighbours[a2].Atoms.push_back(a1);
    neighbours[a2].Bonds.push_back(i);
  }
}

bool getMolAtomPoints(const ROMol &mol, std::vector<RDGeom::Point3D> &atomPoint,
                      bool twod) {
  bool non_zero_z = false;
  atomPoint.resize(mol.getNumAtoms());
  // take X,Y,Z coordinates of each atom
  if (0 != mol.getNumConformers())
    for (auto cnfi = mol.beginConformers(); cnfi != mol.endConformers();
         cnfi++) {
      const Conformer &conf = **cnfi;  // mol.getConformer(confId);
      if (twod || conf.is3D()) {
        for (unsigned i = 0; i < mol.getNumAtoms(); i++) {
          atomPoint[i] = conf.getAtomPos(i);
          if (fabs(atomPoint[i].z) >= 1.e-7) non_zero_z = true;
        }
        break;
      }
    }
  if (atomPoint.empty()) {  // compute XYZ
    // TODO:
    // ???? ..........
  }
  return non_zero_z;
}

typedef std::tuple<std::string, int, int, int> NbrData;

bool lessTuple(const NbrData &left, const NbrData &right) {
  if (std::get<0>(left) < std::get<0>(right)) return true;
  if (std::get<0>(left) > std::get<0>(right)) return false;

  if (std::get<1>(left) < std::get<1>(right)) return true;
  if (std::get<1>(left) > std::get<1>(right)) return false;

  if (std::get<2>(left) < std::get<2>(right)) return true;
  if (std::get<2>(left) > std::get<2>(right)) return false;

  if (std::get<3>(left) < std::get<3>(right)) return true;

  return false;
}

std::string LogNeighbourhood(
    const ROMol &mol, unsigned int idx,
    const std::vector<Neighbourhood> &neighbour_array) {
  std::stringstream oss;
  // FIX ME turn into utility func?
  std::string name("");
  mol.getPropIfPresent(common_properties::_Name, name);
  oss << "atom '" << name << "' idx=" << idx << " AA: ";
  const Atom &atm = *mol.getAtomWithIdx(idx);
  oss << atm.getSymbol();

  if (atm.getFormalCharge())
    oss << (atm.getFormalCharge() > 0 ? "+" : "") << atm.getFormalCharge();

  if (atm.getNumRadicalElectrons()) oss << atm.getNumRadicalElectrons();

  // these neighbors should be sorted properly?
  size_t numNbrs = neighbour_array[idx].Atoms.size();
  // CHECK_INVARIANT(numNBrs == neighbour_array[idx].Bonds.size());

  std::vector<NbrData> nbrs;

  for (size_t i = 0; i < numNbrs; ++i) {
    const Bond *bond = mol.getBondWithIdx(neighbour_array[idx].Bonds[i]);
    const Atom &nbr = *mol.getAtomWithIdx(neighbour_array[idx].Atoms[i]);
    nbrs.push_back(NbrData(nbr.getSymbol(), bond->getBondType(),
                           nbr.getFormalCharge(),
                           nbr.getNumRadicalElectrons()));
  }

  std::sort(nbrs.begin(), nbrs.end(), lessTuple);
  for (auto &nbr : nbrs) {
    std::string bs = "";
    switch (std::get<1>(nbr)) {
      case Bond::SINGLE:
        bs = "-";
        break;
      case Bond::DOUBLE:
        bs = "=";
        break;
      case Bond::TRIPLE:
        bs = "#";
        break;
      case Bond::AROMATIC:
        bs = "~";
        break;
    }
    if (bs.size())
      oss << "(" << bs << std::get<0>(nbr);
    else
      oss << "("
          << "?" << (int)std::get<1>(nbr) << "?" << std::get<0>(nbr);
    if (std::get<2>(nbr))
      oss << (std::get<2>(nbr) > 0 ? "+" : "") << std::get<2>(nbr);

    if (std::get<3>(nbr)) oss << std::get<3>(nbr);
    oss << ")";
  }
  return oss.str();
}

}  // namespace StructureCheck
}  // namespace RDKit
