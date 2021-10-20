//
//  Copyright (C) 2016 Novartis Institutes for BioMedical Research
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

typedef boost::tuple<std::string, int, int, int> NbrData;

bool lessTuple(const NbrData &left, const NbrData &right) {
  if (left.get<0>() < right.get<0>()) return true;
  if (left.get<0>() > right.get<0>()) return false;

  if (left.get<1>() < right.get<1>()) return true;
  if (left.get<1>() > right.get<1>()) return false;

  if (left.get<2>() < right.get<2>()) return true;
  if (left.get<2>() > right.get<2>()) return false;

  if (left.get<3>() < right.get<3>()) return true;

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
    switch (nbr.get<1>()) {
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
      oss << "(" << bs << nbr.get<0>();
    else
      oss << "("
          << "?" << (int)nbr.get<1>() << "?" << nbr.get<0>();
    if (nbr.get<2>()) oss << (nbr.get<2>() > 0 ? "+" : "") << nbr.get<2>();

    if (nbr.get<3>()) oss << nbr.get<3>();
    oss << ")";
  }
  return oss.str();
}

}  // namespace StructureCheck
}  // namespace RDKit
