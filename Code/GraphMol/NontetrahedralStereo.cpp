//
//  Copyright (C) 2022 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/types.h>
#include <RDGeneral/Invariant.h>
#include <RDGeneral/RDLog.h>

#include "Chirality.h"

namespace RDKit {

namespace Chirality {
constexpr unsigned char squareplanar_across[4][4] = {
    {4, 4, 4, 4},  //
    {2, 3, 0, 1},  // SP1
    {1, 0, 3, 2},  // SP2
    {3, 2, 1, 0}   // SP3
};

constexpr unsigned char trigonalbipyramidal_across[21][5] = {
    {5, 5, 5, 5, 5},  //
    {4, 5, 5, 5, 0},  // TB1
    {4, 5, 5, 5, 0},  // TB2
    {3, 5, 5, 0, 5},  // TB3
    {3, 5, 5, 0, 5},  // TB4
    {2, 5, 0, 5, 5},  // TB5
    {2, 5, 0, 5, 5},  // TB6
    {1, 0, 5, 5, 5},  // TB7
    {1, 0, 5, 5, 5},  // TB8
    {5, 4, 5, 5, 1},  // TB9
    {5, 3, 5, 1, 5},  // TB10
    {5, 4, 5, 5, 1},  // TB11
    {5, 3, 5, 1, 5},  // TB12
    {5, 2, 1, 5, 5},  // TB13
    {5, 2, 1, 5, 5},  // TB14
    {5, 5, 4, 5, 2},  // TB15
    {5, 5, 3, 2, 5},  // TB16
    {5, 5, 5, 4, 3},  // TB17
    {5, 5, 5, 4, 3},  // TB18
    {5, 5, 3, 2, 5},  // TB19
    {5, 5, 4, 5, 2}   // TB20
};

constexpr unsigned char octahedral_across[31][6] = {
    {6, 6, 6, 6, 6, 6},  //
    {5, 3, 4, 1, 2, 0},  // OH1
    {5, 3, 4, 1, 2, 0},  // OH2
    {4, 3, 5, 1, 0, 2},  // OH3
    {5, 4, 3, 2, 1, 0},  // OH4
    {4, 5, 3, 2, 0, 1},  // OH5
    {3, 4, 5, 0, 1, 2},  // OH6
    {3, 5, 4, 0, 2, 1},  // OH7
    {5, 2, 1, 4, 3, 0},  // OH8
    {4, 2, 1, 5, 0, 3},  // OH9
    {5, 2, 1, 4, 3, 0},  // OH10
    {4, 2, 1, 5, 0, 3},  // OH11
    {3, 2, 1, 0, 5, 4},  // OH12
    {3, 2, 1, 0, 5, 4},  // OH13
    {5, 4, 3, 2, 1, 0},  // OH14
    {4, 5, 3, 2, 0, 1},  // OH15
    {4, 3, 5, 1, 0, 2},  // OH16
    {3, 5, 4, 0, 2, 1},  // OH17
    {3, 4, 5, 0, 1, 2},  // OH18
    {2, 4, 0, 5, 1, 3},  // OH19
    {2, 5, 0, 4, 3, 1},  // OH20
    {2, 3, 0, 1, 5, 4},  // OH21
    {2, 3, 0, 1, 5, 4},  // OH22
    {2, 5, 0, 4, 3, 1},  // OH23
    {2, 4, 0, 5, 1, 3},  // OH24
    {1, 0, 4, 5, 2, 3},  // OH25
    {1, 0, 5, 4, 3, 2},  // OH26
    {1, 0, 3, 2, 5, 4},  // OH27
    {1, 0, 3, 2, 5, 4},  // OH28
    {1, 0, 5, 4, 3, 2},  // OH29
    {1, 0, 4, 5, 2, 3}   // OH30
};

Bond *getChiralAcrossBond(const Atom *cen, const Bond *qry) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(qry, "bad query pointer");
  PRECONDITION(cen->hasOwningMol() && qry->hasOwningMol(), "no owning mol");
  PRECONDITION(&cen->getOwningMol() == &qry->getOwningMol(),
               "cen and query must come from the same molecule");

  Atom::ChiralType tag = cen->getChiralTag();
  unsigned int ref_max = 0;

  switch (tag) {
    case Atom::ChiralType::CHI_SQUAREPLANAR:
      ref_max = 4;
      break;
    case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
      ref_max = 5;
      break;
    case Atom::ChiralType::CHI_OCTAHEDRAL:
      ref_max = 6;
      break;
    default:
      return nullptr;
  }

  unsigned int perm = 0;
  cen->getPropIfPresent(common_properties::_chiralPermutation, perm);
  if (!perm) {
    return nullptr;
  }

  auto &mol = cen->getOwningMol();
  unsigned int count = 0;
  Bond *ref[6];
  int found = -1;

  for (auto bnd : mol.atomBonds(cen)) {
    if (count == ref_max) {
      return nullptr;
    }
    ref[count] = bnd;
    if (bnd == qry) {
      found = count;
    }
    count++;
  }

  if (found >= 0) {
    switch (tag) {
      case Atom::ChiralType::CHI_SQUAREPLANAR:
        if (perm <= 3) {
          found = squareplanar_across[perm][found];
        } else {
          found = 4;
        }
        break;
      case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
        if (perm <= 20) {
          found = trigonalbipyramidal_across[perm][found];
        } else {
          found = 5;
        }
        break;
      case Atom::ChiralType::CHI_OCTAHEDRAL:
        if (perm <= 30) {
          found = octahedral_across[perm][found];
        } else {
          found = 6;
        }
        break;
      default:
        return nullptr;
    }
    if (static_cast<unsigned int>(found) < count) {
      return ref[found];
    }
  }
  return nullptr;
}

/* Alternate wrappers */
Bond *getChiralAcrossBond(const Atom *cen, const Atom *qry) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(qry, "bad query pointer");
  PRECONDITION(cen->hasOwningMol(), "no owning mol");
  auto bnd =
      cen->getOwningMol().getBondBetweenAtoms(cen->getIdx(), qry->getIdx());
  return getChiralAcrossBond(cen, bnd);
}

Atom *getChiralAcrossAtom(const Atom *cen, const Bond *qry) {
  auto bnd = getChiralAcrossBond(cen, qry);
  return bnd ? bnd->getOtherAtom(cen) : nullptr;
}

Atom *getChiralAcrossAtom(const Atom *cen, const Atom *qry) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(qry, "bad query pointer");
  PRECONDITION(cen->hasOwningMol(), "no owning mol");
  auto bnd =
      cen->getOwningMol().getBondBetweenAtoms(cen->getIdx(), qry->getIdx());
  bnd = getChiralAcrossBond(cen, bnd);
  if (bnd) {
    return bnd->getOtherAtom(cen);
  } else {
    return nullptr;
  }
}
}  // namespace Chirality
}  // namespace RDKit
