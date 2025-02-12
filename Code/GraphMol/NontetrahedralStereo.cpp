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

constexpr unsigned char swap_squareplanar_table[4][6] = {
    // 0  0  0  1  1  2
    // 1  2  3  2  3  3
    {0, 0, 0, 0, 0, 0},
    {3, 1, 2, 2, 1, 3},  // SP1
    {2, 3, 1, 1, 3, 2},  // SP2
    {1, 2, 3, 3, 2, 1}   // SP3
};

constexpr unsigned char swap_trigonalbipyramidal_table[21][10] = {
    // 0   0   0   0   1   1   1   2   2   3
    // 1   2   3   4   2   3   4   3   4   4
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {9, 20, 17, 2, 2, 2, 7, 2, 6, 3},        // TB1
    {11, 15, 18, 1, 1, 1, 8, 1, 5, 4},       // TB2
    {10, 19, 4, 18, 4, 8, 4, 5, 4, 1},       // TB3
    {12, 16, 3, 17, 3, 7, 3, 6, 3, 2},       // TB4
    {13, 6, 16, 20, 7, 6, 6, 3, 2, 6},       // TB5
    {14, 5, 19, 15, 8, 5, 5, 4, 1, 5},       // TB6
    {8, 14, 10, 11, 5, 4, 1, 8, 8, 8},       // TB7
    {7, 13, 12, 9, 6, 3, 2, 7, 7, 7},        // TB8
    {1, 11, 11, 8, 15, 18, 11, 11, 14, 10},  // TB9
    {3, 12, 7, 12, 16, 12, 17, 13, 12, 9},   // TB10
    {2, 9, 9, 7, 20, 17, 9, 9, 13, 12},      // TB11
    {4, 10, 8, 10, 19, 10, 18, 14, 10, 11},  // TB12
    {5, 8, 14, 14, 14, 19, 15, 10, 11, 14},  // TB13
    {6, 7, 13, 13, 13, 16, 20, 12, 9, 13},   // TB14
    {20, 2, 20, 6, 9, 20, 13, 17, 20, 16},   // TB15
    {19, 4, 5, 19, 10, 14, 19, 19, 18, 15},  // TB16
    {18, 18, 1, 4, 18, 11, 10, 15, 19, 18},  // TB17
    {17, 17, 2, 3, 17, 9, 12, 20, 16, 17},   // TB18
    {16, 3, 6, 16, 12, 13, 16, 16, 17, 20},  // TB19
    {15, 1, 15, 5, 11, 15, 14, 18, 15, 19}   // TB20
};

constexpr unsigned char swap_octahedral_table[31][15] = {
    // 0   0   0   0   0   1   1   1   1   2   2   2   3   3   4
    // 1   2   3   4   5   2   3   4   5   3   4   5   4   5   5
    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {17, 16, 30, 21, 2, 14, 2, 10, 25, 8, 2, 22, 4, 7, 3},       // OH1
    {7, 3, 25, 22, 1, 4, 1, 8, 30, 10, 1, 21, 14, 17, 16},       // OH2
    {18, 2, 29, 16, 22, 15, 16, 26, 11, 9, 21, 16, 6, 5, 1},     // OH3
    {15, 18, 19, 28, 14, 2, 8, 14, 27, 14, 10, 24, 1, 6, 5},     // OH4
    {14, 17, 20, 15, 27, 16, 9, 28, 15, 15, 23, 11, 7, 3, 4},    // OH5
    {16, 14, 18, 26, 24, 17, 29, 18, 13, 19, 12, 18, 3, 4, 7},   // OH6
    {2, 15, 17, 23, 25, 18, 30, 12, 17, 20, 17, 13, 5, 1, 6},    // OH7
    {23, 26, 11, 12, 10, 10, 4, 2, 29, 1, 14, 20, 10, 13, 9},    // OH8
    {24, 25, 10, 11, 13, 11, 5, 30, 16, 3, 19, 15, 12, 11, 8},   // OH9
    {20, 29, 9, 13, 8, 8, 14, 1, 26, 2, 4, 23, 8, 12, 11},       // OH10
    {19, 30, 8, 9, 12, 9, 15, 25, 3, 16, 24, 5, 13, 9, 10},      // OH11
    {22, 27, 13, 8, 11, 13, 28, 7, 18, 21, 6, 17, 9, 10, 13},    // OH12
    {21, 28, 12, 10, 9, 12, 27, 17, 6, 22, 18, 7, 11, 8, 12},    // OH13
    {5, 6, 24, 27, 4, 1, 10, 4, 28, 4, 8, 19, 2, 18, 15},        // OH14
    {4, 7, 23, 5, 28, 3, 11, 27, 5, 5, 20, 9, 17, 16, 14},       // OH15
    {6, 1, 26, 3, 21, 5, 3, 29, 9, 11, 22, 3, 18, 15, 2},        // OH16
    {1, 5, 7, 20, 30, 6, 25, 13, 7, 23, 7, 12, 15, 2, 18},       // OH17
    {3, 4, 6, 29, 19, 7, 26, 6, 12, 24, 13, 6, 16, 14, 17},      // OH18
    {11, 24, 4, 30, 18, 25, 23, 24, 22, 6, 9, 14, 21, 24, 20},   // OH19
    {10, 23, 5, 17, 29, 26, 24, 21, 23, 7, 15, 8, 23, 22, 19},   // OH20
    {13, 22, 28, 1, 16, 27, 22, 20, 24, 12, 3, 2, 19, 23, 22},   // OH21
    {12, 21, 27, 2, 3, 28, 21, 23, 19, 13, 16, 1, 24, 20, 21},   // OH22
    {8, 20, 15, 7, 26, 29, 19, 22, 20, 17, 5, 10, 20, 21, 24},   // OH23
    {9, 19, 14, 25, 6, 30, 20, 19, 21, 18, 11, 4, 22, 19, 23},   // OH24
    {30, 9, 2, 24, 7, 19, 17, 11, 1, 29, 30, 28, 27, 30, 26},    // OH25
    {29, 8, 16, 6, 23, 20, 18, 3, 10, 30, 27, 29, 29, 28, 25},   // OH26
    {28, 12, 22, 14, 5, 21, 13, 15, 4, 28, 26, 30, 25, 29, 28},  // OH27
    {27, 13, 21, 4, 15, 22, 12, 5, 14, 27, 29, 25, 30, 26, 27},  // OH28
    {26, 10, 3, 18, 20, 23, 6, 16, 8, 25, 28, 26, 26, 27, 30},   // OH29
    {25, 11, 1, 19, 17, 24, 7, 9, 2, 26, 25, 27, 28, 25, 29}     // OH30
};

static unsigned int swap_squareplanar(unsigned int perm, unsigned int x,
                                      unsigned int y) {
  constexpr unsigned int offset[3] = {0, 2, 3};
  unsigned int swapidx;
  if (x == y) {
    return perm;
  }
  if (x < y) {
    if (y > 3) {
      return 0;
    }
    swapidx = offset[x] + (y - 1);
  } else /* x > y */ {
    if (x > 3) {
      return 0;
    }
    swapidx = offset[y] + (x - 1);
  }
  return perm < 4 ? swap_squareplanar_table[perm][swapidx] : 0;
}

static unsigned int swap_trigonalbipyramidal(unsigned int perm, unsigned int x,
                                             unsigned int y) {
  constexpr unsigned int offset[4] = {0, 3, 5, 6};
  unsigned int swapidx;
  if (x == y) {
    return perm;
  }
  if (x < y) {
    if (y > 4) {
      return 0;
    }
    swapidx = offset[x] + (y - 1);
  } else /* x > y */ {
    if (x > 4) {
      return 0;
    }
    swapidx = offset[y] + (x - 1);
  }
  return perm < 21 ? swap_trigonalbipyramidal_table[perm][swapidx] : 0;
}

static unsigned int swap_octahedral(unsigned int perm, unsigned int x,
                                    unsigned int y) {
  constexpr unsigned int offset[5] = {0, 4, 7, 9, 10};
  unsigned int swapidx;
  if (x == y) {
    return perm;
  }
  if (x < y) {
    if (y > 5) {
      return 0;
    }
    swapidx = offset[x] + (y - 1);
  } else /* x > y */ {
    if (x > 5) {
      return 0;
    }
    swapidx = offset[y] + (x - 1);
  }
  return perm < 31 ? swap_octahedral_table[perm][swapidx] : 0;
}

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

constexpr unsigned char trigonalbipyramidal_axial[21][2] = {
    {5, 5},  //
    {0, 4},  // TB1
    {0, 4},  // TB2
    {0, 3},  // TB3
    {0, 3},  // TB4
    {0, 2},  // TB5
    {0, 2},  // TB6
    {0, 1},  // TB7
    {0, 1},  // TB8
    {1, 4},  // TB9
    {1, 4},  // TB10
    {1, 3},  // TB11
    {1, 3},  // TB12
    {1, 2},  // TB13
    {1, 2},  // TB14
    {2, 4},  // TB15
    {2, 3},  // TB16
    {3, 4},  // TB17
    {3, 4},  // TB18
    {2, 3},  // TB19
    {2, 4}   // TB20
};

int isTrigonalBipyramidalAxialBond(const Atom *cen, const Bond *qry) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(qry, "bad query pointer");
  PRECONDITION(cen->hasOwningMol() && qry->hasOwningMol(), "no owning mol");
  PRECONDITION(&cen->getOwningMol() == &qry->getOwningMol(),
               "center and query must come from the same molecule");
  if (cen->getDegree() > 5 ||
      cen->getChiralTag() != Atom::CHI_TRIGONALBIPYRAMIDAL) {
    return false;
  }
  unsigned int perm = 0;
  cen->getPropIfPresent(RDKit::common_properties::_chiralPermutation, perm);
  if (perm == 0 || perm > 20) return 0;

  unsigned int count = 0;
  for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
    if (bnd == qry) {
      if (count == trigonalbipyramidal_axial[perm][0]) {
        return 1;
      }
      if (count == trigonalbipyramidal_axial[perm][1]) {
        return -1;
      }
      return 0;
    }
    count++;
  }
  return 0;
}

int isTrigonalBipyramidalAxialAtom(const Atom *cen, const Atom *qry) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(qry, "bad query pointer");
  PRECONDITION(cen->hasOwningMol(), "no owning mol");
  auto bnd =
      cen->getOwningMol().getBondBetweenAtoms(cen->getIdx(), qry->getIdx());
  if (!bnd) {
    return 0;
  }
  return isTrigonalBipyramidalAxialBond(cen, bnd);
}

unsigned int getMaxNbors(const Atom::ChiralType tag) {
  switch (tag) {
    case Atom::CHI_TETRAHEDRAL_CW:
    case Atom::CHI_TETRAHEDRAL_CCW:
    case Atom::CHI_TETRAHEDRAL:  // fall through
      return 4;
    case Atom::CHI_ALLENE:
      return 2;  // not used other than SMI/SMA parsers?
    case Atom::CHI_SQUAREPLANAR:
      return 4;
    case Atom::CHI_TRIGONALBIPYRAMIDAL:
      return 5;
    case Atom::CHI_OCTAHEDRAL:
      return 6;
    default:
      BOOST_LOG(rdWarningLog)
          << "Warning: unexpected chiral tag getMaxNbors(): " << tag
          << std::endl;
      return 0;
  }
}

Bond *getChiralAcrossBond(const Atom *cen, const Bond *qry) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(qry, "bad query pointer");
  PRECONDITION(cen->hasOwningMol() && qry->hasOwningMol(), "no owning mol");
  PRECONDITION(&cen->getOwningMol() == &qry->getOwningMol(),
               "center and query must come from the same molecule");

  Atom::ChiralType tag = cen->getChiralTag();
  unsigned int perm = 0;
  cen->getPropIfPresent(common_properties::_chiralPermutation, perm);
  if (!perm) {
    return nullptr;
  }

  auto &mol = cen->getOwningMol();
  unsigned int count = 0;
  Bond *ref[6];
  int found = -1;

  unsigned int ref_max = getMaxNbors(tag);
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
  if (!bnd) {
    return nullptr;
  }
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
  if (!bnd) {
    return nullptr;
  }
  bnd = getChiralAcrossBond(cen, bnd);
  return bnd ? bnd->getOtherAtom(cen) : nullptr;
}

double getIdealAngleBetweenLigands(const Atom *cen, const Atom *lig1,
                                   const Atom *lig2) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(lig1 && lig2, "bad ligand pointer");
  PRECONDITION(
      cen->hasOwningMol() && lig1->hasOwningMol() && lig2->hasOwningMol(),
      "no owning mol");
  PRECONDITION(&cen->getOwningMol() == &lig1->getOwningMol() &&
                   &cen->getOwningMol() == &lig2->getOwningMol(),
               "center and ligands must come from the same molecule");
  auto tag = cen->getChiralTag();
  switch (tag) {
    case Atom::ChiralType::CHI_SQUAREPLANAR:
    case Atom::ChiralType::CHI_OCTAHEDRAL:
      return getChiralAcrossAtom(cen, lig1) == lig2 ? 180 : 90;
      break;
    case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
      if (getChiralAcrossAtom(cen, lig1) == lig2) {
        // both are axial
        return 180;
      } else if (isTrigonalBipyramidalAxialAtom(cen, lig1) ||
                 isTrigonalBipyramidalAxialAtom(cen, lig2)) {
        // one is axial, the other equatorial
        return 90;
      } else {
        // both axial
        return 120;
      }
      break;
    default:
      return 0;
  }
}

Bond *getTrigonalBipyramidalAxialBond(const Atom *cen, int axial) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(cen->hasOwningMol(), "no owning mol");
  if (cen->getChiralTag() != RDKit::Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL ||
      cen->getDegree() > 5)
    return nullptr;

  unsigned int perm = 0;
  cen->getPropIfPresent(RDKit::common_properties::_chiralPermutation, perm);
  if (perm == 0 || perm > 20) return nullptr;

  unsigned int idx = (axial != -1) ? trigonalbipyramidal_axial[perm][0]
                                   : trigonalbipyramidal_axial[perm][1];

  unsigned int count = 0;
  for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
    if (count == idx) return bnd;
    count++;
  }
  return nullptr;
}

Atom *getTrigonalBipyramidalAxialAtom(const Atom *cen, int axial) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(cen->hasOwningMol(), "no owning mol");
  auto bnd = getTrigonalBipyramidalAxialBond(cen, axial);
  return bnd ? bnd->getOtherAtom(cen) : nullptr;
}
bool hasNonTetrahedralStereo(const Atom *cen) {
  PRECONDITION(cen, "bad center pointer");
  if (!cen->hasOwningMol()) {
    return false;
  }
  auto tag = cen->getChiralTag();
  return tag == Atom::ChiralType::CHI_SQUAREPLANAR ||
         tag == Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL ||
         tag == Atom::ChiralType::CHI_OCTAHEDRAL;
}

unsigned int getChiralPermutation(const Atom *cen, const INT_LIST &probe,
                                  bool inverse) {
  PRECONDITION(cen, "bad center pointer");
  PRECONDITION(cen->hasOwningMol(), "no owning mol");

  int perm;
  if (!cen->getPropIfPresent(common_properties::_chiralPermutation, perm) ||
      perm <= 0) {
    return 0;
  }

  decltype(&swap_octahedral) swap_func = nullptr;
  switch (cen->getChiralTag()) {
    case Atom::ChiralType::CHI_OCTAHEDRAL:
      if (probe.size() > 6) {
        return 0;
      }
      swap_func = swap_octahedral;
      break;
    case Atom::ChiralType::CHI_TRIGONALBIPYRAMIDAL:
      if (probe.size() > 5) {
        return 0;
      }
      swap_func = swap_trigonalbipyramidal;
      break;
    case Atom::ChiralType::CHI_SQUAREPLANAR:
      if (probe.size() > 4) {
        return 0;
      }
      swap_func = swap_squareplanar;
      break;
    default:
      break;
  }
  if (!swap_func) {
    return 0;
  }
  unsigned int nbrIdx = 0;
  std::vector<int> order(cen->getOwningMol().getNumBonds(), -1);
  for (const auto bnd : cen->getOwningMol().atomBonds(cen)) {
    order[bnd->getIdx()] = nbrIdx++;
  }

  // nbrPerm maps original index to array position
  std::vector<unsigned int> nbrPerm(nbrIdx);
  std::iota(nbrPerm.begin(), nbrPerm.end(), 0);
  std::vector<unsigned int> probePerm(probe.size());
  nbrIdx = 0;
  for (auto v : probe) {
    probePerm[nbrIdx++] = v < 0 ? -1 : order[v];
  }

  // Missing (implicit) neighbors are at the end when in storage order
  if (nbrPerm.size() < nbrIdx)
    nbrPerm.insert(nbrPerm.end(), nbrIdx - nbrPerm.size(), -1);

  CHECK_INVARIANT(nbrPerm.size() == probePerm.size(),
                  "probe vector size does not match");

  if (inverse) {
    std::swap(nbrPerm, probePerm);
  }

  boost::dynamic_bitset<> swapped(probe.size());
  for (unsigned int i = 0; i < probePerm.size() - 1; ++i) {
    auto pval = probePerm[i];
    if (nbrPerm[i] == pval) {
      continue;
    }
    auto tgt = std::find(nbrPerm.begin() + i, nbrPerm.end(), pval);
    TEST_ASSERT(tgt != nbrPerm.end());
    perm = swap_func(perm, i, tgt - nbrPerm.begin());
    std::swap(*tgt, nbrPerm[i]);
  }

  return perm;
}

}  // namespace Chirality
}  // namespace RDKit
