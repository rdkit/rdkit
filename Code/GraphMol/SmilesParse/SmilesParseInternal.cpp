//
//  Copyright (C) 2025 NVIDIA Corporation & Affiliates and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SmilesParseInternal.h"
#include "SmilesParse.h"

#include <GraphMol/RDMol.h>

#include <vector>
#include <assert.h>
#include <numeric>
#include <memory>

using namespace RDKit;
using BondEnums::BondType;

constexpr static uint32_t INVALID_OPEN_BOND =
    SmilesParseTemp::INVALID_OPEN_BOND;

namespace {

constexpr size_t elementSymbolLookupStart[27] = {
    0,    // Ac Ag Al Am Ar As At Au             ( 8)
    8,    // B  Ba Be Bh Bi Bk Br                ( 7)
    15,   // C  Ca Cd Ce Cf Cl Cm Cn Co Cr Cs Cu (12)
    27,   // Db Ds Dy                            ( 3)
    30,   // Er Es Eu                            ( 3)
    33,   // F  Fe Fl Fm Fr                      ( 5)
    38,   // Ga Gd Ge                            ( 3)
    41,   // H  He Hf Hg Ho Hs                   ( 6)
    47,   // I  In Ir                            ( 3)
    50,   //                                     ( 0)
    50,   // K  Kr                               ( 2)
    52,   // La Li Lr Lu Lv                      ( 5)
    57,   // Mc Md Mg Mn Mo Mt                   ( 6)
    63,   // N  Na Nb Nd Ne Nh Ni No Np          ( 9)
    72,   // O  Og Os                            ( 3)
    75,   // P  Pa Pb Pd Pm Po Pr Pt Pu          ( 9)
    84,   //                                     ( 0)
    84,   // Ra Rb Re Rf Rg Rh Rn Ru             ( 8)
    92,   // S  Sb Sc Se Sg Si Sm Sn Sr          ( 9)
    101,  // Ta Tb Tc Te Th Ti Tl Tm Ts          ( 9)
    110,  // U                                   ( 1)
    111,  // V                                   ( 1)
    112,  // W                                   ( 1)
    113,  // Xe                                  ( 1)
    114,  // Y  Yb                               ( 2)
    116,  // Zn Zr                               ( 2)
    118};

// These are in *alphabetical* order by element symbol
constexpr char elementSymbolSecond[118] = {
    'c', 'g', 'l', 'm', 'r', 's', 't', 'u', 0,   'a', 'e', 'h', 'i',
    'k', 'r', 0,   'a', 'd', 'e', 'f', 'l', 'm', 'n', 'o', 'r', 's',
    'u', 'b', 's', 'y', 'r', 's', 'u', 0,   'e', 'l', 'm', 'r', 'a',
    'd', 'e', 0,   'e', 'f', 'g', 'o', 's', 0,   'n', 'r',

    0,   'r', 'a', 'i', 'r', 'u', 'v', 'c', 'd', 'g', 'n', 'o', 't',
    0,   'a', 'b', 'd', 'e', 'h', 'i', 'o', 'p', 0,   'g', 's', 0,
    'a', 'b', 'd', 'm', 'o', 'r', 't', 'u',

    'a', 'b', 'e', 'f', 'g', 'h', 'n', 'u', 0,   'b', 'c', 'e', 'g',
    'i', 'm', 'n', 'r', 'a', 'b', 'c', 'e', 'h', 'i', 'l', 'm', 's',
    0,   0,   0,   'e', 0,   'b', 'n', 'r'};

// These are in *alphabetical* order by element symbol, same order as
// elementSymbolSecond
constexpr int8_t elementAtomicNum[118] = {
    89, 47,  13,  95,  18,  33,  85,  79,  5,   56,  4,   107, 83,
    97, 35,  6,   20,  48,  58,  98,  17,  96,  112, 27,  24,  55,
    29, 105, 110, 66,  68,  99,  63,  9,   26,  114, 100, 87,  31,
    64, 32,  1,   2,   72,  80,  67,  108, 53,  49,  77,

    19, 36,  57,  3,   103, 71,  116, 115, 101, 12,  25,  42,  109,
    7,  11,  41,  60,  10,  113, 28,  102, 93,  8,   118, 76,  15,
    91, 82,  46,  61,  84,  59,  78,  94,

    88, 37,  75,  104, 111, 45,  86,  44,  16,  51,  21,  34,  106,
    14, 62,  50,  38,  73,  65,  43,  52,  90,  22,  81,  69,  117,
    92, 23,  74,  54,  39,  70,  30,  40};

// Fill vectors representing the molecule's graph in CSR format.
void toCSR(const size_t numAtoms, const std::vector<BondData>& bonds, RDMol& mol) {
  const size_t numBonds = bonds.size();

  mol.getAtomBondStartsVector().resize(numAtoms + 1);
  uint32_t* neighbourStarts = mol.getAtomBondStartsVector().data();
  for (size_t i = 0; i <= numAtoms; ++i) {
    neighbourStarts[i] = 0;
  }

  // First, get atom neighbour counts
  for (size_t i = 0; i < numBonds; ++i) {
    uint32_t a = bonds[i].getBeginAtomIdx();
    uint32_t b = bonds[i].getEndAtomIdx();
    // NOTE: +1 is because first entry will stay zero.
    ++neighbourStarts[a + 1];
    ++neighbourStarts[b + 1];
  }

  // Find the starts by partial-summing the neighbour counts.
  // NOTE: +1 is because first entry will stay zero.
  std::partial_sum(neighbourStarts + 1, neighbourStarts + 1 + numAtoms,
                   neighbourStarts + 1);

  // Fill in the neighbours and bondIndices arrays.
  mol.getOtherAtomIndicesVector().resize(2 * numBonds);
  uint32_t* neighbours = mol.getOtherAtomIndicesVector().data();
  mol.getBondDataIndicesVector().resize(2 * numBonds);
  uint32_t* bondIndices = mol.getBondDataIndicesVector().data();
  for (uint32_t i = 0; i < numBonds; ++i) {
    uint32_t a = bonds[i].getBeginAtomIdx();
    uint32_t b = bonds[i].getEndAtomIdx();

    uint32_t ai = neighbourStarts[a];
    neighbours[ai] = b;
    bondIndices[ai] = i;
    ++neighbourStarts[a];

    uint32_t bi = neighbourStarts[b];
    neighbours[bi] = a;
    bondIndices[bi] = i;
    ++neighbourStarts[b];
  }

  // Shift neighbourStarts forward one after incrementing it.
  uint32_t previous = 0;
  for (size_t i = 0; i < numAtoms; ++i) {
    uint32_t next = neighbourStarts[i];
    neighbourStarts[i] = previous;
    previous = next;
  }
  assert(neighbourStarts[numAtoms] == previous);
}

// text starts out pointing to it
inline int parseInt(const char*& text) {
  char c = *text;
  int v = c - '0';
  ++text;
  c = *text;
  while (c >= '0' && c <= '9') {
    v = v * 10 + (c - '0');
    ++text;
    c = *text;
  }
  return v;
}

static bool parseBracketAtom(const char*& text, AtomData& atom,
                             uint32_t atomIndex, SmilesParseTemp& temp) {
  char c = *text;
  if (c >= '0' && c <= '9') {
    // Isotope
    atom.setIsotope(parseInt(text));
    c = *text;
  }

  if (c >= 'A' && c <= 'Z') {
    // Element symbol
    const size_t firstLetterIndex = c - 'A';
    ++text;
    c = *text;
    char secondLetter = (c >= 'a' && c <= 'z') ? c : 0;
    const size_t firstLetterIndexBegin =
        elementSymbolLookupStart[firstLetterIndex];
    const size_t firstLetterIndexEnd =
        elementSymbolLookupStart[firstLetterIndex + 1];
    size_t elementAlphabeticalIndex = firstLetterIndexBegin;
    while (elementAlphabeticalIndex < firstLetterIndexEnd &&
           secondLetter != elementSymbolSecond[elementAlphabeticalIndex]) {
      ++elementAlphabeticalIndex;
    }
    if (elementAlphabeticalIndex == firstLetterIndexEnd) {
      return false;
    }
    if (secondLetter != 0) {
      ++text;
      c = *text;
    }
    atom.setAtomicNum(elementAtomicNum[elementAlphabeticalIndex]);
  } else if (c >= 'a' && c <= 'z') {
    ++text;
    // Aromatic element symbol
    // More elements are allowed to be aromatic when bracketed than unbracketed.
    uint32_t atomicNum;
    switch (c) {
      case 'a':
        c = *text;
        if (c != 's') {
          return false;
        }
        ++text;
        atomicNum = 33;  // As
        break;
      case 'b':
        atomicNum = 5;  // B
        break;
      case 'c':
        atomicNum = 6;  // C
        break;
      case 'n':
        atomicNum = 7;  // N
        break;
      case 'o':
        atomicNum = 8;  // O
        break;
      case 'p':
        atomicNum = 15;  // P
        break;
      case 's':
        c = *text;
        if (c != 'e') {
          atomicNum = 16;  // S
        } else {
          ++text;
          atomicNum = 34;  // Se
        }
        break;
      default:
        return false;
    }
    atom.setAtomicNum(atomicNum);
    atom.setIsAromatic(true);
    c = *text;
  } else if (c == '*') {
    // Unknown element
    atom.setAtomicNum(0);
    ++text;
    c = *text;
  } else {
    return false;
  }

  if (c == '@') {
    // Chirality
    ++text;
    c = *text;
    if ((c == 'A') || (c == 'O') || (c == 'S') || (c == 'T')) {
      ++text;
      char c2 = *text;
      if (c == 'A') {
        // Only AL is valid "Allenal" (1 or 2)
        if (c2 != 'L') {
          return false;
        }
        atom.setChiralTag(AtomEnums::ChiralType::CHI_ALLENE);
      } else if (c2 == 'O') {
        // Only OH is valid "Octahedral" (up to 2 digits)
        if (c2 != 'H') {
          return false;
        }
        atom.setChiralTag(AtomEnums::ChiralType::CHI_OCTAHEDRAL);
      } else if (c == 'S') {
        // Only SP is valid "Square Planar" (1, 2, or 3)
        if (c2 != 'P') {
          return false;
        }
        atom.setChiralTag(AtomEnums::ChiralType::CHI_SQUAREPLANAR);
      } else if (c == 'T') {
        if (c2 == 'B') {
          // TB "Trigonal Bipyramidal" (up to 2 digits)
          atom.setChiralTag(AtomEnums::ChiralType::CHI_TRIGONALBIPYRAMIDAL);
        } else if (c2 == 'H') {
          // TH "Tetrahedral" (1 or 2)
          atom.setChiralTag(AtomEnums::ChiralType::CHI_TETRAHEDRAL);
        } else {
          return false;
        }
      }
      ++text;
      c = *text;
      if (c < '0' || c > '9') {
        return false;
      }
      int number = (c - '0');
      ++text;
      c = *text;
      if (c >= '0' && c <= '9') {
        number = number * 10 + (c - '0');
        ++text;
        c = *text;
      }
      if (temp.chiralityPermutations.size() <= atomIndex) {
        temp.chiralityPermutations.resize(atomIndex + 1, 0);
      }
      temp.chiralityPermutations[atomIndex] = number;
    } else {
      if (c == '@') {
        ++text;
        c = *text;
        // chirality permutation number 2;
        atom.setChiralTag(AtomEnums::ChiralType::CHI_TETRAHEDRAL_CW);
      } else {
        // chirality permutation number 1;
        atom.setChiralTag(AtomEnums::ChiralType::CHI_TETRAHEDRAL_CCW);
      }
    }
  }

  if (c == 'H') {
    // Explicit hydrogen count
    ++text;
    c = *text;
    if (c >= '0' && c <= '9') {
      atom.setNumExplicitHs(c - '0');
      ++text;
      c = *text;
    } else {
      atom.setNumExplicitHs(1);
    }
  } else {
    // Bracketed atoms need explicit hydrogen counts, else no hydrogens.
    atom.setNumExplicitHs(0);
  }
  atom.setNoImplicit(true);

  if (c == '+' || c == '-') {
    // Charge
    const bool isNegative = (c == '-');
    const char sign = c;
    ++text;
    c = *text;
    int count = 1;
    if (c == sign) {
      count = 2;
      ++text;
      c = *text;
    } else if (c >= '0' && c <= '9') {
      count = (c - '0');
      ++text;
      c = *text;
      // Limit of 2 digits for charge
      if (c >= '0' && c <= '9') {
        count = count * 10 + (c - '0');
        ++text;
        c = *text;
      }
    }
    atom.setFormalCharge(isNegative ? -count : count);
  }

  if (c == ':') {
    // Class
    ++text;
    c = *text;
    if (!(c >= '0' && c <= '9')) {
      return false;
    }
    if (temp.atomMapNumbers.size() <= atomIndex) {
      temp.atomMapNumbers.resize(atomIndex + 1, -1);
    }
    temp.atomMapNumbers[atomIndex] = parseInt(text);
    //atom.atomClass = parseInt(text);
    c = *text;
  }

  if (c == ']') {
    ++text;
    return true;
  }
  return false;
}

static bool parseAtom(const char*& text, AtomData& atom, atomindex_t atomIndex, SmilesParseTemp &temp) {
  char c = *text;
  if (c == '[') {
    ++text;
    return parseBracketAtom(text, atom, atomIndex, temp);
  }

  ++text;
  bool isAromatic = (c & 0x20) != 0;
  c = c & ~0x20;
  int atomicNum;
  switch (c) {
    case 'B':
      if (!isAromatic && *text == 'r') {
        ++text;
        atomicNum = 35;  // Br
      } else {
        atomicNum = 5;  // B
      }
      break;
    case 'C':
      if (!isAromatic && *text == 'l') {
        ++text;
        atomicNum = 17;  // Cl
      } else {
        atomicNum = 6;  // C
      }
      break;
    case 'F':
      if (isAromatic) {
        return false;
      }
      atomicNum = 9;
      break;
    case 'I':
      if (isAromatic) {
        return false;
      }
      atomicNum = 53;
      break;
    case 'N':
      atomicNum = 7;
      break;
    case 'O':
      atomicNum = 8;
      break;
    case 'P':
      atomicNum = 15;
      break;
    case 'S':
      atomicNum = 16;
      break;
    case '*':
      atomicNum = 0;
      // '*' could be aromatic, but that needs to be determined later.
      isAromatic = false;
      // No implicit hydrogens on unknown element.
      atom.setNumExplicitHs(0);
      atom.setNoImplicit(true);
      break;
    default:
      return false;
  }
  atom.setAtomicNum(atomicNum);
  atom.setIsAromatic(isAromatic);
  return true;
}

inline BondType bondTypeFromChar(char c) {
  switch (c) {
    case '.':
      return BondType::ZERO;
    case '-':
    case '/':
    case '\\':
      return BondType::SINGLE;
    case '=':
      return BondType::DOUBLE;
    case '#':
      return BondType::TRIPLE;
    case '$':
      return BondType::QUADRUPLE;
    case ':':
      return BondType::AROMATIC;
  }
  return BondType::OTHER;
}

static void handleAnyRingBonds(const char*& text, uint32_t recentAtomIndex,
                               std::vector<BondData>& bonds,
                               SmilesParseTemp& temp) {
  const char* localText = text;
  char c = *localText;
  while ((c < 'A' || c == '\\') && (c != '(') && (c != ')')) {
    BondType type = bondTypeFromChar(c);
    if (type == BondType::ZERO) {
      // Dot isn't allowed for ring bonds.
      return;
    }
    if (type != BondType::OTHER) {
      ++localText;
      c = *localText;
    } else {
      type = BondType::UNSPECIFIED;
    }
    bool twoDigits = (c == '%');
    if (twoDigits) {
      ++localText;
      c = *localText;
    }
    if (c < '0' || c > '9') {
      // Possibly a regular bond, not a ring bond.
      return;
    }
    size_t index = (c - '0');
    ++localText;
    c = *localText;
    if (twoDigits) {
      if (c < '0' || c > '9') {
        return;
      }
      index = index * 10 + (c - '0');
      ++localText;
      c = *localText;
    }
    if (index >= temp.openBondTable.size()) {
      temp.openBondTable.resize(index + 1, INVALID_OPEN_BOND);
    }
    uint32_t existingAtomIndex = temp.openBondTable[index];
    if (existingAtomIndex == INVALID_OPEN_BOND) {
      // Open a new bond
      temp.openBondTable[index] = recentAtomIndex;
      ++temp.numOpenBonds;
    } else {
      // Complete an open bond
      assert(temp.numOpenBonds > 0);
      assert(existingAtomIndex < recentAtomIndex);
      BondData bond;
      bond.setBeginAtomIdx(existingAtomIndex);
      bond.setEndAtomIdx(recentAtomIndex);
      bond.setBondType(type);
      bonds.push_back(std::move(bond));
      //++numNeighbors[existingAtomIndex];
      //++numNeighbors[recentAtomIndex];
      temp.openBondTable[index] = INVALID_OPEN_BOND;
      --temp.numOpenBonds;
    }

    // Done ring bond, so update text pointer in caller.
    text = localText;
  }

  // Not a bond type or '%' or number, so return
}

}  // namespace

namespace SmilesParseInternal {

bool parseSmiles(const char* text, RDMol &mol, SmilesParseTemp& temp) {
  std::vector<AtomData> &atoms = mol.getAtomDataVector();
  std::vector<BondData> &bonds = mol.getBondDataVector();

  atoms.resize(0);
  bonds.resize(0);
  assert(temp.numOpenBonds == 0);
  temp.atomMapNumbers.resize(0);
  temp.chiralityPermutations.resize(0);

  uint32_t prevAtomIndex = 0;
  bool isBondAllowed = false;
  bool isBranchAllowed = false;
  char c = *text;
  while (c > ' ') {
    // If branching, the bond is after the bracket, not before,
    // so process the bracket first.
    if (isBranchAllowed && (c == '(')) {
      ++text;
      c = *text;
      temp.branchStack.push_back(prevAtomIndex);
      isBondAllowed = true;
      isBranchAllowed = false;
      continue;
    }

    const uint32_t nextAtomIndex = uint32_t(atoms.size());
    if (isBondAllowed) {
      BondData bond;
      bond.setBeginAtomIdx(prevAtomIndex);
      bond.setEndAtomIdx(nextAtomIndex);
      if (c < 'A') {
        auto bondType = bondTypeFromChar(c);
        if (bondType == BondType::OTHER) {
          temp.clearForError();
          return false;
        }
        bond.setBondType(bondType);
        ++text;
        c = *text;
      } else if (c == '\\') {
        bond.setBondType(BondType::SINGLE);
        ++text;
        c = *text;
      } else {
        bond.setBondType(BondType::UNSPECIFIED);
      }
      if (bond.getBondType() != BondType::ZERO) {
        bonds.push_back(bond);
      }
    }

    AtomData atom;
    bool success = parseAtom(text, atom, nextAtomIndex, temp);
    if (!success) {
      temp.clearForError();
      return false;
    }
    atoms.push_back(atom);
    prevAtomIndex = nextAtomIndex;
    handleAnyRingBonds(text, prevAtomIndex, bonds, temp);

    c = *text;
    while (c == ')') {
      if (temp.branchStack.size() == 0) {
        temp.clearForError();
        return false;
      }
      ++text;
      c = *text;
      prevAtomIndex = temp.branchStack.back();
      temp.branchStack.pop_back();
    }

    isBondAllowed = true;
    isBranchAllowed = true;
  }
  if (temp.branchStack.size() != 0) {
    temp.clearForError();
    return false;
  }
  if (temp.numOpenBonds != 0) {
    temp.clearForError();
    return false;
  }

  // Post-processing

  // Prepare CSR graph
  toCSR(atoms.size(), bonds, mol);

  // Change any implicit bonds to single or aromatic.
  for (BondData& bond : bonds) {
    uint32_t atom0 = bond.getBeginAtomIdx();
    uint32_t atom1 = bond.getEndAtomIdx();
    if (bond.getBondType() == BondType::UNSPECIFIED) {
      const bool isAromatic0 = atoms[atom0].getIsAromatic();
      const bool isAromatic1 = atoms[atom1].getIsAromatic();
      if (isAromatic0 == isAromatic1) {
        // Both atoms aromatic or not
        bond.setBondType(isAromatic0 ? BondType::AROMATIC : BondType::SINGLE);
      } else {
#if 0
        const bool isUnknown0 = (atoms[atom0].getAtomicNum() == 0);
        const bool isUnknown1 = (atoms[atom1].getAtomicNum() == 0);
        // Exactly one of the two atoms is marked as aromatic.
        // Heuristic guess that if the bond is in a ring and the non-aromatic
        // atom is unknown, it should be aromatic, but just couldn't be marked
        // aromatic because it was unknown. Unknown atoms are rare, so this
        // almost always assigns a single bond type.
        bond.setBondType(((isUnknown0 || isUnknown1) && bond.isInRing)
                            ? BondType::AROMATIC
                            : BondType::SINGLE);
#endif
        // Because we're not determining bond rings here, just assign SINGLE.
        // FIXME: This should probably be handled properly.
        bond.setBondType(BondType::SINGLE);
      }
    }
  }
  return true;
}

}  // namespace SmilesParseInternal
