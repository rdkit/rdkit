//
//  Copyright (C) 2013-2014 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <cstring>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/FileParsers/FileParserUtils.h>
#include <RDGeneral/LocaleSwitcher.h>
#include "ProximityBonds.h"

#include <GraphMol/MonomerInfo.h>


// PDBWriter support multiple "flavors" of PDB output
// flavor & 1 : Ignore atoms in alternate conformations and dummy atoms
// flavor & 2 : Read each MODEL into a separate molecule.
namespace RDKit {

namespace {

constexpr int BCNAM(char A, char B, char C) { return (A << 16) | (B << 8) | C; }

Atom *PDBAtomFromSymbol(const char *symb) {
  PRECONDITION(symb, "bad char ptr");
  if (symb[0] == 'D' && !symb[1]) {
    auto *result = new Atom(1);
    result->setIsotope(2);
    return result;
  } else if (symb[0] == 'T' && !symb[1]) {
    auto *result = new Atom(1);
    result->setIsotope(3);
    return result;
  }
  int elemno = PeriodicTable::getTable()->getAtomicNumber(symb);
  return elemno > 0 ? new Atom(elemno) : (Atom *)nullptr;
}

void PDBAtomLine(RWMol *mol, const char *ptr, unsigned int len,
                 unsigned int flavor, std::map<int, Atom *> &amap) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");
  std::string tmp;

  if (len < 16) {
    return;
  }

  if ((flavor & 1) == 0) {
    // Ignore alternate locations of atoms.
    if (len >= 17 && ptr[16] != ' ' && ptr[16] != 'A' && ptr[16] != '1') {
      return;
    }
    // Ignore XPLOR pseudo atoms
    if (len >= 54 && !memcmp(ptr + 30, "9999.0009999.0009999.000", 24)) {
      return;
    }
    // Ignore NMR pseudo atoms
    if (ptr[12] == ' ' && ptr[13] == 'Q') {
      return;
    }
    // Ignore PDB dummy residues
    if (len >= 20 && !memcmp(ptr + 18, "DUM", 3)) {
      return;
    }
  }

  int serialno;
  tmp = std::string(ptr + 6, 5);
  try {
    serialno = FileParserUtils::toInt(tmp);
  } catch (boost::bad_lexical_cast &) {
    std::ostringstream errout;
    errout << "Non-integer PDB serial number " << tmp;
    throw FileParseException(errout.str());
  }

  Atom *atom = (Atom *)nullptr;
  char symb[3];

  // Attempt #1:  Atomic Symbol in columns 76 and 77
  if (len >= 78) {
    if (ptr[76] >= 'A' && ptr[76] <= 'Z') {
      symb[0] = ptr[76];
      if (ptr[77] >= 'A' && ptr[77] <= 'Z') {
        symb[1] = ptr[77] + 32;  // tolower
        symb[2] = '\0';
      } else if (ptr[77] >= 'a' && ptr[77] <= 'z') {
        symb[1] = ptr[77];
        symb[2] = '\0';
      } else {
        symb[1] = '\0';
      }
    } else if (ptr[76] == ' ' && ptr[77] >= 'A' && ptr[77] <= 'Z') {
      symb[0] = ptr[77];
      symb[1] = '\0';
    } else {
      symb[0] = '\0';
    }
  } else if (len == 77) {
    if (ptr[76] >= 'A' && ptr[76] <= 'Z') {
      symb[0] = ptr[76];
      symb[1] = '\0';
    } else {
      symb[0] = '\0';
    }
  } else {
    symb[0] = '\0';
  }

  if (symb[0]) {
    atom = PDBAtomFromSymbol(symb);
  }

  if (!atom) {
    // Attempt #2: Atomic Symbol from PDB atom name
    if (ptr[13] >= 'A' && ptr[13] <= 'Z') {
      if (ptr[12] == ' ') {
        symb[0] = ptr[13];
        if (ptr[14] >= 'a' && ptr[14] <= 'z') {
          symb[1] = ptr[14];
          symb[2] = '\0';
        } else {
          symb[1] = '\0';
        }
      } else if (ptr[12] >= 'A' && ptr[12] <= 'Z') {
        symb[0] = ptr[12];
        symb[1] = ptr[13] + 32;  // tolower
        symb[2] = '\0';
        if (ptr[12] == 'H' && ptr[0] == 'A') {
          // No He, Hf, Hg, Ho or Hs in ATOM records
          symb[0] = 'H';
          symb[1] = '\0';
        }
      } else if (ptr[12] >= '0' && ptr[12] <= '9') {
        symb[0] = ptr[13];
        symb[1] = '\0';
      } else {
        symb[0] = '\0';
      }
    } else {
      symb[0] = '\0';
    }

    if (symb[0]) {
      atom = PDBAtomFromSymbol(symb);
    }
  }

  if (!atom) {
    std::ostringstream errout;
    errout << "Cannot determine element for PDB atom #" << serialno;
    throw FileParseException(errout.str());
  }

  mol->addAtom(atom, true, true);
  amap[serialno] = atom;

  if (len >= 38) {
    RDGeom::Point3D pos;
    try {
      pos.x = FileParserUtils::toDouble(std::string(ptr + 30, 8));
      if (len >= 46) {
        pos.y = FileParserUtils::toDouble(std::string(ptr + 38, 8));
      }
      if (len >= 54) {
        pos.z = FileParserUtils::toDouble(std::string(ptr + 46, 8));
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with coordinates for PDB atom #" << serialno;
      throw FileParseException(errout.str());
    }

    Conformer *conf;
    if (!mol->getNumConformers()) {
      conf = new RDKit::Conformer(mol->getNumAtoms());
      conf->set3D(pos.z != 0.0);
      conf->setId(0);
      mol->addConformer(conf, false);
    } else {
      conf = &mol->getConformer();
      if (pos.z != 0.0) {
        conf->set3D(true);
      }
    }
    conf->setAtomPos(atom->getIdx(), pos);
  }

  if (len >= 79) {
    int charge = 0;
    if (ptr[78] >= '1' && ptr[78] <= '9') {
      if (ptr[79] == '-') {
        charge = -(ptr[78] - '0');
      } else if (ptr[79] == '+' || ptr[79] == ' ' || !ptr[79]) {
        charge = ptr[78] - '0';
      }
    } else if (ptr[78] == '+') {
      if (ptr[79] >= '1' && ptr[79] <= '9') {
        charge = ptr[79] - '0';
      } else if (ptr[79] == '+') {
        charge = 2;
      } else if (ptr[79] != '0') {
        charge = 1;
      }
    } else if (ptr[78] == '-') {
      if (ptr[79] >= '1' && ptr[79] <= '9') {
        charge = ptr[79] - '0';
      } else if (ptr[79] == '-') {
        charge = -2;
      } else if (ptr[79] != '0') {
        charge = -1;
      }
    } else if (ptr[78] == ' ') {
      if (ptr[79] >= '1' && ptr[79] <= '9') {
        charge = ptr[79] - '0';
      } else if (ptr[79] == '+') {
        charge = 1;
      } else if (ptr[79] == '-') {
        charge = -1;
      }
    }
    if (charge != 0) {
      atom->setFormalCharge(charge);
    }
  }

  tmp = std::string(ptr + 12, 4);
  AtomPDBResidueInfo *info = new AtomPDBResidueInfo(tmp, serialno);
  atom->setMonomerInfo(info);

  if (len >= 20) {
    tmp = std::string(ptr + 17, 3);
    // boost::trim(tmp);
  } else {
    tmp = "UNL";
  }
  info->setResidueName(tmp);
  if (ptr[0] == 'H') {
    info->setIsHeteroAtom(true);
  }
  if (len >= 17) {
    tmp = std::string(ptr + 16, 1);
  } else {
    tmp = " ";
  }
  info->setAltLoc(tmp);
  if (len >= 22) {
    tmp = std::string(ptr + 21, 1);
  } else {
    tmp = " ";
  }
  info->setChainId(tmp);
  if (len >= 27) {
    tmp = std::string(ptr + 26, 1);
  } else {
    tmp = " ";
  }
  info->setInsertionCode(tmp);

  int resno = 1;
  if (len >= 26) {
    try {
      resno = FileParserUtils::toInt(std::string(ptr + 22, 4));
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with residue number for PDB atom #" << serialno;
      throw FileParseException(errout.str());
    }
  }
  info->setResidueNumber(resno);

  double occup = 1.0;
  if (len >= 60) {
    try {
      occup = FileParserUtils::toDouble(std::string(ptr + 54, 6));
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with occupancy for PDB atom #" << serialno;
      throw FileParseException(errout.str());
    }
  }
  info->setOccupancy(occup);

  double bfactor = 0.0;
  if (len >= 66) {
    try {
      bfactor = FileParserUtils::toDouble(std::string(ptr + 60, 6));
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with temperature factor for PDB atom #" << serialno;
      throw FileParseException(errout.str());
    }
  }
  info->setTempFactor(bfactor);
}

void PDBBondLine(RWMol *mol, const char *ptr, unsigned int len,
                 std::map<int, Atom *> &amap, std::map<Bond *, int> &bmap) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");

  if (len < 16) {
    return;
  }

  std::string tmp(ptr + 6, 5);
  bool fail = false;
  int src, dst;

  try {
    src = FileParserUtils::toInt(tmp);
    if (amap.find(src) == amap.end()) {
      return;
    }
  } catch (boost::bad_lexical_cast &) {
    fail = true;
  }

  if (!fail) {
    if (len > 41) {
      len = 41;
    }
    for (unsigned int pos = 11; pos + 5 <= len; pos += 5) {
      if (!memcmp(ptr + pos, "     ", 5)) {
        break;
      }
      try {
        dst = FileParserUtils::toInt(std::string(ptr + pos, 5));
        if (dst == src || amap.find(dst) == amap.end()) {
          continue;
        }
      } catch (boost::bad_lexical_cast &) {
        fail = true;
      }

      if (!fail) {
        Bond *bond =
            mol->getBondBetweenAtoms(amap[src]->getIdx(), amap[dst]->getIdx());
        if (bond && bond->getBondType() != Bond::ZERO) {
          // Here we use a single byte bitmap to count duplicates
          // Low nibble counts src < dst, high nibble for src > dst
          int seen = bmap[bond];
          if (src < dst) {
            if ((seen & 0x0f) == 0x01) {
              bmap[bond] = seen | 0x02;
              if ((seen & 0x20) == 0) {
                bond->setBondType(Bond::DOUBLE);
              }
            } else if ((seen & 0x0f) == 0x03) {
              bmap[bond] = seen | 0x04;
              if ((seen & 0x40) == 0) {
                bond->setBondType(Bond::TRIPLE);
              }
            } else if ((seen & 0x0f) == 0x07) {
              bmap[bond] = seen | 0x08;
              if ((seen & 0x80) == 0) {
                bond->setBondType(Bond::QUADRUPLE);
              }
            }
          } else /* src < dst */ {
            if ((seen & 0xf0) == 0x10) {
              bmap[bond] = seen | 0x20;
              if ((seen & 0x02) == 0) {
                bond->setBondType(Bond::DOUBLE);
              }
            } else if ((seen & 0xf0) == 0x30) {
              bmap[bond] = seen | 0x40;
              if ((seen & 0x04) == 0) {
                bond->setBondType(Bond::TRIPLE);
              }
            } else if ((seen & 0xf0) == 0x70) {
              bmap[bond] = seen | 0x80;
              if ((seen & 0x08) == 0) {
                bond->setBondType(Bond::QUADRUPLE);
              }
            }
          }
        } else if (!bond) {
          // Bonds in PDB file are explicit
          // if they are not sanitize friendly, set their order to zero
          if (IsBlacklistedPair(amap[src], amap[dst])) {
            bond = new Bond(Bond::ZERO);
          } else {
            bond = new Bond(Bond::SINGLE);
          }
          bond->setOwningMol(mol);
          bond->setBeginAtom(amap[src]);
          bond->setEndAtom(amap[dst]);
          mol->addBond(bond, true);
          bmap[bond] = (src < dst) ? 0x01 : 0x10;
        } else {
          break;
        }
      } else {
        break;
      }
    }
  }

  if (fail) {
    std::ostringstream errout;
    errout << "Problem with CONECT record for PDB atom #" << tmp;
    throw FileParseException(errout.str());
  }
}

void PDBTitleLine(RWMol *mol, const char *ptr, unsigned int len) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");
  std::string title;
  while (ptr[len - 1] == ' ') {
    len--;
  }
  if (ptr[len - 1] == ';') {
    len--;
  }
  if (len > 21 && !strncmp(ptr + 10, " MOLECULE: ", 11)) {
    title = std::string(ptr + 21, len - 21);
  } else if (len > 10) {
    title = std::string(ptr + 10, len - 10);
  }
  if (!title.empty()) {
    mol->setProp(common_properties::_Name, title);
  }
}

void PDBConformerLine(RWMol *mol, const char *ptr, unsigned int len,
                      Conformer *&conf, int &conformer_atmidx) {
  PRECONDITION(mol, "bad mol");
  PRECONDITION(ptr, "bad char ptr");

  if (len >= 38) {
    RDGeom::Point3D pos;
    try {
      pos.x = FileParserUtils::toDouble(std::string(ptr + 30, 8));
      if (len >= 46) {
        pos.y = FileParserUtils::toDouble(std::string(ptr + 38, 8));
      }
      if (len >= 54) {
        pos.z = FileParserUtils::toDouble(std::string(ptr + 46, 8));
      }
    } catch (boost::bad_lexical_cast &) {
      std::ostringstream errout;
      errout << "Problem with multi-conformer coordinates";
      throw FileParseException(errout.str());
    }

    if (conformer_atmidx == 0) {
      conf = new RDKit::Conformer(mol->getNumAtoms());
      conf->setId(mol->getNumConformers());
      conf->set3D(pos.z != 0.0);
      mol->addConformer(conf, false);
    } else if (pos.z != 0.0) {
      conf->set3D(true);
    }

    if (conformer_atmidx < rdcast<int>(mol->getNumAtoms())) {
      conf->setAtomPos(conformer_atmidx, pos);
      conformer_atmidx++;
    }
  }
}

// This function determines whether a standard atom name in
// in a recognized PDB amino acid should be chiral or not.
// This is used to avoid chirality on VAL.CG and LEU.CG.
bool StandardPDBChiralAtom(const char *resnam, const char *atmnam) {
  switch (BCNAM(resnam[0], resnam[1], resnam[2])) {
    case BCNAM('G', 'L', 'Y'):
      return false;
    case BCNAM('I', 'L', 'E'):
    case BCNAM('T', 'H', 'R'):
      // Alpha and beta carbons (" CA " and " CB ").
      return atmnam[0] == ' ' && atmnam[1] == 'C' &&
             (atmnam[2] == 'A' || atmnam[2] == 'B') && atmnam[3] == ' ';
    case BCNAM('A', 'L', 'A'):
    case BCNAM('A', 'R', 'G'):
    case BCNAM('A', 'S', 'N'):
    case BCNAM('A', 'S', 'P'):
    case BCNAM('C', 'Y', 'S'):
    case BCNAM('G', 'L', 'N'):
    case BCNAM('G', 'L', 'U'):
    case BCNAM('H', 'I', 'S'):
    case BCNAM('L', 'E', 'U'):
    case BCNAM('L', 'Y', 'S'):
    case BCNAM('M', 'E', 'T'):
    case BCNAM('P', 'H', 'E'):
    case BCNAM('P', 'R', 'O'):
    case BCNAM('S', 'E', 'R'):
    case BCNAM('T', 'R', 'P'):
    case BCNAM('T', 'Y', 'R'):
    case BCNAM('V', 'A', 'L'):
      return atmnam[0] == ' ' && atmnam[1] == 'C' && atmnam[2] == 'A' &&
             atmnam[3] == ' ';
  }
  return false;
}

void StandardPDBResidueChirality(RWMol *mol) {
  for (ROMol::AtomIterator atomIt = mol->beginAtoms();
       atomIt != mol->endAtoms(); ++atomIt) {
    Atom *atom = *atomIt;
    if (atom->getChiralTag() != Atom::CHI_UNSPECIFIED) {
      auto *info = (AtomPDBResidueInfo *)atom->getMonomerInfo();
      if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE &&
          !info->getIsHeteroAtom() &&
          !StandardPDBChiralAtom(info->getResidueName().c_str(),
                                 info->getName().c_str())) {
        atom->setChiralTag(Atom::CHI_UNSPECIFIED);
        if (atom->hasProp(common_properties::_CIPCode)) {
          atom->clearProp(common_properties::_CIPCode);
        }
      }
    }
  }
}

void BasicPDBCleanup(RWMol &mol) {
  ROMol::VERTEX_ITER atBegin, atEnd;
  boost::tie(atBegin, atEnd) = mol.getVertices();
  while (atBegin != atEnd) {
    Atom *atom = mol[*atBegin];
    atom->calcExplicitValence(false);

    // correct four-valent neutral N -> N+
    // This was github #1029
    if (atom->getAtomicNum() == 7 && atom->getFormalCharge() == 0 &&
        atom->getExplicitValence() == 4) {
      atom->setFormalCharge(1);
    }
    ++atBegin;
  }
}

void parsePdbBlock(RWMol *&mol, const char *str, bool sanitize, bool removeHs,
                   unsigned int flavor, bool proximityBonding) {
  PRECONDITION(str, "bad char ptr");
  std::map<int, Atom *> amap;
  std::map<Bond *, int> bmap;
  Utils::LocaleSwitcher ls;
  bool multi_conformer = false;
  int conformer_atmidx = 0;
  Conformer *conf = nullptr;

  while (*str) {
    unsigned int len;
    const char *next = str + 1;
    for (;;) {
      if (*next == '\r') {
        len = (unsigned int)(next - str);
        if (next[1] == '\n') {
          next += 2;
        } else {
          next++;
        }
        break;
      } else if (*next == '\n') {
        len = (unsigned int)(next - str);
        next++;
        break;
      } else if (*next == '\0') {
        len = (unsigned int)(next - str);
        break;
      }
      next++;
    }

    // ATOM records
    if (str[0] == 'A' && str[1] == 'T' && str[2] == 'O' && str[3] == 'M' &&
        str[4] == ' ' && str[5] == ' ') {
      if (!multi_conformer) {
        if (!mol) {
          mol = new RWMol();
        }
        PDBAtomLine(mol, str, len, flavor, amap);
      } else {
        PDBConformerLine(mol, str, len, conf, conformer_atmidx);
      }
      // HETATM records
    } else if (str[0] == 'H' && str[1] == 'E' && str[2] == 'T' &&
               str[3] == 'A' && str[4] == 'T' && str[5] == 'M') {
      if (!multi_conformer) {
        if (!mol) {
          mol = new RWMol();
        }
        PDBAtomLine(mol, str, len, flavor, amap);
      } else {
        PDBConformerLine(mol, str, len, conf, conformer_atmidx);
      }
      // CONECT records
    } else if (str[0] == 'C' && str[1] == 'O' && str[2] == 'N' &&
               str[3] == 'E' && str[4] == 'C' && str[5] == 'T') {
      if (mol && !multi_conformer) {
        PDBBondLine(mol, str, len, amap, bmap);
      }
      // COMPND records
    } else if (str[0] == 'C' && str[1] == 'O' && str[2] == 'M' &&
               str[3] == 'P' && str[4] == 'N' && str[5] == 'D') {
      if (!mol) {
        mol = new RWMol();
      }
      if (len > 10 &&
          (str[9] == ' ' || !strncmp(str + 9, "2 MOLECULE: ", 12))) {
        PDBTitleLine(mol, str, len);
      }
      // HEADER records
    } else if (str[0] == 'H' && str[1] == 'E' && str[2] == 'A' &&
               str[3] == 'D' && str[4] == 'E' && str[5] == 'R') {
      if (!mol) {
        mol = new RWMol();
      }
      PDBTitleLine(mol, str, len < 50 ? len : 50);
      // ENDMDL records
    } else if (str[0] == 'E' && str[1] == 'N' && str[2] == 'D' &&
               str[3] == 'M' && str[4] == 'D' && str[5] == 'L') {
      if (!mol) {
        break;
      }
      multi_conformer = true;
      conformer_atmidx = 0;
      conf = nullptr;
    }
    str = next;
  }

  if (!mol) {
    return;
  }

  if (proximityBonding) {
    ConnectTheDots(mol, ctdIGNORE_H_H_CONTACTS);
  }
  // flavor & 8 doesn't encode double bonds
  if (proximityBonding || flavor & 8) {
    StandardPDBResidueBondOrders(mol);
  }

  BasicPDBCleanup(*mol);

  if (sanitize) {
    if (removeHs) {
      MolOps::removeHs(*mol, false, false);
    } else {
      MolOps::sanitizeMol(*mol);
    }
  } else {
    // we need some properties for the chiral setup
    mol->updatePropertyCache(false);
  }

  /* Set tetrahedral chirality from 3D co-ordinates */
  MolOps::assignChiralTypesFrom3D(*mol);
  StandardPDBResidueChirality(mol);
}
}  // namespace

RWMol *PDBBlockToMol(const char *str, bool sanitize, bool removeHs,
                     unsigned int flavor, bool proximityBonding) {
  RWMol *mol = nullptr;
  try {
    parsePdbBlock(mol, str, sanitize, removeHs, flavor, proximityBonding);
  } catch (...) {
    delete mol;
    throw;
  }

  return mol;
}

RWMol *PDBBlockToMol(const std::string &str, bool sanitize, bool removeHs,
                     unsigned int flavor, bool proximityBonding) {
  return PDBBlockToMol(str.c_str(), sanitize, removeHs, flavor,
                       proximityBonding);
}

RWMol *PDBDataStreamToMol(std::istream *inStream, bool sanitize, bool removeHs,
                          unsigned int flavor, bool proximityBonding) {
  PRECONDITION(inStream, "bad stream");
  std::string buffer;
  while (!inStream->eof() && !inStream->fail()) {
    std::string line;
    std::getline(*inStream, line);
    buffer += line;
    buffer += '\n';
    auto ptr = line.c_str();
    // Check for END
    if (ptr[0] == 'E' && ptr[1] == 'N' && ptr[2] == 'D' &&
        (ptr[3] == ' ' || ptr[3] == '\r' || ptr[3] == '\n' || !ptr[3])) {
      break;
    }
    // Check for ENDMDL
    if ((flavor & 2) != 0 && ptr[0] == 'E' && ptr[1] == 'N' && ptr[2] == 'D' &&
        ptr[3] == 'M' && ptr[4] == 'D' && ptr[5] == 'L') {
      break;
    }
  }
  return PDBBlockToMol(buffer.c_str(), sanitize, removeHs, flavor,
                       proximityBonding);
}
RWMol *PDBDataStreamToMol(std::istream &inStream, bool sanitize, bool removeHs,
                          unsigned int flavor, bool proximityBonding) {
  return PDBDataStreamToMol(&inStream, sanitize, removeHs, flavor,
                            proximityBonding);
}

RWMol *PDBFileToMol(const std::string &fileName, bool sanitize, bool removeHs,
                    unsigned int flavor, bool proximityBonding) {
  std::ifstream ifs(fileName.c_str(), std::ios_base::binary);
  if (!ifs || ifs.bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << fileName;
    throw BadFileException(errout.str());
  }
  return PDBDataStreamToMol(static_cast<std::istream *>(&ifs), sanitize,
                            removeHs, flavor, proximityBonding);
}

namespace {

constexpr size_t GROUP_PDB=0, ATOMID=1, TYPE_SYMBOL=2, ATOMNAME=3, ALT_ID=4, \
 RESNUM=5, RESNAME=6, CHAINID=7, INSCODE=8, OCCUPANCY=9, TEMPFACTOR=10, \
 CHARGE=11,CARTNX=12, CARTNY=13, CARTNZ=14, AUTH_RESNUM=15;
constexpr size_t N_ATTRS=16;
constexpr int NOT_PRESENT = -1;

bool is_quotation(char c) {
  return c == '"';
}

void parseMmcifAtoms(std::istream &ifs, const std::array<int, N_ATTRS> &name2column,
                     RWMol *&mol, unsigned int flavor) {
  std::vector<std::string> tokens;
  std::string line;

  mol = new RWMol();

  Conformer *conf = nullptr;
  if (name2column[CARTNX] != NOT_PRESENT ||
      name2column[CARTNY] != NOT_PRESENT ||
      name2column[CARTNZ] != NOT_PRESENT) {
    conf = new Conformer();
    conf->set3D(true);
    mol->addConformer(conf, false);
  }

  size_t atomnum = 0;
  while (ifs.peek() != '#') {
    std::getline(ifs, line);
    boost::split(tokens, line, boost::is_any_of(" "),
                 boost::token_compress_on);

    // skip altloc?
    if ((flavor & 1) == 0 && name2column[ALT_ID] != NOT_PRESENT) {
      if (tokens[name2column[ALT_ID]][0] != '.')
      continue;
    }

    std::string tmp = tokens[name2column[TYPE_SYMBOL]];
    // BR -> Br etc
    for (unsigned int i=1; i<tmp.length(); ++i) {
      tmp[i] = tolower(tmp[i]);
    }
    auto a = PDBAtomFromSymbol(tmp.c_str());

    if (!a) {
      std::ostringstream errout;
      errout << "Cannot determine element for PDB atom #" << atomnum;
      throw FileParseException(errout.str());
    }
    mol->addAtom(a, true, true);

    if (name2column[CHARGE] != NOT_PRESENT) {
      tmp = tokens[name2column[CHARGE]];
      int charge = 0;
      if (tmp[0] != '?') {
        try {
          charge = FileParserUtils::toInt(tokens[name2column[CHARGE]]);
        } catch (boost::bad_lexical_cast &) {
          std::ostringstream errout;
          errout << "Problem with formal charge for PDB atom #" << a->getIdx();
          throw FileParseException(errout.str());
        }
        if (charge != 0) {
          a->setFormalCharge(charge);
        }
      }
    }

    if (name2column[ATOMNAME != NOT_PRESENT]) {
      tmp = name2column[ATOMNAME];
      // "O5'" -> O5'
      boost::algorithm::trim_if(tmp, is_quotation);
    } else {
      tmp = "UNL";
    }

    tmp = name2column[ATOMNAME] != NOT_PRESENT ? tokens[name2column[ATOMNAME]] : "UNL";

    int serialno = 0;
    if (name2column[ATOMID] != NOT_PRESENT) {
      // these are usually but not always integers, use 0 if non integer
      try {
        serialno = FileParserUtils::toInt(tokens[name2column[ATOMID]]);
      } catch (boost::bad_lexical_cast &) {
        serialno = 0;
      }
    }

    auto info = new AtomPDBResidueInfo(tmp, serialno);
    a->setMonomerInfo(info);
    if (tokens[name2column[GROUP_PDB]][0] == 'H') {  // HETATM
      info->setIsHeteroAtom(true);
    }

    int resno = 1;
    if (name2column[AUTH_RESNUM] != NOT_PRESENT) {
      try {
        resno = FileParserUtils::toInt(tokens[name2column[AUTH_RESNUM]]);
      } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Problem with residue number for PDB atom #" << a->getIdx();
        throw FileParseException(errout.str());
      }
    } else if (name2column[RESNUM] != NOT_PRESENT) {
      try {
        resno = FileParserUtils::toInt(tokens[name2column[RESNUM]]);
      } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Problem with residue number for PDB atom #" << a->getIdx();
        throw FileParseException(errout.str());
      }
    }
    info->setResidueNumber(resno);

    tmp = name2column[RESNAME] != NOT_PRESENT ? tokens[name2column[RESNAME]] : "UNL";
    info->setResidueName(tmp);

    tmp = name2column[ALT_ID] != NOT_PRESENT ? tokens[name2column[ALT_ID]] : " ";
    if (tmp == ".") { // "." used as non value instead of " "
      tmp = " ";
    }
    info->setAltLoc(tmp);

    tmp = name2column[CHAINID] != NOT_PRESENT ? tokens[name2column[CHAINID]] : " ";
    if (tmp == ".") {
      tmp = " ";
    }
    info->setChainId(tmp);

    tmp = name2column[INSCODE] != NOT_PRESENT ? tokens[name2column[INSCODE]] : " ";
    if (tmp == ".") {
      tmp = " ";
    }
    info->setInsertionCode(tmp);

    double occup = 1.0;
    if (name2column[OCCUPANCY] != NOT_PRESENT) {
      try {
        occup = FileParserUtils::toDouble(tokens[name2column[OCCUPANCY]]);
      } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Problem with occupancy for PDB atom #" << a->getIdx();
        throw FileParseException(errout.str());
      }
    }
    info->setOccupancy(occup);

    double bfactor = 0.0;
    if (name2column[TEMPFACTOR] != NOT_PRESENT) {
      try {
        bfactor = FileParserUtils::toDouble(tokens[name2column[TEMPFACTOR]]);
      } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Problem with temperature factor for PDB atom #" << a->getIdx();
        throw FileParseException(errout.str());
      }
    }
    info->setTempFactor(bfactor);

    if (conf) {
      RDGeom::Point3D pos = {0.0, 0.0, 0.0};

      try {
        if (name2column[CARTNX] != NOT_PRESENT) {
          pos.x = FileParserUtils::toDouble(tokens[name2column[CARTNX]]);
        }
        if (name2column[CARTNY] != NOT_PRESENT) {
          pos.y = FileParserUtils::toDouble(tokens[name2column[CARTNY]]);
        }
        if (name2column[CARTNZ] != NOT_PRESENT) {
          pos.z = FileParserUtils::toDouble(tokens[name2column[CARTNZ]]);
        }
      } catch (boost::bad_lexical_cast &) {
        std::ostringstream errout;
        errout << "Problem with coordinates for PDB atom #" << a->getIdx();
        throw FileParseException(errout.str());
      }

      conf->setAtomPos(a->getIdx(), pos);
    }
    atomnum++;
  }
}


void parseMMCifFile(RWMol *&mol,
                    const std::string &fname, bool sanitize, bool removeHs,
                    unsigned int flavor, bool proximityBonding) {
  PRECONDITION(fname.length(), "no filename");

  // we need to read up to 14 columns, which can be in arbitrary order
  // this keeps track of which column is where, with -1 being column not present
  std::array<int, N_ATTRS> name2column;
  name2column.fill(NOT_PRESENT);

  std::ifstream ifs(fname.c_str(), std::ios_base::binary);
  if (!ifs || ifs.bad()) {
    std::ostringstream errout;
    errout << "Bad input file " << fname;
    throw BadFileException(errout.str());
  }

  // first find header of atom_site section
  std::string line;
  bool found_atom_site = false;
  while (!ifs.eof()) {
    std::getline(ifs, line);
    if (line != "loop_") {
      continue;
    }
    // check what sort of loop we've hit
    std::getline(ifs, line);
    if (line.compare(0, 11, "_atom_site.") == 0) {  // startswith
      found_atom_site = true;
      break;
    }
  }
  if (found_atom_site) {
    int i = 0;
    while (true) {  // sort out which columns we have and where they are
      line.erase(0, 11);  // remove "_atom_site." prefix
      boost::algorithm::trim_right(line);  // sometimes trailing whitespace

      if (line == "group_PDB") {
        name2column[GROUP_PDB] = i;
      } else if (line == "id") {
        name2column[ATOMID] = i;
      } else if (line == "type_symbol") {
        name2column[TYPE_SYMBOL] = i;
      } else if (line == "label_atom_id") {
        name2column[ATOMNAME] = i;
      } else if (line == "label_alt_id") {
        name2column[ALT_ID] = i;
      } else if (line == "label_seq_id") {
        name2column[RESNUM] = i;
      } else if (line == "auth_seq_id") {
        name2column[AUTH_RESNUM] = i;
      } else if (line == "label_comp_id") {
        name2column[RESNAME] = i;
      } else if (line == "label_asym_id") {
        name2column[CHAINID] = i;
      } else if (line == "pdbx_PDB_ins_code") {
        name2column[INSCODE] = i;
      } else if (line == "occupancy") {
        name2column[OCCUPANCY] = i;
      } else if (line == "B_iso_or_equiv") {
        name2column[TEMPFACTOR] = i;
      } else if (line == "pdbx_formal_charge") {
        name2column[CHARGE] = i;
      } else if (line == "Cartn_x") {
        name2column[CARTNX] = i;
      } else if (line == "Cartn_y") {
        name2column[CARTNY] = i;
      } else if (line == "Cartn_z") {
        name2column[CARTNZ] = i;
      }
      if (ifs.peek() != '_') {
        break;
      }
      std::getline(ifs, line);
      i++;
    }

    // slurp all rows from atom_site section
    {
      Utils::LocaleSwitcher ls;
      parseMmcifAtoms(ifs, name2column, mol, flavor);
    }
  }
  if (!mol) {
    return;
  }

  if (proximityBonding) {
    ConnectTheDots(mol, ctdIGNORE_H_H_CONTACTS);
  }
  // flavor & 8 doesn't encode double bonds
  if (proximityBonding || flavor & 8) {
    StandardPDBResidueBondOrders(mol);
  }

  BasicPDBCleanup(*mol);

  if (sanitize) {
    if (removeHs) {
      MolOps::removeHs(*mol, false, false);
    } else {
      MolOps::sanitizeMol(*mol);
    }
  } else {
    // we need some properties for the chiral setup
    mol->updatePropertyCache(false);
  }

  /* Set tetrahedral chirality from 3D co-ordinates */
  MolOps::assignChiralTypesFrom3D(*mol);
  StandardPDBResidueChirality(mol);
}
}  // namespace

RWMol *mmcifFileToMol(const std::string &fname, bool sanitize, bool removeHs,
                      unsigned int flavor, bool proximityBonding) {
  RWMol *mol = nullptr;
  try {
    parseMMCifFile(mol, fname, sanitize, removeHs, flavor, proximityBonding);
  } catch (...) {
    delete mol;
    throw;
  }

  return mol;
}


}  // namespace RDKit
