//
//  Copyright (C) 2013 Greg Landrum and NextMove Software
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>

#ifdef RDK_BUILD_MAEPARSER_SUPPORT
#undef RDK_BUILD_MAEPARSER_SUPPORT
#endif
#include <GraphMol/FileParsers/MolWriters.h>
#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/MonomerInfo.h>
#include <GraphMol/RDKitBase.h>
#include <RDGeneral/BadFileException.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/LocaleSwitcher.h>

// PDBWriter support multiple "flavors" of PDB output
// flavor & 1 : Write MODEL/ENDMDL lines around each record
// flavor & 2 : Don't write any CONECT records
// flavor & 4 : Write CONECT records in both directions
// flavor & 8 : Don't use multiple CONECTs to encode bond order
// flavor & 16 : Write MASTER record
// flavor & 32 : Write TER record

namespace RDKit {

// Get the next atom number of an element, formatted as a 2-letter string.
std::string GetDefaultAtomNumber(const Atom* atom,
        std::map<unsigned int, unsigned int> &elem);

std::string GetPDBAtomLine(const Atom *atom, const Conformer *conf,
                           std::map<unsigned int, unsigned int> &elem) {
  PRECONDITION(atom, "bad atom");
  std::stringstream ss;

  std::string symb = atom->getSymbol();
  char at1, at2;
  switch (symb.length()) {
    case 0:
      at1 = ' ';
      at2 = 'X';
      break;
    case 1:
      at1 = ' ';
      at2 = symb[0];
      break;
    default:
      at1 = symb[0];
      at2 = symb[1];
      if (at2 >= 'a' && at2 <= 'z') {
        at2 -= 32;  // toupper
      }
      break;
  }

  auto *info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
  if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
    ss << (info->getIsHeteroAtom() ? "HETATM" : "ATOM  ");
    ss << std::setw(5) << atom->getIdx() + 1;
    ss << ' ';
    const std::string& name = info->getName();
    if (name.empty()) {
      std::string atnum = GetDefaultAtomNumber(atom, elem);
      ss << at1 << at2 << atnum;
    } else {
      ss << std::setw(4) << name.substr(0, 4);
    }
    const char *ptr = info->getAltLoc().c_str();
    if (*ptr == '\0') {
      ptr = " ";
    }
    ss << *ptr;
    ss << std::setw(3) << info->getResidueName().substr(0, 3);
    ss << ' ';
    ptr = info->getChainId().c_str();
    if (*ptr == '\0') {
      ptr = " ";
    }
    ss << *ptr;
    ss << std::setw(4) << info->getResidueNumber();
    ptr = info->getInsertionCode().c_str();
    if (*ptr == '\0') {
      ptr = " ";
    }
    ss << *ptr;
    ss << "   ";
  } else {
    info = (AtomPDBResidueInfo *)nullptr;
    std::string atnum = GetDefaultAtomNumber(atom, elem);
    ss << "HETATM";
    ss << std::setw(5) << atom->getIdx() + 1;
    ss << ' ';
    ss << at1 << at2 << atnum;
    ss << " UNL     1    ";
  }

  if (conf) {
    const RDGeom::Point3D pos = conf->getAtomPos(atom->getIdx());
    ss << boost::format("%8.3f%8.3f%8.3f") % pos.x % pos.y % pos.z;
  } else {
    ss << "   0.000   0.000   0.000";
  }

  if (info) {
    ss << boost::format("%6.2f%6.2f") % info->getOccupancy() %
              info->getTempFactor();
    ss << "          ";
  } else {
    ss << "  1.00  0.00          ";
  }

  ss << at1;
  ss << at2;
  int charge = atom->getFormalCharge();
  if (charge > 0 && charge < 10) {
    ss << (char)('0' + charge);
    ss << '+';
  } else if (charge < 0 && charge > -10) {
    ss << (char)('0' - charge);
    ss << '-';
  } else {
    ss << "  ";
  }
  return ss.str();
}

std::string GetDefaultAtomNumber(const Atom* atom, std::map<unsigned int, unsigned int> &elem) {
  std::string ret = "  ";
  unsigned int atno = atom->getAtomicNum();
  if (elem.find(atno) == elem.end()) {
    elem[atno] = 1;
    ret[0] = '1';
  } else {
    unsigned int tmp = elem[atno] + 1;
    elem[atno] = tmp;
    if (tmp < 10) {
      ret[0] = tmp + '0';
    } else if (tmp < 100) {
      ret[0] = (tmp / 10) + '0';
      ret[1] = (tmp % 10) + '0';
    } else if (tmp < 360) {
      ret[0] = ((tmp - 100) / 10) + 'A';
      ret[1] = ((tmp - 100) % 10) + '0';
    } else if (tmp < 1036) {
      ret[0] = ((tmp - 360) / 26) + 'A';
      ret[1] = ((tmp - 360) % 26) + 'A';
    }
  }
  return ret;
}

std::string GetPDBBondLines(const Atom *atom, bool all, bool both, bool mult,
                            unsigned int &conect_count) {
  PRECONDITION(atom, "bad atom");
  unsigned int src = atom->getIdx() + 1;
  std::vector<unsigned int> v;

  ROMol *mol = &atom->getOwningMol();
  for (ROMol::OBOND_ITER_PAIR bondIt = mol->getAtomBonds(atom);
       bondIt.first != bondIt.second; ++bondIt.first) {
    Bond *bptr = (*mol)[*bondIt.first];
    Atom *nptr = bptr->getOtherAtom(atom);
    unsigned int dst = nptr->getIdx() + 1;
    if (dst < src && !both) {
      continue;
    }
    Bond::BondType btype = Bond::SINGLE;
    if (mult) {
      btype = bptr->getBondType();
    }
    switch (btype) {
      default:
      case Bond::SINGLE:
      case Bond::AROMATIC:
        if (all) {
          v.push_back(dst);
        }
        break;
      case Bond::QUADRUPLE:
        v.push_back(dst);
        /* FALLTHRU */
      case Bond::TRIPLE:
        v.push_back(dst);
        /* FALLTHRU */
      case Bond::DOUBLE:
        v.push_back(dst);
        v.push_back(dst);
    }
  }

  auto count = rdcast<unsigned int>(v.size());
  if (count == 0) {
    return "";
  }

  std::sort(v.begin(), v.end());
  std::stringstream ss;
  for (unsigned int i = 0; i < count; i++) {
    if ((i & 3) == 0) {
      if (i != 0) {
        ss << '\n';
      }
      ss << "CONECT";
      ss << std::setw(5) << src;
      conect_count++;
    }
    ss << std::setw(5) << v[i];
  }
  ss << '\n';
  return ss.str();
}

std::string MolToPDBBody(const ROMol &mol, const Conformer *conf,
                         unsigned int flavor, unsigned int &atm_count,
                         unsigned int &ter_count, unsigned int &conect_count) {
  std::string res;
  std::string last;
  std::map<unsigned int, unsigned int> elem;
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    last = GetPDBAtomLine(*atomIt, conf, elem);
    res += last;
    res += '\n';
    atm_count++;
  }

  if (ter_count == 0 && atm_count && (flavor & 32)) {
    std::stringstream ss;
    ss << "TER   ";
    ss << std::setw(5) << atm_count + 1;
    if (last.length() >= 27) {
      ss << "      ";
      ss << last.substr(17, 10);
    }
    ss << '\n';
    res += ss.str();
    ter_count = 1;
  }

  bool all = (flavor & 2) == 0;
  bool both = (flavor & 4) != 0;
  bool mult = (flavor & 8) == 0;
  if (all || mult) {
    for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
         atomIt != mol.endAtoms(); ++atomIt) {
      res += GetPDBBondLines(*atomIt, all, both, mult, conect_count);
    }
  }
  return res;
}

std::string MolToPDBBlock(const ROMol &imol, int confId, unsigned int flavor) {
  ROMol mol(imol);
  auto &trwmol = static_cast<RWMol &>(mol);
  MolOps::Kekulize(trwmol);
  Utils::LocaleSwitcher ls;

  std::string res;
  std::string name;
  if (mol.getPropIfPresent(common_properties::_Name, name)) {
    if (!name.empty()) {
      res += "COMPND    ";
      res += name;
      res += '\n';
    }
  }

  unsigned int atm_count = 0;
  unsigned int ter_count = 0;
  unsigned int conect_count = 0;

  const Conformer *conf;
  if (confId < 0 && mol.getNumConformers() > 1) {
    int count = mol.getNumConformers();
    for (confId = 0; confId < count; confId++) {
      conf = &(mol.getConformer(confId));
      std::stringstream ss;
      ss << "MODEL     ";
      ss << std::setw(4) << (confId + 1);
      ss << "\n";
      res += ss.str();
      res +=
          MolToPDBBody(mol, conf, flavor, atm_count, ter_count, conect_count);
      res += "ENDMDL\n";
    }
  } else {
    if (confId < 0 && mol.getNumConformers() == 0) {
      conf = nullptr;
    } else {
      conf = &(mol.getConformer(confId));
    }
    res += MolToPDBBody(mol, conf, flavor, atm_count, ter_count, conect_count);
  }

  if (flavor & 16) {
    std::stringstream ss;
    ss << "MASTER        0    0    0    0    0    0    0    0";
    ss << std::setw(5) << atm_count;
    ss << std::setw(5) << ter_count;
    ss << std::setw(5) << conect_count;
    ss << "    0\n";
    res += ss.str();
  }

  res += "END\n";
  return res;
}

PDBWriter::PDBWriter(const std::string &fileName, unsigned int flavor) {
  if (fileName != "-") {
    auto *tmpStream = new std::ofstream(fileName.c_str());
    df_owner = true;
    if (!(*tmpStream) || (tmpStream->bad())) {
      delete tmpStream;
      std::ostringstream errout;
      errout << "Bad output file " << fileName;
      throw BadFileException(errout.str());
    }
    dp_ostream = static_cast<std::ostream *>(tmpStream);
  } else {
    dp_ostream = static_cast<std::ostream *>(&std::cout);
    df_owner = false;
  }
  d_flavor = flavor;
  d_count = 0;
}

PDBWriter::PDBWriter(std::ostream *outStream, bool takeOwnership,
                     unsigned int flavor) {
  PRECONDITION(outStream, "null stream");
  if (outStream->bad()) {
    throw FileParseException("Bad output stream");
  }
  dp_ostream = outStream;
  df_owner = takeOwnership;
  d_flavor = flavor;
  d_count = 0;
}

PDBWriter::~PDBWriter() {
  // close the writer if it's still open:
  if (dp_ostream != nullptr) {
    close();
  }
}

void PDBWriter::write(const ROMol &mol, int confId) {
  PRECONDITION(dp_ostream, "no output stream");

  d_count++;
  if (d_flavor & 1) {
    std::stringstream ss;
    ss << "MODEL     ";
    ss << std::setw(4) << d_count;
    ss << "\n";
    (*dp_ostream) << ss.str();
  }

  // write the molecule
  (*dp_ostream) << MolToPDBBlock(mol, confId, d_flavor);

  if (d_flavor & 1) {
    (*dp_ostream) << "ENDMDL\n";
  }
}

void MolToPDBFile(const ROMol &mol, const std::string &fname, int confId,
                  unsigned int flavor) {
  PDBWriter w(fname, flavor);
  w.write(mol, confId);
}

}  // namespace RDKit
