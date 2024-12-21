//
//  Copyright (C) 2015,2016 Greg Landrum and NextMove Software
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

#include "SequenceWriters.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/MonomerInfo.h>

namespace RDKit {

static char getOneLetterAACode(const AtomPDBResidueInfo *info) {
  const char *ptr = info->getResidueName().c_str();
  switch (ptr[0]) {
    case 'A':
      if (!strcmp(ptr, "ALA")) {
        return 'A';
      }
      if (!strcmp(ptr, "ARG")) {
        return 'R';
      }
      if (!strcmp(ptr, "ASN")) {
        return 'N';
      }
      if (!strcmp(ptr, "ASP")) {
        return 'D';
      }
      break;
    case 'C':
      if (!strcmp(ptr, "CYS")) {
        return 'C';
      }
      break;
    case 'D':
      if (!strcmp(ptr, "DAL")) {
        return 'a';
      }
      if (!strcmp(ptr, "DAR")) {
        return 'r';
      }
      if (!strcmp(ptr, "DAS")) {
        return 'd';
      }
      if (!strcmp(ptr, "DCY")) {
        return 'c';
      }
      if (!strcmp(ptr, "DGL")) {
        return 'e';
      }
      if (!strcmp(ptr, "DGN")) {
        return 'q';
      }
      if (!strcmp(ptr, "DHI")) {
        return 'h';
      }
      if (!strcmp(ptr, "DIL")) {
        return 'i';
      }
      if (!strcmp(ptr, "DLE")) {
        return 'l';
      }
      if (!strcmp(ptr, "DLY")) {
        return 'k';
      }
      if (!strcmp(ptr, "DPN")) {
        return 'f';
      }
      if (!strcmp(ptr, "DPR")) {
        return 'p';
      }
      if (!strcmp(ptr, "DSG")) {
        return 'n';
      }
      if (!strcmp(ptr, "DSN")) {
        return 's';
      }
      if (!strcmp(ptr, "DTH")) {
        return 't';
      }
      if (!strcmp(ptr, "DTR")) {
        return 'w';
      }
      if (!strcmp(ptr, "DTY")) {
        return 'y';
      }
      if (!strcmp(ptr, "DVA")) {
        return 'v';
      }
      break;
    case 'G':
      if (!strcmp(ptr, "GLU")) {
        return 'E';
      }
      if (!strcmp(ptr, "GLN")) {
        return 'Q';
      }
      if (!strcmp(ptr, "GLY")) {
        return 'G';
      }
      break;
    case 'H':
      if (!strcmp(ptr, "HIS")) {
        return 'H';
      }
      break;
    case 'I':
      if (!strcmp(ptr, "ILE")) {
        return 'I';
      }
      break;
    case 'L':
      if (!strcmp(ptr, "LEU")) {
        return 'L';
      }
      if (!strcmp(ptr, "LYS")) {
        return 'K';
      }
      break;
    case 'M':
      if (!strcmp(ptr, "MET")) {
        return 'M';
      }
      if (!strcmp(ptr, "MED")) {
        return 'm';
      }
      break;
    case 'P':
      if (!strcmp(ptr, "PHE")) {
        return 'F';
      }
      if (!strcmp(ptr, "PRO")) {
        return 'P';
      }
      break;
    case 'S':
      if (!strcmp(ptr, "SER")) {
        return 'S';
      }
      break;
    case 'T':
      if (!strcmp(ptr, "THR")) {
        return 'T';
      }
      if (!strcmp(ptr, "TRP")) {
        return 'W';
      }
      if (!strcmp(ptr, "TYR")) {
        return 'Y';
      }
      break;
    case 'V':
      if (!strcmp(ptr, "VAL")) {
        return 'V';
      }
      break;
  }
  return 'X';
}

static char getOneLetterNACode(const AtomPDBResidueInfo *info) {
  const char *ptr = info->getResidueName().c_str();
  if (ptr[0] == ' ' && (ptr[1] == ' ' || ptr[1] == 'D')) {
    switch (ptr[2]) {
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'U':
        return ptr[2];
    }
  }
  return 'X';
}

std::string MolToSequence(const ROMol &mol) {
  std::string result;
  std::string chain;

  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    const Atom *atom = *atomIt;
    auto *info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
    if (info && info->getMonomerType() == AtomMonomerInfo::PDBRESIDUE) {
      if (info->getName() == " CA ") {
        if (!chain.empty() && chain != info->getChainId()) {
          chain = info->getChainId();
          result += '.';
        }
        result += getOneLetterAACode(info);
      } else if (info->getName() == " C1'") {
        if (!chain.empty() && chain != info->getChainId()) {
          chain = info->getChainId();
          result += '.';
        }
        result += getOneLetterNACode(info);
      }
    }
  }
  return result;
}

std::string MolToFASTA(const ROMol &mol) {
  std::string seq = MolToSequence(mol);
  if (seq.empty()) {
    return "";
  }

  std::string result = ">";
  std::string name;
  if (mol.getPropIfPresent(common_properties::_Name, name)) {
    result += name;
  }
  result += '\n';
  result += seq;
  result += '\n';
  return result;
}

static const char *getHELMAAMonomer(const AtomPDBResidueInfo *info) {
  const char *ptr = info->getResidueName().c_str();
  switch (ptr[0]) {
    case 'A':
      if (!strcmp(ptr, "ALA")) {
        return "A";
      }
      if (!strcmp(ptr, "ARG")) {
        return "R";
      }
      if (!strcmp(ptr, "ASN")) {
        return "N";
      }
      if (!strcmp(ptr, "ASP")) {
        return "D";
      }
      break;
    case 'C':
      if (!strcmp(ptr, "CYS")) {
        return "C";
      }
      break;
    case 'D':
      if (!strcmp(ptr, "DAL")) {
        return "[dA]";
      }
      if (!strcmp(ptr, "DAR")) {
        return "[dR]";
      }
      if (!strcmp(ptr, "DAS")) {
        return "[dD]";
      }
      if (!strcmp(ptr, "DCY")) {
        return "[dC]";
      }
      if (!strcmp(ptr, "DGL")) {
        return "[dE]";
      }
      if (!strcmp(ptr, "DGN")) {
        return "[dQ]";
      }
      if (!strcmp(ptr, "DHI")) {
        return "[dH]";
      }
      if (!strcmp(ptr, "DLE")) {
        return "[dL]";
      }
      if (!strcmp(ptr, "DLY")) {
        return "[dK]";
      }
      if (!strcmp(ptr, "DPN")) {
        return "[dF]";
      }
      if (!strcmp(ptr, "DPR")) {
        return "[dP]";
      }
      if (!strcmp(ptr, "DSG")) {
        return "[dN]";
      }
      if (!strcmp(ptr, "DSN")) {
        return "[dS]";
      }
      if (!strcmp(ptr, "DTR")) {
        return "[dW]";
      }
      if (!strcmp(ptr, "DTY")) {
        return "[dY]";
      }
      if (!strcmp(ptr, "DVA")) {
        return "[dV]";
      }
      break;
    case 'G':
      if (!strcmp(ptr, "GLU")) {
        return "E";
      }
      if (!strcmp(ptr, "GLN")) {
        return "Q";
      }
      if (!strcmp(ptr, "GLY")) {
        return "G";
      }
      break;
    case 'H':
      if (!strcmp(ptr, "HIS")) {
        return "H";
      }
      break;
    case 'I':
      if (!strcmp(ptr, "ILE")) {
        return "I";
      }
      break;
    case 'L':
      if (!strcmp(ptr, "LEU")) {
        return "L";
      }
      if (!strcmp(ptr, "LYS")) {
        return "K";
      }
      break;
    case 'M':
      if (!strcmp(ptr, "MET")) {
        return "M";
      }
      if (!strcmp(ptr, "MED")) {
        return "[dM]";
      }
      if (!strcmp(ptr, "MSE")) {
        return "[seC]";
      }
      break;
    case 'N':
      if (!strcmp(ptr, "NLE")) {
        return "[Nle]";
      }
      if (!strcmp(ptr, "NVA")) {
        return "[Nva]";
      }
      break;
    case 'O':
      if (!strcmp(ptr, "ORN")) {
        return "[Orn]";
      }
      break;
    case 'P':
      if (!strcmp(ptr, "PHE")) {
        return "F";
      }
      if (!strcmp(ptr, "PRO")) {
        return "P";
      }
      break;
    case 'S':
      if (!strcmp(ptr, "SER")) {
        return "S";
      }
      break;
    case 'T':
      if (!strcmp(ptr, "THR")) {
        return "T";
      }
      if (!strcmp(ptr, "TRP")) {
        return "W";
      }
      if (!strcmp(ptr, "TYR")) {
        return "Y";
      }
      break;
    case 'V':
      if (!strcmp(ptr, "VAL")) {
        return "V";
      }
      break;
  }
  return (const char *)nullptr;
}

static bool IsEupeptideBond(AtomPDBResidueInfo *src, AtomPDBResidueInfo *dst) {
  if (src->getChainId() != dst->getChainId()) {
    return false;
  }

  int sresno = src->getResidueNumber();
  int dresno = dst->getResidueNumber();

  // Peptides
  if (src->getName() == " C  " && dst->getName() == " N  ") {
    if (sresno != dresno) {
      return dresno > sresno;
    }
    return dst->getInsertionCode() > src->getInsertionCode();
  }
  if (src->getName() == " N  " && dst->getName() == " C  ") {
    if (sresno != dresno) {
      return dresno < sresno;
    }
    return dst->getInsertionCode() < src->getInsertionCode();
  }

  // Nucleic acids
  if (src->getName() == " O3'" && dst->getName() == " P  ") {
    if (sresno != dresno) {
      return dresno > sresno;
    }
    return dst->getInsertionCode() > src->getInsertionCode();
  }
  if (src->getName() == " P  " && dst->getName() == " O3'") {
    if (sresno != dresno) {
      return dresno > sresno;
    }
    return dst->getInsertionCode() > src->getInsertionCode();
  }
  return false;
}

static bool IsSupportedHELMBond(AtomPDBResidueInfo *src,
                                AtomPDBResidueInfo *dst) {
  if (src->getName() == " SG " && dst->getName() == " SG ") {
    if ((src->getResidueName() == "CYS" || src->getResidueName() == "DCY") &&
        (dst->getResidueName() == "CYS" || dst->getResidueName() == "DCY")) {
      return true;
    }
  }

  if (src->getName() == " N  " && dst->getName() == " C  ") {
    return true;
  }
  if (src->getName() == " C  " && dst->getName() == " N  ") {
    return true;
  }
  return false;
}

static bool FindHELMAtom(std::vector<AtomPDBResidueInfo *> *seq,
                         AtomPDBResidueInfo *info, std::string &id,
                         std::string &pos) {
  char buffer[32];
  char ch;

  const char *ptr = info->getName().c_str();
  if (ptr[0] == ' ' && ptr[1] == 'S' && ptr[2] == 'G' && ptr[3] == ' ') {
    ch = '3';
  } else if (ptr[0] == ' ' && ptr[1] == 'N' && ptr[2] == ' ' && ptr[3] == ' ') {
    ch = '1';
  } else if (ptr[0] == ' ' && ptr[1] == 'C' && ptr[2] == ' ' && ptr[3] == ' ') {
    ch = '2';
  } else {
    return false;
  }

  int resno = info->getResidueNumber();
  for (unsigned int i = 1; i < 10; i++) {
    auto len = (unsigned int)seq[i].size();
    for (unsigned int j = 0; j < len; j++) {
      AtomPDBResidueInfo *targ = seq[i][j];
      if (targ->getResidueNumber() == resno &&
          targ->getResidueName() == info->getResidueName() &&
          targ->getChainId() == info->getChainId() &&
          targ->getInsertionCode() == info->getInsertionCode()) {
        id = "PEPTIDE";
        id += (char)(i + '0');
        sprintf(buffer, "%u:R%c", j + 1, ch);
        pos = buffer;
        return true;
      }
    }
  }
  return false;
}

static std::string NameHELMBond(std::vector<AtomPDBResidueInfo *> *seq,
                                AtomPDBResidueInfo *src,
                                AtomPDBResidueInfo *dst) {
  std::string id1, pos1;
  if (!FindHELMAtom(seq, src, id1, pos1)) {
    return "";
  }
  std::string id2, pos2;
  if (!FindHELMAtom(seq, dst, id2, pos2)) {
    return "";
  }

  std::string result = id1;
  result += ',';
  result += id2;
  result += ',';
  result += pos1;
  result += '-';
  result += pos2;
  return result;
}

static const char *getHELMNAMonomer(const AtomPDBResidueInfo *info) {
  const char *ptr = info->getResidueName().c_str();
  if (ptr[0] == ' ') {
    if (ptr[1] == ' ') {
      switch (ptr[2]) {
        case 'A':
          return "R(A)";
        case 'C':
          return "R(C)";
        case 'G':
          return "R(G)";
        case 'T':
          return "R(T)";
        case 'U':
          return "R(U)";
      }
    } else if (ptr[1] == 'D') {
      switch (ptr[2]) {
        case 'A':
          return "[dR](A)";
        case 'C':
          return "[dR](C)";
        case 'G':
          return "[dR](G)";
        case 'T':
          return "[dR](T)";
        case 'U':
          return "[dR](U)";
      }
    }
  }
  return (const char *)nullptr;
}

// Pstate finite state machine
// 0 start of string	P -> 1     B -> 2	T -> fail
// 1 after 5'-P		P -> fail  B -> 3	T -> fail
// 2 after base		P -> 4     B -> fail	T -> 5
// 3 after 5'-base	P -> 1	   B -> fail	T -> 5
// 4 after base-3'	P -> fail  B -> 2	T -> fail
// 5 after terminal-3'	P -> fail  B -> fail	T -> fail

std::string MolToHELM(const ROMol &mol) {
  std::vector<AtomPDBResidueInfo *> seq[10];
  std::string result;
  unsigned int Pstate = 0;
  bool peptide = true;
  bool first = true;
  std::string chain;
  int id = 1;

  /* First pass: Monomers */
  for (ROMol::ConstAtomIterator atomIt = mol.beginAtoms();
       atomIt != mol.endAtoms(); ++atomIt) {
    const Atom *atom = *atomIt;
    auto *info = (AtomPDBResidueInfo *)(atom->getMonomerInfo());
    // We can only write HELM if all atoms have PDB residue information
    if (!info || info->getMonomerType() != AtomMonomerInfo::PDBRESIDUE) {
      return "";
    }

    if (info->getName() == " CA ") {
      const char *mono = getHELMAAMonomer(info);
      if (!mono) {
        return "";
      }
      if (first) {
        chain = info->getChainId();
        result = "PEPTIDE1{";
        peptide = true;
        first = false;
      } else if (!peptide || info->getChainId() != chain) {
        // Nine chains should be enough?
        if (id == 9) {
          return "";
        }
        id++;
        chain = info->getChainId();
        result += "}|PEPTIDE";
        result += (char)(id + '0');
        result += "{";
        peptide = true;
      } else {
        result += ".";
      }
      result += mono;
      seq[id].push_back(info);
    } else if (info->getResidueName() == "NH2" && info->getName() == " N  ") {
      if (first || !peptide) {
        return "";
      }
      result += ".[am]";
    } else if (info->getResidueName() == "ACE" && info->getName() == " C  ") {
      if (first) {
        chain = info->getChainId();
        result = "PEPTIDE1{[ac]";
        peptide = true;
        first = false;
      } else if (info->getChainId() != chain) {
        // Nine chains should be enough?
        if (id == 9) {
          return "";
        }
        id++;
        chain = info->getChainId();
        result += "}|PEPTIDE";
        result += (char)(id + '0');
        result += "{[ac]";
        peptide = true;
      } else {
        return "";
      }
      seq[id].push_back(info);
    } else if (info->getName() == " C1'") {
      const char *mono = getHELMNAMonomer(info);
      if (!mono) {
        return "";
      }
      if (first) {
        chain = info->getChainId();
        result = "RNA1{";
        peptide = false;
        first = false;
        Pstate = 2;
      } else if (peptide || info->getChainId() != chain) {
        // Nine chains should be enough?
        if (id == 9) {
          return "";
        }
        id++;
        chain = info->getChainId();
        result += "}|RNA";
        result += (char)(id + '0');
        result += "{";
        peptide = false;
        Pstate = 2;
      } else {
        switch (Pstate) {
          case 1:
          case 4:
            result += '.';
            Pstate = 2;
            break;
          default:
            return "";
        }
      }
      result += mono;
    } else if (info->getName() == " P  ") {
      if (getHELMNAMonomer(info)) {
        if (first) {
          chain = info->getChainId();
          result = "RNA1{P";
          peptide = false;
          first = false;
          Pstate = 1;
        } else if (peptide || info->getChainId() != chain) {
          // Nine chains should be enough?
          if (id == 9) {
            return "";
          }
          id++;
          chain = info->getChainId();
          result += "}|RNA";
          result += (char)(id + '0');
          result += "{P";
          peptide = false;
          Pstate = 1;
        } else {
          switch (Pstate) {
            case 2:
              result += 'P';
              Pstate = 4;
              break;
            case 3:
              result += ".P";
              Pstate = 1;
              break;
            default:
              return "";
          }
        }
      } else if (info->getResidueName() == "2PO") {
        if (first || peptide || info->getChainId() != chain) {
          return "";
        }
        switch (Pstate) {
          case 2:
            result += 'P';
            break;
          case 3:
            result += ".P";
            break;
          default:
            return "";
        }
        Pstate = 5;
      }
    }
  }

  if (first) {
    return "";
  }

  result += "}$";

  first = true;
  for (ROMol::ConstBondIterator bondIt = mol.beginBonds();
       bondIt != mol.endBonds(); ++bondIt) {
    const Bond *bond = *bondIt;
    Atom *beg = bond->getBeginAtom();
    Atom *end = bond->getEndAtom();
    if (!beg || !end) {
      continue;
    }
    auto *binfo = (AtomPDBResidueInfo *)(beg->getMonomerInfo());
    auto *einfo = (AtomPDBResidueInfo *)(end->getMonomerInfo());
    if (!binfo || !einfo) {
      continue;
    }
    // Test if this is an uninteresting intra-residue bond
    if (binfo->getResidueNumber() == einfo->getResidueNumber() &&
        binfo->getResidueName() == einfo->getResidueName() &&
        binfo->getChainId() == einfo->getChainId() &&
        binfo->getInsertionCode() == einfo->getInsertionCode()) {
      continue;
    }
    if (bond->getBondType() != Bond::SINGLE) {
      return "";
    }
    if (IsEupeptideBond(binfo, einfo)) {
      continue;
    }
    if (!IsSupportedHELMBond(binfo, einfo)) {
      return "";
    }
    std::string tmp = NameHELMBond(seq, binfo, einfo);
    if (tmp.empty()) {
      return "";
    }
    if (!first) {
      result += "|";
    } else {
      first = false;
    }
    result += tmp;
  }

  result += "$$$";
  return result;
}

}  // namespace RDKit
