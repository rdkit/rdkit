//
//  Copyright (C) 2016 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <iostream>
#include "SmilesParse.h"
#include "SmilesParseOps.h"

namespace SmilesParseOps {
using namespace RDKit;

namespace parser {
template <typename Iterator>
bool read_int(Iterator &first, Iterator last, unsigned int &res) {
  std::string num = "";
  while (first != last && *first >= '0' && *first <= '9') {
    num += *first;
    ++first;
  }
  if (num == "") {
    return false;
  }
  res = boost::lexical_cast<unsigned int>(num);
  return true;
}
template <typename Iterator>
bool read_int_pair(Iterator &first, Iterator last, unsigned int &n1,
                   unsigned int &n2, char sep = '.') {
  if (!read_int(first, last, n1)) return false;
  if (first >= last || *first != sep) return false;
  ++first;
  return read_int(first, last, n2);
}

template <typename Iterator>
std::string read_text_to(Iterator &first, Iterator last, std::string delims) {
  std::string res = "";
  Iterator start = first;
  // EFF: there are certainly faster ways to do this
  while (first != last && delims.find_first_of(*first) == std::string::npos) {
    if (*first == '&' && std::distance(first, last) > 2 &&
        *(first + 1) == '#') {
      // escaped char
      if (start != first) {
        res += std::string(start, first);
      }
      Iterator next = first + 2;
      while (next != last && *next >= '0' && *next <= '9') {
        ++next;
      }
      if (next == last || *next != ';')
        throw RDKit::SmilesParseException(
            "failure parsing CXSMILES extensions: quoted block not terminated "
            "with ';'");
      if (next > first + 2) {
        std::string blk = std::string(first + 2, next);
        res += (char)(boost::lexical_cast<int>(blk));
      }
      first = next + 1;
      start = first;
    } else {
      ++first;
    }
  }
  if (start != first) res += std::string(start, first);
  return res;
}

template <typename Iterator>
bool parse_atom_values(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last || *first != ':') return false;
  ++first;
  unsigned int atIdx = 0;
  while (first != last && *first != '$') {
    std::string tkn = read_text_to(first, last, ";$");
    if (tkn != "") {
      mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
                                         tkn);
    }
    ++atIdx;
    if (first != last && *first != '$') ++first;
  }
  if (first == last || *first != '$') return false;
  ++first;
  return true;
}

template <typename Iterator>
bool parse_atom_props(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last) return false;
  while (first != last && *first != '|' && *first != ',') {
    unsigned int atIdx;
    if (read_int(first, last, atIdx)) {
      if (first == last || *first != '.') return false;
      ++first;
      std::string pname = read_text_to(first, last, ".");
      if (pname != "") {
        if (first == last || *first != '.') return false;
        ++first;
        std::string pval = read_text_to(first, last, ":|,");
        if (pval != "") {
          mol.getAtomWithIdx(atIdx)->setProp(pname, pval);
        }
      }
    }
    if (first != last && *first != '|' && *first != ',') ++first;
  }
  if (first != last && *first != '|' && *first != ',') return false;
  if (*first != '|') ++first;
  return true;
}

template <typename Iterator>
bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last || *first != '$') return false;
  ++first;
  unsigned int atIdx = 0;
  while (first != last && *first != '$') {
    std::string tkn = read_text_to(first, last, ";$");
    if (tkn != "") {
      mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::atomLabel,
                                         tkn);
    }
    ++atIdx;
    if (first != last && *first != '$') ++first;
  }
  if (first == last || *first != '$') return false;
  ++first;
  return true;
}

template <typename Iterator>
bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last || *first != '(') return false;

  auto *conf = new Conformer(mol.getNumAtoms());
  mol.addConformer(conf);
  ++first;
  unsigned int atIdx = 0;
  while (first != last && *first != ')') {
    RDGeom::Point3D pt;
    std::string tkn = read_text_to(first, last, ";)");
    if (tkn != "") {
      std::vector<std::string> tokens;
      boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
      if (tokens.size() >= 1 && tokens[0].size())
        pt.x = boost::lexical_cast<double>(tokens[0]);
      if (tokens.size() >= 2 && tokens[1].size())
        pt.y = boost::lexical_cast<double>(tokens[1]);
      if (tokens.size() >= 3 && tokens[2].size())
        pt.z = boost::lexical_cast<double>(tokens[2]);
    }

    conf->setAtomPos(atIdx, pt);
    ++atIdx;
    if (first != last && *first != ')') ++first;
  }
  if (first == last || *first != ')') return false;
  ++first;
  return true;
}

template <typename Iterator>
bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last || *first != 'C') return false;
  ++first;
  if (first >= last || *first != ':') return false;
  ++first;
  while (first != last && *first >= '0' && *first <= '9') {
    unsigned int aidx;
    unsigned int bidx;
    if (read_int_pair(first, last, aidx, bidx)) {
      Bond *bnd = mol.getBondWithIdx(bidx);
      if (bnd->getBeginAtomIdx() != aidx && bnd->getEndAtomIdx() != aidx) {
        std::cerr << "BOND NOT FOUND! " << bidx << " involving atom " << aidx
                  << std::endl;
        return false;
      }
      bnd->setBondType(Bond::DATIVE);
      if (bnd->getBeginAtomIdx() != aidx) {
        unsigned int tmp = bnd->getBeginAtomIdx();
        bnd->setBeginAtomIdx(aidx);
        bnd->setEndAtomIdx(tmp);
      }
    } else {
      return false;
    }
    if (first < last && *first == ',') ++first;
  }
  return true;
}

template <typename Iterator>
bool processRadicalSection(Iterator &first, Iterator last, RDKit::RWMol &mol,
                           unsigned int numRadicalElectrons) {
  if (first >= last) return false;
  ++first;
  if (first >= last || *first != ':') return false;
  ++first;
  unsigned int atIdx;
  if (!read_int(first, last, atIdx)) return false;
  mol.getAtomWithIdx(atIdx)->setNumRadicalElectrons(numRadicalElectrons);
  while (first < last && *first == ',') {
    ++first;
    if (first < last && (*first < '0' || *first > '9')) return true;
    if (!read_int(first, last, atIdx)) return false;
    mol.getAtomWithIdx(atIdx)->setNumRadicalElectrons(numRadicalElectrons);
  }
  if (first >= last) return false;
  return true;
}

template <typename Iterator>
bool parse_radicals(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last || *first != '^') return false;
  while (*first == '^') {
    ++first;
    if (first >= last) return false;
    if (*first < '1' || *first > '7')
      return false;  // these are the values that are allowed to be there
    switch (*first) {
      case '1':
        if (!processRadicalSection(first, last, mol, 1)) return false;
        break;
      case '2':
      case '3':
      case '4':
        if (!processRadicalSection(first, last, mol, 2)) return false;
        break;
      case '5':
      case '6':
      case '7':
        if (!processRadicalSection(first, last, mol, 3)) return false;
        break;
      default:
        BOOST_LOG(rdWarningLog) << "Radical specification " << *first
                                << " ignored.";
    }
  }
  return true;
}

template <typename Iterator>
bool parse_it(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  if (first >= last || *first != '|') return false;
  ++first;
  while (first < last && *first != '|') {
    typename Iterator::difference_type length = std::distance(first, last);
    if (*first == '(') {
      if (!parse_coords(first, last, mol)) return false;
    } else if (*first == '$') {
      if (length > 4 && *(first + 1) == '_' && *(first + 2) == 'A' &&
          *(first + 3) == 'V' && *(first + 4) == ':') {
        first += 4;
        if (!parse_atom_values(first, last, mol)) return false;
      } else {
        if (!parse_atom_labels(first, last, mol)) return false;
      }
    } else if (length > 9 && std::string(first, first + 9) == "atomProp:") {
      first += 9;
      if (!parse_atom_props(first, last, mol)) return false;
    } else if (*first == 'C') {
      if (!parse_coordinate_bonds(first, last, mol)) return false;
    } else if (*first == '^') {
      if (!parse_radicals(first, last, mol)) return false;
    } else {
      ++first;
    }
    // if(first < last && *first != '|') ++first;
  }
  if (first >= last || *first != '|') return false;
  ++first;  // step past the last '|'
  return true;
}
}  // end of namespace parser

namespace {
template <typename Q>
void addquery(Q *qry, std::string symbol, RDKit::RWMol &mol, unsigned int idx) {
  PRECONDITION(qry, "bad query");
  auto *qa = new QueryAtom(0);
  qa->setQuery(qry);
  qa->setNoImplicit(true);
  mol.replaceAtom(idx, qa);
  if (symbol != "")
    mol.getAtomWithIdx(idx)->setProp(RDKit::common_properties::atomLabel,
                                     symbol);
  delete qa;
}
void processCXSmilesLabels(RDKit::RWMol &mol) {
  for (RDKit::ROMol::AtomIterator atIt = mol.beginAtoms();
       atIt != mol.endAtoms(); ++atIt) {
    std::string symb = "";
    if ((*atIt)->getPropIfPresent(RDKit::common_properties::atomLabel, symb)) {
      if (symb.size() > 3 && symb[0] == '_' && symb[1] == 'A' &&
          symb[2] == 'P') {
        unsigned int mapNum =
            boost::lexical_cast<unsigned int>(symb.substr(3, symb.size() - 3));
        (*atIt)->setAtomMapNum(mapNum);
      } else if (symb == "star_e") {
        /* according to the MDL spec, these match anything, but in MARVIN they
        are "unspecified end groups" for polymers */
        addquery(makeAtomNullQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "Q_e") {
        addquery(makeQAtomQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "QH_p") {
        addquery(makeQHAtomQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "AH_p") {  // this seems wrong...
        /* According to the MARVIN Sketch, AH is "any atom, including H" -
        this would be "*" in SMILES - and "A" is "any atom except H".
        The CXSMILES docs say that "A" can be represented normally in SMILES
        and that "AH" needs to be written out as AH_p. I'm going to assume that
        this is a Marvin internal thing and just parse it as they describe it.
        This means that "*" in the SMILES itself needs to be treated
        differently, which we do below. */
        addquery(makeAHAtomQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "X_p") {
        addquery(makeXAtomQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "XH_p") {
        addquery(makeXHAtomQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "M_p") {
        addquery(makeMAtomQuery(), symb, mol, (*atIt)->getIdx());
      } else if (symb == "MH_p") {
        addquery(makeMHAtomQuery(), symb, mol, (*atIt)->getIdx());
      }
    } else if ((*atIt)->getAtomicNum() == 0 && (*atIt)->getSymbol() == "*") {
      addquery(makeAAtomQuery(), "", mol, (*atIt)->getIdx());
    }
  }
}

}  // end of anonymous namespace

void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                       std::string::const_iterator &first) {
  // std::cerr << "parseCXNExtensions: " << extText << std::endl;
  if (!extText.size() || extText[0] != '|') return;
  first = extText.begin();
  bool ok = parser::parse_it(first, extText.end(), mol);
  if (!ok)
    throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
  if (ok) {
    processCXSmilesLabels(mol);
  }
}
}  // end of namespace SmilesParseOps
