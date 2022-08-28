//
//  Copyright (C) 2016-2021 Greg Landrum and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RDGeneral/BoostStartInclude.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <iostream>
#include <algorithm>
#include "SmilesWrite.h"
#include "SmilesParse.h"
#include "SmilesParseOps.h"
#include <GraphMol/MolEnumerator/LinkNode.h>

namespace SmilesParseOps {
using namespace RDKit;

const std::string cxsmilesindex = "_cxsmilesindex";

// FIX: once this can be automated using constexpr, do so
const std::vector<std::string_view> pseudoatoms{"Pol", "Mod"};
const std::vector<std::string_view> pseudoatoms_p{"Pol_p", "Mod_p"};

std::map<std::string, std::string> sgroupTypemap = {
    {"n", "SRU"},   {"mon", "MON"}, {"mer", "MER"}, {"co", "COP"},
    {"xl", "CRO"},  {"mod", "MOD"}, {"mix", "MIX"}, {"f", "FOR"},
    {"any", "ANY"}, {"gen", "GEN"}, {"c", "COM"},   {"grf", "GRA"},
    {"alt", "COP"}, {"ran", "COP"}, {"blk", "COP"}};

template <typename Q>
void addquery(Q *qry, std::string symbol, RDKit::RWMol &mol, unsigned int idx) {
  PRECONDITION(qry, "bad query");
  auto *qa = new QueryAtom(0);
  qa->setQuery(qry);
  qa->setNoImplicit(true);
  mol.replaceAtom(idx, qa);
  if (symbol != "") {
    mol.getAtomWithIdx(idx)->setProp(RDKit::common_properties::atomLabel,
                                     symbol);
  }
  delete qa;
}

void processCXSmilesLabels(RWMol &mol) {
  if (mol.hasProp("_cxsmilesLabelsProcessed")) {
    return;
  }
  for (auto atom : mol.atoms()) {
    std::string symb = "";
    if (atom->getPropIfPresent(common_properties::atomLabel, symb)) {
      if (symb == "star_e") {
        /* according to the MDL spec, these match anything, but in MARVIN they
        are "unspecified end groups" for polymers */
        addquery(makeAtomNullQuery(), symb, mol, atom->getIdx());
      } else if (symb == "Q_e") {
        addquery(makeQAtomQuery(), symb, mol, atom->getIdx());
      } else if (symb == "QH_p") {
        addquery(makeQHAtomQuery(), symb, mol, atom->getIdx());
      } else if (symb == "AH_p") {  // this seems wrong...
        /* According to the MARVIN Sketch, AH is "any atom, including H" -
        this would be "*" in SMILES - and "A" is "any atom except H".
        The CXSMILES docs say that "A" can be represented normally in SMILES
        and that "AH" needs to be written out as AH_p. I'm going to assume that
        this is a Marvin internal thing and just parse it as they describe it.
        This means that "*" in the SMILES itself needs to be treated
        differently, which we do below. */
        addquery(makeAHAtomQuery(), symb, mol, atom->getIdx());
      } else if (symb == "X_p") {
        addquery(makeXAtomQuery(), symb, mol, atom->getIdx());
      } else if (symb == "XH_p") {
        addquery(makeXHAtomQuery(), symb, mol, atom->getIdx());
      } else if (symb == "M_p") {
        addquery(makeMAtomQuery(), symb, mol, atom->getIdx());
      } else if (symb == "MH_p") {
        addquery(makeMHAtomQuery(), symb, mol, atom->getIdx());
      } else if (std::find(pseudoatoms_p.begin(), pseudoatoms_p.end(), symb) !=
                 pseudoatoms_p.end()) {
        // strip off the "_p":
        atom->setProp(common_properties::dummyLabel,
                      symb.substr(0, symb.size() - 2));
        atom->clearProp(common_properties::atomLabel);
      }
    } else if (atom->getAtomicNum() == 0 && !atom->hasQuery() &&
               atom->getSymbol() == "*") {
      addquery(makeAAtomQuery(), "", mol, atom->getIdx());
    }
  }
  mol.setProp("_cxsmilesLabelsProcessed", 1, true);
}

namespace parser {

const std::string _headCrossings = "_headCrossings";
const std::string _tailCrossings = "_tailCrossings";

template <typename Iterator>
bool read_int(Iterator &first, Iterator last, unsigned int &res) {
  std::string num = "";
  while (first <= last && *first >= '0' && *first <= '9') {
    num += *first;
    ++first;
  }
  if (num.empty()) {
    return false;
  }
  res = boost::lexical_cast<unsigned int>(num);
  return true;
}
template <typename Iterator>
bool read_int_list(Iterator &first, Iterator last,
                   std::vector<unsigned int> &res, char sep = ',') {
  while (1) {
    std::string num = "";
    while (first <= last && *first >= '0' && *first <= '9') {
      num += *first;
      ++first;
    }
    if (!num.empty()) {
      res.push_back(boost::lexical_cast<unsigned int>(num));
    }
    if (first >= last || *first != sep) {
      break;
    }
    ++first;
  }
  return true;
}
template <typename Iterator>
bool read_int_pair(Iterator &first, Iterator last, unsigned int &n1,
                   unsigned int &n2, char sep = '.') {
  if (!read_int(first, last, n1)) {
    return false;
  }
  if (first >= last || *first != sep) {
    return false;
  }
  ++first;
  return read_int(first, last, n2);
}

template <typename Iterator>
std::string read_text_to(Iterator &first, Iterator last, std::string delims) {
  std::string res = "";
  Iterator start = first;
  // EFF: there are certainly faster ways to do this
  while (first <= last && delims.find_first_of(*first) == std::string::npos) {
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
      if (next == last || *next != ';') {
        throw RDKit::SmilesParseException(
            "failure parsing CXSMILES extensions: quoted block not terminated "
            "with ';'");
      }
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
  if (start != first) {
    res += std::string(start, first);
  }
  return res;
}
namespace {

// this is the super fun case where no information about bonds in/out of the
// sgroup is present.
void setupUnmarkedPolymerSGroup(RWMol &mol, SubstanceGroup &sgroup,
                                std::vector<unsigned int> &headCrossings,
                                std::vector<unsigned int> &tailCrossings) {
  const auto &atoms = sgroup.getAtoms();
  if (atoms.empty()) {
    throw SmilesParseException("no atoms in polymer sgroup");
  }
  const auto firstAtom = mol.getAtomWithIdx(atoms.front());
  for (auto nbr : boost::make_iterator_range(mol.getAtomNeighbors(firstAtom))) {
    const auto nbrAtom = mol[nbr];
    if (std::find(atoms.begin(), atoms.end(), nbrAtom->getIdx()) ==
        atoms.end()) {
      // in most cases we just add this to the set of headCrossings.
      // The exception occurs when there's only one atom in the SGroup and
      //  we already have a headCrossing, in which case we may put this one
      //  as a tailCrossing
      if (atoms.size() > 1 || headCrossings.empty()) {
        headCrossings.push_back(
            mol.getBondBetweenAtoms(firstAtom->getIdx(), nbrAtom->getIdx())
                ->getIdx());
      } else if (atoms.size() == 1) {
        if (tailCrossings.empty()) {
          tailCrossings.push_back(
              mol.getBondBetweenAtoms(firstAtom->getIdx(), nbrAtom->getIdx())
                  ->getIdx());
        } else {
          BOOST_LOG(rdWarningLog)
              << " single atom polymer Sgroup has more than two bonds to "
                 "external atoms. Ignoring all bonds after the first two."
              << std::endl;
        }
      }
    }
  }
  if (atoms.size() > 1) {
    const auto lastAtom = mol.getAtomWithIdx(atoms.back());
    for (auto nbr :
         boost::make_iterator_range(mol.getAtomNeighbors(lastAtom))) {
      const auto nbrAtom = mol[nbr];
      if (std::find(atoms.begin(), atoms.end(), nbrAtom->getIdx()) ==
          atoms.end()) {
        tailCrossings.push_back(
            mol.getBondBetweenAtoms(lastAtom->getIdx(), nbrAtom->getIdx())
                ->getIdx());
      }
    }
  }
}

// deal with setting up the crossing bonds, etc.
void finalizePolymerSGroup(RWMol &mol, SubstanceGroup &sgroup) {
  bool isFlipped = false;
  std::string connect = "EU";
  if (sgroup.getPropIfPresent("CONNECT", connect)) {
    if (connect.find(",f") != std::string::npos) {
      isFlipped = true;
      boost::replace_all(connect, ",f", "");
    }
  }
  if (connect == "hh") {
    connect = "HH";
  } else if (connect == "ht") {
    connect = "HT";
  } else if (connect == "eu") {
    connect = "EU";
  } else {
    BOOST_LOG(rdWarningLog) << "unrecognized CXSMILES CONNECT value: '"
                            << connect << "'. Assuming 'eu'" << std::endl;
    connect = "EU";
  }
  sgroup.setProp("CONNECT", connect);

  std::vector<unsigned int> headCrossings;
  std::vector<unsigned int> tailCrossings;
  sgroup.getPropIfPresent(_headCrossings, headCrossings);
  sgroup.clearProp(_headCrossings);
  sgroup.getPropIfPresent(_tailCrossings, tailCrossings);
  sgroup.clearProp(_tailCrossings);
  if (headCrossings.empty() && tailCrossings.empty()) {
    setupUnmarkedPolymerSGroup(mol, sgroup, headCrossings, tailCrossings);
  }
  if (headCrossings.empty() && tailCrossings.empty()) {
    // we tried... nothing more we can do
    return;
  }

  for (auto &bondIdx : headCrossings) {
    sgroup.addBondWithIdx(bondIdx);
  }
  sgroup.setProp("XBHEAD", headCrossings);

  for (auto &bondIdx : tailCrossings) {
    sgroup.addBondWithIdx(bondIdx);
  }

  // now we can setup XBCORR
  std::vector<unsigned int> xbcorr;
  for (unsigned int i = 0;
       i < std::min(headCrossings.size(), tailCrossings.size()); ++i) {
    unsigned headIdx = headCrossings[i];
    unsigned tailIdx = tailCrossings[i];
    if (isFlipped) {
      tailIdx = tailCrossings[tailCrossings.size() - i - 1];
    }
    xbcorr.push_back(headIdx);
    xbcorr.push_back(tailIdx);
  }
  sgroup.setProp("XBCORR", xbcorr);
}

}  // end of anonymous namespace

// we use this pattern a lot and it's a long function call, but a very short
// #define
#define VALID_ATIDX(_atidx_) \
  ((_atidx_) >= startAtomIdx && (_atidx_) < startAtomIdx + mol.getNumAtoms())

template <typename Iterator>
bool parse_atom_values(Iterator &first, Iterator last, RDKit::RWMol &mol,
                       unsigned int startAtomIdx) {
  if (first >= last || *first != ':') {
    return false;
  }
  ++first;
  unsigned int atIdx = 0;
  while (first <= last && *first != '$') {
    std::string tkn = read_text_to(first, last, ";$");
    if (tkn != "" && VALID_ATIDX(atIdx)) {
      mol.getAtomWithIdx(atIdx)->setProp(RDKit::common_properties::molFileValue,
                                         tkn);
    }
    ++atIdx;
    if (first <= last && *first != '$') {
      ++first;
    }
  }
  if (first >= last || *first != '$') {
    return false;
  }
  ++first;
  return true;
}

template <typename Iterator>
bool parse_atom_props(Iterator &first, Iterator last, RDKit::RWMol &mol,
                      unsigned int startAtomIdx) {
  if (first >= last) {
    return false;
  }
  while (first <= last && *first != '|' && *first != ',') {
    unsigned int atIdx;
    if (read_int(first, last, atIdx)) {
      if (first >= last || *first != '.') {
        return false;
      }
      ++first;
      std::string pname = read_text_to(first, last, ".");
      if (!pname.empty()) {
        if (first >= last || *first != '.') {
          return false;
        }
        ++first;
        std::string pval = read_text_to(first, last, ":|,");
        if (VALID_ATIDX(atIdx) && !pval.empty()) {
          mol.getAtomWithIdx(atIdx - startAtomIdx)->setProp(pname, pval);
        }
      }
    }
    if (first <= last && *first != '|' && *first != ',') {
      ++first;
    }
  }
  if (first <= last && *first != '|' && *first != ',') {
    return false;
  }
  if (*first != '|') {
    ++first;
  }
  return true;
}

template <typename Iterator>
bool parse_atom_labels(Iterator &first, Iterator last, RDKit::RWMol &mol,
                       unsigned int startAtomIdx) {
  if (first >= last || *first != '$') {
    return false;
  }
  ++first;
  unsigned int atIdx = 0;
  while (first <= last && *first != '$') {
    std::string tkn = read_text_to(first, last, ";$");
    if (!tkn.empty() && VALID_ATIDX(atIdx)) {
      mol.getAtomWithIdx(atIdx - startAtomIdx)
          ->setProp(RDKit::common_properties::atomLabel, tkn);
    }
    ++atIdx;
    if (first <= last && *first != '$') {
      ++first;
    }
  }
  if (first >= last || *first != '$') {
    return false;
  }
  ++first;
  return true;
}

template <typename Iterator>
bool parse_coords(Iterator &first, Iterator last, RDKit::RWMol &mol,
                  unsigned int startAtomIdx) {
  if (first >= last || *first != '(') {
    return false;
  }

  auto *conf = new Conformer(mol.getNumAtoms());
  mol.addConformer(conf);
  ++first;
  unsigned int atIdx = 0;
  while (first <= last && *first != ')') {
    RDGeom::Point3D pt;
    std::string tkn = read_text_to(first, last, ";)");
    if (VALID_ATIDX(atIdx)) {
      if (!tkn.empty()) {
        std::vector<std::string> tokens;
        boost::split(tokens, tkn, boost::is_any_of(std::string(",")));
        if (tokens.size() >= 1 && tokens[0].size()) {
          pt.x = boost::lexical_cast<double>(tokens[0]);
        }
        if (tokens.size() >= 2 && tokens[1].size()) {
          pt.y = boost::lexical_cast<double>(tokens[1]);
        }
        if (tokens.size() >= 3 && tokens[2].size()) {
          pt.z = boost::lexical_cast<double>(tokens[2]);
        }
      }

      conf->setAtomPos(atIdx - startAtomIdx, pt);
    }
    ++atIdx;
    if (first <= last && *first != ')') {
      ++first;
    }
  }
  if (first >= last || *first != ')') {
    return false;
  }
  ++first;
  return true;
}

template <typename Iterator>
bool parse_coordinate_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                            Bond::BondType typ, unsigned int startAtomIdx,
                            unsigned int startBondIdx) {
  if (first >= last || (*first != 'C' && *first != 'H')) {
    return false;
  }
  ++first;
  if (first >= last || *first != ':') {
    return false;
  }
  ++first;
  while (first <= last && *first >= '0' && *first <= '9') {
    unsigned int aidx;
    unsigned int bidx;
    if (read_int_pair(first, last, aidx, bidx)) {
      if (VALID_ATIDX(aidx) && bidx >= startBondIdx &&
          bidx < startBondIdx + mol.getNumBonds()) {
        Bond *bnd = nullptr;
        for (auto bond : mol.bonds()) {
          unsigned int smilesIdx;
          if (bond->getPropIfPresent("_cxsmilesBondIdx", smilesIdx) &&
              smilesIdx + startBondIdx == bidx) {
            bnd = bond;
            break;
          }
        }
        if (!bnd || (bnd->getBeginAtomIdx() != aidx - startAtomIdx &&
                     bnd->getEndAtomIdx() != aidx - startAtomIdx)) {
          BOOST_LOG(rdWarningLog) << "BOND NOT FOUND! " << bidx
                                  << " involving atom " << aidx << std::endl;
          return false;
        }
        bnd->setBondType(typ);
        if (bnd->getBeginAtomIdx() != aidx - startAtomIdx) {
          unsigned int tmp = bnd->getBeginAtomIdx();
          bnd->setBeginAtomIdx(aidx - startAtomIdx);
          bnd->setEndAtomIdx(tmp);
        }
      }
    } else {
      return false;
    }
    if (first < last && *first == ',') {
      ++first;
    }
  }
  return true;
}

template <typename Iterator>
bool parse_unsaturation(Iterator &first, Iterator last, RDKit::RWMol &mol,
                        unsigned int startAtomIdx) {
  if (first + 1 >= last || *first != 'u') {
    return false;
  }
  ++first;
  if (first >= last || *first != ':') {
    return false;
  }
  ++first;
  while (first < last && *first >= '0' && *first <= '9') {
    unsigned int idx;
    if (!read_int(first, last, idx)) {
      return false;
    }
    if (VALID_ATIDX(idx)) {
      auto atom = mol.getAtomWithIdx(idx - startAtomIdx);
      if (!atom->hasQuery()) {
        atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
      }
      atom->expandQuery(makeAtomUnsaturatedQuery(), Queries::COMPOSITE_AND);
    }
    if (first < last && *first == ',') {
      ++first;
    }
  }
  return true;
}

template <typename Iterator>
bool parse_ring_bonds(Iterator &first, Iterator last, RDKit::RWMol &mol,
                      unsigned int startAtomIdx) {
  if (first >= last || *first != 'r' || first + 1 >= last ||
      *(first + 1) != 'b' || first + 2 >= last || *(first + 2) != ':') {
    return false;
  }
  first += 3;
  while (first < last && *first >= '0' && *first <= '9') {
    unsigned int n1;
    if (!read_int(first, last, n1)) {
      return false;
    }
    // check that we can read at least two more characters:
    if (first + 1 >= last || *first != ':') {
      return false;
    }
    ++first;
    unsigned int n2;
    bool gt = false;
    if (*first == '*') {
      ++first;
      n2 = 0xDEADBEEF;
      if (VALID_ATIDX(n1)) {
        mol.setProp(common_properties::_NeedsQueryScan, 1);
      }
    } else {
      if (!read_int(first, last, n2)) {
        return false;
      }
      switch (n2) {
        case 0:
        case 2:
        case 3:
          break;
        case 4:
          gt = true;
          break;
        default:
          BOOST_LOG(rdWarningLog)
              << "unrecognized rb value: " << n2 << std::endl;
          return false;
      }
    }
    if (VALID_ATIDX(n1)) {
      auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
      if (!atom->hasQuery()) {
        atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
      }
      if (!gt) {
        atom->expandQuery(makeAtomRingBondCountQuery(n2),
                          Queries::COMPOSITE_AND);
      } else {
        auto q = static_cast<ATOM_EQUALS_QUERY *>(new ATOM_LESSEQUAL_QUERY);
        q->setVal(n2);
        q->setDescription("AtomRingBondCount");
        q->setDataFunc(queryAtomRingBondCount);
        atom->expandQuery(q, Queries::COMPOSITE_AND);
      }
    }
    if (first < last && *first == ',') {
      ++first;
    }
  }
  return true;
}

template <typename Iterator>
bool parse_linknodes(Iterator &first, Iterator last, RDKit::RWMol &mol,
                     unsigned int startAtomIdx) {
  // these look like: |LN:1:1.3.2.6,4:1.4.3.6|
  // that's two records:
  //   1:1.3.2.6: 1-3 repeats, atom 1-2, 1-6
  //   4:1.4.3.6: 1-4 repeats, atom 4-3, 4-6
  // which maps to the property value "1 3 2 2 3 2 7|1 4 2 5 4 5 7"
  // If the linking atom only has two neighbors then the outer atom
  // specification (the last two digits) can be left out. So for a molecule
  // where atom 1 has bonds only to atoms 2 and 6 we could have
  // |LN:1:1.3|
  // instead of
  // |LN:1:1.3.2.6|
  if (first >= last || *first != 'L' || first + 1 >= last ||
      *(first + 1) != 'N' || first + 2 >= last || *(first + 2) != ':') {
    return false;
  }
  first += 3;
  std::string accum = "";
  while (first < last && *first >= '0' && *first <= '9') {
    unsigned int atidx;
    if (!read_int(first, last, atidx)) {
      return false;
    }
    // check that we can read at least two more characters:
    if (first + 1 >= last || *first != ':') {
      return false;
    }
    ++first;
    unsigned int startReps;
    if (!read_int(first, last, startReps)) {
      return false;
    }
    if (first + 1 >= last || *first != '.') {
      return false;
    }
    ++first;
    unsigned int endReps;
    if (!read_int(first, last, endReps)) {
      return false;
    }
    unsigned int idx1;
    unsigned int idx2;
    if (first < last && *first == '.') {
      ++first;
      if (!read_int(first, last, idx1)) {
        return false;
      }
      ++first;
      if (!read_int(first, last, idx2)) {
        return false;
      }
    } else if (VALID_ATIDX(atidx) &&
               mol.getAtomWithIdx(atidx - startAtomIdx)->getDegree() == 2) {
      auto nbrs =
          mol.getAtomNeighbors(mol.getAtomWithIdx(atidx - startAtomIdx));
      idx1 = *nbrs.first;
      nbrs.first++;
      idx2 = *nbrs.first;
    } else if (VALID_ATIDX(atidx)) {
      return false;
    }
    if (first < last && *first == ',') {
      ++first;
    }
    if (VALID_ATIDX(atidx)) {
      if (!accum.empty()) {
        accum += "|";
      }
      accum += (boost::format("%d %d 2 %d %d %d %d") % startReps % endReps %
                (atidx - startAtomIdx + 1) % (idx1 - startAtomIdx + 1) %
                (atidx - startAtomIdx + 1) % (idx2 - startAtomIdx + 1))
                   .str();
    }
  }
  if (!accum.empty()) {
    mol.setProp(common_properties::molFileLinkNodes, accum);
  }
  return true;
}

template <typename Iterator>
bool parse_data_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
                       unsigned int startAtomIdx, unsigned int nSGroups) {
  // these look like: |SgD:2,1:FIELD:info::::|
  // example from CXSMILES docs:
  //    SgD:3,2,1,0:name:data:like:unit:t:(1.,1.)
  // the fields are:
  //    SgD:[atom indices]:[field name]:[data value]:[query
  //    operator]:[unit]:[tag]:[coords]
  //   coords are (-1) if atomic coordinates are present
  if (first >= last || *first != 'S' || first + 3 >= last ||
      *(first + 1) != 'g' || *(first + 2) != 'D' || *(first + 3) != ':') {
    return false;
  }
  first += 4;
  std::vector<unsigned int> atoms;
  if (!read_int_list(first, last, atoms)) {
    return false;
  }
  SubstanceGroup sgroup(&mol, std::string("DAT"));
  sgroup.setProp(cxsmilesindex, nSGroups);
  bool keepSGroup = false;
  for (auto idx : atoms) {
    if (VALID_ATIDX(idx)) {
      keepSGroup = true;
      sgroup.addAtomWithIdx(idx - startAtomIdx);
    }
  }
  ++first;
  std::string name = read_text_to(first, last, ":");
  ++first;
  if (keepSGroup && !name.empty()) {
    sgroup.setProp("FIELDNAME", name);
  }
  // FIX:
  if (keepSGroup) {
    sgroup.setProp("FIELDDISP", "    0.0000    0.0000    DR    ALL  0       0");
  }

  std::string data = read_text_to(first, last, ":");
  ++first;
  if (!data.empty() && keepSGroup) {
    std::vector<std::string> dataFields = {data};
    sgroup.setProp("DATAFIELDS", dataFields);
  }

  std::string oper = read_text_to(first, last, ":");
  ++first;
  if (!oper.empty() && keepSGroup) {
    sgroup.setProp("QUERYOP", oper);
  }
  std::string unit = read_text_to(first, last, ":");
  ++first;
  if (!unit.empty() && keepSGroup) {
    sgroup.setProp("FIELDINFO", unit);
  }
  std::string tag = read_text_to(first, last, ":");
  ++first;
  if (!tag.empty() && keepSGroup) {
    // not actually part of what ends up in the output, but
    // it is part of CXSMARTS
    sgroup.setProp("FIELDTAG", tag);
  }
  if (first < last && *first == '(') {
    // FIX
    std::string coords = read_text_to(first, last, ")");
    ++first;
    if (keepSGroup) {
      sgroup.setProp("COORDS", coords);
    }
  }
  // the label processing can destroy sgroup info, so do that now
  // (the function will immediately return if already called)
  if (keepSGroup) {
    processCXSmilesLabels(mol);
    sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);
    addSubstanceGroup(mol, sgroup);
  }
  return true;
}

namespace {
std::vector<RDKit::SubstanceGroup>::iterator find_matching_sgroup(
    std::vector<RDKit::SubstanceGroup> &sgs, unsigned int targetId) {
  return std::find_if(sgs.begin(), sgs.end(), [targetId](const auto &sg) {
    unsigned int pval;
    if (sg.getPropIfPresent(cxsmilesindex, pval)) {
      if (pval == targetId) {
        return true;
      }
    }
    return false;
  });
}
}  // namespace
template <typename Iterator>
bool parse_sgroup_hierarchy(Iterator &first, Iterator last, RDKit::RWMol &mol) {
  // these look like: |SgH:1:0|
  // from CXSMILES docs:
  //    SgH:parentSgroupIndex1:childSgroupIndex1.childSgroupIndex2,parentSgroupIndex2:childSgroupIndex1
  if (first >= last || *first != 'S' || first + 3 >= last ||
      *(first + 1) != 'g' || *(first + 2) != 'H' || *(first + 3) != ':') {
    return false;
  }
  first += 4;
  auto &sgs = getSubstanceGroups(mol);
  while (1) {
    unsigned int parentId;
    if (!read_int(first, last, parentId)) {
      return false;
    }

    bool validParent = true;
    auto psg = find_matching_sgroup(sgs, parentId);
    if (psg == sgs.end()) {
      validParent = false;
    } else {
      psg->getPropIfPresent("index", parentId);
    }
    if (first <= last && *first == ':') {
      ++first;
      std::vector<unsigned int> children;
      if (!read_int_list(first, last, children, '.')) {
        return false;
      }
      if (validParent) {
        for (auto childId : children) {
          if (childId >= sgs.size()) {
            throw SmilesParseException(
                "child id references non-existent SGroup");
          }
          auto csg = find_matching_sgroup(sgs, childId);
          if (csg != sgs.end()) {
            unsigned int cid;
            csg->getProp("index", cid);
            csg->setProp("PARENT", parentId);
          }
        }
      }
      if (first <= last && *first == ',') {
        ++first;
      } else {
        break;
      }
    } else {
      return false;
    }
  }

  return true;
}

template <typename Iterator>
bool parse_polymer_sgroup(Iterator &first, Iterator last, RDKit::RWMol &mol,
                          unsigned int startAtomIdx, unsigned int nSGroups) {
  // these look like:
  //    |Sg:n:6,1,2,4::hh&#44;f:6,0,:4,2,|
  // example from CXSMILES docs:
  // the fields are:
  //    Sg:[type]:[atom indices]:[subscript]:[superscript]:[head crossing
  //    bonds]:[tail crossing bonds]:
  //
  // note that it's legit for empty fields to be completely missing.
  //   for example, this doesn't have any crossing bonds indicated:
  // *-CCCN-* |$star_e;;;;;star_e$,Sg:n:4,1,2,3::hh|
  // this last bit makes the whole thing doubleplusfun to parse

  if (first >= last || *first != 'S' || first + 2 >= last ||
      *(first + 1) != 'g' || *(first + 2) != ':') {
    return false;
  }
  first += 3;

  std::string typ = read_text_to(first, last, ":");
  ++first;
  if (sgroupTypemap.find(typ) == sgroupTypemap.end()) {
    return false;
  }
  bool keepSGroup = false;
  SubstanceGroup sgroup(&mol, sgroupTypemap[typ]);
  sgroup.setProp(cxsmilesindex, nSGroups);
  if (typ == "alt") {
    sgroup.setProp("SUBTYPE", std::string("ALT"));
  } else if (typ == "ran") {
    sgroup.setProp("SUBTYPE", std::string("RAN"));
  } else if (typ == "blk") {
    sgroup.setProp("SUBTYPE", std::string("BLO"));
  }

  std::vector<unsigned int> atoms;
  if (!read_int_list(first, last, atoms)) {
    return false;
  }
  //++first;
  for (auto idx : atoms) {
    if (VALID_ATIDX(idx)) {
      sgroup.addAtomWithIdx(idx - startAtomIdx);
      keepSGroup = true;
    }
  }
  std::vector<unsigned int> headCrossing;
  std::vector<unsigned int> tailCrossing;
  if (first <= last && *first == ':') {
    ++first;
    std::string subscript = read_text_to(first, last, ":|");
    if (keepSGroup && !subscript.empty()) {
      sgroup.setProp("LABEL", subscript);
    }
    if (first <= last && *first == ':') {
      ++first;
      std::string superscript = read_text_to(first, last, ":|,");
      if (keepSGroup && !superscript.empty()) {
        sgroup.setProp("CONNECT", superscript);
      }

      if (first <= last && *first == ':') {
        ++first;
        if (!read_int_list(first, last, headCrossing)) {
          return false;
        }
        if (keepSGroup && !headCrossing.empty()) {
          for (auto &cidx : headCrossing) {
            if (VALID_ATIDX(cidx)) {
              cidx -= startAtomIdx;
            } else {
              keepSGroup = false;
              break;
            }
          }
          sgroup.setProp(_headCrossings, headCrossing, true);
        }
        if (first <= last && *first == ':') {
          ++first;
          if (!read_int_list(first, last, tailCrossing)) {
            return false;
          }
        }
        if (keepSGroup && !tailCrossing.empty()) {
          for (auto &cidx : tailCrossing) {
            if (VALID_ATIDX(cidx)) {
              cidx -= startAtomIdx;
            } else {
              keepSGroup = false;
              break;
            }
          }
          sgroup.setProp("_tailCrossings", tailCrossing, true);
        }
      }
    }
  }
  if (keepSGroup) {  // the label processing can destroy sgroup info, so do that
                     // now (the function will immediately return if already
                     // called)
    processCXSmilesLabels(mol);

    finalizePolymerSGroup(mol, sgroup);
    sgroup.setProp<unsigned int>("index", getSubstanceGroups(mol).size() + 1);

    addSubstanceGroup(mol, sgroup);
  }
  return true;
}

template <typename Iterator>
bool parse_variable_attachments(Iterator &first, Iterator last,
                                RDKit::RWMol &mol, unsigned int startAtomIdx) {
  // these look like: CO*.C1=CC=NC=C1 |m:2:3.5.4|
  // that corresponds to replacing the bond to atom 2 with bonds to atom 3, 5,
  // or 4
  //
  if (first >= last || *first != 'm' || first + 1 >= last ||
      *(first + 1) != ':') {
    return false;
  }
  first += 2;

  while (first < last && *first >= '0' && *first <= '9') {
    unsigned int at1idx;
    if (!read_int(first, last, at1idx)) {
      return false;
    }

    if (VALID_ATIDX(at1idx) &&
        mol.getAtomWithIdx(at1idx - startAtomIdx)->getDegree() != 1) {
      BOOST_LOG(rdWarningLog)
          << "position variation bond to atom with more than one bond"
          << std::endl;
      return false;
    }
    if (first < last && *first == ':') {
      ++first;
    } else {
      BOOST_LOG(rdWarningLog) << "improperly formatted m: block" << std::endl;
      return false;
    }
    std::vector<std::string> others;
    while (first < last && *first >= '0' && *first <= '9') {
      unsigned int aidx;
      if (!read_int(first, last, aidx)) {
        return false;
      }
      if (VALID_ATIDX(aidx)) {
        others.push_back(
            (boost::format("%d") % (aidx - startAtomIdx + 1)).str());
      }
      if (first < last && *first == '.') {
        ++first;
      }
    }
    if (VALID_ATIDX(at1idx)) {
      std::string endPts = (boost::format("(%d") % others.size()).str();
      for (auto idx : others) {
        endPts += " " + idx;
      }
      endPts += ")";

      for (auto nbri : boost::make_iterator_range(
               mol.getAtomBonds(mol.getAtomWithIdx(at1idx - startAtomIdx)))) {
        auto bnd = mol[nbri];
        bnd->setProp(common_properties::_MolFileBondEndPts, endPts);
        bnd->setProp(common_properties::_MolFileBondAttach, std::string("ANY"));
      }
    }
    if (first < last && *first == ',') {
      ++first;
    }
  }
  return true;
}
template <typename Iterator>
bool parse_substitution(Iterator &first, Iterator last, RDKit::RWMol &mol,
                        unsigned int startAtomIdx) {
  if (first >= last || *first != 's' || first + 1 >= last ||
      *(first + 1) != ':') {
    return false;
  }
  first += 2;
  while (first < last && *first >= '0' && *first <= '9') {
    unsigned int n1;
    if (!read_int(first, last, n1)) {
      return false;
    }
    // check that we can read at least two more characters:
    if (first + 1 >= last || *first != ':') {
      return false;
    }
    ++first;
    unsigned int n2;
    if (*first == '*') {
      ++first;
      n2 = 0xDEADBEEF;
      if (VALID_ATIDX(n1)) {
        mol.setProp(common_properties::_NeedsQueryScan, 1);
      }
    } else {
      if (!read_int(first, last, n2)) {
        return false;
      }
    }
    if (VALID_ATIDX(n1)) {
      auto atom = mol.getAtomWithIdx(n1 - startAtomIdx);
      if (!atom->hasQuery()) {
        atom = QueryOps::replaceAtomWithQueryAtom(&mol, atom);
      }
      atom->expandQuery(makeAtomNonHydrogenDegreeQuery(n2),
                        Queries::COMPOSITE_AND);
    }
    if (first < last && *first == ',') {
      ++first;
    }
  }
  return true;
}

template <typename Iterator>
bool processRadicalSection(Iterator &first, Iterator last, RDKit::RWMol &mol,
                           unsigned int numRadicalElectrons,
                           unsigned int startAtomIdx) {
  if (first >= last) {
    return false;
  }
  ++first;
  if (first >= last || *first != ':') {
    return false;
  }
  ++first;
  unsigned int atIdx;
  if (!read_int(first, last, atIdx)) {
    return false;
  }
  if (VALID_ATIDX(atIdx)) {
    mol.getAtomWithIdx(atIdx - startAtomIdx)
        ->setNumRadicalElectrons(numRadicalElectrons);
  }
  while (first < last && *first == ',') {
    ++first;
    if (first < last && (*first < '0' || *first > '9')) {
      return true;
    }
    if (!read_int(first, last, atIdx)) {
      return false;
    }
    if (VALID_ATIDX(atIdx)) {
      mol.getAtomWithIdx(atIdx - startAtomIdx)
          ->setNumRadicalElectrons(numRadicalElectrons);
    }
  }
  return first < last;
}

template <typename Iterator>
bool parse_radicals(Iterator &first, Iterator last, RDKit::RWMol &mol,
                    unsigned int startAtomIdx) {
  if (first >= last || *first != '^') {
    return false;
  }
  while (*first == '^') {
    ++first;
    if (first >= last) {
      return false;
    }
    if (*first < '1' || *first > '7') {
      return false;  // these are the values that are allowed to be there
    }
    switch (*first) {
      case '1':
        if (!processRadicalSection(first, last, mol, 1, startAtomIdx)) {
          return false;
        }
        break;
      case '2':
      case '3':
      case '4':
        if (!processRadicalSection(first, last, mol, 2, startAtomIdx)) {
          return false;
        }
        break;
      case '5':
      case '6':
      case '7':
        if (!processRadicalSection(first, last, mol, 3, startAtomIdx)) {
          return false;
        }
        break;
      default:
        BOOST_LOG(rdWarningLog)
            << "Radical specification " << *first << " ignored.";
    }
  }
  return true;
}

template <typename Iterator>
bool parse_enhanced_stereo(Iterator &first, Iterator last, RDKit::RWMol &mol,
                           unsigned int startAtomIdx) {
  StereoGroupType group_type = StereoGroupType::STEREO_ABSOLUTE;
  if (*first == 'a') {
    group_type = StereoGroupType::STEREO_ABSOLUTE;
  } else if (*first == 'o') {
    group_type = StereoGroupType::STEREO_OR;
  } else if (*first == '&') {
    group_type = StereoGroupType::STEREO_AND;
  }
  ++first;

  // OR and AND groups carry a group number
  if (group_type != StereoGroupType::STEREO_ABSOLUTE) {
    unsigned int group_id = 0;
    read_int(first, last, group_id);
  }

  if (first >= last || *first != ':') {
    return false;
  }
  ++first;

  std::vector<Atom *> atoms;
  while (first <= last && *first >= '0' && *first <= '9') {
    unsigned int aidx;
    if (read_int(first, last, aidx)) {
      if (VALID_ATIDX(aidx)) {
        Atom *atom = mol.getAtomWithIdx(aidx - startAtomIdx);
        if (!atom) {
          BOOST_LOG(rdWarningLog)
              << "Atom " << aidx << " not found!" << std::endl;
          return false;
        }
        atoms.push_back(atom);
      }
    } else {
      return false;
    }

    if (first < last && *first == ',') {
      ++first;
    }
  }
  if (!atoms.empty()) {
    std::vector<StereoGroup> mol_stereo_groups(mol.getStereoGroups());
    mol_stereo_groups.emplace_back(group_type, std::move(atoms));
    mol.setStereoGroups(std::move(mol_stereo_groups));
  }

  return true;
}

template <typename Iterator>
bool parse_it(Iterator &first, Iterator last, RDKit::RWMol &mol,
              unsigned int startAtomIdx, unsigned int startBondIdx) {
  if (first >= last || *first != '|') {
    return false;
  }
  ++first;
  unsigned int nSGroups = 0;
  while (first < last && *first != '|') {
    typename Iterator::difference_type length = std::distance(first, last);
    if (*first == '(') {
      if (!parse_coords(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == '$') {
      if (length > 4 && *(first + 1) == '_' && *(first + 2) == 'A' &&
          *(first + 3) == 'V' && *(first + 4) == ':') {
        first += 4;
        if (!parse_atom_values(first, last, mol, startAtomIdx)) {
          return false;
        }
      } else {
        if (!parse_atom_labels(first, last, mol, startAtomIdx)) {
          return false;
        }
      }
    } else if (length > 9 && std::string(first, first + 9) == "atomProp:") {
      first += 9;
      if (!parse_atom_props(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 'C') {
      if (!parse_coordinate_bonds(first, last, mol, Bond::DATIVE, startAtomIdx,
                                  startBondIdx)) {
        return false;
      }
    } else if (*first == 'H') {
      if (!parse_coordinate_bonds(first, last, mol, Bond::HYDROGEN,
                                  startAtomIdx, startBondIdx)) {
        return false;
      }
    } else if (*first == '^') {
      if (!parse_radicals(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 'a' || *first == 'o' ||
               (*first == '&' && first + 1 < last && first[1] != '#')) {
      if (!parse_enhanced_stereo(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 'r' && first + 1 < last && first[1] == 'b') {
      if (!parse_ring_bonds(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 'L' && first + 1 < last && first[1] == 'N') {
      if (!parse_linknodes(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
               first[2] == 'D') {
      if (!parse_data_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
        return false;
      }
    } else if (*first == 'S' && first + 2 < last && first[1] == 'g' &&
               first[2] == 'H') {
      if (!parse_sgroup_hierarchy(first, last, mol)) {
        return false;
      }
    } else if (*first == 'S' && first + 1 < last && first[1] == 'g') {
      if (!parse_polymer_sgroup(first, last, mol, startAtomIdx, nSGroups++)) {
        return false;
      }
    } else if (*first == 'u') {
      if (!parse_unsaturation(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 's') {
      if (!parse_substitution(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else if (*first == 'm') {
      if (!parse_variable_attachments(first, last, mol, startAtomIdx)) {
        return false;
      }
    } else {
      ++first;
    }
    // if(first < last && *first != '|') ++first;
  }
  if (first >= last || *first != '|') {
    return false;
  }
  ++first;  // step past the last '|'
  return true;
}
}  // namespace parser

void parseCXExtensions(RDKit::RWMol &mol, const std::string &extText,
                       std::string::const_iterator &first,
                       unsigned int startAtomIdx, unsigned int startBondIdx) {
  // BOOST_LOG(rdWarningLog) << "parseCXNExtensions: " << extText << std::endl;
  if (extText.empty()) {
    return;
  }
  if (extText[0] != '|') {
    throw RDKit::SmilesParseException(
        "CXSMILES extension does not start with |");
  }
  first = extText.begin();
  bool ok =
      parser::parse_it(first, extText.end(), mol, startAtomIdx, startBondIdx);
  if (!ok) {
    throw RDKit::SmilesParseException("failure parsing CXSMILES extensions");
  }
  processCXSmilesLabels(mol);
  mol.clearProp("_cxsmilesLabelsProcessed");
}
}  // end of namespace SmilesParseOps

namespace RDKit {
namespace SmilesWrite {
namespace {
std::string quote_string(const std::string &txt) {
  // FIX
  return txt;
}

std::string quote_atomprop_string(const std::string &txt) {
  // at a bare minimum, . needs to be escaped
  std::string res;
  for(auto c : txt) {
    if(c == '.') {
      res += "&#46;";
    } else {
      res += c;
    }
  }
  return res;
}

std::string get_enhanced_stereo_block(
    const ROMol &mol, const std::vector<unsigned int> &atomOrder) {
  if (mol.getStereoGroups().empty()) {
    return "";
  }
  std::stringstream res;
  // we need a map from original atom idx to output idx:
  std::vector<unsigned int> revOrder(mol.getNumAtoms());
  for (unsigned i = 0; i < atomOrder.size(); ++i) {
    revOrder[atomOrder[i]] = i;
  }
  std::vector<unsigned int> absAts;
  std::vector<std::vector<unsigned int>> orGps;
  std::vector<std::vector<unsigned int>> andGps;

  // we want this to be canonical (future proofing)
  for (const auto &sg : mol.getStereoGroups()) {
    std::vector<unsigned int> aids;
    aids.reserve(sg.getAtoms().size());
    for (const auto at : sg.getAtoms()) {
      aids.push_back(revOrder[at->getIdx()]);
    }
    switch (sg.getGroupType()) {
      case StereoGroupType::STEREO_ABSOLUTE:
        absAts.insert(absAts.end(), aids.begin(), aids.end());
        break;
      case StereoGroupType::STEREO_OR:
        std::sort(aids.begin(), aids.end());
        orGps.push_back(aids);
        break;
      case StereoGroupType::STEREO_AND:
        std::sort(aids.begin(), aids.end());
        andGps.push_back(aids);
        break;
    }
  }
  if (!absAts.empty()) {
    res << "a:";
    std::sort(absAts.begin(), absAts.end());
    for (auto aid : absAts) {
      res << aid << ",";
    }
  }
  if (!orGps.empty()) {
    std::sort(orGps.begin(), orGps.end());
    unsigned int gIdx = 1;
    for (const auto &gp : orGps) {
      res << "o" << gIdx++ << ":";
      for (auto aid : gp) {
        res << aid << ",";
      }
    }
  }
  if (!andGps.empty()) {
    std::sort(andGps.begin(), andGps.end());
    unsigned int gIdx = 1;
    for (const auto &gp : andGps) {
      res << "&" << gIdx++ << ":";
      for (auto aid : gp) {
        res << aid << ",";
      }
    }
  }
  std::string resStr = res.str();
  if (!resStr.empty() && resStr.back() == ',') {
    resStr.pop_back();
  }
  return resStr;
}

std::string get_sgroup_hierarchy_block(const ROMol &mol) {
  const auto &sgs = getSubstanceGroups(mol);
  if (sgs.empty()) {
    return "";
  }
  std::stringstream res;
  // we need a map from sgroup index to output index;
  std::map<unsigned int, unsigned int> sgroupOrder;
  bool parentPresent = false;
  for (const auto &sg : sgs) {
    if (sg.hasProp("_cxsmilesOutputIndex")) {
      unsigned int sgidx = sg.getIndexInMol();
      sg.getPropIfPresent("index", sgidx);
      sgroupOrder[sgidx] = sg.getProp<unsigned int>("_cxsmilesOutputIndex");
      sg.clearProp("_cxsmilesOutputIndex");
    }
    if (sg.hasProp("PARENT")) {
      parentPresent = true;
    }
  }

  if (parentPresent) {
    // now loop over them and add the information
    std::map<unsigned int, std::vector<unsigned int>> accum;
    for (const auto &sg : sgs) {
      unsigned pidx;
      if (sg.getPropIfPresent("PARENT", pidx) &&
          sgroupOrder.find(pidx) != sgroupOrder.end()) {
        unsigned int sgidx = sg.getIndexInMol();
        sg.getPropIfPresent("index", sgidx);
        if (sgroupOrder.find(sgidx) != sgroupOrder.end()) {
          accum[sgroupOrder[pidx]].push_back(sgroupOrder[sgidx]);
        }
      }
    }
    if (!accum.empty()) {
      res << "SgH:";
      for (const auto &pr : accum) {
        res << pr.first << ":";
        for (auto v : pr.second) {
          res << v << ".";
        }
        // remove the extra ".":
        res.seekp(-1, res.cur);
        res << ",";
      }
    }
    std::string resStr = res.str();
    while (!resStr.empty() && resStr.back() == ',') {
      resStr.pop_back();
    }
    return resStr;
  } else {
    return "";
  }
}

std::string get_sgroup_polymer_block(
    const ROMol &mol, const std::vector<unsigned int> &atomOrder,
    const std::vector<unsigned int> &bondOrder) {
  const auto &sgs = getSubstanceGroups(mol);
  if (sgs.empty()) {
    return "";
  }
  unsigned int sgroupOutputIndex = 0;
  mol.getPropIfPresent("_cxsmilesOutputIndex", sgroupOutputIndex);
  std::stringstream res;
  // we need a map from original atom idx to output idx:
  std::vector<unsigned int> revAtomOrder(mol.getNumAtoms());
  for (unsigned i = 0; i < atomOrder.size(); ++i) {
    revAtomOrder[atomOrder[i]] = i;
  }
  // we need a map from original bond idx to output idx:
  std::vector<unsigned int> revBondOrder(mol.getNumBonds());
  for (unsigned i = 0; i < bondOrder.size(); ++i) {
    revBondOrder[bondOrder[i]] = i;
  }

  std::map<std::string, std::string> reverseTypemap;
  for (const auto &pr : SmilesParseOps::sgroupTypemap) {
    if (reverseTypemap.find(pr.second) == reverseTypemap.end()) {
      reverseTypemap[pr.second] = pr.first;
    }
  }

  for (const auto &sg : sgs) {
    std::string typ;
    if (sg.getPropIfPresent("TYPE", typ) &&
        reverseTypemap.find(typ) != reverseTypemap.end()) {
      sg.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
      sgroupOutputIndex++;

      res << "Sg:";
      std::string subtype;
      if (typ == "COP" && sg.getPropIfPresent("SUBTYPE", subtype)) {
        if (subtype == "ALT") {
          res << "alt";
        } else if (subtype == "RAN") {
          res << "ran";
        } else if (subtype == "BLO") {
          res << "blk";
        } else {
          res << reverseTypemap["COP"];
        }
      } else {
        res << reverseTypemap[typ];
      }
      res << ":";
      for (const auto oaid : sg.getAtoms()) {
        res << revAtomOrder[oaid] << ",";
      }
      // remove the extra ",":
      res.seekp(-1, res.cur);
      res << ":";
      std::string label;
      if (sg.getPropIfPresent("LABEL", label)) {
        res << label;
      }
      res << ":";
      std::string connect;
      if (sg.getPropIfPresent("CONNECT", connect)) {
        boost::algorithm::to_lower(connect);
        res << connect;
      }
      res << ":";
      std::vector<unsigned int> headCrossings;
      if (sg.getPropIfPresent("XBHEAD", headCrossings) &&
          headCrossings.size() > 1) {
        for (auto v : headCrossings) {
          res << bondOrder[v] << ",";
        }
        // remove the extra ",":
        res.seekp(-1, res.cur);
      }
      res << ":";
      std::vector<unsigned int> tailCrossings;
      if (sg.getPropIfPresent("XBCORR", tailCrossings) &&
          tailCrossings.size() > 2) {
        for (unsigned int i = 1; i < tailCrossings.size(); i += 2) {
          res << bondOrder[tailCrossings[i]] << ",";
        }
        // remove the extra ",":
        res.seekp(-1, res.cur);
      }
      res << ":";
    }
    res << ",";
  }

  std::string resStr = res.str();
  while (!resStr.empty() && resStr.back() == ',') {
    resStr.pop_back();
  }
  mol.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);

  return resStr;
}

std::string get_sgroup_data_block(const ROMol &mol,
                                  const std::vector<unsigned int> &atomOrder) {
  const auto &sgs = getSubstanceGroups(mol);
  if (sgs.empty()) {
    return "";
  }

  unsigned int sgroupOutputIndex = 0;
  mol.getPropIfPresent("_cxsmilesOutputIndex", sgroupOutputIndex);

  std::stringstream res;
  // we need a map from original atom idx to output idx:
  std::vector<unsigned int> revOrder(mol.getNumAtoms());
  for (unsigned i = 0; i < atomOrder.size(); ++i) {
    revOrder[atomOrder[i]] = i;
  }

  for (const auto &sg : sgs) {
    if (sg.hasProp("TYPE") && sg.getProp<std::string>("TYPE") == "DAT") {
      sg.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);
      sgroupOutputIndex++;

      res << "SgD:";
      // we don't attempt to canonicalize the atom order because the user
      // may ascribe some significance to the ordering of the atoms
      for (const auto oaid : sg.getAtoms()) {
        res << revOrder[oaid] << ",";
      }
      // remove the extra ",":
      res.seekp(-1, res.cur);
      res << ":";
      std::string prop;
      if (sg.getPropIfPresent("FIELDNAME", prop) && !prop.empty()) {
        res << prop;
      }
      res << ":";
      std::vector<std::string> vprop;
      if (sg.getPropIfPresent("DATAFIELDS", vprop) && !vprop.empty()) {
        for (const auto &pv : vprop) {
          res << pv << ",";
        }
        // remove the extra ",":
        res.seekp(-1, res.cur);
      }
      res << ":";
      if (sg.getPropIfPresent("QUERYOP", prop) && !prop.empty()) {
        res << prop;
      }
      res << ":";
      if (sg.getPropIfPresent("FIELDINFO", prop) && !prop.empty()) {
        res << prop;
      }
      res << ":";
      if (sg.getPropIfPresent("FIELDTAG", prop) && !prop.empty()) {
        res << prop;
      }
      res << ":";
      // FIX: do something about the coordinates
    }
    res << ",";
  }

  std::string resStr = res.str();
  if (!resStr.empty() && resStr.back() == ',') {
    resStr.pop_back();
  }
  mol.setProp("_cxsmilesOutputIndex", sgroupOutputIndex);

  return resStr;
}

std::string get_atomlabel_block(const ROMol &mol,
                                const std::vector<unsigned int> &atomOrder) {
  std::string res = "";
  for (auto idx : atomOrder) {
    if (idx != atomOrder.front()) {
      res += ";";
    }
    std::string lbl;
    const auto atom = mol.getAtomWithIdx(idx);
    if (atom->getPropIfPresent(common_properties::_QueryAtomGenericLabel,
                               lbl)) {
      res += quote_string(lbl + "_p");
    } else if (!atom->getAtomicNum() &&
               atom->getPropIfPresent(common_properties::dummyLabel, lbl) &&
               std::find(SmilesParseOps::pseudoatoms.begin(),
                         SmilesParseOps::pseudoatoms.end(),
                         lbl) != SmilesParseOps::pseudoatoms.end()) {
      res += quote_string(lbl + "_p");
    } else if (atom->getPropIfPresent(common_properties::atomLabel, lbl)) {
      res += quote_string(lbl);
    }
  }
  // if we didn't find anything return an empty string
  if (std::find_if_not(res.begin(), res.end(),
                       [](const auto c) { return c == ';'; }) == res.end()) {
    res.clear();
  }
  return res;
}

std::string get_value_block(const ROMol &mol,
                            const std::vector<unsigned int> &atomOrder,
                            const std::string &prop) {
  std::string res = "";
  bool first = true;
  for (auto idx : atomOrder) {
    if (!first) {
      res += ";";
    } else {
      first = false;
    }
    std::string lbl;
    if (mol.getAtomWithIdx(idx)->getPropIfPresent(prop, lbl)) {
      res += quote_string(lbl);
    }
  }
  return res;
}
std::string get_radical_block(const ROMol &mol,
                              const std::vector<unsigned int> &atomOrder) {
  std::string res = "";
  std::map<unsigned int, std::vector<unsigned int>> rads;
  for (unsigned int i = 0; i < atomOrder.size(); ++i) {
    auto idx = atomOrder[i];
    auto nrad = mol.getAtomWithIdx(idx)->getNumRadicalElectrons();
    if (nrad) {
      rads[nrad].push_back(i);
    }
  }
  if (rads.size()) {
    for (const auto &pr : rads) {
      switch (pr.first) {
        case 1:
          res += "^1:";
          break;
        case 2:
          res += "^2:";
          break;
        case 3:
          res += "^5:";
          break;
        default:
          BOOST_LOG(rdWarningLog) << "unsupported number of radical electrons "
                                  << pr.first << std::endl;
      }
      for (auto aidx : pr.second) {
        res += boost::str(boost::format("%d,") % aidx);
      }
    }
  }
  return res;
}
double zero_small_vals(double val) {
  if (fabs(val) < 1e-4) {
    return 0.0;
  }
  return val;
}
std::string get_coords_block(const ROMol &mol,
                             const std::vector<unsigned int> &atomOrder) {
  std::string res = "";
  const auto &conf = mol.getConformer();
  bool first = true;
  for (auto idx : atomOrder) {
    const auto &pt = conf.getAtomPos(idx);
    if (!first) {
      res += ";";
    } else {
      first = false;
    }
    res += boost::str(boost::format("%g,%g,") % zero_small_vals(pt.x) %
                      zero_small_vals(pt.y));
    if (conf.is3D()) {
      auto zc = boost::str(boost::format("%g") % zero_small_vals(pt.z));
      if (zc != "0") {
        res += zc;
      }
    }
  }
  return res;
}

std::string get_atom_props_block(const ROMol &mol,
                                 const std::vector<unsigned int> &atomOrder) {
  std::vector<std::string> skip = {common_properties::atomLabel,
                                   common_properties::molFileValue,
                                   common_properties::molParity};
  std::string res = "";
  unsigned int which = 0;
  for (auto idx : atomOrder) {
    const auto atom = mol.getAtomWithIdx(idx);
    bool includePrivate = false, includeComputed = false;
    for (const auto &pn : atom->getPropList(includePrivate, includeComputed)) {
      if (std::find(skip.begin(), skip.end(), pn) == skip.end()) {
        std::string pv = atom->getProp<std::string>(pn);
        if (pn == "dummyLabel" &&
            std::find(SmilesParseOps::pseudoatoms.begin(),
                      SmilesParseOps::pseudoatoms.end(),
                      pv) != SmilesParseOps::pseudoatoms.end()) {
          // it's a pseudoatom, skip it
          continue;
        }
        if (res.size() == 0) {
          res += "atomProp";
        }
        res += boost::str(boost::format(":%d.%s.%s") % which %
                          quote_atomprop_string(pn) % quote_atomprop_string(pv));
      }
    }
    ++which;
  }
  return res;
}

std::string get_linknodes_block(const ROMol &mol,
                                const std::vector<unsigned int> &atomOrder) {
  bool strict = false;
  auto linkNodes = MolEnumerator::utils::getMolLinkNodes(mol, strict);
  if (linkNodes.empty()) {
    return "";
  }
  // we need a map from original atom idx to output idx:
  std::vector<unsigned int> revOrder(mol.getNumAtoms());
  for (unsigned i = 0; i < atomOrder.size(); ++i) {
    revOrder[atomOrder[i]] = i;
  }

  std::stringstream res;
  res << "LN:";
  for (const auto &ln : linkNodes) {
    unsigned int atomIdx = atomOrder[ln.bondAtoms[0].first];
    res << atomIdx << ":" << ln.minRep << "." << ln.maxRep;
    if (mol.getAtomWithIdx(ln.bondAtoms[0].first)->getDegree() > 2) {
      // include the outer atom indices
      res << "." << atomOrder[ln.bondAtoms[0].second] << "."
          << atomOrder[ln.bondAtoms[1].second];
    }
    res << ",";
  }

  std::string resStr = res.str();
  if (!resStr.empty() && resStr.back() == ',') {
    resStr.pop_back();
  }
  return resStr;
}

void appendToCXExtension(const std::string &addition, std::string &base) {
  if (!addition.empty()) {
    if (base.size() > 1) {
      base += ",";
    }
    base += addition;
  }
}

}  // namespace
std::string getCXExtensions(const ROMol &mol, std::uint32_t flags) {
  std::string res = "|";
  // we will need atom and bond orderings. Get them now:
  const std::vector<unsigned int> &atomOrder =
      mol.getProp<std::vector<unsigned int>>(
          common_properties::_smilesAtomOutputOrder);
  const std::vector<unsigned int> &bondOrder =
      mol.getProp<std::vector<unsigned int>>(
          common_properties::_smilesBondOutputOrder);
  bool needLabels = false;
  bool needValues = false;
  for (auto idx : atomOrder) {
    const auto at = mol.getAtomWithIdx(idx);
    if (at->hasProp(common_properties::atomLabel) ||
        at->hasProp(common_properties::_QueryAtomGenericLabel) ||
        at->hasProp(common_properties::dummyLabel)) {
      needLabels = true;
    }
    if (at->hasProp(common_properties::molFileValue)) {
      needValues = true;
    }
  }
  if ((flags & SmilesWrite::CXSmilesFields::CX_COORDS) &&
      mol.getNumConformers()) {
    res += "(" + get_coords_block(mol, atomOrder) + ")";
  }
  if ((flags & SmilesWrite::CXSmilesFields::CX_ATOM_LABELS) && needLabels) {
    auto lbls = get_atomlabel_block(mol, atomOrder);
    if (!lbls.empty()) {
      if (res.size() > 1) {
        res += ",";
      }
      res += "$" + lbls + "$";
    }
  }
  if ((flags & SmilesWrite::CXSmilesFields::CX_MOLFILE_VALUES) && needValues) {
    if (res.size() > 1) {
      res += ",";
    }
    res += "$_AV:" +
           get_value_block(mol, atomOrder, common_properties::molFileValue) +
           "$";
  }
  auto radblock = get_radical_block(mol, atomOrder);
  if ((flags & SmilesWrite::CXSmilesFields::CX_RADICALS) && radblock.size()) {
    if (res.size() > 1) {
      res += ",";
    }
    res += radblock;
    if (res.back() == ',') {
      res.erase(res.size() - 1);
    }
  }

  if (flags & SmilesWrite::CXSmilesFields::CX_ATOM_PROPS) {
    const auto atomblock = get_atom_props_block(mol, atomOrder);
    appendToCXExtension(atomblock, res);
  }

  if (flags & SmilesWrite::CXSmilesFields::CX_LINKNODES) {
    const auto linknodeblock = get_linknodes_block(mol, atomOrder);
    appendToCXExtension(linknodeblock, res);
  }
  if (flags & SmilesWrite::CXSmilesFields::CX_ENHANCEDSTEREO) {
    const auto stereoblock = get_enhanced_stereo_block(mol, atomOrder);
    appendToCXExtension(stereoblock, res);
  }
  if (flags & SmilesWrite::CXSmilesFields::CX_SGROUPS) {
    const auto sgroupdatablock = get_sgroup_data_block(mol, atomOrder);
    appendToCXExtension(sgroupdatablock, res);
  }
  if (flags & SmilesWrite::CXSmilesFields::CX_POLYMER) {
    const auto sgrouppolyblock =
        get_sgroup_polymer_block(mol, atomOrder, bondOrder);
    appendToCXExtension(sgrouppolyblock, res);
  }
  if (flags & (SmilesWrite::CXSmilesFields::CX_SGROUPS |
               SmilesWrite::CXSmilesFields::CX_POLYMER)) {
    const auto sgrouphierarchyblock = get_sgroup_hierarchy_block(mol);
    appendToCXExtension(sgrouphierarchyblock, res);
  }
  mol.clearProp("_cxsmilesOutputIndex");
  if (res.size() > 1) {
    res += "|";
  } else {
    res = "";
  }
  return res;
}
}  // namespace SmilesWrite
}  // namespace RDKit
