//
//  Copyright (C) 2002-2019 Greg Landrum and Rational Discovery LLC
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "SmartsWrite.h"
#include <sstream>
#include <cstdint>
#include <boost/algorithm/string.hpp>
#include "SmilesWrite.h"
#include <GraphMol/RDKitBase.h>
#include <GraphMol/RDKitQueries.h>
#include <GraphMol/Canon.h>
#include <GraphMol/new_canon.h>
#include <RDGeneral/Exceptions.h>
#include <RDGeneral/RDLog.h>

namespace RDKit {
using namespace Canon;

// local utility namespace
namespace {

enum class QueryBoolFeatures {
  HAS_AND = 0x1,
  HAS_LOWAND = 0x2,
  HAS_OR = 0x4,
  HAS_RECURSION = 0x8
};

std::string _recurseGetSmarts(const QueryAtom *qatom,
                              const QueryAtom::QUERYATOM_QUERY *node,
                              bool negate, unsigned int &features);
std::string _recurseBondSmarts(const Bond *bond,
                               const QueryBond::QUERYBOND_QUERY *node,
                               bool negate, int atomToLeftIdx,
                               unsigned int &features);

std::string _combineChildSmarts(std::string cs1, unsigned int features1,
                                std::string cs2, unsigned int features2,
                                std::string descrip, unsigned int &features) {
  std::string res = "";
  if ((descrip.find("Or") > 0) && (descrip.find("Or") < descrip.length())) {
    // if either of child smarts already have a "," and ";" we can't have one
    // more OR here
    if ((features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND) &&
         features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR)) ||
        (features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND) &&
         features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR))) {
      throw ValueErrorException(
          "This is a non-smartable query - OR above and below AND in the "
          "binary tree");
    }
    res += cs1;
    if (!(cs1.empty() || cs2.empty())) {
      res += ",";
    }
    res += cs2;
    features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_OR);
  } else if ((descrip.find("And") > 0) &&
             (descrip.find("And") < descrip.length())) {
    std::string symb;
    if (features1 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR) ||
        features2 & static_cast<unsigned int>(QueryBoolFeatures::HAS_OR)) {
      symb = ";";
      features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_LOWAND);
    } else {
      symb = "&";
      features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_AND);
    }
    res += cs1;
    if (!(cs1.empty() || cs2.empty())) {
      res += symb;
    }
    res += cs2;
  } else {
    std::stringstream err;
    err << "Don't know how to combine using " << descrip;
    throw ValueErrorException(err.str());
  }
  features |= features1;
  features |= features2;

  return res;
}  // namespace

template <typename T>
void describeQuery(const T *query, std::string leader = "\t") {
  // BOOST_LOG(rdInfoLog) << leader << query->getDescription() << std::endl;
  typename T::CHILD_VECT_CI iter;
  for (iter = query->beginChildren(); iter != query->endChildren(); ++iter) {
    describeQuery(iter->get(), leader + "\t");
  }
}

const static std::string _qatomHasStereoSet = "_qatomHasStereoSet";
std::string getAtomSmartsSimple(const QueryAtom *qatom,
                                const ATOM_EQUALS_QUERY *query,
                                bool &needParen) {
  PRECONDITION(query, "bad query");

  std::string descrip = query->getDescription();
  bool hasVal = false;
  enum class Modifiers : std::uint8_t { NONE, RANGE, LESS, GREATER };
  Modifiers mods = Modifiers::NONE;
  if (boost::starts_with(descrip, "range_")) {
    mods = Modifiers::RANGE;
    descrip = descrip.substr(6);
  } else if (boost::starts_with(descrip, "less_")) {
    mods = Modifiers::LESS;
    descrip = descrip.substr(5);
  } else if (boost::starts_with(descrip, "greater_")) {
    mods = Modifiers::GREATER;
    descrip = descrip.substr(8);
  }
  std::stringstream res;
  if (descrip == "AtomImplicitHCount") {
    res << "h";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomHasImplicitH") {
    res << "h";
    needParen = true;
  } else if (descrip == "AtomTotalValence") {
    res << "v";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomAtomicNum") {
    res << "#";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomExplicitDegree") {
    res << "D";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomNonHydrogenDegree") {
    res << "d";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomTotalDegree") {
    res << "X";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomHasRingBond") {
    res << "x";
    needParen = true;
  } else if (descrip == "AtomHCount") {
    res << "H";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomIsAliphatic") {
    res << "A";
    needParen = false;
  } else if (descrip == "AtomIsAromatic") {
    res << "a";
    needParen = false;
  } else if (descrip == "AtomNull") {
    res << "*";
    needParen = false;
  } else if (descrip == "AtomInRing") {
    res << "R";
    needParen = true;
  } else if (descrip == "AtomMinRingSize") {
    res << "r";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomInNRings") {
    res << "R";
    if (mods == Modifiers::NONE && query->getVal() >= 0) {
      hasVal = true;
    }
    needParen = true;
  } else if (descrip == "AtomHasHeteroatomNeighbors") {
    res << "z";
    needParen = true;
  } else if (descrip == "AtomNumHeteroatomNeighbors") {
    res << "z";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomHasAliphaticHeteroatomNeighbors") {
    res << "Z";
    needParen = true;
  } else if (descrip == "AtomNumAliphaticHeteroatomNeighbors") {
    res << "Z";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomFormalCharge") {
    int val = query->getVal();
    if (val < 0) {
      res << "-";
    } else {
      res << "+";
    }
    if (mods == Modifiers::NONE && abs(val) != 1) {
      res << abs(val);
    }
    needParen = true;
  } else if (descrip == "AtomNegativeFormalCharge") {
    int val = query->getVal();
    if (val < 0) {
      res << "+";
    } else {
      res << "-";
    }
    if (mods == Modifiers::NONE && abs(val) != 1) {
      res << abs(val);
    }
    needParen = true;
  } else if (descrip == "AtomHybridization") {
    res << "^";
    switch (query->getVal()) {
      case Atom::S:
        res << "0";
        break;
      case Atom::SP:
        res << "1";
        break;
      case Atom::SP2:
        res << "2";
        break;
      case Atom::SP3:
        res << "3";
        break;
      case Atom::SP3D:
        res << "4";
        break;
      case Atom::SP3D2:
        res << "5";
        break;
    }
    needParen = true;
  } else if (descrip == "AtomMass") {
    res << query->getVal() / massIntegerConversionFactor << "*";
    needParen = true;
  } else if (descrip == "AtomIsotope") {
    res << query->getVal() << "*";
    needParen = true;
  } else if (descrip == "AtomRingBondCount") {
    res << "x";
    hasVal = true;
    needParen = true;
  } else if (descrip == "AtomUnsaturated") {
    res << "$(*=,:,#*)";
    needParen = true;
  } else if (descrip == "AtomType") {
    int atNum;
    bool isAromatic;
    parseAtomType(query->getVal(), atNum, isAromatic);
    std::string symbol = PeriodicTable::getTable()->getElementSymbol(atNum);
    if (isAromatic) {
      symbol[0] += ('a' - 'A');
    }
    res << symbol;
    if (!SmilesWrite::inOrganicSubset(atNum)) {
      needParen = true;
    }
  } else {
    BOOST_LOG(rdWarningLog)
        << "Cannot write SMARTS for query type : " << descrip
        << ". Ignoring it." << std::endl;
    res << "*";
  }

  if (mods != Modifiers::NONE) {
    res << "{";
    switch (mods) {
      case Modifiers::LESS:
        res << ((const ATOM_LESSEQUAL_QUERY *)query)->getVal() << "-";
        break;
      case Modifiers::RANGE:
        res << ((const ATOM_RANGE_QUERY *)query)->getLower() << "-"
            << ((const ATOM_RANGE_QUERY *)query)->getUpper();
        break;
      case Modifiers::GREATER:
        res << "-" << ((const ATOM_GREATEREQUAL_QUERY *)query)->getVal();
        break;
      default:
        break;
    }
    res << "}";
  } else if (hasVal) {
    res << query->getVal();
  }

  // handle atomic stereochemistry
  if (qatom->hasOwningMol() &&
      qatom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    if (qatom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        !qatom->hasProp(_qatomHasStereoSet) &&
        !qatom->hasProp(common_properties::_brokenChirality)) {
      qatom->setProp(_qatomHasStereoSet, 1);
      switch (qatom->getChiralTag()) {
        case Atom::CHI_TETRAHEDRAL_CW:
          res << "@@";
          needParen = true;
          break;
        case Atom::CHI_TETRAHEDRAL_CCW:
          res << "@";
          needParen = true;
          break;
        default:
          break;
      }
    }
  }

  return res.str();
}

std::string getRecursiveStructureQuerySmarts(
    const QueryAtom::QUERYATOM_QUERY *query) {
  PRECONDITION(query, "bad query");
  PRECONDITION(query->getDescription() == "RecursiveStructure", "bad query");
  const auto *rquery = static_cast<const RecursiveStructureQuery *>(query);
  auto *qmol = const_cast<ROMol *>(rquery->getQueryMol());
  std::string res = MolToSmarts(*qmol);
  res = "$(" + res + ")";
  if (rquery->getNegation()) {
    res = "!" + res;
  }
  return res;
}

std::string getBasicBondRepr(Bond::BondType typ, Bond::BondDir dir,
                             bool doIsomericSmiles, bool reverseDative) {
  std::string res;
  switch (typ) {
    case Bond::SINGLE:
      res = "-";
      if (doIsomericSmiles) {
        if (dir == Bond::ENDDOWNRIGHT) {
          res = "\\";
        } else if (dir == Bond::ENDUPRIGHT) {
          res = "/";
        }
      }
      break;
    case Bond::DOUBLE:
      res = "=";
      break;
    case Bond::TRIPLE:
      res = "#";
      break;
    case Bond::AROMATIC:
      res = ":";
      break;
    case Bond::DATIVE:
      if (reverseDative) {
        res = "<-";
      } else {
        res = "->";
      }
      break;
    default:
      res = "";
  }
  return res;
}  // namespace

std::string getBondSmartsSimple(const Bond *bond,
                                const BOND_EQUALS_QUERY *bquery,
                                int atomToLeftIdx) {
  PRECONDITION(bond, "bad bond");

  PRECONDITION(bquery, "bad query");
  std::string descrip = bquery->getDescription();
  std::string res = "";
  if (descrip == "BondNull") {
    res += "~";
  } else if (descrip == "BondInRing") {
    res += "@";
  } else if (descrip == "SingleOrAromaticBond") {
    // don't need to do anything here... :-)
  } else if (descrip == "SingleOrDoubleOrAromaticBond") {
    res += "-,=,:";
  } else if (descrip == "BondDir") {
    int val = bquery->getVal();
    if (val == static_cast<int>(Bond::ENDDOWNRIGHT)) {
      res += "\\";
    } else if (val == static_cast<int>(Bond::ENDUPRIGHT)) {
      res += "/";
    } else {
      throw "Can't write smarts for this bond dir type";
    }
  } else if (descrip == "BondOrder") {
    bool reverseDative =
        (atomToLeftIdx >= 0 &&
         bond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx));
    bool doIsoSmiles =
        !bond->hasOwningMol() ||
        bond->getOwningMol().hasProp(common_properties::_doIsoSmiles);
    res += getBasicBondRepr(static_cast<Bond::BondType>(bquery->getVal()),
                            bond->getBondDir(), doIsoSmiles, reverseDative);
  } else {
    std::stringstream msg;
    msg << "Can't write smarts for this query bond type: " << descrip;
    throw msg.str().c_str();
  }
  return res;
}

std::string _recurseGetSmarts(const QueryAtom *qatom,
                              const QueryAtom::QUERYATOM_QUERY *node,
                              bool negate, unsigned int &features) {
  PRECONDITION(node, "bad node");
  // the algorithm goes something like this
  // - recursively get the smarts for the child queries
  // - combine the child smarts using the following rules:
  //      - if we are currently at an OR query, combine the subqueries with a
  //      ",",
  //        but only if neither of child smarts do not contain "," and ";"
  //        This situation leads to a no smartable situation and throw an
  //        error
  //      - if we are currently at an and query, combine the child smarts with
  //      "&"
  //        if neither of the child smarts contain a "," - otherwise combine
  //        them
  //        the child smarts with a ";"
  //
  // There is an additional complication with composite nodes that carry a
  // negation - in this
  // case we will propagate the negation to the child nodes using the
  // following rules
  //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
  //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))

  auto descrip = node->getDescription();

  unsigned int child1Features = 0;
  unsigned int child2Features = 0;
  auto chi = node->beginChildren();
  auto child1 = chi->get();
  auto dsc1 = child1->getDescription();

  ++chi;
  CHECK_INVARIANT(chi != node->endChildren(),
                  "Not enough children on the query");

  bool needParen;
  std::string csmarts1;
  // deal with the first child
  if (dsc1 == "RecursiveStructure") {
    csmarts1 = getRecursiveStructureQuerySmarts(child1);
    features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_RECURSION);
  } else if ((dsc1 != "AtomOr") && (dsc1 != "AtomAnd")) {
    // child 1 is a simple node
    const auto *tchild = static_cast<const ATOM_EQUALS_QUERY *>(child1);
    csmarts1 = getAtomSmartsSimple(qatom, tchild, needParen);
    bool nneg = (negate) ^ (tchild->getNegation());
    if (nneg) {
      csmarts1 = "!" + csmarts1;
    }
  } else {
    // child 1 is composite node - recurse
    bool nneg = (negate) ^ (child1->getNegation());
    csmarts1 = _recurseGetSmarts(qatom, child1, nneg, child1Features);
  }
  // ok if we have a negation and we have an OR , we have to change to
  // an AND since we propagated the negation
  // i.e NOT (A OR B) = (NOT (A)) AND (NOT(B))
  if (negate) {
    if (descrip == "AtomOr") {
      descrip = "AtomAnd";
    } else if (descrip == "AtomAnd") {
      descrip = "AtomOr";
    }
  }
  auto res = csmarts1;
  while (chi != node->endChildren()) {
    auto child2 = chi->get();
    ++chi;

    auto dsc2 = child2->getDescription();
    std::string csmarts2;

    // deal with the next child
    if (dsc2 == "RecursiveStructure") {
      csmarts2 = getRecursiveStructureQuerySmarts(child2);
      features |= static_cast<unsigned int>(QueryBoolFeatures::HAS_RECURSION);
    } else if ((dsc2 != "AtomOr") && (dsc2 != "AtomAnd")) {
      // child 2 is a simple node
      const auto *tchild = static_cast<const ATOM_EQUALS_QUERY *>(child2);
      csmarts2 = getAtomSmartsSimple(qatom, tchild, needParen);
      bool nneg = (negate) ^ (tchild->getNegation());
      if (nneg) {
        csmarts2 = "!" + csmarts2;
      }
    } else {
      bool nneg = (negate) ^ (child2->getNegation());
      csmarts2 = _recurseGetSmarts(qatom, child2, nneg, child2Features);
    }

    res = _combineChildSmarts(res, child1Features, csmarts2, child2Features,
                              descrip, features);
  }
  return res;
}

std::string _recurseBondSmarts(const Bond *bond,
                               const QueryBond::QUERYBOND_QUERY *node,
                               bool negate, int atomToLeftIdx,
                               unsigned int &features) {
  // the algorithm goes something like this
  // - recursively get the smarts for the child query bonds
  // - combine the child smarts using the following rules:
  //      - if we are currently at an OR query, combine the subqueries with a
  //      ",",
  //        but only if neither of child smarts do not contain "," and ";"
  //        This situation leads to a no smartable situation and throw an
  //        error
  //      - if we are currently at an and query, combine the child smarts with
  //      "&"
  //        if neither of the child smarts contain a "," - otherwise combine
  //        them
  //        the child smarts with a ";"
  //
  // There is an additional complication with composite nodes that carry a
  // negation - in this
  // case we will propagate the negation to the child nodes using the
  // following rules
  //   NOT (a AND b) = ( NOT (a)) AND ( NOT (b))
  //   NOT (a OR b) = ( NOT (a)) OR ( NOT (b))
  PRECONDITION(bond, "bad bond");
  PRECONDITION(node, "bad node");
  std::string descrip = node->getDescription();
  std::string res = "";

  const QueryBond::QUERYBOND_QUERY *child1;
  const QueryBond::QUERYBOND_QUERY *child2;
  unsigned int child1Features = 0;
  unsigned int child2Features = 0;
  QueryBond::QUERYBOND_QUERY::CHILD_VECT_CI chi;

  chi = node->beginChildren();
  child1 = chi->get();
  chi++;
  child2 = chi->get();
  chi++;
  // OK we should be at the end of vector by now - since we can have only two
  // children,
  // well - at least in this case
  CHECK_INVARIANT(chi == node->endChildren(), "Too many children on the query");

  std::string dsc1, dsc2;
  dsc1 = child1->getDescription();
  dsc2 = child2->getDescription();
  std::string csmarts1, csmarts2;

  if ((dsc1 != "BondOr") && (dsc1 != "BondAnd")) {
    // child1 is  simple node get the smarts directly
    const auto *tchild = static_cast<const BOND_EQUALS_QUERY *>(child1);
    csmarts1 = getBondSmartsSimple(bond, tchild, atomToLeftIdx);
    bool nneg = (negate) ^ (tchild->getNegation());
    if (nneg) {
      csmarts1 = "!" + csmarts1;
    }
  } else {
    // child1 is a composite node recurse further
    bool nneg = (negate) ^ (child1->getNegation());
    csmarts1 =
        _recurseBondSmarts(bond, child1, nneg, atomToLeftIdx, child1Features);
  }

  // now deal with the second child
  if ((dsc2 != "BondOr") && (dsc2 != "BondAnd")) {
    // child 2 is a simple node
    const auto *tchild = static_cast<const BOND_EQUALS_QUERY *>(child2);
    csmarts2 = getBondSmartsSimple(bond, tchild, atomToLeftIdx);
    bool nneg = (negate) ^ (tchild->getNegation());
    if (nneg) {
      csmarts2 = "!" + csmarts2;
    }
  } else {
    // child two is a composite node - recurse
    bool nneg = (negate) ^ (child2->getNegation());
    csmarts1 =
        _recurseBondSmarts(bond, child2, nneg, atomToLeftIdx, child2Features);
  }

  // ok if we have a negation and we have to change the underlying logic,
  // since we propagated the negation i.e NOT (A OR B) = (NOT (A)) AND
  // (NOT(B))
  if (negate) {
    if (descrip == "BondOr") {
      descrip = "BondAnd";
    } else if (descrip == "BondAnd") {
      descrip = "BondOr";
    }
  }
  res += _combineChildSmarts(csmarts1, child1Features, csmarts2, child2Features,
                             descrip, features);
  return res;
}

std::string FragmentSmartsConstruct(
    ROMol &mol, unsigned int atomIdx, std::vector<Canon::AtomColors> &colors,
    UINT_VECT &ranks, bool doIsomericSmiles,
    const boost::dynamic_bitset<> *bondsInPlay = nullptr) {
  Canon::MolStack molStack;
  molStack.reserve(mol.getNumAtoms() + mol.getNumBonds());
  std::stringstream res;

  // this is dirty trick get around the fact that canonicalizeFragment
  // thinks we already called findSSSR - to do some atom ranking
  // but for smarts we are going to ignore that part. We will artificially
  // set the "SSSR" property to an empty property
  mol.getRingInfo()->reset();
  mol.getRingInfo()->initialize();
  for (auto &atom : mol.atoms()) {
    atom->updatePropertyCache(false);
  }

  // Another dirty trick to avoid reordering of lead chiral atoms in
  // canonicalizeFragment: since we are writing SMARTS, there won't be a
  // reordering on parsing.
  // The trick is as easy as choosing the first non-chiral atom we can find as
  // root of the string. This should not be a problem, since SMARTS do not get
  // canonicalized.
  if (molStack.empty()) {
    for (const auto atom : mol.atoms()) {
      if (colors[atom->getIdx()] == Canon::WHITE_NODE &&
          atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CCW &&
          atom->getChiralTag() != Atom::CHI_TETRAHEDRAL_CW) {
        atomIdx = atom->getIdx();
        break;
      }
    }
  }

  Canon::canonicalizeFragment(mol, atomIdx, colors, ranks, molStack,
                              bondsInPlay, nullptr, doIsomericSmiles);

  // now clear the "SSSR" property
  mol.getRingInfo()->reset();
  for (auto &msCI : molStack) {
    switch (msCI.type) {
      case Canon::MOL_STACK_ATOM: {
        auto *qatm = static_cast<QueryAtom *>(msCI.obj.atom);
        res << SmartsWrite::GetAtomSmarts(qatm);
        break;
      }
      case Canon::MOL_STACK_BOND: {
        auto *qbnd = static_cast<QueryBond *>(msCI.obj.bond);
        res << SmartsWrite::GetBondSmarts(qbnd, msCI.number);
        break;
      }
      case Canon::MOL_STACK_RING: {
        if (msCI.number < 10) {
          res << msCI.number;
        } else {
          res << "%" << msCI.number;
        }
        break;
      }
      case Canon::MOL_STACK_BRANCH_OPEN: {
        res << "(";
        break;
      }
      case Canon::MOL_STACK_BRANCH_CLOSE: {
        res << ")";
        break;
      }
      default:
        break;
    }
  }
  return res.str();
}

// this is the used when converting a SMILES or
// non-query atom from a mol file into SMARTS.
std::string getNonQueryAtomSmarts(const QueryAtom *qatom) {
  PRECONDITION(qatom, "bad atom");
  PRECONDITION(!qatom->hasQuery(), "atom should not have query");
  std::stringstream res;
  res << "[";

  int isotope = qatom->getIsotope();
  if (isotope) {
    res << isotope;
  }

  if (SmilesWrite::inOrganicSubset(qatom->getAtomicNum())) {
    res << "#" << qatom->getAtomicNum();
  } else {
    res << qatom->getSymbol();
  }

  if (qatom->hasOwningMol() &&
      qatom->getOwningMol().hasProp(common_properties::_doIsoSmiles)) {
    if (qatom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
        !qatom->hasProp(_qatomHasStereoSet) &&
        !qatom->hasProp(common_properties::_brokenChirality)) {
      qatom->setProp(_qatomHasStereoSet, 1);
      switch (qatom->getChiralTag()) {
        case Atom::CHI_TETRAHEDRAL_CW:
          res << "@@";
          break;
        case Atom::CHI_TETRAHEDRAL_CCW:
          res << "@";
          break;
        default:
          break;
      }
    }
  }

  auto hs = qatom->getNumExplicitHs();
  // FIX: probably should be smarter about Hs:
  if (hs) {
    res << "H";
    if (hs > 1) {
      res << hs;
    }
  }
  auto chg = qatom->getFormalCharge();
  if (chg) {
    if (chg == -1) {
      res << "-";
    } else if (chg == 1) {
      res << "+";
    } else if (chg < 0) {
      res << qatom->getFormalCharge();
    } else {
      res << "+" << qatom->getFormalCharge();
    }
  }
  int mapNum;
  if (qatom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    res << ":";
    res << mapNum;
  }
  res << "]";
  return res.str();
}

// this is the used when converting a SMILES or
// non-query bond from a mol file into SMARTS.
std::string getNonQueryBondSmarts(const QueryBond *qbond, int atomToLeftIdx) {
  PRECONDITION(qbond, "bad bond");
  PRECONDITION(!qbond->hasQuery(), "bond should not have query");
  RDUNUSED_PARAM(atomToLeftIdx);
  std::string res;

  if (qbond->getIsAromatic()) {
    res = ":";
  } else {
    bool reverseDative =
        (atomToLeftIdx >= 0 &&
         qbond->getBeginAtomIdx() == static_cast<unsigned int>(atomToLeftIdx));
    bool doIsoSmiles =
        !qbond->hasOwningMol() ||
        qbond->getOwningMol().hasProp(common_properties::_doIsoSmiles);
    res = getBasicBondRepr(qbond->getBondType(), qbond->getBondDir(),
                           doIsoSmiles, reverseDative);
  }

  return res;
}

std::string molToSmarts(const ROMol &inmol, bool doIsomericSmiles,
                        std::vector<AtomColors> &&colors,
                        const boost::dynamic_bitset<> *bondsInPlay = nullptr) {
  ROMol mol(inmol);
  const unsigned int nAtoms = mol.getNumAtoms();
  UINT_VECT ranks;
  ranks.reserve(nAtoms);
  // For smiles writing we would be canonicalizing but we will not do that here.
  // We will simply use the atom indices as the rank
  for (const auto &atom : mol.atoms()) {
    ranks.push_back(atom->getIdx());
  }

  if (doIsomericSmiles) {
    mol.setProp(common_properties::_doIsoSmiles, 1);
  }

  std::string res;
  auto colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
  while (colorIt != colors.end()) {
    unsigned int nextAtomIdx = 0;
    unsigned int nextRank;
    std::string subSmi;
    nextRank = nAtoms + 1;
    for (unsigned int i = 0; i < nAtoms; ++i) {
      if (colors[i] == Canon::WHITE_NODE && ranks[i] < nextRank) {
        nextRank = ranks[i];
        nextAtomIdx = i;
      }
    }

    subSmi = FragmentSmartsConstruct(mol, nextAtomIdx, colors, ranks,
                                     doIsomericSmiles, bondsInPlay);
    res += subSmi;

    colorIt = std::find(colors.begin(), colors.end(), Canon::WHITE_NODE);
    if (colorIt != colors.end()) {
      res += ".";
    }
  }
  return res;
}

}  // namespace

namespace SmartsWrite {
std::string GetAtomSmarts(const QueryAtom *qatom) {
  PRECONDITION(qatom, "bad atom");
  std::string res;
  bool needParen = false;

  // BOOST_LOG(rdInfoLog)<<"Atom: " <<qatom->getIdx()<<std::endl;
  if (!qatom->hasQuery()) {
    res = getNonQueryAtomSmarts(qatom);
    // BOOST_LOG(rdInfoLog)<<"\tno query:" <<res;
    return res;
  }
  const auto query = qatom->getQuery();
  PRECONDITION(query, "atom has no query");
  // describeQuery(query);
  unsigned int queryFeatures = 0;
  std::string descrip = qatom->getQuery()->getDescription();
  if (descrip.empty()) {
    // we have simple atom - just generate the smiles and return
    res = SmilesWrite::GetAtomSmiles(qatom);
    if (res[0] == '[') {
      // chop the brackets off, we'll put them back on later:
      needParen = true;
      res = res.substr(1, res.size() - 2);
    }
  } else if ((descrip == "AtomOr") || (descrip == "AtomAnd")) {
    // we have a composite query
    needParen = true;
    res = _recurseGetSmarts(qatom, query, query->getNegation(), queryFeatures);
    if (res.length() == 1) {  // single atom symbol we don't need parens
      needParen = false;
    }
  } else if (descrip == "RecursiveStructure") {
    // it's a bare recursive structure query:
    res = getRecursiveStructureQuerySmarts(query);
    needParen = true;
  } else {  // we have a simple smarts
    auto *tquery = static_cast<ATOM_EQUALS_QUERY *>(qatom->getQuery());
    res = getAtomSmartsSimple(qatom, tquery, needParen);
    if (tquery->getNegation()) {
      res = "!" + res;
    }
  }
  std::string mapNum;
  if (qatom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
    needParen = true;
    res += ":" + mapNum;
  }
  if (needParen) {
    res = "[" + res + "]";
  }
  return res;
}

std::string GetBondSmarts(const QueryBond *bond, int atomToLeftIdx) {
  PRECONDITION(bond, "bad bond");
  std::string res = "";

  // BOOST_LOG(rdInfoLog) << "bond: " << bond->getIdx() << std::endl;
  ;
  // it is possible that we are regular single bond and we don't need to write
  // anything
  if (!bond->hasQuery()) {
    res = getNonQueryBondSmarts(bond, atomToLeftIdx);
    // BOOST_LOG(rdInfoLog) << "\tno query:" << res << std::endl;
    return res;
  }
  // describeQuery(bond->getQuery());
  if ((typeid(*bond) == typeid(Bond)) &&
      ((bond->getBondType() == Bond::SINGLE) ||
       (bond->getBondType() == Bond::AROMATIC))) {
    BOOST_LOG(rdInfoLog) << "\tbasic:" << res << std::endl;
    return res;
  }
  const auto query = bond->getQuery();
  PRECONDITION(query, "bond has no query");

  unsigned int queryFeatures = 0;
  auto descrip = query->getDescription();
  if ((descrip == "BondAnd") || (descrip == "BondOr")) {
    // composite query
    res = _recurseBondSmarts(bond, query, query->getNegation(), atomToLeftIdx,
                             queryFeatures);
  } else {
    // simple query
    if (query->getNegation()) {
      res = "!";
    }
    const auto *tquery = static_cast<const BOND_EQUALS_QUERY *>(query);
    res += getBondSmartsSimple(bond, tquery, atomToLeftIdx);
  }
  // BOOST_LOG(rdInfoLog) << "\t  query:" << descrip << " " << res << std::endl;
  return res;
}

}  // end of namespace SmartsWrite

std::string MolToSmarts(const ROMol &mol, bool doIsomericSmarts) {
  const unsigned int nAtoms = mol.getNumAtoms();
  if (!nAtoms) {
    return "";
  }

  std::vector<AtomColors> colors(nAtoms, Canon::WHITE_NODE);
  return molToSmarts(mol, doIsomericSmarts, std::move(colors));
}

std::string MolFragmentToSmarts(const ROMol &mol,
                                const std::vector<int> &atomsToUse,
                                const std::vector<int> *bondsToUse,
                                bool doIsomericSmarts) {
  PRECONDITION(!atomsToUse.empty(), "no atoms provided");
  PRECONDITION(!bondsToUse || !bondsToUse->empty(), "no bonds provided");

  auto nAtoms = mol.getNumAtoms();
  if (!nAtoms) {
    return "";
  }

  std::unique_ptr<boost::dynamic_bitset<>> bondsInPlay(nullptr);
  if (bondsToUse != nullptr) {
    bondsInPlay.reset(new boost::dynamic_bitset<>(mol.getNumBonds(), 0));
    for (auto bidx : *bondsToUse) {
      bondsInPlay->set(bidx);
    }
  }

  // Mark all atoms except the ones in atomIndices as already processed.
  // white: unprocessed
  // grey: partial
  // black: complete
  std::vector<AtomColors> colors(nAtoms, Canon::BLACK_NODE);
  for (const auto &idx : atomsToUse) {
    colors[idx] = Canon::WHITE_NODE;
  }

  return molToSmarts(mol, doIsomericSmarts, std::move(colors),
                     bondsInPlay.get());
}

}  // namespace RDKit