//
//  Copyright (c) 2007-2021, Novartis Institutes for BioMedical Research Inc.
//  and other RDKit contributors
//
//  All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Novartis Institutes for BioMedical Research Inc.
//       nor the names of its contributors may be used to endorse or promote
//       products derived from this software without specific prior written
//       permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionPickler.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/QueryOps.h>
#include <boost/dynamic_bitset.hpp>
#include <map>
#include <algorithm>
#include <GraphMol/ChemTransforms/ChemTransforms.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include "GraphMol/ChemReactions/ReactionRunner.h"

namespace RDKit {

std::vector<MOL_SPTR_VECT> ChemicalReaction::runReactants(
    const MOL_SPTR_VECT reactants, unsigned int maxProducts) const {
  return run_Reactants(*this, reactants, maxProducts);
}

std::vector<MOL_SPTR_VECT> ChemicalReaction::runReactant(
    const ROMOL_SPTR reactant, unsigned int reactionTemplateIdx) const {
  return run_Reactant(*this, reactant, reactionTemplateIdx);
}

bool ChemicalReaction::runReactant(RWMol &reactant,
                                   bool removeUnmatchedAtoms) const {
  if (getReactants().size() != 1 || getProducts().size() != 1) {
    throw ChemicalReactionException(
        "Only single reactant - single product reactions can be run in place.");
  }
  return run_Reactant(*this, reactant, removeUnmatchedAtoms);
}

ChemicalReaction::ChemicalReaction(const std::string &pickle) {
  ReactionPickler::reactionFromPickle(pickle, this);
}

void ChemicalReaction::initReactantMatchers(bool silent) {
  unsigned int nWarnings, nErrors;
  if (!this->validate(nWarnings, nErrors, silent)) {
    BOOST_LOG(rdErrorLog) << "initialization failed\n";
    this->df_needsInit = true;
  } else {
    this->df_needsInit = false;
  }
}

bool ChemicalReaction::validate(unsigned int &numWarnings,
                                unsigned int &numErrors, bool silent) const {
  bool res = true;
  numWarnings = 0;
  numErrors = 0;

  if (!this->getNumReactantTemplates()) {
    if (!silent) {
      BOOST_LOG(rdErrorLog) << "reaction has no reactants\n";
    }
    numErrors++;
    res = false;
  }

  if (!this->getNumProductTemplates()) {
    if (!silent) {
      BOOST_LOG(rdErrorLog) << "reaction has no products\n";
    }
    numErrors++;
    res = false;
  }

  std::vector<int> mapNumbersSeen;
  std::map<int, const Atom *> reactingAtoms;
  unsigned int molIdx = 0;
  for (auto molIter = this->beginReactantTemplates();
       molIter != this->endReactantTemplates(); ++molIter) {
    bool thisMolMapped = false;
    for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
         atomIt != (*molIter)->endAtoms(); ++atomIt) {
      int mapNum;
      if ((*atomIt)->getPropIfPresent(common_properties::molAtomMapNumber,
                                      mapNum)) {
        thisMolMapped = true;
        if (std::find(mapNumbersSeen.begin(), mapNumbersSeen.end(), mapNum) !=
            mapNumbersSeen.end()) {
          if (!silent) {
            BOOST_LOG(rdErrorLog) << "reactant atom-mapping number " << mapNum
                                  << " found multiple times.\n";
          }
          numErrors++;
          res = false;
        } else {
          mapNumbersSeen.push_back(mapNum);
          reactingAtoms[mapNum] = *atomIt;
        }
      }
    }
    if (!thisMolMapped) {
      if (!silent) {
        BOOST_LOG(rdWarningLog)
            << "reactant " << molIdx << " has no mapped atoms.\n";
      }
      numWarnings++;
    }
    molIdx++;
  }

  std::vector<int> productNumbersSeen;
  molIdx = 0;
  for (auto molIter = this->beginProductTemplates();
       molIter != this->endProductTemplates(); ++molIter) {
    // clear out some possible cached properties to prevent
    // misleading warnings
    for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
         atomIt != (*molIter)->endAtoms(); ++atomIt) {
      if ((*atomIt)->hasProp(common_properties::_QueryFormalCharge)) {
        (*atomIt)->clearProp(common_properties::_QueryFormalCharge);
      }
      if ((*atomIt)->hasProp(common_properties::_QueryHCount)) {
        (*atomIt)->clearProp(common_properties::_QueryHCount);
      }
      if ((*atomIt)->hasProp(common_properties::_QueryMass)) {
        (*atomIt)->clearProp(common_properties::_QueryMass);
      }
      if ((*atomIt)->hasProp(common_properties::_QueryIsotope)) {
        (*atomIt)->clearProp(common_properties::_QueryIsotope);
      }
    }
    bool thisMolMapped = false;
    for (ROMol::AtomIterator atomIt = (*molIter)->beginAtoms();
         atomIt != (*molIter)->endAtoms(); ++atomIt) {
      int mapNum;
      if ((*atomIt)->getPropIfPresent(common_properties::molAtomMapNumber,
                                      mapNum)) {
        thisMolMapped = true;
        bool seenAlready =
            std::find(productNumbersSeen.begin(), productNumbersSeen.end(),
                      mapNum) != productNumbersSeen.end();
        if (seenAlready) {
          if (!silent) {
            BOOST_LOG(rdWarningLog) << "product atom-mapping number " << mapNum
                                    << " found multiple times.\n";
          }
          numWarnings++;
          // ------------
          //   Always check to see if the atoms connectivity changes independent
          //   if it is mapped multiple times
          // ------------
          const Atom *rAtom = reactingAtoms[mapNum];
          CHECK_INVARIANT(rAtom, "missing atom");
          if (rAtom->getDegree() != (*atomIt)->getDegree()) {
            (*atomIt)->setProp(common_properties::_ReactionDegreeChanged, 1);
          }

        } else {
          productNumbersSeen.push_back(mapNum);
        }
        auto ivIt =
            std::find(mapNumbersSeen.begin(), mapNumbersSeen.end(), mapNum);
        if (ivIt == mapNumbersSeen.end()) {
          if (!seenAlready) {
            if (!silent) {
              BOOST_LOG(rdWarningLog) << "product atom-mapping number "
                                      << mapNum << " not found in reactants.\n";
            }
            numWarnings++;
          }
        } else {
          mapNumbersSeen.erase(ivIt);

          // ------------
          //   The atom is mapped, check to see if its connectivity changes
          // ------------
          const Atom *rAtom = reactingAtoms[mapNum];
          CHECK_INVARIANT(rAtom, "missing atom");
          if (rAtom->getDegree() != (*atomIt)->getDegree()) {
            (*atomIt)->setProp(common_properties::_ReactionDegreeChanged, 1);
          }
        }
      }

      // ------------
      //    Deal with queries
      // ------------
      if ((*atomIt)->hasQuery()) {
        std::list<const Atom::QUERYATOM_QUERY *> queries;
        queries.push_back((*atomIt)->getQuery());
        while (!queries.empty()) {
          const Atom::QUERYATOM_QUERY *query = queries.front();
          queries.pop_front();
          for (auto qIter = query->beginChildren();
               qIter != query->endChildren(); ++qIter) {
            queries.push_back((*qIter).get());
          }
          if (query->getDescription() == "AtomFormalCharge" ||
              query->getDescription() == "AtomNegativeFormalCharge") {
            int qval;
            int neg =
                query->getDescription() == "AtomNegativeFormalCharge" ? -1 : 1;
            if ((*atomIt)->getPropIfPresent(
                    common_properties::_QueryFormalCharge, qval) &&
                (neg * qval) !=
                    static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
              if (!silent) {
                BOOST_LOG(rdWarningLog)
                    << "atom " << (*atomIt)->getIdx() << " in product "
                    << molIdx << " has multiple charge specifications.\n";
              }
              numWarnings++;
            } else {
              int neg = query->getDescription() == "AtomNegativeFormalCharge"
                            ? -1
                            : 1;
              (*atomIt)->setProp(
                  common_properties::_QueryFormalCharge,
                  neg *
                      static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
            }
          } else if (query->getDescription() == "AtomHCount") {
            int qval;
            if ((*atomIt)->getPropIfPresent(common_properties::_QueryHCount,
                                            qval) &&
                qval !=
                    static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
              if (!silent) {
                BOOST_LOG(rdWarningLog)
                    << "atom " << (*atomIt)->getIdx() << " in product "
                    << molIdx << " has multiple H count specifications.\n";
              }
              numWarnings++;
            } else {
              (*atomIt)->setProp(
                  common_properties::_QueryHCount,
                  static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
            }
          } else if (query->getDescription() == "AtomMass") {
            int qval;
            if ((*atomIt)->getPropIfPresent(common_properties::_QueryMass,
                                            qval) &&
                qval !=
                    static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
              if (!silent) {
                BOOST_LOG(rdWarningLog)
                    << "atom " << (*atomIt)->getIdx() << " in product "
                    << molIdx << " has multiple mass specifications.\n";
              }
              numWarnings++;
            } else {
              (*atomIt)->setProp(
                  common_properties::_QueryMass,
                  static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal() /
                      massIntegerConversionFactor);
            }
          } else if (query->getDescription() == "AtomIsotope") {
            int qval;
            if ((*atomIt)->getPropIfPresent(common_properties::_QueryIsotope,
                                            qval) &&
                qval !=
                    static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal()) {
              if (!silent) {
                BOOST_LOG(rdWarningLog)
                    << "atom " << (*atomIt)->getIdx() << " in product "
                    << molIdx << " has multiple isotope specifications.\n";
              }
              numWarnings++;
            } else {
              (*atomIt)->setProp(
                  common_properties::_QueryIsotope,
                  static_cast<const ATOM_EQUALS_QUERY *>(query)->getVal());
            }
          }
        }
      }
    }
    if (!thisMolMapped) {
      if (!silent) {
        BOOST_LOG(rdWarningLog)
            << "product " << molIdx << " has no mapped atoms.\n";
      }
      numWarnings++;
    }
    molIdx++;
  }
  if (!mapNumbersSeen.empty()) {
    if (!silent) {
      std::ostringstream ostr;
      ostr
          << "mapped atoms in the reactants were not mapped in the products.\n";
      ostr << "  unmapped numbers are: ";
      for (std::vector<int>::const_iterator ivIt = mapNumbersSeen.begin();
           ivIt != mapNumbersSeen.end(); ++ivIt) {
        ostr << *ivIt << " ";
      }
      ostr << "\n";
      BOOST_LOG(rdWarningLog) << ostr.str();
    }
    numWarnings++;
  }

  return res;
}

bool isMoleculeReactantOfReaction(const ChemicalReaction &rxn, const ROMol &mol,
                                  std::vector<unsigned int> &which,
                                  bool stopAtFirstMatch) {
  if (!rxn.isInitialized()) {
    throw ChemicalReactionException(
        "initReactantMatchers() must be called first");
  }
  which.clear();
  unsigned int reactant_template_idx = 0;
  for (auto iter = rxn.beginReactantTemplates();
       iter != rxn.endReactantTemplates(); ++iter, ++reactant_template_idx) {
    auto tvect = SubstructMatch(mol, **iter, rxn.getSubstructParams());
    if (!tvect.empty()) {
      which.push_back(reactant_template_idx);
      if (stopAtFirstMatch) {
        return true;
      }
    }
  }
  return !which.empty();
}

bool isMoleculeReactantOfReaction(const ChemicalReaction &rxn, const ROMol &mol,
                                  unsigned int &which) {
  std::vector<unsigned int> matches;
  bool is_reactant = isMoleculeReactantOfReaction(rxn, mol, matches, true);
  if (matches.empty()) {
    which = rxn.getNumReactantTemplates();
  } else {
    which = matches[0];
  }
  return is_reactant;
}

bool isMoleculeReactantOfReaction(const ChemicalReaction &rxn,
                                  const ROMol &mol) {
  unsigned int ignore;
  return isMoleculeReactantOfReaction(rxn, mol, ignore);
}

bool isMoleculeProductOfReaction(const ChemicalReaction &rxn, const ROMol &mol,
                                 std::vector<unsigned int> &which,
                                 bool stopAtFirstMatch) {
  if (!rxn.isInitialized()) {
    throw ChemicalReactionException(
        "initReactantMatchers() must be called first");
  }
  which.clear();
  unsigned int product_template_idx = 0;
  for (auto iter = rxn.beginProductTemplates();
       iter != rxn.endProductTemplates(); ++iter, ++product_template_idx) {
    auto tvect = SubstructMatch(mol, **iter, rxn.getSubstructParams());
    if (!tvect.empty()) {
      which.push_back(product_template_idx);
      if (stopAtFirstMatch) {
        return true;
      }
    }
  }
  return !which.empty();
}

bool isMoleculeProductOfReaction(const ChemicalReaction &rxn, const ROMol &mol,
                                 unsigned int &which) {
  std::vector<unsigned int> matches;
  bool is_product = isMoleculeProductOfReaction(rxn, mol, matches, true);
  if (matches.empty()) {
    which = rxn.getNumProductTemplates();
  } else {
    which = matches[0];
  }
  return is_product;
}

bool isMoleculeProductOfReaction(const ChemicalReaction &rxn,
                                 const ROMol &mol) {
  unsigned int ignore;
  return isMoleculeProductOfReaction(rxn, mol, ignore);
}

bool isMoleculeAgentOfReaction(const ChemicalReaction &rxn, const ROMol &mol,
                               unsigned int &which) {
  if (!rxn.isInitialized()) {
    throw ChemicalReactionException(
        "initReactantMatchers() must be called first");
  }
  which = 0;
  for (auto templ : rxn.getAgents()) {
    if (templ->getNumHeavyAtoms() != mol.getNumHeavyAtoms()) {
      ++which;
      continue;
    }
    if (templ->getNumBonds() != mol.getNumBonds()) {
      ++which;
      continue;
    }
    if (MolOps::getAvgMolWt(*templ) != MolOps::getAvgMolWt(mol)) {
      ++which;
      continue;
    }
    auto tvect = SubstructMatch(mol, *templ, rxn.getSubstructParams());
    if (!tvect.empty()) {
      return true;
    }
    ++which;
  }
  return false;
}

bool isMoleculeAgentOfReaction(const ChemicalReaction &rxn, const ROMol &mol) {
  unsigned int ignore;
  return isMoleculeAgentOfReaction(rxn, mol, ignore);
}

void addRecursiveQueriesToReaction(
    ChemicalReaction &rxn, const std::map<std::string, ROMOL_SPTR> &queries,
    const std::string &propName,
    std::vector<std::vector<std::pair<unsigned int, std::string>>>
        *reactantLabels) {
  if (!rxn.isInitialized()) {
    throw ChemicalReactionException(
        "initReactantMatchers() must be called first");
  }

  if (reactantLabels != nullptr) {
    (*reactantLabels).resize(0);
  }

  for (MOL_SPTR_VECT::const_iterator rIt = rxn.beginReactantTemplates();
       rIt != rxn.endReactantTemplates(); ++rIt) {
    if (reactantLabels != nullptr) {
      std::vector<std::pair<unsigned int, std::string>> labels;
      addRecursiveQueries(**rIt, queries, propName, &labels);
      (*reactantLabels).push_back(labels);
    } else {
      addRecursiveQueries(**rIt, queries, propName);
    }
  }
}

namespace {
// recursively looks for atomic number queries anywhere in this set of children
// or its children
int numComplexQueries(
    Queries::Query<int, Atom const *, true>::CHILD_VECT_CI childIt,
    Queries::Query<int, Atom const *, true>::CHILD_VECT_CI endChildren) {
  int res = 0;
  while (childIt != endChildren) {
    std::string descr = (*childIt)->getDescription();
    if (descr == "AtomAtomicNum" || descr == "AtomNull") {
      ++res;
    } else {
      res += numComplexQueries((*childIt)->beginChildren(),
                               (*childIt)->endChildren());
    }
    ++childIt;
  }
  return res;
}
#if 0
// FIX: this is adapted from Fingerprints.cpp and we really should have code
// like this centralized
bool isComplexQuery(const Atom &a) {
  if (!a.hasQuery()) return false;
  // negated things are always complex:
  if (a.getQuery()->getNegation()) return true;
  std::string descr = a.getQuery()->getDescription();
  if (descr == "AtomAtomicNum") return false;
  if (descr == "AtomOr" || descr == "AtomXor") return true;
  if (descr == "AtomAnd") {
    auto childIt = a.getQuery()->beginChildren();
    int ncq = numComplexQueries(childIt, a.getQuery()->endChildren());
    if (ncq == 1) {
      return false;
    }
  }
  return true;
}
#endif
bool isChangedAtom(const Atom &rAtom, const Atom &pAtom, int mapNum,
                   const std::map<int, const Atom *> &mappedProductAtoms) {
  PRECONDITION(mappedProductAtoms.find(mapNum) != mappedProductAtoms.end(),
               "atom not mapped in products");

  if (rAtom.getAtomicNum() != pAtom.getAtomicNum() &&
      pAtom.getAtomicNum() > 0) {
    // the atomic number changed and the product wasn't a dummy
    return true;
  } else if (rAtom.getDegree() != pAtom.getDegree()) {
    // the degree changed
    return true;
  } else if (pAtom.getAtomicNum() > 0 && isComplexQuery(&rAtom)) {
    // more than a simple query
    return true;
  }

  // now check bond layout:
  std::map<unsigned int, const Bond *> reactantBonds;
  for (const auto &nbrIdx : boost::make_iterator_range(
           rAtom.getOwningMol().getAtomNeighbors(&rAtom))) {
    const Atom *nbr = rAtom.getOwningMol()[nbrIdx];
    int mapNum;
    if (nbr->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
      reactantBonds[mapNum] = rAtom.getOwningMol().getBondBetweenAtoms(
          rAtom.getIdx(), nbr->getIdx());
    } else {
      // if we have an un-mapped neighbor, we are automatically a reacting
      // atom:
      return true;
    }
  }
  for (const auto &nbrIdx : boost::make_iterator_range(
           pAtom.getOwningMol().getAtomNeighbors(&pAtom))) {
    const Atom *nbr = pAtom.getOwningMol()[nbrIdx];
    int mapNum;
    if (nbr->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
      // if we don't have a bond to a similarly mapped atom in the reactant,
      // we're done:
      if (reactantBonds.find(mapNum) == reactantBonds.end()) {
        return true;
      }
      const Bond *rBond = reactantBonds[mapNum];
      const Bond *pBond = pAtom.getOwningMol().getBondBetweenAtoms(
          pAtom.getIdx(), nbr->getIdx());

      // bond comparison logic:
      if (rBond->hasQuery()) {
        if (!pBond->hasQuery()) {
          // reactant query, product not query: always a change
          return true;
        } else {
          if (pBond->getQuery()->getDescription() == "BondNull") {
            // null queries are trump, they match everything
          } else if (rBond->getBondType() == Bond::SINGLE &&
                     pBond->getBondType() == Bond::SINGLE &&
                     ((rBond->getQuery()->getDescription() == "BondOr" &&
                       pBond->getQuery()->getDescription() == "BondOr") ||
                      (rBond->getQuery()->getDescription() ==
                           "SingleOrAromaticBond" &&
                       pBond->getQuery()->getDescription() ==
                           "SingleOrAromaticBond"))) {
            // The SMARTS parser tags unspecified bonds as single, but then
            // adds a query so that they match single or double. these cases
            // match
          } else {
            if (rBond->getBondType() == pBond->getBondType() &&
                rBond->getQuery()->getDescription() == "BondOrder" &&
                pBond->getQuery()->getDescription() == "BondOrder" &&
                static_cast<BOND_EQUALS_QUERY *>(rBond->getQuery())->getVal() ==
                    static_cast<BOND_EQUALS_QUERY *>(pBond->getQuery())
                        ->getVal()) {
              // bond order queries with equal orders also match
            } else {
              // anything else does not match
              return true;
            }
          }
        }
      } else if (pBond->hasQuery()) {
        // reactant not query, product query
        // if product is anything other than the null query
        // it's a change:
        if (pBond->getQuery()->getDescription() != "BondNull") {
          return true;
        }

      } else {
        // neither has a query, just compare the types
        if (rBond->getBondType() != pBond->getBondType()) {
          return true;
        }
      }
    }
  }

  // haven't found anything to say that we are changed, so we must
  // not be
  return false;
}

template <class T>
void getMappedAtoms(T &rIt, std::map<int, const Atom *> &mappedAtoms) {
  for (const auto atom : rIt->atoms()) {
    // we only worry about mapped atoms:
    int mapNum;
    if (atom->getPropIfPresent(common_properties::molAtomMapNumber, mapNum)) {
      mappedAtoms[mapNum] = atom;
    }
  }
}
}  // end of anonymous namespace

VECT_INT_VECT getReactingAtoms(const ChemicalReaction &rxn,
                               bool mappedAtomsOnly) {
  if (!rxn.isInitialized()) {
    throw ChemicalReactionException(
        "initReactantMatchers() must be called first");
  }
  VECT_INT_VECT res;
  res.resize(rxn.getNumReactantTemplates());

  // find mapped atoms in the products :
  std::map<int, const Atom *> mappedProductAtoms;
  for (auto rIt = rxn.beginProductTemplates(); rIt != rxn.endProductTemplates();
       ++rIt) {
    getMappedAtoms(*rIt, mappedProductAtoms);
  }

  // now loop over mapped atoms in the reactants, keeping track of
  // which reactant they are associated with, and check for changes.
  auto resIt = res.begin();
  for (auto rIt = rxn.beginReactantTemplates();
       rIt != rxn.endReactantTemplates(); ++rIt, ++resIt) {
    for (const auto oAtom : (*rIt)->atoms()) {
      // unmapped atoms are definitely changing:
      int mapNum;
      if (!oAtom->getPropIfPresent(common_properties::molAtomMapNumber,
                                   mapNum)) {
        if (!mappedAtomsOnly) {
          resIt->push_back(oAtom->getIdx());
        }
      } else {
        // but mapped ones require more careful consideration
        // if this is found in a reactant:
        if (mappedProductAtoms.find(mapNum) != mappedProductAtoms.end()) {
          if (isChangedAtom(*oAtom, *(mappedProductAtoms[mapNum]), mapNum,
                            mappedProductAtoms)) {
            resIt->push_back(oAtom->getIdx());
          }
        }
      }
    }
  }
  return res;
}

void ChemicalReaction::removeUnmappedReactantTemplates(
    double thresholdUnmappedAtoms, bool moveToAgentTemplates,
    MOL_SPTR_VECT *targetVector) {
  MOL_SPTR_VECT res_reactantTemplates;
  for (auto iter = beginReactantTemplates(); iter != endReactantTemplates();
       ++iter) {
    if (isReactionTemplateMoleculeAgent(*iter->get(), thresholdUnmappedAtoms)) {
      if (moveToAgentTemplates) {
        m_agentTemplates.push_back(*iter);
      }
      if (targetVector) {
        targetVector->push_back(*iter);
      }
    } else {
      res_reactantTemplates.push_back(*iter);
    }
  }
  m_reactantTemplates.clear();
  m_reactantTemplates.insert(m_reactantTemplates.begin(),
                             res_reactantTemplates.begin(),
                             res_reactantTemplates.end());
  res_reactantTemplates.clear();
}

void ChemicalReaction::removeUnmappedProductTemplates(
    double thresholdUnmappedAtoms, bool moveToAgentTemplates,
    MOL_SPTR_VECT *targetVector) {
  MOL_SPTR_VECT res_productTemplates;
  for (auto iter = beginProductTemplates(); iter != endProductTemplates();
       ++iter) {
    if (isReactionTemplateMoleculeAgent(*iter->get(), thresholdUnmappedAtoms)) {
      if (moveToAgentTemplates) {
        m_agentTemplates.push_back(*iter);
      }
      if (targetVector) {
        targetVector->push_back(*iter);
      }
    } else {
      res_productTemplates.push_back(*iter);
    }
  }
  m_productTemplates.clear();
  m_productTemplates.insert(m_productTemplates.begin(),
                            res_productTemplates.begin(),
                            res_productTemplates.end());
  res_productTemplates.clear();
}

void ChemicalReaction::removeAgentTemplates(MOL_SPTR_VECT *targetVector) {
  if (targetVector) {
    for (auto iter = beginAgentTemplates(); iter != endAgentTemplates();
         ++iter) {
      targetVector->push_back(*iter);
    }
  }
  m_agentTemplates.clear();
}

}  // namespace RDKit
