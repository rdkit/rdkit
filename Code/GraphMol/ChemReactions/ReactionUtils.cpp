//
//  Copyright (c) 2014, Novartis Institutes for BioMedical Research Inc.
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
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <cmath>

namespace RDKit {

MOL_SPTR_VECT::const_iterator getStartIterator(const ChemicalReaction &rxn,
                                               ReactionMoleculeType t) {
  MOL_SPTR_VECT::const_iterator begin;
  if (t == Reactant) {
    begin = rxn.beginReactantTemplates();
  }
  if (t == Product) {
    begin = rxn.beginProductTemplates();
    ;
  }
  if (t == Agent) {
    begin = rxn.beginAgentTemplates();
  }
  return begin;
}

MOL_SPTR_VECT::const_iterator getEndIterator(const ChemicalReaction &rxn,
                                             ReactionMoleculeType t) {
  MOL_SPTR_VECT::const_iterator end;
  if (t == Reactant) {
    end = rxn.endReactantTemplates();
  }
  if (t == Product) {
    end = rxn.endProductTemplates();
    ;
  }
  if (t == Agent) {
    end = rxn.endAgentTemplates();
  }
  return end;
}

namespace {

bool hasReactionMoleculeTemplateSubstructMatch(
    const RDKit::ChemicalReaction &rxn,
    const RDKit::ChemicalReaction &query_rxn, RDKit::ReactionMoleculeType t) {
  for (auto begin = getStartIterator(rxn, t); begin != getEndIterator(rxn, t);
       ++begin) {
    for (auto begin_query = getStartIterator(query_rxn, t);
         begin_query != getEndIterator(query_rxn, t); ++begin_query) {
      MatchVectType tvect;
      if (SubstructMatch(*begin->get(), *begin_query->get(), tvect)) {
        return true;
      }
    }
  }
  return false;
}
}  // namespace

bool hasReactantTemplateSubstructMatch(const ChemicalReaction &rxn,
                                       const ChemicalReaction &query_rxn) {
  if (rxn.getNumReactantTemplates() < query_rxn.getNumReactantTemplates()) {
    return false;
  }
  if (query_rxn.getNumReactantTemplates() == 0) {
    return true;
  }
  return hasReactionMoleculeTemplateSubstructMatch(rxn, query_rxn, Reactant);
}

bool hasProductTemplateSubstructMatch(const ChemicalReaction &rxn,
                                      const ChemicalReaction &query_rxn) {
  if (rxn.getNumProductTemplates() < query_rxn.getNumProductTemplates()) {
    return false;
  }
  if (query_rxn.getNumProductTemplates() == 0) {
    return true;
  }
  return hasReactionMoleculeTemplateSubstructMatch(rxn, query_rxn, Product);
}

bool hasAgentTemplateSubstructMatch(const ChemicalReaction &rxn,
                                    const ChemicalReaction &query_rxn) {
  if (rxn.getNumAgentTemplates() < query_rxn.getNumAgentTemplates()) {
    return false;
  }
  if (query_rxn.getNumAgentTemplates() == 0) {
    return true;
  }
  return hasReactionMoleculeTemplateSubstructMatch(rxn, query_rxn, Agent);
}

bool hasReactionSubstructMatch(const ChemicalReaction &rxn,
                               const ChemicalReaction &query_rxn,
                               bool includeAgents) {
  if (includeAgents) {
    return (hasReactantTemplateSubstructMatch(rxn, query_rxn) &&
            hasProductTemplateSubstructMatch(rxn, query_rxn) &&
            hasAgentTemplateSubstructMatch(rxn, query_rxn));
  }
  return (hasReactantTemplateSubstructMatch(rxn, query_rxn) &&
          hasProductTemplateSubstructMatch(rxn, query_rxn));
}

bool hasReactionAtomMapping(const ChemicalReaction &rxn) {
  auto begin = getStartIterator(rxn, Reactant);
  auto end = getEndIterator(rxn, Reactant);
  for (; begin != end; ++begin) {
    const ROMol &reactant = *begin->get();
    if (MolOps::getNumAtomsWithDistinctProperty(
            reactant, common_properties::molAtomMapNumber)) {
      return true;
    }
  }
  begin = getStartIterator(rxn, Product);
  end = getEndIterator(rxn, Product);
  for (; begin != end; ++begin) {
    const ROMol &reactant = *begin->get();
    if (MolOps::getNumAtomsWithDistinctProperty(
            reactant, common_properties::molAtomMapNumber)) {
      return true;
    }
  }
  return false;
}

bool isReactionTemplateMoleculeAgent(const ROMol &mol, double agentThreshold) {
  unsigned numMappedAtoms = MolOps::getNumAtomsWithDistinctProperty(
      mol, common_properties::molAtomMapNumber);
  unsigned numAtoms = mol.getNumHeavyAtoms();
  if (numAtoms > 0 &&
      static_cast<double>(numMappedAtoms) / static_cast<double>(numAtoms) >=
          agentThreshold) {
    return false;
  }
  return true;
}

namespace {

void getMappingNumAtomIdxMapReactants(
    const ChemicalReaction &rxn, std::map<int, Atom *> &reactantAtomMapping) {
  for (auto reactIt = rxn.beginReactantTemplates();
       reactIt != rxn.endReactantTemplates(); ++reactIt) {
    for (const auto atom : (*reactIt)->atoms()) {
      int reactMapNum;
      if (atom->getPropIfPresent(common_properties::molAtomMapNumber,
                                 reactMapNum)) {
        reactantAtomMapping[reactMapNum] = atom;
      }
    }
  }
}
}  // namespace

// returns -1 if we don't find a good match
int countSwapsBetweenReactantAndProduct(const ROMol &product,
                                        const Atom *reactAtom,
                                        const Atom *prodAtom) {
  PRECONDITION(reactAtom, "bad atom");
  PRECONDITION(prodAtom, "bad atom");
  if (reactAtom->getDegree() >= 3 && prodAtom->getDegree() >= 3 &&
      std::abs(static_cast<int>(prodAtom->getDegree()) -
               static_cast<int>(reactAtom->getDegree())) <= 1) {
    std::vector<int> reactOrder;
    unsigned int nReactUnmapped = 0;
    for (const auto &nbri : boost::make_iterator_range(
             reactAtom->getOwningMol().getAtomNeighbors(reactAtom))) {
      const auto &nbrAtom = reactAtom->getOwningMol()[nbri];
      if (nbrAtom->getAtomMapNum() > 0) {
        reactOrder.push_back(nbrAtom->getAtomMapNum());
      } else {
        reactOrder.push_back(-1);
        ++nReactUnmapped;
      }
    }
    if (reactAtom->getDegree() < prodAtom->getDegree()) {
      reactOrder.push_back(-1);
      ++nReactUnmapped;
    }
    if (nReactUnmapped <= 1) {
      std::vector<int> prodOrder;
      unsigned int nProdUnmapped = 0;
      for (const auto &nbri :
           boost::make_iterator_range(product.getAtomNeighbors(prodAtom))) {
        const auto &nbrAtom = product[nbri];
        if (nbrAtom->getAtomMapNum() > 0) {
          prodOrder.push_back(nbrAtom->getAtomMapNum());
        } else {
          prodOrder.push_back(-1);
          ++nProdUnmapped;
        }
      }
      if (prodAtom->getDegree() < reactAtom->getDegree()) {
        prodOrder.push_back(-1);
        ++nProdUnmapped;
      }
      if (nProdUnmapped <= 1) {
        // check that each element of the product mappings is
        // in the reactant mappings
        bool allFound = true;
        for (auto poElem : prodOrder) {
          if (poElem >= 0) {
            if (std::find(reactOrder.begin(), reactOrder.end(), poElem) ==
                reactOrder.end()) {
              // this one was not there, is there an unmapped slot for
              // it (i.e. a -1 value in the reactOrder)?
              if (nReactUnmapped) {
                auto negOne =
                    std::find(reactOrder.begin(), reactOrder.end(), -1);
                if (negOne != reactOrder.end()) {
                  *negOne = poElem;
                  --nReactUnmapped;
                } else {
                  allFound = false;
                  break;
                }
              } else {
                allFound = false;
                break;
              }
            }
          }
        }
        if (allFound) {
          // found a match for all the product atoms, what about all
          // the reactant atoms?
          for (auto roElem : reactOrder) {
            if (std::find(prodOrder.begin(), prodOrder.end(), roElem) ==
                prodOrder.end()) {
              if (nProdUnmapped) {
                auto negOne = std::find(prodOrder.begin(), prodOrder.end(), -1);
                if (negOne != prodOrder.end()) {
                  *negOne = roElem;
                  --nProdUnmapped;
                } else {
                  allFound = false;
                  break;
                }
              } else {
                allFound = false;
                break;
              }
            }
          }
        }
        if (allFound) {
          return countSwapsToInterconvert(reactOrder, prodOrder);
        }
      }
    }
  }
  return -1;
}

void updateProductsStereochem(ChemicalReaction *rxn) {
  std::map<int, Atom *> reactantMapping;
  getMappingNumAtomIdxMapReactants(*rxn, reactantMapping);
  for (MOL_SPTR_VECT::const_iterator prodIt = rxn->beginProductTemplates();
       prodIt != rxn->endProductTemplates(); ++prodIt) {
    for (auto prodAtom : (*prodIt)->atoms()) {
      if (prodAtom->hasProp(common_properties::molInversionFlag)) {
        continue;
      }
      if (!prodAtom->hasProp(common_properties::molAtomMapNumber)) {
        // if we have stereochemistry specified, it's automatically
        // creating stereochem:
        prodAtom->setProp(common_properties::molInversionFlag, 4);
        continue;
      }
      int mapNum;
      prodAtom->getProp(common_properties::molAtomMapNumber, mapNum);
      if (reactantMapping.find(mapNum) != reactantMapping.end()) {
        const auto reactAtom = reactantMapping[mapNum];
        if (prodAtom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
            prodAtom->getChiralTag() != Atom::CHI_OTHER) {
          if (reactAtom->getChiralTag() != Atom::CHI_UNSPECIFIED &&
              reactAtom->getChiralTag() != Atom::CHI_OTHER) {
            // both have stereochem specified, we're either preserving
            // or inverting
            if (reactAtom->getChiralTag() == prodAtom->getChiralTag()) {
              prodAtom->setProp(common_properties::molInversionFlag, 2);
            } else {
              // FIX: this is technically fragile: it should be checking
              // if the atoms both have tetrahedral chirality. However,
              // at the moment that's the only chirality available, so
              // there's no need to go monkeying around.
              prodAtom->setProp(common_properties::molInversionFlag, 1);
            }

            // FIX this should move out into a separate function
            // last thing to check here: if the ordering of the bonds
            // around the atom changed from reactants->products then we
            // may need to adjust the inversion flag
            int nSwaps = countSwapsBetweenReactantAndProduct(
                **prodIt, reactAtom, prodAtom);
            if (nSwaps >= 0 && nSwaps % 2) {
              auto mival =
                  prodAtom->getProp<int>(common_properties::molInversionFlag);
              if (mival == 1) {
                prodAtom->setProp(common_properties::molInversionFlag, 2);
              } else if (mival == 2) {
                prodAtom->setProp(common_properties::molInversionFlag, 1);
              } else {
                CHECK_INVARIANT(false, "inconsistent molInversionFlag");
              }
            }
          } else {
            // stereochem in the product, but not in the reactant
            prodAtom->setProp(common_properties::molInversionFlag, 4);
          }
        } else if (reactantMapping[mapNum]->getChiralTag() !=
                       Atom::CHI_UNSPECIFIED &&
                   reactantMapping[mapNum]->getChiralTag() != Atom::CHI_OTHER) {
          // stereochem in the reactant, but not the product:
          prodAtom->setProp(common_properties::molInversionFlag, 3);
        }
      } else {
        // introduction of new stereocenter by the reaction
        prodAtom->setProp(common_properties::molInversionFlag, 4);
      }
    }
  }
}

namespace {

void removeMappingNumbersFromReactionMoleculeTemplate(
    const MOL_SPTR_VECT &molVec) {
  for (const auto &begin : molVec) {
    ROMol &mol = *begin.get();
    for (ROMol::AtomIterator atomIt = mol.beginAtoms();
         atomIt != mol.endAtoms(); ++atomIt) {
      if ((*atomIt)->hasProp(common_properties::molAtomMapNumber)) {
        (*atomIt)->clearProp(common_properties::molAtomMapNumber);
      }
    }
  }
}

}  // namespace

void removeMappingNumbersFromReactions(const ChemicalReaction &rxn) {
  removeMappingNumbersFromReactionMoleculeTemplate(rxn.getAgents());
  removeMappingNumbersFromReactionMoleculeTemplate(rxn.getProducts());
  removeMappingNumbersFromReactionMoleculeTemplate(rxn.getReactants());
}

}  // namespace RDKit
