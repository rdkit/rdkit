// $Id$
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
}

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

void getMappingNumAtomIdxMapReactants(const ChemicalReaction& rxn, std::map<int,Atom::ChiralType>& reactantMapping){
  for (auto reactIt = rxn.beginReactantTemplates();
       reactIt != rxn.endReactantTemplates(); ++reactIt) {
    for(ROMol::AtomIterator reactAtomIt=(*reactIt)->beginAtoms();
        reactAtomIt!=(*reactIt)->endAtoms();++reactAtomIt){
      int reactMapNum = -1;
      (*reactAtomIt)->getPropIfPresent(common_properties::molAtomMapNumber,reactMapNum);
      if(reactMapNum > -1){
        reactantMapping.insert(std::make_pair(reactMapNum,(*reactAtomIt)->getChiralTag()));
      }
    }
  }
}
}

void updateProductsStereochem(ChemicalReaction *rxn) {
  std::map<int, Atom::ChiralType> reactantMapping;
  getMappingNumAtomIdxMapReactants(*rxn, reactantMapping);
  for (MOL_SPTR_VECT::const_iterator prodIt = rxn->beginProductTemplates();
       prodIt != rxn->endProductTemplates(); ++prodIt) {
    for (ROMol::AtomIterator prodAtomIt = (*prodIt)->beginAtoms();
         prodAtomIt != (*prodIt)->endAtoms(); ++prodAtomIt) {
      if ((*prodAtomIt)->hasProp(common_properties::molInversionFlag)) {
        continue;
      }
      if (!(*prodAtomIt)->hasProp(common_properties::molAtomMapNumber)) {
        // if we have stereochemistry specified, it's automatically creating
        // stereochem:
        (*prodAtomIt)->setProp(common_properties::molInversionFlag, 4);
        continue;
      }
      int mapNum;
      (*prodAtomIt)->getProp(common_properties::molAtomMapNumber, mapNum);
      if (reactantMapping.find(mapNum) != reactantMapping.end()) {
        if ((*prodAtomIt)->getChiralTag() != Atom::CHI_UNSPECIFIED &&
            (*prodAtomIt)->getChiralTag() != Atom::CHI_OTHER) {
          if (reactantMapping[mapNum] != Atom::CHI_UNSPECIFIED &&
              reactantMapping[mapNum] != Atom::CHI_OTHER) {
            // both have stereochem specified:
            if (reactantMapping[mapNum] == (*prodAtomIt)->getChiralTag()) {
              (*prodAtomIt)->setProp(common_properties::molInversionFlag, 2);
            } else {
              // FIX: this is technically fragile: it should be checking
              // if the atoms both have tetrahedral chirality. However,
              // at the moment that's the only chirality available, so there's
              // no need to go monkeying around.
              (*prodAtomIt)->setProp(common_properties::molInversionFlag, 1);
            }
          } else {
            // stereochem in the product, but not in the reactant
            (*prodAtomIt)->setProp(common_properties::molInversionFlag, 4);
          }
        } else if (reactantMapping[mapNum] != Atom::CHI_UNSPECIFIED &&
                   reactantMapping[mapNum] != Atom::CHI_OTHER) {
          // stereochem in the reactant, but not the product:
          (*prodAtomIt)->setProp(common_properties::molInversionFlag, 3);
        }
      } else {
        // introduction of new stereocenter by the reaction
        (*prodAtomIt)->setProp(common_properties::molInversionFlag, 4);
      }
    }
  }
}

namespace {

void removeMappingNumbersFromReactionMoleculeTemplate(const MOL_SPTR_VECT &molVec){
  for (const auto &begin : molVec) {
    ROMol &mol = *begin.get();
    for(ROMol::AtomIterator atomIt=mol.beginAtoms();
        atomIt!=mol.endAtoms();++atomIt){
      if((*atomIt)->hasProp(common_properties::molAtomMapNumber)){
        (*atomIt)->clearProp(common_properties::molAtomMapNumber);
      }
    }
  }
}

}

void removeMappingNumbersFromReactions(const ChemicalReaction &rxn){

  removeMappingNumbersFromReactionMoleculeTemplate(rxn.getAgents());
  removeMappingNumbersFromReactionMoleculeTemplate(rxn.getProducts());
  removeMappingNumbersFromReactionMoleculeTemplate(rxn.getReactants());
}

} // end of RDKit namespace

