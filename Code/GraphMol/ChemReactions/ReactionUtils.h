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
//       products derived from this software without specific prior written permission.
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

#ifndef __RD_REACTION_UTILS_H
#define __RD_REACTION_UTILS_H

#include <GraphMol/ChemReactions/Reaction.h>

namespace RDKit{

  enum ReactionMoleculeType{
  	Reactant,
  	Product,
  	Agent
  };

  MOL_SPTR_VECT::const_iterator getStartIterator(
     const ChemicalReaction &rxn,
     ReactionMoleculeType t);
  MOL_SPTR_VECT::const_iterator getEndIterator(
     const ChemicalReaction &rxn,
     ReactionMoleculeType t);

  unsigned getReactantTemplatesNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy=true);
  unsigned getProductTemplatesNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy=true);
  unsigned getAgentTemplatesNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy=true);
  unsigned getReactionNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy=true, bool includeAgents=false);
  
  unsigned getReactantTemplatesNumBonds(const ChemicalReaction &rxn, bool onlyHeavy=true);
  unsigned getProductTemplatesNumBonds(const ChemicalReaction &rxn, bool onlyHeavy=true);
  unsigned getAgentTemplatesNumBonds(const ChemicalReaction &rxn, bool onlyHeavy=true);
  unsigned getReactionNumBonds(const ChemicalReaction &rxn, bool onlyHeavy=true, bool includeAgents=false);
  
  double calcReactantTemplatesMW(const ChemicalReaction &rxn, bool onlyHeavy=true);
  double calcProductTemplatesMW(const ChemicalReaction &rxn, bool onlyHeavy=true);
  double calcAgentTemplatesMW(const ChemicalReaction &rxn, bool onlyHeavy=true);
  double calcReactionMW(const ChemicalReaction &rxn, bool onlyHeavy=true, bool includeAgents=false);

  unsigned getReactantTemplatesNumRings(const ChemicalReaction &rxn);
  unsigned getProductTemplatesNumRings(const ChemicalReaction &rxn);
  unsigned getAgentTemplatesNumRings(const ChemicalReaction &rxn);
  unsigned getReactionNumRings(const ChemicalReaction &rxn, bool includeAgents=false);

  bool hasReactantTemplateSubstructMatch(
    const ChemicalReaction &rxn,
    const ChemicalReaction &query_rxn);

  bool hasProductTemplateSubstructMatch(
    const ChemicalReaction &rxn,
    const ChemicalReaction &query_rxn);

  bool hasAgentTemplateSubstructMatch(
    const ChemicalReaction &rxn,
    const ChemicalReaction &query_rxn);

  bool hasReactionSubstructMatch(
    const ChemicalReaction &rxn,
    const ChemicalReaction &query_rxn, bool includeAgents=false);

  bool hasReactionAtomMapping(const ChemicalReaction &rxn);

  bool isReactionTemplateMoleculeAgent(const ROMol &mol, double agentThreshold);

} // end of RDKit namespace

#endif
