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

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ChemReactions/ReactionUtils.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Descriptors/MolDescriptors.h>

namespace RDKit {

MOL_SPTR_VECT::const_iterator getStartIterator(
   const ChemicalReaction &rxn,
   ReactionMoleculeType t)
{
  MOL_SPTR_VECT::const_iterator begin;
  if(t == Reactant){
    begin = rxn.beginReactantTemplates();
  }
  if(t == Product){
    begin = rxn.beginProductTemplates();;
  }
  if(t == Agent){
    begin = rxn.beginAgentTemplates();
  }
  return begin;
}

MOL_SPTR_VECT::const_iterator getEndIterator(
   const ChemicalReaction &rxn,
   ReactionMoleculeType t)
{
  MOL_SPTR_VECT::const_iterator end;
  if(t == Reactant){
	end = rxn.endReactantTemplates();
  }
  if(t == Product){
	end = rxn.endProductTemplates();;
  }
  if(t == Agent){
    end = rxn.endAgentTemplates();
  }
  return end;
}


namespace {

unsigned getReactionMoleculeTemplatesNumAtoms(
    const RDKit::ChemicalReaction &rxn,
    RDKit::ReactionMoleculeType t,
    bool onlyHeavy)
{
  unsigned numAtoms=0;
  RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, t);
  RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, t);

  for(; begin != end; ++begin){
	if(onlyHeavy){
      numAtoms += begin->get()->getNumHeavyAtoms();
	}
	else{
      numAtoms += begin->get()->getNumAtoms();
	}
  }
  return numAtoms;
}

unsigned getReactionMoleculeTemplatesNumBonds(
    const RDKit::ChemicalReaction &rxn,
    RDKit::ReactionMoleculeType t,
    bool onlyHeavy)
{
  unsigned numBonds=0;
  RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, t);
  RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, t);

  for(; begin != end; ++begin){
    numBonds += begin->get()->getNumBonds(onlyHeavy);
  }
  return numBonds;
}

double getReactionMoleculeTemplatesMW(
    const RDKit::ChemicalReaction &rxn,
    RDKit::ReactionMoleculeType t,
    bool onlyHeavy)
{
  double mw=0.0;
  RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, t);
  RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, t);

  for(; begin != end; ++begin){
    mw += RDKit::Descriptors::calcAMW(*begin->get(),onlyHeavy);
  }
  return mw;
}

unsigned getReactionMoleculeTemplatesNumRings(
    const RDKit::ChemicalReaction &rxn,
    RDKit::ReactionMoleculeType t)
{
  unsigned numRings=0;
  RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, t);
  RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, t);

  for(; begin != end; ++begin){
	if(begin->get()->getRingInfo()->isInitialized()){
	  numRings += begin->get()->getRingInfo()->numRings();
	}
	else{
	  begin->get()->updatePropertyCache();
	  RDKit::MolOps::findSSSR(*begin->get());
	  numRings += begin->get()->getRingInfo()->numRings();
	}
  }
  return numRings;
}

bool hasReactionMoleculeTemplateSubstructMatch(
  const RDKit::ChemicalReaction &rxn,
  const RDKit::ChemicalReaction &query_rxn,
  RDKit::ReactionMoleculeType t)
{
  for(RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, t);
		  begin != getEndIterator(rxn, t); ++begin){
    for(RDKit::MOL_SPTR_VECT::const_iterator begin_query = getStartIterator(query_rxn, t);
    	  begin_query != getEndIterator(query_rxn, t); ++begin_query){
	  MatchVectType tvect;
	  if(SubstructMatch(*begin->get(), *begin_query->get(),tvect )){
		return true;
	  }
	}
  }
  return false;
}

}

unsigned getReactantTemplatesNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesNumAtoms(rxn, Reactant, onlyHeavy);
}

unsigned getProductTemplatesNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesNumAtoms(rxn, Product, onlyHeavy);
}

unsigned getAgentTemplatesNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesNumAtoms(rxn, Agent, onlyHeavy);
}

unsigned getReactionNumAtoms(const ChemicalReaction &rxn, bool onlyHeavy, bool includeAgents)
{
  if(includeAgents){
    return (getReactantTemplatesNumAtoms(rxn, onlyHeavy) +
    		 getProductTemplatesNumAtoms(rxn, onlyHeavy) +
    		 getAgentTemplatesNumAtoms(rxn, onlyHeavy));
  }
  return (getReactantTemplatesNumAtoms(rxn, onlyHeavy) +
		  getProductTemplatesNumAtoms(rxn, onlyHeavy));

}

unsigned getReactantTemplatesNumBonds(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesNumBonds(rxn, Reactant, onlyHeavy);
}

unsigned getProductTemplatesNumBonds(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesNumBonds(rxn, Product, onlyHeavy);
}

unsigned getAgentTemplatesNumBonds(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesNumBonds(rxn, Agent, onlyHeavy);
}

unsigned getReactionNumBonds(const ChemicalReaction &rxn, bool onlyHeavy, bool includeAgents)
{
  if(includeAgents){
    return (getReactantTemplatesNumBonds(rxn, onlyHeavy) +
    		 getProductTemplatesNumBonds(rxn, onlyHeavy) +
    		 getAgentTemplatesNumBonds(rxn, onlyHeavy));
  }
  return (getReactantTemplatesNumBonds(rxn, onlyHeavy) +
		  getProductTemplatesNumBonds(rxn, onlyHeavy));
}


double calcReactantTemplatesMW(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesMW(rxn, Reactant, onlyHeavy);
}

double calcProductTemplatesMW(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesMW(rxn, Product, onlyHeavy);
}

double calcAgentTemplatesMW(const ChemicalReaction &rxn, bool onlyHeavy)
{
  return getReactionMoleculeTemplatesMW(rxn, Agent, onlyHeavy);
}

double calcReactionMW(const ChemicalReaction &rxn, bool onlyHeavy, bool includeAgents)
{
  if(includeAgents){
    return(calcReactantTemplatesMW(rxn, onlyHeavy) +
    		calcProductTemplatesMW(rxn, onlyHeavy) +
    		calcAgentTemplatesMW(rxn, onlyHeavy));
  }
  return(calcReactantTemplatesMW(rxn, onlyHeavy) +
		 calcProductTemplatesMW(rxn, onlyHeavy));
}

unsigned getReactantTemplatesNumRings(const ChemicalReaction &rxn)
{
  return getReactionMoleculeTemplatesNumRings(rxn, Reactant);
}

unsigned getProductTemplatesNumRings(const ChemicalReaction &rxn)
{
  return getReactionMoleculeTemplatesNumRings(rxn, Product);
}

unsigned getAgentTemplatesNumRings(const ChemicalReaction &rxn)
{
  return getReactionMoleculeTemplatesNumRings(rxn, Agent);
}

unsigned getReactionNumRings(const ChemicalReaction &rxn, bool includeAgents)
{
  if(includeAgents){
    return(getReactantTemplatesNumRings(rxn) +
    		getProductTemplatesNumRings(rxn) +
    		getAgentTemplatesNumRings(rxn));
  }
  return(getReactantTemplatesNumRings(rxn) +
		  getProductTemplatesNumRings(rxn));
}

bool hasReactantTemplateSubstructMatch(
  const ChemicalReaction &rxn,
  const ChemicalReaction &query_rxn)
{
  if(rxn.getNumReactantTemplates() < query_rxn.getNumReactantTemplates()){
    return false;
  }
  if(query_rxn.getNumReactantTemplates() == 0){
  	return true;
  }
  return hasReactionMoleculeTemplateSubstructMatch(rxn, query_rxn, Reactant);
}

bool hasProductTemplateSubstructMatch(
  const ChemicalReaction &rxn,
  const ChemicalReaction &query_rxn)
{
  if(rxn.getNumProductTemplates() < query_rxn.getNumProductTemplates()){
  	return false;
  }
  if(query_rxn.getNumProductTemplates() == 0){
  	return true;
  }
  return hasReactionMoleculeTemplateSubstructMatch(rxn, query_rxn, Product);
}

bool hasAgentTemplateSubstructMatch(
  const ChemicalReaction &rxn,
  const ChemicalReaction &query_rxn)
{
  if(rxn.getNumAgentTemplates() < query_rxn.getNumAgentTemplates()){
	return false;
  }
  if(query_rxn.getNumAgentTemplates() == 0){
	return true;
  }
  return hasReactionMoleculeTemplateSubstructMatch(rxn, query_rxn, Agent);
}

bool hasReactionSubstructMatch(
  const ChemicalReaction &rxn,
  const ChemicalReaction &query_rxn, bool includeAgents)
{
  if(includeAgents){
    return (hasReactantTemplateSubstructMatch(rxn, query_rxn) &&
    	  	  hasProductTemplateSubstructMatch(rxn, query_rxn) &&
  	  	      hasAgentTemplateSubstructMatch(rxn, query_rxn));
  }
  return (hasReactantTemplateSubstructMatch(rxn, query_rxn) &&
	  	  hasProductTemplateSubstructMatch(rxn, query_rxn));
}

bool hasReactionAtomMapping(const ChemicalReaction &rxn)
{
  RDKit::MOL_SPTR_VECT::const_iterator begin = getStartIterator(rxn, Reactant);
  RDKit::MOL_SPTR_VECT::const_iterator end = getEndIterator(rxn, Reactant);
  for(; begin != end; ++begin){
    const ROMol &reactant = *begin->get();
    if(MolOps::getNumAtomsWithDistinctProperty(reactant, "molAtomMapNumber")){
      return true;
    }
  }
  begin = getStartIterator(rxn, Product);
  end = getEndIterator(rxn, Product);
  for(; begin != end; ++begin){
    const ROMol &reactant = *begin->get();
    if(MolOps::getNumAtomsWithDistinctProperty(reactant, "molAtomMapNumber")){
      return true;
    }
  }
  return false;
}

bool isReactionTemplateMoleculeAgent(const ROMol &mol, double agentThreshold)
{
  unsigned numMappedAtoms = MolOps::getNumAtomsWithDistinctProperty(mol, "molAtomMapNumber");
  unsigned numAtoms = mol.getNumHeavyAtoms();
  if(numAtoms > 0 &&
       static_cast<double>(numMappedAtoms)/static_cast<double>(numAtoms) >= agentThreshold){
    return false;
  }
  return true;
}

} // end of RDKit namespace
