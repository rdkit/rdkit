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

#ifndef RD_REACTION_RUNNER_H
#define RD_REACTION_RUNNER_H

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ROMol.h>

namespace RDKit {
//! Runs the reaction on a set of reactants
/*!
  \param rxn:       the template reaction we are interested
  \param reactants: the reactants to be used. The length of this must be equal
  to
                    rxn->getNumReactantTemplates()
                    Caution: The order of the reactant templates determines the
  order of the reactants!

  \return a vector of vectors of products. Each subvector will be
          rxn->getNumProductTemplates() long.

  We return a vector of vectors of products because each individual template may
  map multiple times onto its reactant. This leads to multiple possible result
  sets.
*/
std::vector<MOL_SPTR_VECT> run_Reactants(const ChemicalReaction& rxn,
                                         const MOL_SPTR_VECT& reactants);

//! Runs a single reactant against a single reactant template
/*!
  \param reactant The single reactant to use

  \param reactantTemplateIdx the reactant template to target in the reaction
*/

std::vector<MOL_SPTR_VECT> run_Reactant(const ChemicalReaction& rxn,
                                        const ROMOL_SPTR& reactant,
                                        unsigned int reactantIdx);

//! Reduce the product of a reaction to the sidechains that come from the
//reagents
/*!
  \param addDummyAtoms If true, add dummy atoms to the sidechains for the
      non-reagent parts of the sidechain.  Dummy atoms are annotated with
      the atom maps from the reaction.
      If False, then any sidechain atom where a bond was cleaved is annotated with:
         _rgroupAtomMaps property which indicates the scaffold atommaps that where bonded
         _rgroupBonds property which indicates the bondtype for each atommap bonded
*/

ROMol* reduceProductToSideChains(const ROMOL_SPTR& product,
                                 bool addDummyAtoms = true);
}  // end of RDKit namespace

#endif
