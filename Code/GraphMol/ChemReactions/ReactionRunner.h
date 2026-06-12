//
//  Copyright (c) 2014-2021, Novartis Institutes for BioMedical Research Inc.
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

#include <RDGeneral/export.h>
#ifndef RD_REACTION_RUNNER_H
#define RD_REACTION_RUNNER_H

#include <map>
#include <tuple>

#include <GraphMol/ChemReactions/Reaction.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/Substruct/SubstructMatch.h>

namespace RDKit {
using VectMatchVectType = std::vector<MatchVectType>;
//! Caches reactant-template substructure matches across run_Reactants() calls.
/*!
  Keyed on (reactant template index, reagent pointer identity, effective
  maxMatches, dedupeSymmetricMatches). The ROMol objects whose raw pointers
  are used as keys must outlive this cache and must not be structurally
  modified between calls, otherwise stale match data may be returned.
  Not thread-safe.
*/
using ReactantMatchCache =
    std::map<std::tuple<unsigned int, const ROMol *, unsigned int, bool>,
             VectMatchVectType>;

//! Runs the reaction on a set of reactants
/*!
  \param rxn:       the template reaction we are interested
  \param reactants: the reactants to be used. The length of this must be equal
  to
                    rxn->getNumReactantTemplates()
                    Caution: The order of the reactant templates determines the
  order of the reactants!
  \param maxProducts:  if non zero, the maximum number of products to generate
                       before stopping.  If hit a warning will be generated.
  \return a vector of vectors of products. Each subvector will be
          rxn->getNumProductTemplates() long.

  We return a vector of vectors of products because each individual template may
  map multiple times onto its reactant. This leads to multiple possible result
  sets.
*/
RDKIT_CHEMREACTIONS_EXPORT std::vector<MOL_SPTR_VECT> run_Reactants(
    const ChemicalReaction &rxn, const MOL_SPTR_VECT &reactants,
    unsigned int maxProducts = 1000);

//! Runs a reaction on a set of reactants using a reactant match cache
/*!
  This overload reuses previously computed reactant-template matches for
  repeated reagent/template combinations. The cache key is the reactant
  template index, the reagent pointer identity, the effective
  maxMatches/maxProducts value, and whether symmetry deduplication is on.
*/
RDKIT_CHEMREACTIONS_EXPORT std::vector<MOL_SPTR_VECT> run_Reactants(
    const ChemicalReaction &rxn, const MOL_SPTR_VECT &reactants,
    ReactantMatchCache &cache, bool dedupeSymmetricMatches = false,
    unsigned int maxProducts = 1000);

//! Runs a single reactant against a single reactant template
/*!
  \param reactant The single reactant to use

  \param reactantTemplateIdx the reactant template to target in the reaction

  \return a vector of vectors of products. Each subvector will be
          rxn->getNumProductTemplates() long.

  We return a vector of vectors of products because each individual template may
  map multiple times onto its reactant. This leads to multiple possible result
  sets.

*/

RDKIT_CHEMREACTIONS_EXPORT std::vector<MOL_SPTR_VECT> run_Reactant(
    const ChemicalReaction &rxn, const ROMOL_SPTR &reactant,
    unsigned int reactantIdx);

RDKIT_CHEMREACTIONS_EXPORT bool run_Reactant(const ChemicalReaction &rxn,
                                             RWMol &reactant,
                                             bool removeUnmatchedAtoms = true);

//! Reduce the product generated by run_Reactants or run_Reactant to
/// the sidechains that come from the reagents
//
//  n.b. molecules that might be a product of the given reaction
//       but were not generated by run_Reactant(s) currently
//       produce no sidechains.
/*!
  \param addDummyAtoms If true, add dummy atoms to the sidechains for the
      non-reagent parts of the sidechain.  Dummy atoms are annotated with
      the atom maps from the reaction.
      If False, then any sidechain atom where a bond was cleaved is annotated
  with: _rgroupAtomMaps property which indicates the scaffold atommaps that
  where bonded _rgroupBonds property which indicates the bondtype for each
  atommap bonded
*/

RDKIT_CHEMREACTIONS_EXPORT ROMol *reduceProductToSideChains(
    const ROMOL_SPTR &product, bool addDummyAtoms = true);

namespace ReactionRunnerUtils {
struct ReactantGraft {
  RWMOL_SPTR iso;
  unsigned int templateAtomCount;
  unsigned int templateBondCount;
  std::vector<unsigned int> anchorAtomIdxs;
  std::vector<unsigned int> anchorBondIdxs;
};

RDKIT_CHEMREACTIONS_EXPORT VectMatchVectType getReactantMatchesToTemplate(
  const ROMol &reactant, const ROMol &templ, unsigned int maxMatches,
  const SubstructMatchParameters &ssparams);

RDKIT_CHEMREACTIONS_EXPORT VectMatchVectType dedupeMatchesBySymmetry(
    const ROMol &reactant, const VectMatchVectType &matches);

RDKIT_CHEMREACTIONS_EXPORT MOL_SPTR_VECT generateOneProductSet(
    const ChemicalReaction &rxn, const MOL_SPTR_VECT &reactants,
    const std::vector<MatchVectType> &reactantsMatch);

RDKIT_CHEMREACTIONS_EXPORT RWMOL_SPTR
convertTemplateToMol(const ROMOL_SPTR prodTemplateSptr);

RDKIT_CHEMREACTIONS_EXPORT ReactantGraft extractReactantGraft(
  const ChemicalReaction &rxn, const ROMOL_SPTR &productTemplate,
  const ROMOL_SPTR &reactantTemplate, const ROMOL_SPTR &reactant,
  const MatchVectType &match, unsigned int reactantId, bool doConfs);

RDKIT_CHEMREACTIONS_EXPORT void applyReactantGraft(RWMOL_SPTR product,
                           Conformer *productConf,
                           const ReactantGraft &graft);
}  // namespace ReactionRunnerUtils

}  // namespace RDKit

#endif
