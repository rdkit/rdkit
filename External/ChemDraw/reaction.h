//
//  Copyright (c) 2024, Glysade Inc
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
#ifndef CHEMDRAW_REACTION_H
#define CHEMDRAW_REACTION_H

#include <GraphMol/RDKitBase.h>

#include "ChemDrawStartInclude.h"
#include "chemdraw/CDXStdObjects.h"
#include "ChemDrawEndInclude.h"

#include <map>
#include <vector>

namespace RDKit {

struct ReactionStepInfo {
  // Holds the current reaction step information so that we can convert
  // chemdraw molecules into rdkit reactions
  unsigned int scheme_id;
  unsigned int step_id;
  std::vector<int> ReactionStepProducts;
  std::vector<int> ReactionStepReactants;
  std::vector<int> ReactionStepObjectsAboveArrow;
  std::vector<int> ReactionStepObjectsBelowArrow;
  std::vector<std::pair<int, int>> ReactionStepAtomMap;

  void set_reaction_data(
      std::string type, std::string prop, const std::vector<int> &frag_ids,
      const std::map<unsigned int, size_t> &fragments,
      std::map<unsigned int, std::vector<int>> &grouped_fragments,
      const std::vector<std::unique_ptr<RWMol>> &mols) const;

  void set_reaction_step(
      size_t scheme_id, std::map<unsigned int, Atom *> &atoms,
      const std::map<unsigned int, size_t> &fragments,
      std::map<unsigned int, std::vector<int>> &grouped_fragments,
      const std::vector<std::unique_ptr<RWMol>> &mols) const;
};

class ReactionInfo {
  // Holds the information form the CDX data so that we can convert
  //  the molecules in the file to RDKit Reactions

  std::vector<ReactionStepInfo> steps;
  unsigned int scheme_id;

 public:
  ReactionInfo(CDXReactionScheme &scheme);

  void set_reaction_steps(
      std::map<unsigned int, std::vector<int>> &grouped_fragments,
      const std::vector<std::unique_ptr<RWMol>> &mols) const;
};
}  // namespace RDKit

#endif
