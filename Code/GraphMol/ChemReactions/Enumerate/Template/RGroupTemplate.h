//
//  Copyright (c) 2015, Novartis Institutes for BioMedical Research Inc.
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

#ifndef RDKIT_RGROUP_TEMPLATE
#define RDKIT_RGROUP_TEMPLATE
#include <GraphMol/RDKitBase.h>
#include "../EnumerateBase.h"

#include <string>
#include <vector>
#include <boost/shared_ptr.hpp>

namespace RDKit {
namespace PropTypes {
  const unsigned int MW = 0;
  const unsigned int TPSA = 1;
  const unsigned int ALOGP = 2;
  const unsigned int NUM_ROTORS = 3;
  const unsigned int UNASSIGNED_STEREO = 4;
  const unsigned int MAX = 5;
}

// Helper template to hold implementation details
struct RGroupTemplate {
  ROMOL_SPTR               sideChainMol;
  std::string              smiles;         // A smiles-like template
  std::vector<int>         mappings;
  std::vector<int>         counts;
  unsigned int             reactantIdx;
  unsigned int             bbIdx;
  std::vector<double>      props;
  
 RGroupTemplate() : sideChainMol(), smiles(),
    mappings(), counts(), reactantIdx(-1), bbIdx(-1),
    props(PropTypes::MAX) {}
  
  bool initialize(ROMOL_SPTR &product,
                  unsigned int reactantIdx,
                  unsigned int bbIdx);
  
  int addReactionMapping(unsigned int mapno);
  
  bool isValid() const;       // true of the template is valid
};  

typedef std::vector<std::vector<RGroupTemplate> > VectVectRGroupTemplate;

struct ReactantTemplates {
  std::string               m_scaffoldSmiles;
  std::vector<unsigned int> m_mappings;
  VectVectRGroupTemplate    m_templates;
  
  std::string smiles(const RGROUPS &reactantIds) const;
};

bool ReactantsToTemplates(ReactantTemplates &templates,
                          const ChemicalReaction &rxn,
                          const BBS &bbs);


}
#endif
