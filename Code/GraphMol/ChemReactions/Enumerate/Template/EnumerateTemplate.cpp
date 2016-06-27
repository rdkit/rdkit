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

#include "EnumerateTemplate.h"
#include "../EnumerationStrategyBase.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <string>
#include <vector>

#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit
{

namespace {
  // XXX Fix me
  //  Making this temporary is silly, we should store the template info
  //    on the side-chain molecule as generic prop data

void templatesToBBs(BBS &templateBBs, const ReactantTemplates &m_reactants) {
  templateBBs.clear();
  templateBBs.resize(m_reactants.m_templates.size());
  for(size_t i=0;i<m_reactants.m_templates.size();++i) {
    templateBBs[i].resize(m_reactants.m_templates[i].size());
    for(size_t bbidx=0;bbidx<m_reactants.m_templates[i].size(); ++bbidx) {
      templateBBs[i][bbidx] = m_reactants.m_templates[i][bbidx].sideChainMol;
    }
  }
}
}
  
EnumerateTemplateLibrary::EnumerateTemplateLibrary(const ChemicalReaction &rxn,
                                                   const BBS &bbs) :
    EnumerateLibraryBase(rxn, new CartesianProductStrategy)
{
  bool success = ReactantsToTemplates(m_reactants, m_rxn, bbs);
  PRECONDITION(success, "Failed to create reactant templates");
  BBS templateBBs;
  templatesToBBs(templateBBs, m_reactants);
  m_enumerator->initialize(rxn, templateBBs);
  PRECONDITION(static_cast<bool>(*m_enumerator), "empty enumerator");
}

EnumerateTemplateLibrary::EnumerateTemplateLibrary(
    const ChemicalReaction &rxn,
    const BBS &bbs,
    const EnumerationStrategyBase &enumerator) :
    EnumerateLibraryBase(rxn) {
  m_enumerator.reset( enumerator.Clone() );
  bool success = ReactantsToTemplates(m_reactants, m_rxn, bbs);
  PRECONDITION(success, "Failed to create reactant templates");
  BBS templateBBs;
  templatesToBBs(templateBBs, m_reactants);
  m_enumerator->initialize(rxn, templateBBs);
  PRECONDITION(static_cast<bool>(*m_enumerator), "empty enumerator");
}

EnumerateTemplateLibrary::EnumerateTemplateLibrary(const EnumerateTemplateLibrary &rhs) :
    EnumerateLibraryBase(rhs), 
    m_reactants(rhs.m_reactants) {
}
    
std::vector<MOL_SPTR_VECT> EnumerateTemplateLibrary::next() {
  PRECONDITION( static_cast<bool>(*this), "No more enumerations");
  const bool sanitize=false;
  const std::vector<std::vector<std::string> > &smiles = nextSmiles();
  std::vector<MOL_SPTR_VECT> res(smiles.size());
  for (size_t i=0; i<smiles.size(); ++i) {
    for(size_t j=0; j<smiles[i].size(); ++j) {
      ROMOL_SPTR mol(SmilesToMol(smiles[i][j],0,sanitize));
      res[i].push_back(mol);
    }
  }
  return res;
}

      
std::vector<std::vector<std::string> > EnumerateTemplateLibrary::nextSmiles() {
  PRECONDITION( static_cast<bool>(*this), "No more enumerations");
  const RGROUPS &reactantIndices = m_enumerator->next();
  
  std::vector<std::vector<std::string> > v(1);
  v[0].push_back(m_reactants.smiles( reactantIndices ));
  return v;
}

void EnumerateTemplateLibrary::toStream(std::ostream &ss, bool enumerationStateOnly) const {
  boost::archive::binary_oarchive ar(ss);
  ar << *this;
}
  
void EnumerateTemplateLibrary::initFromStream(std::istream &ss)
{
  boost::archive::binary_iarchive ar(ss);
  ar >> *this;
}
}
