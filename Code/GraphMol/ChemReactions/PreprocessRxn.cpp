//
//  Copyright (c) 2016, Novartis Institutes for BioMedical Research Inc.
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

#include "PreprocessRxn.h"
#include "ReactionParser.h"
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/shared_ptr.hpp>
#include <boost/thread/once.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include <GraphMol/FilterCatalog/FunctionalGroupHierarchy.h>

namespace RDKit {

bool preprocessReaction(ChemicalReaction &rxn,
                        const std::string &propName)
{
  const bool normalized=true;
  return preprocessReaction(rxn,
                           GetFlattenedFunctionalGroupHierarchy(normalized),
                           propName);
}

bool preprocessReaction(ChemicalReaction &rxn,
                        unsigned int &numWarnings,
                        unsigned int &numErrors,
                        std::vector<
                          std::vector<std::pair<unsigned int,std::string> > >&reactantLabels,                        
                        const std::string &propName)
{
  const bool normalized = true;
  return preprocessReaction(rxn,
                            numWarnings,
                            numErrors,
                            reactantLabels,
                            GetFlattenedFunctionalGroupHierarchy(normalized),
                            propName);
}

bool preprocessReaction(ChemicalReaction &rxn,
                        const std::map<std::string, ROMOL_SPTR> &queries,
                        const std::string &propName) {
  unsigned int numWarnings, numErrors;
  std::vector<
    std::vector<std::pair<unsigned int,std::string> > >reactantLabels;

  return preprocessReaction(rxn,
                            numWarnings,
                            numErrors,
                            reactantLabels,
                            queries,
                            propName);
}

bool preprocessReaction(ChemicalReaction &rxn,
                        unsigned int &numWarnings,
                        unsigned int &numErrors,
                        std::vector<
                          std::vector<std::pair<unsigned int,std::string> > >&reactantLabels,
                        const std::map<std::string, ROMOL_SPTR> &queries,
                        const std::string &propName) {
  rxn.setImplicitPropertiesFlag(true);
  rxn.initReactantMatchers();

  if (rxn.validate(numWarnings, numErrors)) {
    addRecursiveQueriesToReaction(rxn,
                                  queries,
                                  propName,
                                  &reactantLabels);
    return true;
  }

  return false;
}

}
