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
#ifndef RDKIT_ENUMERATE_TEMPLATE_H
#define RDKIT_ENUMERATE_TEMPLATE_H
#include "../EnumerateBase.h"
#include "RGroupTemplate.h"
#include <GraphMol/RDKitBase.h>

namespace RDKit
{

class EnumerateTemplateLibrary : public EnumerateLibraryBase
{
public:
  ReactantTemplates m_reactants;
  
public:
EnumerateTemplateLibrary() : EnumerateLibraryBase(), m_reactants() {}
  EnumerateTemplateLibrary(const ChemicalReaction &rxn,
                           const BBS &bbs);
  EnumerateTemplateLibrary(const ChemicalReaction &rxn,
                           const BBS &bbs,
                           const EnumerationStrategyBase &enumerator);
  EnumerateTemplateLibrary(const EnumerateTemplateLibrary &rhs);
  
  virtual std::vector<MOL_SPTR_VECT> next();
  virtual std::vector<std::vector<std::string> > nextSmiles();
  
  void toStream(std::ostream &ss, bool enumerationStateOnly=false) const;
  void initFromStream(std::istream &ss);
  
private:
  friend class boost::serialization::access;
  template<class Archive>    
      void serialize (Archive &ar, const unsigned int /*version*/) {
    ar & boost::serialization::base_object<EnumerateLibraryBase>(*this);
    //ar & m_reactants;
  }
};
}
#endif
