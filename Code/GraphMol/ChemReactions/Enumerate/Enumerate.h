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
#ifndef RDKIT_ENUMERATE_H
#define RDKIT_ENUMERATE_H
#include "EnumerateBase.h"

namespace RDKit {
//! This is a class for running reactions on sets of reagents.
/*!
  This class is a fully self contained reaction engine that can be
  serialized and restarted.  For example, a million products can
  be generated, the engine can be saved for later and reloaded
  to retrieve the next million products.

  basic usage will be something like:
  \verbatim
   ChemicalReaction rxn = ...
   BBS bbs(num_rgroups);
   LoadRGroups(bbs[0]);
   LoadRGroups(bbs[1]..);
   ...
   EnumerateLibrary enumerator(en, bbs);
   for(; (bool)en; ++i) {
     // This is the same as rxn.run_Reactants( reagents );
     std::vector<MOL_SPTR_VECT> products = en.next();
     ...
   }
   \endverbatim

   In general, reactions will enumerate to more products than desired,
   a standard use is:

   \verbatim
   for(int i=0;i<num_samples && (bool)en; ++i) {
     std::vector<MOL_SPTR_VECT> products = en.next();
     ...
   }
   \endverbatim
 */
class EnumerateLibrary : public EnumerateLibraryBase {
  BBS m_bbs;
  
 public:
  EnumerateLibrary() : EnumerateLibraryBase(), m_bbs() {}
  EnumerateLibrary(const ChemicalReaction &rxn, const BBS &bbs);
  EnumerateLibrary(const ChemicalReaction &rxn, const BBS &bbs,
                   const EnumerationStrategyBase &enumerator);
  EnumerateLibrary(const EnumerateLibrary &rhs);

  //! Return the reagents used in the library
  const BBS &getReagents() const { return m_bbs; }

  //! Get the next product set
  std::vector<MOL_SPTR_VECT> next();

  void toStream(std::ostream &ss) const;
  void initFromStream(std::istream &ss);

 private:
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive &ar, const unsigned int /*version*/) const {
    ar &boost::serialization::base_object<EnumerateLibraryBase>(*this);
    size_t sz = m_bbs.size();
    ar &sz;
    
    std::string pickle;
    for (size_t i = 0; i < m_bbs.size(); ++i) {
      sz = m_bbs[i].size();
      ar &sz;
      for (size_t j = 0; j < m_bbs[i].size(); ++j) {
        MolPickler::pickleMol(*m_bbs[i][j], pickle);
        ar &pickle;
      }
    }
  }
  template <class Archive>
  void load(Archive &ar, const unsigned int /*version*/) {
    ar &boost::serialization::base_object<EnumerateLibraryBase>(*this);
    
    size_t sz;
    ar &sz;

    m_bbs.resize(sz);

    for (size_t i = 0; i < m_bbs.size(); ++i) {
      ar &sz;
      m_bbs[i].resize(sz);
      std::string pickle;
      for (size_t j = 0; j < m_bbs[i].size(); ++j) {
        ar &pickle;
        RWMol *mol = new RWMol();
        MolPickler::molFromPickle(pickle, *mol);
        m_bbs[i][j].reset(mol);
      }
    }
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER();
};
}
#endif
