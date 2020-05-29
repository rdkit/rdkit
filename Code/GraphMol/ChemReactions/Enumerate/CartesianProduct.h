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

#include <RDGeneral/export.h>
#ifndef CARTESIANPRODUCT_H
#define CARTESIANPRODUCT_H

#include "EnumerationStrategyBase.h"

namespace RDKit {
//! This is a class for enumerating reagents using Cartesian Products of
// reagents.
/*!
  CartesianProductStrategy produces a  standard walk through all possible
  reagent combinations:

   (0,0,0), (1,0,0), (2,0,0) ...

  basic usage:

  \verbatim
  std::vector<MOL_SPTR_VECT> bbs;
  bbs.push_back( bbs_for_reactants_1 );
  bbs.push_back( bbs_for_reactants_2 );

  RGRUOPS num_bbs;
  num_bbs.push_back(bbs[0].size());
  num_bbs.push_back(bbs[1].size());

  CartesianProductStrategy rgroups(num_bbs);
  for(size_t i=0; i<num_samples && rgroups; ++i) {
    MOL_SPTR_VECT rvect = getReactantsFromRGroups(bbs, rgroups.next());
    std::vector<MOL_SPTR_VECT> lprops = rxn.RunReactants(rvect);
    ...
  }
  \endverbatim

See EnumerationStrategyBase for more details and usage.
*/

class RDKIT_CHEMREACTIONS_EXPORT CartesianProductStrategy
    : public EnumerationStrategyBase {
  size_t m_numPermutationsProcessed{};

 public:
  CartesianProductStrategy()
      : EnumerationStrategyBase() {}

  using EnumerationStrategyBase::initialize;

  virtual void initializeStrategy(const ChemicalReaction &,
                                  const EnumerationTypes::BBS &) {
    m_numPermutationsProcessed = 0;
  }

  virtual const char *type() const { return "CartesianProductStrategy"; }

  //! The current permutation {r1, r2, ...}
  virtual const EnumerationTypes::RGROUPS &next() {
    if (m_numPermutationsProcessed) {
      increment();
    } else
      ++m_numPermutationsProcessed;

    return m_permutation;
  }

  virtual boost::uint64_t getPermutationIdx() const {
    return m_numPermutationsProcessed;
  }

  virtual operator bool() const { return hasNext(); }

  EnumerationStrategyBase *copy() const {
    return new CartesianProductStrategy(*this);
  }

 private:
  void increment() {
    next(0);
    ++m_numPermutationsProcessed;
  }

  bool hasNext() const {
    // Fix me -> use multiprecision int here???
    if (m_numPermutations == EnumerationStrategyBase::EnumerationOverflow ||
        m_numPermutationsProcessed < rdcast<size_t>(m_numPermutations)) {
      return true;
    } else {
      return false;
    }
  }

  void next(size_t rowToIncrement) {
    if (!hasNext()) return;
    m_permutation[rowToIncrement] += 1;
    size_t max_index_of_row = m_permutationSizes[rowToIncrement] - 1;
    if (m_permutation[rowToIncrement] > max_index_of_row) {
      m_permutation[rowToIncrement] = 0;
      next(rowToIncrement + 1);
    }
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
    ar &boost::serialization::base_object<EnumerationStrategyBase>(*this);
    ar &m_numPermutationsProcessed;
  }
#endif
};
}  // namespace RDKit

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::CartesianProductStrategy, 1)
#endif

#endif
