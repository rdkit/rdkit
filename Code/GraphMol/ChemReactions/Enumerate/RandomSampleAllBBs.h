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
#ifndef RGROUP_RANDOM_SAMPLE_ALLBBS_H
#define RGROUP_RANDOM_SAMPLE_ALLBBS_H

#include "EnumerationStrategyBase.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <sstream>

namespace RDKit {
//! RandomSampleAllBBsStrategy
//!  Randomly sample rgroup indices

//! This is a class for randomly enumerating reagents that ensures all reagents
//  are sampled.
/*!
  basic usage:

  \verbatim
  std::vector<MOL_SPTR_VECT> bbs;
  bbs.push_back( bbs_for_reactants_1 );
  bbs.push_back( bbs_for_reactants_2 );

  RandomSampleAllBBsStrategy rgroups;
  rgroups.initialize(rxn, bbs);
  for(size_t i=0; i<num_samples && rgroups; ++i) {
    MOL_SPTR_VECT rvect = getReactantsFromRGroups(bbs, rgroups.next());
    std::vector<MOL_SPTR_VECT> lprops = rxn.RunReactants(rvect);
    ...
  }
  \endverbatim

  See EnumerationStrategyBase for more details and usage.
*/

class RDKIT_CHEMREACTIONS_EXPORT RandomSampleAllBBsStrategy
    : public EnumerationStrategyBase {
  boost::uint64_t m_numPermutationsProcessed{0};
  size_t m_offset{0};
  size_t m_maxoffset{0};

  boost::minstd_rand m_rng;
  std::vector<boost::random::uniform_int_distribution<>> m_distributions;

 public:
  RandomSampleAllBBsStrategy()
      : EnumerationStrategyBase(),
        
        m_rng(),
        m_distributions() {
    for (size_t i = 0; i < m_permutation.size(); ++i) {
      m_distributions.emplace_back(0, m_permutation[i] - 1);
    }
  }
  using EnumerationStrategyBase::initialize;

  void initializeStrategy(const ChemicalReaction &,
                          const EnumerationTypes::BBS &) {
    m_distributions.clear();
    m_permutation.resize(m_permutationSizes.size());
    m_offset = 0;
    m_maxoffset =
        *std::max_element(m_permutationSizes.begin(), m_permutationSizes.end());
    for (size_t i = 0; i < m_permutationSizes.size(); ++i) {
      m_distributions.emplace_back(
          0, m_permutationSizes[i] - 1);
    }

    m_numPermutationsProcessed = 0;
  }

  virtual const char *type() const { return "RandomSampleAllBBsStrategy"; }

  //! The current permutation {r1, r2, ...}
  virtual const EnumerationTypes::RGROUPS &next() {
    if (m_offset >= m_maxoffset) {
      for (size_t i = 0; i < m_permutation.size(); ++i) {
        m_permutation[i] = m_distributions[i](m_rng);
      }
      m_offset = 0;
    } else {
      for (size_t i = 0; i < m_permutation.size(); ++i) {
        m_permutation[i] = (m_permutation[i] + 1) % m_permutationSizes[i];
      }
      ++m_offset;
    }
    ++m_numPermutationsProcessed;

    return m_permutation;
  }

  virtual boost::uint64_t getPermutationIdx() const {
    return m_numPermutationsProcessed;
  }

  virtual operator bool() const { return true; }

  EnumerationStrategyBase *copy() const {
    return new RandomSampleAllBBsStrategy(*this);
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;

  template <class Archive>
  void save(Archive &ar, const unsigned int /*version*/) const {
    // invoke serialization of the base class
    ar << boost::serialization::base_object<const EnumerationStrategyBase>(
        *this);
    ar << m_numPermutationsProcessed;

    std::stringstream random;
    random << m_rng;
    std::string s = random.str();
    ar << s;

    ar << m_offset;
    ar << m_maxoffset;
  }

  template <class Archive>
  void load(Archive &ar, const unsigned int /*version*/) {
    // invoke serialization of the base class
    ar >> boost::serialization::base_object<EnumerationStrategyBase>(*this);
    ar >> m_numPermutationsProcessed;
    std::string s;
    ar >> s;
    std::stringstream random(s);
    random >> m_rng;
    ar >> m_offset;
    ar >> m_maxoffset;

    // reset the uniform distributions
    m_distributions.clear();
    for (size_t i = 0; i < m_permutationSizes.size(); ++i) {
      m_distributions.emplace_back(
          0, m_permutationSizes[i] - 1);
    }
  }

  template <class Archive>
  void serialize(Archive &ar, const unsigned int file_version) {
    boost::serialization::split_member(ar, *this, file_version);
  }
#endif
};
}  // namespace RDKit

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::RandomSampleAllBBsStrategy, 1)
#endif

#endif
