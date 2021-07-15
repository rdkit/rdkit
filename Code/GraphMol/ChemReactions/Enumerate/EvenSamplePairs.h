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

#include <RDGeneral/export.h>
#ifndef RGROUP_EVEN_SAMPLE_H
#define RGROUP_EVEN_SAMPLE_H

#include "EnumerationStrategyBase.h"
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/serialization/set.hpp>
#endif
#include <cstdint>

namespace RDKit {
//! EvenSamplePairsStrategy
/*!  Randomly sample Pairs evenly from a collection of building blocks
      This is a good strategy for choosing a relatively small selection
      of building blocks from a larger set.  As the amount of work needed
      to retrieve the next evenly sample building block grows with the
      number of samples, this method performs progressively worse as the
      number of samples gets larger.

      See EnumeartionStrategyBase for more details.
*/

class RDKIT_CHEMREACTIONS_EXPORT EvenSamplePairsStrategy
    : public EnumerationStrategyBase {
  boost::uint64_t m_numPermutationsProcessed{};

  std::vector<boost::int64_t> used_count;
  std::vector<std::vector<boost::uint64_t>> var_used;
  std::vector<std::vector<boost::uint64_t>> pair_used;
  std::vector<std::vector<boost::uint64_t>> pair_counts;
  std::set<boost::uint64_t> selected;

  boost::uint64_t seed{};         // last seed for permutation (starts at 0)
  boost::uint64_t M{}, a{}, b{};  // random number stuff
  boost::uint64_t nslack{}, min_nslack{};
  boost::uint64_t rejected_period{}, rejected_unique{};
  boost::uint64_t rejected_slack_condition{}, rejected_bb_sampling_condition{};

 public:
  EvenSamplePairsStrategy()
      : EnumerationStrategyBase(),

        used_count(),
        var_used(),
        pair_used(),
        pair_counts(),
        selected() {}

  EvenSamplePairsStrategy(const EvenSamplePairsStrategy &rhs)
      : EnumerationStrategyBase(rhs),
        m_numPermutationsProcessed(rhs.m_numPermutationsProcessed),
        used_count(rhs.used_count),
        var_used(rhs.var_used),
        pair_used(rhs.pair_used),
        pair_counts(rhs.pair_counts),
        selected(rhs.selected),
        seed(rhs.seed),
        M(rhs.M),
        a(rhs.a),
        b(rhs.b),
        nslack(rhs.nslack),
        min_nslack(rhs.min_nslack),
        rejected_period(rhs.rejected_period),
        rejected_unique(rhs.rejected_unique),
        rejected_slack_condition(rhs.rejected_slack_condition),
        rejected_bb_sampling_condition(rhs.rejected_bb_sampling_condition) {}

  const char *type() const override { return "EvenSamplePairsStrategy"; }

  //! This is a class for enumerating RGroups using Cartesian Products of
  //! reagents.
  /*!
    basic usage:

    \verbatim
    std::vector<MOL_SPTR_VECT> bbs;
    bbs.push_back( bbs_for_reactants_1 );
    bbs.push_back( bbs_for_reactants_2 );

    EvenSamplePairsStrategy rgroups;
    rgroups.initialize(rxn, bbs);
    for(boost::uint64_t i=0; i<num_samples && rgroups; ++i) {
      MOL_SPTR_VECT rvect = getReactantsFromRGroups(bbs, rgroups.next());
      std::vector<MOL_SPTR_VECT> lprops = rxn.RunReactants(rvect);
      ...
    }
    \endverbatim
  */
  using EnumerationStrategyBase::initialize;

  void initializeStrategy(const ChemicalReaction &,
                          const EnumerationTypes::BBS &) override;

  //! The current permutation {r1, r2, ...}
  const EnumerationTypes::RGROUPS &next() override;

  boost::uint64_t getPermutationIdx() const override {
    return m_numPermutationsProcessed;
  }

  operator bool() const override { return true; }

  EnumerationStrategyBase *copy() const override {
    return new EvenSamplePairsStrategy(*this);
  }

  std::string stats() const;

 private:
  friend class boost::serialization::access;

  // decode a packed integer into an RGroup selection
  const EnumerationTypes::RGROUPS &decode(boost::uint64_t seed) {
    for (boost::int64_t j = m_permutationSizes.size() - 1; j >= 0; j--) {
      m_permutation[j] = seed % m_permutationSizes[j];
      seed /= m_permutationSizes[j];
    }
    return m_permutation;
  }

  bool try_add(boost::uint64_t seed);

 public:
#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
    // invoke serialization of the base class
    ar &boost::serialization::base_object<EnumerationStrategyBase>(*this);
    ar &m_numPermutationsProcessed;
    ar &used_count;
    ar &var_used;
    ar &pair_used;
    ar &pair_counts;
    ar &selected;

    ar &seed;

    ar &M;
    ar &a;
    ar &b;

    ar &nslack;
    ar &min_nslack;
    ar &rejected_period;
    ar &rejected_unique;
    ar &rejected_slack_condition;
    ar &rejected_bb_sampling_condition;
  }
#endif
};
}  // namespace RDKit

BOOST_CLASS_VERSION(RDKit::EvenSamplePairsStrategy, 1)

#endif
