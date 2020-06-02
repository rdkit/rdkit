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
#ifndef ENUMERATION_STRATEGY_H
#define ENUMERATION_STRATEGY_H

#include "EnumerateTypes.h"
#include "../Reaction.h"
#include <vector>
#include <RDGeneral/BoostStartInclude.h>
#include <cstdint>
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/vector.hpp>
// the next two includes need to be there for boost 1.56
#include <boost/serialization/singleton.hpp>
#include <boost/serialization/extended_type_info.hpp>
#include <boost/serialization/shared_ptr.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

#include <GraphMol/RDKitBase.h>

namespace RDKit {

//! class for flagging enumeration strategy errors
class RDKIT_CHEMREACTIONS_EXPORT EnumerationStrategyException
    : public std::exception {
 public:
  EnumerationStrategyException(const char *msg) : _msg(msg){};
  EnumerationStrategyException(const std::string &msg) : _msg(msg){};
  const char *what() const noexcept override { return _msg.c_str(); };
  ~EnumerationStrategyException() noexcept {};

 private:
  std::string _msg;
};

//! Return the number of elements per input vector
/*!  \param bbs vector<vector<T> >

  \result vector<unint64_t> number of elements in each vector
 */
template <class T>
EnumerationTypes::RGROUPS getSizesFromBBs(
    const std::vector<std::vector<T>> &bbs) {
  EnumerationTypes::RGROUPS sizes;
  for (size_t i = 0; i < bbs.size(); ++i) sizes.push_back(bbs[i].size());
  return sizes;
}

//! getSizesFromReactants
//!  Helper function for enumeration, bbs are stored in a
//!   std::vector< std::vector<boost:shared_ptr<ROMol> >
//
RDKIT_CHEMREACTIONS_EXPORT EnumerationTypes::RGROUPS getSizesFromReactants(
    const std::vector<MOL_SPTR_VECT> &bbs);

//! getReactantsFromRGroups
//!  Helper function for enumeration, bbs are stored in a
//!   std::vector< std::vector<boost:shared_ptr<ROMol> >
//
RDKIT_CHEMREACTIONS_EXPORT MOL_SPTR_VECT
getReactantsFromRGroups(const std::vector<MOL_SPTR_VECT> &bbs,
                        const EnumerationTypes::RGROUPS &rgroups);

//! computeNumProducts
//!  Returns the number of possible product combination from
//!   The given numbers of building blocks for each rgroup
//!   or EnumerationStrategyBase::EnumerationOverflow if the
//!   number will not fit into the machines integer type.
//!   n.b. An overflow simply means there are a lot of products
//!     not that they cannot be enumerated
RDKIT_CHEMREACTIONS_EXPORT boost::uint64_t computeNumProducts(
    const EnumerationTypes::RGROUPS &sizes);

//! Base Class for enumeration strategies
//!  Usage:
//!  EnumerationStrategyBase must be initialized with both a reaction
//!   and the building block (molecule) vector to be sampled.
//!
//!  \verbatim
//!  EnumerationStrategyBase &eb = ...
//!   if(eb) { // can we get another entry
//!    const std::vector<int> &v = eb.next();
//!    v[0] // RGroup 0 position
//!    v[1] // RGroup 1 position...
//!   }
//!  \endverbatim

class RDKIT_CHEMREACTIONS_EXPORT EnumerationStrategyBase {
 protected:
  EnumerationTypes::RGROUPS m_permutation;  // where are we currently?
  EnumerationTypes::RGROUPS
      m_permutationSizes;  // m_permutationSizes num bbs per group
  boost::uint64_t
      m_numPermutations{};  // total number of permutations for this group
                          //  -1 if > ssize_t::max
 public:
  static const boost::uint64_t EnumerationOverflow =
      static_cast<boost::uint64_t>(-1);
  EnumerationStrategyBase()
      : m_permutation(), m_permutationSizes() {}

  virtual ~EnumerationStrategyBase() {}

  virtual const char *type() const { return "EnumerationStrategyBase"; }

  //! Initialize the enumerator based on the reaction and the
  //! supplied building blocks
  //!  This is the standard API point.
  //!  This calls the derived class's initializeStrategy method which must be implemented
  void initialize(const ChemicalReaction &reaction,
                  const EnumerationTypes::BBS &building_blocks) {
    // default initialization, may be overridden (sets the # reactants
    //  and computes the default # of permutations)
    m_permutationSizes = getSizesFromBBs(building_blocks);
    m_permutation.resize(m_permutationSizes.size());

    m_numPermutations = computeNumProducts(m_permutationSizes);
    std::fill(m_permutation.begin(), m_permutation.end(), 0);

    initializeStrategy(reaction, building_blocks);
  }

  // ! Initialize derived class. Must exist.
  // ! EnumerationStrategyBase structures are already initialized:
  // !  m_permutationSizes - [ length of building blocks for each reactant set ]
  // !  m_numPermutations - number of possible permutations ( -1 if not computable )
  // !  m_permutation - the first permutation, always the first supplied reactants
  virtual void initializeStrategy(
      const ChemicalReaction &reaction,
      const EnumerationTypes::BBS &building_blocks) = 0;

  //! returns true if there are more permutations left
  //!  random enumerators may always return true...
  virtual operator bool() const = 0;

  //! The current permutation {r1, r2, ...}
  virtual const EnumerationTypes::RGROUPS &next() = 0;

  //! copy the enumeration strategy complete with current state
  virtual EnumerationStrategyBase *copy() const = 0;

  //! The current position in the enumeration
  const EnumerationTypes::RGROUPS &getPosition() const { return m_permutation; }

  //! a result of EnumerationOverflow indicates that the number of
  //!  permutations is not computable with the current
  //!  rdlonglong size.
  boost::uint64_t getNumPermutations() const { return m_numPermutations; }

  //! Returns how many permutations have been processed by this strategy
  virtual boost::uint64_t getPermutationIdx() const = 0;

  //! Skip the specified number of permutations (useful for
  //!  resetting state to a known position)
  bool skip(boost::uint64_t skipCount) {
    for (boost::uint64_t i = 0; i < skipCount; ++i) next();
    return true;
  }

 protected:
  //! Initialize the internal data structures
  //!  i.e. RGROUPS = {10,40,50};
  void internalInitialize(const EnumerationTypes::RGROUPS &rgroups) {
    m_permutation.resize(rgroups.size());
    m_permutationSizes = rgroups;
    m_numPermutations = computeNumProducts(m_permutationSizes);
    std::fill(m_permutation.begin(), m_permutation.end(), 0);
  }

 private:
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /*version*/) {
    ar &m_permutation;
    ar &m_permutationSizes;
    ar &m_numPermutations;
  }
};
#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_SERIALIZATION_ASSUME_ABSTRACT(EnumerationStrategyBase)
#endif
}  // namespace RDKit

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_CLASS_VERSION(RDKit::EnumerationStrategyBase, 1)
#endif

#endif
