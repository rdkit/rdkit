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

#include "Enumerate.h"
#include "CartesianProduct.h"
#include "RandomSample.h"
#include "RandomSampleAllBBs.h"
#include "EvenSamplePairs.h"
#include "../ReactionPickler.h"
#include <GraphMol/MolPickler.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include "EnumerationPickler.h"

// Since we are exporting the classes for serialization,
//  we should declare the archives types used here
#include <RDGeneral/BoostStartInclude.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>
#include <RDGeneral/BoostEndInclude.h>

BOOST_CLASS_EXPORT(RDKit::EnumerationStrategyBase);
BOOST_CLASS_EXPORT(RDKit::CartesianProductStrategy);
BOOST_CLASS_EXPORT(RDKit::RandomSampleStrategy);
BOOST_CLASS_EXPORT(RDKit::RandomSampleAllBBsStrategy);
BOOST_CLASS_EXPORT(RDKit::EvenSamplePairsStrategy);
BOOST_CLASS_EXPORT(RDKit::EnumerateLibrary);

namespace RDKit {
const RGROUPS &EnumerateLibraryBase::getPosition() const {
  return m_enumerator->currentPosition();
}

std::string EnumerateLibraryBase::getState() const {
  PRECONDITION(m_enumerator.get(), "Null Enumerator");
  std::string state;
  EnumerationStrategyPickler::pickle(m_enumerator, state);
  return state;
}

void EnumerateLibraryBase::setState(const std::string &state) {
  m_enumerator = EnumerationStrategyPickler::fromPickle(state);
}

std::vector<std::vector<std::string> > EnumerateLibraryBase::nextSmiles() {
  std::vector<std::vector<std::string> > result;
  std::vector<MOL_SPTR_VECT> mols = next();
  const bool doisomeric = true;
  result.resize(mols.size());
  for (size_t i = 0; i < mols.size(); ++i) {
    result[i].resize(mols[i].size());
    for (size_t j = 0; j < mols[i].size(); ++j) {
      if (mols[i][j].get()) result[i][j] = MolToSmiles(*mols[i][j], doisomeric);
    }
  }
  return result;
}

EnumerateLibrary::EnumerateLibrary(const ChemicalReaction &rxn, const BBS &bbs)
    : EnumerateLibraryBase(rxn, new CartesianProductStrategy), m_bbs(bbs) {
  m_enumerator->initialize(rxn, bbs);  // getSizesFromBBs(bbs));
  PRECONDITION(static_cast<bool>(*m_enumerator), "Nothing to enumerate.");
}

EnumerateLibrary::EnumerateLibrary(const ChemicalReaction &rxn, const BBS &bbs,
                                   const EnumerationStrategyBase &enumerator)
    : EnumerateLibraryBase(rxn), m_bbs(bbs) {
  m_enumerator.reset(enumerator.Clone());
  m_enumerator->initialize(rxn, bbs);
  PRECONDITION(static_cast<bool>(*m_enumerator), "dkjfdkf");
}

EnumerateLibrary::EnumerateLibrary(const EnumerateLibrary &rhs)
    : EnumerateLibraryBase(rhs), m_bbs(rhs.m_bbs) {}

void EnumerateLibrary::reset() {
    m_enumerator->initialize(m_rxn, m_bbs);
}

std::vector<MOL_SPTR_VECT> EnumerateLibrary::next() {
  PRECONDITION(static_cast<bool>(*this), "No more enumerations");
  const RGROUPS &reactantIndices = m_enumerator->next();
  MOL_SPTR_VECT reactants(m_bbs.size());

  for (size_t i = 0; i < m_bbs.size(); ++i) {
    reactants[i] = m_bbs[i][reactantIndices[i]];
  }

  return m_rxn.runReactants(reactants);
}

void EnumerateLibrary::toStream(std::ostream &ss) const {
  boost::archive::binary_oarchive ar(ss);
  ar << *this;
}

void EnumerateLibrary::initFromStream(std::istream &ss) {
  boost::archive::binary_iarchive ar(ss);
  ar >> *this;
}

ssize_t computeNumProducts(const RGROUPS &sizes) {
  boost::multiprecision::cpp_int myint = 1;

  for (size_t i = 0; i < sizes.size(); ++i) {
    myint *= sizes[i];
  }

  if (myint < std::numeric_limits<ssize_t>::max())
    return myint.convert_to<ssize_t>();
  else
    return EnumerationStrategyBase::EnumerationOverflow;
}

MOL_SPTR_VECT getReactantsFromRGroups(const std::vector<MOL_SPTR_VECT> &bbs,
                                      const RGROUPS &rgroups) {
  PRECONDITION(bbs.size() == rgroups.size(),
               "BBS and RGROUPS must have the same # reactants");
  MOL_SPTR_VECT result;
  result.reserve(bbs.size());
  for (size_t i = 0; i < bbs.size(); ++i) {
    result.push_back(bbs[i][rgroups[i]]);
  }
  return result;
}
}
