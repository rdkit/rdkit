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

#include <RDGeneral/BoostStartInclude.h>
#include <boost/multiprecision/cpp_int.hpp>
#ifdef RDK_USE_BOOST_SERIALIZATION
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/export.hpp>
#endif
#include <RDGeneral/BoostEndInclude.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
// Since we are exporting the classes for serialization,
//  we should declare the archives types used here
BOOST_CLASS_EXPORT(RDKit::EnumerationStrategyBase);
BOOST_CLASS_EXPORT(RDKit::CartesianProductStrategy);
BOOST_CLASS_EXPORT(RDKit::RandomSampleStrategy);
BOOST_CLASS_EXPORT(RDKit::RandomSampleAllBBsStrategy);
BOOST_CLASS_EXPORT(RDKit::EvenSamplePairsStrategy);
BOOST_CLASS_EXPORT(RDKit::EnumerateLibrary);
#endif

namespace RDKit {
using namespace EnumerationTypes;

const RGROUPS &EnumerateLibraryBase::getPosition() const {
  return m_enumerator->getPosition();
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

void EnumerateLibraryBase::resetState() {
  PRECONDITION(m_initialEnumerator.get(), "Unset initial enumerator");
  m_enumerator.reset(m_initialEnumerator->copy());
}

std::vector<std::vector<std::string>> EnumerateLibraryBase::nextSmiles() {
  std::vector<std::vector<std::string>> result;
  std::vector<MOL_SPTR_VECT> mols = next();
  const bool doisomeric = true;
  result.resize(mols.size());
  for (size_t i = 0; i < mols.size(); ++i) {
    result[i].resize(mols[i].size());
    for (size_t j = 0; j < mols[i].size(); ++j) {
      if (mols[i][j].get()) {
        result[i][j] = MolToSmiles(*mols[i][j], doisomeric);
      }
    }
  }
  return result;
}

namespace {
size_t countMatches(const ROMol &bb, const ROMol &query, int maxMatches) {
  std::vector<MatchVectType> matches;
  const bool uniquify = true;
  const bool useChirality = true;
  const bool useQueryQueryMatches = false;

  SubstructMatch(bb, query, matches, uniquify, true, useChirality,
                 useQueryQueryMatches, maxMatches + 1);
  return matches.size();
}
}  // namespace
BBS removeNonmatchingReagents(const ChemicalReaction &rxn, BBS bbs,
                              const EnumerationParams &params) {
  PRECONDITION(bbs.size() <= rxn.getNumReactantTemplates(),
               "Number of Reagents not compatible with reaction templates");
  BBS result;
  result.resize(bbs.size());

  for (size_t reactant_idx = 0; reactant_idx < bbs.size(); ++reactant_idx) {
    size_t removedCount = 0;
    const unsigned int maxMatches =
        (params.reagentMaxMatchCount == INT_MAX)
            ? 0
            : rdcast<unsigned int>(params.reagentMaxMatchCount);

    ROMOL_SPTR reactantTemplate = rxn.getReactants()[reactant_idx];
    for (size_t reagent_idx = 0; reagent_idx < bbs[reactant_idx].size();
         ++reagent_idx) {
      ROMOL_SPTR mol = bbs[reactant_idx][reagent_idx];
      size_t matches =
          countMatches(*mol.get(), *reactantTemplate.get(), maxMatches);

      bool removeReagent = false;
      if (!matches || matches > rdcast<size_t>(params.reagentMaxMatchCount)) {
        removeReagent = true;
      }

      if (!removeReagent && params.sanePartialProducts) {
        // see if we have any sane products in the results
        std::vector<MOL_SPTR_VECT> partialProducts =
            rxn.runReactant(mol, reactant_idx);
        for (auto &partialProduct : partialProducts) {
          int saneProducts = 0;
          for (auto &product_idx : partialProduct) {
            try {
              auto *m = dynamic_cast<RWMol *>(product_idx.get());
              MolOps::sanitizeMol(*m);
              saneProducts++;
            } catch (...) {
            }
          }

          if (!saneProducts) {
            // if any product template has no sane products, we bail
            removeReagent = true;
            break;
          }
        }
      }

      if (removeReagent) {
        removedCount++;
      } else {
        result[reactant_idx].push_back(mol);
      }
    }

    if (removedCount) {
      BOOST_LOG(rdInfoLog) << "Removed " << removedCount
                           << " non matching reagents at template "
                           << reactant_idx << std::endl;
    }
  }
  return result;
}

EnumerateLibrary::EnumerateLibrary(const ChemicalReaction &rxn, const BBS &bbs,
                                   const EnumerationParams &params)
    : EnumerateLibraryBase(rxn, new CartesianProductStrategy),
      m_bbs(removeNonmatchingReagents(m_rxn, bbs, params)) {
  m_enumerator->initialize(m_rxn, m_bbs);  // getSizesFromBBs(bbs));
  m_initialEnumerator.reset(m_enumerator->copy());
}

EnumerateLibrary::EnumerateLibrary(const ChemicalReaction &rxn, const BBS &bbs,
                                   const EnumerationStrategyBase &enumerator,
                                   const EnumerationParams &params)
    : EnumerateLibraryBase(rxn),
      m_bbs(removeNonmatchingReagents(m_rxn, bbs, params)) {
  m_enumerator.reset(enumerator.copy());
  m_enumerator->initialize(m_rxn, m_bbs);
  m_initialEnumerator.reset(m_enumerator->copy());
}

EnumerateLibrary::EnumerateLibrary(const EnumerateLibrary &rhs)
    : EnumerateLibraryBase(rhs), m_bbs(rhs.m_bbs) {}

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
#ifdef RDK_USE_BOOST_SERIALIZATION
  boost::archive::text_oarchive ar(ss);
  ar << *this;
#else
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}

void EnumerateLibrary::initFromStream(std::istream &ss) {
#ifdef RDK_USE_BOOST_SERIALIZATION
  boost::archive::text_iarchive ar(ss);
  ar >> *this;
#else
  PRECONDITION(0, "BOOST SERIALIZATION NOT INSTALLED");
#endif
}

boost::uint64_t computeNumProducts(const RGROUPS &sizes) {
#ifdef RDK_HAVE_MULTIPREC
  boost::multiprecision::cpp_int myint = 1;

  for (boost::uint64_t size : sizes) {
    myint *= size;
  }

  if (myint < std::numeric_limits<boost::uint64_t>::max()) {
    return myint.convert_to<boost::uint64_t>();
  } else {
    return EnumerationStrategyBase::EnumerationOverflow;
  }
#else
  boost::uint64_t myint = 1;

  for (size_t i = 0; i < sizes.size(); ++i) {
    if (sizes[i] &&
        (std::numeric_limits<boost::uint64_t>::max() / sizes[i]) < myint) {
      return EnumerationStrategyBase::EnumerationOverflow;
    }
    myint *= sizes[i];
  }
  return myint;
#endif
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

bool EnumerateLibraryCanSerialize() {
#ifdef RDK_USE_BOOST_SERIALIZATION
  return true;
#else
  return false;
#endif
}
}  // namespace RDKit
