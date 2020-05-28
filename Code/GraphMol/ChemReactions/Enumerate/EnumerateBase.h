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
#ifndef RDKIT_ENUMERATEBASE_H
#define RDKIT_ENUMERATEBASE_H

#include <vector>
#include "EnumerateTypes.h"
#include "../Reaction.h"
#include "EnumerationPickler.h"

#include "EnumerationStrategyBase.h"
#include "CartesianProduct.h"
#include "../ReactionPickler.h"
#include <GraphMol/MolPickler.h>

namespace RDKit {
//! Base class for enumerating chemical reactions from collections of
//  building blocks and reagents.
/*!
  basic usage:

  \verbatim
  EnumerateLibraryBase &enumerator;
  while (enumerator) {
    MOL_SPTR_VECT res = enumerator.next();
    // do something with enumeration products here
  }
  \endverbatim

  See Reaction.h for more details on how ChemicalReactions are
  used.
*/
class RDKIT_CHEMREACTIONS_EXPORT EnumerateLibraryBase {
 protected:
  ChemicalReaction m_rxn;
  boost::shared_ptr<EnumerationStrategyBase> m_enumerator;
  boost::shared_ptr<EnumerationStrategyBase> m_initialEnumerator;

 public:
  //! default constructor
  EnumerateLibraryBase() : m_rxn(), m_enumerator(), m_initialEnumerator() {}

  //! construct with a chemical reaction and an enumeration strategy
  EnumerateLibraryBase(const ChemicalReaction &rxn,
                       EnumerationStrategyBase *enumerator = nullptr)
      : m_rxn(rxn),
        m_enumerator(enumerator ? enumerator : new CartesianProductStrategy),
        m_initialEnumerator(m_enumerator->copy()) {
    m_rxn.initReactantMatchers();
  }

  //! Copy constructor
  EnumerateLibraryBase(const EnumerateLibraryBase &rhs)
      : m_rxn(rhs.m_rxn),
        m_enumerator(rhs.m_enumerator ? rhs.m_enumerator->copy() : nullptr),
        m_initialEnumerator(m_enumerator->copy()) {}

  virtual ~EnumerateLibraryBase() {}

  //! Are there any enumerations left?
  virtual operator bool() const {
    PRECONDITION(m_enumerator.get(), "Null enumeration strategy");
    return static_cast<bool>(*m_enumerator);
  }

  //! reset the enumeration to the beginning.
  void reset() {
    if (m_initialEnumerator.get()) {
      m_enumerator.reset(m_initialEnumerator->copy());
    }
  }

  //! returns the underlying chemical reaction
  const ChemicalReaction &getReaction() const { return m_rxn; }

  //! return the current enumeration strategy
  const EnumerationStrategyBase &getEnumerator() {
    PRECONDITION(m_enumerator.get(), "Null Enumerator");
    return *m_enumerator;
  }

  //! get the next set of products (See run_Reactants) for details
  //  This returns a vector of a vector of molecules.
  //  Each result vector corresponds for a product template.
  //  i.e.
  //    res = library.next();
  //    res[0] are the results for library.getReaction().getProdcts()[0]
  virtual std::vector<MOL_SPTR_VECT> next() = 0;

  //! get the next set of products as smiles
  //  This returns a vector of a vector strings.
  //  Each result vector corresponds for a product template.
  virtual std::vector<std::vector<std::string>> nextSmiles();

  //! Get the current position into the reagent vectors
  //   Use getState/setState to save/restart the enumeration
  //   from this position.
  const EnumerationTypes::RGROUPS &getPosition() const;

  //! Get the current state of the enumerator
  //   This is the position of the enumerator and the enumerators
  //   state that can be used to restart enumerating
  //   from this position.
  std::string getState() const;

  //! Set the current state of the enumerator
  //   Restart the enumerator from this position.
  void setState(const std::string &);

  //! Reset the enumerator to the beginning
  void resetState();

  //! serializes (pickles) to a stream
  virtual void toStream(std::ostream &ss) const = 0;

  //! returns a string with a serialized (pickled) representation
  virtual std::string Serialize() const {
    std::stringstream ss;
    toStream(ss);
    return ss.str();
  }

  //! initializes from a stream pickle
  virtual void initFromStream(std::istream &ss) = 0;

  //! initializes from a string pickle
  virtual void initFromString(const std::string &text) {
    std::stringstream ss(text);
    initFromStream(ss);
  }

 private:
#ifdef RDK_USE_BOOST_SERIALIZATION
  friend class boost::serialization::access;
  template <class Archive>
  void save(Archive &ar, const unsigned int) const {
    std::string pickle;
    ReactionPickler::pickleReaction(m_rxn, pickle);
    ar &pickle;
    ar &m_enumerator;
    // we handle the m_initialEnumerator from a string
    //  for backwards compatibility with a unreleased
    //  version
    EnumerationStrategyPickler::pickle(m_initialEnumerator, pickle);
    ar &pickle;
  }
  template <class Archive>
  void load(Archive &ar, const unsigned int /*version*/) {
    // this should only be called on non-initialized reactions
    if (m_rxn.getNumReactantTemplates() || m_rxn.getNumProductTemplates() ||
        m_rxn.getNumAgentTemplates()) {
      throw ValueErrorException("EnumerateBase already created from archive.");
    }
    std::string pickle;
    ar &pickle;
    ReactionPickler::reactionFromPickle(pickle, m_rxn);
    ar &m_enumerator;
    ar &pickle;
    m_initialEnumerator = EnumerationStrategyPickler::fromPickle(pickle);
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER();
#endif
};

#ifdef RDK_USE_BOOST_SERIALIZATION
BOOST_SERIALIZATION_ASSUME_ABSTRACT(EnumerateLibraryBase)
#endif
}  // namespace RDKit
#endif
