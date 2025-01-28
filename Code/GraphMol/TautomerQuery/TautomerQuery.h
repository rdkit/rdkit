//
// Created by Gareth Jones on 5/7/2020.
//
// Copyright 2020-2022 Schrodinger, Inc and other RDKit contributors
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <RDGeneral/export.h>

#ifndef RDKIT_TAUTOMERQUERY_H
#define RDKIT_TAUTOMERQUERY_H

#include <GraphMol/ROMol.h>
#include <GraphMol/MolPickler.h>
#include <vector>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>

#ifdef RDK_USE_BOOST_SERIALIZATION
#include <RDGeneral/BoostStartInclude.h>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/split_member.hpp>
#include <RDGeneral/BoostEndInclude.h>
#endif

namespace RDKit {

class RWMol;

RDKIT_TAUTOMERQUERY_EXPORT bool TautomerQueryCanSerialize();

class RDKIT_TAUTOMERQUERY_EXPORT TautomerQuery {
 private:
  // Tautomers of the query
  std::vector<ROMOL_SPTR> d_tautomers;
  // Template query for substructure search
  std::unique_ptr<const ROMol> d_templateMolecule;
  // Tautomeric bonds and atoms
  std::vector<size_t> d_modifiedAtoms;
  std::vector<size_t> d_modifiedBonds;

  // tests if a match to the template matches a specific tautomer
  bool matchTautomer(const ROMol &mol, const ROMol &tautomer,
                     const std::vector<unsigned int> &match,
                     const SubstructMatchParameters &params) const;

 public:
  TautomerQuery(std::vector<ROMOL_SPTR> tautomers,
                const ROMol *const templateMolecule,
                std::vector<size_t> modifiedAtoms,
                std::vector<size_t> modifiedBonds);

  //! Copy constructor performs a deep copy
  TautomerQuery(const TautomerQuery &other)
      : d_templateMolecule(other.d_templateMolecule
                               ? new ROMol(*other.d_templateMolecule)
                               : nullptr),
        d_modifiedAtoms(other.d_modifiedAtoms),
        d_modifiedBonds(other.d_modifiedBonds) {
    PRECONDITION(other.d_templateMolecule != nullptr, "Null template");
    for (auto taut : other.d_tautomers) {
      PRECONDITION(taut.get() != nullptr, "Null tautomer");
      d_tautomers.push_back(boost::make_shared<ROMol>(*taut));
    }
  }

  TautomerQuery(const std::string &pickle) { initFromString(pickle); }

  // Factory to build TautomerQuery
  // Caller owns the memory
  static TautomerQuery *fromMol(
      const ROMol &molecule,
      const std::string &tautomerTransformFile = std::string());

  // Substructure search
  std::vector<MatchVectType> substructOf(
      const ROMol &mol,
      const SubstructMatchParameters &params = SubstructMatchParameters(),
      std::vector<ROMOL_SPTR> *matchingTautomers = nullptr) const;

  // SubstructureMatch
  bool isSubstructOf(const ROMol &mol, const SubstructMatchParameters &params =
                                           SubstructMatchParameters());

  // Query fingerprint
  ExplicitBitVect *patternFingerprintTemplate(
      unsigned int fpSize = 2048U) const;
  // Static method to Fingerprint a target
  static ExplicitBitVect *patternFingerprintTarget(const ROMol &target,
                                                   unsigned int fpSize = 2048U);

  // accessors

  // pointer is owned by TautomerQuery
  const ROMol &getTemplateMolecule() const { return *d_templateMolecule; }

  const std::vector<ROMOL_SPTR> getTautomers() const { return d_tautomers; }

  const std::vector<size_t> getModifiedAtoms() const { return d_modifiedAtoms; }

  const std::vector<size_t> getModifiedBonds() const { return d_modifiedBonds; }

  //! serializes (pickles) to a stream
  void toStream(std::ostream &ss) const;
  //! returns a string with a serialized (pickled) representation
  std::string serialize() const;
  //! initializes from a stream pickle
  void initFromStream(std::istream &ss);
  //! initializes from a string pickle
  void initFromString(const std::string &text);

  friend class TautomerQueryMatcher;

#ifdef RDK_USE_BOOST_SERIALIZATION
  template <class Archive>
  void save(Archive &ar, const unsigned int version) const {
    RDUNUSED_PARAM(version);
    std::vector<std::string> pkls;
    for (const auto &taut : d_tautomers) {
      std::string pkl;
      MolPickler::pickleMol(*taut, pkl, PicklerOps::AllProps);
      pkls.push_back(pkl);
    }
    ar << pkls;
    std::string molpkl;
    MolPickler::pickleMol(*d_templateMolecule, molpkl, PicklerOps::AllProps);
    ar << molpkl;
    ar << d_modifiedAtoms;
    ar << d_modifiedBonds;
  }

  template <class Archive>
  void load(Archive &ar, const unsigned int version) {
    RDUNUSED_PARAM(version);

    std::vector<std::string> pkls;
    ar >> pkls;
    d_tautomers.clear();
    for (const auto &pkl : pkls) {
      d_tautomers.push_back(ROMOL_SPTR(new ROMol(pkl)));
    }
    std::string molpkl;
    ar >> molpkl;
    d_templateMolecule.reset(new ROMol(molpkl));

    ar >> d_modifiedAtoms;
    ar >> d_modifiedBonds;
  }
  BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif
};

// so we can use the templates in Code/GraphMol/Substruct/SubstructMatch.h
RDKIT_TAUTOMERQUERY_EXPORT std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const TautomerQuery &query,
    const SubstructMatchParameters &params);

}  // namespace RDKit

#endif  // RDKIT_TAUTOMERQUERY_H
