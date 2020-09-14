//
// Created by Gareth Jones on 5/7/2020.
//
// Copyright 2020 Schrodinger, Inc
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.

#include <RDGeneral/export.h>

#ifndef RDKIT_TAUTOMERQUERY_H
#define RDKIT_TAUTOMERQUERY_H

#include <GraphMol/ROMol.h>
#include <vector>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>


namespace RDKit {

class RWMol;

class RDKIT_TAUTOMERQUERY_EXPORT TautomerQuery {
 private:
  // Tautomers of the query
  const std::vector<ROMOL_SPTR> d_tautomers;
  // Template query for substructure search
  const ROMol *const d_templateMolecule;
  // Tautomeric bonds and atoms
  const std::vector<size_t> d_modifiedAtoms;
  const std::vector<size_t> d_modifiedBonds;

  TautomerQuery(const std::vector<ROMOL_SPTR> &tautomers,
                const ROMol *const templateMolecule,
                const std::vector<size_t> &modifiedAtoms,
                const std::vector<size_t> &modifiedBonds);

  // tests if a match to the template matches a specific tautomer
  bool matchTautomer(const ROMol &mol, const ROMol &tautomer,
                     const std::vector<unsigned int> &match,
                     const SubstructMatchParameters &params) const;

 public:
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
  ExplicitBitVect *patternFingerprintTemplate(unsigned int fpSize = 2048U);

  // Static method to Fingerprint a target
  static ExplicitBitVect *patternFingerprintTarget(const ROMol &target,
                                                   unsigned int fpSize = 2048U);

  // accessors

  // pointer is owned by TautomerQuery
  const ROMol & getTemplateMolecule() const { return *d_templateMolecule; }

  const std::vector<ROMOL_SPTR> getTautomers() const { return d_tautomers; }

  const std::vector<size_t> getModifiedAtoms() const { return d_modifiedAtoms; }

  const std::vector<size_t> getModifiedBonds() const { return d_modifiedBonds; }

  ~TautomerQuery();

  friend class TautomerQueryMatcher;
};

// so we can use the templates in Code/GraphMol/Substruct/SubstructMatch.h
RDKIT_TAUTOMERQUERY_EXPORT std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const TautomerQuery &query,
    const SubstructMatchParameters &params);

}  // namespace RDKit

#endif  // RDKIT_TAUTOMERQUERY_H
