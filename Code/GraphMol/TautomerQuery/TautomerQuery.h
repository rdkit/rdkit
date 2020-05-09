//
// Created by Gareth Jones on 5/7/2020.
//
// Copyright 2020 Schrodinger, Inc
//

#include <RDGeneral/export.h>

#ifndef RDKIT_TAUTOMERQUERY_H
#define RDKIT_TAUTOMERQUERY_H

#include <GraphMol/ROMol.h>
#include <vector>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <DataStructs/ExplicitBitVect.h>

namespace RDKit {

class RWMol;

class TautomerQuery {
 private:
  // Original query- not currently used, but will be useful if the query is not
  // tautomeric
  const ROMol &query;
  // Tautomers of the query
  const std::vector<ROMOL_SPTR> tautomers;
  // Template query for substructure search
  const ROMOL_SPTR templateMolecule;
  // Tautomeric bonds and atoms
  const std::vector<size_t> modifiedAtoms;
  const std::vector<size_t> modifiedBonds;

  TautomerQuery(const ROMol &query, const std::vector<ROMOL_SPTR> tautomers,
                const ROMOL_SPTR queryMolecule,
                const std::vector<size_t> modifiedAtoms,
                const std::vector<size_t> modifiedBonds);

  // tests if a match to the template matches a specific conformer
  bool matchTautomer(const ROMol &mol, const ROMol &tautomer,
                     const MatchVectType &match,
                     const SubstructMatchParameters &params) const;

 public:
  // Factory to build TautomerQuery
  static boost::shared_ptr<TautomerQuery> fromMol(
      const ROMol &molecule,
      const std::string &tautomerTransformFile = std::string());

  // Substructure search
  std::vector<MatchVectType> SubstructMatch(
      const ROMol &mol, const SubstructMatchParameters &params,
      std::vector<ROMOL_SPTR> *matchingTautomers = nullptr) const;

  // Query fingerprint
  ExplicitBitVect *patternFingerprintTemplate(uint fpSize = 2048U);

  // accessors

  const ROMOL_SPTR getTemplateMolecule() const { return templateMolecule; }

  const std::vector<ROMOL_SPTR> getTautomers() const { return tautomers; }

  ~TautomerQuery(){};
};

// so we can use the templates in Code/GraphMol/Substruct/SubstructMatch.h
std::vector<MatchVectType> SubstructMatch(
    const ROMol &mol, const TautomerQuery &query,
    const SubstructMatchParameters &params);

}  // namespace RDKit

#endif  // RDKIT_TAUTOMERQUERY_H
