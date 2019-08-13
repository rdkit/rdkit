//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
/*! \file Charge.h

        \brief Defines the Reionizer class and Uncharger class.

*/
#include <RDGeneral/export.h>
#ifndef RD_CHARGE_H
#define RD_CHARGE_H

#include "MolStandardize.h"
#include <Catalogs/Catalog.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogEntry.h>
#include <GraphMol/MolStandardize/AcidBaseCatalog/AcidBaseCatalogParams.h>

namespace RDKit {
class RWMol;
class ROMol;

namespace MolStandardize {

RDKIT_MOLSTANDARDIZE_EXPORT extern const CleanupParameters
    defaultCleanupParameters;

typedef RDCatalog::HierarchCatalog<AcidBaseCatalogEntry, AcidBaseCatalogParams,
                                   int>
    AcidBaseCatalog;

struct RDKIT_MOLSTANDARDIZE_EXPORT ChargeCorrection {
  std::string Name;
  std::string Smarts;
  int Charge;

  ChargeCorrection(std::string name, std::string smarts, int charge)
      : Name(name), Smarts(smarts), Charge(charge) {}
};

// The default list of ChargeCorrections.
RDKIT_MOLSTANDARDIZE_EXPORT extern std::vector<ChargeCorrection>
    CHARGE_CORRECTIONS;

//! The reionizer class to fix charges and reionize a molecule such that the
// strongest acids ionize first.
/*!

  <b>Notes:</b>
    -
*/

class RDKIT_MOLSTANDARDIZE_EXPORT Reionizer {
 public:
  Reionizer();
  //! construct a Reionizer with a particular acidbaseFile
  Reionizer(const std::string acidbaseFile);
  //! construct a Reionizer with a particular acidbaseFile and charge
  // corrections
  Reionizer(const std::string acidbaseFile,
            const std::vector<ChargeCorrection> ccs);
  //! construct a Reionizer with a particular acidbaseFile and charge
  // corrections
  Reionizer(std::istream &acidbaseStream,
            const std::vector<ChargeCorrection> ccs);
  //! making Reionizer objects non-copyable
  Reionizer(const Reionizer &other) = delete;
  Reionizer &operator=(Reionizer const &) = delete;
  ~Reionizer();

  //! Enforce charges on certain atoms, then perform competitive reionization.
  ROMol *reionize(const ROMol &mol);

 private:
  AcidBaseCatalog *d_abcat;
  std::vector<ChargeCorrection> d_ccs;

  std::pair<unsigned int, std::vector<unsigned int>> *strongestProtonated(
      const ROMol &mol,
      const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &abpairs);
  std::pair<unsigned int, std::vector<unsigned int>> *weakestIonized(
      const ROMol &mol,
      const std::vector<std::pair<ROMOL_SPTR, ROMOL_SPTR>> &abpairs);

};  // Reionizer class

//! The Uncharger class for neutralizing ionized acids and bases.
/*!

  <b>Notes:</b>
    - This class uncharges molecules by adding and/or removing hydrogens.
          - For zwitterions, hydrogens are moved to eliminate charges where
  possible.
          - In cases where there is a positive charge that is not neutralizable,
                an	attempt is made to also preserve the corresponding
  negative charge.

*/

class RDKIT_MOLSTANDARDIZE_EXPORT Uncharger {
 public:
  Uncharger();
  Uncharger(bool canonicalOrdering) : Uncharger() {
    df_canonicalOrdering = canonicalOrdering;
  };
  Uncharger(const Uncharger &other);
  ~Uncharger();

  ROMol *uncharge(const ROMol &mol);

 private:
  bool df_canonicalOrdering = true;
  std::shared_ptr<ROMol> pos_h;
  std::shared_ptr<ROMol> pos_noh;
  std::shared_ptr<ROMol> neg;
  std::shared_ptr<ROMol> neg_acid;
};  // Uncharger class

}  // namespace MolStandardize
}  // namespace RDKit
#endif
