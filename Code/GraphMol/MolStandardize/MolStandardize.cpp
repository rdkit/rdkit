//
//  Copyright (C) 2018 Susan H. Leung
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "MolStandardize.h"
#include "Metal.h"
#include "Normalize.h"
#include "Tautomer.h"
#include "Fragment.h"
#include <GraphMol/RDKitBase.h>
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>
#include "Charge.h"

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
using namespace std;

namespace RDKit {
namespace MolStandardize {
const CleanupParameters defaultCleanupParameters;

RWMol *cleanup(const RWMol &mol, const CleanupParameters &params) {
  RWMol m(mol);
  MolOps::sanitizeMol(m);
  MolOps::removeHs(m);

  MolStandardize::MetalDisconnector md;
  md.disconnect(m);
  RWMOL_SPTR normalized(MolStandardize::normalize(&m, params));
  RWMol *reionized = MolStandardize::reionize(normalized.get(), params);
  MolOps::assignStereochemistry(*reionized);

  return reionized;
}

void tautomerParent(RWMol &mol, const CleanupParameters &params) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(params);
  UNDER_CONSTRUCTION("Not yet implmented");
}

// Return the fragment parent of a given molecule.
// The fragment parent is the largest organic covalent unit in the molecule.
//
RWMol *fragmentParent(const RWMol &mol, const CleanupParameters &params,
                      bool skip_standardize) {
  const RWMol *cleaned = nullptr;

  if (!skip_standardize) {
    cleaned = cleanup(mol, params);
  } else {
    cleaned = &mol;
  }

  LargestFragmentChooser lfragchooser(params.preferOrganic);
  ROMol nm(*cleaned);
  ROMOL_SPTR lfrag(lfragchooser.choose(nm));
  delete cleaned;
  return new RWMol(*lfrag);
}

void stereoParent(RWMol &mol, const CleanupParameters &params) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(params);
  UNDER_CONSTRUCTION("Not yet implmented");
}

void isotopeParent(RWMol &mol, const CleanupParameters &params) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(params);
  UNDER_CONSTRUCTION("Not yet implmented");
}

RWMol *chargeParent(const RWMol &mol, const CleanupParameters &params,
                    bool skip_standardize) {
  // Return the charge parent of a given molecule.
  // The charge parent is the uncharged version of the fragment parent.

  RWMol *m = nullptr;

  if (!skip_standardize) {
    m = cleanup(mol, params);
  }

  RWMOL_SPTR fragparent(fragmentParent(*m, params, true));

  // if fragment...
  ROMol nm(*fragparent);

  Uncharger uncharger;
  ROMOL_SPTR uncharged(uncharger.uncharge(nm));
  RWMol *omol = cleanup(static_cast<RWMol>(*uncharged), params);
  return omol;
}

void superParent(RWMol &mol, const CleanupParameters &params) {
  RDUNUSED_PARAM(mol);
  RDUNUSED_PARAM(params);
  UNDER_CONSTRUCTION("Not yet implmented");
}

RWMol *normalize(const RWMol *mol, const CleanupParameters &params) {
  Normalizer normalizer(params.normalizations, params.maxRestarts);

  ROMol m(*mol);
  ROMol *normalized = normalizer.normalize(m);

  return static_cast<RWMol *>(normalized);
}

RWMol *reionize(const RWMol *mol, const CleanupParameters &params) {
  RDUNUSED_PARAM(params);
  Reionizer reionizer;
  ROMol m(*mol);
  ROMol *reionized = reionizer.reionize(m);

  return static_cast<RWMol *>(reionized);
}

std::string standardizeSmiles(const std::string &smiles) {
  RWMOL_SPTR mol(SmilesToMol(smiles, 0, false));
  if (!mol) {
    std::string message =
        "SMILES Parse Error: syntax error for input: " + smiles;
    throw ValueErrorException(message);
  }

  CleanupParameters params;
  RWMOL_SPTR cleaned(cleanup(*mol, params));
  return MolToSmiles(*cleaned);
}

std::vector<std::string> enumerateTautomerSmiles(
    const std::string &smiles, const CleanupParameters &params) {
  std::shared_ptr<RWMol> mol(SmilesToMol(smiles, 0, false));
  cleanup(*mol, params);
  MolOps::sanitizeMol(*mol);

  auto *tautparams = new TautomerCatalogParams(params.tautomerTransforms);
  //	unsigned int ntautomers = tautparams->getNumTautomers();
  TautomerCatalog tautcat(tautparams);
  TautomerEnumerator te;

  std::vector<ROMOL_SPTR> res =
      te.enumerate(static_cast<ROMol>(*mol), &tautcat);

  std::vector<std::string> tsmiles;
  for (const auto &r : res) {
    tsmiles.push_back(MolToSmiles(*r));
  }

  return tsmiles;
}

}  // end of namespace MolStandardize
}  // end of namespace RDKit
