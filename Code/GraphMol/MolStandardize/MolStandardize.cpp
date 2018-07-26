#include "MolStandardize.h"
#include "Metal.h"
#include "Normalize.h"
#include "Tautomer.h"
#include <GraphMol/RDKitBase.h>
#include <iostream>
#include <GraphMol/ROMol.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/MolStandardize/FragmentCatalog/FragmentRemover.h>
#include <GraphMol/MolStandardize/TransformCatalog/TransformCatalogParams.h>
#include "Charge.h"

#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
using namespace std;

namespace RDKit {
namespace MolStandardize {

RWMol *cleanup(const RWMol &mol, const CleanupParameters &params) {
  //	bool passedOp = false;

  //	auto *newM = new RWMol(mol);
  RWMol m(mol);
  MolOps::sanitizeMol(m);
  //	std::cout << "After sanitizeMol: " << MolToSmiles(mol) << std::endl;
  MolOps::removeHs(m);

  // TODO
  MolStandardize::MetalDisconnector md;
  md.disconnect(m);
  RWMOL_SPTR normalized(MolStandardize::normalize(&m, params));
  RWMol *reionized = MolStandardize::reionize(normalized.get(), params);
  std::cout << "After reionized in cleanup: " << MolToSmiles(*reionized)
            << std::endl;
  MolOps::assignStereochemistry(*reionized);

  return reionized;
}

void tautomerParent(RWMol &mol, const CleanupParameters &params) {}

// Return the fragment parent of a given molecule.
// The fragment parent is the largest organic covalent unit in the molecule.
//
RWMol *fragmentParent(const RWMol &mol, const CleanupParameters &params,
                      bool skip_standardize) {
  const RWMol *cleaned = nullptr;

  if (!skip_standardize) {
    cleaned = cleanup(mol, params);
  } else {
    std::cout << "before copy " << std::endl;
    cleaned = &mol;
    std::cout << "after copy" << MolToSmiles(*cleaned) << std::endl;
  }

  // TODO
  // largest fragment
  LargestFragmentChooser lfragchooser(params.preferOrganic);
  ROMol nm(*cleaned);
  ROMOL_SPTR lfrag = lfragchooser.choose(nm);
  delete cleaned;
  std::cout << "lfrag: " << MolToSmiles(*lfrag) << std::endl;
  return new RWMol(*lfrag);
}

void stereoParent(RWMol &mol, const CleanupParameters &params) {}

void isotopeParent(RWMol &mol, const CleanupParameters &params) {}

RWMol *chargeParent(const RWMol &mol, const CleanupParameters &params,
                    bool skip_standardize) {
  // Return the charge parent of a given molecule.
  // The charge parent is the uncharged version of the fragment parent.

  RWMol *m = nullptr;

  if (!skip_standardize) {
    m = cleanup(mol, params);
    std::cout << "After cleanup in chargeparent: " << MolToSmiles(*m)
              << std::endl;
  }
  std::cout << "After cleanup: " << MolToSmiles(*m) << std::endl;

  RWMOL_SPTR fragparent(fragmentParent(*m, params, true));
  //	delete m;
  std::cout << "After fragmentParent: " << MolToSmiles(*fragparent)
            << std::endl;

  // if fragment...
  ROMol nm(*fragparent);

  Uncharger uncharger;
  ROMOL_SPTR uncharged(uncharger.uncharge(nm));
  std::cout << "After uncharge: " << MolToSmiles(*uncharged) << std::endl;
  RWMol *omol = cleanup(static_cast<RWMol>(*uncharged), params);
  std::cout << "After second cleanup: " << MolToSmiles(*omol) << std::endl;
  //	delete m;
  return omol;
}

void superParent(RWMol &mol, const CleanupParameters &params) {}

RWMol *normalize(const RWMol *mol, const CleanupParameters &params) {
  std::shared_ptr<TransformCatalogParams> tparams(
      new TransformCatalogParams(params.normalizations));
  TransformCatalog tcat(tparams.get());
  Normalizer normalizer;

  ROMol m(*mol);
  ROMol *normalized = normalizer.normalize(m, &tcat);

  std::cout << "normalized: " << MolToSmiles(*normalized) << std::endl;
  //	mol = static_cast<RWMol*>(normalized.get());
  return static_cast<RWMol *>(normalized);
}

RWMol *reionize(const RWMol *mol, const CleanupParameters &params) {
  std::unique_ptr<AcidBaseCatalogParams> abparams(
      new AcidBaseCatalogParams(params.acidbaseFile));
  AcidBaseCatalog abcat(abparams.get());
  Reionizer reionizer;
  ROMol m(*mol);
  ROMol *reionized = reionizer.reionize(m, &abcat);
  std::cout << "After reionizing: " << MolToSmiles(*reionized) << std::endl;
  //	mol = nullptr;

  return static_cast<RWMol *>(reionized);  // new RWMol(*reionized);
}

std::string standardizeSmiles(const std::string &smiles) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smiles, 0, false));
  CleanupParameters params;
  RWMOL_SPTR cleaned(cleanup(*mol, params));
  return MolToSmiles(*cleaned);
}

std::vector<std::string> enumerateTautomerSmiles(
    const std::string &smiles, const CleanupParameters &params) {
  std::shared_ptr<RWMol> mol(SmilesToMol(smiles, 0, false));
  std::cout << "Before cleanup: " << MolToSmiles(*mol) << std::endl;
  cleanup(*mol, params);
  std::cout << "After cleanup: " << MolToSmiles(*mol) << std::endl;
  MolOps::sanitizeMol(*mol);
  std::cout << "After after cleanup: " << MolToSmiles(*mol) << std::endl;

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
