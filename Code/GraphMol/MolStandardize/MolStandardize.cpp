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

#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

using namespace std;
namespace RDKit {
namespace MolStandardize {
const CleanupParameters defaultCleanupParameters;

#define PT_OPT_GET(opt) params.opt = pt.get(#opt, params.opt)
void updateCleanupParamsFromJSON(CleanupParameters &params,
                                 const std::string &json) {
  if (json.empty()) {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(rdbase);
  PT_OPT_GET(normalizations);
  PT_OPT_GET(acidbaseFile);
  PT_OPT_GET(fragmentFile);
  PT_OPT_GET(tautomerTransforms);
  PT_OPT_GET(maxRestarts);
  PT_OPT_GET(preferOrganic);
  PT_OPT_GET(doCanonical);
  PT_OPT_GET(maxTautomers);
  PT_OPT_GET(maxTransforms);
  PT_OPT_GET(tautomerRemoveSp3Stereo);
  PT_OPT_GET(tautomerRemoveBondStereo);
  PT_OPT_GET(tautomerRemoveIsotopicHs);
  PT_OPT_GET(tautomerReassignStereo);
  {
    const auto norm_tfs = pt.get_child_optional("normalizationData");
    if (norm_tfs) {
      for (const auto &entry : *norm_tfs) {
        std::string nm = entry.second.get<std::string>("name", "");
        std::string smarts = entry.second.get<std::string>("smarts", "");
        if (nm.empty() || smarts.empty()) {
          BOOST_LOG(rdWarningLog)
              << " empty transformation name or SMARTS" << std::endl;
          continue;
        }
        params.normalizationData.push_back(std::make_pair(nm, smarts));
      }
    }
  }
  {
    const auto frag_tfs = pt.get_child_optional("fragmentData");
    if (frag_tfs) {
      for (const auto &entry : *frag_tfs) {
        std::string nm = entry.second.get<std::string>("name", "");
        std::string smarts = entry.second.get<std::string>("smarts", "");
        if (nm.empty() || smarts.empty()) {
          BOOST_LOG(rdWarningLog)
              << " empty transformation name or SMARTS" << std::endl;
          continue;
        }
        params.fragmentData.push_back(std::make_pair(nm, smarts));
      }
    }
  }
  {
    const auto ab_data = pt.get_child_optional("acidbaseData");
    if (ab_data) {
      for (const auto &entry : *ab_data) {
        std::string nm = entry.second.get<std::string>("name", "");
        std::string acid = entry.second.get<std::string>("acid", "");
        std::string base = entry.second.get<std::string>("base", "");
        if (nm.empty() || acid.empty() || base.empty()) {
          BOOST_LOG(rdWarningLog)
              << " empty component in acidbaseData" << std::endl;
          continue;
        }
        params.acidbaseData.push_back(std::make_tuple(nm, acid, base));
      }
    }
  }
  {
    const auto taut_data = pt.get_child_optional("tautomerTransformData");
    if (taut_data) {
      for (const auto &entry : *taut_data) {
        std::string nm = entry.second.get<std::string>("name", "");
        std::string smarts = entry.second.get<std::string>("smarts", "");
        std::string bonds = entry.second.get<std::string>("bonds", "");
        std::string charges = entry.second.get<std::string>("charges", "");
        if (nm.empty() || smarts.empty()) {
          BOOST_LOG(rdWarningLog)
              << " empty component in tautomerTransformData" << std::endl;
          continue;
        }
        params.tautomerTransformData.push_back(
            std::make_tuple(nm, smarts, bonds, charges));
      }
    }
  }
}

RWMol *cleanup(const RWMol *mol, const CleanupParameters &params) {
  auto nmol = new RWMol(*mol);
  cleanupInPlace(*nmol, params);
  return nmol;
}
void cleanupInPlace(RWMol &mol, const CleanupParameters &params) {
  MolOps::removeHs(mol);
  MolStandardize::MetalDisconnector md;
  md.disconnectInPlace(mol);
  MolStandardize::normalizeInPlace(mol, params);
  MolStandardize::reionizeInPlace(mol, params);
  bool cleanIt = true;
  bool force = true;
  MolOps::assignStereochemistry(mol, cleanIt, force);
}

RWMol *tautomerParent(const RWMol &mol, const CleanupParameters &params,
                      bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  if (!skip_standardize) {
    cleanupInPlace(*res, params);
  }

  std::unique_ptr<RWMol> ct{canonicalTautomer(res.get(), params)};
  cleanupInPlace(*ct, params);
  return ct.release();
}

// Return the fragment parent of a given molecule.
// The fragment parent is the largest organic covalent unit in the molecule.
//
RWMol *fragmentParent(const RWMol &mol, const CleanupParameters &params,
                      bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  if (!skip_standardize) {
    cleanupInPlace(*res, params);
  }
  LargestFragmentChooser lfragchooser(params.preferOrganic);
  return static_cast<RWMol *>(lfragchooser.choose(*res));
}

RWMol *stereoParent(const RWMol &mol, const CleanupParameters &params,
                    bool skip_standardize) {
  RWMol *res = new RWMol(mol);
  if (!skip_standardize) {
    cleanupInPlace(*res, params);
  }

  MolOps::removeStereochemistry(*res);
  return res;
}

RWMol *isotopeParent(const RWMol &mol, const CleanupParameters &params,
                     bool skip_standardize) {
  RWMol *res = new RWMol(mol);
  if (!skip_standardize) {
    cleanupInPlace(*res, params);
  }

  for (auto atom : res->atoms()) {
    atom->setIsotope(0);
  }
  return res;
}

RWMol *chargeParent(const RWMol &mol, const CleanupParameters &params,
                    bool skip_standardize) {
  // Return the charge parent of a given molecule.
  // The charge parent is the uncharged version of the fragment parent.

  std::unique_ptr<RWMol> fragparent{
      fragmentParent(mol, params, skip_standardize)};

  Uncharger uncharger(params.doCanonical);
  uncharger.unchargeInPlace(*fragparent);
  cleanupInPlace(*fragparent, params);
  return fragparent.release();
}

RWMol *superParent(const RWMol &mol, const CleanupParameters &params,
                   bool skip_standardize) {
  std::unique_ptr<RWMol> res;
  if (!skip_standardize) {
    res.reset(cleanup(mol, params));
  } else {
    res.reset(new RWMol(mol));
  }
  // we can skip fragmentParent since the chargeParent takes care of that
  res.reset(chargeParent(*res, params, true));
  res.reset(isotopeParent(*res, params, true));
  res.reset(stereoParent(*res, params, true));
  res.reset(tautomerParent(*res, params, true));
  return cleanup(*res, params);
}

RWMol *normalize(const RWMol *mol, const CleanupParameters &params) {
  PRECONDITION(mol, "bad molecule");
  std::unique_ptr<Normalizer> normalizer{normalizerFromParams(params)};
  return static_cast<RWMol *>(normalizer->normalize(*mol));
}

RWMol *reionize(const RWMol *mol, const CleanupParameters &params) {
  PRECONDITION(mol, "bad molecule");
  std::unique_ptr<Reionizer> reionizer{reionizerFromParams(params)};
  return static_cast<RWMol *>(reionizer->reionize(*mol));
}

void normalizeInPlace(RWMol &mol, const CleanupParameters &params) {
  std::unique_ptr<Normalizer> normalizer{normalizerFromParams(params)};
  normalizer->normalizeInPlace(mol);
}

void reionizeInPlace(RWMol &mol, const CleanupParameters &params) {
  std::unique_ptr<Reionizer> reionizer{reionizerFromParams(params)};
  reionizer->reionizeInPlace(mol);
}

RWMol *removeFragments(const RWMol *mol, const CleanupParameters &params) {
  PRECONDITION(mol, "bad molecule");
  std::unique_ptr<FragmentRemover> remover{fragmentRemoverFromParams(params)};
  return static_cast<RWMol *>(remover->remove(*mol));
}

void removeFragmentsInPlace(RWMol &mol, const CleanupParameters &params) {
  std::unique_ptr<FragmentRemover> remover{fragmentRemoverFromParams(params)};
  remover->removeInPlace(mol);
}

RWMol *canonicalTautomer(const RWMol *mol, const CleanupParameters &params) {
  PRECONDITION(mol, "bad molecule");
  std::unique_ptr<TautomerEnumerator> te{tautomerEnumeratorFromParams(params)};
  return static_cast<RWMol *>(te->canonicalize(*mol));
}

std::string standardizeSmiles(const std::string &smiles) {
  std::unique_ptr<RWMol> mol{SmilesToMol(smiles, 0, false)};
  if (!mol) {
    std::string message =
        "SMILES Parse Error: syntax error for input: " + smiles;
    throw ValueErrorException(message);
  }

  cleanupInPlace(*mol);
  return MolToSmiles(*mol);
}

std::vector<std::string> enumerateTautomerSmiles(
    const std::string &smiles, const CleanupParameters &params) {
  std::unique_ptr<RWMol> mol(SmilesToMol(smiles, 0, false));
  cleanupInPlace(*mol, params);
  MolOps::sanitizeMol(*mol);

  TautomerEnumerator te(params);

  auto res = te.enumerate(*mol);

  return res.smiles();
}

void disconnectOrganometallics(
    RWMol &mol, RDKit::MolStandardize::MetalDisconnectorOptions mdo) {
  RDKit::MolStandardize::MetalDisconnector md(mdo);
  md.disconnect(mol);
}

ROMol *disconnectOrganometallics(
    const ROMol &mol, RDKit::MolStandardize::MetalDisconnectorOptions mdo) {
  RDKit::MolStandardize::MetalDisconnector md(mdo);
  return md.disconnect(mol);
}

}  // end of namespace MolStandardize
}  // namespace RDKit
