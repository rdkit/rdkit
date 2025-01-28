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
#include <RDGeneral/RDThreads.h>

#ifdef RDK_BUILD_THREADSAFE_SSS
#include <thread>
#endif

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

namespace {
template <typename FuncType>
void standardizeMultipleMolsInPlace(FuncType sfunc, std::vector<RWMol *> &mols,
                                    int numThreads,
                                    const CleanupParameters &params) {
  unsigned int numThreadsToUse = std::min(
      static_cast<unsigned int>(mols.size()), getNumThreadsToUse(numThreads));
  if (numThreadsToUse == 1) {
    for (auto molp : mols) {
      sfunc(*molp, params);
    }
  }
#ifdef RDK_BUILD_THREADSAFE_SSS
  else {
    auto func = [&](unsigned int tidx) {
      for (auto mi = tidx; mi < mols.size(); mi += numThreads) {
        sfunc(*mols[mi], params);
      }
    };
    std::vector<std::thread> threads;
    for (auto tidx = 0u; tidx < numThreadsToUse; ++tidx) {
      threads.emplace_back(func, tidx);
    }
    for (auto &t : threads) {
      if (t.joinable()) {
        t.join();
      }
    }
  }
#endif
}

void throwIfMolPtrListContainsDuplicates(const std::vector<RWMol *> &mols) {
  // we could do this with an unordered set, but that requires memory allocation
  //   and in the "normal" case where all elements are unique we *will* have to
  //   insert all of them.
  // This way is O(N^2) instead of O(NlogN) - actually closer to O(N) with an
  //   unordered_set - but we're doing essentially no work inside the loop.
  // And, when you get down to it, this code is going to be a vanishingly small
  //   part of the runtime of any real standardization function, even for large
  //   N
  for (auto i = 1u; i < mols.size(); ++i) {
    for (auto j = 0u; j < i; ++j) {
      if (mols[i] == mols[j]) {
        throw ValueErrorException("duplicate molecule in input list");
      }
    }
  }
}
}  // namespace

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

void cleanupInPlace(std::vector<RWMol *> &mols, int numThreads,
                    const CleanupParameters &params) {
  throwIfMolPtrListContainsDuplicates(mols);
  standardizeMultipleMolsInPlace(
      static_cast<void (*)(RWMol &, const CleanupParameters &)>(cleanupInPlace),
      mols, numThreads, params);
}

void tautomerParentInPlace(RWMol &mol, const CleanupParameters &params,
                           bool skip_standardize) {
  if (!skip_standardize) {
    cleanupInPlace(mol, params);
  }

  canonicalTautomerInPlace(mol, params);
  cleanupInPlace(mol, params);
}

void tautomerParentInPlace(std::vector<RWMol *> &mols, int numThreads,
                           const CleanupParameters &params,
                           bool skip_standardize) {
  throwIfMolPtrListContainsDuplicates(mols);
  auto sfunc = [skip_standardize](RWMol &m, const CleanupParameters &ps) {
    tautomerParentInPlace(m, ps, skip_standardize);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}

RWMol *tautomerParent(const RWMol &mol, const CleanupParameters &params,
                      bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  tautomerParentInPlace(*res, params, skip_standardize);
  return res.release();
}

void fragmentParentInPlace(std::vector<RWMol *> &mols, int numThreads,
                           const CleanupParameters &params,
                           bool skip_standardize) {
  throwIfMolPtrListContainsDuplicates(mols);
  auto sfunc = [skip_standardize](RWMol &m, const CleanupParameters &ps) {
    fragmentParentInPlace(m, ps, skip_standardize);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}
void fragmentParentInPlace(RWMol &mol, const CleanupParameters &params,
                           bool skip_standardize) {
  if (!skip_standardize) {
    cleanupInPlace(mol, params);
  }
  LargestFragmentChooser lfragchooser(params.preferOrganic);
  lfragchooser.chooseInPlace(mol);
}

// Return the fragment parent of a given molecule.
// The fragment parent is the largest organic covalent unit in the molecule.
//
RWMol *fragmentParent(const RWMol &mol, const CleanupParameters &params,
                      bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  fragmentParentInPlace(*res, params, skip_standardize);
  return res.release();
}

void stereoParentInPlace(std::vector<RWMol *> &mols, int numThreads,
                         const CleanupParameters &params,
                         bool skip_standardize) {
  throwIfMolPtrListContainsDuplicates(mols);
  auto sfunc = [skip_standardize](RWMol &m, const CleanupParameters &ps) {
    stereoParentInPlace(m, ps, skip_standardize);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}
void stereoParentInPlace(RWMol &mol, const CleanupParameters &params,
                         bool skip_standardize) {
  if (!skip_standardize) {
    cleanupInPlace(mol, params);
  }

  MolOps::removeStereochemistry(mol);
}
RWMol *stereoParent(const RWMol &mol, const CleanupParameters &params,
                    bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  stereoParentInPlace(*res, params, skip_standardize);
  return res.release();
}
void isotopeParentInPlace(std::vector<RWMol *> &mols, int numThreads,
                          const CleanupParameters &params,
                          bool skip_standardize) {
  throwIfMolPtrListContainsDuplicates(mols);
  auto sfunc = [skip_standardize](RWMol &m, const CleanupParameters &ps) {
    isotopeParentInPlace(m, ps, skip_standardize);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}

void isotopeParentInPlace(RWMol &mol, const CleanupParameters &params,
                          bool skip_standardize) {
  if (!skip_standardize) {
    cleanupInPlace(mol, params);
  }

  for (auto atom : mol.atoms()) {
    atom->setIsotope(0);
  }
}
RWMol *isotopeParent(const RWMol &mol, const CleanupParameters &params,
                     bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  isotopeParentInPlace(*res, params, skip_standardize);
  return res.release();
}

void chargeParentInPlace(std::vector<RWMol *> &mols, int numThreads,
                         const CleanupParameters &params,
                         bool skip_standardize) {
  throwIfMolPtrListContainsDuplicates(mols);
  auto sfunc = [skip_standardize](RWMol &m, const CleanupParameters &ps) {
    chargeParentInPlace(m, ps, skip_standardize);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}
void chargeParentInPlace(RWMol &mol, const CleanupParameters &params,
                         bool skip_standardize) {
  fragmentParentInPlace(mol, params, skip_standardize);
  Uncharger uncharger(params.doCanonical);
  uncharger.unchargeInPlace(mol);
  cleanupInPlace(mol, params);
}
RWMol *chargeParent(const RWMol &mol, const CleanupParameters &params,
                    bool skip_standardize) {
  // Return the charge parent of a given molecule.
  // The charge parent is the uncharged version of the fragment parent.
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  chargeParentInPlace(*res, params, skip_standardize);
  return res.release();
}

void superParentInPlace(RWMol &mol, const CleanupParameters &params,
                        bool skip_standardize) {
  if (!skip_standardize) {
    cleanupInPlace(mol, params);
  }
  // we can skip fragmentParent since the chargeParent takes care of that
  chargeParentInPlace(mol, params, true);
  isotopeParentInPlace(mol, params, true);
  stereoParentInPlace(mol, params, true);
  tautomerParentInPlace(mol, params, true);
  cleanupInPlace(mol, params);
}

void superParentInPlace(std::vector<RWMol *> &mols, int numThreads,
                        const CleanupParameters &params,
                        bool skip_standardize) {
  throwIfMolPtrListContainsDuplicates(mols);
  auto sfunc = [skip_standardize](RWMol &m, const CleanupParameters &ps) {
    superParentInPlace(m, ps, skip_standardize);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}

RWMol *superParent(const RWMol &mol, const CleanupParameters &params,
                   bool skip_standardize) {
  std::unique_ptr<RWMol> res{new RWMol(mol)};
  superParentInPlace(*res, params, skip_standardize);
  return res.release();
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

void normalizeInPlace(std::vector<RWMol *> &mols, int numThreads,
                      const CleanupParameters &params) {
  throwIfMolPtrListContainsDuplicates(mols);
  std::unique_ptr<Normalizer> normalizer{normalizerFromParams(params)};
  auto sfunc = [&normalizer](RWMol &m, const CleanupParameters &) {
    normalizer->normalizeInPlace(m);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}

void reionizeInPlace(RWMol &mol, const CleanupParameters &params) {
  std::unique_ptr<Reionizer> reionizer{reionizerFromParams(params)};
  reionizer->reionizeInPlace(mol);
}
void reionizeInPlace(std::vector<RWMol *> &mols, int numThreads,
                     const CleanupParameters &params) {
  throwIfMolPtrListContainsDuplicates(mols);
  std::unique_ptr<Reionizer> reionizer{reionizerFromParams(params)};
  auto sfunc = [&reionizer](RWMol &m, const CleanupParameters &) {
    reionizer->reionizeInPlace(m);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
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

void removeFragmentsInPlace(std::vector<RWMol *> &mols, int numThreads,
                            const CleanupParameters &params) {
  throwIfMolPtrListContainsDuplicates(mols);
  std::unique_ptr<FragmentRemover> remover{fragmentRemoverFromParams(params)};
  auto sfunc = [&remover](RWMol &m, const CleanupParameters &) {
    remover->removeInPlace(m);
  };
  standardizeMultipleMolsInPlace(sfunc, mols, numThreads, params);
}

RWMol *canonicalTautomer(const RWMol *mol, const CleanupParameters &params) {
  PRECONDITION(mol, "bad molecule");
  std::unique_ptr<TautomerEnumerator> te{tautomerEnumeratorFromParams(params)};
  return static_cast<RWMol *>(te->canonicalize(*mol));
}
void canonicalTautomerInPlace(RWMol &mol, const CleanupParameters &params) {
  std::unique_ptr<TautomerEnumerator> te{tautomerEnumeratorFromParams(params)};
  te->canonicalizeInPlace(mol);
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

}  // namespace MolStandardize
}  // namespace RDKit
