//  Copyright (c) 2025, Guillaume Godin Osmo Labs, PBC’s and others
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

#include "Osmordred.h"
#include "OsmordredHelpers.h"

#include <boost/functional/hash.hpp>  // For custom hashing of pairs
#include <RDGeneral/RDThreads.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RDKitBase.h>
#include <GraphMol/Substruct/SubstructMatch.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include <GraphMol/Descriptors/BCUT.h>
#include <GraphMol/Descriptors/MolDescriptors.h>
#include <GraphMol/PeriodicTable.h>
#include <GraphMol/Subgraphs/Subgraphs.h>
#include <GraphMol/Subgraphs/SubgraphUtils.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Atom.h>
#include <RDGeneral/export.h>
#include <RDGeneral/types.h>

#include <boost/graph/adjacency_list.hpp>

#include <set>
#include <cmath>  // For M_PI and pow
#include <tuple>
#include <map>
#include <string>
#include <utility>  // for std::pair
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include <future>
#include <chrono>

namespace RDKit {
namespace Descriptors {
namespace Osmordred {
// Fast aggregate: compute all descriptors in C++ in one pass
std::vector<double> calcOsmordred(const ROMol &mol) {
  // Silence RDKit warnings locally
  RDLog::LogStateSetter guard;

  std::vector<double> out;
  out.reserve(4200);  // approximate capacity to avoid many reallocations

  // Always compute the full v2 set; keep signature for ABI stability
  const bool doExEstate = true;

  // Precompute/cached intermediates where safe
  // Note: We do not change algorithms; just reuse intermediates across calls
  std::unique_ptr<RWMol> kekulizedMol(new RWMol(mol));
  try {
    MolOps::Kekulize(*kekulizedMol, false);
  } catch (...) {
    // leave kekulizedMol as-is when kekulization fails
  }

  // No shared matrix caching in baseline fast path

  // Some families already build needed matrices internally; where public
  // helpers exist (e.g., Adj/Dist matrices), we call the descriptor that
  // accepts version flags so we do not duplicate logic.

  // Collect results with minimal inserts
  auto append = [&out](const std::vector<double> &v) {
    out.insert(out.end(), v.begin(), v.end());
  };
  auto appendInt = [&out](const std::vector<int> &v) {
    out.reserve(out.size() + v.size());
    for (int x : v) out.push_back(static_cast<double>(x));
  };

  append(calcABCIndex(mol));
  appendInt(calcAcidBase(mol));
  append(calcAdjMatrixDescsL(mol));
  appendInt(calcAromatic(mol));
  appendInt(calcAtomCounts(mol));
  append(calcAutoCorrelation(mol));
  append(calcBCUTs(mol));
  append(calcBalabanJ(mol));
  // Use L variant (faster), to match Python CalcBaryszMatrix
  append(calcBaryszMatrixDescsL(mol));
  append(calcBertzCT(mol));
  appendInt(calcBondCounts(mol));
  append(calcRNCG_RPCG(mol));
  append(calcCarbonTypes(mol));
  append(calcAllChiDescriptors(mol));
  append(calcConstitutional(mol));
  append(calcDetourMatrixDescsL(mol));
  append(calcDistMatrixDescsL(mol));
  append(calcEStateDescs(mol, doExEstate));
  append(calcEccentricConnectivityIndex(mol));
  append(calcExtendedTopochemicalAtom(mol));
  append(calcFragmentComplexity(mol));
  append(calcFramework(mol));
  append(calcHydrogenBond(mol));
  append(calcLogS(mol));
  append(calcInformationContent(mol, 5));
  append(calcKappaShapeIndex(mol));
  appendInt(calcLipinskiGhose(mol));
  append(calcMcGowanVolume(mol));
  append(calcMoeType(mol));
  append(calcMolecularDistanceEdgeDescs(mol));
  append(calcMolecularId(mol));
  append(calcPathCount(mol));
  append(calcPolarizability(mol));
  appendInt(calcRingDescriptors(mol));
  append(calcRotatableBond(mol));
  append(calcSLogP(mol));
  append(calcTopoPSA(mol));
  append(calcTopologicalChargeDescs(mol));
  append(calcTopologicalIndex(mol));
  append(calcVdwVolumeABC(mol));
  append(calcVertexAdjacencyInformation(mol));
  append(calcWalkCounts(mol));
  append(calcWeight(mol));
  appendInt(calcWienerIndex(mol));
  append(calcZagrebIndex(mol));
  append(calcPol(mol));
  append(calcMR(mol));
  append(calcFlexibility(mol));
  append(calcSchultz(mol));
  append(calcAlphaKappaShapeIndex(mol));
  append(calcHEStateDescs(mol));
  append(calcBEStateDescs(mol));
  append(calcAbrahams(mol));
  append(calcANMat(mol));
  append(calcASMat(mol));
  append(calcAZMat(mol));
  append(calcDSMat(mol));
  append(calcDN2Mat(mol));
  append(calcFrags(mol));
  append(calcAddFeatures(mol));

  return out;
}

// v2.0: Timeout constant for Osmordred computation (60 seconds = 1 minute)
constexpr int OSMORDRED_TIMEOUT_SECONDS = 60;

// v2.0: Single molecule with timeout protection (all-or-nothing)
// Returns NaN vector if computation exceeds timeout_seconds
std::vector<double> calcOsmordredWithTimeout(const ROMol &mol,
                                             int timeout_seconds) {
  auto future =
      std::async(std::launch::async, [&mol]() { return calcOsmordred(mol); });

  int actual_timeout =
      timeout_seconds > 0 ? timeout_seconds : OSMORDRED_TIMEOUT_SECONDS;
  auto status = future.wait_for(std::chrono::seconds(actual_timeout));

  if (status == std::future_status::ready) {
    return future.get();
  } else {
    // Timeout - return NaN vector (3585 NaN values)
    return std::vector<double>(3585, std::numeric_limits<double>::quiet_NaN());
  }
}

// v2.0: Batch version from SMILES: parses each SMILES -> NEW mol (tautomer
// canonical LOST). For tautomer-canonical mols use
// calcOsmordredBatchFromMols(mols) with mols from ToBinary.
std::vector<std::vector<double>> calcOsmordredBatch(
    const std::vector<std::string> &smiles_list, int n_jobs) {
  std::vector<std::vector<double>> results;
  results.reserve(smiles_list.size());

  unsigned int nThreads = getNumThreadsToUse(n_jobs);

  if (nThreads <= 1 || smiles_list.size() < 10) {
    for (const auto &smi : smiles_list) {
      auto future = std::async(std::launch::async, [&smi]() {
        ROMol *mol = SmilesToMol(smi);
        if (mol) {
          auto desc = calcOsmordred(*mol);
          delete mol;
          return desc;
        }
        return std::vector<double>();
      });

      auto status =
          future.wait_for(std::chrono::seconds(OSMORDRED_TIMEOUT_SECONDS));
      if (status == std::future_status::ready) {
        results.push_back(future.get());
      } else {
        results.push_back(std::vector<double>(
            3585, std::numeric_limits<double>::quiet_NaN()));
      }
    }
    return results;
  }

  // Parallel processing using std::async with timeout
  std::vector<std::future<std::vector<double>>> futures;
  futures.reserve(smiles_list.size());

  for (size_t idx = 0; idx < smiles_list.size(); ++idx) {
    const auto &smi = smiles_list[idx];

    futures.emplace_back(std::async(std::launch::async, [smi]() {
      try {
        ROMol *mol = SmilesToMol(smi);
        if (mol) {
          try {
            std::vector<double> descriptors = calcOsmordred(*mol);
            delete mol;
            return descriptors;
          } catch (...) {
            delete mol;
            return std::vector<double>();
          }
        } else {
          return std::vector<double>();
        }
      } catch (...) {
        return std::vector<double>();
      }
    }));
  }

  for (auto &f : futures) {
    auto status = f.wait_for(std::chrono::seconds(OSMORDRED_TIMEOUT_SECONDS));
    if (status == std::future_status::ready) {
      results.push_back(f.get());
    } else {
      results.push_back(
          std::vector<double>(3585, std::numeric_limits<double>::quiet_NaN()));
    }
  }

  return results;
}

// v2.0: Batch version from mol objects: PRESERVES tautomer canonical.
// Python binding uses mol.ToBinary() -> MolPickler::molFromPickle -> these
// mols.
std::vector<std::vector<double>> calcOsmordredBatchFromMols(
    const std::vector<const ROMol *> &mols, int n_jobs) {
  std::vector<std::vector<double>> results;
  results.reserve(mols.size());

  unsigned int nThreads = getNumThreadsToUse(n_jobs);
  const size_t nFeatures = 3585;
  const std::vector<double> nanRow(nFeatures,
                                   std::numeric_limits<double>::quiet_NaN());

  if (nThreads <= 1 || mols.size() < 10) {
    for (const ROMol *mol : mols) {
      if (mol) {
        try {
          results.push_back(calcOsmordred(*mol));
        } catch (...) {
          results.push_back(nanRow);
        }
      } else {
        results.push_back(nanRow);
      }
    }
    return results;
  }

  std::vector<std::future<std::vector<double>>> futures;
  futures.reserve(mols.size());
  for (size_t idx = 0; idx < mols.size(); ++idx) {
    const ROMol *mol = mols[idx];
#if defined(__clang__)
    // lambda capture divergence between clang and msvc
    futures.emplace_back(std::async(std::launch::async, [mol]() {
#else
      futures.emplace_back(std::async(std::launch::async, [mol, nFeatures]() {
#endif
      if (mol) {
        try {
          return calcOsmordred(*mol);
        } catch (...) {
          return std::vector<double>(nFeatures,
                                     std::numeric_limits<double>::quiet_NaN());
        }
      }
      return std::vector<double>(nFeatures,
                                 std::numeric_limits<double>::quiet_NaN());
    }));
  }
  for (auto &f : futures) {
    auto status = f.wait_for(std::chrono::seconds(OSMORDRED_TIMEOUT_SECONDS));
    if (status == std::future_status::ready) {
      results.push_back(f.get());
    } else {
      results.push_back(nanRow);
    }
  }
  return results;
}
      

// v2.0: Get descriptor names in the same order as calcOsmordred returns values
std::vector<std::string> getOsmordredDescriptorNames() {
  std::vector<std::string> names;
  names.reserve(3585);

  auto addNames = [&names](const std::string &baseName, int count) {
    if (count == 1) {
      names.push_back(baseName);
    } else {
      for (int i = 1; i <= count; ++i) {
        names.push_back(baseName + "_" + std::to_string(i));
      }
    }
  };

  // Add names in the exact same order as calcOsmordred appends values
  addNames("ABCIndex", 2);
  addNames("AcidBase", 2);
  addNames("AdjacencyMatrix", 12);
  addNames("Aromatic", 2);
  addNames("AtomCount", 17);
  addNames("Autocorrelation", 606);
  addNames("BCUT", 24);
  addNames("BalabanJ", 1);
  addNames("BaryszMatrix", 104);
  addNames("BertzCT", 1);
  addNames("BondCount", 9);
  addNames("RNCGRPCG", 2);
  addNames("CarbonTypes", 11);
  addNames("Chi", 56);
  addNames("Constitutional", 16);
  addNames("DetourMatrix", 14);
  addNames("DistanceMatrix", 12);
  addNames("EState", 404);
  addNames("EccentricConnectivityIndex", 1);
  addNames("ExtendedTopochemicalAtom", 45);
  addNames("FragmentComplexity", 1);
  addNames("Framework", 1);
  addNames("HydrogenBond", 2);
  addNames("LogS", 1);
  addNames("InformationContent", 42);
  addNames("KappaShapeIndex", 3);
  addNames("Lipinski", 2);
  addNames("McGowanVolume", 1);
  addNames("MoeType", 54);
  addNames("MolecularDistanceEdge", 19);
  addNames("MolecularId", 12);
  addNames("PathCount", 21);
  addNames("Polarizability", 2);
  addNames("RingCount", 138);
  addNames("RotatableBond", 2);
  addNames("SLogP", 2);
  addNames("TopoPSA", 2);
  addNames("TopologicalCharge", 21);
  addNames("TopologicalIndex", 4);
  addNames("VdwVolumeABC", 1);
  addNames("VertexAdjacencyInformation", 1);
  addNames("WalkCount", 21);
  addNames("Weight", 2);
  addNames("WienerIndex", 2);
  addNames("ZagrebIndex", 4);
  addNames("Pol", 1);
  addNames("MR", 1);
  addNames("Flexibility", 1);
  addNames("Schultz", 1);
  addNames("AlphaKappaShapeIndex", 3);
  addNames("HEState", 88);
  addNames("BEState", 1460);
  addNames("Abrahams", 6);
  addNames("ANMat", 25);
  addNames("ASMat", 20);
  addNames("AZMat", 15);
  addNames("DSMat", 20);
  addNames("DN2Mat", 20);
  addNames("Frags", 215);
  addNames("AddFeatures", 7);

  return names;
}

}
}
}  
