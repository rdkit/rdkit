//
//  Copyright (C) 2021-2026 Greg Landrum and other RDKit contributors.
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "Embedder.h"

#include <RDGeneral/BoostStartInclude.h>
#include <boost/lexical_cast.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <RDGeneral/BoostEndInclude.h>

#include <DistGeom/BoundsMatrix.h>
#include <DistGeom/ChiralSet.h>
#include <Geometry/point.h>
#include <string>

namespace RDKit {
namespace DGeomHelpers {

#define EMBED_PARAMS_FIELDS(X)                    \
  X(basinThresh)                                  \
  X(boundsMatForceScaling)                        \
  X(boxSizeMult)                                  \
  X(clearConfs)                                   \
  X(embedForceField)                              \
  X(embedFragmentsSeparately)                     \
  X(enableSequentialRandomSeeds)                  \
  X(enforceChirality)                             \
  X(ETversion)                                    \
  X(forceTransAmides)                             \
  X(ignoreSmoothingFailures)                      \
  X(maxIterations)                                \
  X(numThreads)                                   \
  X(numZeroFail)                                  \
  X(onlyHeavyAtomsForRMS)                         \
  X(optimizerForceTol)                            \
  X(pruneRmsThresh)                               \
  X(randNegEig)                                   \
  X(randomSeed)                                   \
  X(symmetrizeConjugatedTerminalGroupsForPruning) \
  X(timeout)                                      \
  X(trackFailures)                                \
  X(useBasicKnowledge)                            \
  X(useExpTorsionAnglePrefs)                      \
  X(useLegacyImplementation)                      \
  X(useMacrocycle14config)                        \
  X(useMacrocycleTorsions)                        \
  X(useRandomCoords)                              \
  X(useSmallRingTorsions)                         \
  X(useSymmetryForPruning)                        \
  X(verbose)                                      \
  X(onlyInitialEmbedding)

#define PT_OPT_GET(opt) params.opt = pt.get(#opt, params.opt);
#define PT_OPT_PUT(opt) pt.put(#opt, params.opt);

void updateEmbedParametersFromJSON(EmbedParameters &params,
                                   const std::string &json) {
  if (json.empty()) {
    return;
  }
  std::istringstream ss(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);

  EMBED_PARAMS_FIELDS(PT_OPT_GET)

  std::map<int, RDGeom::Point3D> *cmap = nullptr;
  const auto coordMap = pt.get_child_optional("coordMap");
  if (coordMap) {
    // NOTE: this leaks since EmbedParameters uses a naked pointer and we don't
    // have any way to tie the lifetime of the memory we allocate here to the
    // EmbedParameters object itself.
    cmap = new std::map<int, RDGeom::Point3D>();
    for (const auto &entry : *coordMap) {
      RDGeom::Point3D pt;

      auto itm = entry.second.begin();
      pt.x = itm->second.get_value<float>();
      ++itm;
      pt.y = itm->second.get_value<float>();
      ++itm;
      pt.z = itm->second.get_value<float>();

      (*cmap)[boost::lexical_cast<int>(entry.first)] = pt;
    }
    params.coordMap = cmap;
  }

  if (auto opt = pt.get_child_optional("additionalChiralSets")) {
    for (const auto &entry : *opt) {
      auto chiralSet = std::make_shared<DistGeom::ChiralSet>();
      if (const auto val = entry.second.get_optional<unsigned int>("d_idx0")) {
        chiralSet->d_idx0 = *val;
      }
      if (const auto val = entry.second.get_optional<unsigned int>("d_idx1")) {
        chiralSet->d_idx1 = *val;
      }
      if (const auto val = entry.second.get_optional<unsigned int>("d_idx2")) {
        chiralSet->d_idx2 = *val;
      }
      if (const auto val = entry.second.get_optional<unsigned int>("d_idx3")) {
        chiralSet->d_idx3 = *val;
      }
      if (const auto val = entry.second.get_optional<unsigned int>("d_idx4")) {
        chiralSet->d_idx4 = *val;
      }
      if (const auto val =
              entry.second.get_optional<float>("d_volumeLowerBound")) {
        chiralSet->d_volumeLowerBound = *val;
      }
      if (const auto val =
              entry.second.get_optional<float>("d_volumeUpperBound")) {
        chiralSet->d_volumeUpperBound = *val;
      }
      if (const auto val =
              entry.second.get_optional<std::uint64_t>("d_structureFlags")) {
        chiralSet->d_structureFlags = *val;
      }
      // Here you would parse the entry to populate chiralSet
      // This is a placeholder for the actual parsing logic
      params.additionalChiralSets.push_back(chiralSet);
    }
  }
  if (auto opt = pt.get_optional<std::string>("initialEmbeddingMode")) {
    std::stringstream valS{*opt};
    if (!(valS >> params.initialEmbeddingMode)) {
      throw std::invalid_argument(
          "Initial embedding option must be either 'DG_EMBEDDING', 'INTERNAL_COORDINATE_EMBEDDING', or 'RANDOM_COORDINATE_EMBEDDING'");
    }
  }

  if (auto opt = pt.get_optional<unsigned int>("optimizeConfWithId")) {
    params.optimizeConfWithId = *opt;
  }
}

std::string embedParametersToJSON(const EmbedParameters &params) {
  boost::property_tree::ptree pt;

  EMBED_PARAMS_FIELDS(PT_OPT_PUT)
  PT_OPT_PUT(initialEmbeddingMode)

  if (params.optimizeConfWithId.has_value()) {
    pt.put("optimizeConfWithId",
           std::to_string(params.optimizeConfWithId.value()));
  }

  if (params.coordMap) {
    boost::property_tree::ptree coordMapPT;

    for (const auto &kv : *params.coordMap) {
      boost::property_tree::ptree pointPT;
      pointPT.push_back(
          {"", boost::property_tree::ptree(std::to_string(kv.second.x))});
      pointPT.push_back(
          {"", boost::property_tree::ptree(std::to_string(kv.second.y))});
      pointPT.push_back(
          {"", boost::property_tree::ptree(std::to_string(kv.second.z))});

      coordMapPT.add_child(std::to_string(kv.first), pointPT);
    }

    pt.add_child("coordMap", coordMapPT);
  }

  if (params.boundsMat) {
    boost::property_tree::ptree matrixPT;
    const unsigned int N = params.boundsMat->numCols();
    for (unsigned i = 0; i < N; ++i) {
      boost::property_tree::ptree rowPT;

      for (unsigned j = 0; j < N; ++j) {
        boost::property_tree::ptree v;
        v.put("", params.boundsMat->getVal(i, j));
        rowPT.push_back({"", v});
      }

      matrixPT.push_back({"", rowPT});
    }

    pt.add_child("boundsMatrix", matrixPT);
  }

  if (!params.additionalChiralSets.empty()) {
    boost::property_tree::ptree chiralSetsPT;

    for (const auto &chiralSet : params.additionalChiralSets) {
      boost::property_tree::ptree chiralSetPT;
      chiralSetPT.put("d_idx0", chiralSet->d_idx0);
      chiralSetPT.put("d_idx1", chiralSet->d_idx1);
      chiralSetPT.put("d_idx2", chiralSet->d_idx2);
      chiralSetPT.put("d_idx3", chiralSet->d_idx3);
      chiralSetPT.put("d_idx4", chiralSet->d_idx4);
      chiralSetPT.put("d_volumeLowerBound", chiralSet->d_volumeLowerBound);
      chiralSetPT.put("d_volumeUpperBound", chiralSet->d_volumeUpperBound);
      chiralSetPT.put("d_structureFlags", chiralSet->d_structureFlags);

      chiralSetsPT.push_back({"", chiralSetPT});
    }

    pt.add_child("additionalChiralSets", chiralSetsPT);
  }

  std::ostringstream ss;
  boost::property_tree::write_json(ss, pt, false);
  auto str = ss.str();
  boost::algorithm::trim(str);
  return str;
}

}  // namespace DGeomHelpers
}  // namespace RDKit
