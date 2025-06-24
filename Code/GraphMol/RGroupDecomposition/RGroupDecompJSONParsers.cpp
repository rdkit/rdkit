//
//  Copyright (C) 2024 Novartis Biomedical Research and other RDKit contributors
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#define USE_BETTER_ENUMS
#include "RGroupDecompJSONParsers.h"
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

void updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const std::string &details_json) {
  updateRGroupDecompositionParametersFromJSON(params, details_json.c_str());
}

void updateRGroupDecompositionParametersFromJSON(
    RGroupDecompositionParameters &params, const char *details_json) {
  if (details_json && strlen(details_json)) {
    std::istringstream ss;
    boost::property_tree::ptree pt;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);

    std::string labels;
    labels = pt.get<std::string>("labels", labels);
    if (RGroupLabels::_is_valid(labels.c_str())) {
      params.labels = RGroupLabels::_from_string(labels.c_str());
    }

    std::string matchingStrategy;
    matchingStrategy =
        pt.get<std::string>("matchingStrategy", matchingStrategy);
    if (RGroupMatching::_is_valid(matchingStrategy.c_str())) {
      params.matchingStrategy =
          RGroupMatching::_from_string(matchingStrategy.c_str());
    }

    std::string scoreMethod;
    scoreMethod = pt.get<std::string>("scoreMethod", scoreMethod);
    if (RGroupScore::_is_valid(scoreMethod.c_str())) {
      params.scoreMethod = RGroupScore::_from_string(scoreMethod.c_str());
    }

    std::string rgroupLabelling;
    rgroupLabelling = pt.get<std::string>("rgroupLabelling", rgroupLabelling);
    if (RGroupLabelling::_is_valid(rgroupLabelling.c_str())) {
      params.rgroupLabelling =
          RGroupLabelling::_from_string(rgroupLabelling.c_str());
    }

    std::string alignment;
    alignment = pt.get<std::string>("alignment", alignment);
    if (RGroupCoreAlignment::_is_valid(alignment.c_str())) {
      params.alignment = RGroupCoreAlignment::_from_string(alignment.c_str());
    }

    params.chunkSize = pt.get<unsigned int>("chunkSize", params.chunkSize);
    params.onlyMatchAtRGroups =
        pt.get<bool>("onlyMatchAtRGroups", params.onlyMatchAtRGroups);
    params.removeAllHydrogenRGroups = pt.get<bool>(
        "removeAllHydrogenRGroups", params.removeAllHydrogenRGroups);
    params.removeAllHydrogenRGroupsAndLabels =
        pt.get<bool>("removeAllHydrogenRGroupsAndLabels",
                     params.removeAllHydrogenRGroupsAndLabels);
    params.removeHydrogensPostMatch = pt.get<bool>(
        "removeHydrogensPostMatch", params.removeHydrogensPostMatch);
    params.allowNonTerminalRGroups =
        pt.get<bool>("allowNonTerminalRGroups", params.allowNonTerminalRGroups);
    params.allowMultipleRGroupsOnUnlabelled =
        pt.get<bool>("allowMultipleRGroupsOnUnlabelled",
                     params.allowMultipleRGroupsOnUnlabelled);
    params.doTautomers = pt.get<bool>("doTautomers", params.doTautomers);
    params.doEnumeration = pt.get<bool>("doEnumeration", params.doEnumeration);
    params.includeTargetMolInResults = pt.get<bool>(
        "includeTargetMolInResults", params.includeTargetMolInResults);
    params.timeout = pt.get<double>("timeout", params.timeout);
  }
}

}  // end namespace RDKit
