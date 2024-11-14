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
#include "MolFragmenterJSONParser.h"
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {

void parseMolzipParametersJSON(MolzipParams &params, const char *details_json) {
  if (!details_json || !strlen(details_json)) {
    return;
  }
  boost::property_tree::ptree pt;
  std::istringstream ss;
  ss.str(details_json);
  boost::property_tree::read_json(ss, pt);
  std::string label;
  label = pt.get<std::string>("Label", label);
  if (MolzipLabel::_is_valid(label.c_str())) {
    params.label = MolzipLabel::_from_string(label.c_str());
  }
  auto atomSymbolsIt = pt.find("AtomSymbols");
  if (atomSymbolsIt != pt.not_found()) {
    const auto &jsonVect = atomSymbolsIt->second;
    params.atomSymbols.resize(jsonVect.size());
    std::transform(
        jsonVect.begin(), jsonVect.end(), params.atomSymbols.begin(),
        [](const auto &atomSymbolNode) {
          return atomSymbolNode.second.template get_value<std::string>();
        });
  }

  params.atomProperty =
      pt.get<std::string>("AtomProperty", params.atomProperty);
  params.enforceValenceRules =
      pt.get<bool>("EnforceValenceRules", params.enforceValenceRules);
  params.generateCoordinates =
      pt.get<bool>("GenerateCoordinates", params.generateCoordinates);
}

}  // end namespace RDKit
