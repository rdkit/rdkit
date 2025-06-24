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
#include "JSONParsers.h"
#include <GraphMol/MolPickler.h>
#include <GraphMol/FileParsers/PNGParser.h>
#include <GraphMol/SmilesParse/SmilesJSONParsers.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace MinimalLib {

void updatePropertyPickleOptionsFromJSON(unsigned int &propFlags,
                                         const char *details_json) {
  if (details_json && details_json[0]) {
    std::istringstream ss;
    boost::property_tree::ptree pt;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    const auto nodeIt = pt.find("propertyFlags");
    if (nodeIt != pt.not_found()) {
      auto propertyFlagsFromJson =
          (+PicklerOps::PropertyPickleOptions::NoProps)._to_integral();
      for (const auto *key : PicklerOps::PropertyPickleOptions::_names()) {
        if (nodeIt->second.get(key, false)) {
          propertyFlagsFromJson |=
              PicklerOps::PropertyPickleOptions::_from_string(key)
                  ._to_integral();
        }
      }
      propFlags = propertyFlagsFromJson;
    }
  }
}

void updateSanitizeFlagsFromJSON(unsigned int &sanitizeFlags,
                                 const char *details_json) {
  if (details_json && details_json[0]) {
    std::istringstream ss;
    boost::property_tree::ptree pt;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    auto sanitizeFlagsFromJson =
        (+MolOps::SanitizeFlags::SANITIZE_NONE)._to_integral();
    for (const auto *key : MolOps::SanitizeFlags::_names()) {
      if (pt.get(key, false)) {
        sanitizeFlagsFromJson |=
            MolOps::SanitizeFlags::_from_string(key)._to_integral();
      }
    }
    sanitizeFlags = sanitizeFlagsFromJson;
  }
}

void updateRemoveHsParametersFromJSON(MolOps::RemoveHsParameters &ps,
                                      bool &sanitize,
                                      const char *details_json) {
  if (details_json && details_json[0]) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    ps.removeDegreeZero = pt.get("removeDegreeZero", ps.removeDegreeZero);
    ps.removeHigherDegrees =
        pt.get("removeHigherDegrees", ps.removeHigherDegrees);
    ps.removeOnlyHNeighbors =
        pt.get("removeOnlyHNeighbors", ps.removeOnlyHNeighbors);
    ps.removeIsotopes = pt.get("removeIsotopes", ps.removeIsotopes);
    ps.removeAndTrackIsotopes =
        pt.get("removeAndTrackIsotopes", ps.removeAndTrackIsotopes);
    ps.removeDummyNeighbors =
        pt.get("removeDummyNeighbors", ps.removeDummyNeighbors);
    ps.removeDefiningBondStereo =
        pt.get("removeDefiningBondStereo", ps.removeDefiningBondStereo);
    ps.removeWithWedgedBond =
        pt.get("removeWithWedgedBond", ps.removeWithWedgedBond);
    ps.removeWithQuery = pt.get("removeWithQuery", ps.removeWithQuery);
    ps.removeMapped = pt.get("removeMapped", ps.removeMapped);
    ps.removeInSGroups = pt.get("removeInSGroups", ps.removeInSGroups);
    ps.showWarnings = pt.get("showWarnings", ps.showWarnings);
    ps.removeNonimplicit = pt.get("removeNonimplicit", ps.removeNonimplicit);
    ps.updateExplicitCount =
        pt.get("updateExplicitCount", ps.updateExplicitCount);
    ps.removeHydrides = pt.get("removeHydrides", ps.removeHydrides);
    ps.removeNontetrahedralNeighbors = pt.get("removeNontetrahedralNeighbors",
                                              ps.removeNontetrahedralNeighbors);
    sanitize = pt.get("sanitize", sanitize);
  }
}

}  // end namespace MinimalLib
}  // end namespace RDKit
