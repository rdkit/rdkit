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
#include <RDGeneral/JSONHelpers.h>
#include <RDGeneral/BoostStartInclude.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <RDGeneral/BoostEndInclude.h>
#include "SmilesJSONParsers.h"

namespace RDKit {

void updateSmilesWriteParamsFromJSON(SmilesWriteParams &params,
                                     const char *details_json) {
  if (details_json && strlen(details_json)) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    params.doIsomericSmiles =
        pt.get("doIsomericSmiles", params.doIsomericSmiles);
    params.doKekule = pt.get("doKekule", params.doKekule);
    params.rootedAtAtom = pt.get("rootedAtAtom", params.rootedAtAtom);
    params.canonical = pt.get("canonical", params.canonical);
    params.allBondsExplicit =
        pt.get("allBondsExplicit", params.allBondsExplicit);
    params.allHsExplicit = pt.get("allHsExplicit", params.allHsExplicit);
    params.doRandom = pt.get("doRandom", params.doRandom);
  }
}

void updateSmilesWriteParamsFromJSON(SmilesWriteParams &params,
                                     const std::string &details_json) {
  updateSmilesWriteParamsFromJSON(params, details_json.c_str());
}

void updateCXSmilesFieldsFromJSON(std::uint32_t &cxSmilesFields,
                                  unsigned int &restoreBondDirs,
                                  const char *details_json) {
  if (details_json && strlen(details_json)) {
    boost::property_tree::ptree pt;
    std::istringstream ss;
    ss.str(details_json);
    boost::property_tree::read_json(ss, pt);
    bool haveCXSmilesFields = false;
    auto cxSmilesFieldsFromJson = flagsFromJson<SmilesWrite::CXSmilesFields>(pt, &haveCXSmilesFields);
    if (haveCXSmilesFields) {
      cxSmilesFields = cxSmilesFieldsFromJson;
    }
    std::string restoreBondDirOption;
    restoreBondDirOption = pt.get("restoreBondDirOption", restoreBondDirOption);
    if (RestoreBondDirOption::_is_valid(restoreBondDirOption.c_str())) {
      restoreBondDirs =
          RestoreBondDirOption::_from_string(restoreBondDirOption.c_str());
    }
  }
}

void updateCXSmilesFieldsFromJSON(std::uint32_t &cxSmilesFields,
                                  unsigned int &restoreBondDirs,
                                  const std::string &details_json) {
  updateCXSmilesFieldsFromJSON(cxSmilesFields, restoreBondDirs,
                               details_json.c_str());
}

}  // end namespace RDKit
