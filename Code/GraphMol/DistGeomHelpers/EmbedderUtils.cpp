//
//  Copyright (C) 2021 Greg Landrum
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
#include <RDGeneral/BoostEndInclude.h>

namespace RDKit {
namespace DGeomHelpers {

#define PT_OPT_GET(opt) params.opt = pt.get(#opt, params.opt)

void updateEmbedParametersFromJSON(EmbedParameters &params,
                                   const std::string &json) {
  if (json.empty()) {
    return;
  }
  std::istringstream ss;
  ss.str(json);
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(ss, pt);
  PT_OPT_GET(maxIterations);
  PT_OPT_GET(numThreads);
  PT_OPT_GET(randomSeed);
  PT_OPT_GET(clearConfs);
  PT_OPT_GET(useRandomCoords);
  PT_OPT_GET(ignoreSmoothingFailures);
  PT_OPT_GET(enforceChirality);
  PT_OPT_GET(useExpTorsionAnglePrefs);
  PT_OPT_GET(useBasicKnowledge);
  PT_OPT_GET(pruneRmsThresh);
  PT_OPT_GET(onlyHeavyAtomsForRMS);
  PT_OPT_GET(ETversion);
  PT_OPT_GET(embedFragmentsSeparately);
  PT_OPT_GET(useSmallRingTorsions);
  PT_OPT_GET(useMacrocycleTorsions);
  PT_OPT_GET(useMacrocycle14config);
  PT_OPT_GET(forceTransAmides);
  PT_OPT_GET(useSymmetryForPruning);
}

}  // namespace DGeomHelpers
}  // namespace RDKit