//
//  Copyright (C) 2019 Greg Landrum
//
//   @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include <RingDecomposerLib.h>

// We're just building this library to make dependencies brought in via
// cmake's ExternalProject mechanism work properly
// As of this writing (July 2019) I was unable to get dependencies to
// work correctly without adding this bogus library.

RDL_API void do_nothing(unsigned n) {
  RDL_graph *graph = RDL_initNewGraph(n);
  for (auto i = 0; i < n - 1; ++i) {
    RDL_addUEdge(graph, i, i + 1);
  }
  RDL_data *urfdata = RDL_calculate(graph);
  RDL_deleteGraph(graph);
  delete urfdata;
}
