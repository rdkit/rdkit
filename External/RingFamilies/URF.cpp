#include <RingDecomposerLib.h>

void do_nothing(unsigned n) {
  RDL_graph *graph = RDL_initNewGraph(n);
  for (auto i = 0; i < n - 1; ++i) {
    RDL_addUEdge(graph, i, i + 1);
  }
  RDL_data *urfdata = RDL_calculate(graph);
  RDL_deleteGraph(graph);
  delete urfdata;
}