/*--------------------------------------------------------
 * xsubgraph.h
 * Interface of xsubgraph.cc
 * Random extraction of a (possibly) connected subgraph
 * See: argraph.h
 *
 * Author: P. Foggia
 --------------------------------------------------------*/


#ifndef XSUBGRAPH_H


#include "argraph.h"



Graph* ExtractSubgraph(Graph *g, int nodes, bool connected=true);

#endif
