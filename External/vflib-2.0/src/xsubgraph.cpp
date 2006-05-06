/*--------------------------------------------------------
 * xsubgraph.cc
 * Random extraction of a (possibly) connected subgraph
 * See: argraph.h, xsubgraph.h
 *
 * Author: P. Foggia
 --------------------------------------------------------*/

#include <iostream>

#include <assert.h>
#include <stdlib.h>
#include <math.h>

#include "argraph.h"
#include "argedit.h"
#include "xsubgraph.h"
#include "error.h"


inline int irand(int from, int to) 
  { return (int)floor(rand()/(double)(RAND_MAX+1.0)*(to-from))+from;
  }

struct visit_param
  { ARGEdit *ed;
    node_id *map;
  };

typedef Graph::param_type param_type;

static void edge_insert_visitor(Graph *g, node_id n1, node_id n2, void *attr,
                                param_type param);

/*----------------------------------------------------------------
 * Extract a random subgraph from a graph 'g' with a given
 * number of nodes ('nodes'). If 'connected' is true, the 
 * subgraph will be connected. Attributes of the subgraph are
 * shared with the original graph.
 * The subgraph will not inherit the node/edge destroy functions
 * nor the node/edge compatibility predicates.
 * 
 * IMPORTANT
 *   You have to init the random seed by calling srand() before
 *   invoking this function.
 ----------------------------------------------------------------*/
Graph* ExtractSubgraph(Graph *g, int nodes, bool connected)
  { assert(g!=NULL);
    assert(nodes>=0); 
    assert(nodes<=g->NodeCount());

    int i,j;
    int n=g->NodeCount();

    node_id *map=(node_id*)calloc(n, sizeof(node_id));
    if (n>0 && map==NULL)
      OUT_OF_MEMORY();
    for(i=0; i<n; i++)
      map[i]=NULL_NODE;


    ARGEdit ed;

    visit_param param;

    for(i=0; i<nodes; i++)
      { // Choose a node which has not yet been used
        node_id id=irand(0, n-1);
	node_id id0=id;
        
	bool found=false;
        
	do {
	  while (map[id]!=NULL_NODE) 
	    { if (++id == n)
	        id=0;
              if (id==id0)
	        { if (i==0 || !connected)
		    CANT_HAPPEN();
		  else
		    FAIL("Cannot extract a connected subgraph");
		}
       	    }

	  // check for the connectedness of the new node
          if (i>0 && connected)
	    { 
	      for(j=0; j<g->OutEdgeCount(id) && !found; j++)
	        { if (map[g->GetOutEdge(id, j)]!=NULL_NODE)
		    found=true;
		}
	      for(j=0; j<g->InEdgeCount(id) && !found; j++)
	        { if (map[g->GetInEdge(id, j)]!=NULL_NODE)
		    found=true;
		}
	      if (!found)
	        { if (++id == n)
	            id=0;
                  if (id==id0)
	            FAIL("Cannot extract a connected subgraph");
		}
	    }
	  else
	    { found=true;
	    }
	} while (!found); 


	// now add the node to the ARGEdit

	map[id]=i;

	ed.InsertNode(g->GetNodeAttr(id));
      }

    // Now add the edges to the new graph
    param.ed=&ed;
    param.map=map;
    for(i=0; i<n; i++)
      { if (map[i]!=NULL_NODE)
	  g->VisitOutEdges(i, edge_insert_visitor, &param);
      }


    // Construct the graph and return it

    Graph *sub=new Graph(&ed);

    return sub;
  }



/*---------------------------------------------------------
 *  STATIC FUNCTIONS
 --------------------------------------------------------*/
static void edge_insert_visitor(Graph *g, node_id n1, node_id n2, void *attr,
                                param_type param)
  { visit_param *p=(visit_param *)param;
    ARGEdit *ed=p->ed;
    node_id *map=p->map;
    if (map[n1]!=NULL_NODE && map[n2]!=NULL_NODE)
      ed->InsertEdge(map[n1], map[n2], attr);
  }
