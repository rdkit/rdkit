/*----------------------------------------------------
 * sortnodes.cc
 * A function to sort nodes by frequency of in/out
 * degrees.
 ---------------------------------------------------*/

#include <stdlib.h>
#include "sortnodes.h"

#if 0
// Testing main
#include "argedit.h"
#include <stdio.h>
int main(int argc, char *argv[])
  { int n=argc>1? atoi(argv[1]): 10;
    if (n<1) n=10;
    ARGEdit ed;
    int i,j;
    for(i=0; i<n; i++)
      ed.InsertNode(NULL);
    for(i=0; i<n; i++)
      { for(j=1; j<i-(3-i)%4; j++)
          { ed.InsertEdge(i, (i+j)%n, NULL);
	    printf("(%d-%d) ", i, (i+j)%n);
	  }
      }
	printf("\n");
    Graph g(&ed);
    node_id *s=SortNodesByFrequency(&g);
    for(i=0; i<n; i++)
      printf("%d ", s[i]);
    printf("\n");
    return 0;
  }
#endif

typedef int (*compare_fn)(const void *, const void *);

struct NodeInfo
  { node_id id;
    node_id in;
    node_id out;
  };

static int nodeInfoComp1(NodeInfo *a, NodeInfo *b);
static int nodeInfoComp2(NodeInfo *a, NodeInfo *b);

/*----------------------------------------------------
 * Sorts the nodes of a graphs, returning a 
 * heap-allocated vector (using new) with the node ids
 * in the proper orders.
 * The sorting criterion takes into account:
 *    1 - The number of nodes with the same in/out 
 *        degree.
 *    2 - The valence of the nodes.
 * The nodes at the beginning of the vector are
 * the most singular, from which the matching should
 * start.
 *--------------------------------------------------*/
node_id* SortNodesByFrequency(Graph *g)
  { int n=g->NodeCount();
    NodeInfo *vect=new NodeInfo[n];
    int i;
    for(i=0; i<n; i++)
      { vect[i].id=i;
        vect[i].in=g->InEdgeCount(i);
        vect[i].out=g->OutEdgeCount(i);
      }
    
    //for(i=0; i<n; i++)
    //  printf("<%d/%d/%d> ", vect[i].id, vect[i].in, vect[i].out);
    //printf("\n");

    qsort(vect, n, sizeof(vect[0]), (compare_fn)nodeInfoComp1);

    //for(i=0; i<n; i++)
    //  printf("<%d/%d/%d> ", vect[i].id, vect[i].in, vect[i].out);
    //printf("\n");

    int run=1;

    for(i=0; i<n; i+=run)
      { for(run=1; i+run<n && 
                   vect[i+run].in==vect[i].in && 
		   vect[i+run].out==vect[i].out;
            run++)
	  ;
	int j;
	for(j=0; j<run; j++)
	  { vect[i+j].in += vect[i+j].out;
	    vect[i+j].out=run;
	  }
       }
	

    //for(i=0; i<n; i++)
    //  printf("<%d/%d/%d> ", vect[i].id, vect[i].in, vect[i].out);
    //printf("\n");

    qsort(vect, n, sizeof(vect[0]), (compare_fn)nodeInfoComp2);

    //for(i=0; i<n; i++)
    //  printf("<%d/%d/%d> ", vect[i].id, vect[i].in, vect[i].out);
    //printf("\n");

    node_id *nodes=new node_id[n];
    for(i=0; i<n; i++)
      nodes[i]=vect[i].id;
    
    delete[] vect;

    return nodes;
  }

   




/**
 * The ordering by in/out degree
 */
static int nodeInfoComp1(NodeInfo *a, NodeInfo *b)
  { if (a->out < b->out)
      return -1;
    else if (a->out > b->out)
      return +1;
    else if (a->in < b->in)
      return -1;
    else if (a->in > b->in)
      return +1;
    else 
      return 0;
  }

/**
 * The ordering by frequency/valence.
 * The frequency is in the out field, the valence in `in'.
 */
static int nodeInfoComp2(NodeInfo *a, NodeInfo *b)
  { if (a->in==0 && b->in !=0)
      return +1;
    else if (a->in!=0 && b->in==0)
      return -1;
    else if (a->out < b->out)
      return -1;
    else if (a->out > b->out)
      return +1;
    else if (a->in < b->in)
      return -1;
    else if (a->in > b->in)
      return +1;
    else 
      return 0;
  }
