/*--------------------------------------------------
 * gene_mesh.cc
 * Random generation of isomorphic ARGraphs that are
 * almost meshes.
 -------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>

#include "argraph.h"
#include "argedit.h"
#include "error.h"

#include "gene.h"


/*--------------------------------------------------------------------
 * Generate randomly two isomorphic graphs (heap-allocated, and
 * without attributes) with mesh structure.
 * The parameter nodes is rounded to the nearest square number.
 * The parameter sub_nodes allows to generate a graph g1 which is 
 * isomorphic to a subgraph of g2, specifying the desired number
 * of nodes in g1. If sub_nodes<=0, the parameter is ignored.
 *
 * IMPORTANT
 *   You have to init the random number generator by calling srand()
 *   before invoking this function.
 -------------------------------------------------------------------*/
void GenerateMesh(int nodes, int extra_edges, Graph **g1, Graph **g2, 
                  int sub_nodes)
  { 
    int i, j, n1, n2;
    byte **mat;
    short *s;
    int sqnodes=(int)(sqrt((double)nodes)+0.5);
    nodes=sqnodes*sqnodes;
    int edges=extra_edges+ 2*sqnodes*(sqnodes-1);

    assert(nodes>0);
    assert(extra_edges>=0 && edges<=nodes*(nodes-1));
    assert(sub_nodes<=nodes);

    mat=new byte*[nodes];
    for(i=0; i<nodes;i++)
      mat[i]=new byte[nodes];
    s=new short[nodes];

    if (sub_nodes<=0)
      sub_nodes=nodes;

    for(i=0; i<nodes; i++)
      memset(mat[i], '\0', nodes*sizeof(mat[0][0]));

       /* Generate the matrix */
    for(i=0; i<sqnodes; i++)
      for(j=0; j<sqnodes; j++)
        { n1=i*sqnodes+j;
          if (j<sqnodes-1)
            mat[n1][n1+1]=1;
          if (i<sqnodes-1)
            mat[n1][n1+sqnodes]=1;
        }

    for(i=0; i<extra_edges; i++)
          { do {
            n1=(int)floor(rand()/(float)RAND_MAX*nodes);
            n2=(int)floor(rand()/(float)RAND_MAX*nodes);
            } while (n1==n2 || mat[n1][n2]==1);
            mat[n1][n2]=1;
          }



        /* permutation of the matrix */
    for(i=0; i<nodes; i++)
         s[i]=i;
    for(i=0; i<nodes; i++)
          { int tmp;
            j=(int)floor(rand()/(float)RAND_MAX*nodes);
            if (j>=nodes)
              continue;
            tmp=s[j];
            s[j]=s[i];
            s[i]=tmp;
          }
    for(i=0; i<nodes; i++)
          { int tmp;
            j=(int)floor(rand()/(float)RAND_MAX*nodes);
            if (j>=nodes)
              continue;
            tmp=s[j];
            s[j]=s[i];
            s[i]=tmp;
          }

    /* If a strict subgraph is wanted, adjust the permutation */
    if (sub_nodes<nodes)
      { int p,q,r;
        int tmp;

        for(p=1; p<sub_nodes; p++)
          { for(q=p; q<nodes; q++)
              { int connected=0;
                for(r=0; r<p && !connected; r++)
                  { if (mat[s[r]][s[q]] ||
                        mat[s[q]][s[r]])
                      connected=1;
                  }
                if (connected)
                  { tmp=s[p];
                    s[p]=s[q];
                    s[q]=tmp;
                    break;
                  }
              }
            assert(q<nodes);
          }
      }

    /* Generation of the graphs */
    { ARGEdit eg1;

      for(i=0; i<sub_nodes; i++)
        eg1.InsertNode(NULL);
      for(i=0; i<sub_nodes; i++)
        for(j=0; j<sub_nodes; j++)
          if (mat[s[i]][s[j]])
            eg1.InsertEdge(i, j, NULL);
      *g1=new Graph(&eg1);
      if (*g1==NULL)
        error("Out of memory");
    }

        /* permutation of the matrix */
    for(i=0; i<nodes; i++)
         s[i]=i;
    for(i=0; i<nodes; i++)
          { int tmp;
            j=(int)floor(rand()/(float)RAND_MAX*nodes);
            if (j>=nodes)
              continue;
            tmp=s[j];
            s[j]=s[i];
            s[i]=tmp;
          }


    { ARGEdit eg2;

      for(i=0; i<nodes; i++)
        eg2.InsertNode(NULL);
      for(i=0; i<nodes; i++)
        for(j=0; j<nodes; j++)
          if (mat[s[i]][s[j]])
            eg2.InsertEdge(i, j, NULL);
      *g2=new Graph(&eg2);
      if (*g2==NULL)
        error("Out of memory");
    }


   for(i=0; i<nodes;i++)
     delete[] mat[i];
   delete[] mat;
   delete[] s;
            
  }

