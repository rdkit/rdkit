/*--------------------------------------------------
 * gene.cc
 * Random generation of isomorphic ARGraphs
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
 * without attributes).
 * IMPORTANT:
 *   You have to init the random number generator by calling
 *   srand() before using this function
 -------------------------------------------------------------------*/
void Generate(int nodes, int edges, Graph **g1, Graph **g2,
              bool connected)
  { 
    int i, j, n1, n2;
    byte **mat;
    short *s;

    assert(nodes>0);
    assert(edges>=0 && edges<=nodes*(nodes-1));

    mat=new byte*[nodes];
    for(i=0; i<nodes; i++)
      mat[i]=new byte[nodes];
    s=new short[nodes];
     

    for(i=0; i<nodes; i++)
      memset(mat[i], '\0', nodes*sizeof(mat[0][0]));

       /* Generate the matrix */
    for(i=0; i<edges; i++)
          { do {
            n1=(int)floor(rand()/(float)RAND_MAX*nodes);
            n2=(int)floor(rand()/(float)RAND_MAX*nodes);
            } while (n1==n2 || mat[n1][n2]==1);
            mat[n1][n2]=1;
          }

        /* Ensure the graph is connected */
    if (connected)
    { for(i=1; i<nodes; i++)
          { for(j=0; j<i; j++)
              if (mat[i][j] || mat[j][i])
                break;
            if (i==j)
              { j=(int)floor(rand()/(float)RAND_MAX*i);
                if (rand()/(float)RAND_MAX<0.5) mat[i][j]=1;
                else mat[j][i]=1;
              }
          }
    }


        /* permutation of the matrix */
    for(i=0; i<nodes; i++)
         s[i]=i;
    for(i=0; i<nodes; i++)
          { int tmp;
            j=(int)floor(rand()/(float)RAND_MAX*nodes);
            tmp=s[j];
            s[j]=s[i];
            s[i]=tmp;
          }

    /* Generation of the graphs */
    { ARGEdit eg1;

      for(i=0; i<nodes; i++)
        eg1.InsertNode(NULL);
      for(i=0; i<nodes; i++)
        for(j=0; j<nodes; j++)
          if (mat[i][j])
            eg1.InsertEdge(i, j, NULL);
      *g1=new Graph(&eg1);
      if (*g1==NULL)
        error("Out of memory");
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
            
    for(i=0; i<nodes; i++)
      delete[] mat[i];
    delete[] mat;
    delete[] s;
     
  }

